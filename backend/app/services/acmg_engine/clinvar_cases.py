"""ClinVar affected-case counts for PS4 B (proband counting, rulebook ``ps4``).

Reads the ClinVar VCV XML release (``ClinVarVCVRelease_00-latest.xml.gz``) and
keeps, per Variation ID, the submitters that reported the variant in an
affected individual:

* the submission (ClinicalAssertion) is current and not classified benign /
  likely benign by that submitter;
* an ObservedIn sample is germline (not somatic) with ``AffectedStatus`` yes;
* its method is one of ``methods`` (clinical testing, research, provider
  interpretation, case-control by default): literature-only and curation
  submissions repeat published cases and are left out.

PM3 (recessive disorders): in the same pass, submissions (any method, not
benign / likely benign, current) that describe the variant in a homozygous or
compound-heterozygous / in-trans / biallelic state in patients are collected,
from the structured ObservedData counts (``Homozygote``,
``CompoundHeterozygote``) and from the submitter's comment text (sentence
patterns; sentences about population databases, controls, unaffected people
or negations are skipped).  Up to three quoted snippets are kept so the
reviewer can check them.  A sentence that states the trans phase explicitly
("in trans" / "on the opposite allele"), not hedged ("presumed", "phase
unknown", "whether", "in cis"), with a partner named pathogenic or likely
pathogenic, marks the submitter as ``trans`` (schema 3; the partner condition
and ``pm3_independent`` from schema 4).  Submitters whose biallelic sentences
quote a common PMID describe the same published patients and are merged
(``pm3_independent`` groups, ``pm3_trans`` groups with an in-trans report).

PS3 (schema 3): submissions whose comment text reports a functional study
(experimental / in vitro / RNA / enzyme assays …) with an abnormal result for
this variant.  Sentences about in-silico predictions or modelling, other
variants (same codon, different substitution), negative / conflicting /
insufficient results or studies not performed are skipped; PMIDs quoted in
the sentence are kept.

Output ``clinvar_cases.sqlite`` (``cases``: variation_id, submitters,
individuals, submitter_ids, pm3_submitters, pm3_hom, pm3_comphet, pm3_detail,
pm3_trans, ps3_submitters, ps3_detail, pm3_independent; ``meta``).  The engine
caps PS4 from this table at Supporting, scores PM3 with the ClinGen SVI points
(explicit in trans 1, other biallelic reports 0.5 per group, at most Strong)
and caps PS3 at Supporting, because case details, independence and assay
validation cannot be checked.
"""
from __future__ import annotations

import gzip
import json
import os
import re
import shutil
import sqlite3
import subprocess
import time
import xml.etree.ElementTree as ET
from pathlib import Path

CASES_DB_NAME = "clinvar_cases.sqlite"
XML_NAME = "ClinVarVCVRelease_00-latest.xml.gz"
XML_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/" + XML_NAME
DEFAULT_METHODS = ("clinical testing", "research", "provider interpretation", "case-control")


def _local(tag) -> str:
    tag = tag if isinstance(tag, str) else ""
    return tag.rsplit("}", 1)[-1]


def _children(elem, name: str):
    return [c for c in elem if _local(c.tag) == name]


def _descendants(elem, name: str):
    return [c for c in elem.iter() if _local(c.tag) == name]


def _text(elem, name: str) -> str:
    for c in _descendants(elem, name):
        if c.text and c.text.strip():
            return c.text.strip()
    return ""


def _classification(assertion) -> str:
    """Submitter's germline classification (schema 2: Classification/GermlineClassification;
    schema 1: Interpretation/Description)."""
    for c in _children(assertion, "Classification"):
        value = _text(c, "GermlineClassification")
        if value:
            return value
    for c in _children(assertion, "Interpretation"):
        value = _text(c, "Description")
        if value:
            return value
    return ""


def _usable(assertion) -> bool:
    """Current submission that the submitter did not classify benign / likely benign."""
    status = ""
    for c in _children(assertion, "RecordStatus"):
        status = (c.text or "").strip().lower()
    if status and status != "current":
        return False
    cls = _classification(assertion).lower()
    return not ("benign" in cls and "pathogenic" not in cls)


# ── PM3: homozygous / compound heterozygous patients ──────────────────────
_BIALLELIC = re.compile(r"\b(homozyg\w*|compound[\s-]+heterozyg\w*|in[\s-]+trans\b|bi-?allelic|"
                        r"trans\s+configuration)", re.I)
_SKIP = re.compile(r"\b(gnomad|exac|1000\s*genomes|topmed|population|controls?|database|allele\s+frequency|"
                   r"unaffected|healthy|asymptomatic|carriers?\b|no\s+homozyg\w*|not\s+(?:been\s+)?"
                   r"(?:observed|reported|seen|found|detected|identified)|absent)", re.I)
_PERSON = re.compile(r"\b(affected|patients?|individuals?|probands?|cases?|famil(?:y|ies)|siblings?|"
                     r"child(?:ren)?|boys?|girls?|subjects?|infants?|members?|reported|identified|observed|"
                     r"detected|described)\b", re.I)
_SENTENCE = re.compile(r"(?<=[.;!?])\s+")
SNIPPET_CHARS = 240


def _texts(assertion) -> list[str]:
    out = []
    for c in assertion.iter():
        name = _local(c.tag)
        if name == "Comment" or (name == "Attribute" and (c.get("Type") or "") == "Description"):
            if c.text and c.text.strip():
                out.append(" ".join(c.text.split()))
    return out


def _snippet(sentence: str, start: int = 0, end: int = 0) -> str:
    """At most SNIPPET_CHARS of ``sentence``, centred on ``sentence[start:end]``."""
    sentence = sentence.strip()
    if len(sentence) <= SNIPPET_CHARS:
        return sentence
    middle = (start + end) // 2 if end > start else 0
    left = max(0, min(middle - SNIPPET_CHARS // 2, len(sentence) - (SNIPPET_CHARS - 2)))
    right = left + SNIPPET_CHARS - 2
    return ("…" if left else "") + sentence[left:right].strip() + ("…" if right < len(sentence) else "")


def biallelic_sentences(text: str) -> list[tuple[str, str]]:
    """[(kind, sentence)] with kind "hom" / "comphet" for sentences reporting patients in that state."""
    found = []
    for sentence in _SENTENCE.split(text or ""):
        hits = _BIALLELIC.findall(sentence)
        if not hits or _SKIP.search(sentence) or not _PERSON.search(sentence):
            continue
        kinds = {"hom" if h.lower().startswith("homozyg") else "comphet" for h in hits}
        for kind in sorted(kinds):
            found.append((kind, sentence.strip()))
    return found


# explicit trans phase (PM3 Strong, rule v1.13)
_TRANS = re.compile(r"\b(in[\s-]+trans\b|trans\s+configuration|"
                    r"(?:on\s+)?(?:the\s+)?(?:opposite|other|different|separate)\s+(?:parental\s+)?"
                    r"(?:allele|chromosome|haplotype)s?\b)", re.I)
_TRANS_HEDGE = re.compile(
    r"\b(?:presum\w*|possibl\w*|probabl\w*|likely|suspect\w*|assum\w*|thought|believed|inferred|predicted|"
    r"apparently)\s+(?:to\s+be\s+|as\s+)?(?:in[\s-]+trans|trans|on\s+(?:the\s+)?opposite)"
    r"|phase\s+(?:is\s+|was\s+|remains\s+)?(?:unknown|undetermined|unclear|not\s+(?:been\s+)?"
    r"(?:determined|known|established|confirmed|assessed))"
    r"|unknown\s+phase|phase[\s-]+unknown|undetermined\s+phase"
    r"|not\s+(?:been\s+)?(?:determined|confirmed|established|known|tested|available)\b"
    r"|\bwhether\b|\bif\s+(?:the|these|they)\b.{0,40}\bin[\s-]+trans|could\s+not|unable|\bin\s+cis\b"
    r"|parental\s+(?:testing|studies|samples|segregation)\s+(?:was\s+|were\s+)?(?:not|un)", re.I)
_UNCERTAIN_PARTNER = re.compile(r"uncertain\s+(?:clinical\s+)?significance|unknown\s+significance|"
                                r"\bVO?US\b|\bbenign\b|\bnot\s+(?:\w+\s+){0,2}pathogenic", re.I)
# the partner must be named pathogenic or likely pathogenic (rule v1.14)
_PLP_PARTNER = re.compile(r"(?<!non-)(?<!non )(?<!\w)pathogenic\b", re.I)


def explicit_trans(sentence: str, context: str = ""):
    """The match when ``sentence`` states the trans phase explicitly with a
    partner named pathogenic or likely pathogenic; None when hedged, in cis,
    with a VUS / benign partner or a partner whose class is not stated.
    ``context`` (the full stop-delimited sentence around a ``;`` clause) is also
    checked for hedges such as "parental testing was not performed; …"."""
    match = _TRANS.search(sentence or "")
    if not match:
        return None
    for text in (sentence, context):
        if text and (_TRANS_HEDGE.search(text) or _UNCERTAIN_PARTNER.search(text)):
            return None
    if not any(text and _PLP_PARTNER.search(text) for text in (sentence, context)):
        return None
    return match


def _full_sentence(text: str, clause: str) -> str:
    """The stop-delimited sentence of ``text`` that contains ``clause``."""
    at = text.find(clause)
    if at < 0:
        return clause
    start = max(text.rfind(c, 0, at) for c in ".!?") + 1
    stripped = clause.rstrip()
    if stripped.endswith((".", "!", "?")):              # the clause already ends the sentence
        return text[start:at + len(stripped)].strip()
    ends = [i for i in (text.find(c, at + len(clause)) for c in ".!?") if i >= 0]
    return text[start:(min(ends) + 1 if ends else len(text))].strip()


# ── PS3: functional studies (rule v1.13) ───────────────────────────────────
_FUNCTIONAL = re.compile(
    r"\b(functional(?:ly)?\s+(?:stud\w*|assay\w*|analys\w*|experiment\w*|characteri[sz]\w*|evidence|data|"
    r"evaluation|test\w*|expression)|experimental\s+(?:stud\w*|evidence|data)|in\s+vitro|in\s+vivo|"
    r"enzym\w*\s+(?:activit\w*|assay\w*)|minigene|rt-?pcr|(?:m?rna|cdna|transcript)\s+(?:stud\w*|analys\w*|"
    r"sequencing)|splicing\s+(?:assay\w*|stud\w*|analys\w*)|western\s+blot\w*|reporter\s+(?:gene\s+)?assay\w*|"
    r"electrophysiolog\w*|patch[\s-]clamp|(?:patient|cultured)[\s-]+(?:derived\s+)?(?:cells|fibroblasts|"
    r"lymphoblast\w*|lymphocytes|leukocytes)|cell[\s-]based|binding\s+assay\w*|knock-?in|"
    r"studies\s+(?:have\s+)?(?:shown|demonstrated)\s+that\s+this\s+(?:variant|alteration|(?:sequence\s+)?change)\s+"
    r"(?:results\s+in|causes|leads\s+to|disrupts|abolishes|activates|creates))\b", re.I)
_ABNORMAL = re.compile(
    r"\b(damaging|deleterious|deficien\w*|reduc\w*|decreas\w*|loss|lost|abolish\w*|abrogat\w*|impair\w*|"
    r"disrupt\w*|defect\w*|aberrant\w*|abnormal\w*|null|inactiv\w*|diminish\w*|lower|absent|undetectable|"
    r"dominant[\s-]negative|gain[\s-]of[\s-]function|affects?|affected|skip\w*|retention|cryptic|truncat\w*|"
    r"mislocali[sz]\w*|unstable|instability|degrad\w*|misfold\w*|altered|alters|residual|haploinsufficien\w*)\b",
    re.I)
_FUNCTIONAL_SKIP = re.compile(
    r"in\s+silico|predict\w*|computational|algorithm\w*|model(?:l)?ing|expected\s+to|\bsift\b|polyphen|"
    r"\brevel\b|\bcadd\b|spliceai|maxent|conservation|"
    r"different\s+(?:amino|variant|substitution|nucleotide|missense|change|alteration)|"
    r"another\s+(?:variant|missense|alteration|substitution|change|mutation)|other\s+(?:variants|missense|mutations)|"
    r"same\s+(?:codon|residue|position|amino)|this\s+(?:codon|residue|amino\s+acid\s+position)|nearby|neighbo|"
    r"adjacent|similar\s+(?:variant|change|alteration)|knock-?out|"
    r"no\s+(?:significant\s+|appreciable\s+|detectable\s+|substantial\s+|obvious\s+|apparent\s+|measurable\s+)?"
    r"(?:damaging\s+)?(?:effect|impact|difference|change|alteration|reduction|loss|defect)|"
    r"\b(?:does|did|do|was|were|is|has|have)\s+not\b|\bnot\s+(?:substantially\s+|significantly\s+)?"
    r"(?:affect|alter|disrupt|impair|reduce|impact|change)|\bnormal\b|comparable|similar\s+to\s+(?:the\s+)?wild|"
    r"wild[\s-]?type[\s-]like|indistinguishable|retain\w*|conflicting|insufficient|inconclusive|unclear|"
    r"uncertain|no\s+functional|absence\s+of\s+functional|not\s+(?:been\s+)?functionally|"
    r"have\s+not\s+been|has\s+not\s+been|were\s+not\s+performed|unknown", re.I)
_PMID = re.compile(r"\bPMIDs?[:\s]*((?:\d{5,9}[,;\s]*(?:and\s+)?)+)", re.I)


def pmids_in(text: str) -> list[str]:
    """PMIDs quoted as "PMID: 123, 456" in ``text``, in order, without repeats."""
    out: list[str] = []
    for group in _PMID.findall(text or ""):
        out += [p for p in re.findall(r"\d{5,9}", group) if p not in out]
    return out


def functional_sentences(text: str) -> list[tuple[str, re.Match, list[str]]]:
    """[(sentence, functional-term match, PMIDs)] reporting an abnormal functional result for this variant."""
    found = []
    for sentence in _SENTENCE.split(text or ""):
        match = _FUNCTIONAL.search(sentence)
        if not match or not _ABNORMAL.search(sentence) or _FUNCTIONAL_SKIP.search(sentence):
            continue
        found.append((sentence.strip(), match, pmids_in(sentence)))
    return found


def _assertion_functional(assertion) -> dict | None:
    """{"snippet", "pmids"} when the submission's text reports an abnormal functional study."""
    snippet, pmids = "", []
    for text in _texts(assertion):
        for sentence, match, found in functional_sentences(text):
            if not snippet:
                snippet = _snippet(sentence, match.start(), match.end())
            pmids += [p for p in found if p not in pmids]
    return {"snippet": snippet, "pmids": pmids[:10]} if snippet else None


def _structured_biallelic(assertion) -> dict[str, int]:
    """ObservedData counts such as ``Attribute Type="CompoundHeterozygote" integerValue="2"``."""
    out = {"hom": 0, "comphet": 0}
    for observed in _descendants(assertion, "ObservedIn"):
        samples = _children(observed, "Sample")
        if samples and _text(samples[0], "AffectedStatus").lower() == "no":
            continue
        for attr in _descendants(observed, "Attribute"):
            kind = (attr.get("Type") or "").lower().replace("_", "").replace(" ", "")
            raw = attr.get("integerValue") or (attr.text or "").strip()
            value = int(raw) if raw.isdigit() else 0
            if not value:
                continue
            if kind.startswith("compoundheterozygote"):
                out["comphet"] += value
            elif kind.startswith("homozygote"):
                out["hom"] += value
    return out


def _assertion_biallelic(assertion) -> dict | None:
    """{"hom", "comphet", "trans", "snippet", "pmids"} when the submission reports biallelic
    patients; the snippet quotes the explicit in-trans sentence when there is one and
    ``pmids`` are the PMIDs quoted in the biallelic sentences (used to merge submitters
    that describe the same published patients)."""
    counts = _structured_biallelic(assertion)
    snippet, trans_snippet, text_kinds, pmids = "", "", set(), []
    for text in _texts(assertion):
        for kind, sentence in biallelic_sentences(text):
            text_kinds.add(kind)
            context = _full_sentence(text, sentence)
            pmids += [p for p in pmids_in(context) if p not in pmids]
            trans = explicit_trans(sentence, context) if kind == "comphet" else None
            if trans and not trans_snippet:
                trans_snippet = _snippet(sentence, trans.start(), trans.end())
            if not snippet:
                hit = _BIALLELIC.search(sentence)
                snippet = _snippet(sentence, hit.start(), hit.end()) if hit else _snippet(sentence)
    hom = counts["hom"] or ("hom" in text_kinds)
    comphet = counts["comphet"] or ("comphet" in text_kinds)
    if not (hom or comphet):
        return None
    if not snippet:
        parts = [f"{counts['hom']} homozygote(s)" if counts["hom"] else "",
                 f"{counts['comphet']} compound heterozygote(s)" if counts["comphet"] else ""]
        snippet = "ClinVar observation counts: " + ", ".join(p for p in parts if p)
    return {"hom": int(bool(hom)), "comphet": int(bool(comphet)), "trans": int(bool(trans_snippet)),
            "snippet": trans_snippet or snippet, "pmids": pmids[:20]}


def independent_groups(biallelic: dict) -> tuple[int, int]:
    """(groups, groups with an explicit in-trans report): submitters whose biallelic
    sentences quote a common PMID describe the same published patients and count once;
    submitters without a PMID count on their own."""
    parent = {who: who for who in biallelic}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    by_pmid: dict[str, str] = {}
    for who, info in sorted(biallelic.items()):
        for pmid in info.get("pmids") or ():
            if pmid in by_pmid:
                parent[find(who)] = find(by_pmid[pmid])
            else:
                by_pmid[pmid] = who
    groups: dict[str, int] = {}
    for who, info in biallelic.items():
        root = find(who)
        groups[root] = max(groups.get(root, 0), int(info.get("trans") or 0))
    return len(groups), sum(groups.values())


def _assertion_cases(assertion, methods: set) -> int | None:
    """Affected individuals (≥1) reported by one submission, or None if it does not qualify."""
    best = None
    for observed in _descendants(assertion, "ObservedIn"):
        samples = _children(observed, "Sample")
        sample = samples[0] if samples else None
        if sample is None:
            continue
        if _text(sample, "AffectedStatus").lower() != "yes":
            continue
        origin = _text(sample, "Origin").lower()
        if "somatic" in origin:
            continue
        found = {(_text(m, "MethodType") or "").lower() for m in _descendants(observed, "Method")}
        if not found & methods:
            continue
        number = _text(sample, "NumberOfIndividuals")
        count = int(number) if number.isdigit() and int(number) > 0 else 1
        best = max(best or 0, count)
    return best


def iter_records(stream, methods=DEFAULT_METHODS):
    """Yield (variation_id, {submitter: affected individuals}, {submitter: biallelic evidence},
    {submitter: functional-study evidence})."""
    methods = {m.lower() for m in methods}
    variation = ""
    found: dict[str, int] = {}
    biallelic: dict[str, dict] = {}
    functional: dict[str, dict] = {}
    root = None
    for event, elem in ET.iterparse(stream, events=("start", "end")):
        name = _local(elem.tag)
        if event == "start":
            if root is None:
                root = elem
            if name == "VariationArchive":
                variation = elem.get("VariationID") or ""
                found, biallelic, functional = {}, {}, {}
            continue
        if name == "ClinicalAssertion":
            if _usable(elem):
                acc = _children(elem, "ClinVarAccession")
                who = (acc[0].get("OrgID") or acc[0].get("SubmitterName") or acc[0].get("Accession")) if acc else ""
                who = who or f"scv{len(found) + len(biallelic) + len(functional)}"
                cases = _assertion_cases(elem, methods)
                if cases:
                    found[who] = max(found.get(who, 0), cases)
                pm3 = _assertion_biallelic(elem)
                scv = acc[0].get("Accession", "") if acc else ""
                if pm3:
                    pm3["scv"] = scv
                    old = biallelic.get(who)
                    if old is not None:      # several submissions by one submitter: one entry
                        pmids = old["pmids"] + [p for p in pm3["pmids"] if p not in old["pmids"]]
                        if pm3["trans"] <= old["trans"]:
                            pm3 = old
                        pm3["pmids"] = pmids
                    biallelic[who] = pm3
                ps3 = _assertion_functional(elem)
                if ps3 and who not in functional:
                    ps3["scv"] = scv
                    functional[who] = ps3
            elem.clear()
        elif name == "VariationArchive":
            if variation and (found or biallelic or functional):
                yield variation, found, biallelic, functional
            elem.clear()
            if root is not None:
                root.clear()


def iter_cases(stream, methods=DEFAULT_METHODS):
    """Yield (variation_id, {submitter: individuals}) per VariationArchive with affected cases."""
    for variation, found, _, _ in iter_records(stream, methods):
        if found:
            yield variation, found


def open_xml(path: Path):
    """Decompressed stream; uses an external gzip when available (faster)."""
    path = Path(path)
    if str(path).endswith(".gz"):
        gz = shutil.which("pigz") or shutil.which("gzip")
        if gz:
            proc = subprocess.Popen([gz, "-dc", str(path)], stdout=subprocess.PIPE)
            return proc.stdout, proc
        return gzip.open(path, "rb"), None
    return open(path, "rb"), None


def build(xml_path: Path, out: Path, *, release: str = "", methods=DEFAULT_METHODS, progress=print) -> dict:
    """Parse ``xml_path`` and atomically write ``out``; returns stats."""
    started = time.time()
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(f".{out.name}.{os.getpid()}.tmp")
    tmp.unlink(missing_ok=True)
    conn = sqlite3.connect(tmp)
    stream, proc = open_xml(xml_path)
    variants = 0
    try:
        conn.executescript("""
            PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;
            CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE cases (variation_id TEXT PRIMARY KEY, submitters INTEGER NOT NULL,
                                individuals INTEGER NOT NULL, submitter_ids TEXT NOT NULL,
                                pm3_submitters INTEGER NOT NULL DEFAULT 0, pm3_hom INTEGER NOT NULL DEFAULT 0,
                                pm3_comphet INTEGER NOT NULL DEFAULT 0, pm3_detail TEXT NOT NULL DEFAULT '[]',
                                pm3_trans INTEGER NOT NULL DEFAULT 0, ps3_submitters INTEGER NOT NULL DEFAULT 0,
                                ps3_detail TEXT NOT NULL DEFAULT '[]', pm3_independent INTEGER NOT NULL DEFAULT 0)
                                WITHOUT ROWID;
        """)
        batch = []
        pm3_variants = ps3_variants = 0
        sql = "INSERT OR REPLACE INTO cases VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
        for variation, found, biallelic, functional in iter_records(stream, methods):
            ordered = sorted(biallelic.items(), key=lambda kv: (-kv[1]["trans"], kv[0]))
            detail = [{"scv": b["scv"], "snippet": b["snippet"], "trans": bool(b["trans"]),
                       "pmids": b["pmids"][:5]} for _, b in ordered][:3]
            groups, trans_groups = independent_groups(biallelic)
            ps3_detail = [{"scv": f["scv"], "snippet": f["snippet"], "pmids": f["pmids"]}
                          for _, f in sorted(functional.items())][:3]
            batch.append((variation, len(found), sum(found.values()), ",".join(sorted(found)),
                          len(biallelic), sum(b["hom"] for b in biallelic.values()),
                          sum(b["comphet"] for b in biallelic.values()),
                          json.dumps(detail, ensure_ascii=False), trans_groups,
                          len(functional), json.dumps(ps3_detail, ensure_ascii=False), groups))
            variants += 1
            pm3_variants += bool(biallelic)
            ps3_variants += bool(functional)
            if len(batch) >= 5000:
                conn.executemany(sql, batch)
                batch.clear()
                progress(f"  ClinVar 個案資料：{variants:,} 個變異")
        conn.executemany(sql, batch)
        meta = {"schema_version": "4", "release": release, "source": str(xml_path),
                "methods": ",".join(methods), "variants": str(variants), "pm3_variants": str(pm3_variants),
                "ps3_variants": str(ps3_variants),
                "built_at": time.strftime("%Y-%m-%dT%H:%M:%S")}
        conn.executemany("INSERT INTO meta VALUES (?,?)", sorted(meta.items()))
        conn.commit()
    except BaseException:
        conn.close()
        tmp.unlink(missing_ok=True)
        raise
    finally:
        if proc is not None:
            stream.close()
            proc.wait()
    conn.close()
    if proc is not None and proc.returncode not in (0, None):
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"解壓縮 {xml_path} 失敗（exit {proc.returncode}）")
    os.replace(tmp, out)
    return {"variants": variants, "pm3_variants": pm3_variants, "ps3_variants": ps3_variants,
            "seconds": round(time.time() - started, 1), "path": str(out)}
