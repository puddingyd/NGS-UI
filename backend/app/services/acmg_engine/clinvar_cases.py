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
reviewer can check them.

Output ``clinvar_cases.sqlite`` (``cases``: variation_id, submitters,
individuals, submitter_ids, pm3_submitters, pm3_hom, pm3_comphet, pm3_detail;
``meta``).  The engine caps PS4 from this table at Supporting and PM3 at
Moderate because case details, phase and independence cannot be checked.
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
    """{"hom", "comphet", "snippet"} when the submission reports biallelic patients."""
    counts = _structured_biallelic(assertion)
    snippet, text_kinds = "", set()
    for text in _texts(assertion):
        for kind, sentence in biallelic_sentences(text):
            text_kinds.add(kind)
            if not snippet:
                snippet = sentence if len(sentence) <= SNIPPET_CHARS else sentence[:SNIPPET_CHARS - 1] + "…"
    hom = counts["hom"] or ("hom" in text_kinds)
    comphet = counts["comphet"] or ("comphet" in text_kinds)
    if not (hom or comphet):
        return None
    if not snippet:
        parts = [f"{counts['hom']} homozygote(s)" if counts["hom"] else "",
                 f"{counts['comphet']} compound heterozygote(s)" if counts["comphet"] else ""]
        snippet = "ClinVar observation counts: " + ", ".join(p for p in parts if p)
    return {"hom": int(bool(hom)), "comphet": int(bool(comphet)), "snippet": snippet}


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
    """Yield (variation_id, {submitter: affected individuals}, {submitter: biallelic evidence})."""
    methods = {m.lower() for m in methods}
    variation = ""
    found: dict[str, int] = {}
    biallelic: dict[str, dict] = {}
    root = None
    for event, elem in ET.iterparse(stream, events=("start", "end")):
        name = _local(elem.tag)
        if event == "start":
            if root is None:
                root = elem
            if name == "VariationArchive":
                variation = elem.get("VariationID") or ""
                found, biallelic = {}, {}
            continue
        if name == "ClinicalAssertion":
            if _usable(elem):
                acc = _children(elem, "ClinVarAccession")
                who = (acc[0].get("OrgID") or acc[0].get("SubmitterName") or acc[0].get("Accession")) if acc else ""
                who = who or f"scv{len(found) + len(biallelic)}"
                cases = _assertion_cases(elem, methods)
                if cases:
                    found[who] = max(found.get(who, 0), cases)
                pm3 = _assertion_biallelic(elem)
                if pm3 and who not in biallelic:
                    pm3["scv"] = acc[0].get("Accession", "") if acc else ""
                    biallelic[who] = pm3
            elem.clear()
        elif name == "VariationArchive":
            if variation and (found or biallelic):
                yield variation, found, biallelic
            elem.clear()
            if root is not None:
                root.clear()


def iter_cases(stream, methods=DEFAULT_METHODS):
    """Yield (variation_id, {submitter: individuals}) per VariationArchive with affected cases."""
    for variation, found, _ in iter_records(stream, methods):
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
                                pm3_comphet INTEGER NOT NULL DEFAULT 0, pm3_detail TEXT NOT NULL DEFAULT '[]')
                                WITHOUT ROWID;
        """)
        batch = []
        pm3_variants = 0
        sql = "INSERT OR REPLACE INTO cases VALUES (?,?,?,?,?,?,?,?)"
        for variation, found, biallelic in iter_records(stream, methods):
            detail = [{"scv": b["scv"], "snippet": b["snippet"]} for _, b in sorted(biallelic.items())][:3]
            batch.append((variation, len(found), sum(found.values()), ",".join(sorted(found)),
                          len(biallelic), sum(b["hom"] for b in biallelic.values()),
                          sum(b["comphet"] for b in biallelic.values()),
                          json.dumps(detail, ensure_ascii=False)))
            variants += 1
            pm3_variants += bool(biallelic)
            if len(batch) >= 5000:
                conn.executemany(sql, batch)
                batch.clear()
                progress(f"  ClinVar 個案資料：{variants:,} 個變異")
        conn.executemany(sql, batch)
        meta = {"schema_version": "2", "release": release, "source": str(xml_path),
                "methods": ",".join(methods), "variants": str(variants), "pm3_variants": str(pm3_variants),
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
    return {"variants": variants, "pm3_variants": pm3_variants, "seconds": round(time.time() - started, 1),
            "path": str(out)}
