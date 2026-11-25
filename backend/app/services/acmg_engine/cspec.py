"""Turn ClinGen CSpec (VCEP criteria specifications) and UniProt collagen
sequences into the small tables the rule engine reads.

Output of ``build_vcep_specs`` (stored as ``vcep_specs.json``)::

    {"genes": {GENE: {"svi", "version", "vcep", "label", "status", "url",
                      "criteria": {CODE: [allowed strengths]},
                      "thresholds": {"ba1": float, "bs1": float, "pm2": float},
                      "threshold_notes": {...}, "pm1": [{"strength", "ranges"}],
                      "pm1_unparsed": bool}},
     "skipped": {GENE: reason}}

CSpec rules are free text.  Only what can be read reliably is used:

* per criterion, the strengths marked "Applicable" (a criterion with none is
  not applicable for that gene; a stronger strength than allowed is capped);
* allele-frequency thresholds for BA1 / BS1 / PM2: numbers written next to a
  comparator (≥, <, above, below, cutoff …).  When a rule lists several
  values, BA1 takes the largest (fewest benign calls), BS1 the smallest (the
  lower bound of an "≥X and <Y" range, the upper bound being BA1) and PM2 the
  smallest (strictest).  "Absent from controls" with no number means 0;
* PM1 hotspots: protein residue ranges and single residues (``aa 19-43``,
  ``codons 167-931``, ``Tyr138``, ``R391-S392`` …).  Tables only attached as
  files, exon-level and cDNA-level definitions cannot be read; such genes are
  marked ``pm1_unparsed`` and keep the generic ClinVar hotspot rule.
  Genomic coordinates (``17:8002593-8002598``), intra-domain cysteine
  numbering (``Cys2-Cys3``) and disulfide pairs (``C350-C382 disulfide``)
  are not residue ranges.  A few rules the text parser gets wrong are fixed
  in ``PM1_OVERRIDES`` (checked by hand against the CSpec text, applied only
  to the spec version they were written for).

Genes with more than one rule set (different diseases, e.g. RYR1, VWF) are
skipped: the right one cannot be chosen per variant.
"""
from __future__ import annotations

import re

STRENGTHS = {"stand alone": "stand_alone", "very strong": "very_strong", "strong": "strong",
             "moderate": "moderate", "supporting": "supporting"}
STRENGTH_RANK = {"supporting": 1, "moderate": 2, "strong": 3, "very_strong": 4, "stand_alone": 5}
USABLE_STATUS = {"Released", "Approved For Release"}

_AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
        "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
        "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
_AA = r"(?:" + "|".join(_AA3) + r"|[ACDEFGHIKLMNPQRSTVWY])"

# ── text helpers ───────────────────────────────────────────────────────────
def clean(text: str | None) -> str:
    """Drop markdown/HTML reference markers and escapes that break number parsing."""
    t = (text or "").replace("\xa0", " ").replace("\\", "")
    t = re.sub(r"\[<sup>.*?</sup>\]\([^)]*\)", " ", t)          # [<sup>3</sup>](#pmid_…)
    t = re.sub(r"<sup>\s*(?:\d+\s*,?\s*)+</sup>(?!\s*[-−–])", " ", t)  # citation superscripts
    t = re.sub(r"\(#?pmid_?\d+\)|#?pmid_\d+", " ", t, flags=re.I)
    t = re.sub(r"\*\*|__", "", t)
    return t


def short(text: str | None, limit: int = 400) -> str:
    return re.sub(r"\s+", " ", clean(text)).strip()[:limit]


_CMP = (r"(?:≥|>=|=>|>|≤|<=|=<|<|above|below|greater than(?: or equal to)?|greater or equal to|"
        r"less than(?: or equal to)?|exceed(?:s|ing)?|cut-?off(?: of)?|at least|at or above|"
        r"over|under|threshold(?: of)?|frequency of|is|from|between|and|to)")
_NUM = (r"(\d+(?:\.\d+)?)\s*(?:x|×|\*|X)\s*10\s*(?:<sup>|\^)?\s*\(?([-−–]\s*\d+)\)?\s*(?:</sup>)?"
        r"|(\d*\.\d+|\d+)\s*(%)?")
_CMP_NUM = re.compile(_CMP + r"\s*[\s(_*]*(?:" + _NUM + ")", re.I)
_ANY_NUM = re.compile(_NUM, re.I)


def _to_freq(m) -> float | None:
    if m.group(1):
        exp = int(re.sub(r"[−–\s]", "-", m.group(2)).replace("--", "-"))
        value = float(m.group(1)) * 10 ** exp
    else:
        number, pct = m.group(3), m.group(4)
        if pct:
            value = float(number) / 100
        elif "." in number:
            value = float(number)
        else:
            return None                         # integers are allele counts, not frequencies
    return round(value, 12) if 0 < value < 0.3 else None


def frequencies(text: str | None) -> list[float]:
    """Frequencies stated next to a comparator; all frequency-like numbers as fallback."""
    t = clean(text)
    values = {v for m in _CMP_NUM.finditer(t) if (v := _to_freq(m)) is not None}
    if not values:
        values = {v for m in _ANY_NUM.finditer(t) if (v := _to_freq(m)) is not None}
    return sorted(values)


# ── PM1 residues ───────────────────────────────────────────────────────────
_STRIP_FOR_RESIDUES = re.compile(
    r"N[MP]_\d+(?:\.\d+)?|ENS[TP]\d+(?:\.\d+)?|PMIDs?:?\s*[\d,\s]+|\b[cn]\.-?\d+(?:[+-]\d+)?"
    r"|exons?\s+\d+(?:\s*(?:,|and|or|-)\s*\d+)*|\bv\d+(?:\.\d+)*|gnomAD\s*v?\d[\d.]*"
    r"|\b(?:chr)?(?:\d{1,2}|X|Y):\d{4,}(?:\s*[-–]\s*\d+)?(?:\s*\((?:hg|GRCh)\d+\))?", re.I)
_PAIR_WORDS = re.compile(r"^\W{0,3}(?:disulfide|disulphide|bond)", re.I)
_DOMAIN_CYS_MAX = 8          # EGF-like domains number their cysteines Cys1–Cys8
_AA_RANGE = re.compile(r"(?:p\.)?" + _AA + r"-?(\d{1,4})\s*(?:-|–|—|to)\s*(?:p\.)?" + _AA + r"-?(\d{1,4})")
_ITEM = r"(?:p\.)?\d{1,4}(?:\s*(?:-|–|—|to)\s*\d{1,4})?"
_KEY_LIST = re.compile(r"(?:codons?|residues?|amino acids?|positions?|motifs?|\baa(?=[\s\d\[:]))"
                       r"(?:[^\d\n.;]{0,40}?:)?\s*[:\[]?\s*"
                       r"(" + _ITEM + r"(?:\s*(?:,|;|and|or)\s*" + _ITEM + r")*)", re.I)
_AA_SINGLE = re.compile(r"(?<![A-Za-z0-9])(?:p\.)?(" + "|".join(_AA3) + r")-?(\d{1,4})(?![\d])"
                        r"|(?<![A-Za-z0-9.])(?:p\.)?([ACDEFGHIKLMNPQRSTVWY])(\d{2,4})(?:[ACDEFGHIKLMNPQRSTVWY]|\*)?(?![\dA-Za-z])")


def residues(text: str | None) -> list[list[int]]:
    """Protein positions named in a PM1 rule, merged into sorted [start, end] ranges."""
    t = _STRIP_FOR_RESIDUES.sub(" ", clean(text))
    spans: list[tuple[int, int]] = []
    for m in _AA_RANGE.finditer(t):
        a, b = int(m.group(1)), int(m.group(2))
        if max(a, b) <= _DOMAIN_CYS_MAX:          # "between Cys2-Cys3 in cbEGF-like domain"
            continue
        if _PAIR_WORDS.match(t[m.end():]):         # "C350-C382 disulfide": two residues, not a range
            spans.extend([(a, a), (b, b)])
            continue
        spans.append((min(a, b), max(a, b)))
    t2 = _AA_RANGE.sub(" ", t)
    for m in _KEY_LIST.finditer(t2):
        for part in re.split(r"\s*(?:,|;|and|or)\s*", m.group(1)):
            nums = [int(x) for x in re.findall(r"\d+", part)]
            if len(nums) == 1:
                spans.append((nums[0], nums[0]))
            elif len(nums) == 2:
                spans.append((min(nums), max(nums)))
    t3 = _KEY_LIST.sub(" ", t2)
    for m in _AA_SINGLE.finditer(t3):
        pos = int(m.group(2) or m.group(4))
        spans.append((pos, pos))
    spans = sorted({s for s in spans if 0 < s[0] <= s[1] < 40000})
    merged: list[list[int]] = []
    for a, b in spans:
        if merged and a <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return merged


# ── CSpec documents ────────────────────────────────────────────────────────
def _strength_key(label: str | None) -> str | None:
    return STRENGTHS.get(str(label or "").strip().lower())


def _applicable(item: dict) -> bool:
    """"Applicable", "Applicable with VCEP specification", "Applicable as originally described"."""
    return str(item.get("applicability") or "").strip().lower().startswith("applicable")


def rule_sets(doc: dict) -> list[dict]:
    """One entry per rule set of a CSpec SequenceVariantInterpretation document."""
    svi_id = str(doc.get("@id") or "").rsplit("/", 1)[-1]
    status = doc.get("cspecStatus") or doc.get("currentStatus") or doc.get("status") or ""
    vcep = ((doc.get("affiliation") or {}).get("label") or "").strip()
    out = []
    for rs in doc.get("ruleSets") or []:
        genes = [str(g.get("label") or "").strip() for g in rs.get("genes") or [] if g.get("label")]
        criteria: dict[str, dict] = {}
        for cc in rs.get("criteriaCodes") or []:
            code = str(cc.get("label") or "").strip().upper()
            if not code:
                continue
            allowed, texts = [], {}
            for item in cc.get("evidenceStrengths") or []:
                key = _strength_key(item.get("label"))
                if key and _applicable(item):
                    allowed.append(key)
                    texts[key] = item.get("description") or cc.get("description") or ""
            criteria[code] = {"allowed": sorted(set(allowed), key=STRENGTH_RANK.get), "texts": texts}
        out.append({"svi": svi_id, "version": str(doc.get("version") or ""), "status": status,
                    "vcep": vcep, "label": doc.get("label") or "",
                    "url": f"https://cspec.genome.network/cspec/ui/svi/doc/{svi_id}" if svi_id else "",
                    "genes": genes, "criteria": criteria})
    return out


_MOI_WORDS = {"ar": re.compile(r"recessive", re.I), "ad": re.compile(r"dominant", re.I)}


def _moi_split(text: str) -> dict[str, list[float]]:
    """{"ar": [...], "ad": [...]} when a rule states separate values per inheritance."""
    out: dict[str, list[float]] = {}
    for clause in re.split(r"[;,\n]|\.\s", clean(text)):
        kinds = [k for k, rx in _MOI_WORDS.items() if rx.search(clause)]
        if len(kinds) != 1:
            continue
        found = frequencies(clause)
        if found:
            out.setdefault(kinds[0], []).extend(found)
    return out if len(out) == 2 else {}


def thresholds(criteria: dict) -> tuple[dict, dict]:
    """(``{"ba1", "bs1", "pm2"}`` values, notes with the rule text used).

    When a rule gives separate values for recessive and dominant inheritance,
    ``ba1_ar`` / ``ba1_ad`` (etc.) are added; the plain key keeps the
    conservative pick across both for genes whose inheritance is unknown.
    """
    values, notes = {}, {}

    def texts(code, prefer):
        crit = criteria.get(code) or {}
        if not crit.get("allowed"):
            return None
        for key in prefer:
            if key in crit["texts"]:
                return crit["texts"][key]
        return next(iter(crit["texts"].values()), None)

    for name, code, prefer, pick in (("ba1", "BA1", ("stand_alone",), max),
                                     ("bs1", "BS1", ("strong",), min),
                                     ("pm2", "PM2", ("supporting", "moderate"), min)):
        text = texts(code, prefer)
        if text is None:
            continue
        found = frequencies(text)
        if found:
            values[name] = pick(found)
            how = "single" if len(found) == 1 else f"{pick.__name__} of {found}"
            for moi, vals in _moi_split(text).items():
                values[f"{name}_{moi}"] = pick(vals)
                how += f"; {moi.upper()} {pick(vals)}"
            notes[name] = f"{how}: {short(text, 300)}"
        elif name == "pm2" and re.search(r"\babsent\b", clean(text), re.I):
            values[name] = 0.0
            notes[name] = f"absent: {short(text, 300)}"
    return values, notes


def _gene_segments(text: str, genes: list[str]) -> dict[str, str] | None:
    """Split a rule that names several of the rule set's genes into per-gene text."""
    hits = sorted((m.start(), g) for g in genes for m in re.finditer(r"\b" + re.escape(g) + r"\b", text))
    if not hits:
        return None
    out: dict[str, str] = {}
    for i, (pos, gene) in enumerate(hits):
        stop = hits[i + 1][0] if i + 1 < len(hits) else len(text)
        out[gene] = out.get(gene, "") + " " + text[pos:stop]
    return out


def pm1_rules(criteria: dict, gene: str | None = None, genes: list[str] | None = None) -> tuple[list[dict], str]:
    """([{"strength", "ranges"}], state) for PM1.

    ``state``: "" (parsed or PM1 not applicable), "unparsed" (rule text could
    not be read) or "other_gene" (a multi-gene rule set defines PM1 only for
    other genes).
    """
    crit = criteria.get("PM1") or {}
    rules, unparsed, other = [], False, False
    for strength in crit.get("allowed") or []:
        text = crit["texts"].get(strength) or ""
        if gene and genes and len(genes) > 1:
            segments = _gene_segments(clean(text), genes)
            if segments is not None:
                if gene not in segments:
                    other = True
                    continue
                text = segments[gene]
        ranges = residues(text)
        if ranges:
            rules.append({"strength": strength, "ranges": ranges})
        else:
            unparsed = True
    if rules:
        return rules, ""
    return [], "unparsed" if unparsed else ("other_gene" if other else "")


# Hand-checked corrections to PM1 rules the text parser reads wrongly.  Each
# entry applies only to the CSpec version it was written for; a newer version
# keeps the parsed rule and is listed in ``pm1_override_stale`` so it gets
# re-checked.  ``replace``: strength -> rule (ranges, or residue_change "Cys" =
# the reference or the new amino acid is cysteine); ``add``: extra ranges;
# ``not_applicable``: the VCEP defines no PM1 hotspot for the gene.
PM1_OVERRIDES = {
    "FBN1": {"version": "1.0", "replace": {"moderate": {"residue_change": "Cys"}},
             "note": "規格以胺基酸種類定義（EGF／TB／hybrid 區域的 Cys、產生 Cys 的變異；cbEGF 的 Cys 為 Strong），"
                     "平台只用「失去或產生 Cys」給 Moderate；鈣結合序列、羥化位置與 Cys 之間的關鍵 Gly 需人工判斷"},
    "PAH": {"version": "2.1", "add": {"moderate": [[246, 266], [280, 283], [322, 326], [377, 379]]},
            "note": "補上輔因子結合位置 246–266、280–283、322–326、377–379（文字列在 His285 等之後，解析器漏讀）"},
    "HNF4A": {"version": "4.0", "replace": {"supporting": {"ranges": [[37, 113], [180, 220], [300, 350]]}},
              "note": "DNA 結合區從 Leu37 起算（規格註明 33–36 不結合 DNA）"},
    "TNNT2": {"version": "1.0", "replace": {"supporting": {"ranges": [[89, 189]]}},
              "note": "用 MANE 轉錄本 NM_001276345.2 的 codons 89–189（ENST00000367318 的 79–179 是同一段）"},
    "GP1BA": {"version": "1.1", "replace": {"moderate": {"ranges": [[20, 20], [33, 33], [225, 225], [227, 227],
                                                                    [264, 264], [280, 280], [526, 527]]}},
              "note": "規格列出的 Cys 位置 20、33、225、227、264、280、526、527（只寫數字，解析器讀不到）"},
    "GP1BB": {"version": "1.1", "replace": {"moderate": {"ranges": [[93, 93], [95, 95], [118, 118], [141, 141],
                                                                    [147, 147]]}},
              "note": "規格列出的 Cys 位置 93、95、118、141、147（只寫數字，解析器讀不到）"},
    "PAX6": {"version": "1.0", "replace": {"moderate": {"ranges": [[4, 143], [226, 281]]}},
             "note": "規格以 cDNA 位置寫成 paired domain c.10–c.427、homeodomain c.676–c.841，換算為 codons 4–143、226–281"},
    "MYO15A": {"version": "1.0", "not_applicable": True,
               "note": "規格寫明 MYO15A 沒有定義 PM1 熱點，不改用 ClinVar 推算"},
    "OTOF": {"version": "1.0", "not_applicable": True,
             "note": "規格寫明 OTOF 沒有定義 PM1 熱點，不改用 ClinVar 推算"},
}


def _apply_pm1_override(gene: str, version: str, rules: list[dict]) -> tuple[list[dict], str | None, bool]:
    """(rules, note, stale) after the hand-checked override for ``gene``."""
    fix = PM1_OVERRIDES.get(gene)
    if not fix:
        return rules, None, False
    if str(version) != fix["version"]:
        return rules, None, True
    if fix.get("not_applicable"):
        return [], fix["note"], False
    by = {r["strength"]: dict(r) for r in rules}
    for strength, rule in (fix.get("replace") or {}).items():
        by[strength] = {"strength": strength, **rule}
    for strength, extra in (fix.get("add") or {}).items():
        merged = sorted({tuple(x) for x in (by.get(strength, {}).get("ranges") or [])} | {tuple(x) for x in extra})
        out: list[list[int]] = []
        for a, b in merged:
            if out and a <= out[-1][1] + 1:
                out[-1][1] = max(out[-1][1], b)
            else:
                out.append([a, b])
        by[strength] = {**by.get(strength, {}), "strength": strength, "ranges": out}
    order = {k: i for i, k in enumerate(STRENGTH_RANK)}
    return sorted(by.values(), key=lambda r: order.get(r["strength"], 9)), fix["note"], False


def build_vcep_specs(documents) -> dict:
    by_gene: dict[str, list[dict]] = {}
    for doc in documents:
        status = doc.get("cspecStatus") or doc.get("currentStatus") or doc.get("status")
        if status not in USABLE_STATUS:
            continue
        for rs in rule_sets(doc):
            for gene in rs["genes"]:
                by_gene.setdefault(gene, []).append(rs)
    genes, skipped, stale = {}, {}, []
    for gene, sets in sorted(by_gene.items()):
        if len(sets) > 1:
            skipped[gene] = "多個規格（" + "、".join(f"{s['svi']} {s['vcep']}" for s in sets) + "）"
            continue
        rs = sets[0]
        values, notes = thresholds(rs["criteria"])
        pm1, state = pm1_rules(rs["criteria"], gene, rs["genes"])
        pm1, override, is_stale = _apply_pm1_override(gene, rs["version"], pm1)
        if is_stale:
            stale.append(gene)
        if override:
            state = ""
        genes[gene] = {
            "svi": rs["svi"], "version": rs["version"], "status": rs["status"], "vcep": rs["vcep"],
            "label": rs["label"], "url": rs["url"], "rule_set_genes": rs["genes"],
            "criteria": {code: c["allowed"] for code, c in sorted(rs["criteria"].items())},
            "thresholds": values, "threshold_notes": notes,
            "pm1": pm1, "pm1_unparsed": state == "unparsed", "pm1_other_gene": state == "other_gene",
            "pm1_text": {k: short(v, 400) for k, v in (rs["criteria"].get("PM1") or {}).get("texts", {}).items()},
        }
        if override:
            genes[gene]["pm1_override"] = override
            if PM1_OVERRIDES[gene].get("not_applicable"):
                genes[gene]["criteria"]["PM1"] = []
    return {"genes": genes, "skipped": skipped, "pm1_override_stale": stale}


# ── collagen Gly-X-Y ──────────────────────────────────────────────────────
MIN_GXY_RUN = 6     # consecutive Gly-X-Y triplets counted as triple helix


def gly_xy_positions(sequence: str, helix: list | None = None, chain: tuple | None = None) -> list[int]:
    """1-based positions of the glycines in Gly-X-Y runs (≥ MIN_GXY_RUN triplets).

    Restricted to UniProt's annotated triple-helical regions when given, else
    to the mature chain (propeptides excluded) when given.
    """
    hits: set[int] = set()
    n = len(sequence)
    for start in range(3):
        run: list[int] = []
        for i in range(start, n, 3):
            if sequence[i] == "G":
                run.append(i + 1)
            else:
                if len(run) >= MIN_GXY_RUN:
                    hits.update(run)
                run = []
        if len(run) >= MIN_GXY_RUN:
            hits.update(run)
    if helix:
        hits = {p for p in hits if any(a <= p <= b for a, b in helix)}
    elif chain:
        hits = {p for p in hits if chain[0] <= p <= chain[1]}
    return sorted(hits)


def uniprot_collagen_entry(record: dict) -> dict:
    """UniProt REST JSON record -> {"accession", "length", "source", "positions"}."""
    seq = record["sequence"]["value"]
    feats = record.get("features") or []
    def span(f):
        start, end = (f.get("location") or {}).get("start", {}).get("value"), (f.get("location") or {}).get("end", {}).get("value")
        return (start, end) if isinstance(start, int) and isinstance(end, int) else None

    helix = [s for f in feats if f.get("type") == "Region"
             and "triple-helical" in (f.get("description") or "").lower() and (s := span(f))]
    chains = [s for f in feats if f.get("type") == "Chain" and (s := span(f))]
    chain = (min(a for a, _ in chains), max(b for _, b in chains)) if chains else None
    return {"accession": record.get("primaryAccession"), "length": len(seq),
            "source": "UniProt triple-helical region" if helix else "Gly-X-Y run within mature chain",
            "helix": [list(h) for h in helix], "positions": gly_xy_positions(seq, helix, chain)}
