"""Field parsers for the ACMG rule engine.

Every helper is pure and tolerant: TSV cells may be empty, ``.``, ``NA`` or
multi-valued (``&`` / ``,``).  Callers decide what "missing" means; the
helpers only return ``None`` (or an empty collection) for it.
"""
from __future__ import annotations

import re
from urllib.parse import unquote

MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "-"}

IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
TRANSCRIPT_RANK = {"MANE_SELECT": 0, "MANE_PLUS_CLINICAL": 1, "CANONICAL": 2, "APPRIS_P1": 3}
MANE_TYPES = {"MANE_SELECT", "MANE_PLUS_CLINICAL"}

NONSENSE = {"stop_gained"}
FRAMESHIFT = {"frameshift_variant"}
CANONICAL_SPLICE = {"splice_donor_variant", "splice_acceptor_variant"}
START_LOST = {"start_lost"}
LOF_CONSEQUENCES = NONSENSE | FRAMESHIFT | CANONICAL_SPLICE | START_LOST
MISSENSE = {"missense_variant"}
INFRAME = {"inframe_insertion", "inframe_deletion"}
STOP_LOST = {"stop_lost"}
SYNONYMOUS = {"synonymous_variant"}
INTRONIC = {"intron_variant", "splice_region_variant", "splice_donor_5th_base_variant",
            "splice_donor_region_variant", "splice_polypyrimidine_tract_variant"}

_AA3 = {
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys",
    "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Ter", "Sec", "Pyl",
}
_POS_RE = re.compile(r"^([-*]?)(\d+)([+-]\d+)?$")
_HGVS_C_RE = re.compile(r"(?:^|:)c\.(?P<start>[-*]?\d+(?:[+-]\d+)?)(?:_(?P<end>[-*]?\d+(?:[+-]\d+)?))?")
_HGVS_P_RE = re.compile(r"p\.\(?(?P<ref>[A-Z][a-z]{2})(?P<pos>\d+)(?P<alt>[A-Z][a-z]{2}|=)")


def is_missing(value) -> bool:
    return value is None or str(value).strip().upper() in MISSING


def text(value) -> str:
    return "" if is_missing(value) else str(value).strip()


def numbers(value) -> list[float]:
    """All numeric parts of a cell (``&``/``,`` separated); missing parts skipped."""
    if is_missing(value):
        return []
    out: list[float] = []
    for part in re.split(r"[&,|]", str(value)):
        part = part.strip()
        if part.upper() in MISSING:
            continue
        try:
            out.append(float(part))
        except ValueError:
            continue
    return out


def first_number(value) -> float | None:
    vals = numbers(value)
    return vals[0] if vals else None


def max_abs(value) -> float | None:
    vals = numbers(value)
    return max((abs(v) for v in vals), default=None)


def consequences(value) -> set[str]:
    if is_missing(value):
        return set()
    return {part.strip().lower() for part in re.split(r"[&,]", str(value)) if part.strip()}


def variant_key(row: dict) -> str:
    chrom = text(row.get("CHROM"))
    if chrom and not chrom.lower().startswith("chr"):
        chrom = "chr" + chrom
    return f"{chrom}-{text(row.get('POS'))}-{text(row.get('REF'))}-{text(row.get('ALT'))}"


def is_snv(row: dict) -> bool:
    return len(text(row.get("REF"))) == 1 and len(text(row.get("ALT"))) == 1


def fraction(value) -> tuple[int | None, int | None]:
    """Parse ``a/b`` or ``a-c/b`` (protein position, exon) -> (a, b)."""
    raw = text(value)
    if not raw or "/" not in raw:
        return None, None
    left, _, right = raw.partition("/")
    start = re.match(r"\??(\d+)", left.strip())
    total = re.match(r"(\d+)", right.strip())
    return (int(start.group(1)) if start else None, int(total.group(1)) if total else None)


def _parse_pos(token: str):
    m = _POS_RE.match(token)
    if not m:
        return None
    prefix, base, offset = m.groups()
    return prefix, int(base), int(offset) if offset else 0


def hgvs_c_offsets(hgvs_c) -> list[int] | None:
    """Intronic offsets covered by an HGVS c. change (3'-normalised form).

    ``c.1097+6_1097+11del`` -> [6..11]; ``c.456+1G>C`` -> [1];
    ``c.1096_1097+3del`` -> [1, 2, 3] (exon end through +3).  Exonic-only
    changes return ``[0]``.  ``None`` when the notation cannot be parsed.
    """
    raw = text(hgvs_c)
    m = _HGVS_C_RE.search(raw)
    if not m:
        return None
    start = _parse_pos(m.group("start"))
    end = _parse_pos(m.group("end")) if m.group("end") else start
    if start is None or end is None:
        return None
    (_, b1, o1), (_, b2, o2) = start, end
    if b1 == b2 and o1 * o2 >= 0 and (o1 or o2):
        lo, hi = min(o1, o2), max(o1, o2)
        if lo == 0:
            lo = 1 if hi > 0 else lo
        if hi == 0:
            hi = -1 if lo < 0 else hi
        return list(range(lo, hi + 1))
    covered: list[int] = []
    if o1 < 0:
        covered.extend(range(o1, 0))
    if o2 > 0:
        covered.extend(range(1, o2 + 1))
    return covered or [0]


def canonical_splice_hit(hgvs_c) -> bool | None:
    offsets = hgvs_c_offsets(hgvs_c)
    if offsets is None:
        return None
    return any(o in (1, 2, -1, -2) for o in offsets)


def deep_intronic(hgvs_c, donor_min: int, acceptor_max: int) -> bool:
    """True when every covered offset lies beyond +donor_min / -acceptor_max."""
    offsets = hgvs_c_offsets(hgvs_c)
    if not offsets or 0 in offsets:
        return False
    if all(o > 0 for o in offsets):
        return min(offsets) >= donor_min
    if all(o < 0 for o in offsets):
        return max(offsets) <= acceptor_max
    return False


def hgvs_p(hgvs) -> tuple[str, int, str] | None:
    """``ENSP…:p.Arg123Trp`` -> ("Arg", 123, "Trp"); synonymous alt is ``=``."""
    raw = unquote(text(hgvs))
    m = _HGVS_P_RE.search(raw)
    if not m or m.group("ref") not in _AA3:
        return None
    alt = m.group("alt")
    if alt != "=" and alt not in _AA3:
        return None
    return m.group("ref"), int(m.group("pos")), alt


def giab_labels(value) -> set[str]:
    return {part.strip().lower() for part in re.split(r"[,&;|]", text(value)) if part.strip()}


_CLINVAR_EXCLUDE = ("low penetrance", "risk factor", "drug response", "association",
                    "protective", "affects", "other", "not provided", "confers sensitivity")


def clinvar_category(significance) -> str | None:
    """Normalise a ClinVar significance to P / LP / VUS / LB / B / CONFLICT / None."""
    raw = text(significance).replace("_", " ").lower()
    if not raw:
        return None
    first = re.split(r"[|;]", raw)[0].strip()
    if "conflicting" in first:
        return "CONFLICT"
    if any(word in first for word in _CLINVAR_EXCLUDE):
        return None
    compact = first.replace(" ", "")
    table = {
        "pathogenic": "P",
        "likelypathogenic": "LP",
        "pathogenic/likelypathogenic": "LP",
        "uncertainsignificance": "VUS",
        "vus": "VUS",
        "likelybenign": "LB",
        "benign": "B",
        "benign/likelybenign": "LB",
    }
    return table.get(compact)


def moi_codes(value) -> set[str]:
    """OMIM-style inheritance text -> codes (AD, AR, XL, XLD, XLR, YL, MT, ...)."""
    if isinstance(value, (set, list, tuple)):
        items = value
    else:
        items = re.split(r"[,;/&\s]+", text(value))
    out = set()
    for item in items:
        item = str(item).strip().upper()
        if item:
            out.add(item)
    return out
