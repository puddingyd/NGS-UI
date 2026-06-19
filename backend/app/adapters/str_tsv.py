"""STR TSV adapter for the analysis card."""
from __future__ import annotations

import csv
from pathlib import Path

STR_TIERS = ["STR-P", "STR-I", "STR-N"]


def _clean(value) -> str:
    s = "" if value is None else str(value).strip()
    return "" if not s or s.upper() in {"NA", "N/A", "."} else s


def _to_float(value):
    s = _clean(value)
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _tier_for(classification: str) -> str:
    c = (classification or "").strip().lower()
    if c == "pathogenic":
        return "STR-P"
    if c in {"intermediate", "borderline"}:
        return "STR-I"
    return "STR-N"


def _distance_score(row: dict) -> float:
    a1 = _to_float(row.get("REPCN_A1")) or 0.0
    a2 = _to_float(row.get("REPCN_A2")) or 0.0
    max_repeat = max(a1, a2)
    pathogenic_min = _to_float(row.get("PATHOGENIC_MIN"))
    intermediate_min = _to_float(row.get("INTERMEDIATE_MIN"))
    threshold = pathogenic_min or intermediate_min
    if threshold and threshold > 0:
        return max_repeat / threshold
    return max_repeat


def _row_to_variant(row: dict) -> dict:
    a1 = _to_float(row.get("REPCN_A1"))
    a2 = _to_float(row.get("REPCN_A2"))
    max_repeat = max([v for v in (a1, a2) if v is not None], default=None)
    return {
        "CHROM": _clean(row.get("CHROM")),
        "POS": _clean(row.get("POS")),
        "END": _clean(row.get("END")),
        "STR_ID": _clean(row.get("STR_ID")),
        "GENE": _clean(row.get("GENE")),
        "MOTIF": _clean(row.get("MOTIF")),
        "LOCUS_STRUCTURE": _clean(row.get("LOCUS_STRUCTURE")),
        "TYPE": _clean(row.get("TYPE")),
        "REPCN_A1": _clean(row.get("REPCN_A1")),
        "REPCN_A2": _clean(row.get("REPCN_A2")),
        "max_repeat": max_repeat,
        "DP": _clean(row.get("DP")),
        "REPCI": _clean(row.get("REPCI")),
        "BENIGN_MIN": _clean(row.get("BENIGN_MIN")),
        "BENIGN_MAX": _clean(row.get("BENIGN_MAX")),
        "PATHOGENIC_MIN": _clean(row.get("PATHOGENIC_MIN")),
        "PATHOGENIC_MAX": _clean(row.get("PATHOGENIC_MAX")),
        "INTERMEDIATE_MIN": _clean(row.get("INTERMEDIATE_MIN")),
        "INTERMEDIATE_MAX": _clean(row.get("INTERMEDIATE_MAX")),
        "CLASSIFICATION": _clean(row.get("CLASSIFICATION")) or "normal",
        "DISEASE": _clean(row.get("DISEASE")),
        "INHERITANCE": _clean(row.get("INHERITANCE")),
        "PIPELINE": _clean(row.get("PIPELINE")),
        "_sort_score": _distance_score(row),
    }


def load_str_tsv(path: Path) -> tuple[dict, dict]:
    variants: dict[str, dict] = {}
    categories = {tier: [] for tier in STR_TIERS}
    if not path.exists():
        return variants, categories

    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for idx, row in enumerate(reader, start=1):
            variant = _row_to_variant(row)
            gene = variant.get("GENE") or "STR"
            str_id = variant.get("STR_ID") or f"row{idx}"
            chrom = variant.get("CHROM") or ""
            pos = variant.get("POS") or ""
            vid = f"str-{str_id}-{chrom}-{pos}"
            tier = _tier_for(variant.get("CLASSIFICATION") or "")
            variants[vid] = variant
            categories.setdefault(tier, []).append(vid)

    rank = {"STR-P": 0, "STR-I": 1, "STR-N": 2}
    for tier, ids in categories.items():
        ids.sort(
            key=lambda vid: (
                rank.get(tier, 9),
                -float(variants[vid].get("_sort_score") or 0),
                variants[vid].get("GENE") or "",
                variants[vid].get("STR_ID") or "",
            )
        )
        for vid in ids:
            variants[vid].pop("_sort_score", None)
    return variants, categories
