"""Adapter: snv_indel.annotated.tsv → webdata-shape JSON for the UI.

The legacy frontend (ported from vcf-analysis-hg38-R) expects a
sample-level dict with `meta`, `patient_phenotype`, `variants` (keyed by
chr-pos-ref-alt id) and `categories` (keyed by category name → list of
variant ids). This module reads the new tertiary TSV and shapes the data
so the existing render code keeps working with minimal changes.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any

# Tier categories defined in 三級輸出計畫.md §2.3 Page 2.
TIERS = ["1A", "1B", "1C", "2", "3"]

_PLP_SIGS = {
    "Pathogenic",
    "Likely_pathogenic",
    "Pathogenic/Likely_pathogenic",
    "Likely_pathogenic/Pathogenic",
}


def _to_num(v: str):
    if v is None or v == "":
        return None
    try:
        f = float(v)
        return int(f) if f.is_integer() else f
    except ValueError:
        return v


def _to_int(v: str, default: int = 0) -> int:
    if v is None or v == "":
        return default
    try:
        return int(float(v))
    except ValueError:
        return default


def _to_bool(v: str) -> bool:
    return str(v).strip().lower() in ("true", "1", "yes", "y", "t")


def _loftee_hc_call(row: dict) -> bool:
    """Old pipeline emitted LOFTEE_HC ('HC' or ''); new pipeline emits
    a single LOFTEE column ('HC' / 'LC' / '.'). Accept either.
    """
    raw = _coalesce(row.get("LOFTEE_HC"), row.get("LOFTEE"))
    return raw.strip().upper() == "HC"


def classify_tier(row: dict) -> str:
    """Map one TSV row to a tier (1A / 1B / 1C / 2 / 3).

    Per spec:
        1A — ClinVar P/LP ≥ 1★
        1B — Frameshift / nonsense (LOFTEE HC)
        1C — ACMG points ≥ 4 (strong-evidence VUS+)
        2  — ClinVar P/LP 0★ or Conflicting (含 P)
        3  — 其餘 (ACMG points < 4)
    """
    sig = (row.get("CLINVAR_SIG") or "").strip()
    stars = _to_int(row.get("CLINVAR_STARS"), 0)
    loftee_hc = _loftee_hc_call(row)
    is_plp = sig in _PLP_SIGS
    is_conflicting = "Conflicting" in sig

    if is_plp and stars >= 1:
        return "1A"
    if loftee_hc:
        return "1B"
    points = _to_num(row.get("ACMG_POINTS")) or 0
    if isinstance(points, (int, float)) and points >= 4:
        return "1C"
    if (is_plp and stars == 0) or is_conflicting:
        return "2"
    return "3"


def _acmg_to_geno_score(acmg_points) -> int | None:
    """Linear-map ACMG_POINTS (clamped to [-10, 10]) onto 0-100.

    Mirrors the LIRICAL compositeLR-to-pheno-score transform so the
    variant card's "Score" line speaks one consistent 0-100 scale.
    """
    if acmg_points is None:
        return None
    try:
        x = float(acmg_points)
    except (TypeError, ValueError):
        return None
    x = max(-10.0, min(10.0, x))
    return int(round((x + 10.0) / 20.0 * 100.0))


# Canonical 5-tier ACMG class strings the frontend's <select> uses as
# option values. The source TSV is inconsistent ("VUS",
# "Uncertain_significance", lowercase, …), so normalise here — anything
# we don't recognise (e.g. stray evidence-code strings like
# "BP4_Strong|BA1" that leaked into this column upstream) is passed
# through verbatim and the UI just shows it as "—".
_ACMG_CLASS_CANON = {
    "pathogenic":             "Pathogenic",
    "likely pathogenic":      "Likely pathogenic",
    "likely_pathogenic":      "Likely pathogenic",
    "uncertain significance": "Uncertain significance",
    "uncertain_significance": "Uncertain significance",
    "vus":                    "Uncertain significance",
    "likely benign":          "Likely benign",
    "likely_benign":          "Likely benign",
    "benign":                 "Benign",
}


def _normalize_acmg_class(raw: str) -> str:
    key = (raw or "").strip().lower()
    return _ACMG_CLASS_CANON.get(key, (raw or "").strip())


def _coalesce(*vals: str) -> str:
    """First non-blank / non-NA value, '' otherwise."""
    for v in vals:
        s = (v or "").strip()
        if s and s not in (".", "NA", "N/A"):
            return s
    return ""


def _row_to_variant(row: dict) -> dict:
    """Reshape one TSV row into the per-variant dict the frontend expects.

    Field mapping is the inverse of scripts/convert_old_json_to_tertiary_tsv.py
    (with the new spec-only fields added so the UI can display them).
    """
    chrom = row["CHROM"]
    pos = _to_int(row["POS"])
    ref = row["REF"]
    alt = row["ALT"]
    vid = f"{chrom}-{pos}-{ref}-{alt}"

    gene = row.get("GENE", "")
    transcript = row.get("TRANSCRIPT", "")
    hgvs_c = row.get("HGVS_C", "")
    hgvs_p = row.get("HGVS_P", "")
    hgvs_full = ":".join(p for p in (gene, transcript, hgvs_c, hgvs_p) if p)

    try:
        mane_all = json.loads(row.get("MANE_ALL") or "[]")
    except json.JSONDecodeError:
        mane_all = []

    return {
        "id": vid,
        "CHROM": chrom,
        "POS": pos,
        "REF": ref,
        "ALT": alt,
        "gene_symbol": gene,
        "transcript": transcript,
        "transcript_type": row.get("TRANSCRIPT_TYPE", ""),
        "HGVS_C": hgvs_c,
        "HGVS_P": hgvs_p,
        "HGVS": hgvs_full,
        "Consequence": row.get("CONSEQUENCE", ""),
        "MANE_ALL": mane_all,
        "callers": row.get("CALLERS", ""),
        "zygosity": row.get("ZYGOSITY", ""),
        "GT_DV": row.get("GT_DV", ""),
        "GT_HC": row.get("GT_HC", ""),
        "exon":   row.get("EXON", ""),
        "intron": row.get("INTRON", ""),
        # Old pipeline emits single AD/VAF; new pipeline splits per caller
        # (AD_DV/AD_HC, VAF_DV/VAF_HC). DV's VAF is more reliable for
        # heteroplasmy estimation, so prefer it.
        "AD":     _coalesce(row.get("AD"),  row.get("AD_DV"),  row.get("AD_HC")),
        "alt_af": _to_num(
                      _coalesce(row.get("VAF"), row.get("VAF_DV"), row.get("VAF_HC"))
                  ),
        "CLNSIG": row.get("CLINVAR_SIG", ""),
        "clinvar_stars": _to_num(row.get("CLINVAR_STARS")),
        "clinvar_dn": row.get("CLINVAR_DN", ""),
        "CLNSIGCONF": _coalesce(row.get("CLINVAR_CONF"),
                                 row.get("CLINVAR_SIGCONF")),
        "AF": _to_num(row.get("GNOMAD_G_AF")),
        "AF_eas": _to_num(row.get("GNOMAD_G_EAS_AF")),
        "AF_exome": _to_num(row.get("GNOMAD_E_AF")),
        "AF_exome_eas": _to_num(row.get("GNOMAD_E_EAS_AF")),
        "TaiwanBioBank": _to_num(row.get("TWB_AF")),
        "PKNN_LLR": _to_num(row.get("PKNN_LLR")),
        "REVEL": _to_num(row.get("REVEL")),
        "BayesDel": _to_num(_coalesce(row.get("BAYESDEL"),
                                       row.get("BAYESDEL_NOAF"))),
        "AlphaMissense_score": _to_num(row.get("ALPHAMISSENSE")),
        "MetaRNN_score": _to_num(row.get("METARNN")),
        "ESM2_score": _to_num(_coalesce(row.get("ESM2_SCORE"),
                                         row.get("ESM1B"))),
        "Evo2_score": _to_num(row.get("EVO2_SCORE")),
        "SpliceAI_score": _to_num(_coalesce(row.get("SPLICEAI_MAX"),
                                             row.get("PANGOLIN_SCORE"))),
        "CADD_phred": _to_num(row.get("CADD_PHRED")),
        "loftee_hc": _coalesce(row.get("LOFTEE_HC"), row.get("LOFTEE")),
        "loftee_filter": row.get("LOFTEE_FILTER", ""),
        "loftee_flags": row.get("LOFTEE_FLAGS", ""),
        "ACMG_criteria": (row.get("ACMG_EVIDENCE") or "").replace("|", ","),
        "ACMG_score": _to_num(row.get("ACMG_POINTS")),
        "ACMG_classification": _normalize_acmg_class(row.get("ACMG_CLASS", "")),
        # Variant score for the "Score" pill: ACMG_POINTS rescaled 0-100.
        "geno_score": _acmg_to_geno_score(_to_num(row.get("ACMG_POINTS"))),
        "phase_group": row.get("PHASE_GROUP", ""),
        "phase_result": row.get("PHASE_RESULT", ""),
        "in_roh": _to_bool(row.get("IN_ROH", "")),
        "in_panel": _to_bool(row.get("IN_PANEL", "")),
        "in_blacklist": _to_bool(row.get("IN_BLACKLIST", "")),
        "OMIM_link": row.get("OMIM_LINK", ""),
        "gnomAD_link": row.get("GNOMAD_LINK", ""),
        "ClinVar_link": row.get("CLINVAR_LINK", ""),
        "report_class": row.get("REPORT_CLASS", ""),
        "tier": classify_tier(row),
    }


def load_snv_tsv(tsv_path: Path) -> tuple[dict[str, dict], dict[str, list[str]]]:
    """Read snv_indel.annotated.tsv → (variants, categories).

    `variants` is keyed by chr-pos-ref-alt id.
    `categories` is keyed by tier (1A / 1B / 2 / 3 / 4 / 5) → ordered ids
    (within tier 4: ACMG_POINTS desc; within tier 5: ACMG_POINTS desc).
    """
    variants: dict[str, dict] = {}
    by_tier: dict[str, list[tuple[float, str]]] = {t: [] for t in TIERS}

    with tsv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            v = _row_to_variant(row)
            variants[v["id"]] = v
            pts = v.get("ACMG_score")
            sort_key = float(pts) if isinstance(pts, (int, float)) else -999.0
            by_tier[v["tier"]].append((sort_key, v["id"]))

    categories: dict[str, list[str]] = {}
    for t in TIERS:
        by_tier[t].sort(key=lambda kv: (-kv[0], kv[1]))
        categories[t] = [vid for _, vid in by_tier[t]]

    return variants, categories
