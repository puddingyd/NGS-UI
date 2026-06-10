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
import urllib.parse
from pathlib import Path
from typing import Any

from ..services import panel_deadzone


class OldFormatError(ValueError):
    """Raised when an snv_indel.annotated.tsv is in the pre-2026-05 layout
    (no ACMG_CRITERIA column). Caller surfaces a clear "please re-run
    pipeline" message in the UI instead of trying to render a broken card.
    """

# Tier categories defined in 三級輸出計畫.md §2.3 Page 2.
TIERS = ["1A", "1B", "1C", "2", "3"]

_PLP_SIGS = {
    "Pathogenic",
    "Likely_pathogenic",
    "Pathogenic/Likely_pathogenic",
    "Likely_pathogenic/Pathogenic",
}

# Genes whose SNV/Indel calls are NOT clinically meaningful from this
# pipeline — they're either repeat-expansion disease genes (where the
# CAG/CTG/etc. count is what matters, separate STR pipeline) or HLA
# (haplotype-level interpretation needed). Hide them entirely from
# the SNV/Indel cards.
EXCLUDED_GENES = {
    # CAG / SCA repeat-expansion disorders
    "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7",
    "HTT", "TBP", "ZFHX3", "THAP11", "JPH3", "PABPN1",
    # HLA — haplotype-level, separate workflow
    "HLA-A", "HLA-B", "HLA-C",
    "HLA-DQA1", "HLA-DQB1", "HLA-DRB1",
    "HLA-DPA1", "HLA-DPB1",
}

# Conflicting ClinVar interpretations are only worth surfacing in
# tier 2 when at least one submitter called P / LP. CLNSIGCONF like
# "Likely_benign(2)|Uncertain_significance(1)" alone is noise.
_CONFLICTING_PATHO_RE = re.compile(r"\b(?:Likely_)?[Pp]athogenic\b")


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


def _parse_mane_all(raw: str) -> list[dict[str, Any]]:
    """Parse MANE_ALL across pipeline/TSV quoting variants.

    Some historical exports kept CSV-style doubled quotes inside the TSV
    field, and copy/paste from spreadsheets can introduce smart quotes.
    Keep this forgiving so a malformed MANE_ALL cell does not erase the
    entire MANE detail table in the UI.
    """
    text = (raw or "").strip()
    if not text:
        return []
    quote_normalized = (
        text.replace("\u201c", '"')
            .replace("\u201d", '"')
            .replace("\u2018", "'")
            .replace("\u2019", "'")
    )
    candidates = [text]
    if quote_normalized != text:
        candidates.append(quote_normalized)
    for candidate in list(candidates):
        stripped = candidate.strip()
        if len(stripped) >= 2 and stripped[0] == stripped[-1] == '"':
            candidates.append(stripped[1:-1].replace('""', '"'))
        if '""' in stripped:
            candidates.append(stripped.replace('""', '"'))
    for candidate in candidates:
        try:
            value = json.loads(candidate)
        except json.JSONDecodeError:
            continue
        if isinstance(value, list):
            return [r for r in value if isinstance(r, dict)]
    return []


def _mane_impact_rank(impact: str) -> int:
    return {
        "HIGH": 4,
        "MODERATE": 3,
        "LOW": 2,
        "MODIFIER": 1,
    }.get((impact or "").strip().upper(), 0)


def _pick_display_mane(mane_all: list[dict[str, str]], enst_base: str) -> dict[str, str] | None:
    select_rows = [
        m for m in mane_all
        if (m.get("transcript_type") or "").upper() == "MANE_SELECT"
        and re.match(r"^NM_", m.get("transcript") or "", flags=re.I)
    ]
    if not select_rows:
        return None

    def score(row: dict[str, str]) -> tuple[int, int, int, int]:
        row_enst = (row.get("enst") or "").split(".")[0]
        hgvs_signal = int(bool(row.get("hgvs_c") or row.get("hgvs_p")))
        coding_signal = int(_mane_impact_rank(row.get("impact") or "") >= 2)
        return (
            hgvs_signal,
            coding_signal,
            _mane_impact_rank(row.get("impact") or ""),
            int(bool(enst_base and row_enst == enst_base)),
        )

    return max(select_rows, key=score)


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
    # ACMG_SCORE priority mirrors the UI's display priority: prefer
    # the GeneBe second-opinion score when present, fall back to the
    # pipeline's classifier (ACMG_SCORE) and then to the legacy
    # ACMG_POINTS column. Keeps tier classification consistent with
    # what the reviewer sees on the card.
    points = _to_num(_coalesce(row.get("GENEBE_ACMG_SCORE"),
                                row.get("ACMG_SCORE"),
                                row.get("ACMG_POINTS"))) or 0
    if isinstance(points, (int, float)) and points >= 4:
        return "1C"
    if is_plp and stars == 0:
        return "2"
    if is_conflicting:
        # CLNSIGCONF format: "Likely_benign(2)|Uncertain_significance(1)" etc.
        # Only count as tier 2 when at least one submitter weighed in
        # with P / LP — pure VUS+LB conflict isn't actionable.
        sigconf = (row.get("CLINVAR_SIGCONF") or row.get("CLINVAR_CONF") or "")
        if _CONFLICTING_PATHO_RE.search(sigconf):
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


def _max_multi(v) -> float | None:
    """Max numeric value across a `&`-separated multi-value cell.

    VEP emits per-transcript / per-consequence scores joined with `&`
    (e.g. AlphaMissense `.&0.9482&0.9432`). Take the worst-case (max)
    so the card surfaces the most pathogenic prediction for the locus.
    Returns None when no part parses to a number.
    """
    if v is None:
        return None
    s = str(v).strip()
    if not s or s in (".", "NA", "N/A"):
        return None
    best: float | None = None
    for part in s.split("&"):
        p = part.strip()
        if not p or p in (".", "NA", "N/A"):
            continue
        try:
            x = float(p)
        except ValueError:
            continue
        if best is None or x > best:
            best = x
    return best


def _min_multi(v) -> float | None:
    """Min numeric value across a `&`-separated multi-value cell.

    ESM1b uses the opposite direction from the other protein-effect
    predictors: lower scores are more damaging. Take the minimum so
    multi-transcript rows surface the most pathogenic prediction.
    """
    if v is None:
        return None
    s = str(v).strip()
    if not s or s in (".", "NA", "N/A"):
        return None
    best: float | None = None
    for part in s.split("&"):
        p = part.strip()
        if not p or p in (".", "NA", "N/A"):
            continue
        try:
            x = float(p)
        except ValueError:
            continue
        if best is None or x < best:
            best = x
    return best


def _first_str(v) -> str:
    """First non-blank/non-NA part of a `&`-separated cell."""
    if v is None:
        return ""
    s = str(v).strip()
    if not s or s in (".", "NA", "N/A"):
        return ""
    for part in s.split("&"):
        p = part.strip()
        if p and p not in (".", "NA", "N/A"):
            return p
    return ""


def _max_abs_multi(v) -> float | None:
    """Like _max_multi but picks the value with the largest |x| —
    preserves sign. Useful for signed splice scores (Pangolin), where
    negative = splice loss and positive = splice gain, and what matters
    clinically is the magnitude of the predicted splice change.
    """
    if v is None:
        return None
    s = str(v).strip()
    if not s or s in (".", "NA", "N/A"):
        return None
    best: float | None = None
    for part in s.split("&"):
        p = part.strip()
        if not p or p in (".", "NA", "N/A"):
            continue
        try:
            x = float(p)
        except ValueError:
            continue
        if best is None or abs(x) > abs(best):
            best = x
    return best


def _vaf_from_ad(ad: str) -> float | None:
    """Derive first-ALT VAF from an AD string like '4,6' or '4,6,0'.

    HaplotypeCaller VCFs don't ship a FORMAT/VAF field, so HC-only
    calls (or DV-call-with-blank-VAF rows) reach the card with
    VAF_DV/VAF_HC both '.'. AD is always present though — compute
    alt/(ref+alt+other_alts) so reviewers see something instead of '—'.
    """
    if not ad:
        return None
    parts = [p.strip() for p in str(ad).split(",")]
    nums = []
    for p in parts:
        if not p or p == "." or p.upper() in ("NA", "N/A"):
            continue
        try:
            nums.append(int(p))
        except ValueError:
            return None
    if len(nums) < 2:
        return None
    total = sum(nums)
    if total <= 0:
        return None
    return nums[1] / total  # first ALT vs all (ref + ALTs)


def _depth_from_ad(ad: str) -> int:
    """Fallback DP = sum(AD parts) when neither DP nor DP_DV/DP_HC is
    populated. Returns 0 on blank / malformed input.
    """
    if not ad:
        return 0
    total = 0
    for p in str(ad).split(","):
        p = p.strip()
        if not p or p == "." or p.upper() in ("NA", "N/A"):
            continue
        try:
            total += int(float(p))
        except ValueError:
            continue
    return total


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

    hgnc_id = (row.get("HGNC_ID") or "").strip()
    gene, hgnc_id = panel_deadzone.canonical_gene_symbol(row.get("GENE", ""), hgnc_id)
    transcript = row.get("TRANSCRIPT", "")
    # New pipeline ships HGVS_P with URL-encoded characters (e.g.
    # `p.Gly282%3D` for synonymous changes). Decode so the UI renders
    # `p.Gly282=` instead of the percent escape.
    hgvs_c = urllib.parse.unquote(row.get("HGVS_C", ""))
    hgvs_p = urllib.parse.unquote(row.get("HGVS_P", ""))

    def _mane_first(r: dict, *keys: str) -> str:
        for key in keys:
            v = r.get(key)
            if v is not None and str(v).strip():
                return str(v).strip()
        return ""

    # Normalise key names. Different pipeline revisions used tx/refseq/nm
    # for the RefSeq transcript and enst/ensembl_transcript for ENST.
    mane_all = [
        {
            "transcript": _mane_first(
                r, "tx", "transcript", "refseq", "refseq_transcript", "nm", "NM"
            ),
            "enst": _mane_first(
                r, "enst", "ensembl", "ensembl_transcript", "transcript_id"
            ),
            "transcript_type": _mane_first(r, "type", "transcript_type"),
            "consequence": _mane_first(r, "consequence"),
            "hgvs_c": urllib.parse.unquote(_mane_first(r, "hgvsc", "hgvs_c", "HGVS_C")),
            "hgvs_p": urllib.parse.unquote(_mane_first(r, "hgvsp", "hgvs_p", "HGVS_P")),
            "impact": _mane_first(r, "impact"),
        }
        for r in _parse_mane_all(row.get("MANE_ALL") or "[]")
    ]

    # Reviewers prefer RefSeq accessions for continuity with old reports.
    # Prefer the most informative MANE_SELECT row. When VEP picked a
    # downstream/MODIFIER ENST with blank HGVS fields, this lets gene search
    # cards show the clinically useful NM_* missense/coding transcript.
    enst_base = (transcript or "").split(".")[0]
    display_mane = _pick_display_mane(mane_all, enst_base)
    if display_mane:
        refseq_nm = display_mane["transcript"]
        transcript = refseq_nm
        mane_hgvs_c = display_mane.get("hgvs_c") or ""
        mane_hgvs_p = display_mane.get("hgvs_p") or ""
        if mane_hgvs_c:
            hgvs_c = mane_hgvs_c
        if mane_hgvs_p:
            hgvs_p = mane_hgvs_p
        if hgvs_c and ":" in hgvs_c:
            hgvs_c = f"{refseq_nm}:{hgvs_c.split(':', 1)[1]}"
        if hgvs_p and ":" in hgvs_p:
            hgvs_p = hgvs_p.split(":", 1)[1]

    # Build the combined display HGVS. New pipeline ships hgvs_c with
    # its own transcript prefix (`TX:c.xxx`); strip that here so we
    # don't end up with `GENE:NM_xxx:NM_xxx:c.xxx`. Empty / '.' parts
    # are dropped so trailing `:.` doesn't appear when HGVS.p is
    # absent (e.g. UTR / synonymous / non-coding variants).
    def _strip_tx(s):
        if not s:
            return ""
        t = s.strip()
        if not t or t == "." or t.upper() in ("NA", "N/A"):
            return ""
        return t.split(":", 1)[1] if ":" in t else t
    hgvs_full = ":".join(p for p in (gene, transcript,
                                      _strip_tx(hgvs_c),
                                      _strip_tx(hgvs_p)) if p)

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
        # heteroplasmy estimation, so prefer it. HC-only calls (no
        # FORMAT/VAF in HaplotypeCaller) get VAF derived from AD so the
        # card doesn't show '—' just because the column is missing.
        "AD":     _coalesce(row.get("AD"),  row.get("AD_DV"),  row.get("AD_HC")),
        # Total read depth at the position — DV column preferred (more
        # accurate for short reads), HC fallback. Used by the UI for
        # depth-based filtering (drop < 10 on WES) and the low-depth
        # red flag (< 20). Falls back to sum(AD) when neither DP column
        # is set.
        "depth": (
            _to_int(row.get("DP")) or
            _to_int(row.get("DP_DV")) or
            _to_int(row.get("DP_HC")) or
            _depth_from_ad(_coalesce(row.get("AD"),
                                     row.get("AD_DV"),
                                     row.get("AD_HC")))
        ),
        "alt_af": (_to_num(_coalesce(row.get("VAF"), row.get("VAF_DV"), row.get("VAF_HC")))
                   or _vaf_from_ad(_coalesce(row.get("AD"), row.get("AD_DV"), row.get("AD_HC")))),
        "CLNSIG": row.get("CLINVAR_SIG", ""),
        "clinvar_stars": _to_num(row.get("CLINVAR_STARS")),
        "clinvar_dn": row.get("CLINVAR_DN", ""),
        "CLNSIGCONF": _coalesce(row.get("CLINVAR_CONF"),
                                 row.get("CLINVAR_SIGCONF")),
        "AF": _to_num(row.get("GNOMAD_G_AF")),
        "AF_eas": _to_num(row.get("GNOMAD_G_EAS_AF")),
        "AF_exome": _to_num(row.get("GNOMAD_E_AF")),
        "AF_exome_eas": _to_num(row.get("GNOMAD_E_EAS_AF")),
        # New pipeline drops Taiwan Biobank; gnomAD EAS covers the same
        # population well enough.
        "TG_eas_af":           _to_num(row.get("TG_EAS_AF")),
        # In-silico predictors. VEP can emit per-transcript scores
        # joined by '&' (e.g. AlphaMissense '.&0.9482&0.9432') — take
        # the worst case (max). Categorical _PRED columns get the first
        # non-empty value.
        "PKNN_LLR":            _to_num(row.get("PKNN_LLR")),
        "PKNN_evidence":       _first_str(row.get("PKNN_EVIDENCE")),
        "AlphaMissense_score": _max_multi(row.get("ALPHAMISSENSE")),
        "AlphaMissense_pred":  _first_str(row.get("ALPHAMISSENSE_PRED")),
        "Pangolin_score":      _max_abs_multi(row.get("PANGOLIN_SCORE")),
        "Pangolin_detail":     (row.get("PANGOLIN_DETAIL") or "").strip(),
        # ESM1b is the model the new pipeline uses (NOT ESM2 — different
        # protein language model). Payload key renamed to match.
        "ESM1b_score":         _min_multi(row.get("ESM1B")),
        "ESM1b_pred":          _first_str(row.get("ESM1B_PRED")),
        "VARITY_R":            _max_multi(row.get("VARITY_R")),
        "BayesDel":            _max_multi(row.get("BAYESDEL_NOAF")),
        "BayesDel_pred":       _first_str(row.get("BAYESDEL_NOAF_PRED")),
        # SpliceAI + MetaRNN come from annotate_extra_vep stop-gap (new
        # pipeline doesn't ship them).
        "SpliceAI_score":      _max_multi(row.get("SPLICEAI_MAX")),
        "MetaRNN_score":       _max_multi(row.get("METARNN")),
        # Others — under ▾ More on the card.
        "DANN":                _max_multi(row.get("DANN")),
        "PhactBoost":          _max_multi(row.get("PHACTBOOST")),
        "PhyloP":              _max_multi(row.get("PHYLOP100")),
        "GERP":                _max_multi(row.get("GERP")),
        "SIFT_score":          _max_multi(row.get("SIFT")),
        "SIFT_pred":           _first_str(row.get("SIFT_PRED")),
        "LOFTOOL":             _max_multi(row.get("LOFTOOL")),
        "loftee_hc":           _coalesce(row.get("LOFTEE_HC"), row.get("LOFTEE")),
        "loftee_filter":       row.get("LOFTEE_FILTER", ""),
        "loftee_flags":        row.get("LOFTEE_FLAGS", ""),
        # ACMG — pipeline-computed (acmg_classifier.py).
        "ACMG_criteria":       (_coalesce(row.get("ACMG_CRITERIA"),
                                          row.get("ACMG_EVIDENCE"))
                                .replace("|", ",")),
        "ACMG_score":          _to_num(_coalesce(row.get("ACMG_SCORE"),
                                                 row.get("ACMG_POINTS"))),
        "ACMG_classification": _normalize_acmg_class(row.get("ACMG_CLASS", "")),
        "ACMG_notes":          (row.get("ACMG_NOTES") or "").strip(),
        # GeneBe — populated by annotate_acmg_genebe.py as a SECOND
        # opinion (does NOT overwrite the pipeline's ACMG columns).
        # All four blank when the variant isn't in GeneBe.
        "genebe_acmg_class":   _normalize_acmg_class(row.get("GENEBE_ACMG_CLASS", "")),
        "genebe_acmg_score":   _to_num(row.get("GENEBE_ACMG_SCORE")),
        "genebe_acmg_criteria": (row.get("GENEBE_ACMG_CRITERIA") or "").strip(),
        "genebe_acmg_notes":   (row.get("GENEBE_ACMG_NOTES") or "").strip(),
        # "Score" pill: pipeline ACMG_SCORE rescaled to 0-100 (same
        # transform regardless of which source the class came from).
        # Same GeneBe-first cascade as classify_tier above so the Score
        # pill on the card matches the tier the variant lands in.
        "geno_score":          _acmg_to_geno_score(_to_num(_coalesce(
                                   row.get("GENEBE_ACMG_SCORE"),
                                   row.get("ACMG_SCORE"),
                                   row.get("ACMG_POINTS")))),
        # New display fields from the 65-col pipeline output.
        "rs_id":               (row.get("RS_ID") or "").strip(),
        "hgnc_id":             hgnc_id,
        "disease_associated":  panel_deadzone.is_disease_associated_gene(gene, hgnc_id),
        "clinvar_variation_id": (row.get("CLINVAR_VARIATION_ID") or "").strip(),
        "omim_ids":            (row.get("OMIM_IDS") or "").strip(),
        "domains":             (row.get("DOMAINS") or "").strip(),
        "swissprot":           (row.get("SWISSPROT") or "").strip(),
        "impact":              (row.get("IMPACT") or "").strip(),
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


# Depth thresholds applied at the adapter level.
#   WES → hard floor at DP=10 (drop noise outside the capture targets)
#         + red flag at DP<20.
#   WGS → no hard floor (uniform coverage); red flag at DP<10 only.
# Both flags land on the variant as `low_depth: bool`; the UI paints
# the Read-depth cell red when true.
_DP_HARD_FLOOR_WES = 10
_DP_LOW_FLAG_WES   = 20
_DP_LOW_FLAG_WGS   = 10


def load_snv_tsv(tsv_path: Path,
                 *,
                 test_type: str = "WES") -> tuple[dict[str, dict], dict[str, list[str]]]:
    """Read snv_indel.annotated.tsv → (variants, categories).

    `variants` is keyed by chr-pos-ref-alt id.
    `categories` is keyed by tier (1A / 1B / 2 / 3 / 4 / 5) → ordered ids
    (within tier 4: ACMG_POINTS desc; within tier 5: ACMG_POINTS desc).

    `test_type` controls the depth gate:
      WES  → variants with DP < 10 dropped entirely (off-target noise);
             low_depth set on the remainder when DP < 20.
      WGS  → no hard floor; low_depth set when DP < 10.
    """
    variants: dict[str, dict] = {}
    by_tier: dict[str, list[tuple[float, str]]] = {t: [] for t in TIERS}
    is_wes    = (test_type or "").upper() == "WES"
    low_flag  = _DP_LOW_FLAG_WES if is_wes else _DP_LOW_FLAG_WGS

    with tsv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        # New (2026-05+) pipeline emits ACMG_CRITERIA; the legacy
        # version doesn't. Reject up-front so the UI shows a clear
        # "please re-run pipeline" message instead of a half-rendered
        # card from a misread schema.
        if "ACMG_CRITERIA" not in (reader.fieldnames or []):
            raise OldFormatError(
                "snv_indel.annotated.tsv 缺 ACMG_CRITERIA 欄 — "
                "為舊格式，請以新版 pipeline 重跑此樣本。"
            )
        for row in reader:
            canonical_gene, _ = panel_deadzone.canonical_gene_symbol(
                row.get("GENE", ""),
                row.get("HGNC_ID", ""),
            )
            if canonical_gene in EXCLUDED_GENES:
                continue
            v = _row_to_variant(row)
            dp = v.get("depth") or 0
            if is_wes and dp < _DP_HARD_FLOOR_WES:
                continue
            v["low_depth"] = bool(dp and dp < low_flag)
            variants[v["id"]] = v
            pts = v.get("ACMG_score")
            sort_key = float(pts) if isinstance(pts, (int, float)) else -999.0
            by_tier[v["tier"]].append((sort_key, v["id"]))

    categories: dict[str, list[str]] = {}
    for t in TIERS:
        by_tier[t].sort(key=lambda kv: (-kv[0], kv[1]))
        categories[t] = [vid for _, vid in by_tier[t]]

    return variants, categories
