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

from ..services import manual_acmg, panel_deadzone


class OldFormatError(ValueError):
    """Raised when an snv_indel.annotated.tsv is in the pre-2026-05 layout
    (no ACMG_CRITERIA column). Caller surfaces a clear "please re-run
    pipeline" message in the UI instead of trying to render a broken card.
    """

# Reviewer-facing SNV/Indel analysis tiers. Tier 1C is a retrieval bucket:
# predictor triggers surface variants for review but never modify ACMG points
# or classification.
TIERS = ["1A", "1B", "1C", "2"]

PREDICTED_SUSPECT_THRESHOLDS = {
    "pknn": 1.0,
    "alphamissense_moderate": 0.906,
    "alphamissense_supporting": 0.792,
    "bayesdel_moderate": 0.27,
    "bayesdel_supporting": 0.13,
    "pangolin": 0.20,
    "revel_moderate": 0.773,
    "revel_supporting": 0.644,
    "spliceai": 0.20,
}

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

_CONSEQUENCE_RANK = {
    "transcript_ablation": 1,
    "splice_acceptor_variant": 2,
    "splice_donor_variant": 3,
    "stop_gained": 4,
    "frameshift_variant": 5,
    "stop_lost": 6,
    "start_lost": 7,
    "transcript_amplification": 8,
    "inframe_insertion": 9,
    "inframe_deletion": 10,
    "missense_variant": 11,
    "protein_altering_variant": 12,
    "splice_region_variant": 13,
    "splice_donor_5th_base_variant": 14,
    "splice_donor_region_variant": 15,
    "splice_polypyrimidine_tract_variant": 16,
    "incomplete_terminal_codon_variant": 17,
    "stop_retained_variant": 18,
    "synonymous_variant": 19,
    "coding_sequence_variant": 20,
    "mature_miRNA_variant": 21,
    "5_prime_UTR_variant": 22,
    "3_prime_UTR_variant": 23,
    "non_coding_transcript_exon_variant": 24,
    "intron_variant": 25,
    "NMD_transcript_variant": 26,
    "non_coding_transcript_variant": 27,
    "upstream_gene_variant": 28,
    "downstream_gene_variant": 29,
    "intergenic_variant": 38,
}

_TRANSCRIPT_TYPE_RANK = {
    "MANE_SELECT": 0,
    "MANE_PLUS_CLINICAL": 1,
    "CANONICAL": 2,
    "APPRIS_P1": 3,
    "BEST_CONSEQUENCE": 4,
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


def _optional_int(v) -> int | None:
    text = str(v or "").strip()
    if not text or text.upper() in {".", "NA", "N/A"}:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _clean_vep_rank(value: str) -> str:
    text = str(value or "").strip()
    return "" if text.upper() in {"", ".", "-", "NA", "N/A"} else text


def _first_num(v: str):
    """First value of a possibly per-ALT (comma-separated) numeric field, or
    None. Used for the in-house AF columns (Number=A, e.g. '0.44,0.42'); the
    card shows the first ALT, consistent with the first-ALT VAF/AD behaviour.
    Returns None for missing / '.' so absent variants render a clean '—'
    (NOT 0 — that would make the card show '(0/0)')."""
    if v is None or v == "":
        return None
    tok = str(v).split(",", 1)[0].strip()
    if tok in ("", "."):
        return None
    return _to_num(tok)


def _to_bool(v: str) -> bool:
    return str(v).strip().lower() in ("true", "1", "yes", "y", "t")


def _worst_consequence_rank(consequence: str) -> int:
    ranks = [
        _CONSEQUENCE_RANK.get(part.strip(), 99)
        for part in str(consequence or "").split("&")
        if part.strip()
    ]
    return min(ranks) if ranks else 99


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


def _effective_acmg_points(row: dict) -> float | int | None:
    """ACMG score shown by the UI: GeneBe first, then pipeline/legacy."""
    points = _to_num(_coalesce(
        row.get("GENEBE_ACMG_SCORE"),
        row.get("ACMG_SCORE"),
        row.get("ACMG_POINTS"),
    ))
    return points if isinstance(points, (int, float)) else None


def _clingen_vcep_score(criteria_text: str) -> tuple[int | None, list[str]]:
    """Derive reviewer-facing Tavtigian points from ERepo criteria.

    ClinGen ERepo does not publish a score column in the v3.6 TSV contract.
    Keep this explicitly derived from the criteria tokens; the VCEP's supplied
    classification remains authoritative and may differ because a panel can
    apply gene-specific specifications.
    """
    criteria, unknown = manual_acmg.parse_criteria_text(criteria_text)
    if not criteria:
        return None, unknown
    return int(manual_acmg.calculate(criteria)["score"]), unknown


def predicted_suspect_evidence(row: dict) -> dict[str, Any]:
    """Return 1C retrieval triggers without changing ACMG evidence.

    Core predictors exist in the standard tertiary TSV. Extra predictors are
    populated only when Extra VEP had the corresponding local resources.
    Missing scores are treated as not evaluated, never as benign.
    """
    t = PREDICTED_SUSPECT_THRESHOLDS
    points = _effective_acmg_points(row)
    pknn = _max_multi(row.get("PKNN_LLR"))
    alpha = _max_multi(row.get("ALPHAMISSENSE"))
    bayes = _max_multi(row.get("BAYESDEL_NOAF"))
    pangolin = _max_abs_multi(row.get("PANGOLIN_SCORE"))
    revel = _max_multi(row.get("REVEL"))
    spliceai = _max_multi(row.get("SPLICEAI_MAX"))

    acmg_trigger = points is not None and points >= 4
    core_reasons: list[str] = []
    extra_reasons: list[str] = []

    alpha_moderate = alpha is not None and alpha >= t["alphamissense_moderate"]
    bayes_moderate = bayes is not None and bayes >= t["bayesdel_moderate"]
    alpha_supporting = alpha is not None and alpha >= t["alphamissense_supporting"]
    bayes_supporting = bayes is not None and bayes >= t["bayesdel_supporting"]

    if pknn is not None and pknn >= t["pknn"]:
        core_reasons.append(f"P-KNN LLR {pknn:g} ≥ {t['pknn']:g}")
    if alpha_moderate:
        core_reasons.append(
            f"AlphaMissense {alpha:g} ≥ {t['alphamissense_moderate']:g}"
        )
    if bayes_moderate:
        core_reasons.append(f"BayesDel {bayes:g} ≥ {t['bayesdel_moderate']:g}")
    if not alpha_moderate and not bayes_moderate and alpha_supporting and bayes_supporting:
        core_reasons.append(
            f"AlphaMissense {alpha:g} ≥ {t['alphamissense_supporting']:g} + "
            f"BayesDel {bayes:g} ≥ {t['bayesdel_supporting']:g}"
        )
    if pangolin is not None and abs(pangolin) >= t["pangolin"]:
        core_reasons.append(f"|Pangolin {pangolin:g}| ≥ {t['pangolin']:g}")

    if revel is not None and revel >= t["revel_moderate"]:
        extra_reasons.append(f"REVEL {revel:g} ≥ {t['revel_moderate']:g}")
    elif (
        revel is not None
        and revel >= t["revel_supporting"]
        and (alpha_supporting or bayes_supporting)
    ):
        partner = []
        if alpha_supporting:
            partner.append(f"AlphaMissense {alpha:g}")
        if bayes_supporting:
            partner.append(f"BayesDel {bayes:g}")
        extra_reasons.append(
            f"REVEL {revel:g} ≥ {t['revel_supporting']:g} + " + " / ".join(partner)
        )
    if spliceai is not None and spliceai >= t["spliceai"]:
        extra_reasons.append(f"SpliceAI {spliceai:g} ≥ {t['spliceai']:g}")

    reasons: list[str] = []
    if acmg_trigger:
        reasons.append(f"ACMG points {points:g} ≥ 4")
    reasons.extend(f"Core: {reason}" for reason in core_reasons)
    reasons.extend(f"Extra: {reason}" for reason in extra_reasons)
    return {
        "acmg_trigger": acmg_trigger,
        "core_trigger": bool(core_reasons),
        "extra_trigger": bool(extra_reasons),
        "core_reasons": core_reasons,
        "extra_reasons": extra_reasons,
        "reasons": reasons,
    }


def classify_tier(row: dict) -> str:
    """Map one TSV row to a tier (1A / 1B / 1C / 2).

    Per spec:
        1A — ClinVar P/LP ≥ 1★
        1B — Frameshift / nonsense (LOFTEE HC)
        1C — Predicted suspect (ACMG ≥4 or Core/Extra predictor trigger)
        2  — Other

    ClinVar P/LP 0★ and conflicting calls receive no special tier. They may
    still enter 1C through an independent ACMG/predictor trigger; otherwise
    they remain in Other.
    """
    sig = (row.get("CLINVAR_SIG") or "").strip()
    stars = _to_int(row.get("CLINVAR_STARS"), 0)
    loftee_hc = _loftee_hc_call(row)
    is_plp = sig in _PLP_SIGS

    if is_plp and stars >= 1:
        return "1A"
    if loftee_hc:
        return "1B"
    evidence = predicted_suspect_evidence(row)
    if evidence["acmg_trigger"] or evidence["core_trigger"] or evidence["extra_trigger"]:
        return "1C"
    return "2"


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

    ESM1b and SIFT use the opposite direction from most protein-effect
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


def _depth_from_ad_optional(ad: str) -> int | None:
    if not str(ad or "").strip():
        return None
    values: list[int] = []
    for part in str(ad).split(","):
        value = part.strip()
        if not value or value == "." or value.upper() in {"NA", "N/A"}:
            continue
        try:
            values.append(int(float(value)))
        except ValueError:
            return None
    return sum(values) if values else None


def _alt_depth_from_ad(ad: str) -> int | None:
    parts = str(ad or "").split(",")
    if len(parts) < 2:
        return None
    return _optional_int(parts[1])


def _strand_bias_payload(row: dict, ref: str, alt: str) -> dict[str, Any]:
    """Normalize v3.5 STRAND_BIAS while leaving old TSVs badge-free."""
    if "STRAND_BIAS" not in row:
        return {
            "strand_bias_status": "",
            "strand_bias_raw": "",
            "strand_bias_fs": None,
            "strand_bias_sor": None,
            "strand_bias_threshold": "",
        }
    raw = str(row.get("STRAND_BIAS") or "").strip()
    is_indel = len(ref) != len(alt)
    threshold = (
        "Indel: FS>200 or SOR>10.0"
        if is_indel
        else "SNV: FS>60 or SOR>3.0"
    )
    status = "pass" if raw.upper() == "PASS" else "warn" if raw.upper().startswith("WARN") else "manual"
    fs_match = re.search(r"(?:^|[,(])FS=([^,)]*)", raw, flags=re.I)
    sor_match = re.search(r"(?:^|[,(])SOR=([^,)]*)", raw, flags=re.I)
    return {
        "strand_bias_status": status,
        "strand_bias_raw": raw,
        "strand_bias_fs": _to_num(fs_match.group(1)) if fs_match else None,
        "strand_bias_sor": _to_num(sor_match.group(1)) if sor_match else None,
        "strand_bias_threshold": threshold,
    }


def _litvar2_payload(row: dict) -> dict:
    """Return a browser-safe literature payload from post-processing columns."""
    def safe_url(value: object) -> str:
        candidate_url = str(value or "").strip()
        try:
            parsed = urllib.parse.urlparse(candidate_url)
            if (
                parsed.scheme == "https"
                and parsed.hostname == "www.ncbi.nlm.nih.gov"
                and parsed.path == "/research/litvar2/docsum"
            ):
                return candidate_url
        except ValueError:
            pass
        return ""

    def parsed_sources(value: object) -> list[dict]:
        if isinstance(value, str):
            try:
                value = json.loads(value or "[]")
            except (TypeError, json.JSONDecodeError):
                value = []
        sources = []
        seen: set[str] = set()
        for raw_source in value if isinstance(value, list) else []:
            if not isinstance(raw_source, dict):
                continue
            litvar_id = str(
                raw_source.get("litvar_id") or raw_source.get("id") or ""
            ).strip()
            source_url = safe_url(raw_source.get("url"))
            dedupe_key = litvar_id or source_url
            if not dedupe_key or dedupe_key in seen:
                continue
            seen.add(dedupe_key)
            sources.append({
                "id": litvar_id,
                "pmid_count": max(
                    0, int(_to_num(raw_source.get("pmids_count")) or 0),
                ),
                "url": source_url,
            })
        return sources

    pmids: list[str] = []
    for token in re.findall(r"\d+", str(row.get("LITVAR2_PMIDS_TOP5") or "")):
        if token not in pmids:
            pmids.append(token)
        if len(pmids) == 5:
            break
    count = _to_num(row.get("LITVAR2_PMID_COUNT"))
    candidates = []
    try:
        raw_candidates = json.loads(str(row.get("LITVAR2_CANDIDATES_JSON") or "[]"))
    except (TypeError, json.JSONDecodeError):
        raw_candidates = []
    seen_candidates: set[str] = set()
    for raw in raw_candidates if isinstance(raw_candidates, list) else []:
        if not isinstance(raw, dict):
            continue
        litvar_id = str(raw.get("litvar_id") or raw.get("id") or "").strip()
        candidate_url = safe_url(raw.get("url"))
        dedupe_key = litvar_id or candidate_url
        if not dedupe_key or dedupe_key in seen_candidates:
            continue
        seen_candidates.add(dedupe_key)
        candidate_pmids = []
        for value in raw.get("pmids") or []:
            token = str(value).strip()
            if token.isdigit() and token not in candidate_pmids:
                candidate_pmids.append(token)
            if len(candidate_pmids) == 5:
                break
        candidates.append({
            "id": litvar_id,
            "rsid": str(raw.get("rsid") or "").strip(),
            "gene": str(raw.get("gene") or "").strip(),
            "hgvs": str(raw.get("hgvs") or "").strip(),
            "pmid_count": max(0, int(_to_num(raw.get("pmids_count")) or 0)),
            "pmids": candidate_pmids,
            "url": candidate_url,
            "merged_record_count": max(
                1, int(_to_num(raw.get("merged_record_count")) or 1),
            ),
            "source_records": parsed_sources(raw.get("source_records")),
        })
    source_records = parsed_sources(row.get("LITVAR2_SOURCES_JSON") or "[]")
    return {
        "id": str(row.get("LITVAR2_ID") or "").strip(),
        "rsid": str(row.get("LITVAR2_RSID") or "").strip(),
        "pmid_count": max(0, int(count or 0)),
        "pmids": pmids,
        "dataset_date": str(row.get("LITVAR2_DATASET_DATE") or "").strip(),
        "match_method": str(row.get("LITVAR2_MATCH_METHOD") or "").strip(),
        "status": str(row.get("LITVAR2_STATUS") or "").strip(),
        "url": safe_url(row.get("LITVAR2_URL")),
        "merged_record_count": max(1, len(source_records)),
        "source_records": source_records,
        "candidates": candidates,
    }


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
    ad = _coalesce(row.get("AD"), row.get("AD_DV"), row.get("AD_HC"))
    depth = next(
        (
            value
            for value in (
                _optional_int(row.get("DP")),
                _optional_int(row.get("DP_DV")),
                _optional_int(row.get("DP_HC")),
                _depth_from_ad_optional(ad),
            )
            if value is not None
        ),
        None,
    )
    alt_depth = _alt_depth_from_ad(ad)
    vaf_value = _to_num(_coalesce(row.get("VAF"), row.get("VAF_DV"), row.get("VAF_HC")))
    if vaf_value is None:
        vaf_value = _vaf_from_ad(ad)

    hgnc_id = (row.get("HGNC_ID") or "").strip()
    gene, hgnc_id = panel_deadzone.canonical_gene_symbol(row.get("GENE", ""), hgnc_id)
    transcript = row.get("TRANSCRIPT", "")
    ensembl_transcript = transcript
    refseq_transcript = (row.get("REFSEQ_NUC") or "").strip()
    refseq_protein = (row.get("REFSEQ_PROT") or "").strip()
    mane_status = (row.get("MANE_STATUS") or "").strip()
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
        refseq_transcript = refseq_transcript or refseq_nm
        ensembl_transcript = display_mane.get("enst") or ensembl_transcript
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
    display_transcript = refseq_transcript or transcript
    display_hgvs_c = hgvs_c
    if refseq_transcript and display_hgvs_c and ":" in display_hgvs_c:
        display_hgvs_c = f"{refseq_transcript}:{display_hgvs_c.split(':', 1)[1]}"
    hgvs_full = ":".join(p for p in (gene, display_transcript,
                                      _strip_tx(display_hgvs_c),
                                      _strip_tx(hgvs_p)) if p)
    latest_applied = _to_bool(row.get("CLINVAR_LATEST_APPLIED", ""))
    has_separate_latest = "CLINVAR_LATEST_SIG" in row
    if latest_applied and not has_separate_latest:
        # Compatibility with the short-lived overlay format that replaced
        # CLINVAR_* and retained the pipeline values in CLINVAR_BASE_*.
        baseline_sig = _coalesce(row.get("CLINVAR_BASE_SIG"))
        baseline_stars = _to_num(row.get("CLINVAR_BASE_STARS"))
        baseline_dn = _coalesce(row.get("CLINVAR_BASE_DN"))
        baseline_sigconf = _coalesce(row.get("CLINVAR_BASE_SIGCONF"))
        baseline_variation_id = _coalesce(row.get("CLINVAR_BASE_VARIATION_ID"))
        latest_sig = _coalesce(row.get("CLINVAR_SIG"))
        latest_stars = _to_num(row.get("CLINVAR_STARS"))
        latest_dn = _coalesce(row.get("CLINVAR_DN"))
        latest_sigconf = _coalesce(
            row.get("CLINVAR_CONF"), row.get("CLINVAR_SIGCONF")
        )
        latest_variation_id = _coalesce(row.get("CLINVAR_VARIATION_ID"))
        tier_row = dict(row)
        tier_row.update({
            "CLINVAR_SIG": baseline_sig,
            "CLINVAR_STARS": "" if baseline_stars is None else str(baseline_stars),
            "CLINVAR_DN": baseline_dn,
            "CLINVAR_SIGCONF": baseline_sigconf,
            "CLINVAR_VARIATION_ID": baseline_variation_id,
        })
    else:
        baseline_sig = _coalesce(row.get("CLINVAR_SIG"))
        baseline_stars = _to_num(row.get("CLINVAR_STARS"))
        baseline_dn = _coalesce(row.get("CLINVAR_DN"))
        baseline_sigconf = _coalesce(
            row.get("CLINVAR_CONF"), row.get("CLINVAR_SIGCONF")
        )
        baseline_variation_id = _coalesce(row.get("CLINVAR_VARIATION_ID"))
        latest_sig = _coalesce(row.get("CLINVAR_LATEST_SIG"))
        latest_stars = _to_num(row.get("CLINVAR_LATEST_STARS"))
        latest_dn = _coalesce(row.get("CLINVAR_LATEST_DN"))
        latest_sigconf = _coalesce(row.get("CLINVAR_LATEST_SIGCONF"))
        latest_variation_id = _coalesce(row.get("CLINVAR_LATEST_VARIATION_ID"))
        tier_row = row
    predicted_evidence = predicted_suspect_evidence(row)
    vcep_criteria = _coalesce(row.get("CLINGEN_VCEP_CRITERIA"))
    vcep_score, vcep_unknown = _clingen_vcep_score(vcep_criteria)
    variant = {
        "id": vid,
        "CHROM": chrom,
        "POS": pos,
        "REF": ref,
        "ALT": alt,
        "gene_symbol": gene,
        "transcript": transcript,
        "ensembl_transcript": ensembl_transcript,
        "refseq_transcript": refseq_transcript,
        "refseq_protein": refseq_protein,
        "mane_status": mane_status,
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
        "exon":   _clean_vep_rank(row.get("EXON", "")),
        "intron": _clean_vep_rank(row.get("INTRON", "")),
        # Old pipeline emits single AD/VAF; new pipeline splits per caller
        # (AD_DV/AD_HC, VAF_DV/VAF_HC). DV's VAF is more reliable for
        # heteroplasmy estimation, so prefer it. HC-only calls (no
        # FORMAT/VAF in HaplotypeCaller) get VAF derived from AD so the
        # card doesn't show '—' just because the column is missing.
        "AD": ad,
        # Total read depth at the position — DV column preferred (more
        # accurate for short reads), HC fallback. Used by the UI for
        # depth-based filtering (drop < 20 on WES) and the WGS/TITAN
        # low-depth red flag (< 10). Falls back to sum(AD) when neither
        # DP column is set.
        "depth": depth,
        "alt_depth": alt_depth,
        "low_alt_support": bool(alt_depth is not None and alt_depth < 10),
        "alt_af": vaf_value,
        "CLNSIG": baseline_sig,
        "clinvar_stars": baseline_stars,
        "clinvar_dn": baseline_dn,
        "CLNSIGCONF": baseline_sigconf,
        "CLNSIG_old": baseline_sig,
        "CLNSIGCONF_old": baseline_sigconf,
        "clinvar_stars_old": baseline_stars,
        "clinvar_dn_old": baseline_dn,
        "clinvar_variation_id_old": baseline_variation_id,
        "clinvar_latest_sig": latest_sig,
        "clinvar_latest_sigconf": latest_sigconf,
        "clinvar_latest_stars": latest_stars,
        "clinvar_latest_dn": latest_dn,
        "clinvar_latest_variation_id": latest_variation_id,
        "clinvar_latest_review_status": _coalesce(
            row.get("CLINVAR_LATEST_REVIEW_STATUS")
        ),
        "clinvar_latest_applied": latest_applied,
        "clinvar_change": (row.get("CLINVAR_CHANGE") or "").strip(),
        "AF": _to_num(row.get("GNOMAD_G_AF")),
        "AF_eas": _to_num(row.get("GNOMAD_G_EAS_AF")),
        "AF_exome": _to_num(row.get("GNOMAD_E_AF")),
        "AF_exome_eas": _to_num(row.get("GNOMAD_E_EAS_AF")),
        # In-house AF from our own WGS cohort (scripts/annotate_inhouse_af.py
        # writes INHOUSE_AC/AN/AF). The card shows an AF_nckuh row "<AF> (AC/AN)".
        # All None when the variant isn't in the DB (or the DB isn't deployed).
        "inhouse_af":          _first_num(row.get("INHOUSE_AF")),
        "inhouse_ac":          _first_num(row.get("INHOUSE_AC")),
        "inhouse_an":          _first_num(row.get("INHOUSE_AN")),
        # New pipeline drops Taiwan Biobank; gnomAD EAS covers the same
        # population well enough.
        "TG_eas_af":           _to_num(row.get("TG_EAS_AF")),
        # In-silico predictors. VEP can emit per-transcript scores
        # joined by '&' (e.g. AlphaMissense '.&0.9482&0.9432') — take
        # the worst case (max). Categorical _PRED columns get the first
        # non-empty value.
        "PKNN_LLR":            _max_multi(row.get("PKNN_LLR")),
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
        # REVEL / MutPred2 / VEST4 / CADD come from the v3.6 Nextflow
        # dbNSFP branch.  SpliceAI remains a research-only post-processing
        # overlay.  Legacy MetaRNN is retained in the payload for old cases,
        # but is no longer part of the reviewer display order.
        "SpliceAI_score":      _max_multi(row.get("SPLICEAI_MAX")),
        "MetaRNN_score":       _max_multi(row.get("METARNN")),
        "REVEL_score":         _max_multi(row.get("REVEL")),
        "CADD_score":          _max_multi(row.get("CADD_PHRED")),
        "MutPred2_score":      _max_multi(row.get("MUTPRED2")),
        "MutPred2_pred":       _first_str(row.get("MUTPRED2_PRED")),
        "VEST4_score":         _max_multi(row.get("VEST4")),
        # Others — under ▾ More on the card.
        "DANN":                _max_multi(row.get("DANN")),
        "PhactBoost":          _max_multi(row.get("PHACTBOOST")),
        "PhyloP":              _max_multi(row.get("PHYLOP100")),
        "GERP":                _max_multi(row.get("GERP")),
        "SIFT_score":          _min_multi(row.get("SIFT")),
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
        # ClinGen Evidence Repository VCEP expert assessment.  This is a
        # comparison source only; it never enters the effective ACMG cascade
        # unless a reviewer explicitly Applies and saves its criteria.
        "clingen_vcep_class": _normalize_acmg_class(
            _coalesce(row.get("CLINGEN_VCEP_CLASS"))
        ),
        "clingen_vcep_criteria": vcep_criteria.replace("|", ","),
        "clingen_vcep_panel": _coalesce(row.get("CLINGEN_VCEP_PANEL")),
        "clingen_vcep_agreement": _coalesce(row.get("CLINGEN_AGREEMENT")),
        "clingen_vcep_score": vcep_score,
        "clingen_vcep_unknown_criteria": vcep_unknown,
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
        # GIAB genome-stratification labels (annotate_giab_strata.py):
        # comma-separated short labels for difficult regions the variant
        # overlaps (homopolymer / tandem_repeat / segdup / ...). Frontend
        # renders one QC badge per label.
        "giab_strata": [s for s in (row.get("GIAB_STRATA", "") or "").split(",") if s.strip()],
        "OMIM_link": row.get("OMIM_LINK", ""),
        "gnomAD_link": row.get("GNOMAD_LINK", ""),
        "ClinVar_link": row.get("CLINVAR_LINK", ""),
        "litvar2": _litvar2_payload(row),
        "report_class": row.get("REPORT_CLASS", ""),
        "tier": classify_tier(tier_row),
        # Preserve the non-ACMG 1C triggers so a later manual ACMG overlay can
        # recalculate the tier without losing independent predictor evidence.
        "predicted_suspect_non_acmg": bool(
            predicted_evidence["core_trigger"] or predicted_evidence["extra_trigger"]
        ),
    }
    variant.update(_strand_bias_payload(row, ref, alt))
    return variant


def _transcript_option_from_variant(v: dict) -> dict:
    key_parts = [
        v.get("gene_symbol") or "",
        v.get("transcript") or "",
        v.get("refseq_transcript") or "",
        v.get("ensembl_transcript") or "",
        v.get("transcript_type") or "",
        v.get("HGVS_C") or "",
        v.get("HGVS_P") or "",
        v.get("Consequence") or "",
    ]
    key = "|".join(str(part) for part in key_parts)
    return {
        "key": key,
        "gene_symbol": v.get("gene_symbol") or "",
        "transcript": v.get("transcript") or "",
        "ensembl_transcript": v.get("ensembl_transcript") or "",
        "refseq_transcript": v.get("refseq_transcript") or "",
        "refseq_protein": v.get("refseq_protein") or "",
        "mane_status": v.get("mane_status") or "",
        "transcript_type": v.get("transcript_type") or "",
        "HGVS_C": v.get("HGVS_C") or "",
        "HGVS_P": v.get("HGVS_P") or "",
        "HGVS": v.get("HGVS") or "",
        "Consequence": v.get("Consequence") or "",
        "impact": v.get("impact") or "",
        "exon": v.get("exon") or "",
        "intron": v.get("intron") or "",
        "hgnc_id": v.get("hgnc_id") or "",
    }


def _transcript_option_sort_key(opt: dict) -> tuple[int, int, str]:
    return (
        _worst_consequence_rank(opt.get("Consequence") or ""),
        _TRANSCRIPT_TYPE_RANK.get(str(opt.get("transcript_type") or "").upper(), 9),
        str(opt.get("key") or ""),
    )


def _apply_transcript_option(v: dict, opt: dict) -> None:
    for src_key, dst_key in (
        ("gene_symbol", "gene_symbol"),
        ("transcript", "transcript"),
        ("ensembl_transcript", "ensembl_transcript"),
        ("refseq_transcript", "refseq_transcript"),
        ("refseq_protein", "refseq_protein"),
        ("mane_status", "mane_status"),
        ("transcript_type", "transcript_type"),
        ("HGVS_C", "HGVS_C"),
        ("HGVS_P", "HGVS_P"),
        ("HGVS", "HGVS"),
        ("Consequence", "Consequence"),
        ("impact", "impact"),
        ("exon", "exon"),
        ("intron", "intron"),
        ("hgnc_id", "hgnc_id"),
    ):
        v[dst_key] = opt.get(src_key) or ""
    v["selected_transcript_key"] = opt.get("key") or ""
    v["default_transcript_key"] = opt.get("key") or ""


def merge_snv_variant_row(variants: dict[str, dict], v: dict) -> dict:
    """Merge one transcript row into the one-card-per-variant payload."""
    option = _transcript_option_from_variant(v)
    v["transcript_options"] = [option]
    v["selected_transcript_key"] = option["key"]
    v["default_transcript_key"] = option["key"]

    existing = variants.get(v["id"])
    if existing is None:
        variants[v["id"]] = v
        return v

    options = list(existing.get("transcript_options") or [])
    seen = {opt.get("key") for opt in options}
    if option["key"] not in seen:
        options.append(option)
    options.sort(key=_transcript_option_sort_key)
    existing["transcript_options"] = options
    best = options[0]
    if _transcript_option_sort_key(option) < _transcript_option_sort_key({
        "key": existing.get("default_transcript_key") or "",
        "transcript_type": existing.get("transcript_type") or "",
        "Consequence": existing.get("Consequence") or "",
    }):
        for key, value in v.items():
            if key not in {"transcript_options", "selected_transcript_key", "default_transcript_key"}:
                existing[key] = value
    _apply_transcript_option(existing, best)
    return existing


# Depth thresholds applied at the adapter level.
#   WES → hard floor at DP=20 (drop low-confidence WES calls).
#   WGS → no hard floor (uniform coverage); red flag at DP<10 only.
# Both flags land on the variant as `low_depth: bool`; the UI paints
# the Read-depth cell red when true.
_DP_HARD_FLOOR_WES = 20
_DP_LOW_FLAG_WES   = 20
_DP_LOW_FLAG_WGS   = 10


def _is_mito_chrom(chrom: str) -> bool:
    value = (chrom or "").strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value.upper() in {"M", "MT"}


def load_snv_tsv(tsv_path: Path,
                 *,
                 test_type: str = "WES") -> tuple[dict[str, dict], dict[str, list[str]]]:
    """Read snv_indel.annotated.tsv → (variants, categories).

    `variants` is keyed by chr-pos-ref-alt id.
    `categories` is keyed by tier (1A / 1B / 1C / 2) → ordered ids.
    The adapter keeps its established ACMG score/id ordering; sample loading
    re-sorts every tier by the final total score after phenotype enrichment.

    `test_type` controls the depth gate:
      WES  → variants with DP < 20 dropped entirely.
      WGS  → no hard floor; low_depth set when DP < 10.
    """
    variants: dict[str, dict] = {}
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
            if _is_mito_chrom(row.get("CHROM") or ""):
                continue
            canonical_gene, _ = panel_deadzone.canonical_gene_symbol(
                row.get("GENE", ""),
                row.get("HGNC_ID", ""),
            )
            if canonical_gene in EXCLUDED_GENES:
                continue
            v = _row_to_variant(row)
            dp = v.get("depth")
            if is_wes and (dp is None or dp < _DP_HARD_FLOOR_WES):
                continue
            v["low_depth"] = bool(dp is not None and dp < low_flag)
            merge_snv_variant_row(variants, v)

    by_tier: dict[str, list[tuple[float, str]]] = {t: [] for t in TIERS}
    for vid, v in variants.items():
        pts = v.get("ACMG_score")
        sort_key = float(pts) if isinstance(pts, (int, float)) else -999.0
        by_tier.setdefault(v.get("tier"), []).append((sort_key, vid))
    categories: dict[str, list[str]] = {}
    for t in TIERS:
        by_tier[t].sort(key=lambda kv: (-kv[0], kv[1]))
        categories[t] = [vid for _, vid in by_tier[t]]

    return variants, categories
