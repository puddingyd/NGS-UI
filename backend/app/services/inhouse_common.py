"""院內常見 (in-house common) SNVs: marking, hiding and a ranking penalty.

Two in-house sources say a variant is common in our own population:

* the WES recurrence table (``acmg_engine.recurrence``): the fraction of WES
  patients carrying the variant, counted only for WES cases and only once the
  table holds ``WES_MIN_PATIENTS`` patients (the case's own patient is left
  out of both counts);
* the in-house WGS AF (``INHOUSE_AF/AN``, about 1,400 genomes), available for
  every case.

Display: the main screen hides a variant (unless rescued; see ``app.js``) when
``AF_NCKUH ≥ 0.05 & AC ≥ 50`` (unchanged) or when ≥ 5 % of WES patients carry
it.  Ranking: the points used for ``total_score`` drop by 4 (level 1) or 8
(level 2); the ACMG class and the tier are not changed.

    WES patients carrying   ≥ 5 % → 4 points   ≥ 10 % → 8 points
    in-house WGS AF         ≥ 2.5 % → 4 points  ≥ 5 % → 8 points   (AN ≥ 2000)

WGS AF is an allele frequency; heterozygous carriers are about twice the AF,
so 2.5 % / 5 % AF line up with 5 % / 10 % carriers.  The higher of the two
levels is used.  BS1/BA1 already in the effective ACMG points count towards
the target (``extra = max(0, target − BS1 4 / BA1 8 deduction)``), because the
in-house AF is part of the max AF used for BS1/BA1.  Reviewer ACMG (manual)
and ClinVar ≥2★ P/LP are exempt.
"""
from __future__ import annotations

import re

WES_MIN_PATIENTS = 100
WES_LEVELS = ((0.10, 2), (0.05, 1))
WGS_LEVELS = ((0.05, 2), (0.025, 1))
WGS_MIN_AN = 2000
LEVEL_POINTS = {0: 0, 1: 4, 2: 8}
HIDE_WES_FRACTION = 0.05

_SNV_ID = re.compile(r"^(?:chr)?[0-9XYxy]{1,2}-\d+-[ACGTNacgtn]+-[ACGTNacgtn]+$")
_PLP = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic", "likely pathogenic/pathogenic"}


def _num(value):
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if out == out else None


def _recurrence_key(variant_id: str) -> str:
    return variant_id if variant_id.lower().startswith("chr") else "chr" + variant_id


def _level(value, levels) -> int:
    if value is None:
        return 0
    for threshold, level in levels:
        if value >= threshold:
            return level
    return 0


def _wes_counts(variants: dict, sample_id: str, mrn: str) -> dict[str, tuple[int, int]]:
    try:
        from .acmg_engine import runner
        store = runner.shared_recurrence()
        if store is None:
            return {}
        patient = runner.case_patient(store, sample_id, mrn)
        keys = {vid: _recurrence_key(vid) for vid in variants if _SNV_ID.match(str(vid))}
        counts = store.lookup_many(set(keys.values()), exclude_patient=patient or None)
        total = store.patient_count() - (1 if patient and store.has_patient(patient) else 0)
    except Exception as exc:  # noqa: BLE001 — marking never blocks loading
        print(f"[inhouse-common] WES recurrence unavailable sample={sample_id}: {type(exc).__name__}: {exc}",
              flush=True)
        return {}
    if total < WES_MIN_PATIENTS:
        return {}
    return {vid: (counts.get(key, 0), total) for vid, key in keys.items()}


def annotate(variants: dict, *, sample_id: str, test_type: str, mrn: str = "") -> None:
    """Add ``wes_recurrence`` (WES cases) and ``inhouse_common`` (level / basis)."""
    if not variants:
        return
    wes = _wes_counts(variants, sample_id, mrn) if str(test_type or "").upper().startswith("WES") else {}
    for vid, variant in variants.items():
        if not isinstance(variant, dict):
            continue
        variant.pop("wes_recurrence", None)
        variant.pop("inhouse_common", None)
        basis = []
        level = 0
        hit = wes.get(vid)
        if hit:
            count, total = hit
            fraction = count / total if total else 0.0
            variant["wes_recurrence"] = {"count": count, "patients": total, "fraction": round(fraction, 4)}
            wes_level = _level(fraction, WES_LEVELS)
            if wes_level:
                basis.append(f"院內 WES {count}/{total} 位病人（{fraction:.1%}）")
                level = max(level, wes_level)
        af, an = _num(variant.get("inhouse_af")), _num(variant.get("inhouse_an"))
        if af is not None and an is not None and an >= WGS_MIN_AN:
            wgs_level = _level(af, WGS_LEVELS)
            if wgs_level:
                basis.append(f"院內 WGS AF {af:.3g}（AN {int(an)}）")
                level = max(level, wgs_level)
        if level:
            variant["inhouse_common"] = {"level": level, "target": LEVEL_POINTS[level], "basis": "；".join(basis),
                                         "hide_wes": bool(hit and hit[0] / hit[1] >= HIDE_WES_FRACTION)}


def frequency_deduction(criteria_text) -> int:
    """Points BA1/BS1 already take off the effective ACMG score."""
    from . import manual_acmg
    parsed, _unknown = manual_acmg.parse_criteria_text(criteria_text)
    out = 0
    for code in ("BA1", "BS1"):
        evidence = parsed.get(code)
        if evidence:
            out = max(out, manual_acmg._POINTS.get(evidence.get("strength"), 0))
    return out


def _clinvar_protected(variant: dict) -> bool:
    sig = str(variant.get("CLNSIG") or "").strip().lower().replace("_", " ")
    stars = _num(variant.get("clinvar_stars"))
    return sig in _PLP and stars is not None and stars >= 2


def apply_penalty(variant: dict, *, scope: str, score, criteria_text) -> float | None:
    """Ranking points after the penalty; records ``inhouse_common.penalty`` / ``exempt``."""
    info = variant.get("inhouse_common")
    points = _num(score)
    if not isinstance(info, dict):
        return points
    info.pop("exempt", None)
    info["penalty"] = 0
    if scope in ("sample", "global"):
        info["exempt"] = "判讀者已手動判讀 ACMG"
        return points
    if _clinvar_protected(variant):
        info["exempt"] = "ClinVar ≥2★ P/LP"
        return points
    already = frequency_deduction(criteria_text)
    extra = max(0, int(info.get("target") or 0) - already)
    info["already"] = already
    info["penalty"] = extra
    if not extra:
        return points
    return (points or 0.0) - extra
