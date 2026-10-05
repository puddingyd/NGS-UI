"""Normalize displayed NCKUH chrX zygosity from the saved ploidy VCF."""
from __future__ import annotations

import json

from . import ploidy, sample_layout


# GRCh38 pseudoautosomal intervals, 1-based and inclusive.
_X_PAR = ((10_001, 2_781_479), (155_701_383, 156_030_895))


def _is_x_nonpar(chrom: str, pos: str, ref: str) -> bool:
    if chrom.strip().lower() not in {"x", "chrx"}:
        return False
    try:
        start = int(pos)
    except (TypeError, ValueError):
        return False
    if start < 1 or not ref or ref == ".":
        return False
    end = start + len(ref) - 1
    return not any(start <= par_end and end >= par_start for par_start, par_end in _X_PAR)


def _nckuh_xy_sample(sample_id: str) -> bool:
    ploidy_result = ploidy.load_sample_ploidy(sample_id)
    if ploidy_result.get("karyotype") != "XY":
        return False
    source_path = sample_layout.state_file(sample_id, "pipeline_source.json")
    try:
        source = json.loads(source_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        source = {}
    mode = (
        str(source.get("pipeline_type") or "").strip().lower()
        if isinstance(source, dict) else ""
    )
    if mode:
        return mode in {"inhouse", "nckuh"}
    kind = str(ploidy_result.get("pipeline_kind") or "").lower()
    if kind in {"nckuh", "dragen"}:
        return kind == "nckuh"
    return sample_id.lower().endswith(("-nckuh", "-inhouse"))


def normalize_loaded_variants(variants: dict[str, dict], sample_id: str) -> int:
    """Correct existing NCKUH cards/reports on read, without a tertiary rerun."""
    if not variants or not _nckuh_xy_sample(sample_id):
        return 0
    corrected = 0
    for variant in variants.values():
        if (
            not variant.get("somatic")
            and str(variant.get("zygosity") or "").strip().lower()
            in {"hom", "homozygous"}
            and not variant.get("haploid_het")
            and not variant.get("haploid_het_callers")
            and _is_x_nonpar(
                str(variant.get("CHROM") or ""),
                str(variant.get("POS") or ""),
                str(variant.get("REF") or ""),
            )
        ):
            variant["zygosity"] = "hemi"
            corrected += 1
    return corrected
