"""Hydrate a per-patient sample directory the pipeline already produced.

The tertiary pipeline lands variant calls in
    tertiary_output/{LIS_ID}/snv_indel.annotated.tsv
on its own. The 載入新個案 flow attaches reviewer-side info on top:
basic identifiers, an empty default analysis, optionally a parsed
copy of the phenotype.txt. After hydration the directory looks like:

    tertiary_output/{LIS_ID}/
      snv_indel.annotated.tsv  (untouched; pipeline output)
      sample_metadata.json     (basic info + empty reviewer state)
      analyses/default/
        analysis.json          (hpo + selected_panels + note)
        {LIS_ID}_{MRN}_phenotype.txt   (audit copy, when provided)

Refusal cases:
  * lis_id directory missing or no TSV → 404 / 400 from the router
  * sample_metadata.json already present → 409 (already registered)
"""
from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path

from ..config import TERTIARY_OUTPUT_ROOT
from . import analyses_store, phenotype_io


_LIS_ID_RE = re.compile(r"^[A-Za-z0-9_-]{1,32}$")
_TEST_TYPES = {"WES", "WGS"}
_GENOME_BUILDS = {"hg19", "hg38"}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _validate_lis_id(lis_id: str) -> None:
    if not _LIS_ID_RE.match(lis_id or ""):
        raise ValueError(
            "lis_id must match [A-Za-z0-9_-]{1,32} (used as directory name)"
        )


def sample_exists(lis_id: str) -> bool:
    return (TERTIARY_OUTPUT_ROOT / lis_id).is_dir()


def is_registered(lis_id: str) -> bool:
    return (TERTIARY_OUTPUT_ROOT / lis_id / "sample_metadata.json").is_file()


def register(
    *,
    lis_id: str,
    name: str,
    mrn: str,
    test_type: str = "WES",
    genome_build: str = "hg38",
    category: str = "",
    vcf_path: str = "",
    phenotype_text: str = "",
) -> dict:
    """Attach reviewer-side info to a pipeline-produced directory.

    The directory tertiary_output/{lis_id}/ must already exist with at
    least snv_indel.annotated.tsv inside (the pipeline drops it there).
    Refuses if the dir is already registered (sample_metadata.json
    present).
    """
    _validate_lis_id(lis_id)
    if not name:
        raise ValueError("name is required")
    if not mrn:
        raise ValueError("mrn is required")
    if test_type not in _TEST_TYPES:
        raise ValueError(f"test_type must be one of {sorted(_TEST_TYPES)}")
    if genome_build not in _GENOME_BUILDS:
        raise ValueError(f"genome_build must be one of {sorted(_GENOME_BUILDS)}")

    sample_dir = TERTIARY_OUTPUT_ROOT / lis_id
    if not sample_dir.is_dir():
        raise FileNotFoundError(
            f"pipeline directory not found: {sample_dir} "
            "(tertiary pipeline drops the TSV here; nothing to register yet)"
        )
    if not (sample_dir / "snv_indel.annotated.tsv").is_file():
        raise FileNotFoundError(
            f"snv_indel.annotated.tsv missing under {sample_dir}"
        )
    if (sample_dir / "sample_metadata.json").is_file():
        raise FileExistsError(f"sample already registered: {lis_id}")

    # Parse phenotype.txt → hpo + panels for the default analysis.
    hpo, panels = phenotype_io.parse(phenotype_text or "")

    # Seed sample_metadata.json with basic info + empty reviewer state.
    now = _now()
    meta = {
        "sample_id":            lis_id,
        "lis_id":               lis_id,
        "name":                 name,
        "mrn":                  mrn,
        "test_type":            test_type,
        "genome_build":         genome_build,
        "category":             category or "",
        "vcf_path":             vcf_path or "",
        "run_date":             now,
        "active_analysis":      "default",
        "clinical_description": "",
        "comment":              "",
        "tags":                 [],
        "status":               {},
        "edits":                {},
        "panels":               {},
        "manual_variants":      [],
        "created_at":           now,
        "updated_at":           now,
    }
    (sample_dir / "sample_metadata.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    # Default analysis.json + audit copy of the parsed phenotype.txt.
    analyses_store.write_version(lis_id, "default", hpo=hpo, panels=panels)
    if hpo or panels:
        phenotype_io.write(
            hpo, panels,
            analyses_store.version_dir(lis_id, "default")
            / f"{lis_id}_{mrn}_phenotype.txt",
        )

    return meta
