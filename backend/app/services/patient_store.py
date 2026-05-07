"""Create / hydrate the per-patient sample directory.

The per-patient layout (post-Phase-1 migration):

    tertiary_output/{LIS_ID}/
      sample_metadata.json     (basic info + reviewer state)
      snv_indel.annotated.tsv  (variant calls)
      analyses/default/
        analysis.json          (hpo + selected_panels + note)

This module owns the "create a brand-new patient" flow: validate the
identifiers, copy/move the TSV into place, parse the phenotype.txt,
write sample_metadata.json and the default analysis.json. The actual
HTTP shell (multipart vs. path pointer) lives in routers/samples.py.
"""
from __future__ import annotations

import json
import re
import shutil
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


def create_new(
    *,
    lis_id: str,
    name: str,
    mrn: str,
    test_type: str = "WES",
    genome_build: str = "hg38",
    category: str = "",
    vcf_path: str = "",
    tsv_src: Path | None = None,
    tsv_bytes: bytes | None = None,
    phenotype_text: str = "",
) -> dict:
    """Create a new patient directory and seed the default analysis.

    Exactly one of `tsv_src` (server-side path) or `tsv_bytes` (uploaded
    blob) must be provided.

    Returns the freshly-loaded sample_metadata.json dict.
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
    if (tsv_src is None) == (tsv_bytes is None):
        raise ValueError("provide exactly one of tsv_src or tsv_bytes")

    sample_dir = TERTIARY_OUTPUT_ROOT / lis_id
    if sample_dir.exists():
        raise FileExistsError(f"sample already exists: {lis_id}")
    sample_dir.mkdir(parents=True)

    # Write the TSV first; if anything below fails we leave the dir
    # behind but at least the variant data is on disk.
    tsv_dst = sample_dir / "snv_indel.annotated.tsv"
    if tsv_src is not None:
        if not tsv_src.is_file():
            shutil.rmtree(sample_dir, ignore_errors=True)
            raise FileNotFoundError(f"tsv source not found: {tsv_src}")
        shutil.copyfile(tsv_src, tsv_dst)
    else:
        tsv_dst.write_bytes(tsv_bytes or b"")

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
