#!/usr/bin/env python3
"""Repoint every sample's vcf_path to the convention-driven location.

Convention (post-this-change):
    tertiary_output/{LIS_ID}/{LIS_ID}.from_tsv.vcf.gz

For each sample dir under tertiary_output/ with sample_metadata.json:
  1. (re)generate the from-TSV VCF if missing or older than the TSV.
  2. Set sample_metadata.vcf_path to the canonical absolute path.

Idempotent. Safe to run multiple times. Existing original VCFs in
NGS_UI/vcf/ are NOT touched — they sit untouched as a separate
archive; the tertiary tools just stop reading them.

Usage:
    PYTHONPATH=backend NGS_UI_HOME=/home/n102968/NGS_UI \\
        python3 scripts/migrate_vcf_path.py
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

# Ensure the package is importable regardless of cwd.
HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent.parent / "backend"))

from app.config import TERTIARY_OUTPUT_ROOT  # noqa: E402
from app.services import vcf_writer           # noqa: E402


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def migrate_one(sample_dir: Path) -> str:
    sid = sample_dir.name
    meta_path = sample_dir / "sample_metadata.json"
    if not meta_path.is_file():
        return "no-meta"
    if not (sample_dir / "snv_indel.annotated.tsv").is_file():
        return "no-tsv"

    # (Re)generate the VCF if missing or stale.
    if vcf_writer.needs_rebuild(sid):
        vcf_writer.from_tsv(sid)

    canonical = str(vcf_writer.vcf_path_for(sid))

    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    if not isinstance(meta, dict):
        return "meta-bad"
    if meta.get("vcf_path") == canonical:
        return "already-canonical"
    meta["vcf_path"] = canonical
    meta["updated_at"] = _now()
    meta_path.write_text(
        json.dumps(meta, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return "repointed"


def main() -> int:
    if not TERTIARY_OUTPUT_ROOT.is_dir():
        print(f"!! {TERTIARY_OUTPUT_ROOT} not found", file=sys.stderr)
        return 1
    counts: dict[str, int] = {}
    for sub in sorted(TERTIARY_OUTPUT_ROOT.iterdir()):
        if not sub.is_dir() or sub.name.startswith("_"):
            continue
        result = migrate_one(sub)
        counts[result] = counts.get(result, 0) + 1
        print(f"  {result:20} {sub.name}")
    print(f"\nDone: {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
