#!/usr/bin/env python3
"""Mirror legacy source-ID tertiary outputs into suffixed UI output dirs.

Older NGS-UI jobs ran Nextflow with the sequencing source sample ID, so
DRAGEN and in-house output could collide under /home/pipeline/tertiary_output.
This script uses each UI sample's pipeline_source.json to repair completed
cases by copying:

    PIPELINE_OUT_ROOT/{source_sample_id}/ -> PIPELINE_OUT_ROOT/{ui_sample_id}/

Only missing destination directories are created. Existing destinations are
left untouched.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

NGS_UI_HOME = Path(os.environ.get("NGS_UI_HOME") or REPO_ROOT)
DEFAULT_TERTIARY_OUTPUT_ROOT = Path(
    os.environ.get("TERTIARY_OUTPUT_ROOT") or (NGS_UI_HOME / "tertiary_output")
)
DEFAULT_PIPELINE_OUT_ROOT = Path(
    os.environ.get("NGS_UI_PIPELINE_OUT_ROOT") or "/home/pipeline/tertiary_output"
)


def iter_repairs(tertiary_root: Path, pipeline_root: Path):
    for sample_dir in sorted(tertiary_root.iterdir() if tertiary_root.is_dir() else []):
        if not sample_dir.is_dir() or sample_dir.name.startswith("_"):
            continue
        sidecar = sample_dir / "pipeline_source.json"
        if not sidecar.is_file():
            continue
        try:
            data = json.loads(sidecar.read_text(encoding="utf-8")) or {}
        except json.JSONDecodeError:
            continue
        ui_sid = sample_dir.name
        source_sid = str(data.get("source_sample_id") or "").strip()
        pipeline_type = str(data.get("pipeline_type") or "").strip()
        if not source_sid or source_sid == ui_sid:
            continue
        if pipeline_type not in {"dragen", "inhouse"}:
            continue
        src = pipeline_root / source_sid
        dst = pipeline_root / ui_sid
        if src.is_dir() and not dst.exists():
            yield ui_sid, source_sid, pipeline_type, src, dst


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tertiary-root", type=Path, default=DEFAULT_TERTIARY_OUTPUT_ROOT)
    ap.add_argument("--pipeline-root", type=Path, default=DEFAULT_PIPELINE_OUT_ROOT)
    ap.add_argument("--apply", action="store_true", help="perform copies; default is dry-run")
    args = ap.parse_args()

    repairs = list(iter_repairs(args.tertiary_root, args.pipeline_root))
    if not repairs:
        print("No legacy pipeline output directories need repair.")
        return 0
    for ui_sid, source_sid, pipeline_type, src, dst in repairs:
        action = "COPY" if args.apply else "DRY-RUN"
        print(f"{action}\t{pipeline_type}\t{source_sid} -> {ui_sid}\t{src} -> {dst}")
        if args.apply:
            shutil.copytree(src, dst, symlinks=True)
    if not args.apply:
        print("\nRe-run with --apply to create the suffixed pipeline output directories.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
