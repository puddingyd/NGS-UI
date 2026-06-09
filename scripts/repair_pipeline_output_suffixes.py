#!/usr/bin/env python3
"""Mirror legacy source-ID tertiary outputs into suffixed UI output dirs.

Older NGS-UI jobs ran Nextflow with the sequencing source sample ID, so
DRAGEN and in-house output could collide under /home/pipeline/tertiary_output.
This script repairs completed cases by copying either:

    PIPELINE_OUT_ROOT/{source_sample_id}/ -> PIPELINE_OUT_ROOT/{ui_sample_id}/

from each UI sample's pipeline_source.json, or by directly treating legacy
VAL-number pipeline output directories as DRAGEN:

    PIPELINE_OUT_ROOT/VAL-57/ -> PIPELINE_OUT_ROOT/VAL-57-dragen/

Only missing destination directories are created. Existing destinations are
left untouched.
"""
from __future__ import annotations

import argparse
import json
import os
import re
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
DEFAULT_LEGACY_DRAGEN_PATTERN = r"^VAL-\d+$"


def iter_repairs_from_sidecars(tertiary_root: Path, pipeline_root: Path):
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
            yield {
                "ui_sid": ui_sid,
                "source_sid": source_sid,
                "pipeline_type": pipeline_type,
                "src": src,
                "dst": dst,
                "reason": "pipeline_source.json",
            }


def iter_repairs_from_pipeline_dirs(pipeline_root: Path, pattern: str):
    if not pattern or not pipeline_root.is_dir():
        return
    rx = re.compile(pattern)
    for src in sorted(pipeline_root.iterdir()):
        if not src.is_dir() or src.name.startswith("_"):
            continue
        source_sid = src.name
        if not rx.search(source_sid):
            continue
        ui_sid = f"{source_sid}-dragen"
        dst = pipeline_root / ui_sid
        if dst.exists():
            continue
        yield {
            "ui_sid": ui_sid,
            "source_sid": source_sid,
            "pipeline_type": "dragen",
            "src": src,
            "dst": dst,
            "reason": f"pipeline-dir-pattern:{pattern}",
        }


def iter_repairs(tertiary_root: Path, pipeline_root: Path, legacy_dragen_pattern: str):
    seen: set[tuple[Path, Path]] = set()
    for item in iter_repairs_from_sidecars(tertiary_root, pipeline_root):
        key = (item["src"], item["dst"])
        if key in seen:
            continue
        seen.add(key)
        yield item
    for item in iter_repairs_from_pipeline_dirs(pipeline_root, legacy_dragen_pattern):
        key = (item["src"], item["dst"])
        if key in seen:
            continue
        seen.add(key)
        yield item


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tertiary-root", type=Path, default=DEFAULT_TERTIARY_OUTPUT_ROOT)
    ap.add_argument("--pipeline-root", type=Path, default=DEFAULT_PIPELINE_OUT_ROOT)
    ap.add_argument(
        "--legacy-dragen-pattern",
        default=DEFAULT_LEGACY_DRAGEN_PATTERN,
        help=(
            "Regex for old unsuffixed DRAGEN pipeline directories to mirror "
            "to '<sample>-dragen'. Use '' to disable. Default: ^VAL-\\d+$"
        ),
    )
    ap.add_argument("--apply", action="store_true", help="perform copies; default is dry-run")
    args = ap.parse_args()

    repairs = list(iter_repairs(args.tertiary_root, args.pipeline_root, args.legacy_dragen_pattern))
    if not repairs:
        print("No legacy pipeline output directories need repair.")
        return 0
    for item in repairs:
        action = "COPY" if args.apply else "DRY-RUN"
        print(
            f"{action}\t{item['pipeline_type']}\t"
            f"{item['source_sid']} -> {item['ui_sid']}\t"
            f"{item['src']} -> {item['dst']}\t{item['reason']}"
        )
        if args.apply:
            shutil.copytree(item["src"], item["dst"], symlinks=True)
    if not args.apply:
        print("\nRe-run with --apply to create the suffixed pipeline output directories.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
