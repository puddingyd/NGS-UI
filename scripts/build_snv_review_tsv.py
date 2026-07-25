#!/usr/bin/env python3
"""Build the compact main-screen TSV from snv_indel.annotated.tsv."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app.services.snv_review import ensure_review_tsv  # noqa: E402


def _infer_test_type(raw_tsv: Path, output_dir: Path | None = None) -> str:
    directory = output_dir or raw_tsv.parent
    sample_id = directory.parent.name if directory.name == "08_postprocessing" else ""
    candidates = (
        directory / f"{sample_id}.sample_metadata.json",
        directory / "sample_metadata.json",
    ) if sample_id else (directory / "sample_metadata.json",)
    for meta_path in candidates:
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        value = str(meta.get("test_type") or "").upper()
        return value if value in {"WES", "WGS"} else "WES"
    return "WES"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv", required=True, help="complete snv_indel.annotated.tsv")
    ap.add_argument("--output-dir", type=Path, help="directory for derived review TSV")
    ap.add_argument("--output-path", type=Path, help="exact review TSV output path")
    ap.add_argument("--manifest-path", type=Path, help="exact review manifest output path")
    ap.add_argument("--overlay", type=Path, help="sparse SNV annotation overlay SQLite")
    ap.add_argument(
        "--test-type",
        choices=["WES", "WGS"],
        help="Apply WES/WGS-specific review TSV filters. Defaults to sample_metadata.json, then WES.",
    )
    args = ap.parse_args()

    raw_tsv = Path(args.tsv).resolve()
    if not raw_tsv.is_file():
        print(f"ERROR: --tsv 找不到：{raw_tsv}", file=sys.stderr)
        return 2
    output_dir = args.output_dir.resolve() if args.output_dir else None
    output_path = args.output_path.resolve() if args.output_path else None
    manifest_path = args.manifest_path.resolve() if args.manifest_path else None
    if output_dir and output_dir.name == "08_postprocessing":
        sample_id = output_dir.parent.name
        output_path = output_path or output_dir / f"{sample_id}.snv_indel.review.tsv"
        manifest_path = (
            manifest_path
            or output_dir / f"{sample_id}.snv_indel.review.tsv.source.json"
        )
    test_type = args.test_type or _infer_test_type(raw_tsv, output_dir)
    review_tsv = ensure_review_tsv(
        raw_tsv,
        test_type=test_type,
        output_dir=output_dir,
        output_path=output_path,
        manifest_path=manifest_path,
        overlay_path=args.overlay.resolve() if args.overlay else None,
    )
    print(f"[review-tsv] {raw_tsv} → {review_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
