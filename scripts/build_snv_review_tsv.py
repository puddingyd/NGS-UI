#!/usr/bin/env python3
"""Build the compact main-screen TSV from snv_indel.annotated.tsv."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app.services import test_types  # noqa: E402
from app.services.snv_review import ensure_review_tsv  # noqa: E402


def _infer_sample_id(raw_tsv: Path, output_dir: Path | None = None) -> str:
    directory = output_dir or raw_tsv.parent
    if directory.name == "08_postprocessing":
        return directory.parent.name
    if raw_tsv.name == "snv_indel.annotated.tsv":   # legacy UI tertiary_output/<SID>/
        return raw_tsv.parent.name
    return ""


def _read_metadata(raw_tsv: Path, output_dir: Path | None, sample_id: str) -> dict:
    directory = output_dir or raw_tsv.parent
    candidates = (
        directory / f"{sample_id}.sample_metadata.json",
        directory / "sample_metadata.json",
    ) if sample_id else (directory / "sample_metadata.json",)
    for meta_path in candidates:
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if isinstance(meta, dict):
            return meta
    return {}


def resolve_test_type(
    raw_tsv: Path,
    output_dir: Path | None,
    requested: str = "",
    sample_id: str = "",
) -> str:
    """Record the same test type the backend loader will ask for
    (sample_loader._effective_test_type): otherwise the review manifest never
    matches and the UI rebuilds the whole review TSV on first open. Year+T
    LIS IDs (26T...) are TITAN-WGS even when the worker says WGS; filtering is
    identical, only the manifest label differs."""
    sample_id = sample_id or _infer_sample_id(raw_tsv, output_dir)
    meta = _read_metadata(raw_tsv, output_dir, sample_id)
    identity = str(meta.get("lis_id") or meta.get("sample_id") or sample_id)
    value = requested or str(meta.get("test_type") or "")
    return test_types.normalize_test_type(value, sample_id=identity, default="WES")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv", required=True, help="complete snv_indel.annotated.tsv")
    ap.add_argument("--output-dir", type=Path, help="directory for derived review TSV")
    ap.add_argument("--output-path", type=Path, help="exact review TSV output path")
    ap.add_argument("--manifest-path", type=Path, help="exact review manifest output path")
    ap.add_argument("--overlay", type=Path, help="sparse SNV annotation overlay SQLite")
    ap.add_argument(
        "--gpn-msa-db",
        type=Path,
        help="fixed GRCh38 GPN-MSA scores.tsv.bgz (default: NGS_UI_GPN_MSA_DB)",
    )
    ap.add_argument(
        "--require-gpn-msa",
        action="store_true",
        help="fail if the GPN-MSA BGZF, .tbi, or tabix executable is unavailable",
    )
    ap.add_argument(
        "--test-type",
        type=str.upper,
        choices=["WES", "WGS", "TITAN-WGS"],
        help="Apply WES/WGS-specific review TSV filters. Defaults to sample_metadata.json, "
             "then WES. Year+T sample IDs are always recorded as TITAN-WGS.",
    )
    ap.add_argument(
        "--sample",
        default="",
        help="UI sample ID (default: inferred from --output-dir / --tsv path); "
             "used for the TITAN-WGS rule",
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
    test_type = resolve_test_type(raw_tsv, output_dir, args.test_type or "", args.sample)
    review_tsv = ensure_review_tsv(
        raw_tsv,
        test_type=test_type,
        output_dir=output_dir,
        output_path=output_path,
        manifest_path=manifest_path,
        overlay_path=args.overlay.resolve() if args.overlay else None,
        gpn_msa_db=args.gpn_msa_db.resolve() if args.gpn_msa_db else None,
        require_gpn_msa=args.require_gpn_msa,
    )
    print(f"[review-tsv] {raw_tsv} → {review_tsv}")
    status_path = manifest_path or review_tsv.with_suffix(
        review_tsv.suffix + ".source.json"
    )
    try:
        status_payload = json.loads(status_path.read_text(encoding="utf-8"))
        gpn_status = status_payload.get("gpn_msa_annotation") or {}
    except (OSError, json.JSONDecodeError):
        gpn_status = {}
    print(
        "[gpn-msa] "
        f"status={gpn_status.get('status') or 'unknown'} "
        f"annotated_rows={gpn_status.get('annotated_rows', 0)} "
        f"manifest={status_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
