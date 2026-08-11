#!/usr/bin/env python3
"""Backfill best-effort GPN-MSA annotation for existing tertiary cases."""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def _atomic_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    os.replace(tmp, path)


def _read_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _test_type(sample_layout, sample_id: str) -> str:
    metadata = _read_json(sample_layout.state_file(sample_id, "sample_metadata.json"))
    value = str(metadata.get("test_type") or "").upper()
    return value if value in {"WES", "WGS"} else "WES"


def _deduplicate(values: list[str]) -> list[str]:
    seen: set[str] = set()
    output: list[str] = []
    for value in values:
        sample_id = value.strip()
        if sample_id and sample_id not in seen:
            seen.add(sample_id)
            output.append(sample_id)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_ids", nargs="*", metavar="SAMPLE_ID")
    parser.add_argument(
        "--sample",
        action="append",
        default=[],
        dest="selected_samples",
        help="sample ID to backfill; may be repeated",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="process every case resolved by the tertiary layout service",
    )
    parser.add_argument(
        "--root",
        type=Path,
        help="override NGS_UI_TERTIARY_ROOT for this run",
    )
    parser.add_argument(
        "--gpn-msa-db",
        type=Path,
        help="override NGS_UI_GPN_MSA_DB for this run",
    )
    args = parser.parse_args()

    requested = _deduplicate([*args.selected_samples, *args.sample_ids])
    if args.all == bool(requested):
        parser.error("choose exactly one of --all or one/more sample IDs")

    if args.root:
        os.environ["NGS_UI_TERTIARY_ROOT"] = str(args.root.resolve())
    if args.gpn_msa_db:
        os.environ["NGS_UI_GPN_MSA_DB"] = str(args.gpn_msa_db.resolve())

    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root / "backend"))
    from app.services import gpn_msa, sample_layout, snv_review  # noqa: E402

    sample_ids = sorted(sample_layout.iter_sample_ids()) if args.all else requested
    if not sample_ids:
        print("[backfill-gpn-msa] no cases found", file=sys.stderr)
        return 1

    db = Path(gpn_msa.GPN_MSA_DB)
    counts: dict[str, int] = {}
    errors = 0

    for sample_id in sample_ids:
        review = sample_layout.review_tsv(sample_id)
        manifest = sample_layout.review_manifest(sample_id, for_write=True)
        try:
            if review.is_file():
                stats = gpn_msa.annotate_review_tsv(review, db)
                payload = _read_json(manifest)
                payload["gpn_msa"] = gpn_msa.database_signature(db)
                payload["gpn_msa_annotation"] = stats
                _atomic_json(manifest, payload)
            else:
                raw = sample_layout.snv_raw_tsv(sample_id)
                if not raw.is_file():
                    label = "ERROR" if not args.all else "WARNING"
                    print(
                        f"[backfill-gpn-msa] {label}: {sample_id}: "
                        "no review TSV or raw SNV TSV",
                        file=sys.stderr,
                    )
                    errors += not args.all
                    continue
                output_dir = sample_layout.state_dir(sample_id, for_write=True)
                review = snv_review.ensure_review_tsv(
                    raw,
                    test_type=_test_type(sample_layout, sample_id),
                    output_dir=output_dir,
                    output_path=sample_layout.review_tsv(sample_id, for_write=True),
                    manifest_path=manifest,
                    overlay_path=sample_layout.snv_overlay_path(sample_id),
                    gpn_msa_db=db,
                    force_gpn_msa=True,
                )
                payload = _read_json(manifest)
                stats = payload.get("gpn_msa_annotation") or {}

            status = str(stats.get("status") or "unknown")
            counts[status] = counts.get(status, 0) + 1
            print(
                f"[backfill-gpn-msa] {sample_id}: status={status} "
                f"annotated_rows={stats.get('annotated_rows', 0)} "
                f"review={review}"
            )
        except Exception as exc:
            errors += 1
            print(
                f"[backfill-gpn-msa] ERROR: {sample_id}: {exc}",
                file=sys.stderr,
            )

    summary = " ".join(f"{key}={value}" for key, value in sorted(counts.items()))
    print(
        f"[backfill-gpn-msa] done cases={sum(counts.values())} "
        f"errors={errors}{' ' + summary if summary else ''}"
    )
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
