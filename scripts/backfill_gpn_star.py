#!/usr/bin/env python3
"""Backfill GPN-Star on existing compact reviews, without rerunning Nextflow.

Requires --research-only. Default is dry-run. --apply takes the tertiary
sample lock. Existing review rows and all other manifest signatures are preserved (including stale ones).
Cases without a review TSV are skipped; rebuild their post-processing first.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path


def publish_pair(review: Path, manifest: Path, staged_review: Path, payload: dict) -> None:
    """Publish manifest last and restore the review if publication fails."""
    backup = staged_review.with_suffix(".backup")
    shutil.copy2(review, backup)
    fd, name = tempfile.mkstemp(prefix=f".{manifest.name}.", dir=manifest.parent)
    staged_manifest = Path(name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
        os.replace(staged_review, review)
        try:
            os.replace(staged_manifest, manifest)
        except BaseException:
            os.replace(backup, review)
            raise
    finally:
        staged_manifest.unlink(missing_ok=True)
        backup.unlink(missing_ok=True)


def backfill_one(sample: str, root: Path, *, apply: bool, research_only: bool = False) -> dict:
    from app.services import gpn_star, sample_layout

    if not research_only:
        raise ValueError("GPN-Star backfill requires --research-only")

    # Resolve exact existing paths: legacy v2 files remain valid and must not
    # be shadowed by an uncoordinated prefixed copy during a score-only refresh.
    review = sample_layout.review_tsv(sample)
    manifest = sample_layout.review_manifest(sample)
    if not review.is_file():
        return {"status": "skipped_missing_review", "sample": sample}
    if not apply:
        return {"status": "dry_run", "sample": sample, "review": str(review), "database": str(root)}
    from app.workers.dragen_run import _acquire_sample_locks, _release_sample_locks

    locks = _acquire_sample_locks([sample])
    try:
        review = sample_layout.review_tsv(sample)
        manifest = sample_layout.review_manifest(sample)
        try:
            payload = json.loads(manifest.read_text(encoding="utf-8"))
        except FileNotFoundError:
            payload = {}
        if not isinstance(payload, dict):
            raise ValueError("review manifest must be an object")
        with tempfile.TemporaryDirectory(prefix=".gpn-star-", dir=review.parent) as folder:
            staged = Path(folder) / review.name
            shutil.copy2(review, staged)
            results = {}
            payload["gpn_star_research_only"] = True
            for name, model in (("gpn_star", gpn_star.MODEL), ("gpn_star_v", gpn_star.V_MODEL)):
                stats = gpn_star.annotate_review_tsv(staged, root, model=model)
                signature = gpn_star.database_signature(root, model=model)
                payload[name] = signature
                payload[f"{name}_annotation"] = stats
                if stats["status"] == "complete":
                    payload[f"{name}_values"] = {**signature, "revision": stats["revision"]}
                results[name] = stats
            publish_pair(review, manifest, staged, payload)
        # Running API processes observe the changed review mtime/size in their
        # SNV cache keys. No unrelated reviewer state or raw index is touched.
        return {"sample": sample, "models": results,
                "status": "complete" if all(s["status"] == "complete" for s in results.values()) else "incomplete"}
    finally:
        _release_sample_locks(locks)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_ids", nargs="*")
    parser.add_argument("--sample", action="append", default=[])
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--root", type=Path)
    parser.add_argument("--gpn-star-dir", type=Path)
    parser.add_argument("--research-only", action="store_true", help="explicitly opt these cases into research GPN-Star")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--apply", action="store_true")
    mode.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if not args.research_only:
        parser.error("GPN-Star backfill requires --research-only")
    requested = list(dict.fromkeys([*args.sample, *args.sample_ids]))
    if args.all == bool(requested):
        parser.error("choose --all or one/more sample IDs")
    if args.root:
        os.environ["NGS_UI_TERTIARY_ROOT"] = str(args.root.resolve())
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))
    from app.services import gpn_star, sample_layout

    available = sorted(set(sample_layout.iter_sample_ids()))
    sample_ids = available if args.all else []
    errors = 0
    for value in requested:
        matches = [value] if value in available else [
            sample for sample in available if sample_layout.source_sample_id(sample) == value
            or sample in {f"{value}-nckuh", f"{value}-dragen", f"{value}-inhouse"}]
        if not matches:
            print(json.dumps({"sample": value, "status": "not_found"}))
            errors += 1
        sample_ids.extend(matches)
    root = args.gpn_star_dir or gpn_star.GPN_STAR_DIR
    for sample in dict.fromkeys(sample_ids):
        try:
            result = backfill_one(sample, root, apply=args.apply, research_only=args.research_only)
        except Exception as exc:
            result = {"sample": sample, "status": "failed", "message": str(exc)}
        print(json.dumps(result, ensure_ascii=False), flush=True)
        errors += result["status"] not in {"complete", "dry_run"}
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
