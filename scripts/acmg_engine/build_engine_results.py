#!/usr/bin/env python3
"""Run the ACMG rule engine for one or more cases (plan §5.2, step 4).

Writes ``08_postprocessing/{LIS_ID}.snv_acmg_engine.sqlite`` from the case's
review TSV.  The tertiary post-processing calls the first form on the staged
case; the second form backfills cases that are already published::

    # one staged case (post-processing)
    build_engine_results.py --review R.tsv --out O.sqlite --sample SID --test-type WES

    # published cases (backfill); stale/missing only unless --force
    build_engine_results.py --all --dry-run
    build_engine_results.py --all
    build_engine_results.py --sample 24WE0001-nckuh --force

A result is stale when the review TSV, reference tables, engine code or
rulebook changed.  Changes in the WES recurrence table or OMIM inheritance do
not mark it stale; run with ``--force`` after a large recurrence update.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))

from app.services.acmg_engine import runner  # noqa: E402


def _load_json(path: Path) -> dict:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def _case_inputs(sample_id: str) -> dict:
    from app.services import sample_layout, test_types

    review = sample_layout.review_tsv(sample_id)
    manifest = _load_json(sample_layout.review_manifest(sample_id))
    meta = _load_json(sample_layout.state_file(sample_id, "sample_metadata.json"))
    source = _load_json(sample_layout.state_file(sample_id, "pipeline_source.json"))
    identity = str(meta.get("lis_id") or meta.get("sample_id") or sample_id)
    test_type = test_types.normalize_test_type(
        meta.get("test_type") or manifest.get("test_type") or source.get("seq_type") or "",
        sample_id=identity, default="WES")
    return {
        "review": review,
        "out": sample_layout.state_file(sample_id, runner.RESULTS_NAME, for_write=True),
        "existing": sample_layout.state_file(sample_id, runner.RESULTS_NAME),
        "test_type": "WGS" if test_types.is_wgs_type(test_type) else "WES",
        "mrn": str(meta.get("mrn") or "").strip(),
    }


def _count_reasons(skipped: dict) -> dict:
    counts: dict[str, int] = {}
    for reason in skipped.values():
        counts[reason] = counts.get(reason, 0) + 1
    return counts


def backfill(sample_ids, *, force: bool, dry_run: bool, refs, recurrence) -> dict:
    from app.workers.dragen_run import _acquire_sample_locks, _release_sample_locks

    report = {"built": [], "current": [], "skipped": {}, "failed": {}}
    for index, sid in enumerate(sample_ids, 1):
        if index == 1 or index % 25 == 0 or index == len(sample_ids):
            print(f"[{index}/{len(sample_ids)}] {sid}", file=sys.stderr, flush=True)
        case = _case_inputs(sid)
        if not case["review"].is_file():
            report["skipped"][sid] = "沒有 review TSV"
            continue
        if not force and runner.is_current(case["existing"], case["review"], refs):
            report["current"].append(sid)
            continue
        if dry_run:
            report["built"].append(sid)
            continue
        try:
            locks = _acquire_sample_locks([sid])
        except RuntimeError as exc:
            report["skipped"][sid] = str(exc)
            continue
        try:
            stats = runner.build_case(review_tsv=case["review"], out=case["out"], sample_id=sid,
                                      test_type=case["test_type"], refs=refs,
                                      recurrence=recurrence, mrn=case["mrn"])
            report["built"].append(sid)
            print(f"  {sid}: {stats['variants']} variants, {stats['seconds']}s", file=sys.stderr, flush=True)
        except Exception as exc:  # noqa: BLE001 — report and continue with the next case
            report["failed"][sid] = f"{type(exc).__name__}: {exc}"
        finally:
            _release_sample_locks(locks)
    return report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--review", type=Path, help="review TSV of one staged case")
    ap.add_argument("--out", type=Path, help="output sqlite for --review")
    ap.add_argument("--sample", action="append", default=[])
    ap.add_argument("--all", action="store_true", help="every published case")
    ap.add_argument("--test-type", default="", help="WES or WGS (with --review)")
    ap.add_argument("--mrn", default="", help="patient MRN (with --review; used to exclude the patient)")
    ap.add_argument("--refs-dir", type=Path, default=None)
    ap.add_argument("--recurrence-db", type=Path, default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    started = time.time()
    refs = runner.open_references(args.refs_dir)
    recurrence = runner.open_recurrence(args.recurrence_db)
    try:
        if args.review:
            if not args.out or len(args.sample) != 1:
                ap.error("--review needs --out and exactly one --sample")
            test_type = "WGS" if args.test_type.upper() in ("WGS", "TITAN-WGS") else "WES"
            stats = runner.build_case(review_tsv=args.review, out=args.out, sample_id=args.sample[0],
                                      test_type=test_type, refs=refs, recurrence=recurrence, mrn=args.mrn)
            print(f"[acmg-engine] sample={stats['sample_id']} variants={stats['variants']} "
                  f"seconds={stats['seconds']} tables={','.join(stats['reference_tables']) or 'none'} "
                  f"classes={json.dumps(stats['classes'], ensure_ascii=False)}", flush=True)
            return 0
        if not (args.all or args.sample):
            ap.error("give --review, --sample or --all")
        from app.services import sample_layout
        sample_ids = sorted(set(sample_layout.iter_sample_ids())) if args.all else args.sample
        report = backfill(sample_ids, force=args.force, dry_run=args.dry_run, refs=refs, recurrence=recurrence)
    finally:
        if recurrence is not None:
            recurrence.close()
    print(json.dumps({
        "dry_run": args.dry_run,
        "built": len(report["built"]),
        "current": len(report["current"]),
        "skipped": _count_reasons(report["skipped"]),
        "failed": report["failed"],
        "reference_tables": sorted(refs.available),
        "seconds": round(time.time() - started),
    }, ensure_ascii=False, indent=2))
    return 1 if report["failed"] else 0


if __name__ == "__main__":
    sys.exit(main())
