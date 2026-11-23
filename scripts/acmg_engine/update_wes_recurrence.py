#!/usr/bin/env python3
"""Add / refresh / remove WES cases in the ACMG engine's recurrence table.

Each WES case contributes the distinct variants of its 08 review TSV; a
re-run replaces the old contribution.  Registered cases merge by MRN, others
by sample ID without caller and ``_1``…``_N`` suffixes; ``…T`` test samples
are never added.  Safe to run repeatedly — unchanged cases are skipped.

    python3 scripts/acmg_engine/update_wes_recurrence.py --all            # scan every case
    python3 scripts/acmg_engine/update_wes_recurrence.py --sample 24WE0001-nckuh
    python3 scripts/acmg_engine/update_wes_recurrence.py --remove 24WE0001-nckuh
Add ``--dry-run`` to only report what would change.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))

from app import config  # noqa: E402
from app.services.acmg_engine.recurrence import WesRecurrence, load_roster, sync_cases  # noqa: E402


def _progress(index: int, total: int, sample_id: str) -> None:
    if index == 1 or index % 25 == 0 or index == total:
        print(f"[{index}/{total}] {sample_id}", file=sys.stderr, flush=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Maintain the WES recurrence table")
    ap.add_argument("--root", type=Path, default=Path(config.PIPELINE_OUT_ROOT))
    ap.add_argument("--db", type=Path, default=Path(config.ACMG_WES_RECURRENCE_DB))
    ap.add_argument("--sample", action="append", default=[])
    ap.add_argument("--remove", action="append", default=[])
    ap.add_argument("--all", action="store_true", help="scan every case and drop cases that no longer qualify")
    ap.add_argument("--force", action="store_true", help="rebuild even when the review TSV is unchanged")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--verbose", action="store_true", help="list skipped cases and their reasons")
    args = ap.parse_args(argv)
    if not (args.all or args.sample or args.remove):
        ap.error("give --all, --sample or --remove")
    args.db.parent.mkdir(parents=True, exist_ok=True)
    with open(str(args.db) + ".lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        store = WesRecurrence(args.db)
        try:
            removed = []
            for sid in args.remove:
                if args.dry_run or store.remove_sample(sid):
                    removed.append(sid)
            summary = {"removed": removed, "patients": store.patient_count()}
            roster = None
            if args.all or args.sample:
                roster = load_roster() or {}
                summary = sync_cases(store, args.root, None if args.all else args.sample,
                                     force=args.force, prune=args.all, dry_run=args.dry_run, roster=roster,
                                     progress=_progress)
                summary["removed"] = removed + summary["removed"]
        finally:
            store.close()
    skipped = summary.get("skipped", {})
    reasons: dict = {}
    for reason in skipped.values():
        reasons[reason] = reasons.get(reason, 0) + 1
    report = {
        "dry_run": args.dry_run,
        "added": len(summary.get("added", [])),
        "updated": len(summary.get("updated", [])),
        "unchanged": len(summary.get("unchanged", [])),
        "removed": summary.get("removed", []),
        "skipped": reasons,
        "patients": summary.get("patients"),
    }
    if roster is not None:
        report["roster_entries"] = len(roster)
    if args.verbose:
        report["skipped_cases"] = skipped
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
