#!/usr/bin/env python3
"""Download/check the LitVar2 bulk export and atomically rebuild local SQLite."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app.services import litvar2_jobs, litvar2_store  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bulk-file",
        type=Path,
        help="Build from an already downloaded .json.gz instead of network download",
    )
    parser.add_argument(
        "--dataset-date",
        default="",
        help="YYYY-MM-DD version for --bulk-file (default: current UTC date)",
    )
    parser.add_argument(
        "--direct",
        action="store_true",
        help="Run directly without updater state/lock (intended for tests only)",
    )
    args = parser.parse_args()
    if args.bulk_file:
        result = litvar2_store.update_database(
            local_bulk=args.bulk_file,
            dataset_date=args.dataset_date,
            force=True,
        )
        print(result)
        return 0
    if args.direct:
        result = litvar2_store.update_database(force=True)
        print(result)
        return 0
    litvar2_jobs.run_update(trigger="scheduled")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
