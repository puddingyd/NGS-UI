"""Foreground LitVar2 updater used by systemd and the manual UI launcher."""
from __future__ import annotations

import argparse

from ..services import litvar2_jobs


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trigger", choices=["manual", "scheduled"], default="scheduled")
    parser.add_argument("--run-id", default="")
    args = parser.parse_args()
    litvar2_jobs.run_update(trigger=args.trigger, run_id=args.run_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
