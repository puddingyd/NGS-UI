"""Detached HPO update started by the home-page button."""
from __future__ import annotations

import argparse

from ..services import hpo_update


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--release", default="latest")
    args = parser.parse_args()
    return hpo_update.run_job(args.run_id, args.release)


if __name__ == "__main__":
    raise SystemExit(main())
