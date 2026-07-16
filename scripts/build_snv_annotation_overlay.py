#!/usr/bin/env python3
"""Build sparse post-processing annotations without duplicating the raw TSV."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app.services.snv_overlay import build_overlay  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", required=True, type=Path)
    parser.add_argument("--annotated", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    path = build_overlay(args.raw.resolve(), args.annotated.resolve(), args.out.resolve())
    print(f"[snv-overlay] {args.raw} + {args.annotated} → {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
