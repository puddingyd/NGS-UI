#!/usr/bin/env python3
"""Build the compact main-screen TSV from snv_indel.annotated.tsv."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app.services.snv_review import ensure_review_tsv  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv", required=True, help="complete snv_indel.annotated.tsv")
    args = ap.parse_args()

    raw_tsv = Path(args.tsv).resolve()
    if not raw_tsv.is_file():
        print(f"ERROR: --tsv 找不到：{raw_tsv}", file=sys.stderr)
        return 2
    review_tsv = ensure_review_tsv(raw_tsv)
    print(f"[review-tsv] {raw_tsv} → {review_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
