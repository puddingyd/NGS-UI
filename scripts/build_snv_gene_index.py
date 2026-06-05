#!/usr/bin/env python3
"""Build snv_gene_index.sqlite for fast full-TSV gene search."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app.services import snv_gene_index  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv", required=True, help="complete snv_indel.annotated.tsv")
    args = ap.parse_args()

    raw_tsv = Path(args.tsv).resolve()
    if not raw_tsv.is_file():
        print(f"ERROR: --tsv 找不到：{raw_tsv}", file=sys.stderr)
        return 2
    index_path = snv_gene_index.build_index(raw_tsv)
    print(f"[gene-index] {raw_tsv} → {index_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
