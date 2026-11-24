#!/usr/bin/env python3
"""Before/after report of the cleaned HPO → gene table on existing cases (phase 0).

Compares phenotype_to_genes.txt with a snapshot from build_hpo_gene_table.py:
table-level changes, each case's phenotype-related gene set, and the ranking of
reviewer-marked SNVs.  Read-only (does not rewrite pheno_score.tsv)::

    .venv/bin/python scripts/hpo/compare_pheno_tables.py --snapshot ~/NGS_UI/phenotype_data/hpo_kb/2026-09-01_2026-09-01
    .venv/bin/python scripts/hpo/compare_pheno_tables.py --snapshot DIR --out ~/NGS_UI/acmg_validation/pheno0 --sample 26WE0142-nckuh

report.md / summary.json / term_changes.tsv have no sample IDs;
case_changes.tsv and marked_variants.tsv do (keep them on the server).
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))

from app.config import NGS_UI_HOME, PHENO_DATA_DIR  # noqa: E402
from app.services import hpo_gene_table as H  # noqa: E402
from app.services import pheno_table_diff  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--snapshot", type=Path, default=None,
                    help="snapshot dir (default: the active one, else the newest under hpo_kb/)")
    ap.add_argument("--old-table", type=Path, default=PHENO_DATA_DIR / "phenotype_to_genes.txt")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--sample", action="append", default=None, help="limit to these sample IDs (repeatable)")
    args = ap.parse_args(argv)

    snap = args.snapshot
    if snap is None:
        active = H.active_table()
        if active.term_status:
            snap = active.gene_table.parent
        else:
            built = sorted(p for p in H.KB_ROOT.glob("*_*") if (p / H.GENE_TABLE).is_file())
            snap = built[-1] if built else None
    if snap is None or not (Path(snap) / H.GENE_TABLE).is_file():
        print("找不到快照；先跑 build_hpo_gene_table.py", file=sys.stderr)
        return 2
    if not args.old_table.is_file():
        print(f"找不到舊表：{args.old_table}", file=sys.stderr)
        return 2
    out = args.out or (NGS_UI_HOME / "acmg_validation" / f"pheno0_{time.strftime('%Y%m%d_%H%M%S')}")
    report = pheno_table_diff.run(out_dir=out, old_table=args.old_table, new_snapshot=Path(snap),
                                  sample_ids=args.sample, progress=lambda m: print(m, file=sys.stderr))
    print(f"完成（{report['elapsed_seconds']} 秒）：{out}/report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
