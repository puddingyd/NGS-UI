#!/usr/bin/env python3
"""Before/after report for replacing the pipeline ACMG with the rule engine (plan §5.2.2, step 6).

Simulates, for every published case, how the effective ACMG class, tier 1C
and ranking would change if the rule engine took the pipeline ACMG's place::

    .venv/bin/python scripts/acmg_engine/diff_replacement.py
    .venv/bin/python scripts/acmg_engine/diff_replacement.py --out ~/NGS_UI/acmg_validation/diff1

Read-only: no case file is changed.  Run with the ngs-ui service's Python (it
reads reviewer edits, manual ACMG and the active phenotype scores).  Cases
whose stored engine results are stale are skipped; run
build_engine_results.py --all first.

report.md / summary.json hold aggregate numbers only.  marked_variants.tsv
contains sample IDs and is for checking on the server only.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))

from app.services.acmg_engine import replacement_diff, runner  # noqa: E402


def main(argv=None) -> int:
    from app import config

    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=None,
                    help="output directory (default NGS_UI_HOME/acmg_validation/diff-<time>)")
    ap.add_argument("--sample", action="append", default=[], help="limit to these cases")
    ap.add_argument("--refs-dir", type=Path, default=None,
                    help="reference tables used to check that results are current (default: current)")
    args = ap.parse_args(argv)

    out = args.out or Path(config.NGS_UI_HOME) / "acmg_validation" / time.strftime("diff-%Y%m%d-%H%M%S")
    refs = runner.open_references(args.refs_dir)
    report = replacement_diff.run(
        out_dir=out, sample_ids=args.sample or None, refs=refs,
        progress=lambda msg: print(msg, file=sys.stderr, flush=True))
    print(f"cases={report['cases']} seconds={report['seconds']}")
    print(f"report: {out / 'report.md'}")
    print(f"marked variants (contains sample IDs, keep on server): {out / 'marked_variants.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
