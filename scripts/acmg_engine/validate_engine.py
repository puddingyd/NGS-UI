#!/usr/bin/env python3
"""Offline validation of the ACMG rule engine (plan §5.2.10, step 5).

Reads every published case's review TSV and stored engine results and writes
an aggregate report (no sample IDs, patient keys or reviewer names)::

    .venv/bin/python scripts/acmg_engine/validate_engine.py
    .venv/bin/python scripts/acmg_engine/validate_engine.py --out ~/acmg_validation/run1

ClinVar time-split: build reference tables from an older ClinVar
``variant_summary`` into a separate directory, then pass both::

    build_references.py --clinvar OLD.txt.gz --out ~/acmg_refs_t0 --clinvar-release 2024-01
    validate_engine.py --refs-dir ~/acmg_refs_t0 --t0-clinvar OLD.txt.gz

Read-only: no case file is changed.  Run with the ngs-ui service's Python
(it reads OMIM inheritance and the active phenotype scores).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))

from app.services.acmg_engine import runner, validation  # noqa: E402


def main(argv=None) -> int:
    from app import config

    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=None,
                    help="output directory (default NGS_UI_HOME/acmg_validation/<time>)")
    ap.add_argument("--sample", action="append", default=[], help="limit to these cases")
    ap.add_argument("--refs-dir", type=Path, default=None, help="reference tables (default: current)")
    ap.add_argument("--t0-clinvar", type=Path, default=None, help="older variant_summary for the time-split")
    ap.add_argument("--recurrence-db", type=Path, default=Path(config.ACMG_WES_RECURRENCE_DB))
    ap.add_argument("--manual-db", type=Path, default=Path(config.MANUAL_ACMG_DB))
    ap.add_argument("--allow-no-omim", action="store_true")
    args = ap.parse_args(argv)

    if runner.omim_status() != "loaded" and not args.allow_no_omim:
        print(f"讀不到 OMIM 遺傳模式（{sys.executable}）；請用 ngs-ui 服務的 Python 執行，"
              "或加 --allow-no-omim。", file=sys.stderr)
        return 2
    out = args.out or Path(config.NGS_UI_HOME) / "acmg_validation" / time.strftime("%Y%m%d-%H%M%S")
    refs = runner.open_references(args.refs_dir)
    recurrence = runner.open_recurrence(args.recurrence_db)
    try:
        report = validation.run(
            out_dir=out, sample_ids=args.sample or None, refs=refs, recurrence=recurrence,
            recurrence_db=args.recurrence_db, manual_db=args.manual_db, t0_clinvar=args.t0_clinvar,
            progress=lambda msg: print(msg, file=sys.stderr, flush=True))
    finally:
        if recurrence is not None:
            recurrence.close()
    print(json.dumps({"out": str(out), "cases": report["cases"], "seconds": report["seconds"],
                      "time_split": report["time_split"]}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
