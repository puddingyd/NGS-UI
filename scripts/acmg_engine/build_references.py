#!/usr/bin/env python3
"""Build ``acmg_refs.sqlite`` for the ACMG rule engine (plan §5.2.11).

Inputs (defaults in brackets):
  --clinvar       ClinVar variant_summary  [<refs>/raw/variant_summary.txt.gz]
  --mane-gtf      MANE Ensembl GTF         [newest <refs>/raw/MANE.GRCh38.v*.ensembl_genomic.gtf.gz]
  --mane-summary  MANE summary             [$NGS_UI_MANE_SUMMARY or biotools/MANE.GRCh38.v1.5.summary.txt.gz]
  --clingen-hi    ClinGen dosage TSV       [pipeline /home/pipeline/reference/hg38/tertiary/clingen/...]
Output: <--out>/acmg_refs.sqlite (atomic replace) [NGS_UI_ACMG_REFS_DIR].

For time-split validation build a separate directory from an archived
variant_summary and pass it to the engine instead of the live one.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from email.utils import parsedate_to_datetime
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))

from app import config  # noqa: E402
from app.services.acmg_engine import ref_store  # noqa: E402

DEFAULT_CLINGEN = Path("/home/pipeline/reference/hg38/tertiary/clingen/ClinGen_gene_curation_list_GRCh38.tsv")


def _release(path: Path) -> str:
    sidecar = path.with_name(path.name + ".json")
    try:
        modified = json.loads(sidecar.read_text(encoding="utf-8")).get("last_modified") or ""
        if modified:
            return parsedate_to_datetime(modified).date().isoformat()
    except (OSError, ValueError, TypeError):
        pass
    return time.strftime("%Y-%m-%d", time.localtime(path.stat().st_mtime))


def main(argv=None) -> int:
    refs_dir = Path(config.ACMG_REFS_DIR)
    raw = refs_dir / "raw"
    gtfs = sorted(raw.glob("MANE.GRCh38.v*.ensembl_genomic.gtf.gz"))
    mane_summary = Path(os.environ.get("NGS_UI_MANE_SUMMARY")
                        or Path(config.BIOTOOLS_DIR) / "MANE.GRCh38.v1.5.summary.txt.gz")
    ap = argparse.ArgumentParser(description="Build the ACMG rule engine reference database")
    ap.add_argument("--clinvar", type=Path, default=raw / "variant_summary.txt.gz")
    ap.add_argument("--mane-gtf", type=Path, default=gtfs[-1] if gtfs else None)
    ap.add_argument("--mane-summary", type=Path, default=mane_summary)
    ap.add_argument("--clingen-hi", type=Path, default=DEFAULT_CLINGEN)
    ap.add_argument("--out", type=Path, default=refs_dir)
    ap.add_argument("--clinvar-release", default="", help="default: Last-Modified date of the download")
    ap.add_argument("--min-truncating", type=int, default=2)
    args = ap.parse_args(argv)
    for label, path in (("ClinVar variant_summary", args.clinvar), ("MANE summary", args.mane_summary),
                        ("ClinGen HI", args.clingen_hi)):
        if not path or not Path(path).is_file():
            raise SystemExit(f"找不到 {label}：{path}")
    if not args.mane_gtf or not Path(args.mane_gtf).is_file():
        print("警告：沒有 MANE GTF，剪接位 PVS1 會保守給 Strong", file=sys.stderr)
        args.mane_gtf = None
    started = time.time()
    counts = ref_store.build_reference_db(
        args.out, clinvar_summary=args.clinvar, mane_summary=args.mane_summary,
        clingen_hi=args.clingen_hi, mane_gtf=args.mane_gtf, min_truncating=args.min_truncating,
        clinvar_release=args.clinvar_release or _release(args.clinvar))
    print(json.dumps({"out": str(args.out / ref_store.REFS_DB_NAME), "counts": counts,
                      "seconds": round(time.time() - started)}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
