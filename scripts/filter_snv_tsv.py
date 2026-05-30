#!/usr/bin/env python3
"""Pre-filter the unfiltered snv_indel.annotated.tsv from the new
tertiary pipeline (Phase 1 — VEP annotation only).

Keep AF, VAF, and IMPACT filtering in the UI so reviewers can relax
those display filters without rerunning the pipeline.

Always drop:
  * '*'-allele rows (no clinical meaning)
  * non-primary contigs unless --keep-alt-contigs is passed

Updates the input TSV in place unless --out-tsv given.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

# Primary assembly contigs (hg38). Anything else — alt haplotypes
# (_alt), unplaced (chrUn_*), random (_random), patches (_fix), decoy
# — has no clinical reporting value and is dropped by default.
_PRIMARY_CONTIGS = (
    {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM", "chrMT"}
    | {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}
)


def _is_primary_contig(chrom: str) -> bool:
    return chrom.strip() in _PRIMARY_CONTIGS


def filter_tsv(
    in_tsv: Path,
    out_tsv: Path,
    *,
    keep_alt_contigs: bool = False,
) -> dict:
    overwriting = in_tsv.resolve() == out_tsv.resolve()
    target = Path(str(out_tsv) + ".tmp") if overwriting else out_tsv
    target.parent.mkdir(parents=True, exist_ok=True)

    n_in = 0
    n_kept = 0
    n_drop_star = 0
    n_drop_contig = 0
    with open(in_tsv, "r", encoding="utf-8", newline="") as fi:
        reader = csv.DictReader(fi, delimiter="\t")
        fieldnames = reader.fieldnames or []
        with open(target, "w", encoding="utf-8", newline="") as fo:
            writer = csv.DictWriter(fo, fieldnames=fieldnames, delimiter="\t",
                                    extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                n_in += 1
                ref = (row.get("REF") or "").strip()
                alt = (row.get("ALT") or "").strip()
                if "*" in (ref, alt):
                    n_drop_star += 1
                    continue
                if not keep_alt_contigs and not _is_primary_contig(row.get("CHROM", "")):
                    n_drop_contig += 1
                    continue

                writer.writerow(row)
                n_kept += 1

    if overwriting:
        os.replace(target, out_tsv)

    return {
        "n_in":         n_in,
        "n_kept":       n_kept,
        "drop_star":    n_drop_star,
        "drop_contig":  n_drop_contig,
    }


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True,
                    help="snv_indel.annotated.tsv (filtered in place unless "
                         "--out-tsv given)")
    ap.add_argument("--out-tsv",
                    help="write filtered TSV here instead of overwriting --tsv")
    ap.add_argument("--keep-alt-contigs", action="store_true",
                    help="keep variants on alt haplotypes / chrUn / random / "
                         "decoy contigs (default: drop them — only "
                         "chr1-22 + chrX/Y/M are reported)")
    args = ap.parse_args()

    in_tsv = Path(args.tsv).resolve()
    if not in_tsv.is_file():
        print(f"ERROR: --tsv 找不到：{in_tsv}", file=sys.stderr)
        return 2
    out_tsv = Path(args.out_tsv).resolve() if args.out_tsv else in_tsv

    print(f"[filter] in  : {in_tsv}", file=sys.stderr)
    print(f"[filter] out : {out_tsv}", file=sys.stderr)
    print("[filter] rule: drop '*' alleles and non-primary contigs only",
          file=sys.stderr)

    stats = filter_tsv(in_tsv, out_tsv,
                       keep_alt_contigs=args.keep_alt_contigs)
    print(f"[filter] read    {stats['n_in']:>10} rows", file=sys.stderr)
    print(f"[filter] kept    {stats['n_kept']:>10} rows", file=sys.stderr)
    print(f"[filter] drop *  {stats['drop_star']:>10}", file=sys.stderr)
    print(f"[filter] drop ctg{stats['drop_contig']:>9}  "
          f"(alt/random/decoy contigs)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
