#!/usr/bin/env python3
"""Pre-filter the unfiltered snv_indel.annotated.tsv from the new
tertiary pipeline (Phase 1 — VEP annotation only) down to a
clinically-tractable set.

Phase 1 dumps every VEP-annotated variant (~5 M rows for a typical
WGS), including upstream_gene_variant / intergenic noise. The GUI and
GeneBe step both choke on that. This script applies the same broad
rare-variant + protein-impact filter the old R pipeline used to do
upstream of TSV creation.

Default rule — KEEP if EITHER:
  (a) max(--af-cols) <= --max-af  AND  IMPACT in --impact
  (b) CLINVAR_SIG matches /pathogenic|likely_pathogenic/i
DROP everything else. Always drop '*'-allele rows (no clinical meaning).

Defaults:
    --max-af 0.01
    --af-cols GNOMAD_G_AF
    --impact  HIGH,MODERATE

Updates the input TSV in place unless --out-tsv given.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path

PATHO_RE = re.compile(r"\b(?:Likely_)?[Pp]athogenic\b")

# Primary assembly contigs (hg38). Anything else — alt haplotypes
# (_alt), unplaced (chrUn_*), random (_random), patches (_fix), decoy
# — has no clinical reporting value and is dropped by default.
_PRIMARY_CONTIGS = (
    {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM", "chrMT"}
    | {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}
)


def _is_primary_contig(chrom: str) -> bool:
    return chrom.strip() in _PRIMARY_CONTIGS


def _max_af(row: dict, cols: list[str]) -> float:
    m = 0.0
    for c in cols:
        v = row.get(c)
        if v is None:
            continue
        s = str(v).strip()
        if not s or s.upper() in ("NA", "N/A", "."):
            continue
        try:
            f = float(s)
        except ValueError:
            continue
        if f > m:
            m = f
    return m


def _is_clinvar_patho(row: dict) -> bool:
    sig = (row.get("CLINVAR_SIG") or "").strip()
    if not sig or sig in ("NA", "."):
        return False
    return bool(PATHO_RE.search(sig))


def filter_tsv(
    in_tsv: Path,
    out_tsv: Path,
    *,
    max_af: float,
    af_cols: list[str],
    impact_keep: set[str],
    keep_alt_contigs: bool = False,
) -> dict:
    overwriting = in_tsv.resolve() == out_tsv.resolve()
    target = Path(str(out_tsv) + ".tmp") if overwriting else out_tsv
    target.parent.mkdir(parents=True, exist_ok=True)

    n_in = 0
    n_kept = 0
    n_drop_star = 0
    n_drop_contig = 0
    n_drop_af = 0
    n_drop_impact = 0
    n_keep_clinvar_rescue = 0  # would have been dropped, kept by ClinVar rule

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

                clinvar_patho = _is_clinvar_patho(row)
                af = _max_af(row, af_cols)
                impact = (row.get("IMPACT") or "").strip().upper()

                if clinvar_patho:
                    pass  # always keep
                else:
                    if af > max_af:
                        n_drop_af += 1
                        continue
                    if impact_keep and impact not in impact_keep:
                        n_drop_impact += 1
                        continue
                if clinvar_patho and (af > max_af or
                                       (impact_keep and impact not in impact_keep)):
                    n_keep_clinvar_rescue += 1
                writer.writerow(row)
                n_kept += 1

    if overwriting:
        os.replace(target, out_tsv)

    return {
        "n_in":         n_in,
        "n_kept":       n_kept,
        "drop_star":    n_drop_star,
        "drop_contig":  n_drop_contig,
        "drop_af":      n_drop_af,
        "drop_impact":  n_drop_impact,
        "clinvar_rescue": n_keep_clinvar_rescue,
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
    ap.add_argument("--max-af", type=float, default=0.01,
                    help="drop rows whose max(--af-cols) > this (default 0.01; "
                         "use a large number e.g. 1.1 to disable)")
    ap.add_argument("--af-cols", default="GNOMAD_G_AF",
                    help="comma-separated AF columns to check "
                         "(default GNOMAD_G_AF — gnomAD genome 'global' AF)")
    ap.add_argument("--impact", default="HIGH,MODERATE",
                    help="comma-separated IMPACT values to keep "
                         "(default HIGH,MODERATE; pass empty string to disable "
                         "the impact filter)")
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

    af_cols = [c.strip() for c in args.af_cols.split(",") if c.strip()]
    impact_keep = {s.strip().upper() for s in args.impact.split(",") if s.strip()}

    print(f"[filter] in  : {in_tsv}", file=sys.stderr)
    print(f"[filter] out : {out_tsv}", file=sys.stderr)
    print(f"[filter] rule: max({','.join(af_cols)}) <= {args.max_af}  "
          f"AND IMPACT in {sorted(impact_keep) if impact_keep else '(any)'}  "
          f"OR CLINVAR_SIG matches /pathogenic|likely_pathogenic/i",
          file=sys.stderr)

    stats = filter_tsv(in_tsv, out_tsv,
                       max_af=args.max_af,
                       af_cols=af_cols,
                       impact_keep=impact_keep,
                       keep_alt_contigs=args.keep_alt_contigs)
    print(f"[filter] read    {stats['n_in']:>10} rows", file=sys.stderr)
    print(f"[filter] kept    {stats['n_kept']:>10} rows  "
          f"(ClinVar P/LP rescue: {stats['clinvar_rescue']})",
          file=sys.stderr)
    print(f"[filter] drop *  {stats['drop_star']:>10}", file=sys.stderr)
    print(f"[filter] drop ctg{stats['drop_contig']:>9}  "
          f"(alt/random/decoy contigs)", file=sys.stderr)
    print(f"[filter] drop AF {stats['drop_af']:>10}", file=sys.stderr)
    print(f"[filter] drop IMP{stats['drop_impact']:>10}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
