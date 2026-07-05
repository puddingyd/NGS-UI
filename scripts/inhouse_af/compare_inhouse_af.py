#!/usr/bin/env python3
"""Phase C — validate the incremental in-house AF against the GLnexus baseline.

Both inputs are the genotype-stripped INHOUSE_* sites VCF (same contract):
  * incremental : publish_af.py output   (inhouse_af.hg38.vcf.gz)
  * glnexus     : build_inhouse_af.sh output (inhouse_af.glnexus.NNN.vcf.gz)

Both are `bcftools norm -m-` normalized, so a variant key (chrom,pos,ref,alt)
matches 1:1. We load the (smaller) GLnexus sites into a dict, stream the
incremental sites, and at every shared key accumulate the running sums needed
for a Pearson correlation of INHOUSE_AF — stratified by SNP/indel and by an AN
floor (so the low-AN tail, where GLnexus no-calls low-confidence samples, is
separated from the well-called sites the method is judged on).

Pure stdlib (gzip streaming). Prints a table; no plotting.

Usage:
  scripts/inhouse_af/compare_inhouse_af.py \
    --incremental $DB/inhouse_af.hg38.vcf.gz \
    --glnexus     $DB/inhouse_af.glnexus.675.vcf.gz \
    [--an-floors 0,1000,1200,1300]
"""
from __future__ import annotations

import argparse
import gzip
import math
import sys


def info_get(info: str, key: str):
    """Fetch a single INFO value (key=val;...). None if absent."""
    kv = key + "="
    for field in info.split(";"):
        if field.startswith(kv):
            return field[len(kv):]
    return None


def iter_sites(path, af_tag, an_tag):
    """Yield (chrom, pos, ref, alt, af, an) from an INHOUSE_* sites VCF."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            F = line.rstrip("\n").split("\t")
            if len(F) < 8:
                continue
            chrom, pos, ref, alt, info = F[0], F[1], F[3], F[4], F[7]
            af_s = info_get(info, af_tag)
            an_s = info_get(info, an_tag)
            if af_s is None or an_s is None:
                continue
            try:
                af = float(af_s.split(",")[0])
                an = int(an_s.split(",")[0])
            except ValueError:
                continue
            yield chrom, int(pos), ref, alt, af, an


class Accum:
    """Streaming Pearson accumulator."""
    __slots__ = ("n", "sx", "sy", "sxx", "syy", "sxy")

    def __init__(self):
        self.n = self.sx = self.sy = self.sxx = self.syy = self.sxy = 0.0

    def add(self, x, y):
        self.n += 1
        self.sx += x; self.sy += y
        self.sxx += x * x; self.syy += y * y
        self.sxy += x * y

    def pearson(self):
        n = self.n
        if n < 2:
            return float("nan")
        cov = self.sxy - self.sx * self.sy / n
        vx = self.sxx - self.sx * self.sx / n
        vy = self.syy - self.sy * self.sy / n
        if vx <= 0 or vy <= 0:
            return float("nan")
        return cov / math.sqrt(vx * vy)


def is_snp(ref, alt):
    return len(ref) == 1 and len(alt) == 1 and ref != "-" and alt != "-"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--incremental", required=True)
    ap.add_argument("--glnexus", required=True)
    ap.add_argument("--af-tag", default="INHOUSE_AF")
    ap.add_argument("--an-tag", default="INHOUSE_AN")
    ap.add_argument("--an-floors", default="0,1000,1200,1300",
                    help="comma-separated AN floors to stratify by (uses min AN of the two)")
    args = ap.parse_args()
    floors = sorted(int(x) for x in args.an_floors.split(",") if x.strip())

    # 1) load GLnexus baseline into a dict keyed by the normalized variant
    print(f"[compare] loading GLnexus sites: {args.glnexus}", file=sys.stderr)
    gln = {}
    for chrom, pos, ref, alt, af, an in iter_sites(args.glnexus, args.af_tag, args.an_tag):
        gln[(chrom, pos, ref, alt)] = (af, an)
    print(f"[compare]   {len(gln)} GLnexus sites", file=sys.stderr)

    # 2) stream incremental; accumulate at shared keys
    # strata: (variant_class, an_floor) -> Accum  (each site counts once per floor it clears)
    acc = {}
    for cls in ("snp", "indel"):
        for fl in floors:
            acc[(cls, fl)] = Accum()
    n_inc = n_shared = 0
    for chrom, pos, ref, alt, af_i, an_i in iter_sites(args.incremental, args.af_tag, args.an_tag):
        n_inc += 1
        g = gln.get((chrom, pos, ref, alt))
        if g is None:
            continue
        n_shared += 1
        af_g, an_g = g
        an_min = an_i if an_i < an_g else an_g
        cls = "snp" if is_snp(ref, alt) else "indel"
        for fl in floors:
            if an_min >= fl:
                acc[(cls, fl)].add(af_i, af_g)

    # 3) report
    print(f"\n[compare] incremental sites : {n_inc}")
    print(f"[compare] shared with GLnexus: {n_shared}  "
          f"({100.0*n_shared/max(1,n_inc):.1f}% of incremental)\n")
    print(f"{'class':6} {'AN>=':>6} {'n_sites':>12} {'pearson':>9}")
    print("-" * 38)
    for cls in ("snp", "indel"):
        for fl in floors:
            a = acc[(cls, fl)]
            print(f"{cls:6} {fl:>6} {int(a.n):>12} {a.pearson():>9.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
