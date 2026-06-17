#!/usr/bin/env python3
"""Phase B (2/2) — publish the in-house AF sites VCF from the accumulated DB.

Joins counts.sqlite (numerator) with an_track.bg.gz (denominator) and emits the
genotype-stripped sites VCF the UI consumes, identical in contract to the
GLnexus path:

  INFO/INHOUSE_AC    = 2·n_hom + n_het + n_hemi   (chrM: n_mt_hom + n_mt_het)
  INFO/INHOUSE_AN    = AN-track value at the site
  INFO/INHOUSE_AF    = AC / AN
  INFO/INHOUSE_NHOM  = n_hom                       (chrM: n_mt_hom, homoplasmic)
  INFO/INHOUSE_HEMI  = n_hemi
  INFO/INHOUSE_HET_MT= n_mt_het                    (chrM heteroplasmic carriers)

The variant rows (sqlite, ORDER BY chrom,pos) and the AN track (bedGraph, sorted
the same way) are merge-joined in a single streaming pass — no tabix lookups, no
big memory, no extra deps. Sites with AN=0 or AC=0 are dropped.

Usage:
  scripts/inhouse_af/publish_af.py --db-dir $NGS_UI_HOME/biotools/inhouse_af \
    --ref /home/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta
"""
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import sqlite3
import subprocess
import sys

INFO_HEADER = """##INFO=<ID=INHOUSE_AC,Number=A,Type=Integer,Description="In-house alt allele count">
##INFO=<ID=INHOUSE_AN,Number=1,Type=Integer,Description="In-house total called alleles (coverage AN track, DP>=10, ploidy-weighted)">
##INFO=<ID=INHOUSE_AF,Number=A,Type=Float,Description="In-house allele frequency = AC/AN">
##INFO=<ID=INHOUSE_NHOM,Number=A,Type=Integer,Description="In-house homozygous-alt individuals (chrM: homoplasmic carriers)">
##INFO=<ID=INHOUSE_HEMI,Number=A,Type=Integer,Description="In-house hemizygous-alt individuals (male non-PAR X / Y)">
##INFO=<ID=INHOUSE_HET_MT,Number=A,Type=Integer,Description="In-house chrM heteroplasmic carriers">"""


def an_track_iter(path):
    """Yield (chrom, start, end, an) from the bedGraph in sorted order."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            c, s, e, v = line.rstrip("\n").split("\t")
            yield c, int(s), int(e), int(v)


def fai_contigs(ref):
    out = []
    fai = ref + ".fai"
    if os.path.exists(fai):
        with open(fai) as f:
            for line in f:
                p = line.split("\t")
                out.append((p[0], p[1]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db-dir", required=True)
    ap.add_argument("--ref", default="/home/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta")
    ap.add_argument("--out", help="default <db-dir>/inhouse_af.hg38.vcf.gz")
    args = ap.parse_args()

    db_path = os.path.join(args.db_dir, "counts.sqlite")
    an_path = os.path.join(args.db_dir, "an_track.bg.gz")
    if not os.path.exists(an_path):
        an_path = os.path.join(args.db_dir, "an_track.bg")
    out = args.out or os.path.join(args.db_dir, "inhouse_af.hg38.vcf.gz")
    for p in (db_path, an_path):
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing {p} (run accumulate.py first)")

    bgzip = shutil.which("bgzip")
    tmp_vcf = out + (".tmp.vcf" if bgzip else "")
    final_plain = out[:-3] if (out.endswith(".gz") and not bgzip) else out
    write_path = tmp_vcf if bgzip else final_plain

    conn = sqlite3.connect(db_path)
    nsamp = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    cur = conn.execute(
        "SELECT chrom,pos,ref,alt,n_hom,n_het,n_hemi,n_mt_hom,n_mt_het "
        "FROM variant_counts ORDER BY chrom,pos,ref,alt")

    an = an_track_iter(an_path)
    cur_iv = next(an, None)
    n_written = n_drop = 0

    with open(write_path, "w") as o:
        o.write("##fileformat=VCFv4.2\n")
        o.write(f"##source=inhouse_af_incremental;samples={nsamp}\n")
        o.write(INFO_HEADER + "\n")
        for c, ln in fai_contigs(args.ref):
            o.write(f"##contig=<ID={c},length={ln}>\n")
        o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for chrom, pos, ref, alt, n_hom, n_het, n_hemi, n_mt_hom, n_mt_het in cur:
            # advance AN track to the interval covering (chrom,pos)
            while cur_iv and (cur_iv[0] < chrom or
                              (cur_iv[0] == chrom and cur_iv[2] <= pos)):
                cur_iv = next(an, None)
            an_val = (cur_iv[3] if (cur_iv and cur_iv[0] == chrom
                                    and cur_iv[1] <= pos < cur_iv[2]) else 0)
            if chrom == "chrM":
                ac = n_mt_hom + n_mt_het
                nhom, hemi, het_mt = n_mt_hom, 0, n_mt_het
            else:
                ac = 2 * n_hom + n_het + n_hemi
                nhom, hemi, het_mt = n_hom, n_hemi, 0
            if an_val <= 0 or ac <= 0:
                n_drop += 1
                continue
            af = ac / an_val
            info = (f"INHOUSE_AC={ac};INHOUSE_AN={an_val};INHOUSE_AF={af:.6g};"
                    f"INHOUSE_NHOM={nhom};INHOUSE_HEMI={hemi};INHOUSE_HET_MT={het_mt}")
            o.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\n")
            n_written += 1
    conn.close()

    if bgzip:
        with open(out, "wb") as fo:
            subprocess.run([bgzip, "-c", write_path], stdout=fo, check=True)
        os.unlink(write_path)
        if shutil.which("tabix"):
            subprocess.run(["tabix", "-f", "-p", "vcf", out], check=False)
        final = out
    else:
        final = final_plain
        print("[publish] WARN: bgzip not found — wrote plain VCF (no .gz/.tbi)", file=sys.stderr)

    print(f"[publish] {n_written} sites written, {n_drop} dropped (AN=0/AC=0) -> {final}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
