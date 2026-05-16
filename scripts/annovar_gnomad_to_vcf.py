#!/usr/bin/env python3
"""ANNOVAR gnomAD txt → sites-only VCF for `bcftools annotate`.

ANNOVAR's tab-separated format (`Chr Start End Ref Alt AF ...`) drops
the VCF padding base required for indels — '-' is used for the empty
allele. We round-trip back to VCF by looking up the missing base from
the reference FASTA:

  ANNOVAR (deletion)  : chr1 101 103 GT  -    AF=0.42
  VCF                 : chr1 100 .   AGT A    AF=0.42

  ANNOVAR (insertion) : chr1 100 100 -   GT   AF=0.001
  VCF                 : chr1 100 .   A   AGT  AF=0.001

  SNV / MNV: direct copy.

Output is sorted + bgzipped + tabix-indexed so `bcftools annotate -a`
can use it.

Usage:
    scripts/annovar_gnomad_to_vcf.py \\
        --txt $HOME/NGS_UI/biotools/hg38_gnomad41_genome.txt \\
        --ref /home/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta \\
        --out $HOME/NGS_UI/biotools/gnomad/gnomad_af.hg38.vcf.gz

Run inside the tertiary_python container so pysam + bcftools are
available; the companion wrapper `scripts/build_gnomad_af_vcf.sh`
handles that.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

try:
    import pysam
except ImportError:
    print("ERROR: pysam not available. Run inside the tertiary_python "
          "container (or pip install pysam).", file=sys.stderr)
    sys.exit(2)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--txt", required=True, help="ANNOVAR gnomAD txt")
    ap.add_argument("--ref", required=True,
                    help="Reference FASTA (used to look up indel padding bases)")
    ap.add_argument("--out", required=True,
                    help="Output sites VCF.gz path (sorted + indexed)")
    ap.add_argument("--af-col", default="AF",
                    help="AF column name in the ANNOVAR header (default 'AF')")
    ap.add_argument("--min-af", type=float, default=None,
                    help="Skip rows whose AF < this (default: keep all)")
    ap.add_argument("--snv-only", action="store_true",
                    help="Skip indels — much faster but indels won't carry "
                         "gnomAD AF downstream")
    args = ap.parse_args()

    txt_path = Path(args.txt)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_vcf = out_path.with_suffix(out_path.suffix + ".unsorted")

    fasta = pysam.FastaFile(args.ref)

    n_in = n_out = n_skip_af = n_skip_indel = n_skip_invalid = 0
    with open(txt_path, "r", encoding="utf-8") as fi, \
         open(tmp_vcf, "w", encoding="utf-8") as fo:

        header_line = fi.readline().lstrip("#").rstrip("\n")
        header_cols = header_line.split("\t")
        try:
            af_idx = header_cols.index(args.af_col)
        except ValueError:
            print(f"ERROR: --af-col '{args.af_col}' not in header: "
                  f"{header_cols[:10]}…", file=sys.stderr)
            return 2
        print(f"[annovar→vcf] AF column '{args.af_col}' at index {af_idx} "
              f"(of {len(header_cols)} cols)", file=sys.stderr)

        fo.write("##fileformat=VCFv4.2\n")
        fo.write('##INFO=<ID=gnomAD_AF,Number=1,Type=Float,'
                 'Description="gnomAD v4.1 allele frequency (from ANNOVAR txt)">\n')
        fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for line in fi:
            n_in += 1
            cols = line.rstrip("\n").split("\t")
            if len(cols) < max(5, af_idx + 1):
                n_skip_invalid += 1
                continue
            chrom = cols[0].strip()
            start_s = cols[1].strip()
            ref = cols[3].strip()
            alt = cols[4].strip()
            af = cols[af_idx].strip()
            try:
                start = int(start_s)
                af_f = float(af)
            except ValueError:
                n_skip_invalid += 1
                continue
            if args.min_af is not None and af_f < args.min_af:
                n_skip_af += 1
                continue

            is_indel = (ref == "-" or alt == "-" or len(ref) != len(alt))
            if args.snv_only and is_indel:
                n_skip_indel += 1
                continue

            try:
                if ref != "-" and alt != "-" and len(ref) == len(alt):
                    # SNV / MNV — copy through
                    vcf_pos, vcf_ref, vcf_alt = start, ref, alt
                elif alt == "-":
                    # Deletion: VCF needs the base BEFORE Start as padding
                    pad = fasta.fetch(chrom, start - 2, start - 1).upper()
                    if len(pad) != 1:
                        n_skip_invalid += 1
                        continue
                    vcf_pos = start - 1
                    vcf_ref = pad + ref
                    vcf_alt = pad
                elif ref == "-":
                    # Insertion: ANNOVAR Start is the base BEFORE the insertion
                    pad = fasta.fetch(chrom, start - 1, start).upper()
                    if len(pad) != 1:
                        n_skip_invalid += 1
                        continue
                    vcf_pos = start
                    vcf_ref = pad
                    vcf_alt = pad + alt
                else:
                    # Block substitution (uneven length, no '-')
                    vcf_pos, vcf_ref, vcf_alt = start, ref, alt
            except (KeyError, ValueError) as e:
                n_skip_invalid += 1
                continue

            fo.write(f"{chrom}\t{vcf_pos}\t.\t{vcf_ref}\t{vcf_alt}\t.\t.\t"
                     f"gnomAD_AF={af}\n")
            n_out += 1
            if n_in % 1_000_000 == 0:
                print(f"[annovar→vcf] processed {n_in:>11,}  wrote {n_out:>11,}",
                      file=sys.stderr)

    print(f"[annovar→vcf] read   {n_in:>11,} rows", file=sys.stderr)
    print(f"[annovar→vcf] wrote  {n_out:>11,} VCF lines", file=sys.stderr)
    print(f"[annovar→vcf] skip AF{n_skip_af:>11,}", file=sys.stderr)
    print(f"[annovar→vcf] skip indel{n_skip_indel:>8,}", file=sys.stderr)
    print(f"[annovar→vcf] skip invalid{n_skip_invalid:>6,}", file=sys.stderr)

    # Sort + bgzip + tabix via bcftools.
    print(f"[annovar→vcf] sorting + bgzipping → {out_path}", file=sys.stderr)
    subprocess.run(["bcftools", "sort", "--max-mem", "4G",
                    str(tmp_vcf), "-Oz", "-o", str(out_path)], check=True)
    subprocess.run(["bcftools", "index", "-t", "-f", str(out_path)], check=True)
    tmp_vcf.unlink()
    print(f"[annovar→vcf] done → {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
