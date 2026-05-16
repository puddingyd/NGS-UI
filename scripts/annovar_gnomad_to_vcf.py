#!/usr/bin/env python3
"""ANNOVAR gnomAD txt → sites-only VCF for `bcftools annotate`.

ANNOVAR's tab-separated format (`Chr Start End Ref Alt AF ...`) drops
the VCF padding base required for indels — '-' is used for the empty
allele. We round-trip back to VCF by looking up the missing base from
the reference FASTA via its `.fai` index (stdlib only — no pysam):

  ANNOVAR (deletion)  : chr1 101 103 GT  -    AF=0.42
  VCF                 : chr1 100 .   AGT A    AF=0.42

  ANNOVAR (insertion) : chr1 100 100 -   GT   AF=0.001
  VCF                 : chr1 100 .   A   AGT  AF=0.001

  SNV / MNV: direct copy.

Output is an UNSORTED VCF (sort + bgzip + tabix is the caller's job —
see scripts/build_gnomad_af_vcf.sh which pipes the output through
`bcftools sort` in a container).

Usage:
    scripts/annovar_gnomad_to_vcf.py \\
        --txt $HOME/NGS_UI/biotools/hg38_gnomad41_genome.txt \\
        --ref /home/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta \\
        --out gnomad_af.hg38.unsorted.vcf \\
        --af-col gnomad41_genome_AF \\
        --min-af 0.01
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


class FastaRandom:
    """Stdlib-only random-access FASTA reader using the `.fai` index.

    Equivalent to pysam.FastaFile.fetch(chrom, start, end) for our
    purposes (0-based half-open). Memory: just the index dict
    (~tens of KB for hg38). Per-fetch cost: one seek + one read.
    """

    def __init__(self, fasta_path: str | Path):
        self.fasta_path = str(fasta_path)
        self.fh = open(self.fasta_path, "rb")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        fai_path = self.fasta_path + ".fai"
        with open(fai_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                chrom, length, offset, linebases, linewidth = parts[:5]
                self.index[chrom] = (int(length), int(offset),
                                      int(linebases), int(linewidth))

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """0-based half-open. Returns "" on out-of-range / missing chrom."""
        meta = self.index.get(chrom)
        if meta is None:
            return ""
        length, offset, linebases, linewidth = meta
        if start < 0 or end > length or start >= end:
            return ""
        start_line, start_col = divmod(start, linebases)
        end_line,   end_col   = divmod(end,   linebases)
        start_byte = offset + start_line * linewidth + start_col
        end_byte   = offset + end_line   * linewidth + end_col
        self.fh.seek(start_byte)
        chunk = self.fh.read(end_byte - start_byte).decode("ascii")
        return chunk.replace("\n", "").replace("\r", "")

    def close(self):
        try:
            self.fh.close()
        except Exception:
            pass


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--txt", required=True, help="ANNOVAR gnomAD txt")
    ap.add_argument("--ref", required=True,
                    help="Reference FASTA (needs companion .fai index)")
    ap.add_argument("--out", required=True,
                    help="Output unsorted VCF path (caller sorts/bgzips)")
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

    if not Path(args.ref + ".fai").is_file():
        print(f"ERROR: {args.ref}.fai not found. Build with:  samtools faidx {args.ref}",
              file=sys.stderr)
        return 2
    fasta = FastaRandom(args.ref)

    n_in = n_out = n_skip_af = n_skip_indel = n_skip_invalid = 0
    with open(txt_path, "r", encoding="utf-8") as fi, \
         open(out_path, "w", encoding="utf-8") as fo:

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

            fo.write(f"{chrom}\t{vcf_pos}\t.\t{vcf_ref}\t{vcf_alt}\t.\t.\t"
                     f"gnomAD_AF={af}\n")
            n_out += 1
            if n_in % 1_000_000 == 0:
                print(f"[annovar→vcf] processed {n_in:>11,}  wrote {n_out:>11,}",
                      file=sys.stderr)

    fasta.close()
    print(f"[annovar→vcf] read   {n_in:>11,} rows", file=sys.stderr)
    print(f"[annovar→vcf] wrote  {n_out:>11,} VCF lines", file=sys.stderr)
    print(f"[annovar→vcf] skip AF{n_skip_af:>12,}", file=sys.stderr)
    print(f"[annovar→vcf] skip indel{n_skip_indel:>9,}", file=sys.stderr)
    print(f"[annovar→vcf] skip invalid{n_skip_invalid:>7,}", file=sys.stderr)
    print(f"[annovar→vcf] done (unsorted) → {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
