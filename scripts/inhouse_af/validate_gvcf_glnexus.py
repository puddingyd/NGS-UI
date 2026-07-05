#!/usr/bin/env python3
"""Pre-flight: does a DRAGEN gVCF have records GLnexus will reject on bulk load?

GLnexus aborts the ENTIRE joint-genotyping run if any one input has a malformed
record ("gVCF record doesn't have expected # of GT entries"). Finding that out
11 hours into a bulk load is painful, so this streams a gVCF and flags the same
class of defect up front — a diploid/haploid human gVCF genotype must have 1 or
2 allele entries, and every allele index must be within the record's ALT list.

Deliberately stdlib-only (gzip streaming, index-free) so it runs on the
air-gapped DGX with no bcftools dependency. It reads from --gvcf directly, or
from stdin (`--stdin`) so a caller can pipe fast C decompression in:

    bgzip -dc sample.gvcf.gz | validate_gvcf_glnexus.py --stdin --sample-id S

Exit code: 0 = clean, 2 = malformed record(s) found (first few printed to
stderr as `chrom:pos`), 1 = usage/IO error. Use --selftest for unit checks.
"""
from __future__ import annotations

import argparse
import gzip
import sys


def gt_entry_count_ok(gt: str, n_alleles: int):
    """Return (ok, reason). `gt` is the raw GT subfield; n_alleles = 1+len(ALTs).

    GLnexus expects a human gVCF genotype to be haploid or diploid, with every
    numeric allele index < n_alleles. Missing ('.') entries are fine."""
    # strip phasing; the GT is the first FORMAT subfield already isolated by caller
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) not in (1, 2):
        return False, f"GT has {len(alleles)} entries ({gt!r})"
    for a in alleles:
        if a in (".", ""):
            continue
        if not a.isdigit():
            return False, f"non-numeric allele {a!r} in GT {gt!r}"
        if int(a) >= n_alleles:
            return False, f"allele index {a} >= n_alleles {n_alleles} (GT {gt!r})"
    return True, ""


def scan(fh, max_report=5):
    """Scan an open text gVCF stream; return list of (chrom, pos, reason)."""
    bad = []
    gt_idx = None
    for line in fh:
        if line.startswith("#"):
            continue
        F = line.rstrip("\n").split("\t")
        if len(F) < 10:
            continue
        chrom, pos, _id, ref, alt, _q, _f, _info, fmt, sample = F[:10]
        # GT is by spec the first FORMAT field, but locate it defensively.
        keys = fmt.split(":")
        gi = 0 if keys and keys[0] == "GT" else (keys.index("GT") if "GT" in keys else None)
        if gi is None:
            continue
        vals = sample.split(":")
        if gi >= len(vals):
            bad.append((chrom, pos, "sample has no GT value"))
            if len(bad) >= max_report:
                break
            continue
        # n_alleles = REF + real ALTs (ignore gVCF symbolic <NON_REF>/<*>)
        alts = [a for a in alt.split(",") if a not in ("<NON_REF>", "<*>")]
        n_alleles = 1 + len(alts)
        ok, reason = gt_entry_count_ok(vals[gi], n_alleles)
        if not ok:
            bad.append((chrom, pos, reason))
            if len(bad) >= max_report:
                break
    return bad


def selftest():
    assert gt_entry_count_ok("0/1", 2)[0]
    assert gt_entry_count_ok("0|1", 2)[0]
    assert gt_entry_count_ok("1", 2)[0]          # haploid (chrX/Y/M)
    assert gt_entry_count_ok("./.", 2)[0]
    assert gt_entry_count_ok(".", 1)[0]
    assert not gt_entry_count_ok("0/1/1", 2)[0]  # 3 entries -> the GLnexus defect
    assert not gt_entry_count_ok("1/2", 2)[0]    # allele index 2 with only 1 ALT
    assert gt_entry_count_ok("1/2", 3)[0]        # ok when 2 ALTs present
    assert not gt_entry_count_ok("0/x", 2)[0]    # non-numeric
    # stream scan over a tiny synthetic gVCF
    import io
    doc = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "chr1\t100\t.\tA\t<NON_REF>\t.\t.\tEND=200\tGT:MIN_DP\t0/0:30\n"
        "chr1\t250\t.\tA\tG,<NON_REF>\t.\tPASS\t.\tGT:DP\t0/1:40\n"
        "chr1\t300\t.\tA\tG,<NON_REF>\t.\tPASS\t.\tGT:DP\t0/1/1:40\n"  # bad
    )
    bad = scan(io.StringIO(doc))
    assert bad == [("chr1", "300", "GT has 3 entries ('0/1/1')")], bad
    print("selftest OK")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gvcf", help="path to *.hard-filtered.gvcf.gz")
    ap.add_argument("--stdin", action="store_true", help="read an uncompressed gVCF stream from stdin")
    ap.add_argument("--sample-id", default="", help="label for reporting")
    ap.add_argument("--max-report", type=int, default=5)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return selftest()

    sid = args.sample_id
    if args.stdin:
        bad = scan(sys.stdin, args.max_report)
    elif args.gvcf:
        opener = gzip.open if args.gvcf.endswith(".gz") else open
        sid = sid or args.gvcf.split("/")[-1].split(".")[0]
        with opener(args.gvcf, "rt") as f:
            bad = scan(f, args.max_report)
    else:
        ap.error("need --gvcf or --stdin (or --selftest)")

    if bad:
        for chrom, pos, reason in bad:
            print(f"[BAD] {sid}\t{chrom}:{pos}\t{reason}", file=sys.stderr)
        return 2
    print(f"[ok] {sid}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
