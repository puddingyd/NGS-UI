#!/usr/bin/env python3
"""Phase A — per-sample ingest for the incremental in-house AF DB.

From ONE DRAGEN sample (its indexed hard-filtered gVCF + ploidy CSV) produce the
per-sample contributions that Phase B accumulates:

  per_sample/{id}/
    callable.weighted.bed.gz   DP>=10 intervals, ploidy-weighted (-> AN track)
    counts.tsv                 one row per PASS ALT: chrom pos ref alt klass af
    qc.json                    sex call, per-region callable bp, variant tallies

Design: scripts/inhouse_af/DESIGN_incremental.md (§4.1/4.2/4.4/4.6). This script
does NOT accumulate or publish — that is Phase B.

  * Sex/karyotype from DRAGEN `ploidy_estimation_metrics.csv` (`Ploidy
    estimation,XX/XY/...`); coverage fallback not implemented here (all cohort
    samples have the CSV).
  * Callable + AC both come from the single indexed gVCF (ref blocks carry
    MIN_DP; PASS variant rows carry FILTER+GT; chrM included).
  * Variant normalization via `bcftools norm -m- -f ref --check-ref x` so keys
    match the GLnexus path and the annotation TSV.

Stdlib only (gzip streaming); shells out to bcftools for normalization.

Usage:
  scripts/inhouse_af/ingest_sample.py \
    --gvcf  /home/datalake_Raw/Novaseq/<run>/other/<id>/<id>.hard-filtered.gvcf.gz \
    --ploidy-csv .../<id>.ploidy_estimation_metrics.csv \
    --ref   /home/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta \
    --out-dir $NGS_UI_HOME/biotools/inhouse_af/per_sample

  scripts/inhouse_af/ingest_sample.py --selftest   # run unit checks, no I/O
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone, timedelta

# hg38 pseudoautosomal regions (1-based, inclusive)
PAR1 = (10001, 2781479)
PAR2 = (155701383, 156030895)
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
PRIMARY = AUTOSOMES | {"chrX", "chrY", "chrM"}
MT_HOMOPLASMIC_AF = 0.95   # chrM AF>=this -> homoplasmic, else heteroplasmic
TAIPEI = timezone(timedelta(hours=8))

# --------------------------------------------------------------------------
# pure helpers (unit-tested by --selftest)
# --------------------------------------------------------------------------

def parse_karyotype(csv_path: str) -> str:
    """Return the DRAGEN 'Ploidy estimation' value (e.g. 'XX'), or '' if absent."""
    try:
        with open(csv_path, encoding="utf-8") as f:
            kary = ""
            for line in f:
                parts = line.rstrip("\n").split(",")
                if len(parts) >= 4 and parts[2].strip() == "Ploidy estimation":
                    kary = parts[3].strip()
            return kary
    except OSError:
        return ""


def sex_class_of(karyotype: str) -> str:
    """Map a karyotype to female / male / ambiguous.

    Only clean XX/XY get X/Y ploidy weights; X0/XXX/XXY/blank/other are
    'ambiguous' -> autosome+chrM counted normally, X/Y excluded (conservative)."""
    k = (karyotype or "").strip().upper()
    if k == "XX":
        return "female"
    if k == "XY":
        return "male"
    return "ambiguous"


def region_of(chrom: str, pos: int) -> str:
    """auto | xpar | xnonpar | y | mt | other (pos is 1-based)."""
    if chrom in AUTOSOMES:
        return "auto"
    if chrom == "chrX":
        if PAR1[0] <= pos <= PAR1[1] or PAR2[0] <= pos <= PAR2[1]:
            return "xpar"
        return "xnonpar"
    if chrom == "chrY":
        return "y"
    if chrom == "chrM":
        return "mt"
    return "other"


def callable_weight(region: str, sex: str) -> int:
    """Ploidy weight a callable base contributes to AN (0 = skip)."""
    if region == "auto" or region == "xpar":
        return 2
    if region == "xnonpar":
        return 2 if sex == "female" else (1 if sex == "male" else 0)
    if region == "y":
        return 1 if sex == "male" else 0
    if region == "mt":
        return 1
    return 0


def classify_variant(region: str, sex: str, gt: str, af):
    """Return klass for a PASS biallelic ALT, or None to skip.

    nuclear: hom | het | hemi ; chrM: mt_hom (homoplasmic) | mt_het."""
    if region == "mt":
        return "mt_hom" if (af is not None and af >= MT_HOMOPLASMIC_AF) else "mt_het"
    if region in ("xnonpar", "y") and sex == "male":
        return "hemi"
    if region in ("xnonpar", "y", "xpar") and sex == "ambiguous":
        return None
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) == 1:
        return "hemi"
    nalt = sum(1 for a in alleles if a not in ("0", ".", ""))
    return "hom" if nalt >= 2 else "het"


def _fmt_value(fmt: str, sample: str, key: str):
    keys = fmt.split(":")
    if key not in keys:
        return None
    vals = sample.split(":")
    i = keys.index(key)
    return vals[i] if i < len(vals) else None


def _real_alts(alt: str):
    return [a for a in alt.split(",") if a not in ("<NON_REF>", "<*>", ".", "")]

# --------------------------------------------------------------------------
# gVCF streaming pass: callable weighted BED + extract PASS variants
# --------------------------------------------------------------------------

def stream_gvcf(gvcf, sex, min_dp, bed_fh, variants_vcf_fh):
    """One pass: write callable weighted BED + copy header + PASS variant lines
    (primary contigs only) to a temp VCF for normalization.
    Returns dict of callable bp per region group."""
    cbp = {"auto": 0, "x": 0, "y": 0, "mt": 0}
    opener = gzip.open if gvcf.endswith(".gz") else open
    with opener(gvcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                variants_vcf_fh.write(line)
                continue
            c = line.find("\t")
            # cheap field split
            F = line.rstrip("\n").split("\t")
            if len(F) < 10:
                continue
            chrom, pos_s, _id, ref, alt, _qual, filt, info, fmt, sample = F[:10]
            if chrom not in PRIMARY:
                continue
            pos = int(pos_s)
            reals = _real_alts(alt)
            is_var = bool(reals)
            region = region_of(chrom, pos)
            w = callable_weight(region, sex)
            # depth: MIN_DP for ref blocks, DP for variant rows
            if is_var:
                dp = _fmt_value(fmt, sample, "DP")
            else:
                dp = _fmt_value(fmt, sample, "MIN_DP")
            try:
                depth = int(dp) if dp not in (None, ".", "") else None
            except ValueError:
                depth = None
            # callable interval
            if w > 0:
                is_callable = (depth is not None and depth >= min_dp) or (is_var and depth is None)
                if is_callable:
                    if is_var:
                        # single anchor base only. A multi-base REF span (e.g. a
                        # deletion TG>T) would overlap the adjacent ref blocks of
                        # the SAME sample and double-count AN at those bases; we
                        # only ever look up AN at the variant's POS anyway.
                        start, end = pos - 1, pos
                    else:
                        end_s = ""
                        if "END=" in info:
                            end_s = info.split("END=", 1)[1].split(";", 1)[0]
                        end = int(end_s) if end_s else pos
                        start = pos - 1
                    bed_fh.write(f"{chrom}\t{start}\t{end}\t{w}\n")
                    grp = "x" if region in ("xpar", "xnonpar") else ("y" if region == "y" else region)
                    cbp[grp] = cbp.get(grp, 0) + (end - start)
            # stash PASS variant rows for normalization
            if is_var and filt == "PASS":
                variants_vcf_fh.write(line)
    return cbp

# --------------------------------------------------------------------------
# variant counting via bcftools norm
# --------------------------------------------------------------------------

def normalize_and_count(variants_vcf, ref, sex, bcftools_cmd):
    """norm -m- the extracted PASS variants and classify each ALT.
    Returns (rows, tally) where rows = [(chrom,pos,ref,alt,klass,af_str)]."""
    norm = bcftools_cmd + ["norm", "-m-", "-f", ref, "--check-ref", "x", variants_vcf]
    query = bcftools_cmd + ["query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%AF]\n"]
    p1 = subprocess.Popen(norm, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(query, stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
    p1.stdout.close()
    rows, tally = [], {}
    for line in p2.stdout:
        chrom, pos, ref_a, alt, gt, af_s = line.rstrip("\n").split("\t")
        if alt in ("<NON_REF>", "<*>", ".", ""):
            continue
        if chrom not in PRIMARY:
            continue
        try:
            af = float(af_s) if af_s not in (".", "") else None
        except ValueError:
            af = None
        region = region_of(chrom, int(pos))
        klass = classify_variant(region, sex, gt, af)
        if klass is None:
            continue
        rows.append((chrom, pos, ref_a, alt, klass, af_s if region == "mt" else "."))
        tally[klass] = tally.get(klass, 0) + 1
    p2.wait()
    return rows, tally


def resolve_bcftools(args):
    sif = os.environ.get("BCFTOOLS_SIF", args.bcftools_sif)
    if sif:
        binds = ",".join(filter(None, [os.environ.get("APPTAINER_BIND", "/home")]))
        return ["apptainer", "exec", "--bind", binds, sif, "bcftools"]
    binp = os.environ.get("BCFTOOLS_BIN", args.bcftools_bin)
    if shutil.which(binp):
        return [binp]
    sys.exit("ERROR: need bcftools — set BCFTOOLS_SIF or put bcftools on PATH (BCFTOOLS_BIN)")

# --------------------------------------------------------------------------

def run(args):
    sample_id = args.sample_id or os.path.basename(args.gvcf).split(".")[0]
    karyotype = parse_karyotype(args.ploidy_csv) if args.ploidy_csv else ""
    sex = sex_class_of(karyotype)
    bcftools_cmd = resolve_bcftools(args)

    out = os.path.join(args.out_dir, sample_id)
    os.makedirs(out, exist_ok=True)
    bed_path = os.path.join(out, "callable.weighted.bed")
    counts_path = os.path.join(out, "counts.tsv")
    qc_path = os.path.join(out, "qc.json")

    print(f"[ingest] {sample_id}  karyotype={karyotype or 'NA'} -> {sex}", file=sys.stderr)

    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as tv, \
            open(bed_path, "w") as bed:
        tmp_vcf = tv.name
        cbp = stream_gvcf(args.gvcf, sex, args.min_dp, bed, tv)

    rows, tally = normalize_and_count(tmp_vcf, args.ref, sex, bcftools_cmd)
    os.unlink(tmp_vcf)

    with open(counts_path, "w") as cf:
        cf.write("chrom\tpos\tref\talt\tklass\taf\n")
        for r in rows:
            cf.write("\t".join(r) + "\n")

    # bgzip the BED (best-effort; plain BED is fine for Phase B too)
    try:
        subprocess.run(bcftools_cmd[:-1] + ["bgzip", "-f", bed_path]
                       if bcftools_cmd[-1] == "bcftools" else ["bgzip", "-f", bed_path],
                       check=True, stderr=subprocess.DEVNULL)
        bed_out = bed_path + ".gz"
    except Exception:
        bed_out = bed_path  # leave uncompressed

    qc = {
        "sample_id": sample_id,
        "sex_source": "dragen_ploidy" if karyotype else "none",
        "karyotype": karyotype,
        "sex_class": sex,
        "min_dp": args.min_dp,
        "callable_bp": cbp,
        "n_variants": tally,
        "gvcf": os.path.abspath(args.gvcf),
        "ingested_at": datetime.now(TAIPEI).strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(qc_path, "w") as qf:
        json.dump(qc, qf, indent=2)

    print(f"[ingest] {sample_id}  callable_bp={cbp}  variants={tally}", file=sys.stderr)
    print(f"[ingest] -> {bed_out}\n          {counts_path}\n          {qc_path}", file=sys.stderr)
    return 0

# --------------------------------------------------------------------------

def selftest():
    assert region_of("chr7", 100) == "auto"
    assert region_of("chrX", 1_000_000) == "xpar"
    assert region_of("chrX", 60_000_000) == "xnonpar"
    assert region_of("chrX", 155_800_000) == "xpar"
    assert region_of("chrY", 5_000_000) == "y"
    assert region_of("chrM", 73) == "mt"
    assert region_of("chr1_KI270706v1_random", 5) == "other"

    assert callable_weight("auto", "male") == 2
    assert callable_weight("xnonpar", "female") == 2
    assert callable_weight("xnonpar", "male") == 1
    assert callable_weight("xnonpar", "ambiguous") == 0
    assert callable_weight("y", "male") == 1
    assert callable_weight("y", "female") == 0
    assert callable_weight("mt", "female") == 1
    assert callable_weight("other", "male") == 0

    # nuclear classification
    assert classify_variant("auto", "female", "1/1", None) == "hom"
    assert classify_variant("auto", "female", "0/1", None) == "het"
    assert classify_variant("auto", "male", "0|1", None) == "het"
    assert classify_variant("xnonpar", "male", "1/1", None) == "hemi"
    assert classify_variant("xnonpar", "male", "1", None) == "hemi"
    assert classify_variant("y", "male", "1", None) == "hemi"
    assert classify_variant("xnonpar", "female", "1/1", None) == "hom"
    assert classify_variant("xpar", "male", "0/1", None) == "het"
    assert classify_variant("xnonpar", "ambiguous", "1/1", None) is None
    # chrM
    assert classify_variant("mt", "female", "1/1", 1.0) == "mt_hom"
    assert classify_variant("mt", "female", "1/1", 0.99) == "mt_hom"
    assert classify_variant("mt", "female", "0/1", 0.054) == "mt_het"
    assert classify_variant("mt", "male", "0/1", 0.19) == "mt_het"

    # FORMAT field extraction
    assert _fmt_value("GT:AD:DP:GQ:MIN_DP", "0/0:42,2:44:24:41", "MIN_DP") == "41"
    assert _fmt_value("GT:GQ:PS", "0|1:54:2540", "DP") is None
    assert _real_alts("G,<NON_REF>") == ["G"]
    assert _real_alts("<NON_REF>") == []
    assert _real_alts("AC,ACC,<NON_REF>") == ["AC", "ACC"]

    assert sex_class_of("XX") == "female"
    assert sex_class_of("XY") == "male"
    assert sex_class_of("X0") == "ambiguous"
    assert sex_class_of("XXY") == "ambiguous"
    assert sex_class_of("") == "ambiguous"
    print("selftest OK")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gvcf")
    ap.add_argument("--ploidy-csv")
    ap.add_argument("--sample-dir",
                    help="DRAGEN other/{id}/ dir; auto-finds the gVCF + ploidy CSV "
                         "under it (handles the germline_seq/ nesting)")
    ap.add_argument("--ref", default="/home/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta")
    ap.add_argument("--out-dir")
    ap.add_argument("--sample-id")
    ap.add_argument("--min-dp", type=int, default=10)
    ap.add_argument("--bcftools-sif", default="")
    ap.add_argument("--bcftools-bin", default="bcftools")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    if args.sample_dir:
        import glob
        if not args.gvcf:
            hits = sorted(glob.glob(os.path.join(args.sample_dir, "**", "*.hard-filtered.gvcf.gz"),
                                    recursive=True))
            if not hits:
                ap.error(f"no *.hard-filtered.gvcf.gz under {args.sample_dir}")
            args.gvcf = hits[0]
        if not args.ploidy_csv:
            pc = sorted(glob.glob(os.path.join(args.sample_dir, "**", "*.ploidy_estimation_metrics.csv"),
                                  recursive=True))
            args.ploidy_csv = pc[0] if pc else None
    if not (args.gvcf and args.out_dir):
        ap.error("--gvcf (or --sample-dir) and --out-dir are required (or use --selftest)")
    if not os.path.exists(args.gvcf):
        ap.error(f"--gvcf not found: {args.gvcf}")
    if not os.path.exists(args.ref):
        ap.error(f"--ref not found: {args.ref}")
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
