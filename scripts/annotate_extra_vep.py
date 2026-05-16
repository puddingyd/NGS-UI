#!/usr/bin/env python3
"""Re-run VEP on a *filtered* TSV's sites to backfill extra in-silico
columns the main pipeline doesn't extract by default.

The tertiary pipeline runs VEP + dbNSFP but only pulls a subset of
dbNSFP fields (BayesDel_noAF, AlphaMissense, ESM1b, VARITY_R, SIFT,
DANN, PHACTboost, phyloP100way, GERP++_RS, gnomAD_exomes). This
script re-runs VEP on the small filtered set with MetaRNN (and
optionally SpliceAI when a scores VCF is available) added to the
dbNSFP/Plugin extract list, then joins those extra columns back
into the TSV.

Workflow:
    1. Read TSV → build sites VCF (CHROM/POS/REF/ALT, sorted).
    2. Run VEP via apptainer with --plugin dbNSFP,...,MetaRNN_score,
       MetaRNN_pred (+ SpliceAI if --spliceai-snv / --spliceai-indel
       supplied).
    3. Parse VEP VCF CSQ — pick the PICK=1 transcript per variant.
    4. Append METARNN_score / METARNN_pred (and SpliceAI columns if
       run) to the TSV. Existing values aren't touched.

Usage:
    scripts/annotate_extra_vep.py \\
        --tsv tertiary_output/VAL-58-dragen/snv_indel.annotated.tsv \\
        [--spliceai-snv  /path/to/spliceai.snv.vcf.gz \\
         --spliceai-indel /path/to/spliceai.indel.vcf.gz]

By default reads pipeline reference paths from the same locations as
nextflow_tertiary.config; flags below let you override.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# Default pipeline reference paths (mirrors nextflow_tertiary.config dgm profile).
DEFAULT_REF_DIR    = "/home/pipeline/reference/hg38"
DEFAULT_VEP_SIF    = "/home/pipeline/nextflow_containers/vep_115.sif"

# Extra dbNSFP fields we want on top of the pipeline's defaults.
EXTRA_DBNSFP_FIELDS = ["MetaRNN_score", "MetaRNN_pred"]

# TSV column names we'll write. Naming mirrors the pipeline's
# convention so the SNV adapter picks them up automatically.
TSV_COL_METARNN_SCORE = "METARNN"
TSV_COL_METARNN_PRED  = "METARNN_PRED"
TSV_COL_SPLICEAI      = "SPLICEAI_MAX"


def _open_vcf(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def tsv_to_sites(tsv_in: Path, vcf_out: Path) -> int:
    seen: set = set()
    rows: list = []
    with open(tsv_in, "r", encoding="utf-8", newline="") as fi:
        for r in csv.DictReader(fi, delimiter="\t"):
            chrom = (r.get("CHROM") or "").strip()
            pos   = (r.get("POS")   or "").strip()
            ref   = (r.get("REF")   or "").strip()
            alt   = (r.get("ALT")   or "").strip()
            if not (chrom and pos and ref and alt):
                continue
            if "*" in (ref, alt):
                continue
            key = (chrom, pos, ref, alt)
            if key in seen:
                continue
            seen.add(key)
            rows.append((chrom, int(pos), ref, alt))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(vcf_out, "w", encoding="utf-8") as fo:
        fo.write("##fileformat=VCFv4.2\n")
        fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in rows:
            fo.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
    return len(rows)


def run_vep(args, sites: Path, vep_out: Path) -> None:
    """Invoke vep via apptainer with the minimal set of plugins we need."""
    dbnsfp_fields = ",".join(EXTRA_DBNSFP_FIELDS)
    plugin_args = [f"--plugin", f"dbNSFP,{args.dbnsfp},{dbnsfp_fields}"]
    if args.spliceai_snv and args.spliceai_indel:
        plugin_args += ["--plugin",
                        f"SpliceAI,snv={args.spliceai_snv},indel={args.spliceai_indel}"]
    cmd = [
        "apptainer", "exec",
        "--bind", f"{args.ref_dir},{sites.parent},{Path(args.dbnsfp).parent}",
        args.vep_sif,
        "vep",
        "--input_file",  str(sites),
        "--output_file", str(vep_out),
        "--vcf",
        "--compress_output", "bgzip",
        "--offline", "--cache",
        "--dir_cache", args.vep_cache,
        "--dir_plugins", "/opt/vep/Plugins",
        "--assembly", "GRCh38",
        "--fasta", args.ref_fasta,
        "--fork", str(args.fork),
        "--symbol", "--canonical", "--biotype",
        "--mane", "--flag_pick",
        "--pick_order", "mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length",
        "--force_overwrite", "--no_stats", "--safe",
        *plugin_args,
    ]
    print(f"[extra-vep] running vep on {sites.name} …", file=sys.stderr)
    subprocess.run(cmd, check=True)


def _parse_csq_format(vcf_path: Path) -> list[str]:
    """Read the ##INFO=<ID=CSQ,...,Format="A|B|C"> line and return the
    pipe-separated field names so we can index by name."""
    pat = re.compile(r'##INFO=<ID=CSQ,.*Format="([^"]+)"', re.IGNORECASE)
    with _open_vcf(vcf_path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            m = pat.match(line.rstrip("\n"))
            if m:
                return m.group(1).split("|")
    return []


def _spliceai_max_from_csq(csq_field_value: str) -> str:
    """Take the worst-case SpliceAI delta score across the four classes.
    VEP's SpliceAI plugin emits something like 'A|GENE|0.01|0.02|0.03|0.04|...'
    where positions 2-5 are DS_AG / DS_AL / DS_DG / DS_DL — pick max abs."""
    if not csq_field_value:
        return ""
    parts = csq_field_value.split("|")
    nums = []
    for p in parts:
        try:
            nums.append(abs(float(p)))
        except ValueError:
            continue
    return f"{max(nums):.4f}" if nums else ""


def parse_vep_vcf(vep_vcf: Path) -> dict:
    """{(chrom, pos, ref, alt): {METARNN, METARNN_PRED, SPLICEAI_MAX}}.

    Picks the CSQ entry with PICK=1 if multiple transcripts; falls
    back to the first.
    """
    fmt = _parse_csq_format(vep_vcf)
    if not fmt:
        return {}
    idx_metarnn_score = fmt.index("MetaRNN_score") if "MetaRNN_score" in fmt else -1
    idx_metarnn_pred  = fmt.index("MetaRNN_pred")  if "MetaRNN_pred"  in fmt else -1
    idx_spliceai_pred = fmt.index("SpliceAI_pred") if "SpliceAI_pred" in fmt else -1
    idx_pick          = fmt.index("PICK")          if "PICK"          in fmt else -1
    out: dict = {}
    info_csq_pat = re.compile(r"(?:^|;)CSQ=([^;]+)")
    with _open_vcf(vep_vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            m = info_csq_pat.search(cols[7])
            if not m:
                continue
            csq_entries = m.group(1).split(",")
            picked = None
            for e in csq_entries:
                vals = e.split("|")
                if idx_pick >= 0 and idx_pick < len(vals) and vals[idx_pick] == "1":
                    picked = vals
                    break
            if picked is None:
                picked = csq_entries[0].split("|")
            row: dict = {}
            if idx_metarnn_score >= 0 and idx_metarnn_score < len(picked):
                row[TSV_COL_METARNN_SCORE] = picked[idx_metarnn_score]
            if idx_metarnn_pred >= 0 and idx_metarnn_pred < len(picked):
                row[TSV_COL_METARNN_PRED] = picked[idx_metarnn_pred]
            if idx_spliceai_pred >= 0 and idx_spliceai_pred < len(picked):
                row[TSV_COL_SPLICEAI] = _spliceai_max_from_csq(picked[idx_spliceai_pred])
            for alt_one in alt.split(","):
                out[(chrom, pos, ref, alt_one)] = row
    return out


def merge_into_tsv(in_tsv: Path, out_tsv: Path, ann: dict) -> tuple[int, int]:
    overwriting = in_tsv.resolve() == out_tsv.resolve()
    target = Path(str(out_tsv) + ".tmp") if overwriting else out_tsv
    target.parent.mkdir(parents=True, exist_ok=True)
    n_in = n_filled = 0
    with open(in_tsv, "r", encoding="utf-8", newline="") as fi:
        reader = csv.DictReader(fi, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for col in (TSV_COL_METARNN_SCORE, TSV_COL_METARNN_PRED, TSV_COL_SPLICEAI):
            if col not in fieldnames:
                fieldnames.append(col)
        with open(target, "w", encoding="utf-8", newline="") as fo:
            writer = csv.DictWriter(fo, fieldnames=fieldnames, delimiter="\t",
                                    extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                n_in += 1
                key = (
                    (row.get("CHROM") or "").strip(),
                    (row.get("POS") or "").strip(),
                    (row.get("REF") or "").strip(),
                    (row.get("ALT") or "").strip(),
                )
                hit = ann.get(key)
                if hit:
                    changed = False
                    for k, v in hit.items():
                        if v and not (row.get(k) or "").strip():
                            row[k] = v
                            changed = True
                    if changed:
                        n_filled += 1
                writer.writerow(row)
    if overwriting:
        os.replace(target, out_tsv)
    return n_in, n_filled


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--tsv", required=True,
                    help="Filtered SNV TSV to backfill (in-place unless --out-tsv)")
    ap.add_argument("--out-tsv",
                    help="Write merged TSV here instead of overwriting --tsv")
    ap.add_argument("--workdir",
                    help="Keep VEP intermediates here (default: temp dir)")
    ap.add_argument("--ref-dir",   default=DEFAULT_REF_DIR)
    ap.add_argument("--vep-sif",   default=DEFAULT_VEP_SIF)
    ap.add_argument("--vep-cache", default=None,
                    help="default $REF_DIR/tertiary/vep_cache")
    ap.add_argument("--dbnsfp",    default=None,
                    help="default $REF_DIR/tertiary/dbnsfp/dbNSFP4.9c_grch38.gz")
    ap.add_argument("--ref-fasta", default=None,
                    help="default $REF_DIR/Homo_sapiens_assembly38.fasta")
    ap.add_argument("--fork", type=int, default=4)
    ap.add_argument("--spliceai-snv",
                    help="Path to spliceai_scores.raw.snv.hg38.vcf.gz "
                         "(if absent, SpliceAI step is skipped)")
    ap.add_argument("--spliceai-indel",
                    help="Path to spliceai_scores.raw.indel.hg38.vcf.gz")
    args = ap.parse_args()

    if args.vep_cache is None:
        args.vep_cache = f"{args.ref_dir}/tertiary/vep_cache"
    if args.dbnsfp is None:
        args.dbnsfp = f"{args.ref_dir}/tertiary/dbnsfp/dbNSFP4.9c_grch38.gz"
    if args.ref_fasta is None:
        args.ref_fasta = f"{args.ref_dir}/Homo_sapiens_assembly38.fasta"

    in_tsv = Path(args.tsv).resolve()
    if not in_tsv.is_file():
        print(f"ERROR: --tsv 找不到：{in_tsv}", file=sys.stderr)
        return 2
    out_tsv = Path(args.out_tsv).resolve() if args.out_tsv else in_tsv

    if args.workdir:
        wd = Path(args.workdir); wd.mkdir(parents=True, exist_ok=True)
        wd_ctx = None
    else:
        wd_ctx = tempfile.TemporaryDirectory(prefix="extra-vep-")
        wd = Path(wd_ctx.name)

    try:
        sites   = wd / "sites.vcf"
        vep_vcf = wd / "sites.vep.vcf.gz"
        n = tsv_to_sites(in_tsv, sites)
        print(f"[extra-vep] {n} unique sites → {sites}", file=sys.stderr)
        if n == 0:
            print("ERROR: 0 sites", file=sys.stderr)
            return 1
        run_vep(args, sites, vep_vcf)
        ann = parse_vep_vcf(vep_vcf)
        if args.spliceai_snv:
            print(f"[extra-vep] SpliceAI plugin enabled "
                  f"(snv={args.spliceai_snv}, indel={args.spliceai_indel})",
                  file=sys.stderr)
        else:
            print("[extra-vep] SpliceAI skipped (no --spliceai-snv supplied)",
                  file=sys.stderr)
        n_in, n_filled = merge_into_tsv(in_tsv, out_tsv, ann)
        print(f"[extra-vep] backfilled {n_filled}/{n_in} rows from VEP ({len(ann)} sites annotated)",
              file=sys.stderr)
        print(f"[extra-vep] done → {out_tsv}", file=sys.stderr)
    finally:
        if wd_ctx is not None:
            wd_ctx.cleanup()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
