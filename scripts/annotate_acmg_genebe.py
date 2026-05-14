#!/usr/bin/env python3
"""GeneBe ACMG annotation — sites-only wrapper.

Why "sites-only": pygenebe doesn't handle multi-sample VCFs cleanly,
and we don't need sample columns for ACMG anyway. Strip to a one-row-
per-variant sites VCF, send it to pygenebe (via apptainer), parse the
INFO.acmg_score / INFO.acmg_criteria back into a small TSV keyed on
(CHROM, POS, REF, ALT).

Inputs (mutually exclusive):
    --tsv snv_indel.annotated.tsv       use the TSV's CHROM/POS/REF/ALT
    --vcf any.vcf[.gz]                  drop FORMAT + samples to sites-only

Output:
    <outdir>/acmg_genebe.tsv   cols: CHROM POS REF ALT ACMG_score
                                      ACMG_criteria ACMG_class

Credentials (env or flags):
    GENEBE_USER       --username
    GENEBE_API_KEY    --api_key
    GENEBE_SIF        --sif        default $HOME/NGS_UI/biotools/genebe.sif

Score → class mapping mirrors the R workflow:
    score >= 10        Pathogenic
    score 6..9         Likely pathogenic
    score 0..5         Uncertain significance
    score -1..-6       Likely benign
    score <= -7        Benign
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
from pathlib import Path

PAT_SCORE = re.compile(r"(?:^|;)acmg_score=([^;]+)")
PAT_CRIT  = re.compile(r"(?:^|;)acmg_criteria=([^;]+)")


def _open(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "r")


def _row_max_af(row: dict, cols: list[str]) -> float:
    """Largest numeric AF across the listed columns; missing/non-numeric → 0.

    "Missing means rare" matches ACMG BA1/BS1 practice: an AF that
    gnomAD doesn't carry isn't evidence of commonness.
    """
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


def tsv_to_sites(
    tsv_in: str,
    vcf_out: Path,
    *,
    max_af: float | None = None,
    af_cols: list[str] | None = None,
) -> tuple[int, int]:
    """snv_indel.annotated.tsv → minimal sites-only VCF.

    Returns (n_written, n_dropped_by_af). max_af=None disables the AF
    filter.
    """
    af_cols = af_cols or ["GNOMAD_G_AF", "GNOMAD_E_AF"]
    seen: set = set()
    dropped = 0
    with open(tsv_in, "r", encoding="utf-8", newline="") as fi, \
         open(vcf_out, "w", encoding="utf-8") as fo:
        fo.write("##fileformat=VCFv4.2\n")
        fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in csv.DictReader(fi, delimiter="\t"):
            chrom = (row.get("CHROM") or "").strip()
            pos   = (row.get("POS")   or "").strip()
            ref   = (row.get("REF")   or "").strip()
            alt   = (row.get("ALT")   or "").strip()
            if not (chrom and pos and ref and alt):
                continue
            key = (chrom, pos, ref, alt)
            if key in seen:
                continue
            if max_af is not None and _row_max_af(row, af_cols) > max_af:
                dropped += 1
                continue
            seen.add(key)
            fo.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
    return len(seen), dropped


def vcf_to_sites(vcf_in: str, vcf_out: Path) -> int:
    """Strip FORMAT + sample columns; keep ## headers + one row per site."""
    seen: set = set()
    with _open(vcf_in) as fi, open(vcf_out, "w", encoding="utf-8") as fo:
        for line in fi:
            if line.startswith("##"):
                fo.write(line)
                continue
            if line.startswith("#CHROM"):
                fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, flt, info = f[:8]
            key = (chrom, pos, ref, alt)
            if key in seen:
                continue
            seen.add(key)
            fo.write("\t".join([chrom, pos, vid or ".", ref, alt,
                                qual or ".", flt or ".", info or "."]) + "\n")
    return len(seen)


def classify(score: float | None) -> str:
    if score is None:
        return ""
    if score >= 10:  return "Pathogenic"
    if score >=  6:  return "Likely pathogenic"
    if score >=  0:  return "Uncertain significance"
    if score >= -6:  return "Likely benign"
    return "Benign"


def parse_annotated_vcf(annot_vcf: Path, tsv_out: Path) -> int:
    n = 0
    with open(annot_vcf, "r", encoding="utf-8") as fi, \
         open(tsv_out, "w", encoding="utf-8", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["CHROM", "POS", "REF", "ALT",
                    "ACMG_score", "ACMG_criteria", "ACMG_class"])
        for line in fi:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _, ref, alt, _, _, info = parts[:8]
            ms = PAT_SCORE.search(info)
            mc = PAT_CRIT.search(info)
            score_raw = ms.group(1) if ms else ""
            try:
                score = float(score_raw) if score_raw else None
            except ValueError:
                score = None
            crit = mc.group(1) if mc else ""
            w.writerow([chrom, pos, ref, alt, score_raw, crit, classify(score)])
            n += 1
    return n


def ensure_sif(sif: Path) -> None:
    if sif.is_file():
        return
    sif.parent.mkdir(parents=True, exist_ok=True)
    print(f"[genebe] pulling docker://genebe/pygenebe:0.0.18 → {sif}", file=sys.stderr)
    subprocess.run(
        ["apptainer", "pull", "--force", str(sif),
         "docker://genebe/pygenebe:0.0.18"],
        check=True,
    )


def run_genebe(sif: Path, sites: Path, annot: Path, user: str, key: str) -> None:
    cmd = [
        "apptainer", "exec", "--bind", str(sites.parent), str(sif),
        "genebe", "annotate", "--genome", "hg38",
        "--input", str(sites), "--output", str(annot),
        "--username", user, "--api_key", key,
    ]
    print(f"[genebe] running pygenebe (sites={sites.name})", file=sys.stderr)
    subprocess.run(cmd, check=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--tsv", help="snv_indel.annotated.tsv (uses CHROM/POS/REF/ALT)")
    g.add_argument("--vcf", help="any VCF[.gz] — multi-sample OK, samples are dropped")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--username", default=os.environ.get("GENEBE_USER"))
    ap.add_argument("--api_key",  default=os.environ.get("GENEBE_API_KEY"))
    ap.add_argument("--sif",
                    default=os.environ.get(
                        "GENEBE_SIF",
                        str(Path.home() / "NGS_UI" / "biotools" / "genebe.sif")))
    ap.add_argument("--max-af", type=float, default=0.01,
                    help="drop sites where any listed AF column exceeds this "
                         "(default 0.01; use --max-af -1 to disable). "
                         "Only applies to --tsv input.")
    ap.add_argument("--af-cols",
                    default="GNOMAD_G_AF,GNOMAD_E_AF",
                    help="comma-separated AF columns to check for --max-af "
                         "(default GNOMAD_G_AF,GNOMAD_E_AF)")
    args = ap.parse_args()

    if not args.username or not args.api_key:
        print("ERROR: 需要 --username + --api_key（或 env GENEBE_USER/GENEBE_API_KEY）",
              file=sys.stderr)
        return 2

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    sites = outdir / "sites.vcf"
    annot = outdir / "sites.genebe.vcf"
    out   = outdir / "acmg_genebe.tsv"

    if args.tsv:
        max_af = None if args.max_af < 0 else args.max_af
        af_cols = [c.strip() for c in args.af_cols.split(",") if c.strip()]
        n_sites, n_dropped = tsv_to_sites(args.tsv, sites,
                                          max_af=max_af, af_cols=af_cols)
        af_note = (f"  (dropped {n_dropped} with max({','.join(af_cols)}) > {max_af})"
                   if max_af is not None else "  (AF filter off)")
    else:
        n_sites = vcf_to_sites(args.vcf, sites)
        n_dropped = 0
        af_note = ""
    print(f"[genebe] {n_sites} unique sites → {sites}{af_note}", file=sys.stderr)
    if n_sites == 0:
        print("ERROR: 0 sites — 檢查輸入 TSV/VCF 是否有 CHROM/POS/REF/ALT 欄，"
              "或 --max-af 過嚴。",
              file=sys.stderr)
        return 1

    ensure_sif(Path(args.sif))
    run_genebe(Path(args.sif), sites, annot, args.username, args.api_key)
    n_out = parse_annotated_vcf(annot, out)
    print(f"[genebe] done → {out}  ({n_out} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
