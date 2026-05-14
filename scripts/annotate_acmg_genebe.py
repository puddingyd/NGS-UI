#!/usr/bin/env python3
"""GeneBe ACMG annotation — backfill into snv_indel.annotated.tsv.

Reads `--tsv snv_indel.annotated.tsv`, builds a sites-only VCF
(skipping spanning-deletion '*' alleles and high-AF variants per
--max-af), runs pygenebe via apptainer, and writes the GeneBe
ACMG_score / ACMG_criteria / ACMG_class back into the TSV's
existing ACMG_POINTS / ACMG_EVIDENCE / ACMG_CLASS columns. Only
empty cells are filled — pipeline-supplied ACMG always wins, so the
day the tertiary pipeline's ACMG_CLASSIFY step ships you simply stop
running this and the TSV is already self-contained.

By default updates --tsv in place; use --out-tsv to write elsewhere.

Score → class mapping (GeneBe ACMG_score → 5-tier label):
    >= 10        Pathogenic
    6..9         Likely pathogenic
    0..5         Uncertain significance
    -1..-6       Likely benign
    <= -7        Benign

Credentials (env or flags):
    GENEBE_USER       --username
    GENEBE_API_KEY    --api_key
    GENEBE_SIF        --sif (default $HOME/NGS_UI/biotools/genebe.sif)
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

PAT_SCORE = re.compile(r"(?:^|;)acmg_score=([^;]+)")
PAT_CRIT  = re.compile(r"(?:^|;)acmg_criteria=([^;]+)")


def classify(score: float | None) -> str:
    if score is None:
        return ""
    if score >= 10:  return "Pathogenic"
    if score >=  6:  return "Likely pathogenic"
    if score >=  0:  return "Uncertain significance"
    if score >= -6:  return "Likely benign"
    return "Benign"


def _row_max_af(row: dict, cols: list[str]) -> float:
    """Largest numeric AF across the listed columns; missing → 0.

    Matches ACMG BA1/BS1 practice: a missing AF in gnomAD isn't
    evidence of commonness, so it's kept.
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
    tsv_in: Path,
    vcf_out: Path,
    *,
    max_af: float | None,
    af_cols: list[str],
) -> tuple[int, int, int]:
    """snv_indel.annotated.tsv → minimal sites-only VCF.

    Returns (n_written, n_dropped_by_af, n_skipped_star).
    """
    seen: set = set()
    n_af = 0
    n_star = 0
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
            if "*" in (ref, alt):
                n_star += 1
                continue
            key = (chrom, pos, ref, alt)
            if key in seen:
                continue
            if max_af is not None and _row_max_af(row, af_cols) > max_af:
                n_af += 1
                continue
            seen.add(key)
            fo.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
    return len(seen), n_af, n_star


def parse_annotated_vcf(annot_vcf: Path) -> dict[tuple[str, str, str, str],
                                                  tuple[str, str, str]]:
    """{(chr, pos, ref, alt): (score_str, criteria, class_label)}."""
    out: dict = {}
    with open(annot_vcf, "r", encoding="utf-8") as fi:
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
                score_num = float(score_raw) if score_raw else None
            except ValueError:
                score_num = None
            crit = mc.group(1) if mc else ""
            out[(chrom, pos, ref, alt)] = (score_raw, crit, classify(score_num))
    return out


def merge_into_tsv(
    in_tsv: Path,
    out_tsv: Path,
    gb: dict,
) -> tuple[int, int]:
    """Write a new TSV with ACMG_* backfilled. Returns (n_filled, n_total).

    Only blank ACMG cells are populated — any non-empty pipeline value
    is preserved, so when ACMG_CLASSIFY ships nothing changes.
    """
    in_tsv = Path(in_tsv)
    out_tsv = Path(out_tsv)
    overwriting = in_tsv.resolve() == out_tsv.resolve()
    target = Path(str(out_tsv) + ".tmp") if overwriting else out_tsv

    n_filled = 0
    n_total = 0
    with open(in_tsv, "r", encoding="utf-8", newline="") as fi:
        reader = csv.DictReader(fi, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for col in ("ACMG_POINTS", "ACMG_EVIDENCE", "ACMG_CLASS"):
            if col not in fieldnames:
                fieldnames.append(col)
        with open(target, "w", encoding="utf-8", newline="") as fo:
            writer = csv.DictWriter(fo, fieldnames=fieldnames, delimiter="\t",
                                    extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                n_total += 1
                key = (
                    (row.get("CHROM") or "").strip(),
                    (row.get("POS")   or "").strip(),
                    (row.get("REF")   or "").strip(),
                    (row.get("ALT")   or "").strip(),
                )
                gbt = gb.get(key)
                if gbt:
                    score, crit, cls = gbt
                    changed = False
                    if not (row.get("ACMG_POINTS") or "").strip() and score:
                        row["ACMG_POINTS"] = score
                        changed = True
                    if not (row.get("ACMG_EVIDENCE") or "").strip() and crit:
                        row["ACMG_EVIDENCE"] = crit
                        changed = True
                    if not (row.get("ACMG_CLASS") or "").strip() and cls:
                        row["ACMG_CLASS"] = cls
                        changed = True
                    if changed:
                        n_filled += 1
                writer.writerow(row)
    if overwriting:
        os.replace(target, out_tsv)
    return n_filled, n_total


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
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--tsv", required=True,
                    help="snv_indel.annotated.tsv (updated in place unless "
                         "--out-tsv given)")
    ap.add_argument("--out-tsv",
                    help="write merged TSV here instead of overwriting --tsv")
    ap.add_argument("--workdir",
                    help="keep intermediate sites.vcf / sites.genebe.vcf here "
                         "(default: temp dir, removed on success)")
    ap.add_argument("--username", default=os.environ.get("GENEBE_USER"))
    ap.add_argument("--api_key",  default=os.environ.get("GENEBE_API_KEY"))
    ap.add_argument("--sif",
                    default=os.environ.get(
                        "GENEBE_SIF",
                        str(Path.home() / "NGS_UI" / "biotools" / "genebe.sif")))
    ap.add_argument("--max-af", type=float, default=0.01,
                    help="drop sites whose GNOMAD_G_AF > this "
                         "(default 0.01; use -1 to disable)")
    ap.add_argument("--af-cols", default="GNOMAD_G_AF",
                    help="comma-separated AF columns to check "
                         "(default GNOMAD_G_AF — gnomAD genome 'global' AF, "
                         "matches filter_snv_tsv.py)")
    args = ap.parse_args()

    if not args.username or not args.api_key:
        print("ERROR: 需要 --username + --api_key（或 env GENEBE_USER/GENEBE_API_KEY）",
              file=sys.stderr)
        return 2

    in_tsv = Path(args.tsv).resolve()
    if not in_tsv.is_file():
        print(f"ERROR: --tsv 找不到：{in_tsv}", file=sys.stderr)
        return 2
    out_tsv = Path(args.out_tsv).resolve() if args.out_tsv else in_tsv

    if args.workdir:
        wd = Path(args.workdir)
        wd.mkdir(parents=True, exist_ok=True)
        wd_ctx = None
    else:
        wd_ctx = tempfile.TemporaryDirectory(prefix="genebe-")
        wd = Path(wd_ctx.name)

    try:
        sites = wd / "sites.vcf"
        annot = wd / "sites.genebe.vcf"

        max_af = None if args.max_af < 0 else args.max_af
        af_cols = [c.strip() for c in args.af_cols.split(",") if c.strip()]
        n_sites, n_af, n_star = tsv_to_sites(
            in_tsv, sites, max_af=max_af, af_cols=af_cols,
        )
        af_note = (f"dropped {n_af} above AF {max_af}"
                   if max_af is not None else "AF filter off")
        print(f"[genebe] {n_sites} unique sites → sites.vcf  "
              f"({af_note}; skipped {n_star} with '*')", file=sys.stderr)
        if n_sites == 0:
            print("ERROR: 0 sites — TSV 是否有 CHROM/POS/REF/ALT，或 --max-af 過嚴？",
                  file=sys.stderr)
            return 1

        ensure_sif(Path(args.sif))
        run_genebe(Path(args.sif), sites, annot, args.username, args.api_key)
        gb = parse_annotated_vcf(annot)
        n_filled, n_total = merge_into_tsv(in_tsv, out_tsv, gb)
        print(f"[genebe] backfilled ACMG for {n_filled}/{n_total} TSV rows "
              f"(GeneBe annotated {len(gb)} sites)", file=sys.stderr)
        print(f"[genebe] done → {out_tsv}", file=sys.stderr)
    finally:
        if wd_ctx is not None:
            wd_ctx.cleanup()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
