"""Generate a minimal VCF from snv_indel.annotated.tsv for Exomiser/LIRICAL.

Both tools look up gnomAD AFs / pathogenicity scores from their own
databases by genomic coordinates, so a VCF carrying just CHROM /
POS / REF / ALT / GT is enough to drive them. Trade-off: variants
filtered out before the TSV (e.g. AF≥0.05) won't show up here, so
gene-level scoring loses some compound-het fidelity. For the
post-pipeline tertiary stage, that's the lesser of two evils
compared with maintaining a separate VCF path per sample.

Output filename convention:
    tertiary_output/{LIS_ID}/{LIS_ID}.from_tsv.vcf.gz
"""
from __future__ import annotations

import csv
import gzip
from datetime import datetime, timezone
from pathlib import Path

from ..config import TERTIARY_OUTPUT_ROOT


VCF_FILENAME_SUFFIX = ".from_tsv.vcf.gz"

# UCSC-style names; both hg19 and hg38 TSVs in this codebase use them.
CONTIGS = [f"chr{n}" for n in range(1, 23)] + ["chrX", "chrY", "chrM"]


def vcf_path_for(lis_id: str) -> Path:
    return TERTIARY_OUTPUT_ROOT / lis_id / f"{lis_id}{VCF_FILENAME_SUFFIX}"


def needs_rebuild(lis_id: str) -> bool:
    """True if the VCF is missing or older than the source TSV.

    Used by the worker to refresh stale VCFs before invoking
    Exomiser/LIRICAL. Fresh registers don't need to call this — they
    just call from_tsv() unconditionally.
    """
    out = vcf_path_for(lis_id)
    if not out.exists():
        return True
    tsv = TERTIARY_OUTPUT_ROOT / lis_id / "snv_indel.annotated.tsv"
    if not tsv.exists():
        return False
    return out.stat().st_mtime < tsv.stat().st_mtime


def _pick_gt(row: dict) -> str:
    """GT_DV preferred (DeepVariant tends to be cleaner on most loci);
    fall back to GT_HC. Both ./. → skip the variant entirely."""
    gt_dv = (row.get("GT_DV") or "").strip()
    gt_hc = (row.get("GT_HC") or "").strip()
    if gt_dv and gt_dv != "./.":
        return gt_dv
    if gt_hc and gt_hc != "./.":
        return gt_hc
    return ""


def _chrom_sort(chrom: str) -> int:
    c = chrom.replace("chr", "").upper()
    if c == "X":  return 23
    if c == "Y":  return 24
    if c in ("M", "MT"): return 25
    try:
        return int(c)
    except ValueError:
        return 99


def from_tsv(lis_id: str) -> Path:
    """Read the sample's TSV and write a minimal gzipped VCF beside it.

    Returns the output path. Raises FileNotFoundError if the TSV is
    missing.
    """
    sample_dir = TERTIARY_OUTPUT_ROOT / lis_id
    tsv = sample_dir / "snv_indel.annotated.tsv"
    if not tsv.is_file():
        raise FileNotFoundError(f"snv_indel.annotated.tsv missing for {lis_id}")
    out = vcf_path_for(lis_id)
    out.parent.mkdir(parents=True, exist_ok=True)

    rows: list[tuple[str, int, str, str, str]] = []
    with tsv.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            chrom = (r.get("CHROM") or "").strip()
            pos_s = (r.get("POS")   or "").strip()
            ref   = (r.get("REF")   or "").strip()
            alt   = (r.get("ALT")   or "").strip()
            if not all([chrom, pos_s, ref, alt]):
                continue
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            gt = _pick_gt(r)
            if not gt:
                continue
            rows.append((chrom, pos, ref, alt, gt))

    rows.sort(key=lambda x: (_chrom_sort(x[0]), x[1]))

    today = datetime.now(timezone.utc).strftime("%Y%m%d")
    with gzip.open(out, "wt", encoding="utf-8", newline="\n") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##fileDate={today}\n")
        f.write(f"##source=NGS-UI/vcf_writer.from_tsv\n")
        for c in CONTIGS:
            f.write(f"##contig=<ID={c}>\n")
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{lis_id}\n")
        for chrom, pos, ref, alt, gt in rows:
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")

    return out
