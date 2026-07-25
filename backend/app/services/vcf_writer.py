"""Generate a minimal VCF from the immutable 03_acmg TSV for Exomiser/LIRICAL.

Both tools look up gnomAD AFs / pathogenicity scores from their own
databases by genomic coordinates, so a VCF carrying just CHROM /
POS / REF / ALT / GT is enough to drive them. Trade-off: variants
filtered out before the TSV (e.g. AF≥0.05) won't show up here, so
gene-level scoring loses some compound-het fidelity. For the
post-pipeline tertiary stage, that's the lesser of two evils
compared with maintaining a separate VCF path per sample.

Layout v3 prefixes the filename with the LIS ID; layout-v2 names remain
readable through :mod:`sample_layout`.
"""
from __future__ import annotations

import csv
import gzip
import json
from datetime import datetime, timezone
from pathlib import Path

from . import sample_layout
from .snv_rows import is_reportable_raw_row


VCF_FILENAME = "vcf_from_tsv.vcf.gz"
VCF_META_FILENAME = "vcf_from_tsv.vcf.gz.source.json"
WRITER_VERSION = 2

# UCSC-style names; both hg19 and hg38 TSVs in this codebase use them.
CONTIGS = [f"chr{n}" for n in range(1, 23)] + ["chrX", "chrY", "chrM"]


def vcf_path_for(lis_id: str, *, for_write: bool = False) -> Path:
    return sample_layout.state_file(lis_id, VCF_FILENAME, for_write=for_write)


def _meta_path_for(lis_id: str, *, for_write: bool = False) -> Path:
    return sample_layout.state_file(
        lis_id,
        VCF_META_FILENAME,
        for_write=for_write,
    )


def _tsv_signature(path: Path) -> dict:
    st = path.stat()
    return {
        "path": str(path),
        "mtime_ns": st.st_mtime_ns,
        "size": st.st_size,
    }


def needs_rebuild(lis_id: str) -> bool:
    """True if the VCF is missing or older than the source TSV.

    Used by the worker to refresh stale VCFs before invoking
    Exomiser/LIRICAL. Fresh registers don't need to call this — they
    just call from_tsv() unconditionally.
    """
    out = vcf_path_for(lis_id)
    if not out.exists():
        return True
    tsv = sample_layout.snv_raw_tsv(lis_id)
    if not tsv.exists():
        return False
    meta_path = _meta_path_for(lis_id)
    if not meta_path.exists():
        return True
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8")) or {}
    except (OSError, json.JSONDecodeError):
        return True
    if meta.get("writer_version") != WRITER_VERSION:
        return True
    return meta.get("tsv") != _tsv_signature(tsv)


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
    sample_dir = sample_layout.state_dir(lis_id)
    tsv = sample_layout.snv_raw_tsv(lis_id)
    if not tsv.is_file():
        raise FileNotFoundError(f"03_acmg SNV TSV missing for {lis_id}")
    out = vcf_path_for(lis_id, for_write=True)
    out.parent.mkdir(parents=True, exist_ok=True)

    rows_by_key: dict[tuple[str, int, str, str], str] = {}
    source_rows = 0
    with tsv.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            source_rows += 1
            if not is_reportable_raw_row(r):
                continue
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
            rows_by_key.setdefault((chrom, pos, ref, alt), gt)

    rows = sorted(rows_by_key.items(), key=lambda x: (_chrom_sort(x[0][0]), x[0][1], x[0][2], x[0][3]))

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
        for (chrom, pos, ref, alt), gt in rows:
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")

    _meta_path_for(lis_id, for_write=True).write_text(
        json.dumps(
            {
                "writer_version": WRITER_VERSION,
                "tsv": _tsv_signature(tsv),
                "records": len(rows),
                "source_rows": source_rows,
                "written_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    return out
