"""Correct NCKUH diploid chrX calls in the disposable SNV working TSV."""
from __future__ import annotations

import shutil
from pathlib import Path

from .ploidy import read_karyotype


# GRCh38 pseudoautosomal intervals, 1-based and inclusive.
_X_PAR = ((10_001, 2_781_479), (155_701_383, 156_030_895))


def _is_x_nonpar(chrom: str, pos: str, ref: str) -> bool:
    if chrom.strip().lower() not in {"x", "chrx"}:
        return False
    try:
        start = int(pos)
    except (TypeError, ValueError):
        return False
    if start < 1 or not ref or ref == ".":
        return False
    end = start + len(ref) - 1
    return not any(start <= par_end and end >= par_start for par_start, par_end in _X_PAR)


def copy_nckuh_work_tsv(
    raw_tsv: Path,
    work_tsv: Path,
    ploidy_vcf: Path | None,
) -> int:
    """Copy raw TSV, changing XY non-PAR chrX hom to hemi; return row count.

    The ploidy sidecar must be the exact source-sample match selected by the
    tertiary worker. Unknown/aneuploid karyotypes leave the TSV untouched.
    GT and HAPLOID_HET stay as recorded by the upstream callers.
    """
    raw_tsv = Path(raw_tsv)
    work_tsv = Path(work_tsv)
    if raw_tsv.resolve() == work_tsv.resolve():
        raise ValueError("NCKUH SNV working TSV must differ from the immutable raw TSV")
    if ploidy_vcf is None or read_karyotype(ploidy_vcf) != "XY":
        shutil.copyfile(raw_tsv, work_tsv)
        return 0

    corrected = 0
    with raw_tsv.open("r", encoding="utf-8", newline="") as source, \
            work_tsv.open("w", encoding="utf-8", newline="") as target:
        header = source.readline()
        target.write(header)
        columns = header.rstrip("\r\n").split("\t")
        required = ("CHROM", "POS", "REF", "ZYGOSITY", "HAPLOID_HET")
        missing = [name for name in required if name not in columns]
        if missing:
            raise ValueError(f"NCKUH SNV TSV missing columns: {', '.join(missing)}")
        indices = {name: columns.index(name) for name in required}
        for line in source:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(columns):
                raise ValueError("NCKUH SNV TSV row has a different column count")
            if (
                fields[indices["ZYGOSITY"]].strip().lower() in {"hom", "homozygous"}
                and fields[indices["HAPLOID_HET"]].strip() in {"", "."}
                and _is_x_nonpar(
                    fields[indices["CHROM"]],
                    fields[indices["POS"]],
                    fields[indices["REF"]],
                )
            ):
                fields[indices["ZYGOSITY"]] = "hemi"
                newline = "\r\n" if line.endswith("\r\n") else "\n" if line.endswith("\n") else ""
                target.write("\t".join(fields) + newline)
                corrected += 1
            else:
                target.write(line)
    return corrected
