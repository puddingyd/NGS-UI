"""Shared raw SNV row eligibility helpers."""
from __future__ import annotations


PRIMARY_CONTIGS = {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}


def normalized_contig(raw: str) -> str:
    value = str(raw or "").strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value.upper()


def is_primary_contig(raw: str) -> bool:
    return normalized_contig(raw) in PRIMARY_CONTIGS


def is_reportable_raw_row(row: dict[str, str]) -> bool:
    """Match the old in-place filter without mutating the pipeline TSV."""
    alt = str(row.get("ALT") or "").strip()
    return bool(alt and alt != "*" and is_primary_contig(row.get("CHROM") or ""))

