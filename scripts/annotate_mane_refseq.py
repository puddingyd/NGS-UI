#!/usr/bin/env python3
"""Annotate SNV/Indel TSV transcript rows with MANE RefSeq IDs.

The v3.5 tertiary pipeline emits Ensembl transcript IDs in TRANSCRIPT.
This post-processing step maps Ensembl_nuc from the MANE summary table to
RefSeq_nuc / RefSeq_prot so the UI and DOCX can display RefSeq first while
still preserving the original Ensembl transcript for auditability.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from pathlib import Path


OUT_FIELDS = ["REFSEQ_NUC", "REFSEQ_PROT", "MANE_STATUS"]


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def _tx_base(value: str) -> str:
    return str(value or "").strip().split(".", 1)[0].upper()


def _load_mane(path: Path) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            enst = (row.get("Ensembl_nuc") or "").strip()
            if not enst:
                continue
            data = {
                "REFSEQ_NUC": (row.get("RefSeq_nuc") or "").strip(),
                "REFSEQ_PROT": (row.get("RefSeq_prot") or "").strip(),
                "MANE_STATUS": (row.get("MANE_status") or "").strip(),
            }
            mapping[enst.upper()] = data
            mapping.setdefault(_tx_base(enst), data)
    return mapping


def annotate_tsv(tsv: Path, mane: Path) -> tuple[int, int]:
    mapping = _load_mane(mane)
    if not mapping:
        raise RuntimeError(f"MANE summary has no transcript mappings: {mane}")

    with tsv.open("r", encoding="utf-8", newline="") as src:
        reader = csv.DictReader(src, delimiter="\t")
        if not reader.fieldnames:
            raise RuntimeError(f"empty TSV: {tsv}")
        fieldnames = list(reader.fieldnames)
        for field in OUT_FIELDS:
            if field not in fieldnames:
                fieldnames.append(field)

        fd, tmp_name = tempfile.mkstemp(
            prefix=f".{tsv.name}.",
            suffix=".tmp",
            dir=str(tsv.parent),
            text=True,
        )
        total = matched = 0
        try:
            with os.fdopen(fd, "w", encoding="utf-8", newline="") as dst:
                writer = csv.DictWriter(
                    dst,
                    fieldnames=fieldnames,
                    delimiter="\t",
                    extrasaction="ignore",
                    lineterminator="\n",
                )
                writer.writeheader()
                for row in reader:
                    total += 1
                    tx = (row.get("TRANSCRIPT") or "").strip()
                    hit = mapping.get(tx.upper()) or mapping.get(_tx_base(tx))
                    if hit:
                        matched += 1
                        for field, value in hit.items():
                            row[field] = value
                    else:
                        for field in OUT_FIELDS:
                            row.setdefault(field, "")
                    writer.writerow(row)
            os.replace(tmp_name, tsv)
        except Exception:
            try:
                os.unlink(tmp_name)
            except OSError:
                pass
            raise
    return total, matched


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True, type=Path)
    ap.add_argument("--mane-summary", type=Path, default=None)
    args = ap.parse_args()

    mane = args.mane_summary
    if mane is None:
        ngs_home = Path(os.environ.get("NGS_UI_HOME", Path.home() / "NGS_UI"))
        biotools = Path(os.environ.get("NGS_UI_BIOTOOLS_DIR", ngs_home / "biotools"))
        mane = Path(os.environ.get("NGS_UI_MANE_SUMMARY", biotools / "MANE.GRCh38.v1.5.summary.txt.gz"))
    if not mane.is_file():
        print(f"[mane-refseq] MANE summary not found, skipped: {mane}")
        return 0
    total, matched = annotate_tsv(args.tsv, mane)
    print(f"[mane-refseq] rows={total} matched={matched} summary={mane}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
