"""Build the compact SNV/Indel TSV used by the main review screen.

The pipeline-owned snv_indel.annotated.tsv remains the complete source of
truth.  This derived file is intentionally disposable: when the raw TSV
changes, the next sample load rebuilds it atomically.
"""
from __future__ import annotations

import csv
import json
import os
import re
from pathlib import Path

REVIEW_TSV_NAME = "snv_indel.review.tsv"
MAX_GNOMAD_G_AF = 0.05

_PATHOGENIC_RE = re.compile(r"(?:^|[/|,; ])(?:likely[_ ]?)?pathogenic(?:$|[/|,; ])", re.I)


def _to_float(raw: str) -> float | None:
    value = (raw or "").strip()
    if not value or value in {".", "NA", "N/A"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _keep_row(row: dict[str, str]) -> bool:
    """Keep rare / unknown-AF rows and ClinVar P/LP rescue rows."""
    sig = (row.get("CLINVAR_SIG") or "").strip()
    if _PATHOGENIC_RE.search(sig):
        return True
    af = _to_float(row.get("GNOMAD_G_AF") or "")
    return af is None or af < MAX_GNOMAD_G_AF


def ensure_review_tsv(raw_tsv: Path, *, keep_ids: set[str] | None = None) -> Path:
    """Return an up-to-date compact review TSV derived from *raw_tsv*."""
    keep_ids = keep_ids or set()
    review_tsv = raw_tsv.with_name(REVIEW_TSV_NAME)
    manifest = review_tsv.with_suffix(review_tsv.suffix + ".source.json")
    source = {
        "raw_mtime_ns": raw_tsv.stat().st_mtime_ns,
        "raw_size": raw_tsv.stat().st_size,
        "keep_ids": sorted(keep_ids),
        "max_gnomad_g_af": MAX_GNOMAD_G_AF,
    }
    if review_tsv.is_file() and manifest.is_file():
        try:
            if json.loads(manifest.read_text(encoding="utf-8")) == source:
                return review_tsv
        except (OSError, json.JSONDecodeError):
            pass

    tmp = review_tsv.with_suffix(review_tsv.suffix + ".tmp")
    manifest_tmp = manifest.with_suffix(manifest.suffix + ".tmp")
    with raw_tsv.open("r", encoding="utf-8", newline="") as src:
        reader = csv.DictReader(src, delimiter="\t")
        fieldnames = reader.fieldnames or []
        with tmp.open("w", encoding="utf-8", newline="") as dst:
            writer = csv.DictWriter(
                dst, fieldnames=fieldnames, delimiter="\t",
                extrasaction="ignore", lineterminator="\n",
            )
            writer.writeheader()
            for row in reader:
                vid = "-".join(
                    (row.get(k) or "").strip() for k in ("CHROM", "POS", "REF", "ALT")
                )
                if vid in keep_ids or _keep_row(row):
                    writer.writerow(row)
    os.replace(tmp, review_tsv)
    manifest_tmp.write_text(
        json.dumps(source, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    os.replace(manifest_tmp, manifest)
    return review_tsv
