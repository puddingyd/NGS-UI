#!/usr/bin/env python3
"""Annotate review-filtered SNV/Indel rows from the local LitVar2 SQLite.

The complete pipeline TSV remains immutable. ``run_stopgaps.sh`` invokes this
script on its disposable working copy after MANE mapping and before the sparse
overlay is built. Candidates use the exact same WES/WGS filter as
``snv_indel.review.tsv``; transcript rows for one genomic variant are queried
once and receive the same literature payload.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app.config import LITVAR2_DB  # noqa: E402
from app.services import litvar2_store, snv_review  # noqa: E402


OUT_FIELDS = (
    "LITVAR2_ID",
    "LITVAR2_RSID",
    "LITVAR2_PMID_COUNT",
    "LITVAR2_PMIDS_TOP5",
    "LITVAR2_DATASET_DATE",
    "LITVAR2_MATCH_METHOD",
    "LITVAR2_STATUS",
    "LITVAR2_URL",
)
_RSID_RE = re.compile(r"(?i)\brs\d+\b")


def _variant_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    chrom = str(row.get("CHROM") or "").strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return (
        chrom.upper(),
        str(row.get("POS") or "").strip(),
        str(row.get("REF") or "").strip().upper(),
        str(row.get("ALT") or "").strip().upper(),
    )


def _append_unique(values: list[str], raw: str) -> None:
    value = str(raw or "").strip()
    if value and value not in values:
        values.append(value)


def _candidate_groups(
    tsv: Path,
    *,
    test_type: str,
) -> tuple[dict[tuple[str, str, str, str], dict[str, list[str]]], int]:
    bed = snv_review.load_candidate_bed()
    groups: dict[tuple[str, str, str, str], dict[str, list[str]]] = {}
    candidate_rows = 0
    with tsv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if not snv_review.is_review_candidate(
                row,
                test_type=test_type,
                bed=bed,
            ):
                continue
            candidate_rows += 1
            group = groups.setdefault(
                _variant_key(row),
                {"rsids": [], "genes": [], "hgvs": []},
            )
            for rsid in _RSID_RE.findall(str(row.get("RS_ID") or "")):
                _append_unique(group["rsids"], rsid)
            _append_unique(group["genes"], row.get("GENE") or "")
            for field in ("HGVS_C", "HGVS_P"):
                _append_unique(group["hgvs"], row.get(field) or "")
    return groups, candidate_rows


def _payload(result: dict[str, object]) -> dict[str, str]:
    return {
        "LITVAR2_ID": str(result.get("litvar_id") or ""),
        "LITVAR2_RSID": str(result.get("rsid") or ""),
        "LITVAR2_PMID_COUNT": str(result.get("pmids_count") or 0),
        "LITVAR2_PMIDS_TOP5": ",".join(
            str(value) for value in (result.get("pmids") or [])
            if str(value).isdigit()
        ),
        "LITVAR2_DATASET_DATE": str(result.get("dataset_date") or ""),
        "LITVAR2_MATCH_METHOD": str(result.get("match_method") or ""),
        "LITVAR2_STATUS": str(result.get("status") or ""),
        "LITVAR2_URL": str(result.get("url") or ""),
    }


def annotate_tsv(
    tsv: Path,
    db_path: Path,
    *,
    test_type: str,
) -> dict[str, int]:
    tsv = Path(tsv)
    groups, candidate_rows = _candidate_groups(tsv, test_type=test_type)
    results: dict[tuple[str, str, str, str], dict[str, str]] = {}
    counts = {"hit": 0, "no_match": 0, "ambiguous": 0}
    with litvar2_store.open_readonly(db_path) as conn:
        for key, identifiers in groups.items():
            result = litvar2_store.lookup_variant(
                conn,
                rsids=identifiers["rsids"],
                genes=identifiers["genes"],
                hgvs_values=identifiers["hgvs"],
            )
            status = str(result.get("status") or "no_match")
            if status in counts:
                counts[status] += 1
            results[key] = _payload(result)

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
                    payload = results.get(_variant_key(row))
                    for field in OUT_FIELDS:
                        row[field] = payload.get(field, "") if payload else ""
                    writer.writerow(row)
            os.replace(tmp_name, tsv)
        except Exception:
            try:
                os.unlink(tmp_name)
            except OSError:
                pass
            raise
    return {
        "candidate_rows": candidate_rows,
        "candidate_variants": len(groups),
        **counts,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tsv", required=True, type=Path)
    parser.add_argument("--db", type=Path, default=LITVAR2_DB)
    parser.add_argument("--test-type", choices=["WES", "WGS"], default="WES")
    args = parser.parse_args()
    if not args.db.is_file():
        print(f"[litvar2] local database not found, skipped: {args.db}")
        return 0
    try:
        stats = annotate_tsv(
            args.tsv,
            args.db,
            test_type=args.test_type,
        )
    except Exception as exc:
        # Literature lookup is reviewer context, never a release gate.
        print(
            f"[litvar2] WARNING: local annotation failed, skipped: {exc}",
            file=sys.stderr,
        )
        return 0
    print(
        "[litvar2] "
        + " ".join(f"{key}={value}" for key, value in stats.items())
        + f" db={args.db}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
