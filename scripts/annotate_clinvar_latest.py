#!/usr/bin/env python3
"""Compare fixed pipeline ClinVar calls with the latest local snapshot.

The input is a disposable post-processing TSV. Fixed 2026-07-20 pipeline
``CLINVAR_*`` values remain untouched. Latest values are written separately as
``CLINVAR_LATEST_*`` fields only when the snapshot has a confident
allele/Variation ID match. Change arrows are limited to crossings into or out
of the P/LP bucket.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from datetime import date, datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app.services import clinvar_latest_store  # noqa: E402


BASE_FIELDS = {
    "CLINVAR_SIG": "CLINVAR_BASE_SIG",
    "CLINVAR_STARS": "CLINVAR_BASE_STARS",
    "CLINVAR_DN": "CLINVAR_BASE_DN",
    "CLINVAR_SIGCONF": "CLINVAR_BASE_SIGCONF",
    "CLINVAR_VARIATION_ID": "CLINVAR_BASE_VARIATION_ID",
}
LATEST_FIELDS = {
    "CLINVAR_SIG": "CLINVAR_LATEST_SIG",
    "CLINVAR_STARS": "CLINVAR_LATEST_STARS",
    "CLINVAR_DN": "CLINVAR_LATEST_DN",
    "CLINVAR_SIGCONF": "CLINVAR_LATEST_SIGCONF",
    "CLINVAR_VARIATION_ID": "CLINVAR_LATEST_VARIATION_ID",
}
OUTPUT_FIELDS = [
    *BASE_FIELDS.values(),
    *LATEST_FIELDS.values(),
    "CLINVAR_LATEST_REVIEW_STATUS",
    "CLINVAR_LATEST_APPLIED",
    "CLINVAR_CHANGE",
]


def _variation_id(row: dict[str, str]) -> str:
    return clinvar_latest_store.normalize_variation_id(
        row.get("CLINVAR_VARIATION_ID")
        or row.get("VARIATION_ID")
        or ""
    )


def _lookup_key(row: dict[str, str]) -> str:
    try:
        return clinvar_latest_store.variant_key(
            row.get("CHROM", ""),
            row.get("POS", ""),
            row.get("REF", ""),
            row.get("ALT", ""),
        )
    except (TypeError, ValueError):
        return ""


def _latest_values(record: dict[str, object]) -> dict[str, str]:
    return {
        "CLINVAR_SIG": str(record.get("significance") or ""),
        "CLINVAR_STARS": str(record.get("stars") if record.get("stars") is not None else ""),
        "CLINVAR_DN": str(record.get("disease") or ""),
        "CLINVAR_SIGCONF": str(record.get("significance_conflict") or ""),
        "CLINVAR_VARIATION_ID": str(record.get("variation_id") or ""),
    }


def annotate_tsv(
    tsv_path: Path,
    db_path: Path,
    marker_path: Path,
    *,
    baseline_release: str,
) -> dict[str, object]:
    tsv_path = Path(tsv_path)
    db_path = Path(db_path)
    db_meta = clinvar_latest_store.metadata(db_path)
    if db_meta.get("schema_version") != clinvar_latest_store.SCHEMA_VERSION:
        raise RuntimeError(
            "ClinVar database schema is stale; run update_clinvar_latest_db.py "
            f"(found {db_meta.get('schema_version') or 'missing'}, "
            f"expected {clinvar_latest_store.SCHEMA_VERSION})"
        )
    latest_release = db_meta.get("release_date", "")
    if not latest_release:
        raise RuntimeError(f"ClinVar database has no release_date metadata: {db_path}")
    try:
        baseline_date = date.fromisoformat(baseline_release)
        latest_date = date.fromisoformat(latest_release)
    except ValueError as exc:
        raise RuntimeError("ClinVar release dates must use YYYY-MM-DD") from exc
    if latest_date < baseline_date:
        marker = {
            "schema_version": 1,
            "status": "skipped_stale_database",
            "baseline_release": baseline_release,
            "latest_release": latest_release,
            "database": str(db_path),
            "annotated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "stats": {"rows": 0},
        }
        clinvar_latest_store.write_manifest(marker_path, marker)
        return marker

    stats = {
        "rows": 0,
        "matched": 0,
        "matched_by_variation_id": 0,
        "changed_annotation": 0,
        "up_to_plp": 0,
        "down_from_plp": 0,
        "verified_no_record": 0,
        "unmatched_without_id": 0,
    }

    target = Path(str(tsv_path) + ".clinvar.tmp")
    try:
        with tsv_path.open("r", encoding="utf-8", newline="") as source:
            reader = csv.DictReader(source, delimiter="\t")
            fields = list(reader.fieldnames or [])
            for field in OUTPUT_FIELDS:
                if field not in fields:
                    fields.append(field)
            with target.open("w", encoding="utf-8", newline="") as destination:
                writer = csv.DictWriter(
                    destination,
                    fieldnames=fields,
                    delimiter="\t",
                    extrasaction="ignore",
                    lineterminator="\n",
                )
                writer.writeheader()
                batch: list[dict[str, str]] = []
                for row in reader:
                    batch.append(dict(row))
                    if len(batch) >= 5_000:
                        _annotate_batch(batch, db_path, stats)
                        writer.writerows(batch)
                        batch.clear()
                if batch:
                    _annotate_batch(batch, db_path, stats)
                    writer.writerows(batch)
        os.replace(target, tsv_path)
    except BaseException:
        target.unlink(missing_ok=True)
        raise

    marker = {
        "schema_version": 1,
        "status": "complete",
        "baseline_release": baseline_release,
        "latest_release": latest_release,
        "database": str(db_path),
        "database_schema_version": db_meta.get("schema_version", ""),
        "database_record_count": int(db_meta.get("record_count") or 0),
        "annotated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "change_policy": {
            "up": "non-P/LP to P/LP",
            "down": "P/LP to non-P/LP",
            "ignored_examples": ["LP to P", "LB to B", "conflicting to VUS"],
        },
        "stats": stats,
    }
    clinvar_latest_store.write_manifest(marker_path, marker)
    return marker


def _annotate_batch(
    rows: list[dict[str, str]],
    db_path: Path,
    stats: dict[str, int],
) -> None:
    keys = [_lookup_key(row) for row in rows]
    variation_ids = [_variation_id(row) for row in rows]
    by_key, by_id = clinvar_latest_store.lookup_records(db_path, keys, variation_ids)
    stats["rows"] += len(rows)
    for row, key, variation_id in zip(rows, keys, variation_ids):
        record = by_key.get(key)
        if record is None and variation_id:
            record = by_id.get(variation_id)
            if record is not None:
                stats["matched_by_variation_id"] += 1
        old_values = {field: str(row.get(field) or "") for field in BASE_FIELDS}
        if record is not None:
            stats["matched"] += 1
            new_values = _latest_values(record)
        elif variation_id:
            # A previously recorded ClinVar variation that is absent by both
            # normalized allele and ID is meaningful as a current no-record.
            stats["verified_no_record"] += 1
            new_values = {field: "" for field in BASE_FIELDS}
        else:
            # Baseline no-record rows cannot distinguish "still absent" from a
            # representation mismatch. Preserve the original UI annotation.
            stats["unmatched_without_id"] += 1
            row["CLINVAR_LATEST_APPLIED"] = ""
            row["CLINVAR_CHANGE"] = ""
            continue

        if all(old_values[field] == new_values[field] for field in BASE_FIELDS):
            row["CLINVAR_LATEST_APPLIED"] = ""
            row["CLINVAR_CHANGE"] = ""
            continue
        stats["changed_annotation"] += 1
        row["CLINVAR_LATEST_APPLIED"] = "1"
        for live_field, baseline_field in BASE_FIELDS.items():
            row[baseline_field] = old_values[live_field]
            row[LATEST_FIELDS[live_field]] = new_values[live_field]
        row["CLINVAR_LATEST_REVIEW_STATUS"] = str(
            (record or {}).get("review_status") or ""
        )
        direction = clinvar_latest_store.meaningful_change(
            old_values["CLINVAR_SIG"], new_values["CLINVAR_SIG"]
        )
        row["CLINVAR_CHANGE"] = direction
        if direction == "UP_TO_PLP":
            stats["up_to_plp"] += 1
        elif direction == "DOWN_FROM_PLP":
            stats["down_from_plp"] += 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--marker", required=True)
    parser.add_argument("--baseline-release", default="2026-07-20")
    args = parser.parse_args()
    try:
        marker = annotate_tsv(
            Path(args.tsv),
            Path(args.db),
            Path(args.marker),
            baseline_release=args.baseline_release,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(marker["stats"], ensure_ascii=False), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
