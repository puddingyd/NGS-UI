#!/usr/bin/env python3
"""backfill_patient_list_uploads.py — reconstruct uploads.json history
from the archived xlsx files in patient_list/.

Why: the upload-history log only started recording after the
patient-list modal was added. Earlier weekly uploads still have their
raw xlsx archived under
    $NGS_UI_HOME/patient_list/{YYYYMMDDTHHMMSSZ}_{name}.xlsx
but no record in uploads.json — so the modal's "歷次上傳記錄" table
shows nothing from before that change.

What this does: replay every archived xlsx in chronological order
(filename timestamp prefix → primary, mtime → fallback), simulate the
merge into a fresh roster to compute per-upload parsed/added/updated/
total_after counts that match what would have been logged at the
time, and merge the result into uploads.json. Records already present
(matched by archive_name) are skipped — this is idempotent.

The on-disk roster.json is NEVER touched; this is read-only on the
archives + write-only on uploads.json.

Usage:
    PYTHONPATH=backend scripts/backfill_patient_list_uploads.py            # dry-run
    PYTHONPATH=backend scripts/backfill_patient_list_uploads.py --apply    # write
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

# Reuse the production parser + paths so the row count matches whatever
# the running service computes, even if the parser evolves.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "backend"))
from app.config import PATIENT_LIST_DIR  # noqa: E402
from app.services import patient_list_store as pls  # noqa: E402


_TS_FMT = "%Y%m%dT%H%M%SZ"


def _ts_from_filename(name: str) -> str | None:
    """Archive filenames start with the ingest timestamp (UTC, ISO-ish
    compact). Extract → return ISO-8601 with offset, or None when the
    prefix doesn't parse.
    """
    head = name.split("_", 1)[0]
    try:
        dt = datetime.strptime(head, _TS_FMT).replace(tzinfo=timezone.utc)
    except ValueError:
        return None
    return dt.isoformat(timespec="seconds")


def _ts_from_mtime(path: Path) -> str:
    dt = datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc)
    return dt.isoformat(timespec="seconds")


def _original_filename(archive_name: str) -> str:
    """`20251128T032015Z_未完成報告清單.xlsx` → `未完成報告清單.xlsx`."""
    parts = archive_name.split("_", 1)
    return parts[1] if len(parts) == 2 else archive_name


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true",
                    help="actually write to uploads.json (default: dry-run)")
    args = ap.parse_args()

    archives = sorted(p for p in PATIENT_LIST_DIR.glob("*.xlsx") if p.is_file())
    if not archives:
        print(f"(no xlsx archives under {PATIENT_LIST_DIR})")
        return 0

    existing = pls.list_uploads()  # already-recorded; we'll skip these
    known = {r.get("archive_name") for r in existing if r.get("archive_name")}
    print(f"found {len(archives)} archived xlsx; {len(known)} already in uploads.json")

    # Replay every archive in chronological order to recompute the
    # added/updated/total_after counts each upload would have produced.
    # Build the synthetic roster fresh — we never touch the real one.
    synth: dict[str, dict] = {}
    new_records: list[dict] = []
    for path in archives:
        ts = _ts_from_filename(path.name) or _ts_from_mtime(path)
        try:
            parsed = pls.parse_xlsx(path)
        except Exception as exc:  # noqa: BLE001
            print(f"  ! {path.name}: parse failed ({exc}) — skipping")
            continue

        added = updated = 0
        for rec in parsed:
            lid = rec["lis_id"]
            prev = synth.get(lid)
            entry = {
                "lis_id":     lid,
                "specimen":   rec["specimen"],
                "mrn":        rec["mrn"],
                "name":       rec["name"],
                "test_name":  rec["test_name"],
                "test_type":  rec["test_type"],
                "department": rec["department"],
            }
            if prev is None:
                synth[lid] = entry
                added += 1
            else:
                changed = any(
                    prev.get(k) != entry.get(k)
                    for k in ("specimen", "mrn", "name",
                              "test_name", "test_type", "department")
                )
                synth[lid] = entry
                if changed:
                    updated += 1

        if path.name in known:
            print(f"  · {path.name}  (already logged — skipping write)")
            continue

        rec = {
            "uploaded_at":       ts,
            "original_filename": _original_filename(path.name),
            "archive_name":      path.name,
            "parsed":            len(parsed),
            "added":             added,
            "updated":           updated,
            "total_after":       len(synth),
        }
        new_records.append(rec)
        print(f"  + {path.name}  parsed={len(parsed)}  added={added}  "
              f"updated={updated}  total_after={len(synth)}")

    if not new_records:
        print("nothing new to backfill — uploads.json already covers all archives.")
        return 0

    print(f"\nwould add {len(new_records)} record(s) to uploads.json")
    if not args.apply:
        print("(dry-run; pass --apply to write)")
        return 0

    merged = existing + new_records
    out = pls._uploads_path()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(merged, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
