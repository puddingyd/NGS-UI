#!/usr/bin/env python3
"""Freeze legacy SNV OMIM disease ticks into content-keyed picks.

Legacy reviews stored ``edits[vid]["report_diseases"] = {"2": True}``, i.e.
a Disease-slot number that every consumer re-read at export time. This
script resolves those slots once against the OMIM workbook currently in use
(exactly what the platform shows today) and writes
``edits[vid]["report_disease_picks"]`` next to them, so later workbook
reorders or additions can no longer swap the reported disease.

The legacy field is left untouched (rollback-safe); readers prefer the new
field whenever it exists. Variants already carrying the new field are
skipped, so the script is idempotent.

Usage (from the repo root, preferably while nobody is editing reports):

    PYTHONPATH=backend python scripts/migrate_report_disease_picks.py          # dry run
    PYTHONPATH=backend python scripts/migrate_report_disease_picks.py --apply
    PYTHONPATH=backend python scripts/migrate_report_disease_picks.py --apply SAMPLE ...
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app.services import (  # noqa: E402
    disease_picks,
    report_store,
    sample_layout,
    sample_loader,
)


def _legacy_slots(edit: dict) -> list[int]:
    return disease_picks._legacy_slots(edit)


def _pending(meta: dict) -> dict[str, list[int]]:
    edits = meta.get("edits") if isinstance(meta.get("edits"), dict) else {}
    out = {}
    for vid, edit in edits.items():
        if not isinstance(edit, dict) or disease_picks.has_picks(edit):
            continue
        slots = _legacy_slots(edit)
        if slots:
            out[vid] = slots
    return out


def _write_atomic(path: Path, data: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    os.replace(tmp, path)


def migrate_sample(sample_id: str, *, apply: bool) -> tuple[int, int, list[str]]:
    """Return (converted, skipped, messages) for one sample."""
    meta_path = report_store._meta_path(sample_id)
    meta = report_store._read_json(meta_path)
    pending = _pending(meta)
    if not pending:
        return 0, 0, []
    payload = sample_loader.load_sample(sample_id, include_aux=False) or {}
    variants = payload.get("variants") or {}
    converted: dict[str, dict] = {}
    messages: list[str] = []
    for vid, slots in pending.items():
        variant = variants.get(vid)
        if not variant:
            messages.append(f"  skip {vid}: variant not in sample payload (legacy slots kept)")
            continue
        picks = disease_picks.picks_from_slots(variant, slots)
        if not picks:
            messages.append(f"  skip {vid}: slots {slots} are empty in the current workbook")
            continue
        converted[vid] = picks
        labels = "; ".join(p["label"] for p in picks.values())
        messages.append(f"  {vid} slots {slots} -> {labels}")
    if apply and converted:
        with report_store._WRITE_LOCK:
            fresh = report_store._read_json(meta_path)
            edits = dict(fresh.get("edits") or {})
            for vid, picks in converted.items():
                edit = dict(edits.get(vid) or {})
                if disease_picks.has_picks(edit):
                    continue  # reviewer touched it meanwhile
                edit[disease_picks.PICKS_FIELD] = picks
                edits[vid] = edit
            fresh["edits"] = edits
            _write_atomic(meta_path, fresh)
    return len(converted), len(pending) - len(converted), messages


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("samples", nargs="*", help="sample IDs (default: all)")
    ap.add_argument("--apply", action="store_true", help="write changes")
    args = ap.parse_args(argv)
    sample_ids = args.samples or list(sample_layout.iter_sample_ids())
    total_conv = total_skip = 0
    for sample_id in sample_ids:
        try:
            conv, skip, messages = migrate_sample(sample_id, apply=args.apply)
        except Exception as exc:  # keep going; report the failure
            print(f"{sample_id}: ERROR {exc}")
            continue
        if conv or skip:
            print(f"{sample_id}: {conv} converted, {skip} skipped")
            for line in messages:
                print(line)
        total_conv += conv
        total_skip += skip
    mode = "applied" if args.apply else "dry run (use --apply to write)"
    print(f"\n{total_conv} variant(s) converted, {total_skip} skipped — {mode}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
