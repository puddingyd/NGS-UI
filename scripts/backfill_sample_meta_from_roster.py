#!/usr/bin/env python3
"""backfill_sample_meta_from_roster.py — fill in department / physician /
sign_received_at on already-registered samples.

These three fields landed on sample_metadata.json after a bunch of
samples were already registered, so their cards display blank. The
fix: walk every registered sample, look up the roster by LIS_ID
(strip -dragen / -inhouse suffix), and copy over whichever of the
three fields is currently empty. Existing reviewer-edited values are
never overwritten.

Usage:
    PYTHONPATH=backend scripts/backfill_sample_meta_from_roster.py            # dry-run
    PYTHONPATH=backend scripts/backfill_sample_meta_from_roster.py --apply    # write
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "backend"))
from app.config import TERTIARY_OUTPUT_ROOT  # noqa: E402
from app.services import patient_list_store as pls  # noqa: E402


# DRAGEN / in-house mode workers register sample IDs with these
# suffixes so the same LIS_ID can be processed by multiple callers
# without dir collisions. The roster keys don't carry the suffix.
_SUFFIXES = ("-dragen", "-inhouse", "-WES", "-WGS")


def _strip_suffix(sid: str) -> str:
    for suf in _SUFFIXES:
        if sid.endswith(suf):
            return sid[: -len(suf)]
    return sid


def _candidate_lis_ids(meta_lis_id: str, dir_name: str) -> list[str]:
    """Look-up keys to try against the roster, most-specific first.

    Reviewer may have registered with a suffixed lis_id (e.g.
    VAL-58-dragen) or the base lis_id; the on-disk dir name may
    differ from meta.lis_id when the dir was renamed. Try both,
    plus the suffix-stripped version of each, so we always find
    the roster row if it exists.
    """
    seen: list[str] = []
    for candidate in (meta_lis_id, dir_name,
                      _strip_suffix(meta_lis_id), _strip_suffix(dir_name)):
        c = (candidate or "").strip()
        if c and c not in seen:
            seen.append(c)
    return seen


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true",
                    help="actually write to sample_metadata.json (default: dry-run)")
    args = ap.parse_args()

    roster = pls.load_roster()
    if not roster:
        print("(roster is empty — nothing to backfill from)")
        return 0

    root = TERTIARY_OUTPUT_ROOT
    if not root.is_dir():
        print(f"(no tertiary_output dir: {root})")
        return 0

    examined = filled = skipped_complete = no_roster = 0
    pending: list[tuple[Path, dict, dict]] = []  # (meta_path, meta, patch)

    for sub in sorted(root.iterdir()):
        if not sub.is_dir() or sub.name.startswith("_"):
            continue
        meta_path = sub / "sample_metadata.json"
        if not meta_path.is_file():
            continue
        examined += 1
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError) as e:
            print(f"  ! {sub.name}: meta read failed ({e})")
            continue
        if not isinstance(meta, dict):
            continue

        keys = _candidate_lis_ids(meta.get("lis_id", ""), sub.name)
        entry = None
        for k in keys:
            entry = roster.get(k)
            if entry:
                break
        if not entry:
            no_roster += 1
            continue

        patch: dict[str, str] = {}
        for field in ("department", "physician", "sign_received_at"):
            cur = (meta.get(field) or "").strip()
            src = (entry.get(field) or "").strip()
            if not cur and src:
                patch[field] = src
        if not patch:
            skipped_complete += 1
            continue

        pending.append((meta_path, meta, patch))
        filled += 1
        print(f"  + {sub.name}  ← roster[{k}]  "
              + "  ".join(f"{k}={v!r}" for k, v in patch.items()))

    print(f"\nexamined={examined}  to_fill={filled}  "
          f"already_complete={skipped_complete}  no_roster_match={no_roster}")

    if not pending:
        print("nothing to do.")
        return 0

    if not args.apply:
        print("(dry-run; pass --apply to write)")
        return 0

    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    for meta_path, meta, patch in pending:
        meta.update(patch)
        meta["metadata_updated_at"] = now
        meta_path.write_text(
            json.dumps(meta, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
    print(f"wrote {len(pending)} sample_metadata.json file(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
