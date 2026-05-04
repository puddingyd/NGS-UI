"""Per-sample report state (status / edits / panels / comment / tags).

Stored as `data/reports/{sample_id}.json`. Phase 4 will swap this for a
DB-backed store with per-user audit trail; for now a flat JSON file
preserves the old UI's save behaviour.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ..config import REPORTS_DIR


_DEFAULT = {
    "status": {},
    "edits": {},
    "panels": {},
    "tags": [],
    "manual_variants": [],
    "comment": "",
    "clinical_description": "",
    "category": None,
    "yield": 0,
    "updated_at": None,
}


def _path(sample_id: str) -> Path:
    return REPORTS_DIR / f"{sample_id}.json"


def load(sample_id: str) -> dict:
    p = _path(sample_id)
    if not p.exists():
        return dict(_DEFAULT)
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return dict(_DEFAULT)
    out = dict(_DEFAULT)
    out.update(data if isinstance(data, dict) else {})
    return out


def save(sample_id: str, payload: dict) -> dict:
    out = dict(_DEFAULT)
    out.update(payload or {})
    out["updated_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    _path(sample_id).write_text(
        json.dumps(out, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return out
