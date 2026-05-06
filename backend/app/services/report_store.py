"""Reviewer-side state (status / edits / panels / comment / tags …).

Post-migration, these fields all live on the per-patient
sample_metadata.json so they survive across analysis versions. Pre
migration, the loader falls back to data/reports/{sample_id}.json so a
deploy doesn't break the UI before the migration script runs.

Phase 4 will swap this for a DB-backed store with per-user audit trail.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from ..config import REPORTS_DIR, TERTIARY_OUTPUT_ROOT


# Reviewer-editable fields that live on sample_metadata.json. Fields
# outside this list (lis_id, name, vcf_path, …) are pipeline metadata
# and never get clobbered by a save() call.
_REVIEWER_FIELDS = {
    "status",
    "edits",
    "panels",
    "manual_variants",
    "tags",
    "comment",
    "clinical_description",
    "category",
    "yield",
}

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


def _meta_path(sample_id: str) -> Path:
    return TERTIARY_OUTPUT_ROOT / sample_id / "sample_metadata.json"


def _legacy_path(sample_id: str) -> Path:
    return REPORTS_DIR / f"{sample_id}.json"


def _read_json(p: Path) -> dict:
    if not p.exists():
        return {}
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except json.JSONDecodeError:
        return {}


def _write_json(p: Path, data: dict) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(
        json.dumps(data, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def _project_reviewer(meta: dict) -> dict:
    """Pull reviewer-only fields out of sample_metadata.json."""
    out = dict(_DEFAULT)
    for k in _REVIEWER_FIELDS:
        if k in meta:
            out[k] = meta[k]
    out["updated_at"] = meta.get("updated_at")
    return out


def load(sample_id: str) -> dict:
    """Return the reviewer-state dict.

    Reads sample_metadata.json first (canonical post-migration). If it
    doesn't have any reviewer fields yet, falls back to the legacy
    data/reports/{sample_id}.json so old samples keep working.
    """
    meta = _read_json(_meta_path(sample_id))
    if any(k in meta for k in _REVIEWER_FIELDS):
        return _project_reviewer(meta)

    legacy = _read_json(_legacy_path(sample_id))
    if legacy:
        out = dict(_DEFAULT)
        out.update(legacy)
        return out
    return dict(_DEFAULT)


def save(sample_id: str, payload: dict) -> dict:
    """Merge reviewer fields into sample_metadata.json.

    Pipeline-owned fields (lis_id, name, mrn, test_type, vcf_path, …)
    are preserved as-is — the merge only touches keys in
    `_REVIEWER_FIELDS` plus `updated_at`.
    """
    p = _meta_path(sample_id)
    meta = _read_json(p)
    payload = payload or {}
    for k in _REVIEWER_FIELDS:
        if k in payload:
            meta[k] = payload[k]
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    meta["updated_at"] = now
    meta.setdefault("created_at", now)
    _write_json(p, meta)
    return _project_reviewer(meta)
