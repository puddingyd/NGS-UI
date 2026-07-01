"""Reviewer-side state (status / edits / panels / comment / tags …).

These fields live on the per-patient sample_metadata.json so they
survive across analysis versions. Pipeline-owned keys (lis_id, name,
mrn, test_type, vcf_path, …) are preserved untouched on save() — only
the whitelist below gets overwritten.

Phase 4 will swap this for a DB-backed store with per-user audit trail.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from ..config import TERTIARY_OUTPUT_ROOT
from . import clinical_presentation_store


_REVIEWER_FIELDS = {
    "status",
    "edits",
    "panels",
    "secondary_findings",
    "manual_variants",
    "cnv_sv_merges",
    "tags",
    "comment",
    "clinical_description",
    "genetic_counseling",
    "category",
    "sry_confirmed",
    "yield",
}

_DEFAULT = {
    "status": {},
    "edits": {},
    "panels": {},
    "secondary_findings": {},
    "tags": [],
    "manual_variants": [],
    "cnv_sv_merges": [],
    "comment": "",
    "clinical_description": "",
    "genetic_counseling": "",
    "category": None,
    "sry_confirmed": False,
    "yield": 0,
    "updated_at": None,
}


def _meta_path(sample_id: str) -> Path:
    return TERTIARY_OUTPUT_ROOT / sample_id / "sample_metadata.json"


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
    meta = _read_json(_meta_path(sample_id))
    out = _project_reviewer(meta)
    if not str(out.get("clinical_description") or "").strip():
        try:
            sidecar = clinical_presentation_store.load(
                code=meta.get("lis_id") or meta.get("sample_id") or sample_id,
                mrn=meta.get("mrn") or "",
            )
        except ValueError:
            sidecar = {}
        if sidecar:
            out["clinical_description"] = (sidecar.get("content") or "").strip()
    return out


def save(sample_id: str, payload: dict) -> dict:
    """Merge reviewer fields into sample_metadata.json."""
    p = _meta_path(sample_id)
    meta = _read_json(p)
    payload = payload or {}
    for k in _REVIEWER_FIELDS:
        if k in payload:
            meta[k] = payload[k]
    if "clinical_description" in payload:
        code = meta.get("lis_id") or meta.get("sample_id") or sample_id
        mrn = meta.get("mrn") or ""
        if code or mrn:
            try:
                clinical_presentation_store.save(
                    code=code,
                    mrn=mrn,
                    content=str(payload.get("clinical_description") or ""),
                )
            except ValueError:
                pass
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    meta["updated_at"] = now
    meta.setdefault("created_at", now)
    _write_json(p, meta)
    try:
        from . import sample_loader
        sample_loader.update_case_table_row(sample_id)
    except Exception as e:
        print(f"[case-table] report save refresh failed for {sample_id}: {e}", flush=True)
    return _project_reviewer(meta)
