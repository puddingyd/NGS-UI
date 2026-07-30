"""Reviewer-side state (status / edits / panels / comment / tags …).

These fields live on the per-patient sample_metadata.json so they
survive across analysis versions. Pipeline-owned keys (lis_id, name,
mrn, test_type, vcf_path, …) are preserved untouched on save() — only
the whitelist below gets overwritten.

Structured manual ACMG revisions and active cross-case observations live in
``manual_acmg.sqlite``; the sample-specific final snapshot remains here so all
report consumers use the exact result shown on that sample.
"""
from __future__ import annotations

import json
import threading
from datetime import datetime, timezone
from pathlib import Path

from . import clinical_presentation_store, sample_layout


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
_WRITE_LOCK = threading.RLock()


def _meta_path(sample_id: str) -> Path:
    return sample_layout.state_file(sample_id, "sample_metadata.json")


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


def save(sample_id: str, payload: dict, *, user: dict | None = None) -> dict:
    """Merge reviewer fields into sample_metadata.json."""
    with _WRITE_LOCK:
        p = _meta_path(sample_id)
        meta = _read_json(p)
        payload = payload or {}
        now = datetime.now(timezone.utc).isoformat(timespec="seconds")
        if "status" in payload and isinstance(payload.get("status"), dict):
            old_status = meta.get("status") if isinstance(meta.get("status"), dict) else {}
            audit = meta.get("status_audit") if isinstance(meta.get("status_audit"), dict) else {}
            audit = dict(audit)
            for variant_id in set(old_status) | set(payload["status"]):
                old = str(old_status.get(variant_id) or "")
                new = str(payload["status"].get(variant_id) or "")
                if old == new:
                    continue
                if new in {"1", "2"}:
                    audit[variant_id] = {
                        "reviewer_user_id": (user or {}).get("id"),
                        "reviewer_username": (user or {}).get("username") or "",
                        "updated_at": now,
                    }
                else:
                    audit.pop(variant_id, None)
            meta["status_audit"] = audit
        if isinstance(payload.get("edits"), dict):
            incoming_edits = {
                variant_id: dict(edit) if isinstance(edit, dict) else edit
                for variant_id, edit in payload["edits"].items()
            }
            stored_edits = meta.get("edits") if isinstance(meta.get("edits"), dict) else {}
            for variant_id, stored_edit in stored_edits.items():
                if not isinstance(stored_edit, dict):
                    continue
                stored_snapshot = stored_edit.get("manual_acmg")
                if not isinstance(stored_snapshot, dict):
                    continue
                incoming_edit = incoming_edits.get(variant_id)
                if not isinstance(incoming_edit, dict):
                    incoming_edit = {}
                    incoming_edits[variant_id] = incoming_edit
                incoming_snapshot = incoming_edit.get("manual_acmg")
                try:
                    stored_revision = int(stored_snapshot.get("revision_id") or 0)
                    incoming_revision = int(
                        (incoming_snapshot or {}).get("revision_id") or 0
                    )
                except (TypeError, ValueError):
                    stored_revision, incoming_revision = 1, 0
                if not isinstance(incoming_snapshot, dict) or stored_revision > incoming_revision:
                    incoming_edit["manual_acmg"] = stored_snapshot
                    incoming_edit["ACMG_classification"] = stored_snapshot.get(
                        "classification", ""
                    )
                    incoming_edit["ACMG_score"] = stored_snapshot.get("score")
                    incoming_edit["ACMG_criteria"] = stored_snapshot.get(
                        "criteria_text", ""
                    )
            payload = dict(payload)
            payload["edits"] = incoming_edits
        for k in _REVIEWER_FIELDS:
            if k in payload:
                meta[k] = payload[k]
        if user:
            meta["last_reviewer_user_id"] = user.get("id")
            meta["last_reviewer_username"] = user.get("username") or ""
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
        meta["updated_at"] = now
        meta.setdefault("created_at", now)
        _write_json(p, meta)
    try:
        from . import manual_acmg
        manual_acmg.sync_observations(
            meta.get("genome_build") or "hg38",
            sample_id,
            meta.get("status") or {},
            reviewer_user_id=(user or {}).get("id"),
            reviewer_username=(user or {}).get("username") or "",
            updated_at=now,
            status_audit=meta.get("status_audit") or {},
        )
    except Exception as e:
        print(f"[manual-acmg] observation sync failed for {sample_id}: {e}", flush=True)
    try:
        from . import sample_loader
        sample_loader.update_case_table_row(sample_id)
    except Exception as e:
        print(f"[case-table] report save refresh failed for {sample_id}: {e}", flush=True)
    return _project_reviewer(meta)


def save_manual_acmg(
    sample_id: str,
    variant_id: str,
    assertion: dict,
    *,
    user: dict | None = None,
) -> dict:
    """Atomically merge one structured ACMG snapshot into sample edits."""
    with _WRITE_LOCK:
        p = _meta_path(sample_id)
        meta = _read_json(p)
        edits = meta.get("edits") if isinstance(meta.get("edits"), dict) else {}
        edits = dict(edits)
        variant_edits = edits.get(variant_id) if isinstance(edits.get(variant_id), dict) else {}
        variant_edits = dict(variant_edits)
        snapshot = {
            "revision_id": assertion.get("revision_id"),
            "genome_build": assertion.get("genome_build"),
            "variant_id": assertion.get("variant_id"),
            "criteria": assertion.get("criteria") or {},
            "criteria_text": assertion.get("criteria_text") or "",
            "score": assertion.get("score"),
            "classification": assertion.get("classification") or "",
            "reviewer_user_id": assertion.get("reviewer_user_id"),
            "reviewer_username": assertion.get("reviewer_username") or "",
            "source_sample_id": assertion.get("source_sample_id") or sample_id,
            "created_at": assertion.get("created_at"),
        }
        variant_edits["manual_acmg"] = snapshot
        # Keep legacy readers and existing exported views compatible.
        variant_edits["ACMG_classification"] = snapshot["classification"]
        variant_edits["ACMG_score"] = snapshot["score"]
        variant_edits["ACMG_criteria"] = snapshot["criteria_text"]
        edits[variant_id] = variant_edits
        meta["edits"] = edits
        now = datetime.now(timezone.utc).isoformat(timespec="seconds")
        meta["updated_at"] = now
        meta.setdefault("created_at", now)
        if user:
            meta["last_reviewer_user_id"] = user.get("id")
            meta["last_reviewer_username"] = user.get("username") or ""
        _write_json(p, meta)
    try:
        from . import sample_loader
        sample_loader.update_case_table_row(sample_id)
    except Exception as e:
        print(f"[case-table] manual ACMG refresh failed for {sample_id}: {e}", flush=True)
    return snapshot
