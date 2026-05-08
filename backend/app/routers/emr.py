"""EMR fetch endpoints.

  GET  /api/emr/{mrn}                   read-only probe; returns
                                        {phenotype, consultation,
                                         fetched_at} for one MRN.

  POST /api/samples/{id}/sync_emr       re-fetch + merge into the
                                        sample's metadata. Updates
                                        sex (overwrite), dob,
                                        genetic_counseling, plus
                                        emr_synced_at; the response
                                        echoes what changed.

When the EMR client_id env var is missing, both endpoints surface a
503 so the UI can hide the sync button instead of silently doing
nothing.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from fastapi import APIRouter, Depends, HTTPException

from ..auth import current_user
from ..config import TERTIARY_OUTPUT_ROOT
from ..services import emr_client

router = APIRouter(prefix="/api", tags=["emr"], dependencies=[Depends(current_user)])


def _require_enabled() -> None:
    if not emr_client.is_enabled():
        raise HTTPException(503, "EMR client_id not configured (set NGS_UI_EMR_CLIENT_ID)")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@router.get("/emr/enabled")
def emr_enabled():
    """Public probe so the UI can hide EMR controls when client_id
    isn't configured. No data leaks here — just a boolean."""
    return {"enabled": emr_client.is_enabled()}


@router.get("/emr/{mrn}")
def emr_lookup(mrn: str):
    """Read-only fetch. The phenotype API works without the client_id
    (it's behind the same intranet-only HIS host as everything else),
    so allow that half even if the APIM key is missing — the UI can
    still show the EMR-side phenotype, just not the consultation
    record."""
    if not mrn:
        raise HTTPException(400, "mrn required")
    return emr_client.fetch(mrn)


@router.post("/samples/{sample_id}/sync_emr")
def sync_emr(sample_id: str):
    """Re-fetch EMR for the sample's stored MRN and merge into
    sample_metadata.json. Sex from EMR overwrites whatever was there;
    dob + genetic_counseling get filled in even when the reviewer
    hadn't entered anything. Returns the changed fields + the raw
    EMR payload so the UI can show a diff."""
    _require_enabled()
    sub = TERTIARY_OUTPUT_ROOT / sample_id
    meta_path = sub / "sample_metadata.json"
    if not meta_path.is_file():
        raise HTTPException(404, f"sample not found: {sample_id}")
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if not isinstance(meta, dict):
            raise ValueError("malformed sample_metadata.json")
    except (json.JSONDecodeError, ValueError) as e:
        raise HTTPException(500, f"meta load failed: {e}")
    mrn = (meta.get("mrn") or "").strip()
    if not mrn:
        raise HTTPException(400, "sample_metadata.json has no MRN")

    payload = emr_client.fetch(mrn)
    consult = payload.get("consultation") or {}
    pheno = payload.get("phenotype") or {}

    changes: dict[str, dict] = {}
    if consult.get("sex"):
        old = meta.get("sex", "")
        if old != consult["sex"]:
            changes["sex"] = {"from": old, "to": consult["sex"]}
        meta["sex"] = consult["sex"]
    if consult.get("date_of_birth"):
        old = meta.get("date_of_birth", "")
        if old != consult["date_of_birth"]:
            changes["date_of_birth"] = {"from": old, "to": consult["date_of_birth"]}
        meta["date_of_birth"] = consult["date_of_birth"]
    if consult.get("text"):
        old = meta.get("genetic_counseling", "")
        if old != consult["text"]:
            changes["genetic_counseling"] = {"chars": len(consult["text"])}
        meta["genetic_counseling"] = consult["text"]
    meta["emr_synced_at"] = _now()
    meta["updated_at"] = meta["emr_synced_at"]
    meta_path.write_text(
        json.dumps(meta, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return {
        "sample_id": sample_id,
        "mrn":       mrn,
        "changes":   changes,
        "phenotype": pheno,    # for completeness; reviewer txt still wins for HPO
        "consultation_found": bool(consult.get("found")),
    }
