"""Per-sample analysis-version CRUD.

Endpoints:
  GET    /api/samples/{id}/analyses
  POST   /api/samples/{id}/analyses             (create or overwrite)
  DELETE /api/samples/{id}/analyses/{name}
  PUT    /api/samples/{id}/active_analysis
"""
from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException

from ..auth import current_user
from ..config import TERTIARY_OUTPUT_ROOT
from ..services import analyses_store

router = APIRouter(prefix="/api", tags=["analyses"], dependencies=[Depends(current_user)])


def _require_sample(sample_id: str) -> None:
    if not (TERTIARY_OUTPUT_ROOT / sample_id).is_dir():
        raise HTTPException(404, f"sample not found: {sample_id}")


@router.get("/samples/{sample_id}/analyses")
def list_versions(sample_id: str):
    _require_sample(sample_id)
    return {
        "active":   analyses_store.active_version(sample_id),
        "versions": analyses_store.list_versions(sample_id),
    }


@router.post("/samples/{sample_id}/analyses")
def create_or_update_version(sample_id: str, payload: dict):
    """Body: { name, hpo: [...], panels: [...], note?, set_active?, clear_sidecars? }

    `clear_sidecars=true` wipes pheno_score / exomiser / lirical sidecars
    + the analysis_files/ dir before writing analysis.json. Use it when
    overwriting an existing version about to be re-run, so a partial
    failure mid-run can't blend old + new outputs.
    """
    _require_sample(sample_id)
    name = (payload or {}).get("name") or ""
    try:
        analyses_store.validate_name(name)
    except ValueError as e:
        raise HTTPException(400, str(e))

    hpo    = payload.get("hpo") or []
    panels = payload.get("panels") or []
    note   = payload.get("note", "")

    if payload.get("clear_sidecars"):
        analyses_store.clear_sidecars(sample_id, name)

    written = analyses_store.write_version(
        sample_id, name, hpo=hpo, panels=panels, note=note,
    )

    if payload.get("set_active"):
        try:
            analyses_store.set_active(sample_id, name)
        except ValueError as e:
            raise HTTPException(400, str(e))

    return {
        "name":     name,
        **written,
        "is_active": analyses_store.active_version(sample_id) == name,
    }


@router.delete("/samples/{sample_id}/analyses/{name}")
def delete_version(sample_id: str, name: str):
    _require_sample(sample_id)
    try:
        analyses_store.validate_name(name)
    except ValueError as e:
        raise HTTPException(400, str(e))
    try:
        removed = analyses_store.delete_version(sample_id, name)
    except ValueError as e:
        # Reserved name (default).
        raise HTTPException(400, str(e))
    if not removed:
        raise HTTPException(404, f"version not found: {name}")

    # If this was the active version, fall back to default if it
    # exists, else the first remaining version, else None.
    current = analyses_store.active_version(sample_id)
    if current is None:
        return {"deleted": name, "active": None}
    return {"deleted": name, "active": current}


@router.put("/samples/{sample_id}/active_analysis")
def set_active_version(sample_id: str, payload: dict):
    _require_sample(sample_id)
    name = (payload or {}).get("name") or ""
    try:
        analyses_store.set_active(sample_id, name)
    except ValueError as e:
        raise HTTPException(400, str(e))
    return {"active": name}
