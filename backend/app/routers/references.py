"""Reference database versions (home page) and the HPO update button.

  GET  /api/references            versions of the databases in use (no login: no patient data)
  GET  /api/references/hpo/update  state of the HPO update job
  POST /api/references/hpo/update  start downloading / building / activating the latest HPO
"""
from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException

from ..auth import current_user
from ..services import hpo_update, reference_versions

router = APIRouter(prefix="/api/references", tags=["references"])


@router.get("")
def get_versions():
    return reference_versions.collect()


@router.get("/hpo/update")
def get_hpo_update(user: dict = Depends(current_user)):
    return hpo_update.status()


@router.post("/hpo/update")
def post_hpo_update(user: dict = Depends(current_user)):
    try:
        return hpo_update.start_background("latest", user=str(user.get("username") or ""))
    except OSError as exc:
        raise HTTPException(500, f"HPO 更新無法啟動：{exc}") from exc
