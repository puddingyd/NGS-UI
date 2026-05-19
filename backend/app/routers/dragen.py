"""HTTP endpoints for the 三級分析 flow.

  GET  /api/dragen/vcfs            cached index of DRAGEN + in-house VCFs
  POST /api/dragen/index/refresh   rescan both root sets, persist
  POST /api/dragen/jobs            spawn a worker; returns job_id
  GET  /api/dragen/jobs            list known jobs (for topbar status)
  GET  /api/dragen/jobs/{jid}      state + log tail for one job
"""
from __future__ import annotations

import asyncio

from fastapi import APIRouter, BackgroundTasks, Body, Depends, HTTPException

from ..auth import current_user
from ..services import dragen_jobs

router = APIRouter(
    prefix="/api/dragen",
    tags=["dragen"],
    dependencies=[Depends(current_user)],
)


def _meta(idx: dict | None) -> dict:
    if not idx:
        return {
            "updated_at": None,
            "dragen_count": 0,
            "inhouse_count": 0,
            "scan_duration_sec": None,
            "stale": True,
        }
    return {
        "updated_at":        idx.get("updated_at"),
        "dragen_count":      len(idx.get("dragen", [])),
        "inhouse_count":     len(idx.get("inhouse", [])),
        "scan_duration_sec": idx.get("scan_duration_sec"),
        "stale":             dragen_jobs.index_is_stale(idx),
    }


@router.get("/vcfs")
async def get_pipeline_vcfs(background: BackgroundTasks):
    """Return the cached index of DRAGEN + in-house VCFs.

    First call (no file yet) builds the index synchronously so the
    modal has something to show. Subsequent stale calls return the old
    index immediately and schedule a background refresh — the UI will
    pick up the new data on the next open or via a manual refresh.
    """
    idx = dragen_jobs.load_index()
    if idx is None:
        # Cold start: scan synchronously so the user sees data on the
        # very first 三級分析 click. Bounded by find(1) speed.
        idx = await asyncio.to_thread(dragen_jobs.refresh_index)
    elif dragen_jobs.index_is_stale(idx):
        background.add_task(dragen_jobs.refresh_index)
    return {
        "meta":    _meta(idx),
        "dragen":  idx.get("dragen", []),
        "inhouse": idx.get("inhouse", []),
    }


@router.post("/index/refresh")
async def post_refresh_index():
    """Force a rescan now (manual 🔄 button). Synchronous; returns new meta."""
    idx = await asyncio.to_thread(dragen_jobs.refresh_index)
    return {
        "meta":    _meta(idx),
        "dragen":  idx.get("dragen", []),
        "inhouse": idx.get("inhouse", []),
    }


@router.post("/jobs")
def post_dragen_job(payload: dict = Body(...)):
    """Spawn a pipeline worker.

    Body:
      mode:           "dragen" (default) | "inhouse"
      vcf_path:       path to the SNV/Indel VCF (the anchor)
      sample_id:      e.g. VAL-58-dragen / VAL-31-inhouse
      with_extra_vep: bool (default true)
      cnv_vcf:        in-house: gcnv VCF       (ignored for dragen)
      sv_vcf:         in-house: delly VCF      (ignored for dragen)
      mito_vcf:       in-house: mito VCF       (ignored for dragen)
    """
    mode = (payload.get("mode") or "dragen").strip()
    vcf  = (payload.get("vcf_path") or "").strip()
    sid  = (payload.get("sample_id") or "").strip()
    with_extra_vep = bool(payload.get("with_extra_vep", True))
    cnv_vcf  = (payload.get("cnv_vcf")  or "").strip()
    sv_vcf   = (payload.get("sv_vcf")   or "").strip()
    mito_vcf = (payload.get("mito_vcf") or "").strip()
    if not vcf:
        raise HTTPException(400, "vcf_path required")
    if not sid:
        raise HTTPException(400, "sample_id required")
    try:
        job_id = dragen_jobs.start_job(
            vcf, sid,
            mode=mode,
            with_extra_vep=with_extra_vep,
            cnv_vcf=cnv_vcf,
            sv_vcf=sv_vcf,
            mito_vcf=mito_vcf,
        )
    except FileNotFoundError as e:
        raise HTTPException(404, str(e))
    except ValueError as e:
        raise HTTPException(400, str(e))
    return {"job_id": job_id}


@router.get("/jobs")
def get_dragen_jobs():
    """Recent jobs for the topbar status / modal listing."""
    return dragen_jobs.list_jobs()


@router.get("/jobs/{job_id}")
def get_dragen_job(job_id: str):
    state = dragen_jobs.load_state(job_id)
    if state is None:
        raise HTTPException(404, f"job not found: {job_id}")
    state = dict(state)
    state["running"]  = dragen_jobs.is_running(job_id)
    state["log_tail"] = dragen_jobs.tail_log(job_id, n=80)
    return state
