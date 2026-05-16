"""HTTP endpoints for the 三級分析 (DRAGEN VCF → GUI sample) flow.

  GET  /api/dragen/vcfs          list available DRAGEN VCFs to pick from
  POST /api/dragen/jobs          spawn a worker; returns job_id
  GET  /api/dragen/jobs          list known jobs (for topbar status)
  GET  /api/dragen/jobs/{jid}    state + log tail for one job
"""
from __future__ import annotations

from fastapi import APIRouter, Body, Depends, HTTPException

from ..auth import current_user
from ..services import dragen_jobs

router = APIRouter(
    prefix="/api/dragen",
    tags=["dragen"],
    dependencies=[Depends(current_user)],
)


@router.get("/vcfs")
def get_dragen_vcfs():
    """List hard-filtered VCFs under the configured roots."""
    return dragen_jobs.list_dragen_vcfs()


@router.post("/jobs")
def post_dragen_job(payload: dict = Body(...)):
    """Spawn a DRAGEN pipeline worker.

    Body: {"vcf_path": "...", "sample_id": "...", "with_extra_vep": true}
    """
    vcf = (payload.get("vcf_path") or "").strip()
    sid = (payload.get("sample_id") or "").strip()
    with_extra_vep = bool(payload.get("with_extra_vep", True))
    if not vcf:
        raise HTTPException(400, "vcf_path required")
    if not sid:
        raise HTTPException(400, "sample_id required")
    try:
        job_id = dragen_jobs.start_job(vcf, sid, with_extra_vep=with_extra_vep)
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
