"""Authenticated, sample-scoped targeted somatic analysis endpoints."""
from fastapi import APIRouter, Depends, HTTPException

from ..auth import current_user
from ..services import somatic, sample_layout

router = APIRouter(prefix="/api/samples/{sample_id}/somatic", tags=["somatic"],
                   dependencies=[Depends(current_user)])


def check_sample(sample_id):
    try:
        somatic.validate_sid(sample_id)
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc
    if not sample_layout.state_file(sample_id, "sample_metadata.json").is_file():
        raise HTTPException(404, "找不到已載入個案")


def owned_job(sample_id, run_id):
    check_sample(sample_id)
    try:
        job = somatic.read_job(run_id)
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc
    if job.get("sample_id") != sample_id:
        raise HTTPException(404, "找不到分析工作")
    return job


@router.get("")
def status(sample_id: str):
    check_sample(sample_id)
    error = ""
    try:
        somatic.settings()
    except (ValueError, OSError) as exc:
        error = str(exc)
    return {"summary": somatic.summary(sample_id), "jobs": somatic.jobs(sample_id),
            "bams": somatic.sample_bams(sample_id), "configuration_error": error}


@router.post("/preview")
def preview(sample_id: str, payload: dict):
    check_sample(sample_id)
    try:
        return somatic.resolve_targets(payload, somatic.settings())
    except (ValueError, OSError) as exc:
        raise HTTPException(400, str(exc)) from exc


@router.post("/jobs")
def create(sample_id: str, payload: dict, user=Depends(current_user)):
    check_sample(sample_id)
    try:
        return somatic.start(sample_id, payload, user.get("username", ""))
    except (ValueError, OSError) as exc:
        raise HTTPException(400, str(exc)) from exc


@router.get("/jobs/{run_id}")
def detail(sample_id: str, run_id: str):
    job = owned_job(sample_id, run_id)
    log_path = somatic.job_dir(run_id) / "log.txt"
    log = ""
    if log_path.is_file():
        with log_path.open("rb") as handle:
            handle.seek(max(0, log_path.stat().st_size - 64000))
            log = handle.read().decode("utf-8", errors="replace")
    candidates = []
    if job.get("status") == "completed":
        candidates = somatic.read_json(somatic.result_dir(sample_id, run_id) / "candidates.json", [])
    return {"job": job, "log": log, "candidates": candidates}


@router.post("/jobs/{run_id}/cancel")
def cancel(sample_id: str, run_id: str):
    with somatic.submission_lock():
        job = owned_job(sample_id, run_id)
        if job.get("status") not in somatic.ACTIVE:
            raise HTTPException(409, "此工作已結束")
        (somatic.job_dir(run_id) / "cancel").touch()
    return {"cancel_requested": True}


@router.delete("/jobs/{run_id}")
def delete(sample_id: str, run_id: str):
    check_sample(sample_id)
    try:
        with somatic.submission_lock():
            return somatic.delete_run(sample_id, run_id)
    except FileNotFoundError as exc:
        raise HTTPException(404, str(exc)) from exc
    except RuntimeError as exc:
        raise HTTPException(409, str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc


@router.post("/jobs/{run_id}/include")
def include(sample_id: str, run_id: str, payload: dict):
    owned_job(sample_id, run_id)
    try:
        with somatic.submission_lock():
            somatic.select_filtered(sample_id, run_id, str(payload.get("variant_id", "")))
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc
    return {"ok": True}
