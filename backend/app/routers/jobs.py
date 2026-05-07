"""Background job endpoints (Phase C).

Endpoints:
  POST /api/samples/{sample_id}/jobs/exomiser_lirical → enqueue
  GET  /api/samples/{sample_id}/jobs                  → recent jobs
  GET  /api/jobs/{job_id}                             → single job

Live status comes from RQ when reachable; otherwise we fall back to
the JSON sidecar in data/jobs/.
"""
from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException

from ..auth import current_user
from ..config import REDIS_URL, TERTIARY_OUTPUT_ROOT
from ..services import analyses_store, job_store

router = APIRouter(prefix="/api", tags=["jobs"], dependencies=[Depends(current_user)])


def _rq_status(job_id: str) -> dict:
    """Look up the live RQ Job; returns {status, exc_info?} or {} if unreachable."""
    try:
        from redis import Redis
        from rq.job import Job
        conn = Redis.from_url(REDIS_URL)
        job = Job.fetch(job_id, connection=conn)
        out = {"rq_status": job.get_status(refresh=True)}
        if job.exc_info:
            out["rq_error"] = str(job.exc_info)[-1000:]
        return out
    except Exception:
        return {}


def _enqueue(
    sample_id: str,
    kind: str = "exomiser_lirical",
    version: str | None = None,
) -> dict:
    sub = TERTIARY_OUTPUT_ROOT / sample_id
    if not sub.is_dir():
        raise HTTPException(404, f"sample not found: {sample_id}")
    try:
        from redis import Redis
        from rq import Queue
    except ImportError as e:
        raise HTTPException(500, f"RQ not installed: {e}")

    # Worker pulls HPO + sidecar paths from active_version. If the
    # caller asked for a specific version, switch active before
    # enqueueing so the worker lands on the right one. Validates the
    # name exists; falls back to the existing active otherwise.
    if version:
        try:
            analyses_store.set_active(sample_id, version)
        except ValueError as e:
            raise HTTPException(400, str(e))

    job_id = job_store.new_job_id()
    record = job_store.write({
        "job_id":    job_id,
        "sample_id": sample_id,
        "kind":      kind,
        "version":   version or analyses_store.active_version(sample_id),
        "status":    "queued",
        "step":      "queued",
    })
    try:
        conn = Redis.from_url(REDIS_URL)
        q = Queue("ngs-ui", connection=conn, default_timeout=4 * 60 * 60)
        q.enqueue(
            "app.workers.exomiser_lirical.run_exomiser_lirical",
            args=(job_id, sample_id),
            job_id=job_id,
            result_ttl=24 * 3600,
            failure_ttl=24 * 3600,
        )
    except Exception as e:
        job_store.update(job_id, {"status": "failed", "step": "enqueue", "error": str(e)})
        raise HTTPException(500, f"could not enqueue job (is redis up?): {e}")
    return record


@router.post("/samples/{sample_id}/jobs/exomiser_lirical")
def post_exomiser_lirical(sample_id: str, payload: dict | None = None):
    """Body: { version?: "default" | "v2_seizure" | ... }

    Omitted version → use the sample's current active analysis.
    """
    version = (payload or {}).get("version")
    return _enqueue(sample_id, "exomiser_lirical", version=version)


@router.get("/samples/{sample_id}/jobs")
def list_sample_jobs(sample_id: str):
    jobs = job_store.list_for_sample(sample_id, limit=20)
    for j in jobs:
        j.update(_rq_status(j["job_id"]))
    return jobs


@router.get("/jobs/{job_id}")
def get_job(job_id: str):
    j = job_store.read(job_id)
    if not j:
        raise HTTPException(404, "unknown job")
    j.update(_rq_status(job_id))
    return j
