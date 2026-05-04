"""Per-sample job metadata stored as `data/jobs/{job_id}.json`.

RQ holds the live job in Redis (status / result / start time), but
those entries get pruned. We keep our own JSON sidecar so the UI can
show a sample's run history even after Redis evicts the entry.
"""
from __future__ import annotations

import json
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

from ..config import JOBS_DIR


def _path(job_id: str) -> Path:
    return JOBS_DIR / f"{job_id}.json"


def new_job_id() -> str:
    return f"job_{int(time.time())}_{uuid.uuid4().hex[:8]}"


def write(job: dict) -> dict:
    if "job_id" not in job:
        raise ValueError("job dict needs a job_id")
    job["updated_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    _path(job["job_id"]).write_text(
        json.dumps(job, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return job


def update(job_id: str, patch: dict) -> dict:
    j = read(job_id) or {"job_id": job_id}
    j.update(patch)
    return write(j)


def read(job_id: str) -> dict | None:
    p = _path(job_id)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def list_for_sample(sample_id: str, limit: int = 20) -> list[dict]:
    out: list[dict] = []
    for p in sorted(JOBS_DIR.glob("*.json"), reverse=True):
        try:
            j = json.loads(p.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        if j.get("sample_id") == sample_id:
            out.append(j)
            if len(out) >= limit:
                break
    return out
