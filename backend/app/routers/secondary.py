"""HTTP endpoints for secondary-analysis FASTQ samplesheet preparation."""
from __future__ import annotations

import asyncio

from fastapi import APIRouter, BackgroundTasks, Body, Depends, HTTPException

from ..auth import current_user
from ..services import secondary_analysis

router = APIRouter(
    prefix="/api/secondary",
    tags=["secondary"],
    dependencies=[Depends(current_user)],
)


def _meta(idx: dict | None) -> dict:
    if not idx:
        return {
            "updated_at": None,
            "wes_count": 0,
            "wgs_count": 0,
            "wgs_lane_count": 0,
            "scan_duration_sec": None,
            "stale": True,
        }
    return {
        "updated_at": idx.get("updated_at"),
        "wes_count": len(idx.get("wes", [])),
        "wgs_count": len(idx.get("wgs", [])),
        "wgs_lane_count": sum(int(row.get("lane_count") or 1) for row in idx.get("wgs", [])),
        "scan_duration_sec": idx.get("scan_duration_sec"),
        "stale": secondary_analysis.index_is_stale(idx),
    }


@router.get("/fastqs")
async def get_fastqs(background: BackgroundTasks):
    idx = secondary_analysis.load_index()
    if idx is None:
        idx = await asyncio.to_thread(secondary_analysis.refresh_index)
    elif secondary_analysis.index_is_stale(idx):
        background.add_task(secondary_analysis.refresh_index)
    return {
        "meta": _meta(idx),
        "wes": idx.get("wes", []),
        "wgs": idx.get("wgs", []),
    }


@router.post("/index/refresh")
async def post_refresh_index():
    idx = await asyncio.to_thread(secondary_analysis.refresh_index)
    return {
        "meta": _meta(idx),
        "wes": idx.get("wes", []),
        "wgs": idx.get("wgs", []),
    }


@router.get("/nf-work/cleanup-command")
def get_cleanup_nf_work_command():
    try:
        return secondary_analysis.cleanup_nf_work_command()
    except ValueError as e:
        raise HTTPException(400, str(e))


@router.post("/samplesheet")
def post_samplesheet(payload: dict = Body(...)):
    seq_type = (payload.get("seq_type") or "").strip()
    batch_name = (payload.get("batch_name") or "").strip()
    samples = payload.get("samples")
    if not isinstance(samples, list) or not samples:
        raise HTTPException(400, "samples must be a non-empty list")
    try:
        return secondary_analysis.create_samplesheet(seq_type, samples, batch_name=batch_name)
    except FileNotFoundError as e:
        raise HTTPException(404, str(e))
    except ValueError as e:
        raise HTTPException(400, str(e))
    except OSError as e:
        raise HTTPException(500, f"建立 samplesheet 失敗：{e}")
