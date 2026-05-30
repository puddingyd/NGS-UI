"""IGV.js helper endpoints — list candidate BAMs for a sample and
stream BAM/BAI files with HTTP range support so the browser-side
igv.js can read alignments directly.

BAMs live outside the NGS-UI tree, on the pipeline's nextflow output
share. We constrain reads to the configured root(s) and only allow
BAM/BAI/CRAM/CRAI extensions so this can't double as a generic file
proxy. Login required (same as the rest of /api).
"""
from __future__ import annotations

import os
import re
from pathlib import Path

from fastapi import APIRouter, Depends, HTTPException, Query, Request
from fastapi.responses import StreamingResponse

from ..auth import current_user

router = APIRouter(
    prefix="/api/igv",
    tags=["igv"],
    dependencies=[Depends(current_user)],
)

# Where pipeline alignments live. Override via env on the dev box.
_BAM_ROOTS = [
    Path(p).resolve() for p in (
        os.environ.get("NGS_UI_BAM_ROOT")
        or "/home/datalake_Intermediate/pipeline/nextflow_output"
    ).split(":") if p
]

_SID_RE = re.compile(r"^[A-Za-z0-9_-]{1,32}$")
_ALLOWED_SUFFIXES = (".bam", ".bai", ".cram", ".crai")
_CHUNK = 1024 * 1024
_SIBLING_LIMIT = 2


def _validate_sid(sid: str) -> str:
    sid = (sid or "").strip()
    if not _SID_RE.match(sid):
        raise HTTPException(400, "invalid sample id")
    return sid


def _bam_for(sid: str, batch: Path) -> Path | None:
    p = batch / sid / "02_alignment" / f"{sid}.aligned.sorted.bam"
    return p if p.is_file() else None


def _path_ok(p: Path) -> bool:
    try:
        rp = p.resolve()
    except OSError:
        return False
    for root in _BAM_ROOTS:
        try:
            rp.relative_to(root)
            return True
        except ValueError:
            continue
    return False


def _bam_entry(sid: str, batch: Path, path: Path) -> dict:
    return {
        "label":     f"{sid} ({batch.name})",
        "sample_id": sid,
        "batch":     batch.name,
        "path":      str(path),
    }


@router.get("/bams")
def list_bams(sample_id: str = Query(...)):
    """Find the BAM for `sample_id` plus up to 2 sibling BAMs from the
    same batch (for trio / comparison viewing). The frontend can add /
    remove / swap these in the IGV modal."""
    sid = _validate_sid(sample_id)
    primary = None
    siblings: list[dict] = []
    for root in _BAM_ROOTS:
        if not root.is_dir():
            continue
        for batch in sorted(root.iterdir()):
            if not batch.is_dir():
                continue
            p = _bam_for(sid, batch)
            if not p:
                continue
            primary = _bam_entry(sid, batch, p)
            for sib in sorted(batch.iterdir()):
                if not sib.is_dir() or sib.name == sid:
                    continue
                sp = _bam_for(sib.name, batch)
                if sp:
                    siblings.append(_bam_entry(sib.name, batch, sp))
                if len(siblings) >= _SIBLING_LIMIT:
                    break
            break
        if primary:
            break
    return {"primary": primary, "siblings": siblings, "roots": [str(r) for r in _BAM_ROOTS]}


@router.get("/batch-samples")
def list_batch_samples(batch: str = Query(...)):
    """List every sample with a BAM in the given batch — used by the
    "add another BAM" picker so reviewers can pull in any sibling."""
    batch = (batch or "").strip()
    if not re.match(r"^[A-Za-z0-9_.\-]{1,64}$", batch):
        raise HTTPException(400, "invalid batch name")
    out: list[dict] = []
    for root in _BAM_ROOTS:
        bd = root / batch
        if not bd.is_dir():
            continue
        for d in sorted(bd.iterdir()):
            if not d.is_dir():
                continue
            p = _bam_for(d.name, bd)
            if p:
                out.append(_bam_entry(d.name, bd, p))
    return {"batch": batch, "samples": out}


@router.get("/file")
def serve_file(request: Request, path: str = Query(...)):
    """Stream a BAM/BAI file with HTTP range support for igv.js."""
    p = Path(path)
    if not p.is_file():
        raise HTTPException(404, "not found")
    if not _path_ok(p):
        raise HTTPException(403, "forbidden path")
    if p.suffix.lower() not in _ALLOWED_SUFFIXES:
        raise HTTPException(400, "unsupported file type")

    size = p.stat().st_size
    rng = request.headers.get("range")
    if not rng:
        def whole():
            with p.open("rb") as fh:
                while True:
                    b = fh.read(_CHUNK)
                    if not b: break
                    yield b
        return StreamingResponse(
            whole(),
            media_type="application/octet-stream",
            headers={"content-length": str(size), "accept-ranges": "bytes"},
        )

    m = re.match(r"bytes=(\d+)-(\d*)", rng.strip())
    if not m:
        raise HTTPException(416, "invalid range")
    start = int(m.group(1))
    end   = int(m.group(2)) if m.group(2) else size - 1
    end   = min(end, size - 1)
    if start > end:
        raise HTTPException(416, "range out of bounds")
    length = end - start + 1

    def part():
        with p.open("rb") as fh:
            fh.seek(start)
            remaining = length
            while remaining > 0:
                b = fh.read(min(_CHUNK, remaining))
                if not b: break
                remaining -= len(b)
                yield b

    return StreamingResponse(
        part(),
        status_code=206,
        media_type="application/octet-stream",
        headers={
            "content-range":  f"bytes {start}-{end}/{size}",
            "content-length": str(length),
            "accept-ranges":  "bytes",
        },
    )
