"""IGV.js helper endpoints — list candidate BAMs for a sample and
stream BAM/BAI files with HTTP range support so the browser-side
igv.js can read alignments directly.

BAMs live outside the NGS-UI tree, on the pipeline's nextflow output
share. We constrain reads to the configured root(s) and only allow
BAM/BAI/CRAM/CRAI extensions so this can't double as a generic file
proxy. Login required (same as the rest of /api).
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path
from urllib.parse import quote

from fastapi import APIRouter, Depends, HTTPException, Query, Request
from fastapi.responses import StreamingResponse

from ..auth import current_user
from ..config import TERTIARY_OUTPUT_ROOT

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

# Reference genome dir (fasta + .fai). Shipped to igv.js as a custom
# genome config so we can stay on hospital intranet — the default
# igv.js hg38 config points at AWS S3, which is firewalled here.
_REF_DIR = Path(os.environ.get(
    "NGS_UI_IGV_REF_DIR",
    "/home/pipeline/reference/hg38")).resolve()

_SID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
_ALLOWED_SUFFIXES = (".bam", ".bai", ".cram", ".crai",
                     ".fa", ".fai", ".fasta", ".dict", ".gzi", ".2bit")
_CHUNK = 1024 * 1024
_SIBLING_LIMIT = 2
_LEGACY_ALIAS_SUFFIXES = ("-dragen", "-inhouse", "-WES", "-WGS")


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
    for root in (*_BAM_ROOTS, _REF_DIR):
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


def _bam_hits(sid: str) -> list[dict]:
    out: list[dict] = []
    for root in _BAM_ROOTS:
        if not root.is_dir():
            continue
        for batch in sorted(root.iterdir()):
            if not batch.is_dir():
                continue
            p = _bam_for(sid, batch)
            if p:
                out.append(_bam_entry(sid, batch, p))
    return out


def _source_sid_from_sidecar(sid: str) -> str:
    path = TERTIARY_OUTPUT_ROOT / sid / "pipeline_source.json"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return ""
    source_sid = str(data.get("source_sample_id") or "").strip()
    return source_sid if _SID_RE.match(source_sid) else ""


def _legacy_source_sid(sid: str) -> str:
    """Resolve old alias cases only when exactly one BAM matches."""
    candidates = [
        sid.removesuffix(suffix)
        for suffix in _LEGACY_ALIAS_SUFFIXES
        if sid.endswith(suffix)
    ]
    hits = [
        hit
        for candidate in candidates
        for hit in _bam_hits(candidate)
    ]
    return hits[0]["sample_id"] if len(hits) == 1 else ""


def _resolve_primary_bam(sid: str) -> tuple[dict | None, str, str]:
    direct = _bam_hits(sid)
    if direct:
        return direct[0], sid, "exact"
    source_sid = _source_sid_from_sidecar(sid)
    if source_sid:
        hits = _bam_hits(source_sid)
        if hits:
            return hits[0], source_sid, "pipeline_source"
    legacy_sid = _legacy_source_sid(sid)
    if legacy_sid:
        hits = _bam_hits(legacy_sid)
        if hits:
            return hits[0], legacy_sid, "legacy_suffix"
    return None, sid, "not_found"


@router.get("/bams")
def list_bams(sample_id: str = Query(...)):
    """Find the BAM for `sample_id` plus up to 2 sibling BAMs from the
    same batch (for trio / comparison viewing). The frontend can add /
    remove / swap these in the IGV modal."""
    sid = _validate_sid(sample_id)
    primary, resolved_sid, resolution = _resolve_primary_bam(sid)
    siblings: list[dict] = []
    if primary:
        batch = Path(primary["path"]).parents[2]
        for sib in sorted(batch.iterdir()):
            if not sib.is_dir() or sib.name == resolved_sid:
                continue
            sp = _bam_for(sib.name, batch)
            if sp:
                siblings.append(_bam_entry(sib.name, batch, sp))
            if len(siblings) >= _SIBLING_LIMIT:
                break
    return {
        "primary": primary,
        "siblings": siblings,
        "requested_sample_id": sid,
        "resolved_sample_id": resolved_sid,
        "resolution": resolution,
        "roots": [str(r) for r in _BAM_ROOTS],
    }


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


# Common file names we've seen for the hg38 reference. First match wins.
_REF_FASTA_NAMES = (
    "hg38.fa", "hg38.fasta",
    "Homo_sapiens_assembly38.fasta",
    "GRCh38.fa", "GRCh38.fasta",
    "GRCh38.primary_assembly.genome.fa",
)


@router.get("/genome")
def igv_genome(build: str = Query("hg38")):
    """Return an igv.js custom-genome config pointing at our proxied
    reference fasta — hg38 only for now. The default igv.js hg38
    config points at AWS S3, which is blocked on the hospital network.

    Returns {ok: bool, config: {...}|null, fasta_path, fai_path, ...}.
    """
    if build != "hg38":
        raise HTTPException(400, "only hg38 supported")
    if not _REF_DIR.is_dir():
        return {"ok": False, "reason": f"ref dir not found: {_REF_DIR}"}
    fasta = None
    for name in _REF_FASTA_NAMES:
        p = _REF_DIR / name
        if p.is_file():
            fasta = p
            break
    if fasta is None:
        # last-ditch: any .fa / .fasta in the dir
        for p in sorted(_REF_DIR.iterdir()):
            if p.suffix.lower() in (".fa", ".fasta") and p.is_file():
                fasta = p
                break
    if fasta is None:
        return {"ok": False, "reason": f"no fasta under {_REF_DIR}"}
    fai = fasta.with_suffix(fasta.suffix + ".fai")
    if not fai.is_file():
        # Some refs are named hg38.fa + hg38.fai (no double-suffix).
        alt = fasta.with_suffix(".fai")
        if alt.is_file(): fai = alt
    if not fai.is_file():
        return {"ok": False, "reason": f"no .fai next to {fasta.name}"}
    return {
        "ok": True,
        "config": {
            "id":        "hg38",
            "name":      "Human (GRCh38/hg38)",
            "fastaURL":  f"/api/igv/file?path={quote(str(fasta), safe='')}",
            "indexURL":  f"/api/igv/file?path={quote(str(fai),   safe='')}",
        },
        "fasta_path": str(fasta),
        "fai_path":   str(fai),
        "ref_dir":    str(_REF_DIR),
    }


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
