"""IGV.js helper endpoints — list candidate BAMs for a sample and
stream BAM/BAI files with HTTP range support so the browser-side
igv.js can read alignments directly.

BAMs live outside the NGS-UI tree. DRAGEN BAMs use the raw Novaseq
layout `<run>/bam/<sample>.bam`; in-house BAMs use the Nextflow layout
`<batch>/<sample>/02_alignment/<sample>.aligned.sorted.bam`. We
constrain reads to the configured root(s) and only allow BAM/BAI/CRAM
/CRAI extensions so this can't double as a generic file proxy. Login
required (same as the rest of /api).
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
_INHOUSE_BAM_ROOTS = [
    Path(p).resolve() for p in (
        os.environ.get("NGS_UI_INHOUSE_BAM_ROOT")
        or os.environ.get("NGS_UI_BAM_ROOT")
        or "/home/datalake_Intermediate/pipeline/nextflow_output"
    ).split(":") if p
]
_DRAGEN_BAM_ROOTS = [
    Path(p).resolve() for p in (
        os.environ.get("NGS_UI_DRAGEN_BAM_ROOT")
        or "/home/datalake_Raw/Novaseq"
    ).split(":") if p
]
_ALL_BAM_ROOTS = tuple(dict.fromkeys([*_INHOUSE_BAM_ROOTS, *_DRAGEN_BAM_ROOTS]))

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


def _inhouse_bam_for(sid: str, batch: Path) -> Path | None:
    p = batch / sid / "02_alignment" / f"{sid}.aligned.sorted.bam"
    if p.is_file():
        return p
    align_dir = batch / sid / "02_alignment"
    if not align_dir.is_dir():
        return None
    hits = sorted(align_dir.glob("*.aligned.sorted.bam"))
    return hits[0] if len(hits) == 1 else None


def _dragen_bam_for(sid: str, run_dir: Path) -> Path | None:
    bam_dir = run_dir / "bam"
    if not bam_dir.is_dir():
        return None
    exact = bam_dir / f"{sid}.bam"
    if exact.is_file():
        return exact
    hits = [
        p for p in sorted(bam_dir.glob(f"{sid}*.bam"))
        if not p.name.endswith(".repeats.bam")
    ]
    return hits[0] if len(hits) == 1 else None


def _sample_id_from_dragen_bam(path: Path) -> str:
    name = path.name
    return name[:-len(".bam")] if name.endswith(".bam") else path.stem


def _sample_id_from_bam(path: Path) -> str:
    name = path.name
    if name.endswith(".aligned.sorted.bam"):
        return name[:-len(".aligned.sorted.bam")]
    if name.endswith(".bam"):
        return name[:-len(".bam")]
    return path.stem


def _bam_index_for(path: Path) -> Path | None:
    candidates = [
        Path(str(path) + ".bai"),
        path.with_suffix(".bai"),
    ]
    for p in candidates:
        if p.is_file():
            return p
    return None


def _path_ok(p: Path) -> bool:
    try:
        rp = p.resolve()
    except OSError:
        return False
    for root in (*_ALL_BAM_ROOTS, _REF_DIR):
        try:
            rp.relative_to(root)
            return True
        except ValueError:
            continue
    return False


def _dir_ok(p: Path) -> bool:
    try:
        rp = p.resolve()
    except OSError:
        return False
    for root in _ALL_BAM_ROOTS:
        try:
            rp.relative_to(root)
            return True
        except ValueError:
            continue
    return False


def _bam_entry(sid: str, batch: Path, path: Path, *, source: str) -> dict:
    index = _bam_index_for(path)
    return {
        "label":     f"{sid} ({batch.name})",
        "sample_id": sid,
        "batch":     batch.name,
        "path":      str(path),
        "index_path": str(index) if index else "",
        "source":    source,
    }


def _manual_bam_entry(path: Path, folder: Path) -> dict:
    sid = _sample_id_from_bam(path)
    index = _bam_index_for(path)
    batch = folder.parent.name if folder.name == "bam" else folder.name
    return {
        "label":     f"{sid} ({batch})",
        "sample_id": sid,
        "batch":     batch,
        "path":      str(path),
        "index_path": str(index) if index else "",
        "source":    "manual",
    }


def _bam_mtime(hit: dict) -> float:
    try:
        return Path(hit["path"]).stat().st_mtime
    except OSError:
        return 0.0


def _bam_hits(sid: str, *, source: str = "") -> list[dict]:
    out: list[dict] = []
    if source in ("", "inhouse", "nckuh"):
        for root in _INHOUSE_BAM_ROOTS:
            if not root.is_dir():
                continue
            for batch in sorted(root.iterdir()):
                if not batch.is_dir():
                    continue
                p = _inhouse_bam_for(sid, batch)
                if p:
                    out.append(_bam_entry(sid, batch, p, source="inhouse"))
    if source in ("", "dragen"):
        for root in _DRAGEN_BAM_ROOTS:
            if not root.is_dir():
                continue
            for run_dir in sorted(root.iterdir()):
                if not run_dir.is_dir():
                    continue
                p = _dragen_bam_for(sid, run_dir)
                if p:
                    out.append(_bam_entry(sid, run_dir, p, source="dragen"))
    out.sort(key=_bam_mtime, reverse=True)
    return out


def _sidecar_for(sid: str) -> dict:
    path = TERTIARY_OUTPUT_ROOT / sid / "pipeline_source.json"
    try:
        data = json.loads(path.read_text(encoding="utf-8")) or {}
    except (OSError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def _source_sid_from_sidecar(sid: str) -> str:
    data = _sidecar_for(sid)
    source_sid = str(data.get("source_sample_id") or "").strip()
    return source_sid if _SID_RE.match(source_sid) else ""


def _pipeline_type_from_sidecar(sid: str) -> str:
    raw = str(_sidecar_for(sid).get("pipeline_type") or "").strip().lower()
    if raw == "dragen":
        return "dragen"
    if raw in ("nckuh", "inhouse"):
        return "inhouse"
    return ""


def _legacy_source_sid(sid: str, *, source: str = "") -> str:
    """Resolve UI sample aliases such as VAL-60-dragen → VAL-60."""
    candidates = [
        sid.removesuffix(suffix)
        for suffix in _LEGACY_ALIAS_SUFFIXES
        if sid.endswith(suffix)
    ]
    hits = [
        hit
        for candidate in candidates
        for hit in _bam_hits(candidate, source=source)
    ]
    return hits[0]["sample_id"] if hits else ""


def _resolve_primary_bam(sid: str) -> tuple[dict | None, str, str]:
    preferred_source = _pipeline_type_from_sidecar(sid)
    direct = _bam_hits(sid, source=preferred_source)
    if direct:
        return direct[0], sid, "exact"
    source_sid = _source_sid_from_sidecar(sid)
    if source_sid:
        hits = _bam_hits(source_sid, source=preferred_source)
        if hits:
            return hits[0], source_sid, "pipeline_source"
    legacy_sid = _legacy_source_sid(sid, source=preferred_source)
    if legacy_sid:
        hits = _bam_hits(legacy_sid, source=preferred_source)
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
        if primary.get("source") == "dragen":
            run_dir = Path(primary["path"]).parents[1]
            bam_dir = run_dir / "bam"
            for sp in sorted(bam_dir.glob("*.bam")) if bam_dir.is_dir() else []:
                if sp.name.endswith(".repeats.bam"):
                    continue
                sib_sid = _sample_id_from_dragen_bam(sp)
                if sib_sid == resolved_sid:
                    continue
                siblings.append(_bam_entry(sib_sid, run_dir, sp, source="dragen"))
                if len(siblings) >= _SIBLING_LIMIT:
                    break
        else:
            batch = Path(primary["path"]).parents[2]
            for sib in sorted(batch.iterdir()):
                if not sib.is_dir() or sib.name == resolved_sid:
                    continue
                sp = _inhouse_bam_for(sib.name, batch)
                if sp:
                    siblings.append(_bam_entry(sib.name, batch, sp, source="inhouse"))
                if len(siblings) >= _SIBLING_LIMIT:
                    break
    return {
        "primary": primary,
        "siblings": siblings,
        "requested_sample_id": sid,
        "resolved_sample_id": resolved_sid,
        "resolution": resolution,
        "roots": [str(r) for r in _ALL_BAM_ROOTS],
        "dragen_roots": [str(r) for r in _DRAGEN_BAM_ROOTS],
        "inhouse_roots": [str(r) for r in _INHOUSE_BAM_ROOTS],
    }


@router.get("/batch-samples")
def list_batch_samples(batch: str = Query(...)):
    """List every sample with a BAM in the given batch — used by the
    "add another BAM" picker so reviewers can pull in any sibling."""
    batch = (batch or "").strip()
    if not re.match(r"^[A-Za-z0-9_.\-]{1,64}$", batch):
        raise HTTPException(400, "invalid batch name")
    out: list[dict] = []
    for root in _INHOUSE_BAM_ROOTS:
        bd = root / batch
        if not bd.is_dir():
            continue
        for d in sorted(bd.iterdir()):
            if not d.is_dir():
                continue
            p = _inhouse_bam_for(d.name, bd)
            if p:
                out.append(_bam_entry(d.name, bd, p, source="inhouse"))
    for root in _DRAGEN_BAM_ROOTS:
        run_dir = root / batch
        bam_dir = run_dir / "bam"
        if not bam_dir.is_dir():
            continue
        for p in sorted(bam_dir.glob("*.bam")):
            if p.name.endswith(".repeats.bam"):
                continue
            sid = _sample_id_from_dragen_bam(p)
            out.append(_bam_entry(sid, run_dir, p, source="dragen"))
    return {"batch": batch, "samples": out}


@router.get("/bam-folder")
def list_bam_folder(dir: str = Query(...)):
    """List BAMs directly under a reviewer-provided server directory.

    The directory must still live under one of the configured BAM roots.
    We intentionally do not recurse; reviewers usually point this at a
    concrete run `bam/` folder, and non-recursive listing avoids
    accidentally sweeping unrelated data.
    """
    folder = Path(dir).expanduser()
    if not folder.is_dir():
        raise HTTPException(404, "folder not found")
    if not _dir_ok(folder):
        raise HTTPException(403, "folder is outside configured BAM roots")
    out = []
    for p in sorted(folder.glob("*.bam")):
        if p.name.endswith(".repeats.bam"):
            continue
        if p.is_file():
            out.append(_manual_bam_entry(p, folder))
    out.sort(key=_bam_mtime, reverse=True)
    return {
        "dir": str(folder),
        "samples": out,
        "roots": [str(r) for r in _ALL_BAM_ROOTS],
    }


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
