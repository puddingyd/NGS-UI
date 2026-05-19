"""DRAGEN pipeline job management.

State lives under DRAGEN_JOBS_DIR/{job_id}/:
    state.json     job metadata + current step (atomically rewritten)
    log.txt        combined stdout/stderr from the chain
    pid            spawned worker PID (for `is_running` check)

Jobs are spawned via subprocess.Popen with start_new_session=True so
they survive a uvicorn reload / restart; we never wait on them
inside the request handler. The frontend polls /api/dragen/jobs/{id}
every few seconds.
"""
from __future__ import annotations

import json
import os
import re
import signal
import subprocess
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

from ..config import (DRAGEN_JOBS_DIR, DRAGEN_VCF_ROOTS, INHOUSE_VCF_ROOTS,
                       PIPELINE_VCF_INDEX_PATH,
                       PIPELINE_VCF_INDEX_TTL_HOURS, REPO_ROOT)

# Final pipeline steps, in order — the worker writes the current one
# into state.json so the UI can show progress.
PIPELINE_STEPS = [
    "queued",
    "mito",
    "stage",
    "nextflow",
    "stop-gaps",
    "done",
]


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# ── VCF discovery ──────────────────────────────────────────────────

_DRAGEN_VCF_GLOBS = [
    "*hard-filtered.vcf.gz",
    "*/vcf.gz/*hard-filtered.vcf.gz",
    "*/*/*hard-filtered.vcf.gz",
]
_SUFFIX_RE = re.compile(r"\.hard-filtered\.vcf\.gz$", re.IGNORECASE)


def list_dragen_vcfs() -> list[dict]:
    """Scan every configured DRAGEN_VCF_ROOTS for hard-filtered VCFs.

    Returns most-recent-first list of
        {path, sample_id, run, size, mtime}.
    `sample_id` is the basename minus the `.hard-filtered.vcf.gz`
    suffix; `run` is the closest parent directory that looks like a
    sequencing-run folder (basename of the dirname containing
    `vcf.gz/` if any, else the immediate parent).
    """
    seen: set[str] = set()
    out: list[dict] = []
    for root in DRAGEN_VCF_ROOTS:
        if not root.exists():
            continue
        for pat in _DRAGEN_VCF_GLOBS:
            for p in root.glob(pat):
                if not p.is_file():
                    continue
                sp = str(p)
                if sp in seen:
                    continue
                seen.add(sp)
                sid = _SUFFIX_RE.sub("", p.name)
                # Locate the run folder: e.g. /datalake/Novaseq/20260428_LH00873/vcf.gz/sample.vcf.gz
                # → run = "20260428_LH00873"
                run = ""
                for parent in p.parents:
                    if parent == root:
                        break
                    if parent.name == "vcf.gz":
                        continue
                    run = parent.name
                    break
                try:
                    st = p.stat()
                except OSError:
                    continue
                out.append({
                    "path": sp,
                    "sample_id": sid,
                    "run": run,
                    "size": st.st_size,
                    "mtime": st.st_mtime,
                })
    out.sort(key=lambda r: r["mtime"], reverse=True)
    return out


_INHOUSE_SNV_REL = "04_snv_indel"
_INHOUSE_SUFFIX_RE = re.compile(r"\.ensemble\.fixed\.vcf\.gz$", re.IGNORECASE)


def _find_inhouse_snv_vcfs(root: Path) -> list[Path]:
    """find(1) is 10-30x faster than pathlib.rglob across the datalake.

    We only care about the SNV anchor; siblings are derived from its
    path. Falls back to Python glob if `find` is missing.
    """
    if not root.is_dir():
        return []
    try:
        proc = subprocess.run(
            ["find", str(root), "-type", "f",
             "-path", f"*/{_INHOUSE_SNV_REL}/*.ensemble.fixed.vcf.gz"],
            check=False, capture_output=True, text=True,
        )
        if proc.returncode == 0:
            return [Path(line) for line in proc.stdout.splitlines() if line]
    except OSError:
        pass
    return list(root.glob(f"**/{_INHOUSE_SNV_REL}/*.ensemble.fixed.vcf.gz"))


def list_inhouse_vcfs() -> list[dict]:
    """Scan INHOUSE_VCF_ROOTS for in-house ensemble Nextflow outputs.

    Anchor on the SNV/Indel VCF, then discover three siblings under the
    sample dir. `run` is the parent-of-sample-dir basename (often a
    batch / study id). Returns most-recent-first.
    """
    seen: set[str] = set()
    out: list[dict] = []
    for root in INHOUSE_VCF_ROOTS:
        for snv in _find_inhouse_snv_vcfs(root):
            sp = str(snv)
            if sp in seen or not snv.is_file():
                continue
            seen.add(sp)
            sid = _INHOUSE_SUFFIX_RE.sub("", snv.name)
            # snv = <root>/.../<run>/<SID>/04_snv_indel/<SID>.ensemble.fixed.vcf.gz
            #       parents:  [0]=04_snv_indel  [1]=<SID>  [2]=<run>
            sample_dir = snv.parents[1] if len(snv.parents) >= 2 else snv.parent
            run = snv.parents[2].name if len(snv.parents) >= 3 else ""

            def sib(rel: str) -> str:
                p = sample_dir / rel
                return str(p) if p.is_file() else ""

            cnv  = sib(f"05_cnv_sv/{sid}.gcnv.vcf.gz")
            sv   = sib(f"05_cnv_sv/{sid}.delly.vcf.gz")
            mito = sib(f"07_mitochondria/{sid}.mito.vcf.gz")

            try:
                st = snv.stat()
            except OSError:
                continue
            out.append({
                "path":        sp,
                "sample_id":   sid,
                "sample_dir":  str(sample_dir),
                "run":         run,
                "cnv_vcf":     cnv,
                "sv_vcf":      sv,
                "mito_vcf":    mito,
                "size":        st.st_size,
                "mtime":       st.st_mtime,
            })
    out.sort(key=lambda r: r["mtime"], reverse=True)
    return out


# ── Pipeline VCF index ─────────────────────────────────────────────

def load_index() -> dict | None:
    p = PIPELINE_VCF_INDEX_PATH
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None


def save_index(idx: dict) -> None:
    p = PIPELINE_VCF_INDEX_PATH
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(idx, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    tmp.replace(p)


def refresh_index() -> dict:
    """Run both scans, persist results, return the new index."""
    t0 = time.time()
    dragen  = list_dragen_vcfs()
    inhouse = list_inhouse_vcfs()
    idx = {
        "updated_at":        _now(),
        "scan_duration_sec": round(time.time() - t0, 2),
        "dragen":            dragen,
        "inhouse":           inhouse,
    }
    save_index(idx)
    return idx


def index_is_stale(idx: dict | None) -> bool:
    if not idx or not idx.get("updated_at"):
        return True
    try:
        ts = datetime.fromisoformat(idx["updated_at"].replace("Z", "+00:00"))
    except ValueError:
        return True
    age_h = (datetime.now(timezone.utc) - ts).total_seconds() / 3600
    return age_h >= PIPELINE_VCF_INDEX_TTL_HOURS


# ── Job state I/O ──────────────────────────────────────────────────

def _job_dir(job_id: str) -> Path:
    return DRAGEN_JOBS_DIR / job_id


def _state_path(job_id: str) -> Path:
    return _job_dir(job_id) / "state.json"


def _log_path(job_id: str) -> Path:
    return _job_dir(job_id) / "log.txt"


def _pid_path(job_id: str) -> Path:
    return _job_dir(job_id) / "pid"


def load_state(job_id: str) -> dict | None:
    p = _state_path(job_id)
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def save_state(job_id: str, state: dict) -> None:
    p = _state_path(job_id)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(state, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    tmp.replace(p)


def _pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except (ProcessLookupError, PermissionError, OSError):
        return False
    return True


def is_running(job_id: str) -> bool:
    st = load_state(job_id)
    if st is None or st.get("state") in ("done", "failed", "cancelled"):
        return False
    pid_file = _pid_path(job_id)
    if not pid_file.is_file():
        return False
    try:
        pid = int(pid_file.read_text().strip())
    except (ValueError, OSError):
        return False
    return _pid_alive(pid)


def tail_log(job_id: str, n: int = 50) -> str:
    p = _log_path(job_id)
    if not p.is_file():
        return ""
    try:
        with p.open("rb") as f:
            try:
                f.seek(-min(p.stat().st_size, 32 * 1024), os.SEEK_END)
            except OSError:
                f.seek(0)
            data = f.read().decode("utf-8", errors="replace")
        lines = data.splitlines()
        return "\n".join(lines[-n:])
    except OSError:
        return ""


def list_jobs(limit: int = 50) -> list[dict]:
    if not DRAGEN_JOBS_DIR.is_dir():
        return []
    jobs: list[dict] = []
    for child in DRAGEN_JOBS_DIR.iterdir():
        if not child.is_dir():
            continue
        st = load_state(child.name)
        if not st:
            continue
        jobs.append(st)
    jobs.sort(key=lambda j: j.get("created_at", ""), reverse=True)
    return jobs[:limit]


# ── Job spawn ──────────────────────────────────────────────────────

def start_job(
    vcf_path: str,
    sample_id: str,
    *,
    mode: str = "dragen",
    with_extra_vep: bool = True,
    cnv_vcf: str = "",
    sv_vcf: str = "",
    mito_vcf: str = "",
) -> str:
    """Spawn a detached worker that runs the chosen pipeline chain.

    mode = "dragen" → DRAGEN germline (single hard-filtered VCF; siblings
                       cnv.vcf.gz / cnv_sv.vcf.gz auto-discovered).
    mode = "inhouse" → in-house ensemble Nextflow output (vcf_path is the
                       ensemble.fixed.vcf.gz; cnv_vcf / sv_vcf / mito_vcf
                       are the explicit sibling paths from the index).

    Returns the job_id. The worker writes state.json + log.txt under
    DRAGEN_JOBS_DIR/<job_id>/; the route polls.
    """
    if mode not in ("dragen", "inhouse"):
        raise ValueError(f"unknown mode: {mode}")
    vcf = Path(vcf_path)
    if not vcf.is_file():
        raise FileNotFoundError(f"VCF not found: {vcf_path}")
    if not sample_id:
        raise ValueError("sample_id required")

    job_id = f"{int(time.time())}-{uuid.uuid4().hex[:8]}"
    jdir = _job_dir(job_id)
    jdir.mkdir(parents=True, exist_ok=True)

    save_state(job_id, {
        "job_id":         job_id,
        "mode":           mode,
        "vcf_path":       str(vcf),
        "sample_id":      sample_id,
        "with_extra_vep": with_extra_vep,
        "cnv_vcf":        cnv_vcf,
        "sv_vcf":         sv_vcf,
        "mito_vcf":       mito_vcf,
        "state":          "queued",
        "step":           "queued",
        "created_at":     _now(),
        "started_at":     None,
        "finished_at":    None,
        "error":          None,
    })

    log_fh = _log_path(job_id).open("w", buffering=1)
    cmd = [
        "python3", "-m", "app.workers.dragen_run",
        "--job-id",  job_id,
        "--vcf",     str(vcf),
        "--sample",  sample_id,
        "--mode",    mode,
    ]
    if with_extra_vep:
        cmd.append("--with-extra-vep")
    if cnv_vcf:  cmd += ["--cnv-vcf",  cnv_vcf]
    if sv_vcf:   cmd += ["--sv-vcf",   sv_vcf]
    if mito_vcf: cmd += ["--mito-vcf", mito_vcf]

    env = os.environ.copy()
    env.setdefault("PYTHONPATH", str(REPO_ROOT / "backend"))

    proc = subprocess.Popen(
        cmd,
        stdout=log_fh, stderr=subprocess.STDOUT,
        cwd=str(REPO_ROOT),
        env=env,
        start_new_session=True,
    )
    _pid_path(job_id).write_text(str(proc.pid))
    return job_id


def cancel_job(job_id: str) -> bool:
    """Best-effort: send SIGTERM to the worker process group."""
    pid_file = _pid_path(job_id)
    if not pid_file.is_file():
        return False
    try:
        pid = int(pid_file.read_text().strip())
    except (ValueError, OSError):
        return False
    try:
        os.killpg(os.getpgid(pid), signal.SIGTERM)
    except (ProcessLookupError, PermissionError, OSError):
        return False
    st = load_state(job_id) or {}
    st.update({"state": "cancelled", "finished_at": _now()})
    save_state(job_id, st)
    return True
