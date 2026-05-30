"""DRAGEN pipeline job management.

State lives under TERTIARY_JOBS_DIR/{job_id}/:
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
import shutil
import signal
import subprocess
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

from ..config import (DRAGEN_VCF_ROOTS, INHOUSE_VCF_ROOTS,
                       LEGACY_DRAGEN_JOBS_DIR, TERTIARY_JOBS_DIR,
                       PIPELINE_VCF_INDEX_PATH,
                       PIPELINE_VCF_INDEX_TTL_HOURS, PIPELINE_OUT_ROOT,
                       REPO_ROOT)

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
_SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _validate_sample_id(sample_id: str) -> str:
    if not _SAMPLE_ID_RE.match(sample_id or ""):
        raise ValueError(
            "sample_id must match [A-Za-z0-9][A-Za-z0-9._-]{0,127}"
        )
    return sample_id


# ── VCF discovery ──────────────────────────────────────────────────

_DRAGEN_VCF_GLOBS = [
    # The only layout we care about now: <root>/<run>/vcf.gz/*hard-
    # filtered.vcf.gz. Restricting to this one pattern keeps the index
    # tight (no random *.hard-filtered.vcf.gz dropped one level up).
    "*/vcf.gz/*hard-filtered.vcf.gz",
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
    current = TERTIARY_JOBS_DIR / job_id
    if current.exists():
        return current
    legacy = LEGACY_DRAGEN_JOBS_DIR / job_id
    return legacy if legacy.exists() else current


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
    jobs: list[dict] = []
    seen: set[str] = set()
    for root in (TERTIARY_JOBS_DIR, LEGACY_DRAGEN_JOBS_DIR):
        if not root.is_dir():
            continue
        for child in root.iterdir():
            if not child.is_dir() or child.name in seen:
                continue
            seen.add(child.name)
            st = load_state(child.name)
            if st:
                jobs.append(st)
    jobs.sort(key=lambda j: j.get("created_at", ""), reverse=True)
    return jobs[:limit]


# ── Finished pipeline output management ───────────────────────────

def _latest_job_for_sample(sample_id: str) -> dict | None:
    jobs = [j for j in list_jobs(limit=1000) if j.get("sample_id") == sample_id]
    return jobs[0] if jobs else None


def list_pipeline_outputs() -> list[dict]:
    """List sample directories under /home/pipeline/tertiary_output."""
    out: list[dict] = []
    if not PIPELINE_OUT_ROOT.is_dir():
        return out
    latest_jobs: dict[str, dict] = {}
    for job in list_jobs(limit=1000):
        sample_id = job.get("sample_id", "")
        if sample_id and sample_id not in latest_jobs:
            latest_jobs[sample_id] = job
    for child in PIPELINE_OUT_ROOT.iterdir():
        if not child.is_dir() or child.name.startswith("_"):
            continue
        try:
            _validate_sample_id(child.name)
            st = child.stat()
        except (ValueError, OSError):
            continue
        acmg_dir = child / "03_acmg"
        has_acmg = acmg_dir.is_dir() and any(acmg_dir.glob("*.snv_indel.acmg.tsv"))
        job = latest_jobs.get(child.name)
        job_id = (job or {}).get("job_id", "")
        out.append({
            "sample_id":     child.name,
            "mtime":         st.st_mtime,
            "has_acmg":      has_acmg,
            "job_id":        job_id,
            "job_state":     (job or {}).get("state", ""),
            "log_available": bool(job_id and _log_path(job_id).is_file()),
        })
    out.sort(key=lambda row: row["mtime"], reverse=True)
    return out


def get_pipeline_output_log(sample_id: str, n: int = 400) -> dict:
    """Return the most recent NGS-UI job log associated with a sample."""
    _validate_sample_id(sample_id)
    sample_dir = PIPELINE_OUT_ROOT / sample_id
    if not sample_dir.is_dir():
        raise FileNotFoundError(f"pipeline output not found: {sample_id}")
    job = _latest_job_for_sample(sample_id)
    if not job:
        return {"sample_id": sample_id, "job_id": "", "log": ""}
    job_id = job.get("job_id", "")
    return {
        "sample_id": sample_id,
        "job_id": job_id,
        "job_state": job.get("state", ""),
        "log": tail_log(job_id, n=max(1, min(n, 2000))),
    }


def delete_pipeline_output(sample_id: str) -> dict:
    """Delete one pipeline output directory unless its job is active."""
    _validate_sample_id(sample_id)
    sample_dir = PIPELINE_OUT_ROOT / sample_id
    if not sample_dir.is_dir():
        raise FileNotFoundError(f"pipeline output not found: {sample_id}")
    job = _latest_job_for_sample(sample_id)
    if job and job.get("state") in ("queued", "running") and is_running(job.get("job_id", "")):
        raise RuntimeError(f"pipeline output is still in use: {sample_id}")
    shutil.rmtree(sample_dir)
    return {"sample_id": sample_id, "deleted": str(sample_dir)}


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
    TERTIARY_JOBS_DIR/<job_id>/; the route polls.
    """
    if mode not in ("dragen", "inhouse"):
        raise ValueError(f"unknown mode: {mode}")
    vcf = Path(vcf_path)
    if not vcf.is_file():
        raise FileNotFoundError(f"VCF not found: {vcf_path}")
    if not sample_id:
        raise ValueError("sample_id required")
    _validate_sample_id(sample_id)

    job_id = f"{int(time.time())}-{uuid.uuid4().hex[:8]}"
    jdir = _job_dir(job_id)
    jdir.mkdir(parents=True, exist_ok=True)

    created_at = _now()
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
        "created_at":     created_at,
        "started_at":     None,
        "finished_at":    None,
        "error":          None,
        "step_started_at": created_at,
        "step_history":   [{"step": "queued", "started_at": created_at}],
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
