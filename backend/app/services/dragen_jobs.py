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
                       REPO_ROOT, TERTIARY_OUTPUT_ROOT)

# Final pipeline steps, in order — the worker writes the current one
# into state.json so the UI can show progress.
PIPELINE_STEPS = [
    "queued",
    "detect-pipeline-output",
    "samplesheet",
    "stage",  # legacy fallback only
    "nextflow",
    "copy-pipeline-tsv",
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


def infer_source_sample_id(vcf: Path, mode: str) -> str:
    """Return the sequencing sample ID encoded in the selected VCF name."""
    name = vcf.name
    suffix = (".ensemble.fixed.vcf.gz" if mode == "inhouse"
              else ".hard-filtered.vcf.gz")
    sid = name.removesuffix(suffix)
    return _validate_sample_id(sid)


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
                input_dir = str(p.parent.parent if p.parent.name == "vcf.gz" else p.parent)
                try:
                    st = p.stat()
                except OSError:
                    continue
                out.append({
                    "path": sp,
                    "sample_id": sid,
                    "run": run,
                    "input_dir": input_dir,
                    "size": st.st_size,
                    "mtime": st.st_mtime,
                })
    out.sort(key=lambda r: r["mtime"], reverse=True)
    return out


_INHOUSE_SNV_REL = "04_snv_indel"
_INHOUSE_SUFFIX_RE = re.compile(r"\.ensemble\.fixed\.vcf\.gz$", re.IGNORECASE)
_INHOUSE_WGS_SIZE_THRESHOLD = 100 * 1024 * 1024


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


def _infer_seq_type_from_size(size: int | None, *, default: str = "WES") -> str:
    """Use the agreed UI heuristic: in-house VCF >=100 MB defaults to WGS."""
    if size is not None and size >= _INHOUSE_WGS_SIZE_THRESHOLD:
        return "WGS"
    return default


def _normalize_seq_type(value: str | None, *, default: str) -> str:
    v = (value or "").strip().upper()
    if v in {"WES", "WGS"}:
        return v
    return default


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
            seq_type = _infer_seq_type_from_size(st.st_size)
            out.append({
                "path":        sp,
                "sample_id":   sid,
                "sample_dir":  str(sample_dir),
                "input_dir":   str(sample_dir),
                "run":         run,
                "cnv_vcf":     cnv,
                "sv_vcf":      sv,
                "mito_vcf":    mito,
                "seq_type":    seq_type,
                "seq_type_inferred_by": "size",
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

def _job_samples(job: dict) -> list[dict]:
    samples = job.get("samples")
    if isinstance(samples, list) and samples:
        out = []
        for sample in samples:
            if not isinstance(sample, dict):
                continue
            sample_id = sample.get("sample_id", "")
            source_sample_id = sample.get("source_sample_id", "")
            if sample_id:
                out.append({
                    "sample_id": sample_id,
                    "source_sample_id": source_sample_id or sample_id,
                })
        if out:
            return out
    sample_id = job.get("sample_id", "")
    if not sample_id:
        return []
    return [{
        "sample_id": sample_id,
        "source_sample_id": job.get("source_sample_id", sample_id) or sample_id,
    }]


def _latest_job_for_sample(sample_id: str) -> dict | None:
    jobs = [
        j for j in list_jobs(limit=1000)
        if any(
            s["sample_id"] == sample_id or s["source_sample_id"] == sample_id
            for s in _job_samples(j)
        )
    ]
    return jobs[0] if jobs else None


def _job_mtime(job: dict) -> float:
    for key in ("finished_at", "started_at", "created_at"):
        raw = job.get(key)
        if not raw:
            continue
        try:
            return datetime.fromisoformat(str(raw).replace("Z", "+00:00")).timestamp()
        except ValueError:
            continue
    return 0.0


def list_pipeline_outputs() -> list[dict]:
    """List pipeline sample directories plus NGS-UI jobs without output."""
    out: list[dict] = []
    latest_jobs: dict[str, dict] = {}
    jobs_by_source: dict[str, tuple[dict, str]] = {}
    for job in list_jobs(limit=1000):
        for sample in _job_samples(job):
            sample_id = sample["sample_id"]
            source_sample_id = sample["source_sample_id"]
            if sample_id and sample_id not in latest_jobs:
                latest_jobs[sample_id] = job
            if source_sample_id and source_sample_id not in jobs_by_source:
                jobs_by_source[source_sample_id] = (job, sample_id)

    sample_dirs: dict[str, tuple[Path, str]] = {}
    if PIPELINE_OUT_ROOT.is_dir():
        for child in PIPELINE_OUT_ROOT.iterdir():
            if not child.is_dir() or child.name.startswith("_"):
                continue
            try:
                _validate_sample_id(child.name)
            except ValueError:
                continue
            job_entry = jobs_by_source.get(child.name)
            ui_sample_id = job_entry[1] if job_entry else child.name
            sample_dirs[ui_sample_id] = (child, child.name)

    for sample_id in set(sample_dirs) | set(latest_jobs):
        sample_entry = sample_dirs.get(sample_id)
        child = sample_entry[0] if sample_entry else None
        pipeline_sample_id = sample_entry[1] if sample_entry else ""
        job = latest_jobs.get(sample_id)
        try:
            output_mtime = child.stat().st_mtime if child else 0.0
        except OSError:
            output_mtime = 0.0
        acmg_dir = child / "03_acmg" if child else None
        has_acmg = bool(
            acmg_dir and acmg_dir.is_dir()
            and any(acmg_dir.glob("*.snv_indel.acmg.tsv"))
        )
        job_id = (job or {}).get("job_id", "")
        job_sample = None
        if job:
            job_sample = next(
                (s for s in _job_samples(job) if s["sample_id"] == sample_id),
                None,
            )
        out.append({
            "sample_id":     sample_id,
            "source_sample_id": (job_sample or {}).get("source_sample_id", pipeline_sample_id),
            "pipeline_sample_id": pipeline_sample_id,
            "mtime":         max(output_mtime, _job_mtime(job or {})),
            "has_output":    child is not None,
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
    job = _latest_job_for_sample(sample_id)
    if not job:
        sample_dir = PIPELINE_OUT_ROOT / sample_id
        if not sample_dir.is_dir():
            raise FileNotFoundError(f"pipeline output or job not found: {sample_id}")
        return {"sample_id": sample_id, "job_id": "", "log": ""}
    job_id = job.get("job_id", "")
    return {
        "sample_id": sample_id,
        "job_id": job_id,
        "job_state": job.get("state", ""),
        "log": tail_log(job_id, n=max(1, min(n, 2000))),
    }


def delete_pipeline_output(sample_id: str) -> dict:
    """Delete UI/pipeline output and all NGS-UI job logs unless one is active."""
    _validate_sample_id(sample_id)
    jobs = [
        j for j in list_jobs(limit=1000)
        if any(
            s["sample_id"] == sample_id or s["source_sample_id"] == sample_id
            for s in _job_samples(j)
        )
    ]
    primary_job = jobs[0] if jobs else None
    matched_sample = None
    if primary_job:
        matched_sample = next(
            (
                s for s in _job_samples(primary_job)
                if s["sample_id"] == sample_id or s["source_sample_id"] == sample_id
            ),
            None,
        )
    ui_sample_id = (matched_sample or {}).get("sample_id") or sample_id
    pipeline_sample_id = (matched_sample or {}).get("source_sample_id") or sample_id
    sample_dir = PIPELINE_OUT_ROOT / pipeline_sample_id
    ui_dir = TERTIARY_OUTPUT_ROOT / ui_sample_id
    if not sample_dir.is_dir() and not ui_dir.is_dir() and not jobs:
        raise FileNotFoundError(f"pipeline output or job not found: {sample_id}")
    for job in jobs:
        job_id = job.get("job_id", "")
        if job.get("state") in ("queued", "running") and is_running(job_id):
            raise RuntimeError(f"pipeline output is still in use: {sample_id}")
    deleted: list[str] = []
    if sample_dir.is_dir():
        shutil.rmtree(sample_dir)
        deleted.append(str(sample_dir))
    if ui_dir.is_dir() and ui_dir.resolve() != sample_dir.resolve():
        shutil.rmtree(ui_dir)
        deleted.append(str(ui_dir))
    for job in jobs:
        if (job.get("sample_count") or 1) > 1:
            remaining = [
                sample for sample in (job.get("samples") or [])
                if not (
                    sample.get("sample_id") == ui_sample_id
                    or sample.get("source_sample_id") == pipeline_sample_id
                )
            ]
            job_id = job.get("job_id", "")
            if remaining:
                job = dict(job)
                job["samples"] = remaining
                job["sample_ids"] = [s.get("sample_id", "") for s in remaining if s.get("sample_id")]
                job["source_sample_ids"] = [
                    s.get("source_sample_id", "") for s in remaining if s.get("source_sample_id")
                ]
                job["sample_count"] = len(remaining)
                job["sample_id"] = job["sample_ids"][0] if job["sample_ids"] else ""
                job["source_sample_id"] = (
                    job["source_sample_ids"][0] if job["source_sample_ids"] else job["sample_id"]
                )
                save_state(job_id, job)
                deleted.append(f"{_job_dir(job_id)}/state.json sample entry")
            else:
                job_dir = _job_dir(job_id)
                if job_dir.is_dir():
                    shutil.rmtree(job_dir)
                    deleted.append(str(job_dir))
            continue
        job_dir = _job_dir(job.get("job_id", ""))
        if job_dir.is_dir():
            shutil.rmtree(job_dir)
            deleted.append(str(job_dir))
    from . import sample_loader
    sample_loader.invalidate_sample_cache(ui_dir)
    return {"sample_id": ui_sample_id, "source_sample_id": pipeline_sample_id, "deleted": deleted}


# ── Job spawn ──────────────────────────────────────────────────────

def start_job(
    vcf_path: str,
    sample_id: str,
    *,
    source_sample_id: str = "",
    mode: str = "dragen",
    seq_type: str = "",
    with_extra_vep: bool = True,
    cnv_vcf: str = "",
    sv_vcf: str = "",
    mito_vcf: str = "",
    samples: list[dict] | None = None,
) -> str:
    """Spawn a detached worker that runs the chosen pipeline chain.

    mode = "dragen" → DRAGEN germline; v3.1 sample sheet uses
                       {input_dir}/vcf.gz/{source_sample_id}.hard-filtered.vcf.gz.
    mode = "inhouse" → NCKUH ensemble output; v3.x sample sheet uses
                       {input_dir}/04_snv_indel/{source_sample_id}.ensemble.fixed.vcf.gz.
                       seq_type should be WES or WGS. cnv_vcf / sv_vcf
                       are still passed for NGS-UI AnnotSV fallback until
                       pipeline CNV/SV lands. mito_vcf is
                       retained in job metadata only; v3.2 pipeline output
                       04_mito/{sample}.mito.tsv is copied back to the UI.

    Returns the job_id. The worker writes state.json + log.txt under
    TERTIARY_JOBS_DIR/<job_id>/; the route polls.
    """
    batch_samples: list[dict] = []
    if samples is not None:
        if not samples:
            raise ValueError("samples must not be empty")
        for item in samples:
            if not isinstance(item, dict):
                raise ValueError("each sample must be an object")
            sample_mode = (item.get("mode") or mode or "dragen").strip()
            if sample_mode not in ("dragen", "inhouse"):
                raise ValueError(f"unknown mode: {sample_mode}")
            sample_vcf = Path((item.get("vcf_path") or "").strip())
            if not sample_vcf.is_file():
                raise FileNotFoundError(f"VCF not found: {sample_vcf}")
            sample_id_i = (item.get("sample_id") or "").strip()
            if not sample_id_i:
                raise ValueError("sample_id required")
            _validate_sample_id(sample_id_i)
            source_id_i = (item.get("source_sample_id") or "").strip() or infer_source_sample_id(sample_vcf, sample_mode)
            _validate_sample_id(source_id_i)
            default_seq = "WGS" if sample_mode == "dragen" else _infer_seq_type_from_size(sample_vcf.stat().st_size)
            batch_samples.append({
                "mode": sample_mode,
                "vcf_path": str(sample_vcf),
                "sample_id": sample_id_i,
                "source_sample_id": source_id_i,
                "seq_type": _normalize_seq_type(item.get("seq_type"), default=default_seq),
                "cnv_vcf": (item.get("cnv_vcf") or "").strip(),
                "sv_vcf": (item.get("sv_vcf") or "").strip(),
                "mito_vcf": (item.get("mito_vcf") or "").strip(),
            })
        modes = {s["mode"] for s in batch_samples}
        if len(modes) != 1:
            raise ValueError("samples in one job must all use the same mode")
        mode = batch_samples[0]["mode"]
        vcf = Path(batch_samples[0]["vcf_path"])
        sample_id = batch_samples[0]["sample_id"]
        source_sample_id = batch_samples[0]["source_sample_id"]
    else:
        if mode not in ("dragen", "inhouse"):
            raise ValueError(f"unknown mode: {mode}")
        vcf = Path(vcf_path)
        if not vcf.is_file():
            raise FileNotFoundError(f"VCF not found: {vcf_path}")
        if not sample_id:
            raise ValueError("sample_id required")
        _validate_sample_id(sample_id)
        source_sample_id = source_sample_id or infer_source_sample_id(vcf, mode)
        _validate_sample_id(source_sample_id)
        default_seq = "WGS" if mode == "dragen" else _infer_seq_type_from_size(vcf.stat().st_size)
        batch_samples = [{
            "mode": mode,
            "vcf_path": str(vcf),
            "sample_id": sample_id,
            "source_sample_id": source_sample_id,
            "seq_type": _normalize_seq_type(seq_type, default=default_seq),
            "cnv_vcf": cnv_vcf,
            "sv_vcf": sv_vcf,
            "mito_vcf": mito_vcf,
        }]

    job_id = f"{int(time.time())}-{uuid.uuid4().hex[:8]}"
    jdir = _job_dir(job_id)
    jdir.mkdir(parents=True, exist_ok=True)

    created_at = _now()
    sample_ids = [s["sample_id"] for s in batch_samples]
    source_sample_ids = [s["source_sample_id"] for s in batch_samples]
    save_state(job_id, {
        "job_id":         job_id,
        "mode":           mode,
        "vcf_path":       str(vcf),
        "sample_id":      sample_id,
        "source_sample_id": source_sample_id,
        "sample_ids":     sample_ids,
        "source_sample_ids": source_sample_ids,
        "sample_count":   len(batch_samples),
        "samples":        batch_samples,
        "pipeline_type":  "dragen" if mode == "dragen" else "nckuh",
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
        "nextflow_step_history": [],
    })

    log_fh = _log_path(job_id).open("w", buffering=1)
    batch_path = _job_dir(job_id) / "batch.json"
    batch_path.write_text(
        json.dumps({"samples": batch_samples}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    cmd = [
        "python3", "-m", "app.workers.dragen_run",
        "--job-id",  job_id,
        "--batch-json", str(batch_path),
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
