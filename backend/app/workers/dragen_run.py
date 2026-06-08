"""Pipeline worker — spawned by dragen_jobs.start_job().

Runs the tertiary chain end-to-end on either a DRAGEN hard-filtered
VCF (`--mode dragen`) or an in-house ensemble Nextflow output
(`--mode inhouse`). Steps:

    1. detect existing pipeline output
       ├ HIT  → skip Nextflow, reuse the pre-existing .acmg.tsv
       └ MISS → run 2+3 below
    2. write a v3.x samplesheet       → data/jobs/tertiary/<job>/samplesheet.csv
       (legacy env fallback can still stage into nf_stage/<SID>/04_snv_indel)
    3. nextflow main_tertiary.nf      → /home/pipeline/tertiary_output/<source SID>/
                                          03_acmg/<source SID>.snv_indel.acmg.tsv
                                          04_mito/<source SID>.mito.tsv
    4. copy pipeline outputs          → NGS_UI/tertiary_output/<SID>/
                                          snv_indel.annotated.tsv
                                          mito.annotated.tsv
                                          + pipeline_source.json (audit)
    5. run_stopgaps.sh                → filter / GeneBe / extra-VEP / CNV-AnnotSV
                                          + pre-build snv_indel.review.tsv
                                          (ClinVar removed — pipeline already
                                           does it; GeneBe writes a SECOND
                                           opinion to GENEBE_* columns)

Mode differences:
  dragen  — sample sheet uses pipeline_type=dragen and the DRAGEN run
            folder as input_dir; CNV/SV AnnotSV auto-discovers sibling
            <SID>.cnv.vcf.gz + <SID>.sv.vcf.gz from the same dir.
  inhouse — sample sheet uses pipeline_type=nckuh and the in-house sample
            output folder as input_dir; CNV/SV AnnotSV runs separately on
            gcnv + delly VCFs.

Existing-output detection (step 2): tries
    /home/pipeline/tertiary_output/<SID>/03_acmg/<SID>.snv_indel.acmg.tsv
    /home/pipeline/tertiary_output/<stripped SID>/03_acmg/<stripped SID>.snv_indel.acmg.tsv
in that order, where <stripped SID> drops the -dragen / -inhouse
suffix. WES/WGS suffixes are test-type identifiers and must not be
stripped, otherwise a WGS run such as VAL-24-WGS could accidentally
reuse an older VAL-24 WES TSV.

Started by `python3 -m app.workers.dragen_run --job-id … --vcf …`.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
import traceback
from datetime import datetime, timedelta, timezone
from pathlib import Path

from ..config import (NGS_UI_HOME, PIPELINE_OUT_ROOT, REPO_ROOT,
                       TERTIARY_JOBS_DIR, TERTIARY_OUTPUT_ROOT)
from ..services import dragen_jobs, mitomap_mito

TERTIARY_NEXTFLOW_CONFIG = Path(os.environ.get(
    "NGS_UI_TERTIARY_CONFIG",
    "/home/n102968/NGS_UI/nextflow_tertiary_no_pgx.config",
))


def _strip_sid_suffix(sid: str) -> str:
    """Drop the -dragen / -inhouse caller suffix the GUI adds for
    directory disambiguation. Returns sid unchanged when no suffix
    matches."""
    for suf in ("-dragen", "-inhouse"):
        if sid.endswith(suf):
            return sid[: -len(suf)]
    return sid


def _find_pipeline_acmg_tsv(*sample_ids: str) -> Path | None:
    """Look for a v3.1 <SID>.snv_indel.acmg.tsv under production output.

    The UI sample ID may intentionally differ from the source VCF sample
    ID (for example, adding "-dragen" to avoid collisions). The pipeline
    itself must run under the source sample ID because v3.1 composes input
    paths from sample_id + input_dir, so try source IDs first and then the
    UI ID / stripped legacy aliases.
    """
    candidates: list[str] = []
    for sid in sample_ids:
        if not sid:
            continue
        for s in (sid, _strip_sid_suffix(sid)):
            if s and s not in candidates:
                candidates.append(s)
    for s in candidates:
        p = PIPELINE_OUT_ROOT / s / "03_acmg" / f"{s}.snv_indel.acmg.tsv"
        if p.is_file():
            return p
    return None


def _find_pipeline_mito_tsv(*sample_ids: str) -> Path | None:
    """Look for a v3.2 <SID>.mito.tsv under production output."""
    candidates: list[str] = []
    for sid in sample_ids:
        if not sid:
            continue
        for s in (sid, _strip_sid_suffix(sid)):
            if s and s not in candidates:
                candidates.append(s)
    for s in candidates:
        p = PIPELINE_OUT_ROOT / s / "04_mito" / f"{s}.mito.tsv"
        if p.is_file():
            return p
    return None


def _find_pipeline_annotsv_tsv(kind: str, *sample_ids: str) -> Path | None:
    """Look for v3.2 06_cnv_sv/<SID>.<kind>.annotated.tsv output."""
    candidates: list[str] = []
    for sid in sample_ids:
        if not sid:
            continue
        for s in (sid, _strip_sid_suffix(sid)):
            if s and s not in candidates:
                candidates.append(s)
    for s in candidates:
        p = PIPELINE_OUT_ROOT / s / "06_cnv_sv" / f"{s}.{kind}.annotated.tsv"
        if p.is_file():
            return p
    return None


def _pipeline_input_dir(vcf: Path, mode: str) -> Path:
    """Return the sample-sheet input_dir for the v3.1 pipeline."""
    if mode == "dragen":
        # guide v3.1: {input_dir}/vcf.gz/{sample_id}.hard-filtered.vcf.gz
        return vcf.parent.parent if vcf.parent.name == "vcf.gz" else vcf.parent
    # guide v3.1: {input_dir}/04_snv_indel/{sample_id}.ensemble.fixed.vcf.gz
    if vcf.parent.name == "04_snv_indel" and len(vcf.parents) >= 2:
        return vcf.parents[1]
    return vcf.parent


def _write_samplesheet(
    job_id: str,
    *,
    samples: list[dict],
) -> Path:
    """Create a v3.1 sample sheet for the selected source VCFs."""
    path = TERTIARY_JOBS_DIR / job_id / "samplesheet.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["sample_id", "pipeline_type", "input_dir", "seq_type", "hpo"],
        )
        writer.writeheader()
        for sample in samples:
            mode = sample["mode"]
            writer.writerow({
                "sample_id": sample["source_sample_id"],
                "pipeline_type": "dragen" if mode == "dragen" else "nckuh",
                "input_dir": str(_pipeline_input_dir(Path(sample["vcf_path"]), mode)),
                "seq_type": sample["seq_type"],
                "hpo": "",
            })
    return path


def _validate_acmg_tsv(path: Path, *, strict_v31: bool) -> None:
    """Catch stale/partial tertiary outputs before copying them into UI."""
    try:
        header = path.open(encoding="utf-8").readline().rstrip("\n").split("\t")
    except OSError as e:
        raise RuntimeError(f"cannot read pipeline TSV: {path}") from e
    required = {"HGNC_ID", "ACMG_CRITERIA", "ACMG_SCORE", "ACMG_CLASS", "ACMG_NOTES"}
    missing = sorted(required - set(header))
    if missing:
        raise RuntimeError(
            f"pipeline TSV missing v3.1 required columns: {', '.join(missing)}"
        )
    if strict_v31 and len(header) < 65:
        raise RuntimeError(
            f"pipeline TSV has {len(header)} columns; v3.1 guide expects at least 65"
        )
    if len(header) != 65:
        _log(f"[validate] warning: pipeline TSV has {len(header)} columns (v3.1 guide: 65)")


def _sample_from_args(args: argparse.Namespace) -> dict:
    mode = args.mode
    return {
        "mode": mode,
        "vcf_path": args.vcf,
        "sample_id": args.sample,
        "source_sample_id": args.source_sample,
        "seq_type": "WGS" if mode == "dragen" else "WES",
        "cnv_vcf": args.cnv_vcf,
        "sv_vcf": args.sv_vcf,
        "mito_vcf": args.mito_vcf,
    }


def _normalize_seq_type(value: str | None, *, default: str) -> str:
    v = (value or "").strip().upper()
    if v in {"WES", "WGS"}:
        return v
    return default


def _load_samples(args: argparse.Namespace) -> list[dict]:
    if args.batch_json:
        raw = json.loads(Path(args.batch_json).read_text(encoding="utf-8"))
        samples = raw.get("samples") if isinstance(raw, dict) else raw
        if not isinstance(samples, list) or not samples:
            raise ValueError("--batch-json must contain a non-empty samples list")
    else:
        samples = [_sample_from_args(args)]

    out: list[dict] = []
    for sample in samples:
        if not isinstance(sample, dict):
            raise ValueError("batch sample must be an object")
        mode = (sample.get("mode") or args.mode or "dragen").strip()
        if mode not in ("dragen", "inhouse"):
            raise ValueError(f"unknown mode in batch: {mode}")
        sample_id = (sample.get("sample_id") or "").strip()
        source_sample_id = (sample.get("source_sample_id") or "").strip()
        vcf_path = (sample.get("vcf_path") or sample.get("vcf") or "").strip()
        if not sample_id or not source_sample_id or not vcf_path:
            raise ValueError("batch sample requires sample_id, source_sample_id, and vcf_path")
        out.append({
            "mode": mode,
            "vcf_path": vcf_path,
            "sample_id": sample_id,
            "source_sample_id": source_sample_id,
            "seq_type": _normalize_seq_type(
                sample.get("seq_type"),
                default=("WGS" if mode == "dragen" else "WES"),
            ),
            "cnv_vcf": (sample.get("cnv_vcf") or "").strip(),
            "sv_vcf": (sample.get("sv_vcf") or "").strip(),
            "mito_vcf": (sample.get("mito_vcf") or "").strip(),
        })
    modes = {s["mode"] for s in out}
    if len(modes) != 1:
        raise ValueError("batch samples must all use the same mode/pipeline_type")
    return out


def _track_pipeline_source(
    sample_dir: Path,
    source: Path,
    *,
    source_sample_id: str,
    source_vcf_path: str,
    pipeline_type: str,
) -> None:
    """Write a small audit record so the reviewer (and a future
    re-sync endpoint) can tell where the SNV TSV originated.
    Lives alongside sample_metadata.json so register() doesn't
    accidentally clobber it.
    """
    try:
        mtime = source.stat().st_mtime
    except OSError:
        mtime = None
    rec = {
        "source_path":  str(source),
        "source_mtime": mtime,
        "source_sample_id": source_sample_id,
        "source_vcf_path": source_vcf_path,
        "pipeline_type": pipeline_type,
        "copied_at":    _now(),
    }
    (sample_dir / "pipeline_source.json").write_text(
        json.dumps(rec, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def _symlink_inhouse_into_nf_stage(vcf: Path, sid: str, nf_stage: Path) -> Path:
    """In-house ensemble VCFs are already 2-sample (DV+HC) and don't
    need staging — but Nextflow's input_dir convention forces a layout
    `<input_dir>/04_snv_indel/<sample_id>.ensemble.fixed.vcf.gz`. Drop
    a symlink at that path so Nextflow sees the original VCF
    unmodified (no gnomAD/BED pre-filter — let the pipeline see every
    variant). Returns the staging directory.
    """
    stage_snv = nf_stage / "04_snv_indel"
    stage_snv.mkdir(parents=True, exist_ok=True)
    stage_vcf = stage_snv / f"{sid}.ensemble.fixed.vcf.gz"
    stage_tbi = stage_snv / f"{sid}.ensemble.fixed.vcf.gz.tbi"
    # Idempotent: remove stale links from a prior run.
    for p in (stage_vcf, stage_tbi):
        if p.is_symlink() or p.exists():
            p.unlink()
    stage_vcf.symlink_to(vcf)
    tbi = vcf.with_suffix(vcf.suffix + ".tbi")  # .vcf.gz.tbi
    if tbi.is_file():
        stage_tbi.symlink_to(tbi)
    return nf_stage


TAIPEI_TZ = timezone(timedelta(hours=8))


def _now():
    return datetime.now(TAIPEI_TZ).isoformat(timespec="seconds")


def _load_secrets() -> None:
    """Populate os.environ from $NGS_UI_HOME/secrets.env if present.

    uvicorn runs under systemd and doesn't inherit interactive shell
    `export`s, so subprocess steps that need GENEBE_USER / GENEBE_API_KEY
    fail unless they come from somewhere outside the repo. The file is
    plain KEY=VAL lines (no quoting, no expansion), git-ignored, mode
    0600 — populated once by the operator. Values already in
    os.environ win (systemd Environment= can still override).
    """
    path = NGS_UI_HOME / "secrets.env"
    if not path.is_file():
        return
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        k = k.strip()
        v = v.strip().strip('"').strip("'")
        if k and k not in os.environ:
            os.environ[k] = v


def _update(job_id: str, **kw) -> None:
    st = dragen_jobs.load_state(job_id) or {}
    st.update(kw)
    dragen_jobs.save_state(job_id, st)


def _log(message: str = "") -> None:
    """Write one worker-owned log line with an ISO timestamp."""
    prefix = f"[{_now()}]"
    print(f"{prefix} {message}" if message else prefix, flush=True)


def _set_step(job_id: str, step: str, **kw) -> None:
    """Persist and log step transitions for post-run timing analysis."""
    now = _now()
    st = dragen_jobs.load_state(job_id) or {}
    history = list(st.get("step_history") or [])
    history.append({"step": step, "started_at": now})
    st.update(kw)
    st.update(step=step, step_started_at=now, step_history=history)
    dragen_jobs.save_state(job_id, st)
    _log(f"[step] {step}")


def _record_nextflow_step(
    job_id: str,
    slug: str,
    process: str,
    event: str,
    *,
    elapsed: float | None = None,
) -> None:
    now = _now()
    st = dragen_jobs.load_state(job_id) or {}
    history = list(st.get("nextflow_step_history") or [])
    item = {
        "step": slug,
        "process": process,
        "event": event,
        "at": now,
    }
    if elapsed is not None:
        item["elapsed_seconds"] = round(elapsed, 1)
    history.append(item)
    st["nextflow_step_history"] = history
    dragen_jobs.save_state(job_id, st)
    if elapsed is None:
        _log(f"[nextflow-step] {slug} {event} process={process}")
    else:
        _log(f"[nextflow-step] {slug} {event} process={process} elapsed={elapsed:.1f}s")


def _run(cmd: list[str], *, label: str, on_line=None) -> None:
    """Stream a subprocess's stdout/stderr into this worker's stdout
    (which is already redirected to log.txt by dragen_jobs.start_job).
    Raises on non-zero exit so the outer try/except records failure.
    """
    started = time.monotonic()
    _log()
    _log(f"========================= [{label}] =========================")
    _log("$ " + " ".join(cmd))
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
        bufsize=1,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        print(line, end="", flush=True)
        if on_line is not None:
            on_line(line)
    proc.wait()
    elapsed = time.monotonic() - started
    _log(f"[command] {label} finished exit={proc.returncode} elapsed={elapsed:.1f}s")
    if proc.returncode != 0:
        raise RuntimeError(f"{label} failed (exit {proc.returncode})")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--job-id", required=True)
    ap.add_argument("--vcf",    default="")
    ap.add_argument("--sample", default="")
    ap.add_argument("--source-sample", default="")
    ap.add_argument("--batch-json", default="")
    ap.add_argument("--mode",   default="dragen", choices=["dragen", "inhouse"])
    ap.add_argument("--with-extra-vep", action="store_true")
    # In-house only — explicit sibling VCF paths from the index.
    ap.add_argument("--cnv-vcf",  default="")
    ap.add_argument("--sv-vcf",   default="")
    ap.add_argument("--mito-vcf", default="")
    args = ap.parse_args()

    _load_secrets()

    job_id = args.job_id
    samples = _load_samples(args)
    mode = samples[0]["mode"]
    sample_ids = [sample["sample_id"] for sample in samples]

    scripts = REPO_ROOT / "scripts"
    batch_slug = sample_ids[0] if len(sample_ids) == 1 else f"batch-{job_id}"
    nf_work  = NGS_UI_HOME / "nf_work" / batch_slug
    pipeline_type = "dragen" if mode == "dragen" else "nckuh"
    legacy_staging = os.environ.get("NGS_UI_TERTIARY_LEGACY_STAGING", "").strip().lower() in {
        "1", "true", "yes", "y", "on"
    }

    started_at = _now()
    _set_step(job_id, "detect-pipeline-output", state="running", started_at=started_at)
    try:
        # 1. Reuse existing pipeline output if the production pipeline
        # has already processed each source sample.
        existing_by_sid: dict[str, Path] = {}
        pending_samples: list[dict] = []
        for sample in samples:
            sid = sample["sample_id"]
            source_sid = sample["source_sample_id"]
            existing = _find_pipeline_acmg_tsv(source_sid, sid)
            if existing is not None:
                existing_by_sid[sid] = existing
                _log(f"[detect] {sid}: reusing existing pipeline TSV: {existing}")
            else:
                pending_samples.append(sample)
                _log(f"[detect] {sid}: no existing pipeline TSV")

        if pending_samples:
            _log(f"[detect] running nextflow for {len(pending_samples)} sample(s)")

            # 2. v3.x pipeline input. The official pipeline now owns
            # PREPARE_VCF / PREPARE_VCF_DRAGEN, so the default path is
            # a one-row sample sheet. Keep the old staging route behind
            # an env switch for deployments that have not upgraded yet.
            if legacy_staging:
                if len(pending_samples) != 1:
                    raise RuntimeError(
                        "batch jobs require the v3.x samplesheet pipeline; "
                        "unset NGS_UI_TERTIARY_LEGACY_STAGING or run one sample at a time"
                    )
                sample = pending_samples[0]
                sid = sample["sample_id"]
                source_sid = sample["source_sample_id"]
                vcf = sample["vcf_path"]
                vcf_path = Path(vcf)
                nf_stage = NGS_UI_HOME / "nf_stage" / sid
                _set_step(job_id, "stage")
                _log("[stage] legacy staging enabled by NGS_UI_TERTIARY_LEGACY_STAGING")
                if mode == "inhouse":
                    if source_sid == sid:
                        _symlink_inhouse_into_nf_stage(vcf_path, sid, nf_stage)
                        _log(f"[stage] in-house symlink → {nf_stage}/04_snv_indel/"
                             f"{sid}.ensemble.fixed.vcf.gz (no filter)")
                    else:
                        _run([str(scripts / "stage_dragen_for_tertiary.sh"),
                              "--in",         vcf,
                              "--sample",     sid,
                              "--skip-norm",
                              "--skip-bed",
                              "--skip-gnomad",
                              "--keep-chrm"],
                             label="2a/4 stage in-house alias")
                else:
                    _run([str(scripts / "stage_dragen_for_tertiary.sh"),
                          "--in",     vcf,
                          "--sample", sid],
                         label="2a/4 stage")
            else:
                _set_step(job_id, "samplesheet")
                samplesheet = _write_samplesheet(
                    job_id,
                    samples=pending_samples,
                )
                _log(f"[samplesheet] {samplesheet}")
                for sample in pending_samples:
                    _log(
                        "[samplesheet] "
                        f"source_sample_id={sample['source_sample_id']} "
                        f"ui_sample_id={sample['sample_id']} "
                        f"pipeline_type={pipeline_type} "
                        f"input_dir={_pipeline_input_dir(Path(sample['vcf_path']), mode)}"
                    )

            # 3. Nextflow → /home/pipeline/tertiary_output/<source SID>/...
            _set_step(job_id, "nextflow")
            nextflow_stages = [
                ("prepare-vcf-dragen", "PREPARE_VCF_DRAGEN", 0),
                ("prepare-vcf", "PREPARE_VCF", 0),
                ("add-callers-tag", "ADD_CALLERS_TAG", 0),
                ("filter-for-annotation", "FILTER_FOR_ANNOTATION", 1),
                ("mito-vep", "MITO_VEP", None),
                ("mito-parse", "MITO_PARSE", None),
                ("str-parse-dragen", "STR_PARSE_DRAGEN", None),
                ("str-parse", "STR_PARSE", None),
                ("vep-annotate", "VEP_ANNOTATE", 2),
                ("pangolin-score", "PANGOLIN_SCORE", 3),
                ("parse-csq", "PARSE_CSQ", 4),
                ("acmg-classify", "ACMG_CLASSIFY", 5),
                ("prepare-cnv-dragen", "PREPARE_CNV_DRAGEN", None),
                ("annotsv-cnv-dragen", "ANNOTSV_CNV_DRAGEN", None),
                ("prepare-sv-dragen", "PREPARE_SV_DRAGEN", None),
                ("annotsv-sv-dragen", "ANNOTSV_SV_DRAGEN", None),
            ]
            nextflow_progress_rank = -1
            nextflow_running: dict[str, float] = {}
            nextflow_seen_events: set[tuple[str, str]] = set()

            def track_nextflow(line: str) -> None:
                nonlocal nextflow_progress_rank
                if not re.search(r"\|\s*\d+\s+of\s+\d+", line):
                    return
                match = re.search(r"\|\s*(\d+)\s+of\s+(\d+)", line)
                if not match:
                    return
                done = int(match.group(1))
                total = int(match.group(2))
                for slug, token, progress_rank in nextflow_stages:
                    if token not in line:
                        continue
                    if done < total:
                        event = "start"
                    else:
                        event = "done"
                    key = (slug, event)
                    if key in nextflow_seen_events:
                        return
                    nextflow_seen_events.add(key)
                    process = token
                    if event == "start":
                        nextflow_running[slug] = time.monotonic()
                        _record_nextflow_step(job_id, slug, process, "start")
                        if progress_rank is not None and progress_rank > nextflow_progress_rank:
                            nextflow_progress_rank = progress_rank
                            _set_step(job_id, f"nextflow:{slug}")
                    else:
                        started = nextflow_running.get(slug)
                        elapsed = (time.monotonic() - started) if started is not None else None
                        _record_nextflow_step(job_id, slug, process, "done", elapsed=elapsed)
                        if slug == "vep-annotate" and nextflow_progress_rank < 6:
                            nextflow_progress_rank = 6
                            _set_step(job_id, "nextflow:vep-annotate-done")
                    return

            if legacy_staging:
                sample = pending_samples[0]
                sid = sample["sample_id"]
                nextflow_cmd = [
                    "nextflow",
                    "-c", str(TERTIARY_NEXTFLOW_CONFIG),
                    "run", "/home/pipeline/tertiary_code/main_tertiary.nf",
                    "-profile", "dgm",
                    "-work-dir", str(nf_work),
                    "--sample_id", sid,
                    "--input_dir", str(nf_stage),
                    "--seq_type",  sample["seq_type"],
                    "--out_dir",   str(PIPELINE_OUT_ROOT),
                    "-resume",
                ]
                _run(nextflow_cmd, label="2b/4 nextflow legacy", on_line=track_nextflow)
            else:
                samplesheet = TERTIARY_JOBS_DIR / job_id / "samplesheet.csv"
                inner = " ".join(shlex.quote(part) for part in [
                    "nextflow",
                    "-c", str(TERTIARY_NEXTFLOW_CONFIG),
                    "run", "/home/pipeline/tertiary_code/main_tertiary.nf",
                    "-profile", "dgm",
                    "-work-dir", str(nf_work),
                    "--pipeline_type", pipeline_type,
                    "--samplesheet", str(samplesheet),
                    "--out_dir", str(PIPELINE_OUT_ROOT),
                    "-resume",
                ])
                env_script = os.environ.get(
                    "NGS_UI_TERTIARY_ENV_SCRIPT",
                    "/home/pipeline/pipeline_code/NGS2ndAnalysis_env.sh",
                )
                shell_cmd = (
                    f"if [ -f {shlex.quote(env_script)} ]; then "
                    f"source {shlex.quote(env_script)}; "
                    "else echo '[nextflow] env script not found; using current environment'; fi; "
                    f"{inner}"
                )
                _run([
                    "bash",
                    "-lc",
                    shell_cmd,
                ], label="2b/4 nextflow v3.x", on_line=track_nextflow)

            for sample in pending_samples:
                sid = sample["sample_id"]
                source_sid = sample["source_sample_id"]
                existing = _find_pipeline_acmg_tsv(source_sid, sid)
                if existing is None:
                    raise RuntimeError(
                        "nextflow finished but expected acmg.tsv not found under "
                        f"{PIPELINE_OUT_ROOT}/{source_sid}/03_acmg/"
                    )
                existing_by_sid[sid] = existing

        # 4. Copy pipeline outputs → NGS-UI side; record source audit.
        _set_step(job_id, "copy-pipeline-tsv")
        pipeline_annotsv_copied: dict[str, set[str]] = {}
        for sample in samples:
            sid = sample["sample_id"]
            source_sid = sample["source_sample_id"]
            source_vcf = sample["vcf_path"]
            sample_dir = TERTIARY_OUTPUT_ROOT / sid
            sample_dir.mkdir(parents=True, exist_ok=True)
            gui_tsv = sample_dir / "snv_indel.annotated.tsv"
            existing = existing_by_sid.get(sid)
            if existing is None:
                raise RuntimeError(f"internal error: no pipeline TSV for {sid}")
            _validate_acmg_tsv(existing, strict_v31=not legacy_staging)
            shutil.copyfile(existing, gui_tsv)
            _track_pipeline_source(
                sample_dir,
                existing,
                source_sample_id=source_sid,
                source_vcf_path=source_vcf,
                pipeline_type=mode,
            )
            _log(f"[copy] {sid}: {existing} → {gui_tsv}")
            mito_src = _find_pipeline_mito_tsv(source_sid, sid)
            if mito_src is None:
                _log(
                    f"[copy] {sid}: mito output not found under "
                    f"{PIPELINE_OUT_ROOT}/{source_sid}/04_mito/"
                )
            else:
                mito_dst = sample_dir / "mito.annotated.tsv"
                shutil.copyfile(mito_src, mito_dst)
                try:
                    mitomap_mito.annotate_mito_tsv(mito_dst)
                except Exception:
                    pass
                _log(f"[copy] {sid}: {mito_src} → {mito_dst}")
            for kind in ("cnv", "sv"):
                annotsv_src = _find_pipeline_annotsv_tsv(kind, source_sid, sid)
                if annotsv_src is None:
                    _log(
                        f"[copy] {sid}: {kind.upper()} AnnotSV output not found under "
                        f"{PIPELINE_OUT_ROOT}/{source_sid}/06_cnv_sv/"
                    )
                    continue
                annotsv_dst = sample_dir / f"{kind}.annotated.tsv"
                shutil.copyfile(annotsv_src, annotsv_dst)
                pipeline_annotsv_copied.setdefault(sid, set()).add(kind)
                _log(f"[copy] {sid}: {annotsv_src} → {annotsv_dst}")

        # 5. Stop-gap chain (no ClinVar; pipeline already populates it).
        _set_step(job_id, "stop-gaps")
        def track_stopgaps(line: str) -> None:
            match = re.search(r"\[(stopgaps-step|sample-step)]\s+([a-z0-9-]+)\s+start", line)
            if match:
                group = "stop-gaps" if match.group(1) == "stopgaps-step" else "sample-step"
                _set_step(job_id, f"{group}:{match.group(2)}")

        stopgap_count = len(samples)
        for stopgap_index, sample in enumerate(samples):
            sid = sample["sample_id"]
            _update(
                job_id,
                stopgap_sample_index=stopgap_index,
                stopgap_sample_count=stopgap_count,
                stopgap_sample_id=sid,
            )
            gui_tsv = TERTIARY_OUTPUT_ROOT / sid / "snv_indel.annotated.tsv"
            stop_args = [str(scripts / "run_stopgaps.sh"),
                         "--tsv",    str(gui_tsv),
                         "--sample", sid]
            if pipeline_annotsv_copied.get(sid) == {"cnv", "sv"}:
                stop_args += ["--skip-cnv"]
                _log(f"[stop-gaps] {sid}: skip AnnotSV fallback; pipeline CNV/SV already copied")
            elif mode == "dragen":
                stop_args += ["--dragen-cnv-source", sample["vcf_path"]]
            elif mode == "inhouse":
                if sample["cnv_vcf"]:
                    stop_args += ["--inhouse-cnv-vcf", sample["cnv_vcf"]]
                if sample["sv_vcf"]:
                    stop_args += ["--inhouse-sv-vcf",  sample["sv_vcf"]]
            if not args.with_extra_vep:
                stop_args.append("--skip-extra-vep")
            _run(stop_args, label=f"stop-gaps {sid}", on_line=track_stopgaps)

        finished_at = _now()
        _set_step(job_id, "done", state="done", finished_at=finished_at)
        _log(f"[tertiary_run] DONE. samples={','.join(sample_ids)}")
        return 0

    except Exception as e:
        traceback.print_exc()
        _update(job_id,
                state="failed",
                error=str(e),
                finished_at=_now())
        _log(f"[tertiary_run] FAILED: {e}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
