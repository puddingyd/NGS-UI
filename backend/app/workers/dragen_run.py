"""Pipeline worker — spawned by dragen_jobs.start_job().

Runs the 4-script chain end-to-end on either a DRAGEN hard-filtered
VCF (`--mode dragen`) or an in-house ensemble Nextflow output
(`--mode inhouse`). Steps:

    1. annotate_mito_vcf.sh           → NGS_UI/tertiary_output/<SID>/mito.annotated.tsv
    2. detect existing pipeline output
       ├ HIT  → skip 3+4, reuse the pre-existing .acmg.tsv
       └ MISS → run 3+4 below
    3. stage_dragen_for_tertiary.sh   → nf_stage/<SID>/04_snv_indel/...
    4. nextflow main_tertiary.nf      → /home/pipeline/tertiary_output/<SID>/
                                          03_acmg/<SID>.snv_indel.acmg.tsv
    5. copy pipeline TSV              → NGS_UI/tertiary_output/<SID>/
                                          snv_indel.annotated.tsv
                                          + pipeline_source.json (audit)
    6. run_stopgaps.sh                → filter / GeneBe / extra-VEP / CNV-AnnotSV
                                          + pre-build snv_indel.review.tsv
                                          (ClinVar removed — pipeline already
                                           does it; GeneBe writes a SECOND
                                           opinion to GENEBE_* columns)

Mode differences:
  dragen  — step 1 reads the hard-filtered VCF (extracts chrM);
            step 6's CNV/SV branch auto-discovers sibling
            <SID>.cnv.vcf.gz + <SID>.sv.vcf.gz from the same dir.
  inhouse — step 1 reads the explicit <SID>.mito.vcf.gz;
            step 6 runs AnnotSV separately on gcnv + delly.

Existing-output detection (step 2): tries
    /home/pipeline/tertiary_output/<SID>/03_acmg/<SID>.snv_indel.acmg.tsv
    /home/pipeline/tertiary_output/<stripped SID>/03_acmg/<stripped SID>.snv_indel.acmg.tsv
in that order, where <stripped SID> drops the -dragen / -inhouse
suffix. Production runs by the pipeline team don't carry our
suffix, so the second candidate matches their output and we reuse
it (saving ~30 min per WES sample).

Started by `python3 -m app.workers.dragen_run --job-id … --vcf …`.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

from ..config import (NGS_UI_HOME, PIPELINE_OUT_ROOT, REPO_ROOT,
                       TERTIARY_OUTPUT_ROOT)
from ..services import dragen_jobs


def _strip_sid_suffix(sid: str) -> str:
    """Drop the -dragen / -inhouse caller suffix the GUI adds for
    directory disambiguation. Returns sid unchanged when no suffix
    matches."""
    for suf in ("-dragen", "-inhouse", "-WES", "-WGS"):
        if sid.endswith(suf):
            return sid[: -len(suf)]
    return sid


def _find_pipeline_acmg_tsv(sid: str) -> Path | None:
    """Look for an existing <SID>.snv_indel.acmg.tsv under the pipeline
    production output root. Tries the full (suffixed) SID first; falls
    back to the stripped base SID so production runs by the pipeline
    team (no caller suffix) are reused.
    """
    candidates = [sid]
    base = _strip_sid_suffix(sid)
    if base != sid:
        candidates.append(base)
    for s in candidates:
        p = PIPELINE_OUT_ROOT / s / "03_acmg" / f"{s}.snv_indel.acmg.tsv"
        if p.is_file():
            return p
    return None


def _track_pipeline_source(sample_dir: Path, source: Path) -> None:
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


def _now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


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
    ap.add_argument("--vcf",    required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--mode",   default="dragen", choices=["dragen", "inhouse"])
    ap.add_argument("--with-extra-vep", action="store_true")
    # In-house only — explicit sibling VCF paths from the index.
    ap.add_argument("--cnv-vcf",  default="")
    ap.add_argument("--sv-vcf",   default="")
    ap.add_argument("--mito-vcf", default="")
    args = ap.parse_args()

    _load_secrets()

    job_id = args.job_id
    vcf    = args.vcf
    sid    = args.sample
    mode   = args.mode
    sample_dir = TERTIARY_OUTPUT_ROOT / sid
    sample_dir.mkdir(parents=True, exist_ok=True)
    # NGS-UI side TSV path the adapter reads (no SID prefix).
    gui_tsv = sample_dir / "snv_indel.annotated.tsv"

    scripts = REPO_ROOT / "scripts"
    nf_work  = NGS_UI_HOME / "nf_work" / sid
    nf_stage = NGS_UI_HOME / "nf_stage" / sid

    # In-house mode reads chrM from the explicit mito VCF; DRAGEN mode
    # extracts chrM from the hard-filtered VCF (the worker's --vcf is
    # the whole-genome input there).
    mito_in = args.mito_vcf if (mode == "inhouse" and args.mito_vcf) else vcf
    seq_type = "WGS" if mode == "dragen" else "WES"

    started_at = _now()
    _set_step(job_id, "mito", state="running", started_at=started_at)
    try:
        # 1. Mito
        if mode == "inhouse" and not args.mito_vcf:
            _log("[mito] skipped — in-house mode but --mito-vcf empty")
        else:
            _run([str(scripts / "annotate_mito_vcf.sh"),
                  "--in",     mito_in,
                  "--sample", sid,
                  "--outdir", str(sample_dir)],
                 label="1/4 mito")

        # 2. Reuse existing pipeline output if the production pipeline
        # has already processed this sample (sid or its stripped base).
        _set_step(job_id, "detect-pipeline-output")
        existing = _find_pipeline_acmg_tsv(sid)
        if existing is not None:
            _log(f"[detect] reusing existing pipeline TSV: {existing}")
        else:
            _log("[detect] no existing pipeline TSV — running nextflow")

            # 3. Stage. DRAGEN's 5M-row hard-filter VCF still goes
            # through the full stager (gnomAD AF<0.01 + gene-body BED
            # filter + DV+HC synthesis) so Pangolin doesn't segfault
            # on the splice candidate set. In-house ensemble VCFs are
            # already small (~40k rows) AND already 2-sample, so we
            # just symlink them straight into the nf_stage layout
            # without any filtering — Nextflow sees every variant.
            _set_step(job_id, "stage")
            if mode == "inhouse":
                source_sid = Path(vcf).name.removesuffix(".ensemble.fixed.vcf.gz")
                if source_sid == sid:
                    _symlink_inhouse_into_nf_stage(Path(vcf), sid, nf_stage)
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

            # 4. Nextflow → /home/pipeline/tertiary_output/<SID>/...
            _set_step(job_id, "nextflow")
            nextflow_stages = [
                ("add-callers-tag", "ADD_CALLERS_TAG"),
                ("filter-for-annotation", "FILTER_FOR_ANNOTATION"),
                ("vep-annotate", "VEP_ANNOTATE"),
                ("pangolin-score", "PANGOLIN_SCORE"),
                ("parse-csq", "PARSE_CSQ"),
                ("acmg-classify", "ACMG_CLASSIFY"),
            ]
            nextflow_stage_rank = -1

            def track_nextflow(line: str) -> None:
                nonlocal nextflow_stage_rank
                if not re.search(r"\|\s*[01]\s+of\s+1", line):
                    return
                for rank, (slug, token) in enumerate(nextflow_stages):
                    if token in line and rank > nextflow_stage_rank:
                        nextflow_stage_rank = rank
                        _set_step(job_id, f"nextflow:{slug}")
                        return

            _run([
                "nextflow",
                "-c", "/home/pipeline/tertiary_code/nextflow_tertiary.config",
                "run", "/home/pipeline/tertiary_code/main_tertiary.nf",
                "-profile", "dgm",
                "-work-dir", str(nf_work),
                "--sample_id", sid,
                "--input_dir", str(nf_stage),
                "--seq_type",  seq_type,
                "--out_dir",   str(PIPELINE_OUT_ROOT),
            ], label="2b/4 nextflow", on_line=track_nextflow)

            existing = _find_pipeline_acmg_tsv(sid)
            if existing is None:
                raise RuntimeError(
                    "nextflow finished but expected acmg.tsv not found under "
                    f"{PIPELINE_OUT_ROOT}/{sid}/03_acmg/"
                )

        # 5. Copy pipeline TSV → NGS-UI side; record source audit.
        _set_step(job_id, "copy-pipeline-tsv")
        shutil.copyfile(existing, gui_tsv)
        _track_pipeline_source(sample_dir, existing)
        _log(f"[copy] {existing} → {gui_tsv}")

        # 6. Stop-gap chain (no ClinVar; pipeline already populates it).
        _set_step(job_id, "stop-gaps")
        stop_args = [str(scripts / "run_stopgaps.sh"),
                     "--tsv",    str(gui_tsv),
                     "--sample", sid]
        if mode == "dragen":
            stop_args += ["--dragen-cnv-source", vcf]
        elif mode == "inhouse":
            if args.cnv_vcf:
                stop_args += ["--inhouse-cnv-vcf", args.cnv_vcf]
            if args.sv_vcf:
                stop_args += ["--inhouse-sv-vcf",  args.sv_vcf]
        if not args.with_extra_vep:
            stop_args.append("--skip-extra-vep")
        def track_stopgaps(line: str) -> None:
            match = re.search(r"\[stopgaps-step]\s+([a-z0-9-]+)\s+start", line)
            if match:
                _set_step(job_id, f"stop-gaps:{match.group(1)}")

        _run(stop_args, label="3/4 stop-gaps", on_line=track_stopgaps)

        finished_at = _now()
        _set_step(job_id, "done", state="done", finished_at=finished_at)
        _log("[tertiary_run] DONE.")
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
