"""Pipeline worker — spawned by dragen_jobs.start_job().

Runs the 4-script chain end-to-end on either a DRAGEN hard-filtered
VCF (`--mode dragen`) or an in-house ensemble Nextflow output
(`--mode inhouse`). Steps (same shape for both modes; only the inputs
differ):

    1. annotate_mito_vcf.sh         → tertiary_output/<SID>/mito.annotated.tsv
    2. stage_dragen_for_tertiary.sh → nf_stage/<SID>/04_snv_indel/...
    3. nextflow main_tertiary.nf    → tertiary_output/<SID>/<SID>.snv_indel.annotated.tsv
    4. run_stopgaps.sh              → ClinVar / filter / GeneBe / extra-VEP / CNV-AnnotSV

Mode differences:
  dragen  — step 1 reads the hard-filtered VCF (extracts chrM);
            step 5 (inside run_stopgaps) auto-discovers sibling
            <SID>.cnv.vcf.gz + cnv_sv.vcf.gz and subtracts overlap.
  inhouse — step 1 reads the explicit <SID>.mito.vcf.gz;
            step 5 runs AnnotSV separately on gcnv + delly (no
            subtraction — different callers, no event overlap).

Started by `python3 -m app.workers.dragen_run --job-id … --vcf …`.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

from ..config import (DRAGEN_JOBS_DIR, NGS_UI_HOME, REPO_ROOT,
                       TERTIARY_OUTPUT_ROOT)
from ..services import dragen_jobs


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


def _run(cmd: list[str], *, label: str) -> None:
    """Stream a subprocess's stdout/stderr into this worker's stdout
    (which is already redirected to log.txt by dragen_jobs.start_job).
    Raises on non-zero exit so the outer try/except records failure.
    """
    print(f"\n========================= [{label}] =========================",
          flush=True)
    print("$", " ".join(cmd), flush=True)
    proc = subprocess.run(cmd, check=False)
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
    tsv = sample_dir / f"{sid}.snv_indel.annotated.tsv"

    scripts = REPO_ROOT / "scripts"
    nf_work  = NGS_UI_HOME / "nf_work" / sid
    nf_stage = NGS_UI_HOME / "nf_stage" / sid

    # In-house mode reads chrM from the explicit mito VCF; DRAGEN mode
    # extracts chrM from the hard-filtered VCF (the worker's --vcf is
    # the whole-genome input there).
    mito_in = args.mito_vcf if (mode == "inhouse" and args.mito_vcf) else vcf

    _update(job_id, state="running", step="mito", started_at=_now())
    try:
        # 1. Mito
        if mode == "inhouse" and not args.mito_vcf:
            print("[1/4 mito] skipped — in-house mode but --mito-vcf empty",
                  flush=True)
        else:
            _run([str(scripts / "annotate_mito_vcf.sh"),
                  "--in",     mito_in,
                  "--sample", sid,
                  "--outdir", str(sample_dir)],
                 label="1/4 mito")

        # 2. Stage (gnomAD AF + gene-body BED filter, caller-agnostic)
        _update(job_id, step="stage")
        _run([str(scripts / "stage_dragen_for_tertiary.sh"),
              "--in",     vcf,
              "--sample", sid],
             label="2/4 stage")

        # 3. Nextflow tertiary pipeline
        _update(job_id, step="nextflow")
        _run([
            "nextflow",
            "-c", "/home/pipeline/tertiary_code/nextflow_tertiary.config",
            "run", "/home/pipeline/tertiary_code/main_tertiary.nf",
            "-profile", "dgm",
            "-work-dir", str(nf_work),
            "--sample_id", sid,
            "--input_dir", str(nf_stage),
            "--seq_type",  "WGS",
            "--out_dir",   str(TERTIARY_OUTPUT_ROOT),
        ], label="3/4 nextflow")

        if not tsv.is_file():
            raise RuntimeError(f"nextflow finished but TSV not found: {tsv}")

        # 4. Stop-gap chain. CNV/SV branch depends on mode:
        #    dragen  → sibling discovery + cnv↔cnv_sv subtraction
        #    inhouse → explicit gcnv + delly VCFs, no subtraction
        _update(job_id, step="stop-gaps")
        stop_args = [str(scripts / "run_stopgaps.sh"),
                     "--tsv",    str(tsv),
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
        _run(stop_args, label="4/4 stop-gaps")

        _update(job_id, state="done", step="done", finished_at=_now())
        print("\n[dragen_run] DONE.", flush=True)
        return 0

    except Exception as e:
        traceback.print_exc()
        _update(job_id,
                state="failed",
                error=str(e),
                finished_at=_now())
        print(f"\n[dragen_run] FAILED: {e}", flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
