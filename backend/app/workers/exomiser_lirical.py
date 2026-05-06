"""RQ worker entry point for re-running Exomiser + LIRICAL on a sample.

Steps (mirrors vcf-analysis-hg19.R §1303–§1510, ported to Python):

  1. Read sample_metadata.json → grab HPO, panels, vcf_path, build.
  2. Render `analysis_files/exomiser_input.yml` from the template, with
     vcf / hpoIds / genomeAssembly / outputDirectory filled in.
  3. java -jar exomiser-cli.jar --analysis exomiser_input.yml
       → produces exomiser_result.variants.tsv
  4. Render `analysis_files/lirical_input.yaml` similarly and run LIRICAL.
       → produces lirical_result.variant.tsv
  5. Parse both outputs and write the spec-compliant
     `exomiser_results.tsv` and `lirical_results.tsv` next to
     snv_indel.annotated.tsv.

All long steps update job_store every chunk so the UI poller sees
progress.
"""
from __future__ import annotations

import json
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import yaml

from ..config import (
    EXOMISER_DATA_HG19,
    EXOMISER_DATA_HG38,
    EXOMISER_HOME,
    EXOMISER_JAR,
    EXOMISER_PROPS,
    JAVA_BIN,
    JAVA_OPTS,
    LIRICAL_JAR,
    REPO_ROOT,
    TERTIARY_OUTPUT_ROOT,
)
from ..services import analyses_store, job_store
from .results_parser import parse_exomiser_variants_tsv, parse_lirical_variant_tsv

EXOMISER_TEMPLATE = REPO_ROOT / "phenotype_reference" / "exomiser_input.yml"
LIRICAL_TEMPLATE  = REPO_ROOT / "phenotype_reference" / "lirical_input.yaml"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _read_yaml(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _write_yaml(path: Path, doc: dict) -> None:
    path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")


def _hpo_ids(hpo_list: list) -> list[str]:
    return [
        h["phenotype"] for h in (hpo_list or [])
        if isinstance(h, dict)
        and isinstance(h.get("phenotype"), str)
        and h["phenotype"].startswith("HP:")
    ]


def _render_exomiser_yml(template: Path, work_dir: Path, vcf: Path,
                         hpo_ids: list[str], assembly: str) -> Path:
    doc = _read_yaml(template)
    analysis = doc.setdefault("analysis", {})
    analysis["genomeAssembly"] = assembly
    analysis["vcf"]    = str(vcf)
    analysis["hpoIds"] = hpo_ids
    output = doc.setdefault("outputOptions", {})
    output["outputDirectory"] = str(work_dir / "analysis_files")
    output["outputFileName"]  = "exomiser_result"
    output["outputFormats"]   = ["TSV_VARIANT"]
    out = work_dir / "analysis_files" / "exomiser_input.yml"
    out.parent.mkdir(parents=True, exist_ok=True)
    _write_yaml(out, doc)
    return out


def _render_lirical_yaml(template: Path, work_dir: Path, sample_id: str,
                         vcf: Path, hpo_ids: list[str]) -> Path:
    doc = _read_yaml(template)
    doc["sampleId"]       = sample_id
    doc["hpoIds"]         = hpo_ids
    doc["negatedHpoIds"]  = doc.get("negatedHpoIds") or []
    doc["vcf"]            = str(vcf)
    out = work_dir / "analysis_files" / "lirical_input.yaml"
    out.parent.mkdir(parents=True, exist_ok=True)
    _write_yaml(out, doc)
    return out


def _run(cmd: list[str], log_path: Path, cwd: Path | None = None,
         timeout: int = 60 * 60) -> int:
    """Run a subprocess; tee stdout/stderr to a log file."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as logf:
        logf.write(f"\n# {_now()} $ {' '.join(shlex.quote(c) for c in cmd)}\n")
        logf.flush()
        proc = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=logf,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
    return proc.returncode


def run_exomiser_lirical(job_id: str, sample_id: str) -> dict:
    """RQ worker entry. Updates job_store as it progresses.

    Reads VCF + genome_build from the patient's sample_metadata.json,
    reads HPO from the active analysis version's analysis.json, and
    writes every intermediate (analysis_files/) + sidecar
    (exomiser_results.tsv, lirical_results.tsv) into the version
    directory. Pre-migration samples (no analyses/ dir) fall back to
    the sample root so old layouts keep working.
    """
    sample_root = TERTIARY_OUTPUT_ROOT / sample_id
    if not sample_root.is_dir():
        raise RuntimeError(f"sample dir not found: {sample_root}")

    meta_path = sample_root / "sample_metadata.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    vcf  = Path(meta.get("vcf_path") or "")
    if not vcf.is_file():
        raise RuntimeError(f"vcf not found: {vcf}")

    assembly = (meta.get("genome_build") or "hg38").lower()
    if assembly not in ("hg19", "hg38"):
        raise RuntimeError(f"unsupported genome_build: {assembly!r}")

    # HPO comes from the active analysis version; fall back to legacy
    # meta.hpo for pre-migration samples that haven't been touched yet.
    active = analyses_store.active_version(sample_id)
    if active:
        analysis = analyses_store.read_version(sample_id, active) or {}
        hpo_list = analysis.get("hpo") or []
    else:
        hpo_list = meta.get("hpo") or []
    hpo_ids = _hpo_ids(hpo_list)
    if not hpo_ids:
        raise RuntimeError("no HP: terms in active analysis (or legacy meta.hpo)")

    work_dir = analyses_store.active_version_dir(sample_id)
    work_dir.mkdir(parents=True, exist_ok=True)
    log_path = work_dir / "analysis_files" / "rerun.log"

    # ---- Exomiser ----
    job_store.update(job_id, {
        "status": "running",
        "step": "exomiser:render",
        "started_at": _now(),
    })
    exo_yml = _render_exomiser_yml(
        EXOMISER_TEMPLATE, work_dir, vcf, hpo_ids, assembly,
    )
    job_store.update(job_id, {"step": "exomiser:run"})
    exo_cmd = [
        JAVA_BIN, *shlex.split(JAVA_OPTS),
        f"-Dspring.config.location={EXOMISER_PROPS}",
        "-jar", str(EXOMISER_JAR),
        "--analysis", str(exo_yml),
    ]
    rc = _run(exo_cmd, log_path, cwd=EXOMISER_HOME, timeout=2 * 60 * 60)
    if rc != 0:
        job_store.update(job_id, {"status": "failed", "step": "exomiser:run", "exit_code": rc})
        raise RuntimeError(f"Exomiser exited {rc}; see {log_path}")

    # ---- LIRICAL ----
    job_store.update(job_id, {"step": "lirical:render"})
    lir_yaml = _render_lirical_yaml(LIRICAL_TEMPLATE, work_dir, sample_id, vcf, hpo_ids)

    job_store.update(job_id, {"step": "lirical:run"})
    ed_dir = EXOMISER_DATA_HG38 if assembly == "hg38" else EXOMISER_DATA_HG19
    ed_flag = "-ed38" if assembly == "hg38" else "-ed19"
    lir_out_dir = work_dir / "analysis_files" / "lirical"
    lir_out_dir.mkdir(parents=True, exist_ok=True)
    lir_cmd = [
        JAVA_BIN, *shlex.split(JAVA_OPTS),
        "-jar", str(LIRICAL_JAR),
        "yaml",
        "-y", str(lir_yaml),
        "--assembly", assembly,
        ed_flag, str(ed_dir),
        "-x", "lirical_result.variant",
        "-o", str(lir_out_dir),
        "-f", "tsv",
    ]
    rc = _run(lir_cmd, log_path, timeout=60 * 60)
    if rc != 0:
        job_store.update(job_id, {"status": "failed", "step": "lirical:run", "exit_code": rc})
        raise RuntimeError(f"LIRICAL exited {rc}; see {log_path}")

    # ---- Parse results into spec-compliant sidecars ----
    job_store.update(job_id, {"step": "parse"})
    exo_var_tsv = work_dir / "analysis_files" / "exomiser_result.variants.tsv"
    lir_var_tsv = lir_out_dir / "lirical_result.variant.tsv"
    n_exo = parse_exomiser_variants_tsv(exo_var_tsv, work_dir / "exomiser_results.tsv")
    n_lir = parse_lirical_variant_tsv(lir_var_tsv, work_dir / "lirical_results.tsv")

    result = {
        "status":      "succeeded",
        "step":        "done",
        "finished_at": _now(),
        "n_exomiser_variants": n_exo,
        "n_lirical_variants":  n_lir,
    }
    job_store.update(job_id, result)
    return result
