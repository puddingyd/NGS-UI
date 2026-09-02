"""Pipeline worker — spawned by dragen_jobs.start_job().

Runs the tertiary chain end-to-end on either a DRAGEN hard-filtered
VCF (`--mode dragen`) or an in-house ensemble Nextflow output
(`--mode inhouse`). Steps:

    1. acquire per-sample locks
    2. write a v3.x samplesheet       → data/jobs/tertiary/<job>/samplesheet.csv
       (legacy env fallback can still stage into nf_stage/<SID>/04_snv_indel)
    3. always run Nextflow -resume    → job-private staging/<source SID>/
    4. validate staged 00-07 outputs
    5. prepare staged post-processing → audit + optional ploidy/MITOMAP output
    6. post-processing                → disposable SNV working TSV → sparse
                                          overlay + review TSV + gene index;
    7. transactionally promote 00-07 and derived files into the live sample,
       preserving reviewer-owned 08_postprocessing state and publishing
       layout.json last
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

Nextflow owns cache/re-run decisions.  Existing published files are never
used to skip the workflow; one shared launch/work lineage per pipeline mode
lets -resume reuse valid task outputs across different batch compositions.
Publishing goes to a fresh job staging root so resumed publishDir rules
cannot leave stale live files.

Started by `python3 -m app.workers.dragen_run --job-id … --vcf …`.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import json
import os
import re
import signal
import shlex
import shutil
import sqlite3
import subprocess
import sys
import time
import traceback
from datetime import datetime, timedelta, timezone
from pathlib import Path

from ..config import (NGS_UI_HOME, PIPELINE_OUT_ROOT, REPO_ROOT,
                       TERTIARY_JOBS_DIR, TERTIARY_NF_WORK_ROOT)
from ..services import dragen_jobs, mitomap_mito, sample_layout

TERTIARY_NEXTFLOW_CONFIG = Path(os.environ.get(
    "NGS_UI_TERTIARY_CONFIG",
    "/home/pipeline/tertiary_code/nextflow_tertiary.config",
))
PIPELINE_CLINVAR_RELEASE = "2026-07-20"

PIPELINE_STAGE_DIRS = (
    "00_prepare",
    "01_vep",
    "02_pangolin",
    "03_acmg",
    "04_mito",
    "05_str",
    "06_cnv_sv",
    "07_pgx",
)
REQUIRED_PIPELINE_STAGE_DIRS = PIPELINE_STAGE_DIRS[:-1]
MANAGED_POSTPROCESSING_NAMES = {
    "layout.json",
    "pipeline_source.json",
    "snv_annotations.sqlite",
    "snv_indel.review.tsv",
    "snv_indel.review.tsv.source.json",
    "snv_gene_index.sqlite",
    "litvar2_annotation.json",
    "clinvar_comparison.json",
    "mito.annotated.tsv",
    "ploidy.vcf.gz",
    "ploidy_qc.txt",
    "cnv.annotated.tsv",
    "sv.annotated.tsv",
}


def _strip_sid_suffix(sid: str) -> str:
    """Drop the -dragen / -nckuh / -inhouse caller suffix the GUI adds for
    directory disambiguation. Returns sid unchanged when no suffix
    matches."""
    for suf in ("-dragen", "-nckuh", "-inhouse"):
        if sid.endswith(suf):
            return sid[: -len(suf)]
    return sid


def _default_ui_sample_id(sample_id: str, mode: str) -> str:
    sid = (sample_id or "").strip()
    if mode == "dragen":
        return sid if sid.lower().endswith("-dragen") else f"{sid}-dragen"
    if mode == "inhouse":
        sid_l = sid.lower()
        if sid_l.endswith("-nckuh") or sid_l.endswith("-inhouse"):
            return sid
        return f"{sid}-nckuh"
    return sid


def _pipeline_candidate_ids(
    sample_id: str,
    source_sample_id: str = "",
    *,
    include_legacy_aliases: bool = True,
) -> list[str]:
    candidates: list[str] = []
    search_ids = [sample_id]
    if include_legacy_aliases:
        search_ids.append(source_sample_id)
    for sid in search_ids:
        sid = (sid or "").strip()
        if sid and sid not in candidates:
            candidates.append(sid)
    if include_legacy_aliases:
        stripped = _strip_sid_suffix(sample_id or "")
        if stripped and stripped not in candidates:
            candidates.append(stripped)
    return candidates


def _find_pipeline_acmg_tsv(
    *sample_ids: str,
    root: Path | None = None,
) -> Path | None:
    """Look for a v3.1 <SID>.snv_indel.acmg.tsv under production output.

    The UI sample ID may intentionally differ from the source VCF sample
    ID (for example, adding "-dragen" to avoid collisions). The pipeline
    itself must run under the source sample ID because v3.1 composes input
    paths from sample_id + input_dir, so try source IDs first and then the
    UI ID / stripped legacy aliases.
    """
    output_root = Path(root) if root is not None else PIPELINE_OUT_ROOT
    candidates: list[str] = []
    for sid in sample_ids:
        if not sid:
            continue
        if sid not in candidates:
            candidates.append(sid)
    for s in candidates:
        p = output_root / s / "03_acmg" / f"{s}.snv_indel.acmg.tsv"
        if p.is_file():
            return p
        d = output_root / s / "03_acmg"
        if d.is_dir():
            hit = next(d.glob("*.snv_indel.acmg.tsv"), None)
            if hit is not None and hit.is_file():
                return hit
    return None


def _find_pipeline_mito_tsv(
    *sample_ids: str,
    root: Path | None = None,
) -> Path | None:
    """Look for a v3.2 <SID>.mito.tsv under production output."""
    output_root = Path(root) if root is not None else PIPELINE_OUT_ROOT
    candidates: list[str] = []
    for sid in sample_ids:
        if not sid:
            continue
        if sid not in candidates:
            candidates.append(sid)
    for s in candidates:
        p = output_root / s / "04_mito" / f"{s}.mito.tsv"
        if p.is_file():
            return p
        d = output_root / s / "04_mito"
        if d.is_dir():
            hit = next(d.glob("*.mito.tsv"), None)
            if hit is not None and hit.is_file():
                return hit
    return None


def _find_pipeline_str_tsv(
    *sample_ids: str,
    root: Path | None = None,
) -> Path | None:
    """Look for v3.x 05_str/<SID>.str.tsv output."""
    output_root = Path(root) if root is not None else PIPELINE_OUT_ROOT
    candidates: list[str] = []
    for sid in sample_ids:
        if not sid:
            continue
        if sid not in candidates:
            candidates.append(sid)
    for s in candidates:
        p = output_root / s / "05_str" / f"{s}.str.tsv"
        if p.is_file():
            return p
        d = output_root / s / "05_str"
        if d.is_dir():
            hit = next(d.glob("*.str.tsv"), None)
            if hit is not None and hit.is_file():
                return hit
    return None


def _find_pipeline_annotsv_tsv(
    kind: str,
    *sample_ids: str,
    root: Path | None = None,
) -> Path | None:
    """Look for v3.2 06_cnv_sv/<SID>.<kind>.annotated.tsv output."""
    output_root = Path(root) if root is not None else PIPELINE_OUT_ROOT
    candidates: list[str] = []
    for sid in sample_ids:
        if not sid:
            continue
        if sid not in candidates:
            candidates.append(sid)
    for s in candidates:
        p = output_root / s / "06_cnv_sv" / f"{s}.{kind}.annotated.tsv"
        if p.is_file():
            return p
        d = output_root / s / "06_cnv_sv"
        if d.is_dir():
            hit = next(d.glob(f"*.{kind}.annotated.tsv"), None)
            if hit is not None and hit.is_file():
                return hit
    return None


def _find_pipeline_pgx_output(
    *sample_ids: str,
    root: Path | None = None,
) -> Path | None:
    """Return any PGx/PharmCAT-looking output file for a sample.

    The formal PGx output contract is still moving, so this intentionally
    accepts the known plan (`pgx.tsv`) plus common PharmCAT / PGx names.
    """
    output_root = Path(root) if root is not None else PIPELINE_OUT_ROOT
    names = [
        "pgx.tsv",
        "pharmcat.json",
        "pharmcat.tsv",
        "pharmcat.html",
    ]
    globs = [
        "*pgx*.tsv", "*PGx*.tsv", "*pharmcat*.json", "*pharmcat*.tsv",
        "*pharmcat*.html", "*PharmCAT*",
    ]
    subdirs = ["", "07_pgx", "pgx", "PGx", "pharmcat", "PharmCAT"]
    for sid in sample_ids:
        if not sid:
            continue
        sample_root = output_root / sid
        if not sample_root.is_dir():
            continue
        for sub in subdirs:
            d = sample_root / sub if sub else sample_root
            if not d.is_dir():
                continue
            for name in names:
                p = d / name
                if p.is_file():
                    return p
            for pat in globs:
                hit = next((p for p in d.glob(pat) if p.is_file()), None)
                if hit is not None:
                    return hit
    return None


def _find_pipeline_pgx_files(
    *sample_ids: str,
    root: Path | None = None,
) -> dict[str, Path]:
    """Return known 07_pgx files keyed by their UI-side destination name."""
    output_root = Path(root) if root is not None else PIPELINE_OUT_ROOT
    wanted = {
        "pgx.tsv": "*.pgx.tsv",
        "pharmcat.report.json": "*.pharmcat.report.json",
        "outside_calls.tsv": "*.outside_calls.tsv",
        "stellarpgx.tsv": "*.stellarpgx.tsv",
        "optitype.tsv": "*.optitype.tsv",
    }
    found: dict[str, Path] = {}
    for sid in sample_ids:
        if not sid:
            continue
        d = output_root / sid / "07_pgx"
        if not d.is_dir():
            continue
        for dst_name, pattern in wanted.items():
            if dst_name in found:
                continue
            exact = d / f"{sid}.{dst_name}"
            if exact.is_file():
                found[dst_name] = exact
                continue
            hit = next(d.glob(pattern), None)
            if hit is not None and hit.is_file():
                found[dst_name] = hit
    return found


def _find_dragen_ploidy_vcf(source_vcf: str | Path, source_sample_id: str) -> Path | None:
    """Find the DRAGEN ploidy VCF corresponding to one source VCF.

    Candidate names are exact so a directory containing multiple samples
    cannot accidentally donate another sample's karyotype.
    """
    source = Path(source_vcf)
    name = source.name
    stems: list[str] = []
    for suffix in (".hard-filtered.vcf.gz", ".vcf.gz"):
        if name.endswith(suffix):
            stems.append(name[:-len(suffix)])
            break
    if source_sample_id:
        stems.append(str(source_sample_id).strip())
    for stem in dict.fromkeys(value for value in stems if value):
        candidate = source.with_name(f"{stem}.ploidy.vcf.gz")
        if candidate.is_file():
            return candidate
    return None


def _find_nckuh_ploidy_vcf(source_vcf: str | Path, source_sample_id: str) -> Path | None:
    """Find 03_alignment_qc/<source>.ploidy.vcf.gz for an NCKUH ensemble VCF."""
    source = Path(source_vcf)
    sample_root = source.parent.parent if source.parent.name == "04_snv_indel" else source.parent
    qc_dir = sample_root / "03_alignment_qc"
    name = source.name
    stems: list[str] = []
    for suffix in (".ensemble.fixed.vcf.gz", ".vcf.gz"):
        if name.endswith(suffix):
            stems.append(name[:-len(suffix)])
            break
    if source_sample_id:
        stems.append(str(source_sample_id).strip())
    for stem in dict.fromkeys(value for value in stems if value):
        candidate = qc_dir / f"{stem}.ploidy.vcf.gz"
        if candidate.is_file():
            return candidate
    return None


def _find_ploidy_qc_txt(
    ploidy_vcf: Path | None,
    *,
    sample_id: str,
    source_sample_id: str,
) -> Path | None:
    """Locate the future-use ploidy QC text without consuming it in the UI."""
    ids = list(dict.fromkeys(
        value for value in (source_sample_id.strip(), sample_id.strip()) if value
    ))
    if ploidy_vcf is not None:
        directory = ploidy_vcf.parent
        stem = ploidy_vcf.name.removesuffix(".ploidy.vcf.gz")
        for name in (
            f"{stem}.ploidy_qc.txt",
            *(f"{sid}.ploidy_qc.txt" for sid in ids),
            "ploidy_qc.txt",
        ):
            candidate = directory / name
            if candidate.is_file():
                return candidate
    for sid_dir in ids:
        directory = PIPELINE_OUT_ROOT / sid_dir / "00_prepare"
        for source_id in ids:
            candidate = directory / f"{source_id}.ploidy_qc.txt"
            if candidate.is_file():
                return candidate
    return None


def _copy_dragen_ploidy_vcf(
    source_vcf: str | Path,
    source_sample_id: str,
    sample_dir: Path,
    sample_id: str | None = None,
) -> tuple[Path, Path] | None:
    source = _find_dragen_ploidy_vcf(source_vcf, source_sample_id)
    sid = sample_id or (
        sample_dir.parent.name
        if sample_dir.name == sample_layout.POSTPROCESSING_DIRNAME
        else sample_dir.name
    )
    destination = sample_layout.scoped_file(
        sample_dir,
        sid,
        "ploidy.vcf.gz",
        for_write=True,
        force_prefixed=True,
    )
    if source is None:
        # Do not retain a stale karyotype when the same LIS ID is re-run
        # against a DRAGEN source that has no matching ploidy sidecar.
        destination.unlink(missing_ok=True)
        return None
    shutil.copyfile(source, destination)
    return source, destination


def _copy_ploidy_artifacts(
    *,
    mode: str,
    source_vcf: str | Path,
    source_sample_id: str,
    sample_id: str,
    post_dir: Path,
) -> tuple[tuple[Path, Path] | None, tuple[Path, Path] | None]:
    source = (
        _find_dragen_ploidy_vcf(source_vcf, source_sample_id)
        if mode == "dragen"
        else _find_nckuh_ploidy_vcf(source_vcf, source_sample_id)
    )
    vcf_destination = sample_layout.scoped_file(
        post_dir,
        sample_id,
        "ploidy.vcf.gz",
        for_write=True,
        force_prefixed=True,
    )
    copied_vcf = None
    if source is not None:
        shutil.copyfile(source, vcf_destination)
        copied_vcf = (source, vcf_destination)
    else:
        vcf_destination.unlink(missing_ok=True)

    qc_source = _find_ploidy_qc_txt(
        source,
        sample_id=sample_id,
        source_sample_id=source_sample_id,
    )
    qc_destination = sample_layout.scoped_file(
        post_dir,
        sample_id,
        "ploidy_qc.txt",
        for_write=True,
        force_prefixed=True,
    )
    copied_qc = None
    if qc_source is not None:
        shutil.copyfile(qc_source, qc_destination)
        copied_qc = (qc_source, qc_destination)
    else:
        qc_destination.unlink(missing_ok=True)
    return copied_vcf, copied_qc


def _pipeline_outputs_for(
    sample_id: str,
    source_sample_id: str,
    *,
    require_pgx: bool,
    require_str: bool = True,
    include_legacy_aliases: bool = True,
    root: Path | None = None,
) -> tuple[dict[str, Path], list[str], str]:
    output_root = Path(root) if root is not None else PIPELINE_OUT_ROOT
    outputs: dict[str, Path] = {}
    found_under = ""
    candidate_ids = _pipeline_candidate_ids(
        sample_id,
        source_sample_id,
        include_legacy_aliases=include_legacy_aliases,
    )
    for sid in candidate_ids:
        checks: list[tuple[str, Path | None]] = [
            ("snv_indel.acmg.tsv", _find_pipeline_acmg_tsv(sid, root=output_root)),
            ("mito.tsv", _find_pipeline_mito_tsv(sid, root=output_root)),
            (
                "cnv.annotated.tsv",
                _find_pipeline_annotsv_tsv("cnv", sid, root=output_root),
            ),
            (
                "sv.annotated.tsv",
                _find_pipeline_annotsv_tsv("sv", sid, root=output_root),
            ),
        ]
        if require_str:
            checks.append(("str.tsv", _find_pipeline_str_tsv(sid, root=output_root)))
        if require_pgx:
            checks.append((
                "PGx/PharmCAT",
                _find_pipeline_pgx_output(sid, root=output_root),
            ))
        present = {name: path for name, path in checks if path is not None}
        if not found_under and present:
            found_under = sid
        missing = [name for name, path in checks if path is None]
        if not missing:
            return present, [], sid
    all_names = [
        "snv_indel.acmg.tsv",
        "mito.tsv",
        "cnv.annotated.tsv",
        "sv.annotated.tsv",
    ]
    if require_str:
        all_names.append("str.tsv")
    if require_pgx:
        all_names.append("PGx/PharmCAT")
    missing = []
    for name in all_names:
        if name == "snv_indel.acmg.tsv":
            path = _find_pipeline_acmg_tsv(*candidate_ids, root=output_root)
        elif name == "mito.tsv":
            path = _find_pipeline_mito_tsv(*candidate_ids, root=output_root)
        elif name == "cnv.annotated.tsv":
            path = _find_pipeline_annotsv_tsv(
                "cnv", *candidate_ids, root=output_root
            )
        elif name == "sv.annotated.tsv":
            path = _find_pipeline_annotsv_tsv(
                "sv", *candidate_ids, root=output_root
            )
        elif name == "str.tsv":
            path = _find_pipeline_str_tsv(*candidate_ids, root=output_root)
        else:
            path = _find_pipeline_pgx_output(*candidate_ids, root=output_root)
        if path is None:
            missing.append(name)
        else:
            outputs[name] = path
    return outputs, missing, found_under


def _nextflow_context(samples: list[dict]) -> tuple[Path, Path]:
    """Return the shared cache lineage for one tertiary pipeline mode.

    A Nextflow task hash includes the session ID.  Keeping one launch/session
    lineage per mode therefore lets the same sample/input task resume even
    when the surrounding batch changes.  A mode-wide advisory lock serializes
    access to the session's LevelDB cache.  A migrated lineage keeps using the
    adopted session's original work root because cached task metadata and work
    outputs are both required for resume.
    """
    mode = str(samples[0].get("mode") or "") if samples else ""
    scope = "dragen" if mode == "dragen" else "nckuh"
    context = TERTIARY_NF_WORK_ROOT / "contexts" / f"shared-{scope}"
    launch_dir = context / "launch"
    work_dir = context / "work"
    pointer = context / "work-dir.txt"
    try:
        pointed_work = _validated_nf_work_path(
            pointer.read_text(encoding="utf-8").strip(),
            require_exists=True,
        )
    except OSError:
        pointed_work = None
    if pointed_work is not None:
        work_dir = pointed_work
    return launch_dir, work_dir


def _command_option(command: list[str], name: str) -> str:
    for index, token in enumerate(command):
        if token == name and index + 1 < len(command):
            return command[index + 1]
        if token.startswith(f"{name}="):
            return token.split("=", 1)[1]
    return ""


def _history_command_samples(command: list[str]) -> set[str]:
    sample_ids: set[str] = set()
    direct = _command_option(command, "--sample_id")
    if direct:
        sample_ids.add(direct)
    samplesheet_value = _command_option(command, "--samplesheet")
    if not samplesheet_value:
        return sample_ids
    samplesheet = Path(samplesheet_value)
    try:
        with samplesheet.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                sample_id = str(row.get("sample_id") or "").strip()
                if sample_id:
                    sample_ids.add(sample_id)
    except OSError:
        pass
    return sample_ids


def _validated_nf_work_path(
    value: str | Path,
    *,
    relative_to: Path | None = None,
    require_exists: bool,
) -> Path | None:
    raw = str(value).strip()
    if not raw:
        return None
    path = Path(raw)
    if not path.is_absolute():
        path = (relative_to or TERTIARY_NF_WORK_ROOT) / path
    try:
        resolved = path.resolve()
        work_root = TERTIARY_NF_WORK_ROOT.resolve()
        if not resolved.is_relative_to(work_root):
            return None
        if require_exists and not resolved.is_dir():
            return None
    except OSError:
        return None
    return resolved


def _write_shared_work_pointer(launch_dir: Path, work_dir: Path) -> None:
    context = Path(launch_dir).parent
    validated = _validated_nf_work_path(work_dir, require_exists=True)
    if validated is None:
        raise RuntimeError(f"unsafe or missing shared Nextflow work root: {work_dir}")
    context.mkdir(parents=True, exist_ok=True)
    pointer = context / "work-dir.txt"
    temporary = context / f".work-dir.{os.getpid()}.tmp"
    temporary.write_text(str(validated) + "\n", encoding="utf-8")
    temporary.replace(pointer)


def _seed_shared_nextflow_context(
    samples: list[dict],
    *,
    launch_dir: Path,
) -> tuple[str, Path] | None:
    """Adopt one previous session when creating a mode-wide cache lineage.

    Cache entries from different session IDs cannot be merged because the
    session ID participates in every task hash.  Choose the prior session with
    the greatest overlap with the requested samples, preferring a successful
    and then newer run, and preserve its original work directories.
    """
    target_nf = Path(launch_dir) / ".nextflow"
    if (target_nf / "history").is_file():
        return None
    mode = str(samples[0].get("mode") or "") if samples else ""
    expected_pipeline_type = "dragen" if mode == "dragen" else "nckuh"
    requested_ids = {
        str(value)
        for sample in samples
        for value in (sample.get("sample_id"), sample.get("source_sample_id"))
        if value
    }
    candidate_nf_roots: list[Path] = [REPO_ROOT / ".nextflow"]
    contexts = TERTIARY_NF_WORK_ROOT / "contexts"
    if contexts.is_dir():
        for candidate_launch in sorted(contexts.glob("*/launch")):
            candidate_nf = candidate_launch / ".nextflow"
            if candidate_nf != target_nf:
                candidate_nf_roots.append(candidate_nf)

    candidates: list[tuple[int, int, str, str, str, Path, Path]] = []
    for source_nf in candidate_nf_roots:
        history = source_nf / "history"
        try:
            lines = history.read_text(encoding="utf-8").splitlines()
        except OSError:
            continue
        for line in lines:
            fields = line.split("\t", 6)
            if len(fields) != 7 or fields[3] not in {"OK", "ERR"}:
                continue
            session_id = fields[5].strip()
            source_cache = source_nf / "cache" / session_id
            if not session_id or not source_cache.is_dir():
                continue
            try:
                command = shlex.split(fields[6])
            except ValueError:
                continue
            source_work = _validated_nf_work_path(
                _command_option(command, "-work-dir"),
                relative_to=source_nf.parent,
                require_exists=True,
            )
            if source_work is None:
                continue
            pipeline_type = _command_option(command, "--pipeline_type").lower()
            prior_samples = _history_command_samples(command)
            overlap = len(requested_ids & prior_samples)
            if pipeline_type:
                accepted_pipeline_types = {expected_pipeline_type}
                if expected_pipeline_type == "nckuh":
                    accepted_pipeline_types.add("inhouse")
                if pipeline_type not in accepted_pipeline_types:
                    continue
            elif overlap == 0:
                # Legacy --sample_id runs did not record pipeline_type. Only
                # adopt one when its sample identity matches this request.
                continue
            candidates.append((
                overlap,
                1 if fields[3] == "OK" else 0,
                fields[0],
                line,
                session_id,
                source_nf,
                source_work,
            ))
    if not candidates:
        return None

    (
        _overlap,
        _success,
        _timestamp,
        line,
        session_id,
        source_nf,
        source_work,
    ) = max(
        candidates,
        key=lambda item: (item[0], item[1], item[2]),
    )
    source_cache = source_nf / "cache" / session_id
    target_cache = target_nf / "cache" / session_id
    try:
        target_cache.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source_cache, target_cache)
        (target_nf / "history").write_text(line + "\n", encoding="utf-8")
    except OSError as exc:
        shutil.rmtree(target_nf, ignore_errors=True)
        _log(f"[nextflow] shared cache seed skipped: {exc}")
        return None
    return session_id, source_work


def _acquire_nextflow_cache_lock(
    mode: str,
    *,
    blocking: bool = True,
) -> object:
    """Serialize access to one mode-wide Nextflow session cache."""
    lock_dir = TERTIARY_JOBS_DIR / ".nextflow_locks"
    lock_dir.mkdir(parents=True, exist_ok=True)
    scope = "dragen" if mode == "dragen" else "nckuh"
    handle = (lock_dir / f"{scope}.lock").open("a+", encoding="utf-8")
    operation = fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB)
    try:
        fcntl.flock(handle.fileno(), operation)
    except BaseException:
        handle.close()
        raise
    handle.seek(0)
    handle.truncate()
    handle.write(f"pid={os.getpid()} acquired_at={_now()}\n")
    handle.flush()
    return handle


def _release_nextflow_cache_lock(handle: object | None) -> None:
    if handle is None:
        return
    try:
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    except OSError:
        pass
    try:
        handle.close()
    except OSError:
        pass


def _acquire_sample_locks(sample_ids: list[str]) -> list[object]:
    """Hold non-blocking advisory locks for every UI sample in a job."""
    lock_dir = TERTIARY_JOBS_DIR / ".sample_locks"
    lock_dir.mkdir(parents=True, exist_ok=True)
    handles: list[object] = []
    try:
        for sample_id in sorted(set(sample_ids)):
            handle = (lock_dir / f"{sample_id}.lock").open("a+", encoding="utf-8")
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                handle.close()
                raise RuntimeError(
                    f"tertiary analysis is already running for {sample_id}"
                ) from exc
            handle.seek(0)
            handle.truncate()
            handle.write(f"pid={os.getpid()} acquired_at={_now()}\n")
            handle.flush()
            handles.append(handle)
        return handles
    except BaseException:
        _release_sample_locks(handles)
        raise


def _release_sample_locks(handles: list[object]) -> None:
    for handle in reversed(handles):
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        except OSError:
            pass
        try:
            handle.close()
        except OSError:
            pass


def _prepare_job_staging(job_id: str) -> tuple[Path, Path]:
    staging_parent = PIPELINE_OUT_ROOT / ".staging"
    rollback_parent = PIPELINE_OUT_ROOT / ".rollback"
    staging_root = staging_parent / job_id
    rollback_root = rollback_parent / job_id
    if staging_root.exists() or rollback_root.exists():
        raise RuntimeError(f"job staging already exists: {job_id}")
    staging_root.mkdir(parents=True, exist_ok=False)
    try:
        rollback_root.mkdir(parents=True, exist_ok=False)
    except BaseException:
        shutil.rmtree(staging_root, ignore_errors=True)
        raise
    return staging_root, rollback_root


def _validate_staged_sample(
    staging_root: Path,
    source_sample_id: str,
    *,
    require_pgx: bool,
) -> dict[str, Path]:
    sample_dir = staging_root / source_sample_id
    required_dirs = list(REQUIRED_PIPELINE_STAGE_DIRS)
    if require_pgx:
        required_dirs.append("07_pgx")
    missing_dirs = [name for name in required_dirs if not (sample_dir / name).is_dir()]
    if missing_dirs:
        raise RuntimeError(
            f"nextflow staging output for {source_sample_id} is missing director"
            f"{'y' if len(missing_dirs) == 1 else 'ies'}: {', '.join(missing_dirs)}"
        )
    outputs, missing, found_under = _pipeline_outputs_for(
        source_sample_id,
        source_sample_id,
        require_pgx=require_pgx,
        require_str=True,
        include_legacy_aliases=False,
        root=staging_root,
    )
    if missing or found_under != source_sample_id:
        detail = ", ".join(missing) if missing else "unexpected sample directory"
        raise RuntimeError(
            f"nextflow staging output validation failed for {source_sample_id}: {detail}"
        )
    return outputs


def _sqlite_set_meta(path: Path, key: str, value: str) -> None:
    if not path.is_file():
        return
    with sqlite3.connect(path) as conn:
        conn.execute(
            "INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)",
            (key, value),
        )
        conn.commit()


def _rebase_staged_derived_paths(
    *,
    sample_id: str,
    stage_post_dir: Path,
    final_raw_tsv: Path,
    final_post_dir: Path,
) -> None:
    """Make staged SQLite/JSON signatures valid after live promotion."""
    final_raw = str(final_raw_tsv.absolute())
    overlay = stage_post_dir / f"{sample_id}.snv_annotations.sqlite"
    gene_index = stage_post_dir / f"{sample_id}.snv_gene_index.sqlite"
    _sqlite_set_meta(overlay, "source_path", final_raw)
    _sqlite_set_meta(gene_index, "source_path", final_raw)

    manifest = stage_post_dir / f"{sample_id}.snv_indel.review.tsv.source.json"
    if manifest.is_file():
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        overlay_sig = payload.get("overlay")
        if isinstance(overlay_sig, dict) and overlay_sig.get("exists"):
            overlay_stat = overlay.stat()
            overlay_sig["path"] = str(
                (final_post_dir / f"{sample_id}.snv_annotations.sqlite").absolute()
            )
            overlay_sig["mtime_ns"] = overlay_stat.st_mtime_ns
            overlay_sig["size"] = overlay_stat.st_size
            overlay_sig["current"] = True
        tmp = manifest.with_suffix(manifest.suffix + ".tmp")
        tmp.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        os.replace(tmp, manifest)


def _path_exists(path: Path) -> bool:
    return path.exists() or path.is_symlink()


def _promote_staged_generation(
    *,
    sample_id: str,
    staged_sample_dir: Path,
    live_sample_dir: Path,
    rollback_sample_dir: Path,
) -> list[tuple[Path | None, Path, Path]]:
    """Promote validated 00-07 + derived state with rollback on failure.

    Reviewer-owned files stay in the existing live 08_postprocessing
    directory.  Only worker-managed prefixed artifacts are replaced, and the
    new layout marker is the final rename in the transaction.
    """
    staged_post = staged_sample_dir / sample_layout.POSTPROCESSING_DIRNAME
    live_post = live_sample_dir / sample_layout.POSTPROCESSING_DIRNAME
    backup_post = rollback_sample_dir / sample_layout.POSTPROCESSING_DIRNAME
    live_sample_dir.mkdir(parents=True, exist_ok=True)
    live_post.mkdir(parents=True, exist_ok=True)
    rollback_sample_dir.mkdir(parents=True, exist_ok=True)

    promoted_legacy = sample_layout.promote_state_tree_in_directory_to_v3(
        live_post,
        sample_id,
        exclude_names=MANAGED_POSTPROCESSING_NAMES,
    )
    if promoted_legacy:
        _log(
            f"[layout] {sample_id}: copied {len(promoted_legacy)} legacy "
            "reviewer/analysis state file(s) to v3 names"
        )

    operations: list[tuple[Path | None, Path, Path]] = []

    def replace(new_path: Path | None, live_path: Path, backup_path: Path) -> None:
        backup_path.parent.mkdir(parents=True, exist_ok=True)
        if _path_exists(live_path):
            os.replace(live_path, backup_path)
        operations.append((new_path, live_path, backup_path))
        if new_path is not None:
            if not _path_exists(new_path):
                raise RuntimeError(f"staged promotion source disappeared: {new_path}")
            live_path.parent.mkdir(parents=True, exist_ok=True)
            os.replace(new_path, live_path)

    marker_name = f"{sample_id}.layout.json"
    try:
        for name in PIPELINE_STAGE_DIRS:
            staged = staged_sample_dir / name
            replace(
                staged if staged.is_dir() else None,
                live_sample_dir / name,
                rollback_sample_dir / name,
            )

        for name in sorted(MANAGED_POSTPROCESSING_NAMES - {"layout.json"}):
            staged = staged_post / f"{sample_id}.{name}"
            replace(
                staged if staged.is_file() else None,
                live_post / f"{sample_id}.{name}",
                backup_post / f"{sample_id}.{name}",
            )

        # Activation marker is deliberately last.
        replace(
            staged_post / marker_name,
            live_post / marker_name,
            backup_post / marker_name,
        )
    except BaseException:
        _rollback_promotion_operations(operations)
        raise
    return operations


def _rollback_promotion_operations(
    operations: list[tuple[Path | None, Path, Path]],
) -> None:
    for new_path, live_path, backup_path in reversed(operations):
        try:
            if _path_exists(live_path):
                if new_path is not None:
                    new_path.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(live_path, new_path)
                elif live_path.is_dir() and not live_path.is_symlink():
                    shutil.rmtree(live_path)
                else:
                    live_path.unlink(missing_ok=True)
            if _path_exists(backup_path):
                live_path.parent.mkdir(parents=True, exist_ok=True)
                os.replace(backup_path, live_path)
        except OSError:
            pass


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


def _validate_acmg_tsv(
    path: Path,
    *,
    strict_v31: bool,
    expect_academic_dbnsfp: bool | None = None,
) -> None:
    """Catch stale/partial tertiary outputs before copying them into UI."""
    try:
        header = path.open(encoding="utf-8").readline().rstrip("\n").split("\t")
    except OSError as e:
        raise RuntimeError(f"cannot read pipeline TSV: {path}") from e
    required = {
        "CHROM", "POS", "REF", "ALT",
        "GENE", "TRANSCRIPT", "TRANSCRIPT_TYPE", "HGVS_C", "HGVS_P",
        "CONSEQUENCE", "IMPACT",
        "HGNC_ID", "ACMG_CRITERIA", "ACMG_SCORE", "ACMG_CLASS", "ACMG_NOTES",
    }
    missing = sorted(required - set(header))
    if missing:
        raise RuntimeError(
            f"pipeline TSV missing required columns: {', '.join(missing)}"
        )
    is_v35_transcript_schema = "MANE_ALL" not in header
    expected_cols = 81
    schema_label = "v3.6 transcript schema" if is_v35_transcript_schema else "legacy MANE_ALL schema"
    if strict_v31 and is_v35_transcript_schema and "STRAND_BIAS" not in header:
        raise RuntimeError("pipeline TSV v3.6 schema is missing STRAND_BIAS")
    v36_required = {
        "CLINGEN_VCEP_CLASS", "CLINGEN_VCEP_CRITERIA", "CLINGEN_VCEP_PANEL",
        "REVEL", "MUTPRED2", "MUTPRED2_PRED", "VEST4", "CADD_PHRED",
        "DBNSFP_VERSION", "CLINGEN_AGREEMENT", "PVS1_STRENGTH", "PVS1_REASON",
    }
    if strict_v31:
        missing_v36 = sorted(v36_required - set(header))
        if missing_v36:
            raise RuntimeError(
                "pipeline TSV is not v3.6; missing columns: "
                + ", ".join(missing_v36)
            )
    if strict_v31 and len(header) < expected_cols:
        raise RuntimeError(
            f"pipeline TSV has {len(header)} columns; {schema_label} expects at least {expected_cols}"
        )
    if len(header) != expected_cols:
        _log(
            f"[validate] warning: pipeline TSV has {len(header)} columns "
            f"({schema_label}: {expected_cols})"
        )
    if strict_v31 and expect_academic_dbnsfp is not None:
        try:
            with path.open(encoding="utf-8") as handle:
                handle.readline()
                first = handle.readline().rstrip("\n").split("\t")
            if not first or first == [""]:
                raise IndexError
        except (OSError, IndexError) as exc:
            raise RuntimeError(f"pipeline TSV has no data rows: {path}") from exc
        version_index = header.index("DBNSFP_VERSION")
        actual = first[version_index].strip() if version_index < len(first) else ""
        expected = "5.3a" if expect_academic_dbnsfp else "4.9c"
        if actual != expected:
            raise RuntimeError(
                f"pipeline DBNSFP_VERSION={actual or '(blank)'}; expected {expected} "
                f"for research_only={expect_academic_dbnsfp}"
            )


def _ensure_v36_annotation_versions(path: Path) -> Path:
    """Record the fixed ClinVar release for v3.6 Nextflow output."""
    suffix = ".snv_indel.acmg.tsv"
    source_name = path.name[:-len(suffix)] if path.name.endswith(suffix) else path.stem
    sidecar = path.with_name(f"{source_name}.annotation_versions.json")
    try:
        existing = json.loads(sidecar.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        existing = {}
    payload = existing if isinstance(existing, dict) else {}
    payload["schema_version"] = payload.get("schema_version") or 1
    payload["pipeline_schema"] = "v3.6"
    databases = payload.setdefault("databases", {})
    if not isinstance(databases, dict):
        databases = payload["databases"] = {}
    clinvar = databases.setdefault("clinvar", {})
    if not isinstance(clinvar, dict):
        clinvar = databases["clinvar"] = {}
    clinvar["release_date"] = PIPELINE_CLINVAR_RELEASE
    clinvar["source"] = "Nextflow v3.6 fixed release"
    tmp = sidecar.with_suffix(sidecar.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    os.replace(tmp, sidecar)
    return sidecar


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
        sample_id = _default_ui_sample_id((sample.get("sample_id") or "").strip(), mode)
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
    sample_id: str,
    sample_dir: Path,
    source: Path,
    *,
    source_sample_id: str,
    source_vcf_path: str,
    pipeline_type: str,
    published_source: Path | None = None,
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
        "source_path":  str(published_source or source),
        "source_mtime": mtime,
        "source_sample_id": source_sample_id,
        "source_vcf_path": source_vcf_path,
        "pipeline_type": pipeline_type,
        "annotated_at": _now(),
    }
    sample_layout.scoped_file(
        sample_dir,
        sample_id,
        "pipeline_source.json",
        for_write=True,
        force_prefixed=True,
    ).write_text(
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
    return datetime.now(TAIPEI_TZ).strftime("%Y-%m-%d %H:%M:%S")


def _load_secrets() -> None:
    """Populate os.environ from $NGS_UI_HOME/secrets.env if present.

    uvicorn runs under systemd and doesn't inherit interactive shell
    `export`s. Optional subprocess integrations such as the review-filtered
    GeneBe live fallback therefore load credentials from this file. It uses
    plain KEY=VAL lines (no quoting, no expansion), is git-ignored and should
    be mode 0600. Values already in os.environ win, so systemd Environment=
    can still override them.
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
    if "nextflow_progress_pct" in kw:
        kw["nextflow_progress_pct"] = max(
            float(st.get("nextflow_progress_pct") or 0),
            float(kw.get("nextflow_progress_pct") or 0),
        )
    st.update(kw)
    dragen_jobs.save_state(job_id, st)


def _log(message: str = "") -> None:
    """Write one worker-owned log line with a Taipei timestamp suffix."""
    suffix = f"[{_now()}]"
    print(f"{message} {suffix}" if message else suffix, flush=True)


def _set_step(job_id: str, step: str, **kw) -> None:
    """Persist and log step transitions for post-run timing analysis."""
    now = _now()
    st = dragen_jobs.load_state(job_id) or {}
    if "nextflow_progress_pct" in kw:
        kw["nextflow_progress_pct"] = max(
            float(st.get("nextflow_progress_pct") or 0),
            float(kw.get("nextflow_progress_pct") or 0),
        )
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
    done: int | None = None,
    total: int | None = None,
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
    if done is not None:
        item["done"] = done
    if total is not None:
        item["total"] = total
    if elapsed is not None:
        item["elapsed_seconds"] = round(elapsed, 1)
    history.append(item)
    st["nextflow_step_history"] = history
    dragen_jobs.save_state(job_id, st)


def _elapsed_minutes_label(elapsed: float | None) -> str:
    if elapsed is None:
        return ""
    return f"elapsed={elapsed / 60:.1f}m"


_ANSI_ESCAPE_RE = re.compile(r"\x1b(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")


def _parse_nextflow_progress_line(
    line: str,
    tokens: tuple[str, ...],
    *,
    expected_total: int,
) -> dict[str, object] | None:
    """Parse one Nextflow process status line for batch progress.

    Nextflow may initially print ``1 of 1`` and then grow the denominator as
    more channel items arrive.  Treating that first line as 100% makes a
    long-running per-sample process such as VEP consume its whole weight at
    once and leaves the UI progress bar apparently frozen.  Use the known
    batch size as the denominator while the process is active, and only
    consider the process terminal when Nextflow marks it complete (or when
    the reported total has reached the whole batch).
    """
    norm = " ".join(_ANSI_ESCAPE_RE.sub("", line).strip().split())
    token = next((candidate for candidate in tokens if candidate in norm), "")
    if not token:
        return None

    count_match = re.search(r"\|\s*(\d+)\s+of\s+(\d+)", norm)
    done = int(count_match.group(1)) if count_match else 0
    reported_total = int(count_match.group(2)) if count_match else 0
    progress_total = max(1, expected_total, reported_total)
    visibly_complete = "✔" in norm or "✓" in norm

    if "[-" in norm:
        event = "queued"
    elif visibly_complete or (
        reported_total >= max(1, expected_total)
        and done >= reported_total
    ):
        event = "done"
    else:
        event = "start"

    if event == "done":
        fraction = 1.0
    elif count_match:
        fraction = max(0.0, min(1.0, done / progress_total))
    else:
        fraction = 0.0

    return {
        "norm": norm,
        "process": token,
        "event": event,
        "done": done,
        "reported_total": reported_total,
        "progress_total": progress_total,
        "fraction": fraction,
        "has_count": bool(count_match),
    }


def _run(
    cmd: list[str],
    *,
    label: str,
    on_line=None,
    display_cmd: list[str] | None = None,
    cwd: Path | None = None,
) -> None:
    """Stream a subprocess's stdout/stderr into this worker's stdout
    (which is already redirected to log.txt by dragen_jobs.start_job).
    Raises on non-zero exit so the outer try/except records failure.
    """
    started = time.monotonic()
    _log()
    _log(f"========================= [{label}] =========================")
    _log("$ " + " ".join(display_cmd or cmd))
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
        bufsize=1,
        cwd=str(cwd) if cwd is not None else None,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        suffix = ""
        if on_line is not None:
            suffix = on_line(line) or ""
        if suffix:
            print(line.rstrip("\n") + suffix, flush=True)
        else:
            print(line, end="", flush=True)
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
    ap.add_argument(
        "--research-only", "--with-extra-vep",
        dest="research_only", action="store_true",
    )
    ap.add_argument("--without-pgx", action="store_true")
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
    nf_launch, nf_work = _nextflow_context(samples)
    pipeline_type = "dragen" if mode == "dragen" else "nckuh"
    legacy_staging = os.environ.get("NGS_UI_TERTIARY_LEGACY_STAGING", "").strip().lower() in {
        "1", "true", "yes", "y", "on"
    }
    lock_handles: list[object] = []
    nextflow_cache_lock: object | None = None
    staging_root: Path | None = None
    rollback_root: Path | None = None
    seeded_resume_session: str | None = None
    previous_main_sigterm = signal.getsignal(signal.SIGTERM)

    def cleanup_on_sigterm(_signum, _frame):
        raise SystemExit(143)

    signal.signal(signal.SIGTERM, cleanup_on_sigterm)

    started_at = _now()
    _set_step(job_id, "detect-pipeline-output", state="running", started_at=started_at)
    try:
        lock_handles = _acquire_sample_locks(
            sample_ids + [sample["source_sample_id"] for sample in samples]
        )
        staging_root, rollback_root = _prepare_job_staging(job_id)
        nf_launch.mkdir(parents=True, exist_ok=True)
        nf_work.mkdir(parents=True, exist_ok=True)
        _update(
            job_id,
            staging_root=str(staging_root),
            nextflow_launch_dir=str(nf_launch),
            nextflow_work_dir=str(nf_work),
            nextflow_cache_scope=pipeline_type,
        )

        # Existing live output is recorded for operator context only.
        # It never decides whether Nextflow runs; dependency/cache decisions
        # belong to Nextflow -resume.
        existing_by_sid: dict[str, dict[str, Path]] = {}
        pending_samples: list[dict] = list(samples)
        for sample in samples:
            sid = sample["sample_id"]
            source_sid = sample["source_sample_id"]
            outputs, missing, found_under = _pipeline_outputs_for(
                sid,
                source_sid,
                require_pgx=not args.without_pgx,
                require_str=True,
            )
            present = ", ".join(sorted(outputs)) or "none"
            status = "complete" if not missing else f"missing {', '.join(missing)}"
            _log(
                f"[detect] {sid}: live output {status} under "
                f"{PIPELINE_OUT_ROOT / (found_under or sid)} "
                f"(present: {present}); Nextflow -resume will run"
            )

        if pending_samples:
            _log(f"[detect] running nextflow for {len(pending_samples)} sample(s)")

            if pending_samples:
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

                # 3. Nextflow → job-private staging/<source SID>/...
                _set_step(job_id, "waiting-nextflow-cache")
                nextflow_stages = [
                    ("prepare-vcf-dragen-add-tag", ("PREPARE_VCF_DRAGEN:ADD_DRAGEN_TAG",), 0.5),
                    ("prepare-vcf-dragen", ("PREPARE_VCF_DRAGEN",), 0.5),
                    ("prepare-vcf", ("PREPARE_VCF",), 1.0),
                    ("mito-vep", ("MITO_ANNOTATE:MITO_VEP",), 1.5),
                    ("mito-parse", ("MITO_ANNOTATE:MITO_PARSE",), 0.8),
                    ("str-parse-dragen", ("STR_PARSE_DRAGEN",), 0.7),
                    ("prepare-cnv-dragen", ("PREPARE_CNV_DRAGEN", "PREPARE_CNV_NCKUH"), 0.7),
                    ("annotsv-cnv-dragen", ("ANNOTSV_CNV_DRAGEN", "ANNOTSV_CNV_NCKUH"), 3.0),
                    ("prepare-sv-dragen", ("PREPARE_SV_DRAGEN", "PREPARE_SV_NCKUH"), 0.7),
                    ("annotsv-sv-dragen", ("ANNOTSV_SV_DRAGEN", "ANNOTSV_SV_NCKUH"), 8.0),
                    ("vep-annotate", ("SNV_ANNOTATE:VEP_ANNOTATE",), 32.0),
                    ("pangolin-score", ("SNV_ANNOTATE:PANGOLIN_SCORE",), 8.0),
                    ("parse-csq", ("PARSE_VEP_CSQ:PARSE_CSQ",), 16.0),
                    ("acmg-classify", ("ACMG_CLASSIFY",), 8.0),
                    ("pgx-stellarpgx", ("PGX_ANNOTATE:PGX_STELLARPGX",), 1.0),
                    ("pgx-hla-extract", ("PGX_ANNOTATE:PGX_HLA_EXTRACT",), 1.0),
                    ("pgx-optitype", ("PGX_ANNOTATE:PGX_OPTITYPE",), 8.0),
                    ("pgx-pharmcat", ("PGX_ANNOTATE:PGX_PHARMCAT",), 4.0),
                    ("pgx-parse", ("PGX_ANNOTATE:PGX_PARSE",), 2.0),
                ]
                if args.without_pgx:
                    nextflow_stages = [row for row in nextflow_stages if not row[0].startswith("pgx-")]
                nextflow_stage_index = {slug: idx for idx, (slug, _tokens, _weight) in enumerate(nextflow_stages)}
                nextflow_progress_start = 3.0
                nextflow_progress_end = 82.0
                nextflow_total_weight = max(1.0, sum(weight for _slug, _tokens, weight in nextflow_stages))
                nextflow_stage_fraction: dict[str, float] = {
                    slug: 0.0 for slug, _tokens, _weight in nextflow_stages
                }
                nextflow_progress_rank = -1
                nextflow_running: dict[str, float] = {}
                nextflow_seen_events: set[tuple[str, str, int, int]] = set()

                def track_nextflow(line: str) -> None:
                    nonlocal nextflow_progress_rank
                    if "[" not in line or "]" not in line:
                        return
                    for slug, tokens, weight in nextflow_stages:
                        parsed = _parse_nextflow_progress_line(
                            line,
                            tokens,
                            expected_total=len(pending_samples),
                        )
                        if parsed is None:
                            continue
                        process = str(parsed["process"])
                        event = str(parsed["event"])
                        done = int(parsed["done"])
                        reported_total = int(parsed["reported_total"])
                        progress_total = int(parsed["progress_total"])
                        fraction = float(parsed["fraction"])
                        has_count = bool(parsed["has_count"])
                        key = (slug, event, done, reported_total)
                        if key in nextflow_seen_events:
                            return
                        nextflow_seen_events.add(key)
                        if event == "start":
                            nextflow_running.setdefault(slug, time.monotonic())
                        elapsed = None
                        if event == "done":
                            started = nextflow_running.get(slug)
                            elapsed = (time.monotonic() - started) if started is not None else None
                        _record_nextflow_step(
                            job_id,
                            slug,
                            process,
                            event,
                            elapsed=elapsed,
                            done=done or None,
                            total=reported_total or None,
                        )
                        suffix = ""
                        elapsed_label = _elapsed_minutes_label(elapsed)
                        if elapsed_label:
                            suffix = f"  {elapsed_label} [{_now()}]"
                        if event == "queued":
                            return
                        stage_idx = nextflow_stage_index.get(slug, 0)
                        nextflow_stage_fraction[slug] = max(
                            nextflow_stage_fraction.get(slug, 0.0),
                            max(0.0, min(1.0, fraction)),
                        )
                        weighted_done = sum(
                            nextflow_stage_fraction.get(stage_slug, 0.0) * stage_weight
                            for stage_slug, _stage_tokens, stage_weight in nextflow_stages
                        )
                        pct = nextflow_progress_start + (
                            weighted_done
                            / nextflow_total_weight
                        ) * (nextflow_progress_end - nextflow_progress_start)
                        if stage_idx > nextflow_progress_rank or event == "done" or has_count:
                            nextflow_progress_rank = max(nextflow_progress_rank, stage_idx)
                            _update(
                                job_id,
                                nextflow_progress_pct=round(pct, 1),
                                nextflow_current={
                                    "step": slug,
                                    "process": process,
                                    "event": event,
                                    "done": done,
                                    "total": progress_total,
                                    "reported_total": reported_total,
                                },
                            )
                        return suffix

                _log(
                    f"[nextflow] waiting for shared {pipeline_type} cache lineage"
                )
                nextflow_cache_lock = _acquire_nextflow_cache_lock(mode)
                _log(f"[nextflow] acquired shared {pipeline_type} cache lock")
                seeded_context = _seed_shared_nextflow_context(
                    samples,
                    launch_dir=nf_launch,
                )
                if seeded_context is not None:
                    seeded_resume_session, nf_work = seeded_context
                    _write_shared_work_pointer(nf_launch, nf_work)
                else:
                    # A job may have waited while another worker initialized
                    # this lineage. Re-read its persisted work root after the
                    # mode lock is acquired instead of using the stale path
                    # calculated before the wait.
                    _same_launch, nf_work = _nextflow_context(samples)
                    nf_work.mkdir(parents=True, exist_ok=True)
                    _write_shared_work_pointer(nf_launch, nf_work)
                _update(
                    job_id,
                    shared_resume_session=seeded_resume_session,
                    nextflow_work_dir=str(nf_work),
                )
                if seeded_resume_session:
                    _log(
                        "[nextflow] adopted previous cache session "
                        f"{seeded_resume_session} into shared {pipeline_type} lineage"
                    )
                _set_step(job_id, "nextflow")

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
                        "--out_dir",   str(staging_root),
                    ]
                    if args.without_pgx:
                        nextflow_cmd += ["--run_pgx", "false"]
                    if args.research_only:
                        nextflow_cmd += ["--academic_dbnsfp", "true"]
                    nextflow_cmd.append("-resume")
                    if seeded_resume_session:
                        nextflow_cmd.append(seeded_resume_session)
                    nextflow_run_cmd = nextflow_cmd
                    nextflow_run_label = "2b/4 nextflow legacy"
                else:
                    samplesheet = TERTIARY_JOBS_DIR / job_id / "samplesheet.csv"
                    nextflow_cmd = [
                        "nextflow",
                        "-c", str(TERTIARY_NEXTFLOW_CONFIG),
                        "run", "/home/pipeline/tertiary_code/main_tertiary.nf",
                        "-profile", "dgm",
                        "-work-dir", str(nf_work),
                        "--pipeline_type", pipeline_type,
                        "--samplesheet", str(samplesheet),
                        "--out_dir", str(staging_root),
                    ]
                    if args.without_pgx:
                        nextflow_cmd += ["--run_pgx", "false"]
                    if args.research_only:
                        nextflow_cmd += ["--academic_dbnsfp", "true"]
                    nextflow_cmd.append("-resume")
                    if seeded_resume_session:
                        nextflow_cmd.append(seeded_resume_session)
                    inner = " ".join(shlex.quote(part) for part in nextflow_cmd)
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
                    nextflow_run_cmd = [
                        "bash",
                        "-lc",
                        shell_cmd,
                    ]
                    nextflow_run_label = "2b/4 nextflow v3.x"

                try:
                    _run(
                        nextflow_run_cmd,
                        label=nextflow_run_label,
                        on_line=track_nextflow,
                        cwd=nf_launch,
                    )
                finally:
                    _release_nextflow_cache_lock(nextflow_cache_lock)
                    nextflow_cache_lock = None

                for sample in pending_samples:
                    sid = sample["sample_id"]
                    source_sid = sample["source_sample_id"]
                    staged_pipeline_id = sid if legacy_staging else source_sid
                    outputs = _validate_staged_sample(
                        staging_root,
                        staged_pipeline_id,
                        require_pgx=not args.without_pgx,
                    )
                    existing_by_sid[sid] = outputs
        # 4. Prepare job-private 08_postprocessing beside staged 00-07.
        _set_step(job_id, "prepare-postprocessing")
        pipeline_annotsv_available: dict[str, set[str]] = {}
        raw_tsv_by_sid: dict[str, Path] = {}
        staged_sample_dir_by_sid: dict[str, Path] = {}
        final_raw_tsv_by_sid: dict[str, Path] = {}
        for sample in samples:
            sid = sample["sample_id"]
            source_sid = sample["source_sample_id"]
            source_vcf = sample["vcf_path"]
            staged_pipeline_id = sid if legacy_staging else source_sid
            staged_sample_dir = staging_root / staged_pipeline_id
            staged_sample_dir_by_sid[sid] = staged_sample_dir
            post_dir = staged_sample_dir / sample_layout.POSTPROCESSING_DIRNAME
            post_dir.mkdir(parents=True, exist_ok=True)
            outputs = existing_by_sid.get(sid) or {}
            existing = outputs.get("snv_indel.acmg.tsv")
            if existing is None:
                raise RuntimeError(f"internal error: no pipeline TSV for {sid}")
            _validate_acmg_tsv(
                existing,
                strict_v31=not legacy_staging,
                expect_academic_dbnsfp=(args.research_only if not legacy_staging else None),
            )
            if not legacy_staging:
                sidecar = _ensure_v36_annotation_versions(existing)
                _log(f"[source] {sid}: ClinVar baseline metadata {sidecar}")
            raw_tsv_by_sid[sid] = existing
            final_raw_tsv = (
                sample_layout.unified_sample_dir(sid)
                / existing.relative_to(staged_sample_dir)
            )
            final_raw_tsv_by_sid[sid] = final_raw_tsv
            _track_pipeline_source(
                sid,
                post_dir,
                existing,
                source_sample_id=source_sid,
                source_vcf_path=source_vcf,
                pipeline_type=mode,
                published_source=final_raw_tsv,
            )
            _log(f"[source] {sid}: staged immutable SNV TSV {existing}")
            ploidy_copy, ploidy_qc_copy = _copy_ploidy_artifacts(
                mode=mode,
                source_vcf=source_vcf,
                source_sample_id=source_sid,
                sample_id=sid,
                post_dir=post_dir,
            )
            if ploidy_copy is None:
                _log(f"[copy] {sid}: matching {mode} ploidy VCF not found for {source_vcf}")
            else:
                ploidy_src, ploidy_dst = ploidy_copy
                _log(f"[copy] {sid}: {ploidy_src} → {ploidy_dst}")
            if ploidy_qc_copy is not None:
                ploidy_qc_src, ploidy_qc_dst = ploidy_qc_copy
                _log(f"[copy] {sid}: future-use ploidy QC {ploidy_qc_src} → {ploidy_qc_dst}")
            mito_src = outputs.get("mito.tsv")
            if mito_src is None:
                _log(
                    f"[copy] {sid}: mito output not found under staged "
                    f"{staged_sample_dir}/04_mito/"
                )
            else:
                mito_dst = sample_layout.scoped_file(
                    post_dir,
                    sid,
                    "mito.annotated.tsv",
                    for_write=True,
                    force_prefixed=True,
                )
                shutil.copyfile(mito_src, mito_dst)
                try:
                    mito_changed = mitomap_mito.annotate_mito_tsv(mito_dst)
                except Exception:
                    mito_changed = False
                if mito_changed:
                    _log(f"[derived] {sid}: MITOMAP-enriched mito TSV → {mito_dst}")
                else:
                    mito_dst.unlink(missing_ok=True)
                    _log(f"[source] {sid}: use pipeline mito TSV directly ({mito_src})")
            str_src = outputs.get("str.tsv")
            if str_src is None:
                _log(
                    f"[copy] {sid}: STR output not found under staged "
                    f"{staged_sample_dir}/05_str/"
                )
            else:
                _log(f"[source] {sid}: use pipeline STR TSV directly ({str_src})")
            pgx_files = _find_pipeline_pgx_files(
                staged_pipeline_id,
                root=staging_root,
            )
            if not pgx_files and not args.without_pgx:
                _log(
                    f"[copy] {sid}: PGx outputs not found under staged "
                    f"{staged_sample_dir}/07_pgx/"
                )
            for dst_name, pgx_src in sorted(pgx_files.items()):
                _log(f"[source] {sid}: use pipeline {dst_name} directly ({pgx_src})")
            for kind in ("cnv", "sv"):
                key = f"{kind}.annotated.tsv"
                annotsv_src = outputs.get(key)
                if annotsv_src is None:
                    _log(
                        f"[copy] {sid}: {kind.upper()} AnnotSV output not found "
                        f"under staged {staged_sample_dir}/06_cnv_sv/"
                    )
                    continue
                pipeline_annotsv_available.setdefault(sid, set()).add(kind)
                _log(f"[source] {sid}: use pipeline {kind.upper()} TSV directly ({annotsv_src})")

        # 5. Post-processing chain. ClinVar is compared against the weekly UI
        # snapshot here; the fixed pipeline annotation remains immutable.
        _set_step(job_id, "post-processing")
        def track_post_processing(line: str) -> None:
            match = re.search(r"\[(post-processing-step|sample-step)]\s+([a-z0-9-]+)\s+start", line)
            if match:
                group = "post-processing" if match.group(1) == "post-processing-step" else "sample-step"
                _set_step(job_id, f"{group}:{match.group(2)}")

        post_processing_count = len(samples)
        for post_processing_index, sample in enumerate(samples):
            sid = sample["sample_id"]
            _update(
                job_id,
                post_processing_sample_index=post_processing_index,
                post_processing_sample_count=post_processing_count,
                post_processing_sample_id=sid,
                # Backward-compatible keys for older frontends/state readers.
                stopgap_sample_index=post_processing_index,
                stopgap_sample_count=post_processing_count,
                stopgap_sample_id=sid,
            )
            raw_tsv = raw_tsv_by_sid[sid]
            staged_sample_dir = staged_sample_dir_by_sid[sid]
            post_dir = staged_sample_dir / sample_layout.POSTPROCESSING_DIRNAME
            for stale_work in post_dir.glob(".snv_indel.*.working.tsv"):
                stale_work.unlink(missing_ok=True)
            work_tsv = post_dir / f".snv_indel.{job_id}.working.tsv"
            previous_sigterm = signal.getsignal(signal.SIGTERM)

            def cleanup_work_on_sigterm(_signum, _frame):
                work_tsv.unlink(missing_ok=True)
                raise SystemExit(143)

            signal.signal(signal.SIGTERM, cleanup_work_on_sigterm)
            try:
                shutil.copyfile(raw_tsv, work_tsv)
            except BaseException:
                work_tsv.unlink(missing_ok=True)
                signal.signal(signal.SIGTERM, previous_sigterm)
                raise
            stop_args = [str(scripts / "run_stopgaps.sh"),
                         "--work-tsv", str(work_tsv),
                         "--raw-tsv", str(raw_tsv),
                         "--post-dir", str(post_dir),
                         "--sample", sid,
                         "--seq-type", sample["seq_type"]]
            if pipeline_annotsv_available.get(sid) == {"cnv", "sv"}:
                stop_args += ["--skip-cnv"]
                _log(f"[post-processing] {sid}: skip AnnotSV fallback; pipeline CNV/SV are available")
            elif mode == "dragen":
                stop_args += ["--dragen-cnv-source", sample["vcf_path"]]
            elif mode == "inhouse":
                if sample["cnv_vcf"]:
                    stop_args += ["--inhouse-cnv-vcf", sample["cnv_vcf"]]
                if sample["sv_vcf"]:
                    stop_args += ["--inhouse-sv-vcf",  sample["sv_vcf"]]
            if not args.research_only:
                stop_args.append("--skip-spliceai")
            display_stop_args = [
                "post-processing",
                "--raw-tsv", str(raw_tsv),
                "--post-dir", str(post_dir),
                "--sample", sid,
                "--seq-type", sample["seq_type"],
            ]
            if "--skip-cnv" in stop_args:
                display_stop_args.append("--skip-cnv")
            if "--skip-spliceai" in stop_args:
                display_stop_args.append("--skip-spliceai")
            try:
                _run(
                    stop_args,
                    label=f"post-processing {sid}",
                    on_line=track_post_processing,
                    display_cmd=display_stop_args,
                )
                _rebase_staged_derived_paths(
                    sample_id=sid,
                    stage_post_dir=post_dir,
                    final_raw_tsv=final_raw_tsv_by_sid[sid],
                    final_post_dir=sample_layout.unified_postprocessing_dir(sid),
                )
                sample_layout.write_layout_marker_in_sample_dir(
                    staged_sample_dir,
                    sid,
                    source_id=sample["source_sample_id"],
                    raw_tsv=raw_tsv,
                )
            finally:
                work_tsv.unlink(missing_ok=True)
                signal.signal(signal.SIGTERM, previous_sigterm)

        # 6. Promote every fully prepared generation. Keep every backup until
        # the whole batch has switched so a later sample failure can roll back
        # earlier samples from the same job.
        _set_step(job_id, "promote-output")
        promoted_batches: list[list[tuple[Path | None, Path, Path]]] = []
        try:
            for sample in samples:
                sid = sample["sample_id"]
                operations = _promote_staged_generation(
                    sample_id=sid,
                    staged_sample_dir=staged_sample_dir_by_sid[sid],
                    live_sample_dir=sample_layout.unified_sample_dir(sid),
                    rollback_sample_dir=rollback_root / sid,
                )
                promoted_batches.append(operations)
                live_post = sample_layout.unified_postprocessing_dir(sid)
                for stale_work in live_post.glob(".snv_indel.*.working.tsv"):
                    stale_work.unlink(missing_ok=True)
                _log(f"[promote] {sid}: staged generation activated")
        except BaseException:
            for operations in reversed(promoted_batches):
                _rollback_promotion_operations(operations)
            raise

        from ..services import sample_loader
        for sid in sample_ids:
            try:
                sample_loader.invalidate_sample_cache(
                    sample_layout.unified_postprocessing_dir(sid)
                )
                sample_loader.update_case_table_row(sid)
            except Exception as cache_error:
                _log(f"[cache] {sid}: refresh failed: {cache_error}")

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
    finally:
        signal.signal(signal.SIGTERM, previous_main_sigterm)
        _release_nextflow_cache_lock(nextflow_cache_lock)
        _release_sample_locks(lock_handles)
        for path in (staging_root, rollback_root):
            if path is not None and path.is_dir():
                shutil.rmtree(path, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
