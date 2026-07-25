"""Hydrate a per-patient sample directory the pipeline already produced.

The tertiary pipeline lands immutable calls in 03_acmg and the UI attaches
reviewer-side info below 08_postprocessing:
basic identifiers, an empty default analysis, optionally a parsed
copy of the phenotype.txt. After hydration the directory looks like:

    tertiary_output/{LIS_ID}/08_postprocessing/
      (03_acmg raw TSV lives at the sample root)
      sample_metadata.json     (basic info + empty reviewer state)
      analyses/default/
        analysis.json          (hpo + selected_panels + note)
        {LIS_ID}_{MRN}_phenotype.txt   (audit copy, when provided)

Refusal cases:
  * lis_id directory missing or no TSV → 404 / 400 from the router
  * sample_metadata.json already present → 409 (already registered)
"""
from __future__ import annotations

import json
import re
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path

from . import analyses_store, emr_client, phenotype_io, sample_layout, test_types, vcf_writer


_LIS_ID_RE = re.compile(r"^[A-Za-z0-9_-]{1,32}$")
_TEST_TYPES = test_types.VALID_TEST_TYPES
_GENOME_BUILDS = {"hg19", "hg38"}


def _log_perf(event: str, started: float, **fields) -> None:
    parts = [f"[perf] {event}", f"elapsed={time.perf_counter() - started:.3f}s"]
    parts.extend(f"{key}={value}" for key, value in fields.items())
    print(" ".join(parts), flush=True)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _validate_lis_id(lis_id: str) -> None:
    if not _LIS_ID_RE.match(lis_id or ""):
        raise ValueError(
            "lis_id must match [A-Za-z0-9_-]{1,32} (used as directory name)"
        )


def sample_exists(lis_id: str) -> bool:
    return sample_layout.is_ui_ready(lis_id)


def is_registered(lis_id: str) -> bool:
    return sample_layout.state_file(lis_id, "sample_metadata.json").is_file()


def delete(lis_id: str, *, delete_pipeline_output: bool = False) -> dict:
    """Delete or unregister one sample.

    By default this only removes reviewer-side registration/report state,
    preserving pipeline outputs so the sample appears in the unregistered
    list again. Passing delete_pipeline_output keeps the older destructive
    behavior used by pipeline-output management.
    """
    _validate_lis_id(lis_id)
    ui_dir = sample_layout.state_dir(lis_id)
    if not ui_dir.is_dir():
        raise FileNotFoundError(f"sample not found: {lis_id}")

    if delete_pipeline_output:
        from . import dragen_jobs
        result = dragen_jobs.delete_pipeline_output(lis_id)
        return {
            **result,
            "pipeline_output_requested": True,
            "pipeline_output_deleted": bool(result["deleted"]),
            "pipeline_output_error": "",
        }

    deleted = []
    for name in ("sample_metadata.json", "case_summary.json"):
        for path in sample_layout.state_file_candidates(lis_id, name):
            if path.exists():
                path.unlink()
                deleted.append(str(path))

    analyses_dir = ui_dir / "analyses"
    if analyses_dir.is_dir():
        shutil.rmtree(analyses_dir)
        deleted.append(str(analyses_dir))

    from . import sample_loader
    sample_loader.invalidate_sample_cache(ui_dir)
    sample_loader.remove_case_table_row(lis_id)
    return {
        "sample_id": lis_id,
        "deleted": deleted,
        "unregistered": True,
        "pipeline_output_requested": False,
        "pipeline_output_deleted": False,
        "pipeline_output_error": "",
    }


def register(
    *,
    lis_id: str,
    name: str,
    mrn: str,
    sex: str = "",
    test_type: str = "WES",
    genome_build: str = "hg38",
    category: str = "",
    department: str = "",
    physician: str = "",
    sign_received_at: str = "",
    phenotype_text: str = "",
    clinical_description: str = "",
    hpo: list | None = None,
    panels: list | None = None,
) -> dict:
    """Attach reviewer-side info to a pipeline-produced directory.

    Phenotype precedence:
      1. Explicit `hpo` / `panels` lists (frontend-edited chips win;
         already include any EMR-sourced terms the reviewer kept).
      2. Else parse `phenotype_text` (the reviewer's
         <LIS>_<MRN>_phenotype.txt content).
      3. Else fall back to the EMR phenotype API.
    """
    """Attach reviewer-side info to a pipeline-produced directory.

    The unified layout marker + 03_acmg raw TSV must already exist.
    Refuses if the dir is already registered (sample_metadata.json
    present).
    """
    started = time.perf_counter()
    _validate_lis_id(lis_id)
    test_type = test_types.normalize_test_type(test_type, sample_id=lis_id)
    if not name:
        raise ValueError("name is required")
    if not mrn:
        raise ValueError("mrn is required")
    if test_type not in _TEST_TYPES:
        raise ValueError(f"test_type must be one of {sorted(_TEST_TYPES)}")
    if genome_build not in _GENOME_BUILDS:
        raise ValueError(f"genome_build must be one of {sorted(_GENOME_BUILDS)}")

    sample_dir = sample_layout.state_dir(lis_id)
    raw_tsv = sample_layout.snv_raw_tsv(lis_id)
    if not sample_layout.is_ui_ready(lis_id) or not raw_tsv.is_file():
        raise FileNotFoundError(
            f"pipeline SNV TSV not found for {lis_id} "
            "(tertiary pipeline drops the TSV here; nothing to register yet)"
        )
    sample_dir.mkdir(parents=True, exist_ok=True)
    if sample_layout.state_file(lis_id, "sample_metadata.json").is_file():
        raise FileExistsError(f"sample already registered: {lis_id}")

    # Parse the reviewer-curated phenotype.txt first; if it had any
    # content treat that as authoritative. Otherwise fall back to the
    # EMR's GetPhenotypeList output (best-effort: reviewer txt wins
    # per the system convention). Frontend-edited chips override
    # both — they were derived from one of these sources and may
    # have been edited.
    pheno_started = time.perf_counter()
    if hpo is not None or panels is not None:
        hpo = list(hpo or [])
        panels = list(panels or [])
        emr_payload = emr_client.fetch(mrn) if mrn else {}
    else:
        hpo, panels = phenotype_io.parse(phenotype_text or "")
        emr_payload = emr_client.fetch(mrn) if mrn else {}
        if not hpo and not panels:
            emr_pheno = emr_payload.get("phenotype") or {}
            if emr_pheno.get("found"):
                hpo = emr_pheno.get("hpo") or []
    _log_perf(
        "patient_store.register.phenotype_emr",
        pheno_started,
        sample=lis_id,
        hpo=len(hpo or []),
        panels=len(panels or []),
    )

    # Sex / dob / genetic_counseling come from the consultation API.
    # Sex from EMR overwrites whatever the reviewer typed (per spec);
    # genetic_counseling lands as-is. Failures are silent — feature
    # disabled / empty consultation just means these fields stay blank.
    consult = emr_payload.get("consultation") or {}
    if consult.get("sex"):
        sex = consult["sex"]                           # overwrite
    if consult.get("date_of_birth"):
        dob_from_emr = consult["date_of_birth"]
    else:
        dob_from_emr = ""
    genetic_counseling = consult.get("text", "") or ""

    # Exomiser/LIRICAL uses a convention-driven VCF path under the
    # sample directory. Building it scans the full raw TSV on WGS, so
    # registration only records an already-existing VCF. The worker
    # creates or refreshes it in the background before analysis.
    vcf_started = time.perf_counter()
    vcf_out = vcf_writer.vcf_path_for(lis_id)
    vcf_path = str(vcf_out) if vcf_out.is_file() else ""
    _log_perf(
        "patient_store.register.vcf_path",
        vcf_started,
        sample=lis_id,
        exists=int(bool(vcf_path)),
    )

    # Seed sample_metadata.json with basic info + empty reviewer state.
    now = _now()
    meta = {
        "sample_id":            lis_id,
        "lis_id":               lis_id,
        "name":                 name,
        "mrn":                  mrn,
        "sex":                  (sex or "").upper() if sex else "",
        "date_of_birth":        dob_from_emr,
        "test_type":            test_type,
        "genome_build":         genome_build,
        "category":             category or "",
        "department":           department or "",
        "physician":            physician or "",
        "sign_received_at":     sign_received_at or "",
        "vcf_path":             vcf_path,
        "run_date":             now,
        "active_analysis":      "default",
        "clinical_description": clinical_description or "",
        "genetic_counseling":   genetic_counseling,
        "emr_synced_at":        now if (consult.get("found") or emr_payload.get("phenotype", {}).get("found")) else "",
        "comment":              "",
        "tags":                 [],
        "status":               {},
        "edits":                {},
        "panels":               {},
        "manual_variants":      [],
        "created_at":           now,
        "updated_at":           now,
    }
    meta_started = time.perf_counter()
    sample_layout.state_file(
        lis_id,
        "sample_metadata.json",
        for_write=True,
    ).write_text(
        json.dumps(meta, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    _log_perf("patient_store.register.metadata", meta_started, sample=lis_id)

    # Default analysis.json + audit copy of the parsed phenotype.txt.
    # write_version side-effects pheno_score.tsv into the version dir,
    # so the freshly-registered sample is immediately ready for the
    # Clinical-block / pheno-score lookups (no need to wait for the
    # reviewer to hit "save" in the analysis page).
    version_started = time.perf_counter()
    analyses_store.write_version(lis_id, "default", hpo=hpo, panels=panels)
    _log_perf(
        "patient_store.register.write_version",
        version_started,
        sample=lis_id,
        hpo=len(hpo or []),
        panels=len(panels or []),
    )
    if hpo or panels:
        phenotype_io.write(
            hpo, panels,
            analyses_store.version_dir(lis_id, "default")
            / f"{lis_id}_{mrn}_phenotype.txt",
        )
        # pheno_score.tsv is the source of truth for in-panel state.
        # SNV loads and gene search apply it dynamically, so registration
        # does not rewrite the large raw TSV.
        from . import phenotype_scorer
        score_started = time.perf_counter()
        scores = phenotype_scorer.compute_pheno_score(hpo or [], panels or [])
        _log_perf(
            "patient_store.register.compute_pheno_score",
            score_started,
            sample=lis_id,
            genes=len(scores),
        )

    _log_perf("patient_store.register.total", started, sample=lis_id)
    try:
        from . import sample_loader
        sample_loader.update_case_table_row(lis_id)
    except Exception as e:
        print(f"[case-table] register refresh failed for {lis_id}: {e}", flush=True)
    return meta
