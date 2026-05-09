"""Hydrate a per-patient sample directory the pipeline already produced.

The tertiary pipeline lands variant calls in
    tertiary_output/{LIS_ID}/snv_indel.annotated.tsv
on its own. The 載入新個案 flow attaches reviewer-side info on top:
basic identifiers, an empty default analysis, optionally a parsed
copy of the phenotype.txt. After hydration the directory looks like:

    tertiary_output/{LIS_ID}/
      snv_indel.annotated.tsv  (untouched; pipeline output)
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
from datetime import datetime, timezone
from pathlib import Path

from ..config import TERTIARY_OUTPUT_ROOT
from . import analyses_store, emr_client, phenotype_io, vcf_writer


_LIS_ID_RE = re.compile(r"^[A-Za-z0-9_-]{1,32}$")
_TEST_TYPES = {"WES", "WGS"}
_GENOME_BUILDS = {"hg19", "hg38"}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _validate_lis_id(lis_id: str) -> None:
    if not _LIS_ID_RE.match(lis_id or ""):
        raise ValueError(
            "lis_id must match [A-Za-z0-9_-]{1,32} (used as directory name)"
        )


def sample_exists(lis_id: str) -> bool:
    return (TERTIARY_OUTPUT_ROOT / lis_id).is_dir()


def is_registered(lis_id: str) -> bool:
    return (TERTIARY_OUTPUT_ROOT / lis_id / "sample_metadata.json").is_file()


def register(
    *,
    lis_id: str,
    name: str,
    mrn: str,
    sex: str = "",
    test_type: str = "WES",
    genome_build: str = "hg38",
    category: str = "",
    phenotype_text: str = "",
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

    The directory tertiary_output/{lis_id}/ must already exist with at
    least snv_indel.annotated.tsv inside (the pipeline drops it there).
    Refuses if the dir is already registered (sample_metadata.json
    present).
    """
    _validate_lis_id(lis_id)
    if not name:
        raise ValueError("name is required")
    if not mrn:
        raise ValueError("mrn is required")
    if test_type not in _TEST_TYPES:
        raise ValueError(f"test_type must be one of {sorted(_TEST_TYPES)}")
    if genome_build not in _GENOME_BUILDS:
        raise ValueError(f"genome_build must be one of {sorted(_GENOME_BUILDS)}")

    sample_dir = TERTIARY_OUTPUT_ROOT / lis_id
    if not sample_dir.is_dir():
        raise FileNotFoundError(
            f"pipeline directory not found: {sample_dir} "
            "(tertiary pipeline drops the TSV here; nothing to register yet)"
        )
    if not (sample_dir / "snv_indel.annotated.tsv").is_file():
        raise FileNotFoundError(
            f"snv_indel.annotated.tsv missing under {sample_dir}"
        )
    if (sample_dir / "sample_metadata.json").is_file():
        raise FileExistsError(f"sample already registered: {lis_id}")

    # Parse the reviewer-curated phenotype.txt first; if it had any
    # content treat that as authoritative. Otherwise fall back to the
    # EMR's GetPhenotypeList output (best-effort: reviewer txt wins
    # per the system convention). Frontend-edited chips override
    # both — they were derived from one of these sources and may
    # have been edited.
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

    # Generate the minimal VCF Exomiser/LIRICAL will consume. The path
    # is convention-driven (tertiary_output/{lis_id}/vcf_from_tsv.vcf.gz)
    # so we don't ask the reviewer to fill it in.
    vcf_out = vcf_writer.from_tsv(lis_id)

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
        "vcf_path":             str(vcf_out),
        "run_date":             now,
        "active_analysis":      "default",
        "clinical_description": "",
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
    (sample_dir / "sample_metadata.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    # Default analysis.json + audit copy of the parsed phenotype.txt.
    # write_version side-effects pheno_score.tsv into the version dir,
    # so the freshly-registered sample is immediately ready for the
    # Clinical-block / pheno-score lookups (no need to wait for the
    # reviewer to hit "save" in the analysis page).
    analyses_store.write_version(lis_id, "default", hpo=hpo, panels=panels)
    if hpo or panels:
        phenotype_io.write(
            hpo, panels,
            analyses_store.version_dir(lis_id, "default")
            / f"{lis_id}_{mrn}_phenotype.txt",
        )
        # Reflect the just-computed pheno set onto the SNV TSV's
        # IN_PANEL column so per-sample loads see the right markers
        # immediately. (No-op when HPO+panels are both empty.)
        from . import phenotype_scorer
        scores = phenotype_scorer.compute_pheno_score(hpo or [], panels or [])
        phenotype_scorer.update_in_panel_column(
            lis_id, {g for g, s in scores.items() if s > 0}
        )

    return meta
