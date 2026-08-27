"""HPO/panel editor + Python pheno_score recompute (Phase A).

Endpoints:
  GET  /api/hpo/search?q=...&limit=20
  GET  /api/panels
  POST /api/samples/{sample_id}/phenotype
       body: {"hpo": [{"phenotype": "HP:0001250", "label": "...", "weight": 2}, ...],
              "panels": ["HIE", "Marfan_panel", ...]}
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from fastapi import APIRouter, Depends, HTTPException, Query

from ..auth import current_user
from ..services import (
    hpo_ontology,
    patient_phenotype_store,
    phenotype_scorer,
    sample_layout,
    sample_loader,
)

router = APIRouter(prefix="/api", tags=["phenotype"], dependencies=[Depends(current_user)])


@router.get("/hpo/search")
def hpo_search(q: str = Query(""), limit: int = Query(20, ge=1, le=100)):
    results = hpo_ontology.search(q, limit=limit)
    # Annotate with the per-term gene count from phenotype_to_genes.txt so
    # the picker can show "Seizure (84 genes)" without a second round-trip.
    for r in results:
        r["gene_count"] = phenotype_scorer.gene_count(r["hpo_id"])
        for parent in r.get("parents", []):
            parent["gene_count"] = phenotype_scorer.gene_count(parent["hpo_id"])
    return results


@router.get("/hpo/{hpo_id:path}")
def hpo_get(hpo_id: str):
    t = hpo_ontology.get(hpo_id)
    if t is None:
        raise HTTPException(404, f"unknown HPO term: {hpo_id}")
    return t.to_dict()


@router.get("/panels")
def panels_list():
    return phenotype_scorer.list_panels()


@router.post("/samples/{sample_id}/phenotype")
def update_phenotype(sample_id: str, payload: dict):
    sub = sample_layout.state_dir(sample_id)
    if not sub.is_dir():
        raise HTTPException(404, f"sample not found: {sample_id}")

    hpo_in = payload.get("hpo") or []
    panels_in = payload.get("panels") or []
    # Caller can target a specific version; otherwise we land on the
    # currently-active version, creating 'default' on the fly for
    # un-migrated samples that have nothing yet.
    target_version = payload.get("version")

    from ..services import analyses_store
    if target_version:
        analyses_store.validate_name(target_version)
    else:
        target_version = analyses_store.active_version(sample_id) or "default"

    meta_path = sample_layout.state_file(sample_id, "sample_metadata.json")
    if not meta_path.is_file():
        raise HTTPException(404, f"registered sample not found: {sample_id}")
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        raise HTTPException(500, "sample metadata could not be read")
    if not isinstance(meta, dict):
        raise HTTPException(500, "sample metadata is malformed")
    mrn = str(meta.get("mrn") or "").strip()
    try:
        patient_phenotype_store.check_token("MRN", mrn, required=True)
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc

    # 1. Persist hpo/panels into analyses/{version}/analysis.json.
    analyses_store.write_version(
        sample_id, target_version,
        hpo=hpo_in, panels=panels_in,
        note=payload.get("note", ""),
    )

    # Update sample_metadata.json's active_analysis pointer + clean up
    # any legacy `hpo` / `selected_panels` left over from before the
    # migration so the loader can stop reading them on next load.
    meta.pop("hpo", None)
    meta.pop("patient_phenotype", None)
    meta.pop("selected_panels", None)
    meta["active_analysis"] = target_version
    meta["phenotype_updated_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")

    # 2. Only the default analysis is the patient's reusable phenotype.
    # Follow-up analysis versions are sample-owned experiments/combinations
    # and must not overwrite the MRN-level snapshot used by future LIS IDs.
    # Empty default input intentionally becomes a header-only file, preventing
    # stale LIS-specific legacy files from resurfacing.
    patient_snapshot = None
    if target_version == "default":
        try:
            patient_snapshot = patient_phenotype_store.save(
                mrn=mrn,
                code=str(meta.get("lis_id") or sample_id),
                hpo=hpo_in,
                panels=panels_in,
            )
        except ValueError as exc:
            raise HTTPException(400, str(exc)) from exc
        except OSError as exc:
            raise HTTPException(500, f"病人 phenotype 儲存失敗：{exc}") from exc

    # 3. Compute pheno_score (analyses_store.write_version already
    # wrote pheno_score.tsv as a side effect; recompute here only to
    # produce response stats). SNV loads and gene search apply
    # in-panel state dynamically from pheno_score.tsv; do not rewrite
    # the large raw TSV here.
    scores = phenotype_scorer.compute_pheno_score(hpo_in, panels_in)
    in_panel_genes = {g for g, s in scores.items() if s > 0}

    # 5. Stats for UI
    top10 = sorted(scores.items(), key=lambda kv: -kv[1])[:10]
    sample_loader.update_case_table_row(sample_id)
    return {
        "sample_id":         sample_id,
        "n_hpo":             len(hpo_in),
        "n_panels":          len(panels_in),
        "n_genes_scored":    len(scores),
        "n_in_panel_genes":  len(in_panel_genes),
        "top_score":         max(scores.values(), default=0.0),
        "top10":             [{"gene": g, "score": round(s, 2)} for g, s in top10],
        "updated_at":        meta["phenotype_updated_at"],
        "patient_phenotype": patient_snapshot,
        "patient_phenotype_synced": patient_snapshot is not None,
    }
