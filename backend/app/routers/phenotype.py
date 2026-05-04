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
from pathlib import Path

from fastapi import APIRouter, HTTPException, Query

from ..config import TERTIARY_OUTPUT_ROOT
from ..services import hpo_ontology, phenotype_scorer

router = APIRouter(prefix="/api", tags=["phenotype"])


@router.get("/hpo/search")
def hpo_search(q: str = Query(""), limit: int = Query(20, ge=1, le=100)):
    return hpo_ontology.search(q, limit=limit)


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
    sub = TERTIARY_OUTPUT_ROOT / sample_id
    if not sub.is_dir():
        raise HTTPException(404, f"sample not found: {sample_id}")

    hpo_in = payload.get("hpo") or []
    panels_in = payload.get("panels") or []

    # 1. Persist into sample_metadata.json
    meta_path = sub / "sample_metadata.json"
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            meta = {}
    else:
        meta = {}
    if not isinstance(meta, dict):
        meta = {}
    meta["hpo"] = hpo_in
    meta["selected_panels"] = panels_in
    meta["phenotype_updated_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")

    # 2. Compute pheno_score
    scores = phenotype_scorer.compute_pheno_score(hpo_in, panels_in)

    # 3. Persist gene → score sidecar
    phenotype_scorer.write_pheno_table(sample_id, scores)

    # 4. Rewrite IN_PANEL column (in_panel iff score > 0)
    in_panel_genes = {g for g, s in scores.items() if s > 0}
    n_updated = phenotype_scorer.update_in_panel_column(sample_id, in_panel_genes)

    # 5. Stats for UI
    top10 = sorted(scores.items(), key=lambda kv: -kv[1])[:10]
    return {
        "sample_id":         sample_id,
        "n_hpo":             len(hpo_in),
        "n_panels":          len(panels_in),
        "n_genes_scored":    len(scores),
        "n_in_panel_genes":  len(in_panel_genes),
        "top_score":         max(scores.values(), default=0.0),
        "top10":             [{"gene": g, "score": round(s, 2)} for g, s in top10],
        "tsv_rows_updated":  n_updated,
        "updated_at":        meta["phenotype_updated_at"],
    }
