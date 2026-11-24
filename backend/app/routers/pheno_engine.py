"""Phenotype prediction from HPO terms (plan §4.3–4.8, phase 2 engine).

  POST /api/pheno-engine/predict   {"hpo": [{"phenotype": "HP:…", "weight": 1}, …], "top": 10}

No login (like the /phenotype/ tool): the request holds only HPO IDs and the
response only knowledge-base rankings; nothing is stored.
"""
from __future__ import annotations

import re

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

router = APIRouter(prefix="/api/pheno-engine", tags=["pheno-engine"])

MAX_TERMS = 200
_HP = re.compile(r"^HP:\d{7}$")


class PredictRequest(BaseModel):
    hpo: list = Field(default_factory=list)
    top: int = 10


def _items(raw: list) -> list[tuple[str, float]]:
    out = []
    for entry in raw[:MAX_TERMS]:
        if isinstance(entry, dict):
            hid = str(entry.get("phenotype") or entry.get("hpo_id") or "").strip()
            weight = entry.get("weight", 1)
        else:
            hid, weight = str(entry).strip(), 1
        hid = hid.upper()
        if _HP.match(hid):
            out.append((hid, weight))
    return out


def predict(items: list, top: int = 10, engine=None) -> dict:
    try:
        from ..services.pheno_engine import engine as E
    except ImportError as exc:                       # numpy missing on this server
        raise HTTPException(503, f"表型引擎無法載入（{exc}）") from exc
    engine = engine or E.load_engine()
    if engine is None:
        raise HTTPException(503, "尚未建立表型知識庫：請先在首頁「資料庫版本」按「更新 HPO」")
    top = max(1, min(int(top or 10), 50))
    result = engine.score(items)
    kb = engine.kb
    diseases = []
    for d in result.top_diseases(top):
        di = engine.did[d["disease_id"]]
        d["n_annotations"] = int(kb.dir_ptr[di + 1] - kb.dir_ptr[di])
        d["explain"] = [{k: r[k] for k in ("hpo_id", "hpo_name", "weight", "annotation", "annotation_name",
                                            "relation", "mica_name", "s", "mu", "net")}
                        for r in result.explain(di)]
        for r in d["explain"]:
            r["relation_zh"] = E.RELATION_ZH[r["relation"]]
        diseases.append(d)
    genes = result.top_genes(top)
    return {
        "snapshot": engine.snapshot,
        "hp_version": engine.manifest.get("hp_version", ""),
        "n_diseases": kb.N,
        "n_genes": len(engine.genes),
        "terms": [{k: t[k] for k in ("hpo_id", "name", "input_id", "how", "weight")} for t in result.pre["terms"]],
        "removed": [dict(r, reason_zh=E.REMOVED_ZH.get(r["reason"], r["reason"])) for r in result.pre["removed"]],
        "genes": genes,
        "diseases": diseases,
        "lead_nats": result.lead_nats(),
        "calibration": "表型點數 P 為文獻病例估計的暫定值",
    }


@router.post("/predict")
def post_predict(req: PredictRequest):
    items = _items(req.hpo or [])
    if not items:
        raise HTTPException(400, "沒有可用的 HPO term（格式 HP:0000000）")
    return predict(items, req.top)
