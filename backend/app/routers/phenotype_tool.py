"""Public endpoints for the standalone 輸入臨床表徵 tool (served at
/phenotype/).

This router has NO auth dependency on purpose — the tool runs on the
hospital intranet and clinicians use it without logging in (only the
analysis app gates behind auth). Writes are constrained: MRN / LIS_ID
must be short alphanumeric tokens, the filename is built from those
validated parts (so no path traversal), and the body size is capped.

Output txt lands in PHENOTYPE_DIR as either
    {lis_id}_{mrn}_phenotype.txt   (when a LIS_ID was given)
    {mrn}_phenotype.txt            (MRN only)
so sample_loader.list_unregistered picks it up automatically — it
looks for {lis_id}_*_phenotype.txt first, then {mrn}_phenotype.txt.
"""
from __future__ import annotations

import re

from fastapi import APIRouter, HTTPException, Query

from ..config import PHENOTYPE_DIR
from ..services import hpo_ontology, panel_deadzone, phenotype_scorer

router = APIRouter(prefix="/api/phenotype-tool", tags=["phenotype-tool"])

_TOKEN_RE = re.compile(r"^[A-Za-z0-9_-]{1,32}$")
_MAX_CONTENT_BYTES = 64 * 1024     # phenotype txt files are tiny


def _check_token(name: str, value: str, *, required: bool) -> str:
    v = (value or "").strip()
    if not v:
        if required:
            raise HTTPException(400, f"{name} 為必填")
        return ""
    if not _TOKEN_RE.match(v):
        raise HTTPException(400, f"{name} 只能是英數 / - / _（最多 32 字）")
    return v


@router.get("/panels")
def list_panels_public():
    """Gene-panel list for the tool's panel autocomplete (read-only)."""
    return phenotype_scorer.list_panels()


@router.get("/fixed-panels")
def list_fixed_panels():
    """Index of the WES-I / WES-II / WGS fixed panels for the new
    Gene-panels UI tabs. Built by scripts/import_fixed_panels.py into
    the repo's phenotype_data/fixed_panels/index.json; we just stream it.
    Returns {"series": []} on a clean dev box (no panels imported).
    """
    import json
    from ..config import FIXED_PANELS_DIR
    idx = FIXED_PANELS_DIR / "index.json"
    if not idx.is_file():
        return {"series": []}
    try:
        return json.loads(idx.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {"series": []}


@router.get("/gene-list")
def gene_list(kind: str = Query("", pattern="^(|hpo|panel)$"), key: str = Query(..., min_length=1)):
    """Return canonical gene symbols for one HPO term or panel key."""
    data = phenotype_scorer.genes_for_key(key, kind)
    if not data.get("gene_count"):
        raise HTTPException(404, "找不到這個 HPO term / panel 的 gene list")
    return data


@router.get("/gene-memberships")
def gene_memberships(gene: str = Query(..., min_length=1, max_length=64)):
    """Find HPO terms and panels that contain one canonical gene symbol."""
    data = phenotype_scorer.memberships_for_gene(gene)
    for item in data.get("hpo", []):
        term = hpo_ontology.get(item.get("id", ""))
        if term:
            item["name"] = term.name
        else:
            item["name"] = ""
    return data


_GENE_SPLIT_RE = re.compile(r"[,\s]+")
_MAX_PANEL_GENES = 5000
_MAX_DEAD_ZONE_GENES = 10000


def _sort_dead_zone_entries(entries: list[dict]) -> list[dict]:
    def bucket(pct: float) -> int:
        if pct >= 70:
            return 0
        if pct >= 50:
            return 1
        if pct >= 30:
            return 2
        return 3

    return sorted(
        entries,
        key=lambda e: (
            bucket(float(e.get("cds_dead_pct") or 0)),
            -float(e.get("cds_dead_pct") or 0),
            str(e.get("gene") or ""),
        ),
    )


@router.post("/dead-zone")
def dead_zone_for_gene_list(payload: dict):
    """Return WES/WGS cohort dead-zone rows for a supplied gene list.

    The phenotype drawer already has the canonical gene list, so this endpoint
    only runs when the reviewer explicitly asks for dead-zone rows. The large
    WES/WGS tables stay behind panel_deadzone's process LRU cache.
    """
    test_type = str((payload or {}).get("test_type") or "").upper()
    if test_type not in {"WES", "WGS"}:
        raise HTTPException(400, "test_type 必須是 WES 或 WGS")
    raw_genes = (payload or {}).get("genes") or []
    if isinstance(raw_genes, str):
        genes = [g for g in _GENE_SPLIT_RE.split(raw_genes) if g]
    elif isinstance(raw_genes, (list, tuple)):
        genes = [str(g or "").strip() for g in raw_genes if str(g or "").strip()]
    else:
        genes = []
    if len(genes) > _MAX_DEAD_ZONE_GENES:
        raise HTTPException(400, f"基因數過多（上限 {_MAX_DEAD_ZONE_GENES}）")
    hits = panel_deadzone.dead_zone_for_genes(test_type, genes)
    entries = _sort_dead_zone_entries([dict(v) for v in hits.values()])
    return {
        "test_type": test_type,
        "threshold": panel_deadzone.dead_zone_threshold(test_type),
        "gene_count": len(set(genes)),
        "entries": entries,
    }


@router.post("/custom-panel")
def create_custom_panel(payload: dict):
    """Create a reusable gene panel from a user-supplied gene list.

    Body: {name, genes}. `genes` may be a comma/whitespace-separated
    string or a list of strings. The name is sanitised server-side
    (non [A-Za-z0-9_-] runs → '_'); a collision with an existing panel
    is a 409. On success the panel file is written and the in-memory
    tables updated, so analysis can use it immediately. Returns
    {"name": <sanitised>, "n_genes": int}.
    """
    raw_genes = (payload or {}).get("genes", "")
    if isinstance(raw_genes, str):
        genes = [g for g in _GENE_SPLIT_RE.split(raw_genes) if g]
    elif isinstance(raw_genes, (list, tuple)):
        genes = []
        for chunk in raw_genes:
            genes.extend(g for g in _GENE_SPLIT_RE.split(str(chunk or "")) if g)
    else:
        genes = []
    if len(genes) > _MAX_PANEL_GENES:
        raise HTTPException(400, f"基因數過多（上限 {_MAX_PANEL_GENES}）")
    try:
        return phenotype_scorer.register_custom_panel(
            (payload or {}).get("name", ""),
            genes,
            source=(payload or {}).get("source", ""),
        )
    except ValueError as e:
        msg = str(e)
        # "已存在" → 409 Conflict; everything else is a bad request.
        raise HTTPException(409 if "已存在" in msg else 400, msg)


@router.post("/save")
def save_phenotype_file(payload: dict):
    """Write the tool's generated phenotype.txt into PHENOTYPE_DIR.

    Body: {mrn, code (LIS_ID), content}. At least one of mrn / code
    is required. Filename: {code}_{mrn}_ when both given, {code}_
    when only LIS_ID, {mrn}_ when only MRN.
    """
    mrn  = _check_token("MRN",    (payload or {}).get("mrn", ""),  required=False)
    code = _check_token("LIS_ID", (payload or {}).get("code", ""), required=False)
    if not mrn and not code:
        raise HTTPException(400, "請至少提供 MRN 或 LIS_ID")
    content = (payload or {}).get("content", "")
    if not isinstance(content, str) or not content.strip():
        raise HTTPException(400, "內容為空")
    if len(content.encode("utf-8")) > _MAX_CONTENT_BYTES:
        raise HTTPException(400, "內容過大")

    if code and mrn:
        fname = f"{code}_{mrn}_phenotype.txt"
    elif code:
        fname = f"{code}_phenotype.txt"
    else:
        fname = f"{mrn}_phenotype.txt"
    PHENOTYPE_DIR.mkdir(parents=True, exist_ok=True)
    out = PHENOTYPE_DIR / fname
    # Defence in depth: the name is built from validated tokens, but
    # resolve and re-check it sits directly under PHENOTYPE_DIR.
    if out.resolve().parent != PHENOTYPE_DIR.resolve():
        raise HTTPException(400, "檔名不合法")
    out.write_text(content if content.endswith("\n") else content + "\n", encoding="utf-8")
    return {"path": str(out), "filename": fname, "mrn": mrn, "code": code}


@router.get("/load")
def load_phenotype_file(
    mrn:  str | None = Query(None),
    code: str | None = Query(None),
):
    """Return an existing phenotype.txt for editing.

    Lookup order: {code}_*_phenotype.txt (if code given) → exact
    {mrn}_phenotype.txt → *_{mrn}_phenotype.txt. 404 when nothing
    matches.
    """
    code = _check_token("LIS_ID", code or "", required=False)
    mrn  = _check_token("MRN",    mrn  or "", required=False)
    if not code and not mrn:
        raise HTTPException(400, "請提供 LIS_ID 或 MRN")
    if not PHENOTYPE_DIR.is_dir():
        raise HTTPException(404, "尚無任何 phenotype 檔")

    candidate = None
    if code:
        # exact "{code}_{mrn}" if both given, then code-only, then any MRN
        if mrn:
            p = PHENOTYPE_DIR / f"{code}_{mrn}_phenotype.txt"
            if p.is_file():
                candidate = p
        if candidate is None:
            p = PHENOTYPE_DIR / f"{code}_phenotype.txt"
            if p.is_file():
                candidate = p
        if candidate is None:
            matches = sorted(PHENOTYPE_DIR.glob(f"{code}_*_phenotype.txt"))
            if matches:
                candidate = matches[0]
    if candidate is None and mrn:
        exact = PHENOTYPE_DIR / f"{mrn}_phenotype.txt"
        if exact.is_file():
            candidate = exact
        else:
            matches = sorted(PHENOTYPE_DIR.glob(f"*_{mrn}_phenotype.txt"))
            if matches:
                candidate = matches[0]
    if candidate is None:
        raise HTTPException(404, "找不到對應的 phenotype 檔")

    try:
        text = candidate.read_text(encoding="utf-8")
    except OSError:
        raise HTTPException(404, "讀取失敗")
    # Best-effort recover {code, mrn} from the filename. Three shapes:
    #   {code}_{mrn}_phenotype.txt   {code}_phenotype.txt   {mrn}_phenotype.txt
    stem = candidate.stem
    core = stem[:-len("_phenotype")] if stem.endswith("_phenotype") else stem
    if code and core == code:
        parsed_code, parsed_mrn = code, ""
    elif mrn and core == mrn:
        parsed_code, parsed_mrn = "", mrn
    elif code and core.startswith(code + "_"):
        parsed_code, parsed_mrn = code, core[len(code) + 1:]
    else:
        parts = core.split("_")
        if len(parts) >= 2:
            parsed_code, parsed_mrn = parts[0], parts[1]
        else:
            parsed_code, parsed_mrn = "", core
    return {
        "filename": candidate.name,
        "content":  text,
        "code":     parsed_code,
        "mrn":      parsed_mrn,
    }
