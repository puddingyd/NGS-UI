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
from ..services import phenotype_scorer

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


@router.post("/save")
def save_phenotype_file(payload: dict):
    """Write the tool's generated phenotype.txt into PHENOTYPE_DIR.

    Body: {mrn (required), code (optional LIS_ID), content}.
    """
    mrn  = _check_token("MRN",     (payload or {}).get("mrn", ""),  required=True)
    code = _check_token("LIS_ID",  (payload or {}).get("code", ""), required=False)
    content = (payload or {}).get("content", "")
    if not isinstance(content, str) or not content.strip():
        raise HTTPException(400, "內容為空")
    if len(content.encode("utf-8")) > _MAX_CONTENT_BYTES:
        raise HTTPException(400, "內容過大")

    fname = f"{code}_{mrn}_phenotype.txt" if code else f"{mrn}_phenotype.txt"
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
    # Best-effort parse of {code}_{mrn}_phenotype or {mrn}_phenotype.
    stem = candidate.stem
    if stem.endswith("_phenotype"):
        core = stem[:-len("_phenotype")]
    else:
        core = stem
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
