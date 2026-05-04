from fastapi import APIRouter, HTTPException

from ..services import report_store, sample_loader

router = APIRouter(prefix="/api", tags=["samples"])


@router.get("/healthz")
def healthz():
    return {"ok": True}


@router.get("/samples")
def list_samples():
    return sample_loader.list_index()


@router.get("/samples/{sample_id}")
def get_sample(sample_id: str):
    payload = sample_loader.load_sample(sample_id)
    if payload is None:
        raise HTTPException(404, f"sample not found: {sample_id}")
    return payload


@router.get("/samples/{sample_id}/report")
def get_report(sample_id: str):
    return report_store.load(sample_id)


@router.put("/samples/{sample_id}/report")
def put_report(sample_id: str, payload: dict):
    return report_store.save(sample_id, payload)


@router.get("/options")
def get_options():
    """Category / tag dropdown suggestions; reads _options.json if present."""
    from ..config import TERTIARY_OUTPUT_ROOT
    import json as _json
    p = TERTIARY_OUTPUT_ROOT / "_options.json"
    if not p.exists():
        return {"category_options": [], "tag_suggestions": []}
    try:
        return _json.loads(p.read_text(encoding="utf-8"))
    except _json.JSONDecodeError:
        return {"category_options": [], "tag_suggestions": []}
