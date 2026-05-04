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


@router.put("/samples/{sample_id}/metadata")
def put_sample_metadata(sample_id: str, payload: dict):
    """Edit a small whitelist of sample_metadata.json fields from the UI.

    Only the operator-facing identifiers + sequencing/build live here.
    HPO + selected_panels go via /api/samples/{id}/phenotype.
    """
    import json as _json
    from datetime import datetime, timezone

    from ..config import TERTIARY_OUTPUT_ROOT
    sub = TERTIARY_OUTPUT_ROOT / sample_id
    if not sub.is_dir():
        raise HTTPException(404, f"sample not found: {sample_id}")
    meta_path = sub / "sample_metadata.json"
    meta = {}
    if meta_path.exists():
        try:
            meta = _json.loads(meta_path.read_text(encoding="utf-8"))
        except _json.JSONDecodeError:
            meta = {}
    if not isinstance(meta, dict):
        meta = {}
    EDITABLE = {"name", "mrn", "lis_id", "test_type", "category",
                "genome_build", "vcf_path", "tags", "run_date"}
    for k, v in (payload or {}).items():
        if k in EDITABLE:
            meta[k] = v
    meta["metadata_updated_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    meta_path.write_text(_json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    return meta


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
