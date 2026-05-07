from urllib.parse import quote

from fastapi import APIRouter, Depends, Form, HTTPException, UploadFile
from fastapi.responses import Response

from ..auth import current_user
from ..services import docx_export, patient_store, report_store, sample_loader

router = APIRouter(prefix="/api", tags=["samples"], dependencies=[Depends(current_user)])


@router.get("/samples/{sample_id}/report.docx")
def get_report_docx(sample_id: str):
    try:
        blob = docx_export.build_diagnosis_docx(sample_id)
    except FileNotFoundError as e:
        raise HTTPException(404, str(e))
    fname = quote(f"{sample_id}_diagnosis.docx")
    return Response(
        content=blob,
        media_type="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
        headers={"Content-Disposition": f"attachment; filename*=UTF-8''{fname}"},
    )


@router.get("/samples")
def list_samples():
    return sample_loader.list_index()


@router.get("/samples/unregistered")
def list_unregistered_samples():
    """Pipeline-dropped directories not yet attached to reviewer info.

    The 載入新個案 modal calls this to populate the LIS_ID dropdown so
    reviewers don't have to retype an ID that already lives on disk.
    """
    return sample_loader.list_unregistered()


@router.post("/samples")
async def register_sample(
    lis_id:         str = Form(...),
    name:           str = Form(...),
    mrn:            str = Form(...),
    test_type:      str = Form("WES"),
    genome_build:   str = Form("hg38"),
    category:       str = Form(""),
    vcf_path:       str = Form(""),
    phenotype_path: str = Form(""),
    phenotype_file: UploadFile | None = None,
):
    """Attach reviewer-side info to a pipeline-produced directory.

    The TSV must already live at
        tertiary_output/{lis_id}/snv_indel.annotated.tsv
    (the pipeline puts it there). This endpoint only writes
    sample_metadata.json + analyses/default/analysis.json on top.
    Refuses with 404 if the directory is missing, 409 if it already has
    a sample_metadata.json.

    Phenotype is optional. Accepts either an uploaded blob
    (phenotype_file) or a server-side path (phenotype_path).
    """
    # Phenotype: pick whichever form the form supplied; both empty is OK.
    phenotype_text = ""
    if phenotype_file is not None and phenotype_file.filename:
        phenotype_text = (await phenotype_file.read()).decode("utf-8", errors="replace")
    elif phenotype_path:
        from pathlib import Path
        p = Path(phenotype_path)
        if not p.is_file():
            raise HTTPException(400, f"phenotype_path not found: {phenotype_path}")
        phenotype_text = p.read_text(encoding="utf-8")

    try:
        meta = patient_store.register(
            lis_id=lis_id, name=name, mrn=mrn,
            test_type=test_type, genome_build=genome_build,
            category=category, vcf_path=vcf_path,
            phenotype_text=phenotype_text,
        )
    except FileExistsError as e:
        raise HTTPException(409, str(e))
    except FileNotFoundError as e:
        raise HTTPException(404, str(e))
    except ValueError as e:
        raise HTTPException(400, str(e))
    return {"sample_id": lis_id, "meta": meta}


@router.get("/samples/{sample_id}")
def get_sample(sample_id: str, version: str | None = None):
    payload = sample_loader.load_sample(sample_id, version=version)
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
