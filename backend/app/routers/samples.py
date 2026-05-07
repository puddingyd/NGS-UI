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


@router.post("/samples")
async def create_sample(
    lis_id:        str = Form(...),
    name:          str = Form(...),
    mrn:           str = Form(...),
    test_type:     str = Form("WES"),
    genome_build:  str = Form("hg38"),
    category:      str = Form(""),
    vcf_path:      str = Form(""),
    tsv_path:      str = Form(""),
    phenotype_path:str = Form(""),
    tsv_file:      UploadFile | None = None,
    phenotype_file:UploadFile | None = None,
):
    """Create a new patient directory.

    Accepts both browser multipart uploads (tsv_file/phenotype_file)
    AND server-side paths (tsv_path/phenotype_path). For each pair you
    can pick either form per file independently. The function refuses
    to create when:

      * lis_id is already a sample dir (UI offers an 'open existing'
        button rather than overwriting),
      * neither tsv_file nor tsv_path is supplied,
      * the supplied tsv_path doesn't resolve to a file.

    Phenotype is optional; with neither file nor path the default
    analysis lands with empty hpo/panels.
    """
    if patient_store.sample_exists(lis_id):
        raise HTTPException(409, f"sample already exists: {lis_id}")

    # TSV: prefer uploaded blob, fall back to path pointer.
    tsv_bytes = None
    tsv_src   = None
    if tsv_file is not None and tsv_file.filename:
        tsv_bytes = await tsv_file.read()
    elif tsv_path:
        from pathlib import Path
        tsv_src = Path(tsv_path)
    else:
        raise HTTPException(400, "either tsv_file or tsv_path is required")

    # Phenotype: same dual-input shape, but optional.
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
        meta = patient_store.create_new(
            lis_id=lis_id, name=name, mrn=mrn,
            test_type=test_type, genome_build=genome_build,
            category=category, vcf_path=vcf_path,
            tsv_src=tsv_src, tsv_bytes=tsv_bytes,
            phenotype_text=phenotype_text,
        )
    except FileExistsError as e:
        raise HTTPException(409, str(e))
    except (FileNotFoundError, ValueError) as e:
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
