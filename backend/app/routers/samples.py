from urllib.parse import quote

from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile
from fastapi.responses import Response

from ..auth import current_user
from ..config import PHENOTYPE_DIR
from ..services import (
    docx_export,
    patient_list_store,
    patient_store,
    report_store,
    sample_loader,
)

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


@router.post("/patient_list")
async def upload_patient_list(file: UploadFile = File(...)):
    """Ingest a 未完成報告清單 xlsx → archive it + merge into roster.json.

    The roster maps LIS_ID → {mrn, name, test_type, department}; the
    載入新個案 modal reads it to auto-fill those fields when the
    reviewer picks a pipeline-dropped TSV.
    """
    content = await file.read()
    if not content:
        raise HTTPException(400, "空白檔案")
    try:
        return patient_list_store.ingest_xlsx(content, file.filename or "upload.xlsx")
    except ValueError as e:
        raise HTTPException(400, str(e))


@router.get("/patient_list")
def get_patient_list():
    """Current merged roster ({lis_id: {...}}). For debugging / UI display."""
    return patient_list_store.load_roster()


@router.get("/patient_list/uploads")
def get_patient_list_uploads():
    """History of every successful 個案清單 ingest, latest first.

    Each entry: {uploaded_at, original_filename, archive_name, parsed,
    added, updated, total_after}. Powers the 上傳記錄 modal.
    """
    return patient_list_store.list_uploads()


@router.post("/samples")
def register_sample(
    lis_id:        str = Form(...),
    name:          str = Form(...),
    mrn:           str = Form(...),
    sex:           str = Form(""),
    test_type:     str = Form("WES"),
    genome_build:  str = Form("hg38"),
    category:      str = Form(""),
    hpo_json:      str = Form(""),
    panels_json:   str = Form(""),
    run_analysis:  str = Form(""),
):
    """Attach reviewer-side info to a pipeline-produced directory.

    The TSV must already live at
        tertiary_output/{lis_id}/snv_indel.annotated.tsv
    (the pipeline puts it there). register() also generates the
    minimal Exomiser/LIRICAL-input VCF beside it
    ({lis_id}.from_tsv.vcf.gz) so the operator never has to point at a
    VCF manually.

    Phenotype is auto-loaded from
        NGS_UI/patient_phenotype/{lis_id}_{mrn}_phenotype.txt
    when the file exists. If it doesn't, the sample registers with
    empty hpo/panels and the response includes phenotype_loaded=false
    so the UI can hint about it.
    """
    # Frontend-edited chips arrive as JSON strings; an empty string
    # means "no override → fall back to file/EMR" (handled inside
    # register()). Bad JSON falls back too rather than 4xxing.
    import json as _json
    hpo_payload = panels_payload = None
    if hpo_json or panels_json:
        try:
            hpo_payload    = _json.loads(hpo_json)    if hpo_json    else []
            panels_payload = _json.loads(panels_json) if panels_json else []
        except _json.JSONDecodeError:
            hpo_payload = panels_payload = None

    pheno_path = PHENOTYPE_DIR / f"{lis_id}_{mrn}_phenotype.txt"
    phenotype_text = ""
    phenotype_loaded = False
    if pheno_path.is_file():
        phenotype_text = pheno_path.read_text(encoding="utf-8")
        phenotype_loaded = True

    try:
        meta = patient_store.register(
            lis_id=lis_id, name=name, mrn=mrn, sex=sex,
            test_type=test_type, genome_build=genome_build,
            category=category,
            phenotype_text=phenotype_text,
            hpo=hpo_payload, panels=panels_payload,
        )
    except FileExistsError as e:
        raise HTTPException(409, str(e))
    except FileNotFoundError as e:
        raise HTTPException(404, str(e))
    except ValueError as e:
        raise HTTPException(400, str(e))
    # If the reviewer asked us to run analysis on register (or if any
    # phenotype data ended up populated), enqueue exomiser/lirical so
    # the cards have scores by the time they finish reading the
    # patient's basic info. Failure to enqueue is non-fatal — the
    # sample is registered, reviewer can hit ▶ 開始分析 manually.
    job_id = None
    should_run = (str(run_analysis).lower() in ("1", "true", "yes", "on")
                  or bool(hpo_payload) or bool(panels_payload))
    if should_run:
        try:
            from . import jobs as _jobs
            rec = _jobs._enqueue(lis_id, "exomiser_lirical", version="default")
            job_id = rec.get("job_id")
        except Exception:
            job_id = None
    return {
        "sample_id": lis_id,
        "meta": meta,
        "phenotype_loaded": phenotype_loaded,
        "phenotype_path":   str(pheno_path),
        "job_id":           job_id,
    }


@router.get("/samples/{sample_id}")
def get_sample(sample_id: str, version: str | None = None):
    # include_aux=False → core payload only (meta + SNV/Indel + report
    # state); the frontend pulls CNV/SV and Mito separately so the
    # SNV/Indel view shows up first (staged loading).
    payload = sample_loader.load_sample(sample_id, version=version, include_aux=False)
    if payload is None:
        raise HTTPException(404, f"sample not found: {sample_id}")
    return payload


@router.get("/samples/{sample_id}/cnv-sv")
def get_sample_cnv_sv(sample_id: str, version: str | None = None):
    """CNV/SV side-channels for the staged loader."""
    payload = sample_loader.load_sample_cnv_sv(sample_id, version=version)
    if payload is None:
        raise HTTPException(404, f"sample not found: {sample_id}")
    return payload


@router.get("/samples/{sample_id}/mito")
def get_sample_mito(sample_id: str, version: str | None = None):
    """Mitochondria side-channel for the staged loader."""
    payload = sample_loader.load_sample_mito(sample_id, version=version)
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
    EDITABLE = {"name", "mrn", "lis_id", "sex", "test_type", "category",
                "genome_build", "tags", "run_date"}
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


# Canonical category list — drives both the load-new-case modal
# dropdown AND the editable Category select on the sample card so the
# values stay in sync. Order matches the reviewer-requested ordering.
_CATEGORY_OPTIONS = [
    "Neurology", "Endocrinology", "MCA", "Nephrology", "GI", "Metabolism",
    "AIR", "Hematology", "Oncology", "Ophthalmology", "Musculoskeletal",
    "Dermatology", "CV", "ENT", "Asymptomatic",
]


@router.get("/options")
def get_options():
    """Category list + (optional) tag suggestions.

    Categories are hard-coded server-side so adding a new one is a
    one-line edit + restart, not a config file the operator has to
    remember to update. Tag suggestions still come from _options.json
    when present so reviewers can keep iterating on the tag vocabulary
    without a deploy.
    """
    from ..config import TERTIARY_OUTPUT_ROOT
    import json as _json
    payload = {"category_options": list(_CATEGORY_OPTIONS), "tag_suggestions": []}
    p = TERTIARY_OUTPUT_ROOT / "_options.json"
    if p.exists():
        try:
            extra = _json.loads(p.read_text(encoding="utf-8"))
            if isinstance(extra, dict) and isinstance(extra.get("tag_suggestions"), list):
                payload["tag_suggestions"] = extra["tag_suggestions"]
        except _json.JSONDecodeError:
            pass
    return payload
