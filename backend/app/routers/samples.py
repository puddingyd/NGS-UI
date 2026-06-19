import time
from urllib.parse import quote

from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile
from fastapi.responses import Response

from ..auth import current_user
from ..config import PHENOTYPE_DIR, TERTIARY_OUTPUT_ROOT
from ..services import (
    clinical_presentation_store,
    docx_export,
    patient_list_store,
    patient_store,
    report_store,
    sample_loader,
)

router = APIRouter(prefix="/api", tags=["samples"], dependencies=[Depends(current_user)])


def _log_perf(event: str, started: float, **fields) -> None:
    parts = [f"[perf] {event}", f"elapsed={time.perf_counter() - started:.3f}s"]
    parts.extend(f"{key}={value}" for key, value in fields.items())
    print(" ".join(parts), flush=True)


@router.get("/samples/{sample_id}/report.docx")
def get_report_docx(sample_id: str, gene_list_mode: str = "grouped"):
    """Render the diagnostic report DOCX, save a copy to
    NGS_UI/report/, and stream it back as a download.

    `gene_list_mode` controls §五.4 「本次檢測基因包括」:
      grouped (default) → one paragraph per HPO term / panel
      merged            → single deduped flat list
    """
    if gene_list_mode not in ("grouped", "merged"):
        gene_list_mode = "grouped"
    try:
        blob = docx_export.build_diagnosis_docx(
            sample_id, gene_list_mode=gene_list_mode,
        )
    except FileNotFoundError as e:
        raise HTTPException(404, str(e))

    # Archive: NGS_UI/report/{SID}_diagnosis_{YYYYMMDDTHHMMSSZ}.docx —
    # reviewers can audit/redownload past renders without re-clicking
    # 匯出. The latest copy is also dropped at {SID}_diagnosis.docx for
    # easy linking from external tools.
    from datetime import datetime, timezone
    from ..config import REPORT_OUTPUT_DIR
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    archive = REPORT_OUTPUT_DIR / f"{sample_id}_diagnosis_{ts}.docx"
    latest  = REPORT_OUTPUT_DIR / f"{sample_id}_diagnosis.docx"
    try:
        archive.write_bytes(blob)
        latest.write_bytes(blob)
    except OSError:
        # Disk write failure shouldn't prevent the download; the
        # browser still gets the docx. Reviewer just loses the archive.
        pass

    fname = quote(f"{sample_id}_diagnosis.docx")
    return Response(
        content=blob,
        media_type="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
        headers={"Content-Disposition": f"attachment; filename*=UTF-8''{fname}"},
    )


@router.get("/samples")
def list_samples():
    return sample_loader.list_index()


@router.get("/samples/case-summary")
def list_case_summaries():
    return sample_loader.list_case_summaries()


@router.get("/samples/unregistered")
def list_unregistered_samples():
    """Pipeline-dropped directories not yet attached to reviewer info.

    The 載入新個案 modal calls this to populate the LIS_ID dropdown so
    reviewers don't have to retype an ID that already lives on disk.
    """
    return sample_loader.list_unregistered()


@router.delete("/samples/{sample_id}")
def delete_sample(sample_id: str, delete_pipeline_output: bool = False):
    """Remove one UI sample directory and optionally its pipeline output."""
    try:
        return patient_store.delete(
            sample_id, delete_pipeline_output=delete_pipeline_output,
        )
    except ValueError as e:
        raise HTTPException(400, str(e))
    except FileNotFoundError as e:
        raise HTTPException(404, f"sample not found: {sample_id}")
    except RuntimeError as e:
        raise HTTPException(409, str(e))
    except OSError as e:
        raise HTTPException(500, f"刪除失敗：{e}")


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


@router.get("/patient_list/options")
def get_patient_list_options():
    """Distinct 科別 / 開單醫師 values from the roster — fills the
    editable datalists on the sample card so reviewers can pick or
    type a new value.
    """
    return patient_list_store.list_options()


@router.post("/samples")
def register_sample(
    lis_id:        str = Form(...),
    name:          str = Form(...),
    mrn:           str = Form(...),
    sex:              str = Form(""),
    test_type:        str = Form("WES"),
    genome_build:     str = Form("hg38"),
    category:         str = Form(""),
    department:       str = Form(""),
    physician:        str = Form(""),
    sign_received_at: str = Form(""),
    hpo_json:         str = Form(""),
    panels_json:      str = Form(""),
    run_analysis:     str = Form(""),
):
    """Attach reviewer-side info to a pipeline-produced directory.

    The TSV must already live at
        tertiary_output/{lis_id}/snv_indel.annotated.tsv
    (the pipeline puts it there). Exomiser/LIRICAL use the conventional
    vcf_from_tsv.vcf.gz path; register records an existing VCF, and the
    background worker builds it when missing or stale.

    Phenotype is auto-loaded from
        NGS_UI/patient_phenotype/{lis_id}_{mrn}_phenotype.txt
    when the file exists. If it doesn't, the sample registers with
    empty hpo/panels and the response includes phenotype_loaded=false
    so the UI can hint about it.
    """
    started = time.perf_counter()
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

    source_sample_id = ""
    pipeline_source = TERTIARY_OUTPUT_ROOT / lis_id / "pipeline_source.json"
    if pipeline_source.is_file():
        try:
            source_info = _json.loads(pipeline_source.read_text(encoding="utf-8")) or {}
            source_sample_id = str(source_info.get("source_sample_id") or "")
        except (OSError, _json.JSONDecodeError):
            source_sample_id = ""

    roster_entry, roster_lis_id = patient_list_store.lookup_with_key(
        lis_id,
        source_sample_id,
    )
    pheno_path = PHENOTYPE_DIR / f"{lis_id}_{mrn}_phenotype.txt"
    phenotype_text = ""
    phenotype_loaded = False
    pheno_lis_candidates = patient_list_store.lookup_candidates(lis_id, roster_lis_id, source_sample_id)
    for pheno_lis_id in patient_list_store.lookup_candidates(lis_id, roster_lis_id, source_sample_id):
        candidate = PHENOTYPE_DIR / f"{pheno_lis_id}_{mrn}_phenotype.txt"
        if candidate.is_file():
            pheno_path = candidate
            phenotype_text = candidate.read_text(encoding="utf-8")
            phenotype_loaded = True
            break
    clinical_description = ""
    try:
        clinical_data = clinical_presentation_store.load(
            code=lis_id,
            mrn=mrn,
            code_candidates=pheno_lis_candidates,
        )
        clinical_description = (clinical_data.get("content") or "").strip()
    except ValueError:
        clinical_description = ""

    # Fall back to the latest roster entry for fields the reviewer
    # didn't explicitly provide. Lets the load-new-case modal stay
    # minimal — only LIS_ID + name + mrn are mandatory; 科別 /
    # 開單醫師 / 簽收時間 ride along from the xlsx upload.
    if not (department and physician and sign_received_at):
        roster_entry = roster_entry or {}
        department       = department or roster_entry.get("department", "")
        physician        = physician  or roster_entry.get("physician", "")
        sign_received_at = sign_received_at or roster_entry.get("sign_received_at", "")

    try:
        register_started = time.perf_counter()
        meta = patient_store.register(
            lis_id=lis_id, name=name, mrn=mrn, sex=sex,
            test_type=test_type, genome_build=genome_build,
            category=category,
            department=department,
            physician=physician,
            sign_received_at=sign_received_at,
            phenotype_text=phenotype_text,
            clinical_description=clinical_description,
            hpo=hpo_payload, panels=panels_payload,
        )
        _log_perf("router.samples.post.register", register_started, sample=lis_id)
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
            enqueue_started = time.perf_counter()
            from . import jobs as _jobs
            rec = _jobs._enqueue(lis_id, "exomiser_lirical", version="default")
            job_id = rec.get("job_id")
            _log_perf(
                "router.samples.post.enqueue",
                enqueue_started,
                sample=lis_id,
                job_id=job_id or "",
            )
        except Exception:
            job_id = None
    _log_perf(
        "router.samples.post.total",
        started,
        sample=lis_id,
        should_run=int(bool(should_run)),
        phenotype_loaded=int(bool(phenotype_loaded)),
    )
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


@router.get("/samples/{sample_id}/cnv")
def get_sample_cnv(sample_id: str, version: str | None = None):
    """CNV side-channel for the staged loader."""
    payload = sample_loader.load_sample_cnv(sample_id, version=version)
    if payload is None:
        raise HTTPException(404, f"sample not found: {sample_id}")
    return payload


@router.get("/samples/{sample_id}/sv")
def get_sample_sv(sample_id: str, version: str | None = None):
    """SV side-channel for the staged loader."""
    payload = sample_loader.load_sample_sv(sample_id, version=version)
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


@router.get("/samples/{sample_id}/str")
def get_sample_str(sample_id: str, version: str | None = None):
    """STR side-channel for the staged loader."""
    payload = sample_loader.load_sample_str(sample_id, version=version)
    if payload is None:
        raise HTTPException(404, f"sample not found: {sample_id}")
    return payload


@router.get("/samples/{sample_id}/pgx")
def get_sample_pgx(sample_id: str, version: str | None = None):
    """PGx / PharmCAT side-channel for the staged loader."""
    payload = sample_loader.load_sample_pgx(sample_id, version=version)
    if payload is None:
        raise HTTPException(404, f"sample not found: {sample_id}")
    return payload


@router.get("/samples/{sample_id}/secondary-snv")
def get_sample_secondary_snv(sample_id: str, version: str | None = None):
    """ACMG SF / Proactive / Carrier SNV side-channel for staged loading."""
    payload = sample_loader.load_sample_secondary_snv(sample_id, version=version)
    if payload is None:
        raise HTTPException(404, f"sample not found: {sample_id}")
    return payload


@router.get("/samples/{sample_id}/report-gene-list")
def get_sample_report_gene_list(sample_id: str, version: str | None = None):
    payload = sample_loader.load_report_gene_list(sample_id, version=version)
    if payload is None:
        raise HTTPException(404, f"sample not found: {sample_id}")
    return payload


@router.get("/samples/{sample_id}/snv-search")
def search_sample_snv(sample_id: str, genes: str, version: str | None = None):
    """Search the complete source TSV, not the compact main-screen TSV."""
    gene_list = [g.strip() for g in genes.split(",") if g.strip()]
    if not gene_list:
        raise HTTPException(400, "genes is required")
    if len(gene_list) > 100:
        raise HTTPException(400, "最多一次搜尋 100 個 genes")
    payload = sample_loader.search_snv_by_genes(
        sample_id, gene_list, version=version,
    )
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
                "genome_build", "tags", "run_date",
                "department", "physician", "sign_received_at"}
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
