import sqlite3
import time
from urllib.parse import quote

from fastapi import APIRouter, Body, Depends, File, Form, HTTPException, UploadFile
from fastapi.responses import Response

from ..auth import current_user
from ..services import (
    analyses_store,
    clinical_presentation_store,
    docx_export,
    litvar2_on_demand,
    patient_documents,
    patient_list_store,
    patient_phenotype_store,
    patient_store,
    report_store,
    sample_layout,
    sample_loader,
    test_types,
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


@router.get("/samples/{sample_id}/health-report.docx")
def get_health_report_docx(sample_id: str, sections: str = "acmg_sf,pgx"):
    """Render the health-screening DOCX with selected secondary findings
    sections. `sections` is a comma-separated list of panel keys plus pgx.
    """
    selected = [s.strip() for s in (sections or "").split(",") if s.strip()]
    try:
        blob = docx_export.build_health_docx(sample_id, sections=selected)
    except FileNotFoundError as e:
        raise HTTPException(404, str(e))

    from datetime import datetime, timezone
    from ..config import REPORT_OUTPUT_DIR
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    archive = REPORT_OUTPUT_DIR / f"{sample_id}_health_{ts}.docx"
    latest  = REPORT_OUTPUT_DIR / f"{sample_id}_health.docx"
    try:
        archive.write_bytes(blob)
        latest.write_bytes(blob)
    except OSError:
        pass

    fname = quote(f"{sample_id}_health.docx")
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
    """Unregister one sample, or optionally delete its pipeline output."""
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
    phenotype_explicit: bool = Form(False),
):
    """Attach reviewer-side info to a pipeline-produced directory.

    The unified layout marker and 03_acmg TSV must already exist.
    Exomiser/LIRICAL use the conventional 08_postprocessing
    vcf_from_tsv.vcf.gz path; register records an existing VCF, and the
    background worker builds it when missing or stale.

    Phenotype is auto-loaded from the patient-level
        NGS_UI/patient_phenotype/{mrn}_phenotype.txt
    with legacy LIS-specific filenames as fallbacks. If no file exists,
    registration can still use explicitly submitted chips or EMR fallback.
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
    pipeline_source = sample_layout.state_file(lis_id, "pipeline_source.json")
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
    pheno_lis_candidates = patient_list_store.lookup_candidates(lis_id, roster_lis_id, source_sample_id)
    try:
        phenotype_data = patient_phenotype_store.load(
            mrn=mrn,
            code=lis_id,
            code_candidates=pheno_lis_candidates,
        )
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc
    phenotype_text = str(phenotype_data.get("content") or "")
    phenotype_loaded = bool(phenotype_data)
    pheno_path = str(phenotype_data.get("path") or "")

    # The frontend normally sends its visible chips. Empty untouched chips,
    # however, mean "no browser override" so a just-entered MRN can still
    # recover the patient snapshot (or fall back to EMR). Explicitly cleared
    # chips remain authoritative via phenotype_explicit=true.
    if not phenotype_explicit and not (hpo_payload or panels_payload):
        if phenotype_data:
            hpo_payload = phenotype_data.get("hpo") or []
            panels_payload = phenotype_data.get("panels") or []
        else:
            hpo_payload = panels_payload = None
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
    # Exomiser/LIRICAL require at least one HPO term. Panels alone still
    # produce pheno_score.tsv during registration, but must not enqueue
    # a worker that cannot render valid Exomiser/LIRICAL inputs.
    job_id = None
    default_analysis = analyses_store.read_version(lis_id, "default") or {}
    should_run = bool(default_analysis.get("hpo"))
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
        "phenotype_path":   pheno_path,
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
    """Secondary-finding SNV side-channel for staged loading."""
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


@router.post("/samples/{sample_id}/litvar2/lookup")
def lookup_sample_litvar2(sample_id: str, payload: dict = Body(...)):
    """Batch lookup selected SNVs against the current local LitVar2 index."""
    try:
        return litvar2_on_demand.lookup_variants(
            sample_id,
            payload.get("variant_ids") or [],
            trigger=str(payload.get("trigger") or "gene_search"),
            force=bool(payload.get("force")),
        )
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc
    except FileNotFoundError as exc:
        raise HTTPException(404, str(exc)) from exc
    except litvar2_on_demand.Litvar2LookupError as exc:
        raise HTTPException(503, str(exc)) from exc
    except (OSError, sqlite3.Error) as exc:
        raise HTTPException(500, f"LitVar2 本地查詢失敗：{exc}") from exc


@router.put("/samples/{sample_id}/metadata")
def put_sample_metadata(
    sample_id: str,
    payload: dict,
    user: dict = Depends(current_user),
):
    """Edit a small whitelist of sample_metadata.json fields from the UI.

    Only the operator-facing identifiers + sequencing/build live here.
    HPO + selected_panels go via /api/samples/{id}/phenotype.
    """
    import json as _json
    from datetime import datetime, timezone

    sub = sample_layout.state_dir(sample_id)
    if not sub.is_dir():
        raise HTTPException(404, f"sample not found: {sample_id}")
    meta_path = sample_layout.state_file(sample_id, "sample_metadata.json")
    meta = {}
    if meta_path.exists():
        try:
            meta = _json.loads(meta_path.read_text(encoding="utf-8"))
        except _json.JSONDecodeError:
            meta = {}
    if not isinstance(meta, dict):
        meta = {}
    old_mrn = str(meta.get("mrn") or "").strip()
    requested_mrn = None
    if "mrn" in (payload or {}):
        try:
            requested_mrn = patient_documents.validate_mrn(
                str((payload or {}).get("mrn") or "")
            )
        except patient_documents.InvalidDocument as exc:
            raise HTTPException(400, str(exc)) from exc

    mrn_migration = None
    if requested_mrn and old_mrn and requested_mrn != old_mrn:
        # A shared old MRN makes a silent one-sample reassignment unsafe: it
        # could detach the other registered cases from their patient-level
        # documents. Single-patient moves remain automatic.
        if patient_documents.has_patient_data(old_mrn):
            shared_by = []
            for other_id in sample_layout.iter_sample_ids():
                if other_id == sample_id:
                    continue
                other_path = sample_layout.state_file(other_id, "sample_metadata.json")
                if not other_path.is_file():
                    continue
                try:
                    other_meta = _json.loads(other_path.read_text(encoding="utf-8")) or {}
                except (_json.JSONDecodeError, OSError):
                    continue
                if str(other_meta.get("mrn") or "").strip() == old_mrn:
                    shared_by.append(other_id)
            if shared_by:
                raise HTTPException(
                    409,
                    "舊病歷號仍被其他個案使用，為避免混合病人資料，不能自動搬移："
                    + ", ".join(sorted(shared_by)[:10]),
                )
        try:
            mrn_migration = patient_documents.move_mrn(
                old_mrn,
                requested_mrn,
                user=user,
            )
        except patient_documents.InvalidDocument as exc:
            raise HTTPException(400, str(exc)) from exc
        except patient_documents.DocumentConflict as exc:
            raise HTTPException(409, str(exc)) from exc
        except OSError as exc:
            raise HTTPException(500, f"病歷號資料搬移失敗：{exc}") from exc
    EDITABLE = {"name", "mrn", "lis_id", "sex", "test_type", "category",
                "genome_build", "tags", "run_date",
                "department", "physician", "sign_received_at"}
    for k, v in (payload or {}).items():
        if k in EDITABLE:
            meta[k] = (
                requested_mrn
                if k == "mrn"
                else
                test_types.normalize_test_type(
                    str(v or ""),
                    sample_id=str(meta.get("lis_id") or meta.get("sample_id") or sample_id),
                )
                if k == "test_type"
                else v
            )
    meta["metadata_updated_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    tmp_meta = meta_path.with_name(meta_path.name + ".tmp")
    try:
        tmp_meta.write_text(
            _json.dumps(meta, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        tmp_meta.replace(meta_path)
    except OSError as exc:
        try:
            tmp_meta.unlink(missing_ok=True)
        except OSError:
            pass
        if mrn_migration and (
            mrn_migration.get("documents") or mrn_migration.get("sidecars")
        ):
            try:
                patient_documents.move_mrn(requested_mrn, old_mrn, user=user)
            except Exception:
                pass
        raise HTTPException(500, f"個案資料儲存失敗：{exc}") from exc
    sample_loader.update_case_table_row(sample_id)
    return meta


@router.get("/samples/{sample_id}/report")
def get_report(sample_id: str):
    return report_store.load(sample_id)


@router.put("/samples/{sample_id}/report")
def put_report(sample_id: str, payload: dict, user: dict = Depends(current_user)):
    return report_store.save(sample_id, payload, user=user)


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
    import json as _json
    payload = {"category_options": list(_CATEGORY_OPTIONS), "tag_suggestions": []}
    p = sample_layout.global_cache_path("_options.json")
    if not p.exists():
        p = sample_layout.legacy_ui_root() / "_options.json"
    if p.exists():
        try:
            extra = _json.loads(p.read_text(encoding="utf-8"))
            if isinstance(extra, dict) and isinstance(extra.get("tag_suggestions"), list):
                payload["tag_suggestions"] = extra["tag_suggestions"]
        except _json.JSONDecodeError:
            pass
    return payload
