import logging
import threading

from fastapi import FastAPI
from fastapi.responses import FileResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from starlette.middleware.gzip import GZipMiddleware
from starlette.middleware.sessions import SessionMiddleware

from .config import FRONTEND_DIR, PHENO_DATA_DIR
from .routers import acmg, analyses, auth, documents, dragen, emr, files, igv, jobs, omim_capture, phenotype, phenotype_tool, pheno_engine, prs, references, samples, secondary, somatic
from .services import (
    clinvar_mito,
    hpo_ontology,
    hpo_update,
    inhouse_af_mito,
    manual_acmg,
    omim_store,
    patient_phenotype_store,
    phenotype_scorer,
    users,
    unregistered_samples,
)

app = FastAPI(title="NGS-UI", version="0.1.0")
logger = logging.getLogger(__name__)

# Large SNV JSON payloads compress well because field names repeat for
# every variant. Browsers negotiate and decompress gzip automatically.
app.add_middleware(GZipMiddleware, minimum_size=5000)

# 8 h session cookie; SameSite=Lax keeps third-party sites from
# silently impersonating the user while still letting the browser
# send the cookie on top-level GETs from internal links.
app.add_middleware(
    SessionMiddleware,
    secret_key=users.session_secret(),
    session_cookie="ngs_session",
    max_age=8 * 60 * 60,
    same_site="lax",
    https_only=False,  # internal hospital network may not be HTTPS yet
)


def _warm_caches() -> None:
    # Repair missing patient files from already-saved registrations before
    # slow reference warm-up. The phenotype load endpoint can also recover
    # one MRN immediately if it arrives before this background pass finishes.
    try:
        result = patient_phenotype_store.backfill_registered()
        logger.info("patient phenotype backfill: created=%d errors=%d", result["created"], result["errors"])
    except Exception:
        logger.exception("patient phenotype backfill failed; lookup will retry")
    # Rebuild the small active Causative/Other registry from sample JSON before
    # the slower phenotype reference warm-up starts.
    try:
        manual_acmg.backfill_observations_from_samples()
    except Exception:
        pass
    # Parse hp.obo (~17 k terms) and load phenotype_to_genes.txt (~1 M
    # rows) once so subsequent requests don't pay the I/O cost.
    hpo_ontology.load()
    phenotype_scorer.load()
    # OMIM.xlsx (~17 k rows via openpyxl) takes ~2 s on first parse;
    # warming it here moves the cost off the first /samples/{id} call.
    # Failures are silent so a missing/misconfigured xlsx doesn't
    # block startup — sample loads will degrade to empty Disease cols.
    try:
        omim_store._ensure_loaded()
    except Exception:
        pass
    # Phenotype engine (knowledge base + match kernel, a few seconds): load
    # it here so the first "預測" click does not wait for it.
    try:
        from .services.pheno_engine import engine as pheno_engine_service
        pheno_engine_service.load_engine()
    except Exception:
        logger.exception("phenotype engine warm-up failed; it will load on first use")
    # Cache the chrM-only ClinVar VCF (~3 k rows) so the first Mito
    # variant load doesn't pay the ~100 ms parse cost. Silent on
    # failure — variants just get blank ClinVar fields until the
    # operator drops the VCF in and restarts.
    try:
        clinvar_mito._load()
    except Exception:
        pass
    # The whole-genome in-house AF VCF is large; use its tabix index to
    # cache only chrM records for Mitochondria cards. Missing DB/index/tools
    # silently disable this optional annotation.
    try:
        inhouse_af_mito._load()
    except Exception:
        pass


@app.on_event("startup")
def _start_cache_warm() -> None:
    # Keep HTTP startup fast: phenotype_to_genes.txt can take more than
    # a minute to parse on the deployment host. Loaders are thread-safe
    # and still fall back to synchronous loading if a request arrives
    # before this best-effort warm-up finishes.
    unregistered_samples.start_refresh()
    thread = threading.Thread(
        target=_warm_caches,
        name="startup-cache-warm",
        daemon=True,
    )
    thread.start()
    logger.info("started background cache warm-up")


@app.get("/api/healthz")
def healthz():
    # Public, so monitoring can probe the service before login.
    return {"ok": True}


app.include_router(auth.router)
app.include_router(samples.router)
app.include_router(acmg.router)
app.include_router(analyses.router)
app.include_router(phenotype.router)
app.include_router(phenotype_tool.router)
app.include_router(documents.router)
app.include_router(jobs.router)
app.include_router(emr.router)
app.include_router(dragen.router)
app.include_router(secondary.router)
app.include_router(files.router)
app.include_router(igv.router)
app.include_router(somatic.router)
app.include_router(prs.router)
app.include_router(prs.documentation_router)
app.include_router(omim_capture.router)
app.include_router(references.router)
app.include_router(pheno_engine.router)


# Friendly direct URL for the standalone HPO/panel tool: /phenotype
# (no trailing slash) → /phenotype/ which the static mount resolves
# to frontend/phenotype/index.html. Registered before the catch-all
# StaticFiles mount so it wins.
@app.get("/phenotype")
def _phenotype_tool_redirect():
    return RedirectResponse("/phenotype/")


@app.get("/phenotype/hpo_data.json")
def _phenotype_tool_hpo_data():
    """The /phenotype/ term list: the copy rebuilt by "更新 HPO" when present."""
    updated = PHENO_DATA_DIR / hpo_update.PICKER_DATA
    path = updated if updated.is_file() else FRONTEND_DIR / "phenotype" / "hpo_data.json"
    return FileResponse(path, media_type="application/json", headers={"Cache-Control": "no-cache"})


@app.get("/omim")
def _omim_page_redirect():
    return RedirectResponse("/omim/")


# The OMIM page used to live at /omim-capture/; keep old links working.
@app.get("/omim-capture")
@app.get("/omim-capture/{rest:path}")
def _omim_capture_old_url(rest: str = ""):
    return RedirectResponse(f"/omim/{rest}", status_code=301)


if FRONTEND_DIR.exists():
    app.mount("/", StaticFiles(directory=str(FRONTEND_DIR), html=True), name="frontend")
