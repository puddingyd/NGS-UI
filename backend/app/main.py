import logging
import threading

from fastapi import FastAPI
from fastapi.responses import RedirectResponse
from fastapi.staticfiles import StaticFiles
from starlette.middleware.gzip import GZipMiddleware
from starlette.middleware.sessions import SessionMiddleware

from .config import FRONTEND_DIR
from .routers import acmg, analyses, auth, dragen, emr, igv, jobs, phenotype, phenotype_tool, samples, secondary
from .services import (
    clinvar_mito,
    hpo_ontology,
    inhouse_af_mito,
    manual_acmg,
    omim_store,
    phenotype_scorer,
    users,
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
app.include_router(jobs.router)
app.include_router(emr.router)
app.include_router(dragen.router)
app.include_router(secondary.router)
app.include_router(igv.router)


# Friendly direct URL for the standalone HPO/panel tool: /phenotype
# (no trailing slash) → /phenotype/ which the static mount resolves
# to frontend/phenotype/index.html. Registered before the catch-all
# StaticFiles mount so it wins.
@app.get("/phenotype")
def _phenotype_tool_redirect():
    return RedirectResponse("/phenotype/")


if FRONTEND_DIR.exists():
    app.mount("/", StaticFiles(directory=str(FRONTEND_DIR), html=True), name="frontend")
