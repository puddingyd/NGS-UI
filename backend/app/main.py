from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from starlette.middleware.sessions import SessionMiddleware

from .config import FRONTEND_DIR
from .routers import analyses, auth, jobs, phenotype, samples
from .services import hpo_ontology, phenotype_scorer, users

app = FastAPI(title="NGS-UI", version="0.1.0")

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


@app.on_event("startup")
def _warm_caches():
    # Parse hp.obo (~17 k terms) and load phenotype_to_genes.txt (~1 M
    # rows) once so subsequent requests don't pay the I/O cost.
    hpo_ontology.load()
    phenotype_scorer.load()


@app.get("/api/healthz")
def healthz():
    # Public, so monitoring can probe the service before login.
    return {"ok": True}


app.include_router(auth.router)
app.include_router(samples.router)
app.include_router(analyses.router)
app.include_router(phenotype.router)
app.include_router(jobs.router)

if FRONTEND_DIR.exists():
    app.mount("/", StaticFiles(directory=str(FRONTEND_DIR), html=True), name="frontend")
