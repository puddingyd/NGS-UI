from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from .config import FRONTEND_DIR
from .routers import jobs, phenotype, samples
from .services import hpo_ontology, phenotype_scorer

app = FastAPI(title="NGS-UI", version="0.1.0")


@app.on_event("startup")
def _warm_caches():
    # Parse hp.obo (~17 k terms) and load phenotype_to_genes.txt (~1 M
    # rows) once so subsequent requests don't pay the I/O cost.
    hpo_ontology.load()
    phenotype_scorer.load()


app.include_router(samples.router)
app.include_router(phenotype.router)
app.include_router(jobs.router)

if FRONTEND_DIR.exists():
    app.mount("/", StaticFiles(directory=str(FRONTEND_DIR), html=True), name="frontend")
