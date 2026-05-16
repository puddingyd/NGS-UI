"""Path / env config.

Layout (production):
    NGS_UI/                    ← NGS_UI_HOME
    ├── NGS-UI/                ← REPO_ROOT (this git checkout)
    ├── biotools/              ← Exomiser + LIRICAL CLIs
    ├── vcf/                   ← per-sample VCFs
    ├── tertiary_output/       ← per-sample TSV + sidecars (NOT in git)
    │   └── _index.json        ← optional sample-list cache (lives next to samples)
    └── data/                  ← server runtime state (users.db, jobs/, …)

Every path is derived from NGS_UI_HOME so the whole tree can be moved
by changing one env var (or the default below). Each piece can still be
overridden individually with its own env var when the layout differs.
"""
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# Default: parent of the repo (NGS_UI/NGS-UI/ → NGS_UI/). Falls back to
# the repo itself when NGS_UI_HOME is unset and no parent layout exists,
# so dev checkouts keep working with everything inside the repo.
_default_home = REPO_ROOT.parent if (REPO_ROOT.parent / "NGS-UI").exists() else REPO_ROOT
NGS_UI_HOME = Path(os.environ.get("NGS_UI_HOME", _default_home))

TERTIARY_OUTPUT_ROOT = Path(os.environ.get(
    "TERTIARY_OUTPUT_ROOT",
    NGS_UI_HOME / "tertiary_output",
))

DATA_ROOT = Path(os.environ.get(
    "NGS_UI_DATA_ROOT",
    NGS_UI_HOME / "data",
))

INDEX_PATH = Path(os.environ.get(
    "NGS_UI_INDEX_PATH",
    NGS_UI_HOME / "tertiary_output" / "_index.json",
))

VCF_DIR = Path(os.environ.get(
    "NGS_UI_VCF_DIR",
    NGS_UI_HOME / "vcf",
))

PHENOTYPE_DIR = Path(os.environ.get(
    "NGS_UI_PHENOTYPE_DIR",
    NGS_UI_HOME / "patient_phenotype",
))

# Uploaded "未完成報告清單" xlsx files + the derived roster.json that
# maps LIS_ID → {mrn, name, test_name, department}. The 載入新個案
# modal reads the roster to auto-fill MRN / 姓名 / Test type.
PATIENT_LIST_DIR = Path(os.environ.get(
    "NGS_UI_PATIENT_LIST_DIR",
    NGS_UI_HOME / "patient_list",
))

# Reference data for phenotype scoring + HPO search: hp.obo,
# phenotype_to_genes.txt, gene_panels/*.txt (incl. user-created custom
# panels). Lives under NGS_UI_HOME so it can be swapped without a
# redeploy.
PHENO_DATA_DIR = Path(os.environ.get(
    "NGS_UI_PHENO_DATA_DIR",
    NGS_UI_HOME / "phenotype_data",
))
GENE_PANELS_DIR = Path(os.environ.get("NGS_UI_GENE_PANELS_DIR", PHENO_DATA_DIR / "gene_panels"))

BIOTOOLS_DIR = Path(os.environ.get(
    "NGS_UI_BIOTOOLS_DIR",
    NGS_UI_HOME / "biotools",
))

# OMIM annotation table (xlsx). Loaded once at first use and lazily
# reloaded when the file mtime changes. Empty value or missing file
# disables OMIM annotation; variants render with empty Disease lists.
OMIM_XLSX = Path(os.environ.get(
    "NGS_UI_OMIM_XLSX",
    NGS_UI_HOME / "OMIM" / "OMIM.xlsx",
))

JOBS_DIR = DATA_ROOT / "jobs"
JOBS_DIR.mkdir(parents=True, exist_ok=True)

FRONTEND_DIR = Path(os.environ.get(
    "FRONTEND_DIR",
    REPO_ROOT / "frontend",
))

# ---- Bioinformatics tool paths (override via env on the server) ----

EXOMISER_HOME = Path(os.environ.get(
    "EXOMISER_HOME",
    BIOTOOLS_DIR / "exomiser-cli-14.1.0",
))
EXOMISER_JAR  = Path(os.environ.get(
    "EXOMISER_JAR",
    EXOMISER_HOME / "exomiser-cli-14.1.0.jar",
))
EXOMISER_PROPS = Path(os.environ.get(
    "EXOMISER_PROPS",
    EXOMISER_HOME / "application.properties",
))
EXOMISER_DATA_HG38 = Path(os.environ.get(
    "EXOMISER_DATA_HG38",
    EXOMISER_HOME / "data" / "2508_hg38",
))
EXOMISER_DATA_HG19 = Path(os.environ.get(
    "EXOMISER_DATA_HG19",
    EXOMISER_HOME / "data" / "2508_hg19",
))

LIRICAL_HOME = Path(os.environ.get(
    "LIRICAL_HOME",
    BIOTOOLS_DIR / "lirical-cli-2.2.1",
))
LIRICAL_JAR  = Path(os.environ.get(
    "LIRICAL_JAR",
    LIRICAL_HOME / "lirical-cli-2.2.1.jar",
))

JAVA_BIN  = os.environ.get("JAVA_BIN",  "java")
JAVA_OPTS = os.environ.get("JAVA_OPTS", "-Xms4g -Xmx16g")

REDIS_URL = os.environ.get("REDIS_URL", "redis://127.0.0.1:6379/0")

# NCKU intranet EMR (HIS phenotype list + APIM consultation gateway).
# Empty value disables every EMR-related code path so the UI can run
# off-network without raising. The legacy VIP_API.sh hard-coded this
# id; we now read it from systemd's Environment= so it lives outside
# the repo. EMR endpoints are intranet-only so this is currently
# secondary defense, not auth.
EMR_CLIENT_ID = os.environ.get("NGS_UI_EMR_CLIENT_ID", "")

# DRAGEN VCF roots scanned by /api/dragen/vcfs. Each entry is a path
# (with shell globs) — files matching `*hard-filtered.vcf.gz` under any
# of them show up as candidates in the 三級分析 modal. Override via env
# (`:`-separated list of paths/globs) when DRAGEN deposits land
# elsewhere.
DRAGEN_VCF_ROOTS = [
    Path(p) for p in os.environ.get(
        "NGS_UI_DRAGEN_VCF_ROOTS",
        "/home/datalake_Raw/Novaseq",
    ).split(":") if p
]
DRAGEN_JOBS_DIR = DATA_ROOT / "jobs" / "dragen"
DRAGEN_JOBS_DIR.mkdir(parents=True, exist_ok=True)
