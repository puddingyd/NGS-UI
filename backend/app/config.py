"""Path / env config.

Layout (production):
    NGS_UI/                    ← NGS_UI_HOME
    ├── NGS-UI/                ← REPO_ROOT (this git checkout)
    ├── biotools/              ← Exomiser + LIRICAL CLIs
    ├── vcf/                   ← per-sample VCFs
    ├── tertiary_output/       ← per-sample TSV + sidecars (NOT in git)
    ├── data/                  ← server runtime state (users.db, jobs/, …)
    └── _index.json            ← optional sample-list cache

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
    NGS_UI_HOME / "_index.json",
))

VCF_DIR = Path(os.environ.get(
    "NGS_UI_VCF_DIR",
    NGS_UI_HOME / "vcf",
))

BIOTOOLS_DIR = Path(os.environ.get(
    "NGS_UI_BIOTOOLS_DIR",
    NGS_UI_HOME / "biotools",
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
