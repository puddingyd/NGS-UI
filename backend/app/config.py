import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

TERTIARY_OUTPUT_ROOT = Path(os.environ.get(
    "TERTIARY_OUTPUT_ROOT",
    REPO_ROOT / "tertiary_output",
))

DATA_ROOT = Path(os.environ.get(
    "NGS_UI_DATA_ROOT",
    REPO_ROOT / "data",
))

REPORTS_DIR = DATA_ROOT / "reports"
JOBS_DIR    = DATA_ROOT / "jobs"
REPORTS_DIR.mkdir(parents=True, exist_ok=True)
JOBS_DIR.mkdir(parents=True, exist_ok=True)

FRONTEND_DIR = Path(os.environ.get(
    "FRONTEND_DIR",
    REPO_ROOT / "frontend",
))

# ---- Bioinformatics tool paths (override via env on the server) ----

EXOMISER_HOME = Path(os.environ.get(
    "EXOMISER_HOME",
    "/home/n102968/biotools/exomiser-cli-14.1.0",
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
    "/home/n102968/biotools/lirical-cli-2.2.1",
))
LIRICAL_JAR  = Path(os.environ.get(
    "LIRICAL_JAR",
    LIRICAL_HOME / "lirical-cli-2.2.1.jar",
))

JAVA_BIN  = os.environ.get("JAVA_BIN",  "java")
JAVA_OPTS = os.environ.get("JAVA_OPTS", "-Xms4g -Xmx16g")

REDIS_URL = os.environ.get("REDIS_URL", "redis://127.0.0.1:6379/0")
