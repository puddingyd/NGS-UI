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
REPORTS_DIR.mkdir(parents=True, exist_ok=True)

FRONTEND_DIR = Path(os.environ.get(
    "FRONTEND_DIR",
    REPO_ROOT / "frontend",
))
