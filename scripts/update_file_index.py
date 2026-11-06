#!/usr/bin/env python3
"""Use the running web service's exact path overrides for the daily index.

Read only the listed discovery settings, not credentials. Fail closed if the
service environment cannot be read rather than silently scanning other roots.
"""
import os
from pathlib import Path
import subprocess
import sys

PATH_KEYS = {
    "NGS_UI_HOME", "NGS_UI_DATA_ROOT", "NGS_UI_SECONDARY_WES_FASTQ_ROOTS",
    "NGS_UI_SECONDARY_WGS_FASTQ_ROOTS", "NGS_UI_SECONDARY_FASTQ_INDEX_TTL_HOURS",
    "NGS_UI_DRAGEN_VCF_ROOTS", "NGS_UI_INHOUSE_VCF_ROOTS", "NGS_UI_DRAGEN_BAM_ROOT",
    "NGS_UI_INHOUSE_BAM_ROOT", "NGS_UI_BAM_ROOT", "NGS_UI_TERTIARY_ROOT",
    "NGS_UI_PIPELINE_OUT_ROOT", "NGS_UI_LEGACY_PIPELINE_OUT_ROOT",
    "NGS_UI_LEGACY_TERTIARY_OUTPUT_ROOT", "TERTIARY_OUTPUT_ROOT",
}


def main():
    pid = subprocess.check_output(["systemctl", "show", "ngs-ui.service", "-p", "MainPID", "--value"], text=True).strip()
    if not pid.isdigit() or int(pid) <= 0:
        raise RuntimeError("ngs-ui.service must be running to copy its discovery paths")
    environment = Path(f"/proc/{pid}/environ").read_bytes()
    for key in PATH_KEYS:
        os.environ.pop(key, None)
    for entry in environment.split(b"\0"):
        name, sep, value = entry.partition(b"=")
        if sep and name.decode(errors="replace") in PATH_KEYS:
            os.environ[name.decode()] = os.fsdecode(value)
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))
    from app.workers.file_index_update import main as update
    return update()


if __name__ == "__main__":
    raise SystemExit(main())
