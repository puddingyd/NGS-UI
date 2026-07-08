"""DRAGEN sex-karyotype sidecar helpers."""
from __future__ import annotations

import gzip
import re
from pathlib import Path


def read_karyotype(path: Path) -> str:
    """Read DRAGEN's estimatedSexKaryotype header, preserving calls such as XXY or X."""
    if not path.is_file():
        return ""
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("##estimatedSexKaryotype="):
                    value = line.split("=", 1)[1].strip().upper()
                    return re.sub(r"[^A-Z0-9+-]", "", value)
                if not line.startswith("##"):
                    break
    except OSError:
        return ""
    return ""


def load_sample_ploidy(sample_dir: Path) -> dict:
    path = sample_dir / "ploidy.vcf.gz"
    karyotype = read_karyotype(path)
    return {
        "karyotype": karyotype,
        "source": path.name if karyotype else "",
    }
