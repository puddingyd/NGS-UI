"""Mito ClinVar lookup — runtime join from chrM-only ClinVar VCF.

Sits on top of `scripts/build_mito_clinvar_vcf.sh`'s output (default
$NGS_UI_HOME/biotools/clinvar/clinvar_mito.vcf.gz; override via env
CLINVAR_MITO_VCF). Loaded once on first call into a (POS, REF, ALT)-
keyed dict and shared across requests.

This is intentionally decoupled from the 三級分析 pipeline: the
pipeline writes pure MITOMAP-driven mito.annotated.tsv and the UI /
docx-export layer joins ClinVar at load time. Updating the slice
doesn't require re-running the pipeline; restart the uvicorn process
(or call reload()) and the new file is picked up.
"""
from __future__ import annotations

import gzip
import os
from pathlib import Path

from ..config import NGS_UI_HOME

_DEFAULT_VCF = NGS_UI_HOME / "biotools" / "clinvar" / "clinvar_mito.vcf.gz"

# Map ClinVar CLNREVSTAT → review-status star count (NCBI convention).
# Same table as scripts/annotate_clinvar.py; keep in sync if NCBI adds
# a new status.
_STAR_MAP: dict[str, int] = {
    "practice_guideline":                                   4,
    "reviewed_by_expert_panel":                             3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter":                  1,
    "criteria_provided,_conflicting_classifications":       1,
    "criteria_provided,_conflicting_interpretations":       1,
    "no_assertion_criteria_provided":                       0,
    "no_assertion_provided":                                0,
    "no_classification_provided":                           0,
}


def _stars(revstat: str) -> int:
    return _STAR_MAP.get((revstat or "").strip(), 0)


def _parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
        else:
            out[kv] = ""
    return out


def _vcf_path() -> Path:
    return Path(os.environ.get("CLINVAR_MITO_VCF", str(_DEFAULT_VCF)))


# (pos, ref, alt) → {CLNSIG, stars, CLNDN, CLNSIGCONF}
_TABLE: dict[tuple[int, str, str], dict] | None = None
_LOADED_FROM: str = ""


def _load(force: bool = False) -> dict[tuple[int, str, str], dict]:
    """Parse the mito ClinVar VCF once (lazy). Returns the cached
    lookup table. Force=True invalidates and reloads."""
    global _TABLE, _LOADED_FROM
    if _TABLE is not None and not force:
        return _TABLE

    out: dict[tuple[int, str, str], dict] = {}
    path = _vcf_path()
    if not path.is_file():
        # Cache the empty table so we don't re-stat the missing file
        # on every variant; user restarts the service after dropping
        # the VCF in.
        _TABLE = out
        _LOADED_FROM = ""
        return out

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if not raw or raw[0] == "#":
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            _chrom, pos_s, _id, ref, alt, _qual, _filter, info = parts[:8]
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            info_d = _parse_info(info)
            sig    = info_d.get("CLNSIG", "")
            star_n = _stars(info_d.get("CLNREVSTAT", ""))
            dn     = info_d.get("CLNDN", "")
            sconf  = info_d.get("CLNSIGCONF", "")
            for a in alt.split(","):
                out[(pos, ref, a)] = {
                    "CLNSIG":      sig,
                    "stars":       star_n,
                    "CLNDN":       dn,
                    "CLNSIGCONF":  sconf,
                }
    _TABLE = out
    _LOADED_FROM = str(path)
    return out


_EMPTY = {"CLNSIG": "", "stars": 0, "CLNDN": "", "CLNSIGCONF": ""}


def lookup(pos: int, ref: str, alt: str) -> dict:
    """Return the ClinVar record for one mito variant, or a blank dict
    when the variant isn't in ClinVar (most polymorphisms).
    """
    return _load().get((pos, ref, alt), _EMPTY)


def reload() -> int:
    """Force a reload from disk (e.g. after operator drops a new
    clinvar_mito.vcf.gz in place). Returns the number of records."""
    _load(force=True)
    return len(_TABLE or {})


def status() -> dict:
    """Diagnostics for the /api/healthz-ish endpoints if we ever
    surface them. Cheap; doesn't trigger a load."""
    return {
        "loaded": _TABLE is not None,
        "path":   _LOADED_FROM,
        "count":  len(_TABLE) if _TABLE is not None else 0,
    }
