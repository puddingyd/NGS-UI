"""Runtime lookup of chrM records in the in-house AF sites VCF.

The published in-house AF database is a whole-genome, bgzip-compressed sites
VCF, so scanning it on every Mitochondria request would be prohibitively
expensive.  This service uses the VCF's tabix index to extract only chrM once,
then keeps a small (POS, REF, ALT)-keyed table in memory.

chrM uses a carrier-frequency contract rather than diploid allele arithmetic:

  INHOUSE_AC      samples carrying the ALT
  INHOUSE_AN      callable mitochondrial genomes (one per sample)
  INHOUSE_AF      INHOUSE_AC / INHOUSE_AN
  INHOUSE_NHOM    homoplasmic carriers
  INHOUSE_HET_MT  heteroplasmic carriers

Missing database/index/query tools silently produce an empty lookup.  The
Mitochondria adapter treats that as annotation unavailable and does not use
these fields for tiering or ACMG classification.
"""
from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from pathlib import Path
from typing import Iterable

from ..config import INHOUSE_AF_DB

_CONTIG_ALIASES = ("chrM", "MT", "M", "chrMT")
_DEFAULT_BCFTOOLS_SIF = Path(
    "/home/pipeline/nextflow_containers/bcftools_1.23.1.sif"
)
_MAX_UNINDEXED_SCAN_BYTES = 64 * 1024 * 1024


def _parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for field in (info or "").split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            out[key] = value
    return out


def _per_alt(raw: str, index: int) -> str:
    values = (raw or "").split(",")
    if not values:
        return ""
    if index < len(values):
        return values[index]
    # Number=1 fields (notably INHOUSE_AN) apply to every ALT.
    return values[0] if len(values) == 1 else ""


def _to_float(raw: str):
    try:
        value = str(raw or "").strip()
        return None if value in ("", ".") else float(value)
    except (TypeError, ValueError):
        return None


def _to_int(raw: str):
    value = _to_float(raw)
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _minimal_repr(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Trim shared suffix/prefix without requiring a reference FASTA."""
    ref = (ref or "").upper()
    alt = (alt or "").upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def _table_from_lines(lines: Iterable[str]) -> dict[tuple[int, str, str], dict]:
    out: dict[tuple[int, str, str], dict] = {}
    for raw in lines:
        if not raw or raw[0] == "#":
            continue
        parts = raw.rstrip("\n").split("\t")
        if len(parts) < 8 or parts[0] not in _CONTIG_ALIASES:
            continue
        try:
            pos = int(parts[1])
        except ValueError:
            continue
        ref = parts[3].upper()
        info = _parse_info(parts[7])
        for alt_index, alt_raw in enumerate(parts[4].split(",")):
            alt = alt_raw.upper()
            if alt in ("", ".", "*", "<NON_REF>", "<*>"):
                continue
            ac = _to_int(_per_alt(info.get("INHOUSE_AC", ""), alt_index))
            an = _to_int(_per_alt(info.get("INHOUSE_AN", ""), alt_index))
            af = _to_float(_per_alt(info.get("INHOUSE_AF", ""), alt_index))
            if ac is None and an is None and af is None:
                continue
            key = _minimal_repr(pos, ref, alt)
            out[key] = {
                "inhouse_ac": ac,
                "inhouse_an": an,
                "inhouse_af": af,
                "inhouse_nhom": _to_int(
                    _per_alt(info.get("INHOUSE_NHOM", ""), alt_index)
                ),
                "inhouse_het_mt": _to_int(
                    _per_alt(info.get("INHOUSE_HET_MT", ""), alt_index)
                ),
            }
    return out


def _bcftools_command() -> list[str] | None:
    binary = os.environ.get("BCFTOOLS_BIN", "bcftools")
    if shutil.which(binary):
        return [binary]

    sif_raw = os.environ.get("BCFTOOLS_SIF", "")
    sif = Path(sif_raw) if sif_raw else _DEFAULT_BCFTOOLS_SIF
    if sif.is_file() and shutil.which("apptainer"):
        binds = os.environ.get("APPTAINER_BIND", "/home")
        return ["apptainer", "exec", "--bind", binds, str(sif), "bcftools"]
    return None


def _indexed_region_lines(path: Path) -> tuple[list[str], str]:
    """Return chrM VCF records plus the query method used."""
    bcftools = _bcftools_command()
    if bcftools is not None:
        for contig in _CONTIG_ALIASES:
            try:
                proc = subprocess.run(
                    bcftools + ["view", "-H", "-r", contig, str(path)],
                    capture_output=True,
                    text=True,
                    timeout=60,
                    check=False,
                )
            except (OSError, subprocess.TimeoutExpired):
                continue
            lines = proc.stdout.splitlines()
            if lines:
                return lines, f"bcftools:{contig}"

    tabix = shutil.which(os.environ.get("TABIX_BIN", "tabix"))
    if tabix:
        for contig in _CONTIG_ALIASES:
            try:
                proc = subprocess.run(
                    [tabix, str(path), contig],
                    capture_output=True,
                    text=True,
                    timeout=60,
                    check=False,
                )
            except (OSError, subprocess.TimeoutExpired):
                continue
            lines = proc.stdout.splitlines()
            if lines:
                return lines, f"tabix:{contig}"
    return [], ""


def _small_file_lines(path: Path) -> tuple[list[str], str]:
    """Development/test fallback; never stream a production-sized WGS DB."""
    try:
        if path.stat().st_size > _MAX_UNINDEXED_SCAN_BYTES:
            return [], ""
    except OSError:
        return [], ""
    opener = gzip.open if str(path).endswith(".gz") else open
    lines: list[str] = []
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                if raw.startswith(_CONTIG_ALIASES):
                    lines.append(raw.rstrip("\n"))
    except OSError:
        return [], ""
    return lines, "small-file-scan"


_TABLE: dict[tuple[int, str, str], dict] | None = None
_LOADED_FROM = ""
_QUERY_METHOD = ""


def _load(force: bool = False) -> dict[tuple[int, str, str], dict]:
    global _TABLE, _LOADED_FROM, _QUERY_METHOD
    if _TABLE is not None and not force:
        return _TABLE

    path = Path(INHOUSE_AF_DB)
    if not path.is_file():
        _TABLE = {}
        _LOADED_FROM = ""
        _QUERY_METHOD = ""
        return _TABLE

    lines, method = _indexed_region_lines(path)
    if not lines:
        lines, method = _small_file_lines(path)
    _TABLE = _table_from_lines(lines)
    _LOADED_FROM = str(path) if _TABLE else ""
    _QUERY_METHOD = method if _TABLE else ""
    return _TABLE


_EMPTY = {
    "inhouse_ac": None,
    "inhouse_an": None,
    "inhouse_af": None,
    "inhouse_nhom": None,
    "inhouse_het_mt": None,
}


def lookup(pos: int, ref: str, alt: str) -> dict:
    """Return chrM carrier frequency/counts for one normalized variant."""
    try:
        key = _minimal_repr(int(pos), ref, alt)
    except (TypeError, ValueError):
        return _EMPTY
    try:
        return _load().get(key, _EMPTY)
    except Exception:
        # Optional annotation must never make the Mitochondria endpoint fail.
        return _EMPTY


def reload() -> int:
    """Reload the chrM slice after the database is atomically updated."""
    _load(force=True)
    return len(_TABLE or {})


def status() -> dict:
    return {
        "loaded": _TABLE is not None,
        "path": _LOADED_FROM,
        "query_method": _QUERY_METHOD,
        "count": len(_TABLE) if _TABLE is not None else 0,
    }
