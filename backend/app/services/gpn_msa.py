"""Join pre-computed GRCh38 GPN-MSA scores onto a compact review TSV.

The upstream table contains all possible single-nucleotide substitutions and
is queried through its tabix index.  GPN-MSA is deliberately a review-layer
annotation: this module never writes the immutable 03_acmg TSV or the sparse
annotation overlay used by raw gene search.
"""
from __future__ import annotations

import csv
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from ..config import GPN_MSA_DB

SCORE_COLUMN = "GPN_MSA_SCORE"
QUERY_BATCH_SIZE = 1000
STATUS_COMPLETE = "complete"
STATUS_MISSING_DB = "skipped_missing_db"
STATUS_MISSING_INDEX = "skipped_missing_index"
STATUS_MISSING_TABIX = "skipped_missing_tabix"
STATUS_FAILED = "failed"
_PRIMARY_CONTIGS = {str(value) for value in range(1, 23)} | {"X", "Y"}
_BASES = {"A", "C", "G", "T"}


def _file_signature(path: Path) -> dict[str, object]:
    try:
        stat = path.stat()
    except OSError:
        return {"path": str(path), "exists": False}
    return {
        "path": str(path),
        "exists": True,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def database_signature(db_path: Path | None = None) -> dict[str, object]:
    """Return the data/index/tabix identity recorded in review manifests."""
    db = Path(db_path or GPN_MSA_DB)
    requested_tabix = os.environ.get("TABIX_BIN") or "tabix"
    resolved_tabix = shutil.which(requested_tabix)
    return {
        "database": _file_signature(db),
        "index": _file_signature(Path(f"{db}.tbi")),
        "tabix": {
            "requested": requested_tabix,
            "path": resolved_tabix or "",
            "exists": bool(resolved_tabix),
        },
    }


def validate_database(db_path: Path | None = None) -> tuple[Path, str]:
    """Validate the fixed score table and return ``(path, tabix_binary)``."""
    db = Path(db_path or GPN_MSA_DB)
    missing = [path for path in (db, Path(f"{db}.tbi")) if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "required GPN-MSA input not found: " + ", ".join(map(str, missing))
        )
    requested_tabix = os.environ.get("TABIX_BIN") or "tabix"
    tabix = shutil.which(requested_tabix)
    if not tabix:
        raise RuntimeError("tabix is required to query the GPN-MSA score table")
    return db, tabix


def _status(
    status: str,
    message: str,
    db: Path,
    *,
    rows: int = 0,
    queryable_variants: int = 0,
    annotated_variants: int = 0,
    annotated_rows: int = 0,
) -> dict[str, object]:
    return {
        "status": status,
        "message": message,
        "database": str(db),
        "rows": rows,
        "queryable_variants": queryable_variants,
        "annotated_variants": annotated_variants,
        "annotated_rows": annotated_rows,
        "database_available": status not in {
            STATUS_MISSING_DB,
            STATUS_MISSING_INDEX,
            STATUS_MISSING_TABIX,
        },
    }


def _warn(message: str) -> None:
    print(f"[gpn-msa] WARNING: {message}", file=sys.stderr, flush=True)


def _unavailable_status(db: Path, exc: Exception) -> tuple[str, str]:
    index = Path(f"{db}.tbi")
    if not db.is_file():
        return STATUS_MISSING_DB, f"GPN-MSA DB not found: {db}"
    if not index.is_file():
        return STATUS_MISSING_INDEX, f"GPN-MSA index not found: {index}"
    if isinstance(exc, RuntimeError) and "tabix is required" in str(exc):
        return STATUS_MISSING_TABIX, str(exc)
    return STATUS_FAILED, f"GPN-MSA setup validation failed: {exc}"


def _normalise_key(row: dict[str, str]) -> tuple[str, int, str, str] | None:
    chrom = str(row.get("CHROM") or "").strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    chrom = chrom.upper()
    ref = str(row.get("REF") or "").strip().upper()
    alt = str(row.get("ALT") or "").strip().upper()
    if chrom not in _PRIMARY_CONTIGS or ref not in _BASES or alt not in _BASES:
        return None
    try:
        pos = int(str(row.get("POS") or "").strip())
    except ValueError:
        return None
    if pos < 1:
        return None
    return chrom, pos, ref, alt


def query_scores(
    db_path: Path,
    keys: set[tuple[str, int, str, str]],
    *,
    tabix_bin: str,
) -> dict[tuple[str, int, str, str], str]:
    """Query exact alleles in bounded tabix batches.

    Regions are deduplicated by position because the table returns all three
    substitutions at each queried base; exact REF/ALT matching happens after
    tabix returns the five-column rows.
    """
    wanted = set(keys)
    positions = sorted(
        {(chrom, pos) for chrom, pos, _ref, _alt in wanted},
        key=lambda item: (
            int(item[0]) if item[0].isdigit() else 23 if item[0] == "X" else 24,
            item[1],
        ),
    )
    found: dict[tuple[str, int, str, str], str] = {}
    for start in range(0, len(positions), QUERY_BATCH_SIZE):
        regions = [
            f"{chrom}:{pos}-{pos}"
            for chrom, pos in positions[start:start + QUERY_BATCH_SIZE]
        ]
        proc = subprocess.run(
            [tabix_bin, str(db_path), *regions],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if proc.returncode != 0:
            detail = proc.stderr.strip() or f"exit status {proc.returncode}"
            raise RuntimeError(f"GPN-MSA tabix query failed: {detail}")
        for line in proc.stdout.splitlines():
            if not line or line.startswith("#"):
                continue
            columns = line.split("\t")
            if len(columns) < 5:
                continue
            parsed = _normalise_key({
                "CHROM": columns[0], "POS": columns[1],
                "REF": columns[2], "ALT": columns[3],
            })
            if parsed not in wanted:
                continue
            try:
                score = float(columns[4])
            except ValueError:
                continue
            if not math.isfinite(score):
                continue
            found[parsed] = f"{score:g}"
    return found


def annotate_review_tsv(
    tsv_path: Path,
    db_path: Path | None = None,
    *,
    required: bool = False,
) -> dict[str, object]:
    """Atomically add ``GPN_MSA_SCORE`` to a compact review TSV.

    Normal post-processing is best-effort: unavailable deployment data and
    query/write failures return an explicit status without changing the input
    TSV. ``required=True`` remains available for deployment diagnostics.
    """
    tsv = Path(tsv_path)
    db = Path(db_path or GPN_MSA_DB)
    try:
        db, tabix = validate_database(db)
    except (FileNotFoundError, RuntimeError, OSError) as exc:
        if required:
            raise
        status, message = _unavailable_status(db, exc)
        _warn(message)
        return _status(status, message, db)

    keys: set[tuple[str, int, str, str]] = set()
    with tsv.open("r", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if not reader.fieldnames:
            raise RuntimeError(f"empty review TSV: {tsv}")
        for row in reader:
            key = _normalise_key(row)
            if key is not None:
                keys.add(key)

    try:
        scores = query_scores(db, keys, tabix_bin=tabix)
        with tsv.open("r", encoding="utf-8", newline="") as source:
            reader = csv.DictReader(source, delimiter="\t")
            fields = list(reader.fieldnames or [])
            if SCORE_COLUMN not in fields:
                fields.append(SCORE_COLUMN)
            fd, tmp_name = tempfile.mkstemp(
                prefix=f".{tsv.name}.", suffix=".tmp", dir=str(tsv.parent), text=True
            )
            rows = annotated = 0
            with os.fdopen(fd, "w", encoding="utf-8", newline="") as destination:
                writer = csv.DictWriter(
                    destination,
                    fieldnames=fields,
                    delimiter="\t",
                    extrasaction="ignore",
                    lineterminator="\n",
                )
                writer.writeheader()
                for row in reader:
                    rows += 1
                    value = scores.get(_normalise_key(row), "")
                    row[SCORE_COLUMN] = value
                    annotated += bool(value)
                    writer.writerow(row)
            os.replace(tmp_name, tsv)
    except Exception as exc:
        if "tmp_name" in locals():
            try:
                os.unlink(tmp_name)
            except OSError:
                pass
        if required:
            raise
        message = f"GPN-MSA annotation failed for {tsv}: {exc}"
        _warn(message)
        return _status(
            STATUS_FAILED,
            message,
            db,
            queryable_variants=len(keys),
        )
    return _status(
        STATUS_COMPLETE,
        "",
        db,
        rows=rows,
        queryable_variants=len(keys),
        annotated_variants=len(scores),
        annotated_rows=annotated,
    )
