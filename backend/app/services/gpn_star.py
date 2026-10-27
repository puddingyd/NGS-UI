"""Local, allele-exact GPN-Star V/M annotation of compact review TSVs only.

Canonical hg38 Parquet uses 1-based positions. No network access, inference,
raw/overlay writes, threshold, or ACMG/tier logic belongs in this module.
"""
from __future__ import annotations

import csv
import importlib.metadata
import json
import math
import os
import re
import sys
import tempfile
from pathlib import Path

from ..config import GPN_STAR_DIR

MODEL = "gpn-star-hg38-m447-200m"
SCORE_FIELD = "llr_calibrated"
# Keep the original M column and manifest names for existing installations.
SCORE_COLUMN = "GPN_STAR_SCORE"
V_MODEL = "gpn-star-hg38-v100-200m"
MODELS = (MODEL, V_MODEL)
SCORE_COLUMNS = {MODEL: SCORE_COLUMN, V_MODEL: "GPN_STAR_V_SCORE"}
PRIMARY_CONTIGS = {str(n) for n in range(1, 23)} | {"X", "Y"}
KEY_COLUMNS = ["chrom", "pos", "ref", "alt"]


def manifest_path(root: Path, model: str = MODEL) -> Path:
    if model not in MODELS:
        raise ValueError(f"unsupported GPN-Star model: {model}")
    return root / ("manifest.json" if model == MODEL else "manifest.v.json")


def _file_signature(path: Path) -> dict:
    try:
        st = path.stat()
        return {"path": str(path), "size": st.st_size, "mtime_ns": st.st_mtime_ns}
    except OSError:
        return {"path": str(path), "exists": False}


def database_signature(db_path: Path | None = None, *, model: str = MODEL) -> dict:
    root = Path(db_path or GPN_STAR_DIR)
    try:
        version = importlib.metadata.version("polars")
    except importlib.metadata.PackageNotFoundError:
        version = ""
    return {
        "model": model, "assembly": "hg38", "score_field": SCORE_FIELD,
        "schema_version": 1, "polars": version,
        "manifest": _file_signature(manifest_path(root, model)),
        "shards": {chrom: _file_signature(shard_path(root, chrom, model))
                   for chrom in sorted(PRIMARY_CONTIGS)},
    }


def shard_path(root: Path, chrom: str, model: str = MODEL) -> Path:
    manifest_path(root, model)  # validate model before constructing a path
    return root / model / "llr" / f"llr_chr{chrom}.parquet"


def validate_database(db_path: Path | None = None, *, model: str = MODEL) -> dict:
    root = Path(db_path or GPN_STAR_DIR)
    manifest = json.loads(manifest_path(root, model).read_text(encoding="utf-8"))
    if (manifest.get("model") != model or manifest.get("assembly") != "hg38"
            or manifest.get("score_field") != SCORE_FIELD
            or manifest.get("schema_version") != 1
            or not re.fullmatch(r"[0-9a-f]{40}", manifest.get("revision", ""))):
        raise ValueError("GPN-Star manifest model/assembly/score/revision mismatch")
    files = manifest.get("files", [])
    if len(files) != 24 or {f["chrom"] for f in files} != PRIMARY_CONTIGS:
        raise ValueError("GPN-Star manifest must contain chromosomes 1-22, X, Y")
    for item in files:
        path = shard_path(root, item["chrom"], model)
        if item["path"] != path.relative_to(root).as_posix():
            raise ValueError("unexpected GPN-Star shard path")
        if not re.fullmatch(r"[0-9a-f]{64}", item.get("sha256", "")):
            raise ValueError("missing GPN-Star shard checksum")
        if path.stat().st_size != item["size"]:
            raise ValueError(f"GPN-Star shard size mismatch: {path}")
    return manifest


def _normalise_key(row: dict) -> tuple[str, int, str, str] | None:
    chrom = str(row.get("CHROM") or "").strip().upper().removeprefix("CHR")
    ref = str(row.get("REF") or "").strip().upper()
    alt = str(row.get("ALT") or "").strip().upper()
    if chrom not in PRIMARY_CONTIGS or ref not in {"A", "C", "G", "T"} or alt not in {"A", "C", "G", "T"} or ref == alt:
        return None
    try:
        pos = int(str(row.get("POS") or "").strip())
    except ValueError:
        return None
    return (chrom, pos, ref, alt) if pos > 0 else None


def query_scores(root: Path, keys: set[tuple], model: str = MODEL) -> dict[tuple, str]:
    import polars as pl

    found = {}
    for chrom in sorted({key[0] for key in keys}):
        wanted = {key for key in keys if key[0] == chrom}
        positions = sorted({key[1] for key in wanted})
        # Projection + predicate pushdown limits reads to relevant row groups.
        # One query per chromosome avoids repeatedly scanning a large shard.
        scores = (pl.scan_parquet(shard_path(root, chrom, model))
                  .filter(pl.col("pos").is_in(positions))
                  .select(*KEY_COLUMNS, SCORE_FIELD)
                  .collect(engine="streaming"))
        for c, pos, ref, alt, value in scores.iter_rows():
            key = (c, pos, ref, alt)
            if key not in wanted or value is None:
                continue
            score = float(value)
            if not math.isfinite(score):
                continue
            if key in found:
                raise ValueError(f"duplicate GPN-Star allele: {key}")
            # Nine significant digits round-trip the canonical Float32 value.
            found[key] = format(score, ".9g")
    return found


def annotate_review_tsv(tsv_path: Path, db_path: Path | None = None, *, required: bool = False, model: str = MODEL) -> dict:
    root = Path(db_path or GPN_STAR_DIR)
    column = SCORE_COLUMNS[model]
    tsv = Path(tsv_path)
    stats = {"status": "failed", "model": model, "assembly": "hg38",
             "score_field": SCORE_FIELD, "database": str(root), "rows": 0,
             "queryable_variants": 0, "annotated_variants": 0, "annotated_rows": 0,
             "unsupported_rows": 0, "unmatched_variants": 0}
    tmp_name = None
    try:
        manifest = validate_database(root, model=model)
        stats["revision"] = manifest["revision"]
        with tsv.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or not {"CHROM", "POS", "REF", "ALT"}.issubset(reader.fieldnames):
                raise ValueError("review TSV lacks CHROM/POS/REF/ALT")
            keys = set()
            for row in reader:
                key = _normalise_key(row)
                if key is not None:
                    keys.add(key)
                else:
                    stats["unsupported_rows"] += 1
        stats["queryable_variants"] = len(keys)
        scores = query_scores(root, keys, model)
        stats["annotated_variants"] = len(scores)
        stats["unmatched_variants"] = len(keys) - len(scores)
        with tsv.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames)
            if column not in fields:
                fields.append(column)
            fd, tmp_name = tempfile.mkstemp(prefix=f".{tsv.name}.", suffix=".tmp", dir=tsv.parent)
            with os.fdopen(fd, "w", encoding="utf-8", newline="") as output:
                writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for row in reader:
                    value = scores.get(_normalise_key(row), "")
                    row[column] = value
                    writer.writerow(row)
                    stats["rows"] += 1
                    stats["annotated_rows"] += value != ""
            os.replace(tmp_name, tsv)
            tmp_name = None
        stats["status"] = "complete"
    except Exception as exc:
        if required:
            raise
        if isinstance(exc, FileNotFoundError):
            stats["status"] = "skipped_missing_db" if not manifest_path(root, model).is_file() else "skipped_missing_shard"
        elif isinstance(exc, ImportError):
            stats["status"] = "skipped_missing_polars"
        stats["message"] = str(exc)
        print(f"[gpn-star {model}] WARNING: {stats['status']}: {exc}", file=sys.stderr, flush=True)
    finally:
        if tmp_name is not None:
            Path(tmp_name).unlink(missing_ok=True)
    return stats
