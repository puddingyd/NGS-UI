"""Versioned local LitVar2 lookups for variants outside review.tsv.

Post-processing remains authoritative for the compact main-screen candidate
set. Gene search can fill previously unannotated variants in background, while
an explicit card refresh may update any resolvable variant against the current
local bulk SQLite. Results live in a separate per-sample cache so reviewer
lookups never mutate the immutable pipeline TSV or its sparse annotation
overlay.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import sqlite3
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from . import litvar2_store, sample_layout, snv_gene_index


CACHE_SCHEMA_VERSION = "1"
MARKER_SCHEMA_VERSION = 1
LOOKUP_RESULT_VERSION = "2"
MAX_BATCH_VARIANTS = 1000
MAX_RAW_FALLBACK_BYTES = 100 * 1024 * 1024
_RSID_RE = re.compile(r"(?i)\brs\d+\b")
_write_lock = threading.Lock()


class Litvar2LookupError(RuntimeError):
    pass


def _read_json(path: Path) -> dict:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _marker(sample_id: str) -> dict:
    value = _read_json(sample_layout.litvar2_marker_path(sample_id))
    valid = (
        value.get("schema_version") == MARKER_SCHEMA_VERSION
        and value.get("status") == "complete"
        and value.get("scope") == "review_candidates"
    )
    return value if valid else {}


def _database_fingerprint(metadata: dict, db_path: Path) -> str:
    sha = str(metadata.get("source_sha256") or "").strip()
    dataset_date = str(metadata.get("dataset_date") or "").strip()
    schema_version = str(metadata.get("schema_version") or "legacy").strip()
    result_prefix = f"result:{LOOKUP_RESULT_VERSION}:schema:{schema_version}:"
    if sha:
        # The displayed LitVar2 version is the bulk dataset date. Include it
        # even when identical compressed content was republished with a new
        # release date, otherwise a cached payload would keep showing the old
        # date after the atomic DB switch.
        return f"{result_prefix}sha256:{sha}:date:{dataset_date}"
    try:
        stat = Path(db_path).stat()
        fallback = f"{dataset_date}:{stat.st_mtime_ns}:{stat.st_size}"
    except OSError:
        fallback = dataset_date or "unavailable"
    return result_prefix + "file:" + hashlib.sha256(
        fallback.encode("utf-8")
    ).hexdigest()


def sample_status(sample_id: str) -> dict[str, object]:
    marker = _marker(sample_id)
    metadata = litvar2_store.database_metadata(litvar2_store.LITVAR2_DB)
    available = bool(metadata.get("available"))
    return {
        "postprocessed": bool(marker),
        "automatic_enabled": bool(marker) and available,
        "db_available": available,
        "dataset_date": str(metadata.get("dataset_date") or ""),
        "marker_dataset_date": str(marker.get("dataset_date") or ""),
        "fingerprint": (
            _database_fingerprint(metadata, litvar2_store.LITVAR2_DB)
            if available else ""
        ),
    }


def _clean_variant_ids(values: Iterable[object]) -> list[str]:
    out: list[str] = []
    for raw in values:
        value = str(raw or "").strip()
        if not value or len(value) > 512 or any(ch in value for ch in "\r\n\x00"):
            continue
        if value not in out:
            out.append(value)
    return out


def _variant_id(row: dict[str, str]) -> str:
    return "-".join(str(row.get(key) or "").strip() for key in ("CHROM", "POS", "REF", "ALT"))


def _scan_rows(path: Path, wanted: set[str]) -> list[dict[str, str]]:
    if not path.is_file() or not wanted:
        return []
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            variant_id = _variant_id(row)
            if variant_id not in wanted:
                continue
            rows.append(row)
    return rows


def _resolve_rows(sample_id: str, variant_ids: list[str]) -> tuple[list[dict[str, str]], list[str]]:
    raw_tsv = sample_layout.snv_raw_tsv(sample_id)
    if not raw_tsv.is_file():
        raise FileNotFoundError(f"SNV source not found for {sample_id}")
    wanted = set(variant_ids)
    rows = snv_gene_index.query_rows_by_ids(
        raw_tsv,
        wanted,
        sample_layout.snv_gene_index_path(sample_id),
    )
    if rows is None:
        rows = _scan_rows(sample_layout.review_tsv(sample_id), wanted)
        found = {_variant_id(row) for row in rows}
        remaining = wanted - found
        try:
            raw_size = raw_tsv.stat().st_size
        except OSError:
            raw_size = 0
        if remaining and raw_size and raw_size <= MAX_RAW_FALLBACK_BYTES:
            rows.extend(_scan_rows(raw_tsv, remaining))
    found = {_variant_id(row) for row in rows}
    return rows, sorted(wanted - found)


def _append_unique(values: list[str], raw: object) -> None:
    value = str(raw or "").strip()
    if value and value not in values:
        values.append(value)


def _identifier_groups(rows: list[dict[str, str]]) -> dict[str, dict[str, list[str]]]:
    groups: dict[str, dict[str, list[str]]] = {}
    for row in rows:
        variant_id = _variant_id(row)
        if not variant_id:
            continue
        group = groups.setdefault(variant_id, {"rsids": [], "genes": [], "hgvs": []})
        for rsid in _RSID_RE.findall(str(row.get("RS_ID") or "")):
            _append_unique(group["rsids"], rsid)
        _append_unique(group["genes"], row.get("GENE") or "")
        _append_unique(group["hgvs"], row.get("HGVS_C") or "")
        _append_unique(group["hgvs"], row.get("HGVS_P") or "")
    return groups


def _identifier_fingerprint(identifiers: dict[str, list[str]]) -> str:
    normalized = {
        key: sorted({str(value).strip() for value in values if str(value).strip()})
        for key, values in identifiers.items()
    }
    raw = json.dumps(normalized, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _browser_payload(result: dict[str, object]) -> dict[str, object]:
    def source_records(payload: dict[str, object]) -> list[dict[str, object]]:
        out = []
        for raw_source in payload.get("source_records") or []:
            if not isinstance(raw_source, dict):
                continue
            out.append({
                "id": str(raw_source.get("litvar_id") or ""),
                "pmid_count": max(0, int(raw_source.get("pmids_count") or 0)),
                "url": str(raw_source.get("url") or ""),
            })
        return out

    candidates = []
    for raw in result.get("candidates") or []:
        if not isinstance(raw, dict):
            continue
        candidates.append({
            "id": str(raw.get("litvar_id") or ""),
            "rsid": str(raw.get("rsid") or ""),
            "gene": str(raw.get("gene") or ""),
            "hgvs": str(raw.get("hgvs") or ""),
            "pmid_count": max(0, int(raw.get("pmids_count") or 0)),
            "pmids": [
                str(value) for value in (raw.get("pmids") or [])
                if str(value).isdigit()
            ][:5],
            "url": str(raw.get("url") or ""),
            "merged_record_count": max(
                1, int(raw.get("merged_record_count") or 1),
            ),
            "source_records": source_records(raw),
        })
    return {
        "id": str(result.get("litvar_id") or ""),
        "rsid": str(result.get("rsid") or ""),
        "pmid_count": max(0, int(result.get("pmids_count") or 0)),
        "pmids": [
            str(value) for value in (result.get("pmids") or [])
            if str(value).isdigit()
        ][:5],
        "dataset_date": str(result.get("dataset_date") or ""),
        "match_method": str(result.get("match_method") or ""),
        "status": str(result.get("status") or "no_match"),
        "url": str(result.get("url") or ""),
        "merged_record_count": max(
            1, int(result.get("merged_record_count") or 1),
        ),
        "source_records": source_records(result),
        "candidates": candidates,
    }


def _ensure_cache(conn: sqlite3.Connection) -> None:
    conn.execute("PRAGMA journal_mode=DELETE")
    conn.execute("CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
    conn.execute(
        "CREATE TABLE IF NOT EXISTS results ("
        "variant_id TEXT NOT NULL,"
        "db_fingerprint TEXT NOT NULL,"
        "identifiers_fingerprint TEXT NOT NULL,"
        "payload_json TEXT NOT NULL,"
        "trigger TEXT NOT NULL,"
        "looked_up_at TEXT NOT NULL,"
        "PRIMARY KEY (variant_id, db_fingerprint, identifiers_fingerprint))"
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS results_current_idx "
        "ON results(db_fingerprint, variant_id, looked_up_at)"
    )
    conn.execute(
        "INSERT OR REPLACE INTO meta(key, value) VALUES ('schema_version', ?)",
        (CACHE_SCHEMA_VERSION,),
    )


def _cached_rows(
    cache_path: Path,
    variant_ids: list[str],
    db_fingerprint: str,
) -> dict[str, dict[str, object]]:
    if not cache_path.is_file() or not variant_ids:
        return {}
    placeholders = ",".join("?" for _ in variant_ids)
    try:
        with sqlite3.connect(cache_path) as conn:
            rows = conn.execute(
                f"SELECT variant_id, payload_json FROM results "
                f"WHERE db_fingerprint = ? AND variant_id IN ({placeholders}) "
                "ORDER BY looked_up_at DESC",
                [db_fingerprint, *variant_ids],
            ).fetchall()
    except sqlite3.Error:
        return {}
    out: dict[str, dict[str, object]] = {}
    for variant_id, raw in rows:
        if str(variant_id) in out:
            continue
        try:
            payload = json.loads(raw)
        except (TypeError, json.JSONDecodeError):
            continue
        if isinstance(payload, dict):
            # Cache rows created before ambiguous candidates were exposed do
            # not contain enough information for the new expandable UI.
            # Treat only those legacy ambiguous rows as a miss; hit/no-match
            # cache entries remain valid for the same DB fingerprint.
            if payload.get("status") == "ambiguous" and not payload.get("candidates"):
                continue
            out[str(variant_id)] = payload
    return out


def _write_results(
    cache_path: Path,
    rows: list[tuple[str, str, str, dict[str, object], str]],
) -> None:
    if not rows:
        return
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    with _write_lock, sqlite3.connect(cache_path, timeout=10) as conn:
        _ensure_cache(conn)
        conn.executemany(
            "INSERT OR REPLACE INTO results("
            "variant_id, db_fingerprint, identifiers_fingerprint, payload_json, trigger, looked_up_at"
            ") VALUES (?, ?, ?, ?, ?, ?)",
            [
                (variant_id, db_fp, identifiers_fp, json.dumps(payload, ensure_ascii=False), trigger, now)
                for variant_id, db_fp, identifiers_fp, payload, trigger in rows
            ],
        )
        conn.commit()


def apply_cached(variants: dict[str, dict], sample_id: str) -> int:
    if not variants:
        return 0
    metadata = litvar2_store.database_metadata(litvar2_store.LITVAR2_DB)
    if not metadata.get("available"):
        return 0
    db_fp = _database_fingerprint(metadata, litvar2_store.LITVAR2_DB)
    cached = _cached_rows(
        sample_layout.litvar2_on_demand_path(sample_id),
        list(variants),
        db_fp,
    )
    for variant_id, payload in cached.items():
        if variant_id in variants:
            variants[variant_id]["litvar2"] = dict(payload)
    return len(cached)


def lookup_variants(
    sample_id: str,
    variant_ids: Iterable[object],
    *,
    trigger: str,
    force: bool = False,
) -> dict[str, object]:
    ids = _clean_variant_ids(variant_ids)
    if not ids:
        raise ValueError("variant_ids is required")
    if len(ids) > MAX_BATCH_VARIANTS:
        raise ValueError(f"最多一次查詢 {MAX_BATCH_VARIANTS} 個 variants")
    if trigger not in {"gene_search", "manual"}:
        raise ValueError("invalid LitVar2 lookup trigger")
    status = sample_status(sample_id)
    if trigger == "gene_search" and not status["automatic_enabled"]:
        return {"results": {}, "missing": [], "eligible": False, "status": status}
    if not status["db_available"]:
        raise Litvar2LookupError("LitVar2 本地資料庫尚未建立")

    rows, unresolved = _resolve_rows(sample_id, ids)
    groups = _identifier_groups(rows)
    metadata = litvar2_store.database_metadata(litvar2_store.LITVAR2_DB)
    db_fp = _database_fingerprint(metadata, litvar2_store.LITVAR2_DB)
    cache_path = sample_layout.litvar2_on_demand_path(sample_id, for_write=True)
    cached = {} if force else _cached_rows(cache_path, list(groups), db_fp)
    results = dict(cached)
    writes: list[tuple[str, str, str, dict[str, object], str]] = []
    with litvar2_store.open_readonly(litvar2_store.LITVAR2_DB) as conn:
        for variant_id, identifiers in groups.items():
            if variant_id in results:
                continue
            result = litvar2_store.lookup_variant(
                conn,
                rsids=identifiers["rsids"],
                genes=identifiers["genes"],
                hgvs_values=identifiers["hgvs"],
            )
            payload = _browser_payload(result)
            results[variant_id] = payload
            writes.append((
                variant_id,
                db_fp,
                _identifier_fingerprint(identifiers),
                payload,
                trigger,
            ))
    _write_results(cache_path, writes)
    missing = sorted(set(unresolved) | (set(ids) - set(groups)))
    return {
        "results": results,
        "missing": missing,
        "eligible": True,
        "status": sample_status(sample_id),
        "cached": len(cached),
        "queried": len(writes),
    }
