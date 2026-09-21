"""Shared, persistent pending-sample index; patient details are loaded on selection."""
from __future__ import annotations

import fcntl
import json
import logging
import os
import threading
import time
from contextlib import contextmanager
from pathlib import Path

from .. import config
from . import patient_list_store, patient_phenotype_store, sample_layout, test_types

logger = logging.getLogger(__name__)
CACHE_VERSION = 1
MAX_AGE_SECONDS = 300
_REFRESH_LOCK = threading.Lock()
_last_error_at = 0.0


def _path() -> Path:
    return Path(config.DATA_ROOT) / "unregistered_samples.json"


def _roots() -> list[str]:
    return [str(p) for p in (sample_layout.unified_root(), sample_layout.legacy_ui_root(), sample_layout.legacy_pipeline_root())]


def _json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def _snapshot() -> dict:
    data = _json(_path())
    if (data.get("version") != CACHE_VERSION or data.get("roots") != _roots()
            or not isinstance(data.get("items"), list)
            or not isinstance(data.get("scanned_at"), (int, float))
            or not isinstance(data.get("updated_at"), (int, float))):
        return {}
    return data


@contextmanager
def _write_lock():
    path = _path()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.with_suffix(".lock").open("a") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _publish(items: list[dict], scanned_at: float) -> dict:
    data = {"version": CACHE_VERSION, "roots": _roots(), "scanned_at": scanned_at,
            "updated_at": time.time(), "items": sorted(items, key=lambda r: r["mtime"], reverse=True)}
    path = _path()
    temporary = path.with_suffix(".json.tmp")
    try:
        temporary.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return data


def describe(sample_id: str) -> dict | None:
    """Read each layout/source JSON once; never parse phenotype or stat source VCF."""
    patient_phenotype_store.check_token("LIS_ID", sample_id, required=True)
    sample = sample_layout.unified_sample_dir(sample_id)
    post = sample / sample_layout.POSTPROCESSING_DIRNAME
    marker = _json(sample_layout.scoped_file(post, sample_id, "layout.json"))
    unified = marker.get("layout_version") in sample_layout.SUPPORTED_LAYOUT_VERSIONS
    state = post if unified else sample_layout.legacy_ui_root() / sample_id
    # Reject registered cases before looking for their annotation inputs.
    if sample_layout.scoped_file(state, sample_id, "sample_metadata.json").is_file():
        return None
    legacy_raw = state / "snv_indel.annotated.tsv"
    if not unified and not legacy_raw.is_file():
        return None
    source = _json(sample_layout.scoped_file(state, sample_id, "pipeline_source.json"))
    source_id = str(source.get("source_sample_id") or sample_id)
    raw = legacy_raw
    if unified:
        raw_value = str(marker.get("raw_tsv") or source.get("source_path") or "")
        raw = Path(raw_value) if raw_value else sample / "03_acmg" / f"{source_id}.snv_indel.acmg.tsv"
        if raw_value and not raw.is_absolute():
            raw = sample / raw
        if not raw.is_file():
            raw = sample_layout._first_file(sample / "03_acmg", f"{source_id}.snv_indel.acmg.tsv", "*.snv_indel.acmg.tsv")
    if not raw.is_file():
        return None
    pipeline = sample if sample.is_dir() else sample_layout.legacy_pipeline_root() / sample_id
    mtime = pipeline.stat().st_mtime if pipeline.is_dir() else state.stat().st_mtime
    return {"lis_id": sample_id, "source_sample_id": str(source.get("source_sample_id") or ""),
            "pipeline_type": str(source.get("pipeline_type") or ""),
            "source_vcf_path": str(source.get("source_vcf_path") or ""),
            "tsv_size": raw.stat().st_size, "mtime": mtime}


def refresh(*, force: bool = False) -> dict:
    """Single writer across API processes and pipeline workers; atomic publication."""
    started = time.perf_counter()
    with _write_lock():
        current = _snapshot()
        if not force and current and time.time() - current["scanned_at"] < MAX_AGE_SECONDS:
            return current
        rows = []
        seen = 0
        for sample_id in sample_layout.iter_sample_ids():
            seen += 1
            try:
                row = describe(sample_id)
            except (FileNotFoundError, ValueError):
                continue  # concurrently removed or non-sample directory
            if row:
                rows.append(row)
        data = _publish(rows, time.time())
    logger.info("pending-index scan seconds=%.3f directories=%d pending=%d", time.perf_counter() - started, seen, len(rows))
    return data


def start_refresh() -> bool:
    """Warm on startup / serve stale data while one background refresh runs."""
    if time.time() - _last_error_at < 30:
        return False
    if not _REFRESH_LOCK.acquire(blocking=False):
        return True

    def run():
        global _last_error_at
        try:
            refresh()
            _last_error_at = 0.0
        except Exception:
            _last_error_at = time.time()
            logger.exception("pending-index background refresh failed")
        finally:
            _REFRESH_LOCK.release()

    threading.Thread(target=run, name="pending-sample-index", daemon=True).start()
    return True


def refresh_samples(sample_ids: list[str]) -> None:
    """Update only changed entries; do not postpone the next full reconciliation."""
    try:
        with _write_lock():
            current = _snapshot()
            if not current:
                return  # first startup/request will build a complete index
            rows = {row["lis_id"]: row for row in current["items"]}
            for sample_id in set(sample_ids):
                try:
                    row = describe(sample_id)
                except (FileNotFoundError, ValueError):
                    row = None
                rows.pop(sample_id, None)
                if row:
                    rows[sample_id] = row
            _publish(list(rows.values()), current["scanned_at"])
    except Exception:
        logger.exception("pending-index incremental refresh failed")


def _enrich(row: dict, roster: dict) -> dict:
    entry, key = patient_list_store.lookup_with_key(row["lis_id"], row.get("source_sample_id", ""), roster=roster)
    fields = ("mrn", "name", "department", "physician", "sign_received_at")
    summary = {name: entry.get(name, "") for name in fields} if entry else None
    if summary is not None:
        summary["test_type"] = test_types.normalize_test_type(entry.get("test_type", ""), sample_id=row["lis_id"], default="")
    return {**row, "roster": summary, "roster_lis_id": key}


def listing(*, force: bool = False) -> dict:
    from . import dragen_jobs
    started = time.perf_counter()
    data = refresh(force=True) if force else _snapshot()
    stale = not data or time.time() - data["scanned_at"] >= MAX_AGE_SECONDS
    refreshing = start_refresh() if stale else False
    roster_started = time.perf_counter()
    roster = patient_list_store.load_roster()
    active = dragen_jobs.active_sample_ids()
    items = [_enrich(row, roster) for row in data.get("items", []) if row["lis_id"] not in active]
    logger.info("pending-index response seconds=%.3f roster_jobs_seconds=%.3f pending=%d refreshing=%s", time.perf_counter() - started, time.perf_counter() - roster_started, len(items), refreshing)
    return {"items": items, "updated_at": data.get("updated_at"), "refreshing": refreshing,
            "error": "清單更新失敗，請按更新清單重試。" if stale and not refreshing else ""}


def detail(sample_id: str, *, mrn: str = "") -> dict:
    """Fresh validation and phenotype lookup for exactly one selected sample."""
    from . import dragen_jobs
    started = time.perf_counter()
    row = describe(sample_id)
    if not row:
        raise FileNotFoundError("此個案已登錄、尚未完成或已移除，請更新清單。")
    if sample_id in dragen_jobs.active_sample_ids():
        raise RuntimeError("此個案正在三級分析，請完成後再載入。")
    row = _enrich(row, patient_list_store.load_roster())
    mrn = patient_phenotype_store.check_token("MRN", mrn, required=False)
    patient_mrn = mrn or (row["roster"] or {}).get("mrn", "")
    candidates = patient_list_store.lookup_candidates(sample_id, row["roster_lis_id"], row["source_sample_id"])
    # An explicitly entered MRN must not fall back to a different patient's LIS.
    pheno = patient_phenotype_store.load(mrn=patient_mrn, **({} if mrn else {"code": sample_id, "code_candidates": candidates}))
    row["phenotype"] = {key: pheno.get(key) for key in ("path", "mrn", "hpo", "panels")} if pheno else None
    row["source_vcf_size"] = 0
    if row["source_vcf_path"]:
        try:
            row["source_vcf_size"] = Path(row["source_vcf_path"]).stat().st_size
        except OSError:
            pass
    logger.info("pending-index detail seconds=%.3f", time.perf_counter() - started)
    return row
