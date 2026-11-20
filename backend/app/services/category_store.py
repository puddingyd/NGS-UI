"""Shared category vocabulary and recoverable updates of registered cases."""
from __future__ import annotations

import fcntl
import hashlib
import inspect
import json
import os
import stat
import tempfile
import threading
from contextlib import contextmanager
from functools import wraps
from pathlib import Path

from fastapi import HTTPException

from .. import config
from . import sample_layout

DEFAULT_CATEGORIES = [
    "Neurology", "Endocrinology", "MCA", "Nephrology", "GI", "Metabolism",
    "AIR", "Hematology", "Oncology", "Ophthalmology", "Musculoskeletal",
    "Dermatology", "CV", "ENT", "Asymptomatic",
]
_LOCK = threading.RLock()
_LOCAL = threading.local()


def _path(name="categories.json") -> Path:
    return Path(config.DATA_ROOT) / name


def _read(path):
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Invalid JSON object: {path.name}")
    return value


def _write(path, value):
    """Keep complete old/new files even on interrupted writes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            if path.exists():
                os.fchmod(stream.fileno(), stat.S_IMODE(path.stat().st_mode))
            json.dump(value, stream, ensure_ascii=False, indent=2)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(name, path)
        _sync_directory(path.parent)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def _sync_directory(directory):
    fd = os.open(directory, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _finish(transaction):
    rollback = transaction.get("rollback", False)
    for change in transaction["changes"]:
        path = sample_layout.state_file(change["sample_id"], "sample_metadata.json")
        meta = _read(path)
        meta["category"] = change["old"] if rollback else change["new"]
        _write(path, meta)
    settings = transaction["before"] if rollback else transaction["after"]
    if settings is None:
        _path().unlink(missing_ok=True)
    else:
        _write(_path(), settings)
    # Derived rows contain metadata signatures. Removing the cache makes the
    # next case-list request rebuild it, including in another web worker.
    sample_layout.global_cache_path("_case_table.json").unlink(missing_ok=True)
    _path("categories.pending.json").unlink()
    _sync_directory(_path().parent)


@contextmanager
def _sample_locks(changes):
    # Pipeline publication replaces whole state directories. Share its locks
    # so a completed analysis cannot publish an older copy of the category.
    if not changes:
        yield
        return
    from ..workers.dragen_run import _acquire_sample_locks, _release_sample_locks
    ids = {change["sample_id"] for change in changes}
    ids.update(sample_layout.source_sample_id(sid) for sid in list(ids))
    try:
        handles = _acquire_sample_locks(list(ids))
    except RuntimeError as exc:
        raise HTTPException(409, "部分受影響個案正在分析或更新，請稍後再儲存分類") from exc
    try:
        yield
    finally:
        _release_sample_locks(handles)


@contextmanager
def locked():
    """Serialize category transactions with all metadata read/modify/writes."""
    with _LOCK:
        if getattr(_LOCAL, "active", False):
            yield
            return
        _path().parent.mkdir(parents=True, exist_ok=True)
        with _path("categories.lock").open("a") as stream:
            fcntl.flock(stream, fcntl.LOCK_EX)
            _LOCAL.active = True
            try:
                pending = _path("categories.pending.json")
                if pending.exists():
                    transaction = _read(pending)
                    with _sample_locks(transaction["changes"]):
                        _finish(transaction)
                yield
            finally:
                _LOCAL.active = False
                fcntl.flock(stream, fcntl.LOCK_UN)


def serialized(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        with locked():
            return function(*args, **kwargs)
    # FastAPI must resolve future annotations against the original module.
    wrapped.__signature__ = inspect.signature(function, eval_str=True)
    return wrapped


def _snapshot():
    settings = _read(_path()) if _path().exists() else None
    names = list(settings["names"] if settings is not None else DEFAULT_CATEGORIES)
    samples = []
    counts = {}
    for sid in sorted(sample_layout.iter_sample_ids()):
        path = sample_layout.state_file(sid, "sample_metadata.json")
        if not path.exists():
            continue
        meta = _read(path)  # An unreadable registered case must not be skipped.
        name = meta.get("category") or ""
        if not isinstance(name, str):
            raise ValueError("個案分類格式錯誤")
        samples.append((sid, meta))
        if name:
            counts[name] = counts.get(name, 0) + 1
            if name not in names:
                names.append(name)
    revision = hashlib.sha256(json.dumps(names, ensure_ascii=False).encode()).hexdigest()
    ordered = sorted(names, key=lambda name: -counts.get(name, 0))
    return settings, samples, names, {
        "category_options": ordered,
        "category_counts": {name: counts.get(name, 0) for name in ordered},
        "category_revision": revision,
    }


@serialized
def options():
    return _snapshot()[3]


def validate_selection(value, revision=None, *, previous=""):
    """Called under locked(); stale tabs cannot undo a rename or removal."""
    if value is not None and not isinstance(value, str):
        raise HTTPException(400, "分類名稱格式錯誤")
    name = (value or "").strip()
    # Existing clients/scripts may omit a revision. Before the first managed
    # edit preserve the historical free-text input behavior.
    if not _path().exists() and (not revision or not name or name in DEFAULT_CATEGORIES or name == (previous or "")):
        return name
    settings = _read(_path()) if _path().exists() else {"names": _snapshot()[2]}
    current_revision = settings.get("revision")
    if revision and current_revision and revision != current_revision and name != (previous or ""):
        raise HTTPException(409, "分類清單已更新，請重新載入個案後再儲存")
    if name and name not in settings["names"] and name != (previous or ""):
        raise HTTPException(409, "此分類已改名或移除，請重新載入分類清單")
    return name


def revision():
    if _path().exists():
        return _read(_path()).get("revision")
    return hashlib.sha256(json.dumps(DEFAULT_CATEGORIES, ensure_ascii=False).encode()).hexdigest()


@serialized
def save(payload, *, user=None):
    before, samples, names, snapshot = _snapshot()
    if payload.get("category_revision") != snapshot["category_revision"]:
        raise HTTPException(409, "分類清單已被更新，請關閉後重新開啟管理分類")
    entries = payload.get("categories")
    if not isinstance(entries, list) or len(entries) > 200:
        raise HTTPException(400, "分類清單格式錯誤，最多可設定 200 個分類")
    result, used, mapping = [], set(), {name: "" for name in names}
    for entry in entries:
        if not isinstance(entry, dict) or not isinstance(entry.get("name"), str):
            raise HTTPException(400, "請填寫分類名稱")
        name, original = entry["name"].strip(), entry.get("original")
        if not name or len(name) > 100 or any(ord(c) < 32 for c in name):
            raise HTTPException(400, "分類名稱須為 1–100 字且不可包含換行或控制字元")
        if name.casefold() in used:
            raise HTTPException(400, "分類名稱不可重複（不分大小寫）")
        used.add(name.casefold())
        if original is not None:
            if not isinstance(original, str) or original not in names or original in result:
                raise HTTPException(400, "原分類資料錯誤")
            result.append(original)
            mapping[original] = name
    # Keep the canonical tie order; the frequency-sorted editor does not
    # silently change the order when counts tie again later.
    new_names = [mapping[name] for name in names if mapping[name]]
    new_names.extend(entry["name"].strip() for entry in entries if entry.get("original") is None)
    revision = hashlib.sha256(json.dumps(new_names, ensure_ascii=False).encode()).hexdigest()
    after = {"names": new_names, "revision": revision, "updated_by": (user or {}).get("username", "")}
    changes = [{"sample_id": sid, "old": meta.get("category"), "new": mapping[meta["category"]]}
               for sid, meta in samples if meta.get("category") and mapping[meta["category"]] != meta["category"]]
    transaction = {"before": before, "after": after, "changes": changes}
    with _sample_locks(changes):
        _write(_path("categories.pending.json"), transaction)
        try:
            _finish(transaction)
        except OSError:
            transaction["rollback"] = True
            _write(_path("categories.pending.json"), transaction)
            _finish(transaction)
            raise
    result = _snapshot()[3]
    result.update(category_changes=mapping, updated_samples=len(changes))
    return result
