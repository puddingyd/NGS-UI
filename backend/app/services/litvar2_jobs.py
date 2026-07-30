"""Background update state shared by systemd and the tertiary modal."""
from __future__ import annotations

import fcntl
import json
import os
import subprocess
import sys
import traceback
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from ..config import (
    LITVAR2_DIR,
    LITVAR2_UPDATE_LOG_PATH,
    LITVAR2_UPDATE_STATE_PATH,
    REPO_ROOT,
)
from . import litvar2_store


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _atomic_state(payload: dict) -> None:
    LITVAR2_DIR.mkdir(parents=True, exist_ok=True)
    tmp = LITVAR2_UPDATE_STATE_PATH.with_suffix(".json.tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    os.replace(tmp, LITVAR2_UPDATE_STATE_PATH)


def _read_state() -> dict:
    try:
        value = json.loads(LITVAR2_UPDATE_STATE_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _pid_alive(pid) -> bool:
    try:
        value = int(pid)
        os.kill(value, 0)
        return True
    except (TypeError, ValueError, OSError):
        return False


@contextmanager
def _exclusive_lock(path: Path, *, blocking: bool):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = path.open("a+")
    flags = fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB)
    try:
        fcntl.flock(handle.fileno(), flags)
        yield handle
    finally:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        finally:
            handle.close()


def status() -> dict:
    state = _read_state()
    running = (
        state.get("state") in {"queued", "running"}
        and _pid_alive(state.get("pid"))
    )
    if state.get("state") in {"queued", "running"} and not running:
        state = {
            **state,
            "state": "failed",
            "step": "worker-exited",
            "finished_at": _now(),
            "error": state.get("error") or "LitVar2 updater process is no longer running",
        }
        _atomic_state(state)
    db = litvar2_store.database_metadata()
    return {
        **state,
        "running": running,
        "database": db,
        "dataset_date": db.get("dataset_date", ""),
    }


def start_background_update() -> dict:
    """Spawn the same updater used by the monthly systemd timer."""
    LITVAR2_DIR.mkdir(parents=True, exist_ok=True)
    spawn_lock = LITVAR2_DIR / ".spawn.lock"
    with _exclusive_lock(spawn_lock, blocking=True):
        current = status()
        if current.get("running"):
            return {**current, "already_running": True}
        run_id = uuid.uuid4().hex
        queued = {
            "run_id": run_id,
            "trigger": "manual",
            "state": "queued",
            "step": "queued",
            "created_at": _now(),
            "started_at": None,
            "finished_at": None,
            "pid": None,
            "error": None,
        }
        _atomic_state(queued)
        env = os.environ.copy()
        backend_path = str(REPO_ROOT / "backend")
        existing_pythonpath = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = (
            backend_path
            if not existing_pythonpath
            else backend_path + os.pathsep + existing_pythonpath
        )
        log_handle = LITVAR2_UPDATE_LOG_PATH.open("a", buffering=1, encoding="utf-8")
        try:
            proc = subprocess.Popen(
                [
                    sys.executable,
                    "-m",
                    "app.workers.litvar2_update",
                    "--trigger",
                    "manual",
                    "--run-id",
                    run_id,
                ],
                cwd=str(REPO_ROOT),
                env=env,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
        finally:
            log_handle.close()
        latest = _read_state()
        if latest.get("run_id") == run_id and latest.get("state") == "queued":
            latest["pid"] = proc.pid
            _atomic_state(latest)
        return status()


def run_update(*, trigger: str, run_id: str = "") -> dict:
    """Run an update in the foreground; safe for systemd or detached worker."""
    LITVAR2_DIR.mkdir(parents=True, exist_ok=True)
    lock_path = LITVAR2_DIR / ".update.lock"
    try:
        lock_context = _exclusive_lock(lock_path, blocking=False)
        with lock_context:
            rid = run_id or uuid.uuid4().hex
            state = {
                "run_id": rid,
                "trigger": trigger,
                "state": "running",
                "step": "checking-source",
                "created_at": _read_state().get("created_at") or _now(),
                "started_at": _now(),
                "finished_at": None,
                "pid": os.getpid(),
                "error": None,
            }
            _atomic_state(state)

            def progress(step: str, fields: dict[str, object]) -> None:
                current = _read_state()
                if current.get("run_id") != rid:
                    current = state.copy()
                current.update(fields)
                current.update(
                    run_id=rid,
                    trigger=trigger,
                    state="running",
                    step=step,
                    pid=os.getpid(),
                )
                _atomic_state(current)
                print(
                    f"[litvar2-update] step={step} "
                    + " ".join(f"{key}={value}" for key, value in fields.items()),
                    flush=True,
                )

            result = litvar2_store.update_database(progress=progress)
            done = {
                **_read_state(),
                "run_id": rid,
                "trigger": trigger,
                "state": "done",
                "step": result.get("action") or "done",
                "finished_at": _now(),
                "pid": os.getpid(),
                "error": None,
                "result": result,
            }
            _atomic_state(done)
            print(
                f"[litvar2-update] done action={result.get('action')} "
                f"dataset_date={result.get('dataset_date', '')}",
                flush=True,
            )
            return done
    except BlockingIOError:
        current = status()
        print("[litvar2-update] another updater already holds the lock", flush=True)
        return {**current, "already_running": True}
    except Exception as exc:
        failed = {
            **_read_state(),
            "run_id": run_id or _read_state().get("run_id") or uuid.uuid4().hex,
            "trigger": trigger,
            "state": "failed",
            "step": "failed",
            "finished_at": _now(),
            "pid": os.getpid(),
            "error": str(exc),
        }
        _atomic_state(failed)
        traceback.print_exc()
        raise
