"""Run per-sample TSV work in forked worker processes.

Batch post-processing resolves each resource once in the parent and then
reads/rewrites every sample's working TSV. Those per-sample passes are
independent and CPU-bound (csv parsing), so they are spread over a small
process pool. ``fork`` lets workers inherit the parent's already-resolved
lookup tables without pickling them; only the per-sample return values travel
back. Forking a multi-threaded process is unsafe, so ``fork_map`` quietly
runs serially whenever other threads exist (results are identical).
"""
from __future__ import annotations

import multiprocessing
import os
import sys

_FN = None
_ITEMS: list = []


def default_workers() -> int:
    try:
        value = int(os.environ.get("NGS_UI_POSTPROC_WORKERS", "") or 0)
    except ValueError:
        value = 0
    # Post-processing runs after Nextflow has released its CPUs.
    return max(1, value or min(8, os.cpu_count() or 1))


def _single_threaded() -> bool:
    """True when forking is safe: no other (Python or native) threads."""
    try:
        return len(os.listdir("/proc/self/task")) <= 1
    except OSError:
        import threading
        return threading.active_count() <= 1


def _invoke(index: int):
    result = _FN(_ITEMS[index])
    sys.stdout.flush()
    sys.stderr.flush()
    return result


def fork_map(fn, items, workers: int | None = None) -> list:
    """``[fn(item) for item in items]``, in parallel when worthwhile.

    Results keep the input order. Any worker exception is re-raised.
    """
    global _FN, _ITEMS
    items = list(items)
    workers = default_workers() if workers is None else max(1, int(workers))
    if (workers == 1 or len(items) <= 1
            or "fork" not in multiprocessing.get_all_start_methods()
            or not _single_threaded()):
        return [fn(item) for item in items]
    sys.stdout.flush()
    sys.stderr.flush()
    _FN, _ITEMS = fn, items
    try:
        context = multiprocessing.get_context("fork")
        with context.Pool(min(workers, len(items))) as pool:
            return pool.map(_invoke, range(len(items)), chunksize=1)
    finally:
        _FN, _ITEMS = None, []
