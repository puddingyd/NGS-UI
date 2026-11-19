"""Time-weighted progress for tertiary analysis jobs.

The progress bar used to be driven by fixed per-process weights matched
against Nextflow's ANSI status lines.  Nextflow truncates long process names
from the left (``ANN…NNOTATE:VEP_ANNOTATE (80)``), so the heaviest processes
(VEP, Pangolin, PARSE_CSQ) were never matched and the bar sat at ~16% for most
of the run, then jumped to 82% when Nextflow exited.

This module replaces that with:

* Nextflow run with ``-ansi-log false -with-trace <job>/trace.txt``.  Task
  submissions come from the plain stdout lines and completions from the trace
  file; both carry the full, untruncated process name.
* Work estimates in *seconds*: every process (and every post-processing step)
  has an expected per-task duration.  Progress = finished work / expected
  work, with partial credit for running tasks.
* A timing profile learned from previous successful jobs, keyed by
  ``<pipeline_type>:<seq_type>`` (``nckuh:WES``, ``dragen:WGS`` …).  It holds
  per-process task durations, per-step post-processing durations and the
  per-sample wall time of the Nextflow and post-processing phases, which
  decides how much of the bar each phase gets.  Built-in defaults (measured
  on a 133-sample NCKUH WES batch) are used until a profile exists.
"""
from __future__ import annotations

import json
import math
import os
import re
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

PROFILE_VERSION = 1
EMA_ALPHA = 0.5
RUNNING_CREDIT_CAP = 0.9

# ---- defaults -----------------------------------------------------------

# Expected seconds per task.  Relative size matters more than absolute value;
# the learned profile replaces these after the first successful job.
DEFAULT_PROCESS_SECONDS: dict[str, float] = {
    "ADD_DRAGEN_TAG": 30.0,
    "ADD_CALLERS_TAG": 30.0,
    "FILTER_FOR_ANNOTATION": 60.0,
    "MITO_VEP": 90.0,
    "MITO_PARSE": 45.0,
    "STR_PARSE_DRAGEN": 40.0,
    "STR_PREPARE_NCKUH": 40.0,
    "STR_PARSE_NCKUH": 40.0,
    "PREPARE_CNV_DRAGEN": 40.0,
    "PREPARE_CNV_NCKUH": 40.0,
    "CNVKIT_TO_BED": 40.0,
    "ANNOTSV_CNV_DRAGEN": 180.0,
    "ANNOTSV_CNV_NCKUH": 180.0,
    "ANNOTSV_CNV_NCKUH_WES": 180.0,
    "ANNOTSV_CNV_NCKUH_WGS": 180.0,
    "PREPARE_SV_DRAGEN": 40.0,
    "PREPARE_SV_NCKUH": 40.0,
    "ANNOTSV_SV_DRAGEN": 480.0,
    "ANNOTSV_SV_NCKUH": 480.0,
    "VEP_ANNOTATE": 1920.0,
    "PANGOLIN_SCORE": 1200.0,
    "PARSE_CSQ": 960.0,
    "ACMG_CLASSIFY": 480.0,
    "PGX_STELLARPGX": 60.0,
    "PGX_HLA_EXTRACT": 60.0,
    "PGX_OPTITYPE": 480.0,
    "PGX_PHARMCAT": 240.0,
    "PGX_PARSE": 120.0,
}
UNKNOWN_PROCESS_SECONDS = 60.0

# Processes that only exist for one pipeline type / sequencing type.  Used to
# keep the default (no-profile) denominator close to what will really run.
_DRAGEN_ONLY = re.compile(r"(?:_DRAGEN$|^ADD_DRAGEN_TAG$)")
_NCKUH_ONLY = re.compile(r"(?:_NCKUH(?:_WES|_WGS)?$|^ADD_CALLERS_TAG$|^CNVKIT_TO_BED$)")
_GENERIC_CNV = {"ANNOTSV_CNV_NCKUH"}  # superseded by the _WES/_WGS variants

# Average seconds per post-processing step (NCKUH WES, 133-sample batch).
DEFAULT_POSTPROCESSING_STEP_SECONDS: dict[str, float] = {
    "clinvar-latest": 10.0,
    "genebe": 240.0,
    "spliceai": 97.0,
    "giab-strata": 10.0,
    "inhouse-af": 31.0,
    "mane-refseq": 1.0,
    "litvar2": 22.0,
    "snv-overlay": 2.0,
    "review-tsv": 79.0,
    "gene-index": 1.0,
    "cnv-rescue": 60.0,
}
# Wall seconds per sample for each phase (same batch: Nextflow 64,189 s and
# post-processing 65,436 s for 133 samples).
DEFAULT_PHASE_SECONDS_PER_SAMPLE = {"nextflow": 483.0, "postprocessing": 492.0}

# ---- small helpers ------------------------------------------------------

_ANSI_ESCAPE_RE = re.compile(r"\x1b(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
_STDOUT_TASK_RE = re.compile(
    r"\[(?P<hash>[0-9a-f]{2}/[0-9a-f]{6})\]\s+"
    r"(?P<event>Submitted|Cached|Re-submitted)\s+process\s+>\s+(?P<name>.+?)\s*$"
)
_DURATION_PART_RE = re.compile(r"(\d+(?:\.\d+)?)\s*(ms|d|h|m|s)")


def simple_process_name(name: str) -> str:
    """``ANNOTATE_SNV:SNV_ANNOTATE:VEP_ANNOTATE (12)`` → ``VEP_ANNOTATE``."""
    text = _ANSI_ESCAPE_RE.sub("", str(name or "")).strip()
    text = re.sub(r"\s*\(.*\)\s*$", "", text)
    return text.rsplit(":", 1)[-1].strip().lstrip("…").strip()


def kind_key(pipeline_type: str, seq_type: str) -> str:
    return f"{(pipeline_type or '').lower()}:{(seq_type or '').upper()}"


def parse_duration(value: str) -> float | None:
    """Parse Nextflow trace durations: ``1h 2m 3s``, ``4.5s``, ``350ms``."""
    text = str(value or "").strip()
    if not text or text == "-":
        return None
    try:
        return float(text) / 1000.0  # raw trace (milliseconds)
    except ValueError:
        pass
    total = 0.0
    found = False
    for number, unit in _DURATION_PART_RE.findall(text):
        found = True
        n = float(number)
        total += {"ms": n / 1000.0, "s": n, "m": n * 60, "h": n * 3600, "d": n * 86400}[unit]
    return total if found else None


def parse_stdout_task(line: str) -> tuple[str, str, str] | None:
    """Return ``(hash, event, process)`` for ``-ansi-log false`` task lines."""
    match = _STDOUT_TASK_RE.search(_ANSI_ESCAPE_RE.sub("", line or ""))
    if not match:
        return None
    event = "submitted" if match["event"] != "Cached" else "cached"
    return match["hash"], event, simple_process_name(match["name"])


def read_trace(path: Path, offset: int = 0) -> tuple[list[dict], int]:
    """Read complete trace rows appended after ``offset``.

    Returns the parsed rows and the new byte offset.  A trailing partial line
    is left for the next call.  The header line is recognised and skipped.
    """
    try:
        with Path(path).open("rb") as handle:
            handle.seek(offset)
            data = handle.read()
    except OSError:
        return [], offset
    end = data.rfind(b"\n")
    if end < 0:
        return [], offset
    chunk = data[: end + 1].decode("utf-8", errors="replace")
    rows: list[dict] = []
    header = None
    if offset:
        header = read_trace.headers.get(str(path))
    for line in chunk.splitlines():
        if not line.strip():
            continue
        cols = line.split("\t")
        if header is None or cols[:2] == ["task_id", "hash"]:
            header = cols
            read_trace.headers[str(path)] = header
            continue
        rows.append(dict(zip(header, cols)))
    return rows, offset + end + 1


read_trace.headers = {}  # type: ignore[attr-defined]


# ---- profile ------------------------------------------------------------

def load_profile(path: Path) -> dict:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {"version": PROFILE_VERSION, "kinds": {}}
    if not isinstance(data, dict) or data.get("version") != PROFILE_VERSION:
        return {"version": PROFILE_VERSION, "kinds": {}}
    data.setdefault("kinds", {})
    return data


def save_profile(path: Path, profile: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(profile, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def _ema(old: float | None, new: float) -> float:
    if old is None or not math.isfinite(old) or old <= 0:
        return float(new)
    return (1 - EMA_ALPHA) * float(old) + EMA_ALPHA * float(new)


def phase_seconds_per_sample(profile: dict, kind: str, phase: str) -> float:
    value = (((profile.get("kinds") or {}).get(kind) or {}).get("phases") or {}).get(phase)
    try:
        value = float(value)
    except (TypeError, ValueError):
        value = 0.0
    return value if value > 0 else DEFAULT_PHASE_SECONDS_PER_SAMPLE[phase]


def postprocessing_step_seconds(profile: dict, kind: str) -> dict[str, float]:
    learned = (((profile.get("kinds") or {}).get(kind) or {}).get("postprocessing_steps") or {})
    out = dict(DEFAULT_POSTPROCESSING_STEP_SECONDS)
    for step, seconds in learned.items():
        try:
            if float(seconds) > 0:
                out[step] = float(seconds)
        except (TypeError, ValueError):
            continue
    return out


def update_profile(
    profile: dict,
    kind: str,
    *,
    trace_rows: list[dict] | None = None,
    sample_count: int = 0,
    nextflow_seconds: float | None = None,
    postprocessing_seconds: float | None = None,
    step_seconds: dict[str, list[float]] | None = None,
) -> dict:
    """Fold one successful job's measurements into ``profile`` (in place)."""
    entry = profile.setdefault("kinds", {}).setdefault(kind, {})
    if trace_rows and sample_count > 0:
        durations: dict[str, list[float]] = {}
        counts: dict[str, int] = {}
        for row in trace_rows:
            status = str(row.get("status") or "").upper()
            process = simple_process_name(row.get("name", ""))
            if not process or status not in {"COMPLETED", "CACHED"}:
                continue
            counts[process] = counts.get(process, 0) + 1
            if status != "COMPLETED":
                continue
            seconds = parse_duration(row.get("realtime", "")) or parse_duration(row.get("duration", ""))
            if seconds is not None and seconds > 0:
                durations.setdefault(process, []).append(seconds)
        processes = entry.setdefault("processes", {})
        for process, count in counts.items():
            item = processes.setdefault(process, {})
            values = durations.get(process)
            if values:
                item["seconds"] = round(_ema(item.get("seconds"), sum(values) / len(values)), 2)
            item["tasks_per_sample"] = round(_ema(item.get("tasks_per_sample"), count / sample_count), 3)
    if sample_count > 0:
        phases = entry.setdefault("phases", {})
        for phase, seconds in (("nextflow", nextflow_seconds), ("postprocessing", postprocessing_seconds)):
            # Sub-second phases (dry runs, everything cached) carry no signal.
            if seconds and seconds / sample_count >= 1.0:
                phases[phase] = round(_ema(phases.get(phase), seconds / sample_count), 2)
    if step_seconds:
        steps = entry.setdefault("postprocessing_steps", {})
        for step, values in step_seconds.items():
            if values:
                steps[step] = round(_ema(steps.get(step), sum(values) / len(values)), 2)
    entry["updated_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    entry["jobs"] = int(entry.get("jobs") or 0) + 1
    return profile


# ---- phase bands --------------------------------------------------------

def phase_bands(profile: dict, kinds: list[str], *, with_prs: bool) -> dict[str, tuple[float, float]]:
    """Split 0-100% into phases by expected wall time.

    ``kinds`` has one entry per sample.  Fixed bands: start 0-2%, and the
    promote / PRS tail.  Nextflow and post-processing share the rest in
    proportion to their learned per-sample wall time.
    """
    nf = sum(phase_seconds_per_sample(profile, kind, "nextflow") for kind in kinds) or 1.0
    pp = sum(phase_seconds_per_sample(profile, kind, "postprocessing") for kind in kinds) or 1.0
    start = 2.0
    work_end = 94.0 if with_prs else 97.0
    split = start + (work_end - start) * nf / (nf + pp)
    bands = {
        "pre": (0.0, start),
        "nextflow": (start, split),
        "postprocessing": (split, work_end),
        "promote": (work_end, 95.0 if with_prs else 99.0),
    }
    if with_prs:
        bands["prs"] = (95.0, 99.0)
    return bands


def band_pct(bands: dict[str, tuple[float, float]], phase: str, fraction: float) -> float:
    lo, hi = bands[phase]
    fraction = max(0.0, min(1.0, float(fraction)))
    return round(lo + (hi - lo) * fraction, 1)


# ---- Nextflow tracker ---------------------------------------------------

@dataclass
class _Process:
    seconds: float
    expected: float  # expected task count for this batch
    submitted: int = 0
    completed: int = 0
    cached: int = 0
    failed: int = 0
    first_seen: float | None = None
    done_logged: bool = False


@dataclass
class NextflowProgress:
    """Track one Nextflow run; thread-safe."""

    kinds: list[str]
    profile: dict
    with_pgx: bool = True
    _processes: dict[str, _Process] = field(default_factory=dict)
    _running: dict[str, tuple[str, float]] = field(default_factory=dict)  # hash → (process, start)
    _finished: set[str] = field(default_factory=set)
    _lock: threading.RLock = field(default_factory=threading.RLock)
    _trace_offset: int = 0

    def __post_init__(self) -> None:
        expected: dict[str, float] = {}
        seconds: dict[str, list[float]] = {}
        for kind in self.kinds:
            pipeline_type, _sep, seq_type = kind.partition(":")
            learned = (((self.profile.get("kinds") or {}).get(kind) or {}).get("processes") or {})
            if learned:
                table = {
                    name: (float(item.get("seconds") or DEFAULT_PROCESS_SECONDS.get(name, UNKNOWN_PROCESS_SECONDS)),
                           float(item.get("tasks_per_sample") or 1.0))
                    for name, item in learned.items()
                }
            else:
                table = {
                    name: (secs, 1.0)
                    for name, secs in DEFAULT_PROCESS_SECONDS.items()
                    if self._default_applies(name, pipeline_type, seq_type)
                }
            for name, (secs, per_sample) in table.items():
                if not self.with_pgx and name.startswith("PGX_"):
                    continue
                expected[name] = expected.get(name, 0.0) + per_sample
                seconds.setdefault(name, []).append(secs)
        for name, count in expected.items():
            values = seconds[name]
            self._processes[name] = _Process(seconds=sum(values) / len(values), expected=count)

    @staticmethod
    def _default_applies(name: str, pipeline_type: str, seq_type: str) -> bool:
        if name in _GENERIC_CNV:
            return False
        if pipeline_type == "dragen" and _NCKUH_ONLY.search(name):
            return False
        if pipeline_type != "dragen" and _DRAGEN_ONLY.search(name):
            return False
        if name.endswith("_WES") and seq_type != "WES":
            return False
        if name.endswith("_WGS") and seq_type != "WGS":
            return False
        if name == "CNVKIT_TO_BED" and seq_type != "WGS":
            return False
        return True

    def _proc(self, name: str) -> _Process:
        proc = self._processes.get(name)
        if proc is None:
            proc = _Process(seconds=DEFAULT_PROCESS_SECONDS.get(name, UNKNOWN_PROCESS_SECONDS), expected=0.0)
            self._processes[name] = proc
        return proc

    # -- events --
    def on_stdout(self, line: str, now: float | None = None) -> list[dict]:
        """Feed one Nextflow stdout line; return process start/done events."""
        parsed = parse_stdout_task(line)
        if parsed is None:
            return []
        task_hash, event, name = parsed
        now = time.monotonic() if now is None else now
        events: list[dict] = []
        with self._lock:
            if task_hash in self._finished:
                return []
            proc = self._proc(name)
            if proc.first_seen is None:
                proc.first_seen = now
                events.append({"process": name, "event": "start"})
            if event == "cached":
                self._finished.add(task_hash)
                self._running.pop(task_hash, None)
                proc.cached += 1
                events.extend(self._maybe_done(name, now))
            elif task_hash not in self._running:
                proc.submitted += 1
                self._running[task_hash] = (name, now)
        return events

    def on_trace_row(self, row: dict, now: float | None = None) -> list[dict]:
        task_hash = str(row.get("hash") or "")
        status = str(row.get("status") or "").upper()
        name = simple_process_name(row.get("name", ""))
        if not task_hash or not name or status not in {"COMPLETED", "CACHED", "FAILED", "ABORTED"}:
            return []
        now = time.monotonic() if now is None else now
        with self._lock:
            if task_hash in self._finished:
                return []
            self._finished.add(task_hash)
            started = self._running.pop(task_hash, None)
            proc = self._proc(name)
            if proc.first_seen is None:
                proc.first_seen = started[1] if started else now
            if status == "CACHED":
                proc.cached += 1
            elif status == "COMPLETED":
                proc.completed += 1
            else:
                proc.failed += 1
            return self._maybe_done(name, now)

    def poll_trace(self, path: Path, now: float | None = None) -> list[dict]:
        with self._lock:
            rows, self._trace_offset = read_trace(path, self._trace_offset)
            events: list[dict] = []
            for row in rows:
                events.extend(self.on_trace_row(row, now))
            return events

    def _maybe_done(self, name: str, now: float) -> list[dict]:
        proc = self._processes[name]
        finished = proc.completed + proc.cached
        target = max(proc.expected, proc.submitted + proc.cached)
        if proc.done_logged or finished < max(1.0, round(target)) or any(
            p == name for p, _t in self._running.values()
        ):
            return []
        proc.done_logged = True
        elapsed = now - (proc.first_seen if proc.first_seen is not None else now)
        return [{"process": name, "event": "done", "done": finished,
                 "total": int(round(target)), "elapsed": elapsed}]

    # -- queries --
    def fraction(self, now: float | None = None) -> float:
        now = time.monotonic() if now is None else now
        with self._lock:
            total = 0.0
            done = 0.0
            running_by_proc: dict[str, list[float]] = {}
            for name, started in self._running.values():
                running_by_proc.setdefault(name, []).append(now - started)
            for name, proc in self._processes.items():
                expected = max(proc.expected, proc.submitted + proc.cached + proc.failed)
                # Cached tasks cost nothing in this run: drop them from both sides.
                total += max(0.0, expected - proc.cached) * proc.seconds
                done += (proc.completed + proc.failed) * proc.seconds
                for elapsed in running_by_proc.get(name, []):
                    done += min(RUNNING_CREDIT_CAP, elapsed / max(1.0, proc.seconds)) * proc.seconds
            if total <= 0:
                return 0.0
            return max(0.0, min(0.99, done / total))

    def current(self) -> dict | None:
        """The process with the most remaining work, for the step label."""
        with self._lock:
            best: dict | None = None
            best_remaining = -1.0
            for name, proc in self._processes.items():
                if proc.first_seen is None:
                    continue
                expected = max(proc.expected, proc.submitted + proc.cached)
                finished = proc.completed + proc.cached
                remaining = max(0.0, expected - finished) * proc.seconds
                if remaining > best_remaining:
                    best_remaining = remaining
                    best = {
                        "step": name.lower().replace("_", "-"),
                        "process": name,
                        "done": finished,
                        "total": int(round(expected)),
                    }
            return best

    def finish_events(self, now: float | None = None) -> list[dict]:
        """Close out processes whose done event never fired (e.g. fewer tasks)."""
        now = time.monotonic() if now is None else now
        events: list[dict] = []
        with self._lock:
            for name, proc in self._processes.items():
                if proc.first_seen is None or proc.done_logged:
                    continue
                proc.done_logged = True
                finished = proc.completed + proc.cached
                elapsed = now - proc.first_seen
                events.append({"process": name, "event": "done", "done": finished,
                               "total": finished, "elapsed": elapsed})
        return events


# ---- post-processing tracker -------------------------------------------

@dataclass
class PostprocessingProgress:
    """Sample-by-sample post-processing progress weighted by step time."""

    sample_kinds: list[str]
    profile: dict
    research_only: bool = False
    dragen: bool = False
    _index: int = 0
    _step: str | None = None
    _step_started: float | None = None
    _done_steps: list[str] = field(default_factory=list)
    step_seconds: dict[str, list[float]] = field(default_factory=dict)
    _lock: threading.RLock = field(default_factory=threading.RLock)

    def _weights(self, kind: str) -> dict[str, float]:
        weights = postprocessing_step_seconds(self.profile, kind)
        if not self.research_only:
            weights.pop("spliceai", None)
        if not self.dragen:
            weights.pop("cnv-rescue", None)
        return weights

    def start_sample(self, index: int) -> None:
        with self._lock:
            self._index = index
            self._step = None
            self._step_started = None
            self._done_steps = []

    def on_marker(self, step: str, state: str, elapsed: float | None = None, now: float | None = None) -> None:
        now = time.monotonic() if now is None else now
        with self._lock:
            if state == "start":
                self._step = step
                self._step_started = now
            elif state == "done":
                if elapsed is None and self._step == step and self._step_started is not None:
                    elapsed = now - self._step_started
                if elapsed is not None:
                    self.step_seconds.setdefault(step, []).append(float(elapsed))
                if step not in self._done_steps:
                    self._done_steps.append(step)
                if self._step == step:
                    self._step = None
                    self._step_started = None

    def fraction(self, now: float | None = None) -> float:
        now = time.monotonic() if now is None else now
        with self._lock:
            total_samples = max(1, len(self.sample_kinds))
            kind = self.sample_kinds[min(self._index, total_samples - 1)] if self.sample_kinds else ""
            weights = self._weights(kind)
            total = sum(weights.values()) or 1.0
            done = sum(weights.get(step, 0.0) for step in self._done_steps)
            if self._step and self._step_started is not None:
                expected = weights.get(self._step, 0.0)
                if expected > 0:
                    done += min(RUNNING_CREDIT_CAP, (now - self._step_started) / expected) * expected
            within = max(0.0, min(1.0, done / total))
            return max(0.0, min(1.0, (self._index + within) / total_samples))


# ---- ticker -------------------------------------------------------------

class Ticker:
    """Call ``fn`` every ``interval`` seconds on a daemon thread."""

    def __init__(self, fn, interval: float = 15.0):
        self._fn = fn
        self._interval = interval
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._loop, name="tertiary-progress", daemon=True)

    def _loop(self) -> None:
        while not self._stop.wait(self._interval):
            try:
                self._fn()
            except Exception:  # progress must never break the job
                pass

    def start(self) -> "Ticker":
        self._thread.start()
        return self

    def stop(self) -> None:
        self._stop.set()
        if self._thread.is_alive():
            self._thread.join(timeout=self._interval + 5)

    def __enter__(self):
        return self.start()

    def __exit__(self, *_exc):
        self.stop()
        return False


def eta_seconds(started: float, pct: float, now: float | None = None) -> int | None:
    """Remaining seconds from elapsed time and a time-calibrated percentage."""
    now = time.monotonic() if now is None else now
    if pct < 5 or pct >= 100:
        return None
    elapsed = now - started
    return int(elapsed * (100.0 - pct) / pct)
