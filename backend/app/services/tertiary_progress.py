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
    # With the pipeline's Pangolin score cache (2026-11-21) nearly every
    # splice candidate is a lookup: ~30-80 s per sample instead of ~8 min.
    "PANGOLIN_SCORE": 60.0,
    "PARSE_CSQ": 960.0,
    "ACMG_CLASSIFY": 480.0,
    "PGX_STELLARPGX": 60.0,
    "PGX_HLA_EXTRACT": 60.0,
    "PGX_OPTITYPE": 480.0,
    "PGX_PHARMCAT": 240.0,
    "PGX_PARSE": 120.0,
}
UNKNOWN_PROCESS_SECONDS = 60.0

# WGS differs mostly in size: ~6M SNV rows per sample instead of ~0.1M.
# Measured on a 9-sample DRAGEN WGS Research-only batch (2026-11-22); used
# until the profile has learned the WGS kind.
DEFAULT_PROCESS_SECONDS_BY_SEQ: dict[str, dict[str, float]] = {
    "WGS": {
        "VEP_ANNOTATE": 3000.0,
        "PANGOLIN_SCORE": 70.0,
        "PARSE_CSQ": 510.0,
        "ACMG_CLASSIFY": 62.0,
        "PGX_OPTITYPE": 1100.0,
    },
}


def _seq_of(kind: str) -> str:
    return (kind or "").partition(":")[2].upper()


def default_process_seconds(name: str, seq_type: str = "") -> float:
    by_seq = DEFAULT_PROCESS_SECONDS_BY_SEQ.get((seq_type or "").upper(), {})
    if name in by_seq:
        return by_seq[name]
    return DEFAULT_PROCESS_SECONDS.get(name, UNKNOWN_PROCESS_SECONDS)

# Processes that only exist for one pipeline type / sequencing type.  Used to
# keep the default (no-profile) denominator close to what will really run.
_DRAGEN_ONLY = re.compile(r"(?:_DRAGEN$|^ADD_DRAGEN_TAG$)")
_NCKUH_ONLY = re.compile(r"(?:_NCKUH(?:_WES|_WGS)?$|^ADD_CALLERS_TAG$|^CNVKIT_TO_BED$)")
_GENERIC_CNV = {"ANNOTSV_CNV_NCKUH"}  # superseded by the _WES/_WGS variants

# Wall time of a phase or a batch-wide post-processing step is modelled as
# ``fixed + per_sample * n``: each resource has a one-off cost (start VEP,
# open/scan a database) shared by every sample of the batch, plus a
# per-sample part. Defaults are fitted from two NCKUH WES Research-only runs
# (1 sample and 36 samples, 2026-11-20); the learned profile replaces them.
BATCH_STEPS = ("clinvar-latest", "genebe", "spliceai", "giab-strata", "inhouse-af",
               "mane-refseq", "litvar2", "gpn-prefetch")
SAMPLE_STEPS = ("snv-overlay", "review-tsv", "gene-index", "acmg-engine", "cnv-rescue")
DEFAULT_LINEAR: dict[str, tuple[float, float]] = {
    "phase:nextflow": (114.0, 473.0),
    "phase:postprocessing": (650.0, 60.0),
    "phase:promote": (5.0, 10.0),
    "step:clinvar-latest": (30.0, 0.95),
    "step:genebe": (315.0, 15.0),
    "step:spliceai": (98.0, 10.0),
    "step:giab-strata": (10.0, 0.2),
    "step:inhouse-af": (31.0, 0.25),
    "step:mane-refseq": (1.0, 0.12),
    "step:litvar2": (20.0, 10.5),
    "step:gpn-prefetch": (141.0, 23.0),
}
# Per-sample steps (each sample on its own; run in parallel workers).
DEFAULT_SAMPLE_STEP_SECONDS: dict[str, float] = {
    "snv-overlay": 2.0,
    "review-tsv": 2.0,
    "gene-index": 1.0,
    "acmg-engine": 3.0,
    "cnv-rescue": 60.0,
}
# WGS starting point (same 2026-11-22 batch; post-processing figures assume
# the sorted GeneBe lookup and streaming GIAB pass). Replaced by learning.
DEFAULT_LINEAR_BY_SEQ: dict[str, dict[str, tuple[float, float]]] = {
    "WGS": {
        "phase:nextflow": (600.0, 3150.0),
        "phase:postprocessing": (1000.0, 330.0),
        "phase:promote": (5.0, 60.0),
        "step:clinvar-latest": (30.0, 39.0),
        "step:genebe": (120.0, 90.0),
        "step:spliceai": (98.0, 53.0),
        "step:giab-strata": (10.0, 20.0),
        "step:inhouse-af": (300.0, 30.0),
        "step:mane-refseq": (1.0, 15.0),
        "step:litvar2": (20.0, 15.0),
        "step:gpn-prefetch": (141.0, 40.0),
    },
}
DEFAULT_SAMPLE_STEP_SECONDS_BY_SEQ: dict[str, dict[str, float]] = {
    "WGS": {"snv-overlay": 120.0, "review-tsv": 180.0, "gene-index": 90.0, "acmg-engine": 60.0,
            "cnv-rescue": 60.0},
}


def default_linear(kind: str, key: str) -> tuple[float, float]:
    return DEFAULT_LINEAR_BY_SEQ.get(_seq_of(kind), {}).get(key, DEFAULT_LINEAR[key])


MAX_OBSERVATIONS = 12

# ---- small helpers ------------------------------------------------------

_ANSI_ESCAPE_RE = re.compile(r"\x1b(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
_STDOUT_TASK_RE = re.compile(
    r"\[(?P<hash>[0-9a-f]{2}/[0-9a-f]{6})\]\s+"
    r"(?P<event>Submitted|Cached|Re-submitted)\s+process\s+>\s+(?P<name>.+?)\s*$"
)
_DURATION_PART_RE = re.compile(r"(\d+(?:\.\d+)?)\s*(ms|d|h|m|s)")


def _process_path(name: str) -> str:
    """``ANNOTATE_SNV:SNV_ANNOTATE:VEP_ANNOTATE (12)`` → without the task tag."""
    text = _ANSI_ESCAPE_RE.sub("", str(name or "")).strip()
    return re.sub(r"\s*\(.*\)\s*$", "", text)


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


def fit_linear(observations, default: tuple[float, float]) -> tuple[float, float]:
    """Fit ``seconds = fixed + per_sample * n`` to ``[[n, seconds], ...]``.

    With observations at a single batch size the default's shape is kept and
    scaled to match; with two or more sizes it is a least-squares fit, with
    negative coefficients clamped to zero.
    """
    obs = []
    for item in observations or []:
        try:
            n, seconds = float(item[0]), float(item[1])
        except (TypeError, ValueError, IndexError):
            continue
        if n > 0 and seconds >= 0 and math.isfinite(seconds):
            obs.append((n, seconds))
    a0, b0 = default
    if not obs:
        return a0, b0
    sizes = {round(n, 3) for n, _s in obs}
    if len(sizes) < 2:
        n_mean = sum(n for n, _s in obs) / len(obs)
        s_mean = sum(sec for _n, sec in obs) / len(obs)
        base = a0 + b0 * n_mean
        scale = s_mean / base if base > 0 else 1.0
        return a0 * scale, b0 * scale
    count = len(obs)
    mean_n = sum(n for n, _s in obs) / count
    mean_s = sum(sec for _n, sec in obs) / count
    var = sum((n - mean_n) ** 2 for n, _s in obs)
    b = sum((n - mean_n) * (sec - mean_s) for n, sec in obs) / var if var > 0 else 0.0
    a = mean_s - b * mean_n
    if b < 0:
        return mean_s, 0.0
    if a < 0:
        return 0.0, sum(n * sec for n, sec in obs) / sum(n * n for n, _s in obs)
    return a, b


def _entry(profile: dict, kind: str) -> dict:
    return ((profile.get("kinds") or {}).get(kind) or {})


def estimate_seconds(profile: dict, kind: str, key: str, n: float) -> float:
    """Expected wall seconds of ``key`` (``phase:…`` / ``step:…``) for n samples."""
    if n <= 0:
        return 0.0
    observations = (_entry(profile, kind).get("observations") or {}).get(key)
    a, b = fit_linear(observations, default_linear(kind, key))
    return max(0.0, a + b * n)


def sample_step_seconds(profile: dict, kind: str) -> dict[str, float]:
    out = dict(DEFAULT_SAMPLE_STEP_SECONDS)
    out.update(DEFAULT_SAMPLE_STEP_SECONDS_BY_SEQ.get(_seq_of(kind), {}))
    for step, seconds in (_entry(profile, kind).get("sample_steps") or {}).items():
        try:
            if float(seconds) > 0:
                out[step] = float(seconds)
        except (TypeError, ValueError):
            continue
    return out


def process_seconds(profile: dict, kind: str, process: str) -> float:
    item = (_entry(profile, kind).get("processes") or {}).get(process) or {}
    try:
        value = float(item.get("seconds") or 0)
    except (TypeError, ValueError):
        value = 0.0
    return value if value > 0 else default_process_seconds(process, _seq_of(kind))


def noncached_share(profile: dict, kind: str, trace_rows: list[dict]) -> float:
    """Share of this run's Nextflow work that actually executed (by time).

    Counting tasks is misleading: dozens of quick cached tasks can hide that
    the expensive VEP/Pangolin chain re-ran. Each task is weighted by its
    process's expected duration instead.
    """
    ran = cached = 0.0
    for row in trace_rows:
        status = str(row.get("status") or "").upper()
        weight = process_seconds(profile, kind, simple_process_name(row.get("name", "")))
        if status == "CACHED":
            cached += weight
        elif status in {"COMPLETED", "FAILED"}:
            ran += weight
    total = ran + cached
    return ran / total if total > 0 else 0.0


def _observe(entry: dict, key: str, n: float, seconds: float | None) -> None:
    if seconds is None or n <= 0 or not math.isfinite(seconds) or seconds < 0:
        return
    observations = entry.setdefault("observations", {}).setdefault(key, [])
    observations.append([round(float(n), 3), round(float(seconds), 1)])
    del observations[:-MAX_OBSERVATIONS]


def update_profile(
    profile: dict,
    kind: str,
    *,
    trace_rows: list[dict] | None = None,
    sample_count: int = 0,
    nextflow_seconds: float | None = None,
    postprocessing_seconds: float | None = None,
    promote_seconds: float | None = None,
    step_observations: list[tuple[str, float, float]] | None = None,
    sample_step_seconds_seen: dict[str, list[float]] | None = None,
) -> dict:
    """Fold one successful job's measurements into ``profile`` (in place).

    Per-task process durations are smoothed (EMA). Phase and batch-step wall
    times are stored as ``[n, seconds]`` observations for the
    ``fixed + per_sample * n`` fit; the Nextflow phase uses the effective
    number of samples that really ran (time-weighted, cached work excluded).
    """
    entry = profile.setdefault("kinds", {}).setdefault(kind, {})
    share = noncached_share(profile, kind, trace_rows or []) if trace_rows else 1.0
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
        effective = sample_count * share
        # A (nearly) fully cached run says nothing about real Nextflow time.
        if nextflow_seconds and effective >= 0.5:
            _observe(entry, "phase:nextflow", effective, nextflow_seconds)
        if postprocessing_seconds and postprocessing_seconds >= 1:
            _observe(entry, "phase:postprocessing", sample_count, postprocessing_seconds)
        if promote_seconds is not None:
            _observe(entry, "phase:promote", sample_count, promote_seconds)
    for step, n, seconds in step_observations or []:
        if f"step:{step}" in DEFAULT_LINEAR:
            _observe(entry, f"step:{step}", n, seconds)
    if sample_step_seconds_seen:
        steps = entry.setdefault("sample_steps", {})
        for step, values in sample_step_seconds_seen.items():
            if values:
                steps[step] = round(_ema(steps.get(step), sum(values) / len(values)), 2)
    # Superseded per-sample averages from the first implementation.
    entry.pop("phases", None)
    entry.pop("postprocessing_batch_steps", None)
    entry["updated_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    entry["jobs"] = int(entry.get("jobs") or 0) + 1
    return profile


# ---- phase bands --------------------------------------------------------

def _kind_counts(kinds: list[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for kind in kinds:
        counts[kind] = counts.get(kind, 0) + 1
    return counts


def phase_bands(profile: dict, kinds: list[str], *, with_prs: bool) -> dict[str, tuple[float, float]]:
    """Split 0-100% into phases by expected wall time for this batch size.

    ``kinds`` has one entry per sample. 0-2% covers start-up and the PRS tail
    keeps 95-99%; Nextflow, post-processing and promotion share the rest in
    proportion to their ``fixed + per_sample * n`` estimates, so a single
    sample (fixed costs dominate) and a large batch get different splits.
    """
    counts = _kind_counts(kinds)
    est = {
        phase: sum(estimate_seconds(profile, kind, f"phase:{phase}", n) for kind, n in counts.items()) or 1.0
        for phase in ("nextflow", "postprocessing", "promote")
    }
    start = 2.0
    end = 95.0 if with_prs else 99.0
    total = sum(est.values())
    nf_end = start + (end - start) * est["nextflow"] / total
    pp_end = nf_end + (end - start) * est["postprocessing"] / total
    bands = {
        "pre": (0.0, start),
        "nextflow": (start, nf_end),
        "postprocessing": (nf_end, pp_end),
        "promote": (pp_end, end),
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
    full_name: str = ""
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
                    name: (float(item.get("seconds") or default_process_seconds(name, seq_type)),
                           float(item.get("tasks_per_sample") or 1.0))
                    for name, item in learned.items()
                }
            else:
                table = {
                    name: (default_process_seconds(name, seq_type), 1.0)
                    for name in DEFAULT_PROCESS_SECONDS
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
        full = _STDOUT_TASK_RE.search(_ANSI_ESCAPE_RE.sub("", line))["name"]
        now = time.monotonic() if now is None else now
        events: list[dict] = []
        with self._lock:
            if task_hash in self._finished:
                return []
            proc = self._proc(name)
            proc.full_name = proc.full_name or _process_path(full)
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
            proc.full_name = proc.full_name or _process_path(row.get("name", ""))
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
                # Progress is "how much of this batch's work is done": a task
                # Nextflow reuses from an earlier run (-resume) counts as done.
                # A re-run after a failure therefore starts where the cache
                # left off, for whatever subset/superset of samples is run.
                total += expected * proc.seconds
                done += (proc.completed + proc.failed + proc.cached) * proc.seconds
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

    def snapshot(self, *, finished: bool = False) -> list[dict]:
        """Per-process status rows, in pipeline order, for the live block.

        Mirrors Nextflow's console block: one row per process with
        ``done of total`` (cached included), running and failed counts.
        Processes expected for this batch but not started yet appear as
        pending; after the run (``finished``) those never started are dropped.
        """
        order = {name: index for index, name in enumerate(DEFAULT_PROCESS_SECONDS)}
        with self._lock:
            running: dict[str, int] = {}
            for name, _started in self._running.values():
                running[name] = running.get(name, 0) + 1
            rows = []
            for seen_index, (name, proc) in enumerate(self._processes.items()):
                if proc.first_seen is None and (finished or proc.expected <= 0):
                    continue
                done = proc.completed + proc.cached
                total = int(round(max(proc.expected, proc.submitted + proc.cached + proc.failed, done)))
                if proc.failed:
                    status = "failed"
                elif proc.first_seen is None:
                    status = "pending"
                elif done >= total and not running.get(name):
                    status = "done"
                else:
                    status = "running"
                rows.append({
                    "process": name,
                    "name": proc.full_name or name,
                    "done": done,
                    "total": total,
                    "cached": proc.cached,
                    "running": running.get(name, 0),
                    "failed": proc.failed,
                    "status": status,
                    "_order": (order.get(name, 1000 + seen_index)),
                })
            rows.sort(key=lambda row: row.pop("_order"))
            return rows

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
    """Batch post-processing progress weighted by expected step time.

    ``chunks`` holds the sample kinds of each chunk run by one orchestrator
    call. A batch-wide step (``[post-processing-step]``) costs
    ``fixed + per_sample * n`` for the chunk; per-sample steps
    (``[sample-step] … sample=SID``) run ``workers`` at a time, so the chunk
    spends ``n * seconds / workers`` of wall time on them.
    """

    chunks: list[list[str]]
    profile: dict
    research_only: bool = False
    dragen: bool = False
    workers: int = 1
    step_observations: list[tuple[str, float, float]] = field(default_factory=list)
    sample_seconds: dict[str, list[float]] = field(default_factory=dict)
    _chunk: int = 0
    _done: float = 0.0
    _running: dict[tuple[str, str], float] = field(default_factory=dict)
    _lock: threading.RLock = field(default_factory=threading.RLock)

    def _steps(self) -> tuple[list[str], list[str]]:
        batch = [s for s in BATCH_STEPS if self.research_only or s != "spliceai"]
        sample = [s for s in SAMPLE_STEPS if self.dragen or s != "cnv-rescue"]
        return batch, sample

    def _chunk_plan(self, index: int) -> dict[str, tuple[float, float]]:
        """{step: (wall seconds for the whole chunk, seconds per sample run)}."""
        kinds = self.chunks[index] if index < len(self.chunks) else []
        counts = _kind_counts(kinds)
        size = max(1, len(kinds))
        parallel = max(1, min(self.workers, size))
        batch, sample = self._steps()
        plan: dict[str, tuple[float, float]] = {}
        for step in batch:
            wall = sum(estimate_seconds(self.profile, kind, f"step:{step}", n) for kind, n in counts.items())
            plan[step] = (wall, wall)
        for step in sample:
            per = sum(sample_step_seconds(self.profile, kind)[step] * n for kind, n in counts.items()) / size
            plan[step] = (per * size / parallel, per)
        return plan

    def _size(self, index: int) -> int:
        return max(1, len(self.chunks[index]) if index < len(self.chunks) else 1)

    def _total(self) -> float:
        return sum(sum(w for w, _p in self._chunk_plan(i).values()) for i in range(len(self.chunks))) or 1.0

    def start_chunk(self, index: int) -> None:
        with self._lock:
            self._chunk = index
            self._running.clear()

    def finish_chunk(self, index: int) -> None:
        """Credit whatever a chunk did not report (skipped/absent steps)."""
        with self._lock:
            before = sum(sum(w for w, _p in self._chunk_plan(i).values()) for i in range(index + 1))
            self._done = max(self._done, before)
            self._running.clear()

    def on_marker(self, step: str, state: str, elapsed: float | None = None,
                  sample: str | None = None, now: float | None = None) -> None:
        now = time.monotonic() if now is None else now
        with self._lock:
            key = (step, sample or "")
            if state == "start":
                self._running[key] = now
                return
            started = self._running.pop(key, None)
            if elapsed is None and started is not None:
                elapsed = now - started
            wall, _per = self._chunk_plan(self._chunk).get(step, (0.0, 0.0))
            size = self._size(self._chunk)
            if sample:
                self._done += wall / size
                if elapsed is not None:
                    self.sample_seconds.setdefault(step, []).append(float(elapsed))
            else:
                self._done += wall
                if elapsed is not None:
                    self.step_observations.append((step, float(size), float(elapsed)))

    def fraction(self, now: float | None = None) -> float:
        now = time.monotonic() if now is None else now
        with self._lock:
            plan = self._chunk_plan(self._chunk)
            size = self._size(self._chunk)
            done = self._done
            for (step, sample), started in self._running.items():
                wall, per = plan.get(step, (0.0, 0.0))
                expected = per if sample else wall
                credit = wall / size if sample else wall
                if expected > 0:
                    done += min(RUNNING_CREDIT_CAP, (now - started) / expected) * credit
            return max(0.0, min(1.0, done / self._total()))


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


def eta_seconds(started: float, pct: float, now: float | None = None,
                base_pct: float = 0.0) -> int | None:
    """Remaining seconds from the progress rate since a baseline.

    ``started``/``base_pct`` mark where real work began: job start, or the
    moment cached Nextflow tasks were last credited (a resumed run jumps
    ahead instantly; that jump must not be read as speed).
    """
    now = time.monotonic() if now is None else now
    if pct < 5 or pct >= 100 or pct - base_pct < 2:
        return None
    elapsed = now - started
    return int(elapsed * (100.0 - pct) / (pct - base_pct))
