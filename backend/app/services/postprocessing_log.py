"""Compact, consistent lifecycle events; complete tool output is kept separately."""
from collections import deque
from datetime import datetime, timedelta, timezone
import re
import time


def event(sample, step, status, *, elapsed=None, detail=''):
    value = f'[postprocessing] sample={sample} step={step} status={status}'
    if elapsed is not None:
        value += f' elapsed={elapsed:.1f}s'
    return value + (f' | {detail}' if detail else '')


def timestamped(message, at=None):
    tz = timezone(timedelta(hours=8))
    now = datetime.fromtimestamp(at, tz) if at is not None else datetime.now(tz)
    return f'[{now:%Y-%m-%d %H:%M:%S}] {message}'


def print_event(message):
    print(timestamped(message), flush=True)


def failure_detail(log, fallback):
    """Prefer the actual tool failure to its later wrapper traceback."""
    lines = log.splitlines()
    for pattern in (r'FATAL:', r'(?:OSError|ValueError|FileNotFoundError):',
                    r'\[PRS] container temporary directory check failed:', r'(?:ERROR|Error|error):'):
        found = [line.strip() for line in lines if re.search(pattern, line)]
        if found:
            return found[-1][:1500]
    return fallback


def read_tail(path, size=50000):
    try:
        with path.open('rb') as stream:
            stream.seek(max(0, path.stat().st_size-size))
            return stream.read().decode('utf-8', errors='replace')
    except OSError:
        return ''


# "skipped" as a status word only: counters such as GeneBe's
# "'*'-skipped 0" are informational and must not mark the step skipped.
_SKIPPED_RE = re.compile(r"(?<![\w'*-])skipped\b(?!\s*\d)")


class Summary:
    def __init__(self, sample, step, logger=print_event):
        self.sample, self.step, self.logger = sample, step, logger
        self.started = {}
        self.tail = deque(maxlen=20)
        self.errors = deque(maxlen=10)
        self.warnings = set()

    def start(self, step=None, detail='', sample=None):
        step = step or self.step
        key = (sample or self.sample, step)
        if key not in self.started:
            self.started[key] = time.monotonic()
            self.logger(event(key[0], step, 'start', detail=detail))

    def done(self, step, *, elapsed=None, status='done', detail='', sample=None):
        key = (sample or self.sample, step)
        start = self.started.pop(key, time.monotonic())
        self.logger(event(key[0], step, status,
                          elapsed=elapsed if elapsed is not None else time.monotonic()-start, detail=detail))

    def _current(self):
        """(sample, step) that most recently started and is still open."""
        return next(reversed(self.started), (self.sample, self.step))

    def feed(self, line):
        text = line.strip()
        self.tail.append(text)
        # Batch post-processing tags per-sample steps with ``sample=SID``;
        # batch-wide steps belong to this run's label.
        marker = re.search(r'\[(?:post-processing-step|sample-step)]\s+([a-z0-9-]+)\s+(start|done)'
                           r'(?:\s+sample=(\S+))?', text)
        if marker:
            step, state, sample = marker.groups()
            if state == 'start': self.start(step, sample=sample)
            else:
                elapsed = re.search(r'elapsed=([\d.]+)s', text)
                self.done(step, elapsed=float(elapsed[1]) if elapsed else None, sample=sample)
            return
        gpn = re.search(r'\[perf] snv_review\.gpn_(start|done).*provider=(\S+)', text)
        if gpn:
            state, provider = gpn.groups()
            if state == 'start': self.start(provider)
            else:
                elapsed = re.search(r'elapsed=([\d.]+)s', text)
                status = re.search(r'\bstatus=(\S+)', text)
                self.done(provider, elapsed=float(elapsed[1]) if elapsed else None,
                          detail=status[0] if status else '')
            return
        if re.search(r'FATAL:|(?:ERROR|Error|OSError|ValueError|FileNotFoundError):', text):
            self.errors.append(text)
            self.logger(event(*self._current(), 'error', detail=text[:1500]))
        elif re.search(r'\bWARNING\b|\bwarning\b', text) or _SKIPPED_RE.search(text):
            # Repeated bind warnings can flood stderr; retain every raw line.
            if text not in self.warnings:
                self.warnings.add(text)
                status = 'warning' if re.search(r'\bWARNING\b|\bwarning\b', text) else 'skipped'
                self.logger(event(*self._current(), status, detail=text[:500]))
        elif re.search(r'\bbackfilled\b|candidate_rows=|integrated_PASS=|\[gpn-[^]]+] status=', text):
            self.logger(event(*self._current(), 'info', detail=text.split(' manifest=')[0][:500]))

    def finish(self, code):
        detail = failure_detail('\n'.join(self.errors or self.tail), f'exit={code}') if code else ''
        for sample, step in reversed(list(self.started)):
            self.done(step, status='failed' if code else 'done', detail=detail if code else '', sample=sample)
        return detail
