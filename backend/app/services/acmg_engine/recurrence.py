"""Incremental WES recurrence table (plan §5.2.5 step 0, §5.2.11).

Counts in how many WES *patients* each variant appears, from every WES case's
08 review TSV (rare candidate variants only).  Each case's contribution is
stored separately so a re-run replaces it and a deleted case can be removed;
old WES cases can be added later as their tertiary analysis is re-run.

Patients: registered cases merge by MRN; otherwise the sample ID with its
caller suffix and ``_1``…``_N`` replicate suffix removed.  Sample IDs ending in
``T`` (test samples) are never added.  WGS cases use the in-house AF instead.
"""
from __future__ import annotations

import csv
import json
import re
import sqlite3
import time
from pathlib import Path

from . import parse as P

_CALLER_SUFFIX = re.compile(r"-(nckuh|dragen|inhouse)$", re.IGNORECASE)
_REPLICATE_SUFFIX = re.compile(r"_\d+$")
_TEST_SAMPLE = re.compile(r"\d+T$")


def base_sample_id(sample_id: str) -> str:
    return _CALLER_SUFFIX.sub("", str(sample_id or "").strip())


def is_test_sample(sample_id: str) -> bool:
    """``23WE0044T`` style IDs are test samples (TITAN ``26T00970`` is not)."""
    return bool(_TEST_SAMPLE.search(base_sample_id(sample_id)))


def patient_key(sample_id: str, mrn: str = "") -> str:
    mrn = str(mrn or "").strip()
    if mrn:
        return f"mrn:{mrn}"
    return "id:" + _REPLICATE_SUFFIX.sub("", base_sample_id(sample_id))


def review_variant_keys(review_tsv: Path) -> set[str]:
    """Distinct CHROM-POS-REF-ALT keys of a review TSV (mitochondrial rows skipped)."""
    keys: set[str] = set()
    with open(review_tsv, encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = P.text(row.get("CHROM")).upper().replace("CHR", "")
            if chrom in ("M", "MT") or not chrom:
                continue
            keys.add(P.variant_key(row))
    return keys


class WesRecurrence:
    """SQLite store; safe for concurrent readers and short writer transactions."""

    def __init__(self, path: Path, *, readonly: bool = False):
        self.path = Path(path)
        if readonly:
            self.conn = sqlite3.connect(f"file:{self.path}?mode=ro", uri=True, timeout=30,
                                        check_same_thread=False)
        else:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.conn = sqlite3.connect(self.path, timeout=30, check_same_thread=False)
            self.conn.execute("PRAGMA journal_mode=WAL")
            self.conn.executescript(
                """
                CREATE TABLE IF NOT EXISTS samples (
                    sample_id TEXT PRIMARY KEY, patient_key TEXT NOT NULL, signature TEXT NOT NULL,
                    variant_count INTEGER NOT NULL, updated_at TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS sample_variants (
                    sample_id TEXT NOT NULL, variant_key TEXT NOT NULL,
                    PRIMARY KEY (sample_id, variant_key)) WITHOUT ROWID;
                CREATE INDEX IF NOT EXISTS sample_variants_key ON sample_variants(variant_key);
                """
            )
        self._patients: int | None = None

    def close(self) -> None:
        self.conn.close()

    # writes
    def upsert_sample(self, sample_id: str, patient: str, signature: dict, variant_keys) -> int:
        keys = sorted(set(variant_keys))
        sig = json.dumps(signature, sort_keys=True)
        with self.conn:
            self.conn.execute("DELETE FROM sample_variants WHERE sample_id=?", (sample_id,))
            self.conn.executemany("INSERT INTO sample_variants VALUES (?,?)", [(sample_id, k) for k in keys])
            self.conn.execute(
                "INSERT INTO samples VALUES (?,?,?,?,?) ON CONFLICT(sample_id) DO UPDATE SET "
                "patient_key=excluded.patient_key, signature=excluded.signature, "
                "variant_count=excluded.variant_count, updated_at=excluded.updated_at",
                (sample_id, patient, sig, len(keys), time.strftime("%Y-%m-%dT%H:%M:%S%z")))
        self._patients = None
        return len(keys)

    def remove_sample(self, sample_id: str) -> bool:
        with self.conn:
            self.conn.execute("DELETE FROM sample_variants WHERE sample_id=?", (sample_id,))
            removed = self.conn.execute("DELETE FROM samples WHERE sample_id=?", (sample_id,)).rowcount
        self._patients = None
        return bool(removed)

    # reads
    def stored_signature(self, sample_id: str) -> dict | None:
        hit = self.conn.execute("SELECT signature FROM samples WHERE sample_id=?", (sample_id,)).fetchone()
        return json.loads(hit[0]) if hit else None

    def sample_ids(self) -> set[str]:
        return {r[0] for r in self.conn.execute("SELECT sample_id FROM samples")}

    def patient_count(self) -> int:
        if self._patients is None:
            self._patients = self.conn.execute("SELECT COUNT(DISTINCT patient_key) FROM samples").fetchone()[0]
        return self._patients

    def lookup(self, variant_key: str) -> tuple[int, int]:
        """(patients carrying the variant, patients in the table)."""
        hit = self.conn.execute(
            "SELECT COUNT(DISTINCT s.patient_key) FROM sample_variants v JOIN samples s USING (sample_id) "
            "WHERE v.variant_key=?", (variant_key,)).fetchone()
        return int(hit[0] or 0), self.patient_count()

    def signature(self) -> dict:
        row = self.conn.execute("SELECT COUNT(*), MAX(updated_at) FROM samples").fetchone()
        return {"path": str(self.path), "samples": row[0], "patients": self.patient_count(),
                "updated_at": row[1] or ""}


# ── case scanning (shared by the CLI and post-processing) ─────────────────
def _scoped(directory: Path, sample_id: str, name: str) -> Path:
    prefixed = directory / f"{sample_id}.{name}"
    return prefixed if prefixed.exists() else directory / name


def _load_json(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def load_roster() -> dict | None:
    """Uploaded case roster (LIS ID -> MRN) or None outside the app."""
    try:
        from .. import patient_list_store
        return patient_list_store.load_roster()
    except Exception:
        return None


def _roster_mrn(roster: dict | None, sample_id: str, source_id: str) -> str:
    if not roster:
        return ""
    try:
        from .. import patient_list_store
        entry, _ = patient_list_store.lookup_with_key(sample_id, source_id, roster=roster)
    except Exception:
        return ""
    return str((entry or {}).get("mrn") or "").strip()


def case_contribution(sample_dir: Path, sample_id: str, roster: dict | None = None) -> dict:
    """Describe one case's WES contribution or why it is skipped.

    The patient is the registered MRN, else the roster MRN, else the base ID.

    Returns ``{"status": "ok", "patient", "signature", "review"}`` or
    ``{"status": "skip", "reason"}``.
    """
    from ..test_types import normalize_test_type  # local import: keeps engine core standalone

    post = Path(sample_dir) / "08_postprocessing"
    if is_test_sample(sample_id):
        return {"status": "skip", "reason": "測試檢體（T）"}
    manifest = _load_json(_scoped(post, sample_id, "snv_indel.review.tsv.source.json"))
    source = _load_json(_scoped(post, sample_id, "pipeline_source.json"))
    test_type = normalize_test_type(str(manifest.get("test_type") or source.get("seq_type") or ""),
                                    sample_id=sample_id, default="")
    if test_type != "WES":
        return {"status": "skip", "reason": f"非 WES（{test_type or '未知'}）"}
    review = _scoped(post, sample_id, "snv_indel.review.tsv")
    if not review.is_file():
        return {"status": "skip", "reason": "沒有 review TSV"}
    meta = _load_json(_scoped(post, sample_id, "sample_metadata.json"))
    mrn = (str(meta.get("mrn") or "").strip()
           or _roster_mrn(roster, sample_id, str(source.get("source_sample_id") or "")))
    st = review.stat()
    return {
        "status": "ok",
        "patient": patient_key(sample_id, mrn),
        "signature": {"review_size": st.st_size, "review_mtime_ns": st.st_mtime_ns, "mrn": mrn},
        "review": review,
    }


def sync_cases(store: WesRecurrence, root: Path, sample_ids=None, *, force: bool = False,
               prune: bool = False, dry_run: bool = False, roster: dict | None = None,
               progress=None) -> dict:
    """Add/replace the given WES cases (all directories under ``root`` when None)."""
    root = Path(root)
    if roster is None:
        roster = load_roster()
    if sample_ids is None:
        sample_ids = sorted(p.name for p in root.iterdir()
                            if p.is_dir() and not p.name.startswith((".", "_")))
    summary = {"added": [], "updated": [], "unchanged": [], "skipped": {}, "removed": []}
    seen = set()
    for index, sid in enumerate(sample_ids, 1):
        if progress:
            progress(index, len(sample_ids), sid)
        info = case_contribution(root / sid, sid, roster)
        if info["status"] != "ok":
            summary["skipped"][sid] = info["reason"]
            if store.stored_signature(sid) is not None:
                if not dry_run:
                    store.remove_sample(sid)
                summary["removed"].append(sid)
            continue
        seen.add(sid)
        old = store.stored_signature(sid)
        if old == info["signature"] and not force:
            summary["unchanged"].append(sid)
            continue
        if not dry_run:
            store.upsert_sample(sid, info["patient"], info["signature"], review_variant_keys(info["review"]))
        summary["updated" if old else "added"].append(sid)
    if prune:
        for sid in sorted(store.sample_ids() - seen):
            if not dry_run:
                store.remove_sample(sid)
            summary["removed"].append(sid)
    summary["patients"] = store.patient_count()
    return summary
