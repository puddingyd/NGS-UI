"""Run the rule engine over one case and store the results (plan §5.2, step 4).

Input is the case's 08 review TSV (raw 03_acmg rows already merged with the
post-processing overlay), so every candidate shown on the main screen gets an
engine result.  Gene-search variants outside the review TSV are evaluated on
demand by the loader with :func:`evaluate_rows`.

Results go to ``{LIS_ID}.snv_acmg_engine.sqlite`` (one row per genomic
variant, keyed like the UI variant id ``chr1-123-A-G``) with a manifest.  The
manifest's *core* signature (review TSV, reference tables, engine code and
rulebook) decides whether the file is current.  The WES recurrence table and
OMIM inheritance are recorded but do not invalidate the file: they change
whenever any case or OMIM entry is updated, and a later backfill picks the
changes up.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sqlite3
import time
import zlib
from functools import lru_cache
from pathlib import Path

from . import parse as P
from .engine import Context, evaluate_variant, load_rulebook
from .ref_store import SqliteReferences, load_ba1_exceptions
from .refs import References

RESULTS_NAME = "snv_acmg_engine.sqlite"
FORMAT_VERSION = 4          # 2: CaseRecurrence delegates the missense window / gene count lookups
                            # 3: CaseRecurrence delegates the collagen Gly-X-Y lookup
                            # 4: CaseRecurrence answers PS4 lookups (in-house probands, ClinVar cases)
_PACKAGE_DIR = Path(__file__).resolve().parent


_RESULT_CODE_FILES = ("engine.py", "parse.py", "refs.py", "ref_store.py")


@lru_cache(maxsize=1)
def code_signature() -> str:
    """Hash of the files that decide results: engine, parsers, reference lookups,
    rulebook and BA1 exception list.

    Tooling in the package (this runner's storage/attach code, recurrence
    maintenance, validation) is left out so editing it does not mark every
    stored result stale.  A change to ``evaluate_rows``, ``CaseRecurrence`` or
    ``omim_moi`` here that alters results must bump ``FORMAT_VERSION``.
    """
    digest = hashlib.sha256()
    for path in sorted(_PACKAGE_DIR.iterdir()):
        if path.is_file() and (path.name in _RESULT_CODE_FILES or path.suffix in (".json", ".tsv")):
            digest.update(path.name.encode())
            digest.update(path.read_bytes())
    return digest.hexdigest()[:16]


def ui_variant_id(row: dict) -> str:
    """Same id as ``adapters.snv_tsv._row_to_variant`` (raw CHROM, integer POS)."""
    pos = P.text(row.get("POS"))
    try:
        pos = str(int(float(pos)))
    except ValueError:
        pass
    return f"{P.text(row.get('CHROM'))}-{pos}-{P.text(row.get('REF'))}-{P.text(row.get('ALT'))}"


def group_variant_rows(rows) -> dict[str, list[dict]]:
    """Group transcript rows by genomic variant (mitochondrial rows skipped)."""
    groups: dict[str, list[dict]] = {}
    for row in rows:
        chrom = P.text(row.get("CHROM")).upper().replace("CHR", "")
        if not chrom or chrom in ("M", "MT") or not P.text(row.get("POS")):
            continue
        groups.setdefault(ui_variant_id(row), []).append(row)
    return groups


# ── PS4 A: in-house probands ───────────────────────────────────────────────
def _sample_mrn(sample_id: str) -> str:
    try:
        from .. import sample_layout
        meta = json.loads(sample_layout.state_file(sample_id, "sample_metadata.json").read_text(encoding="utf-8"))
        return str(meta.get("mrn") or "").strip() if isinstance(meta, dict) else ""
    except Exception:  # noqa: BLE001
        return ""


class InhouseProbands:
    """Variants marked Causative (status 1) in other patients' cases.

    Read from the observation registry in ``manual_acmg.sqlite`` (Causative /
    Other statuses that reviewers save).  Samples of the same patient (same
    MRN, or the same base sample ID) count once, and the case being scored is
    never counted."""

    def __init__(self, by_variant: dict[str, set], patient_of=None, signature: str = ""):
        self.by_variant = by_variant
        self._patient_of = patient_of
        self._patients: dict[str, str] = {}
        self._signature = signature

    @classmethod
    def from_db(cls, db: Path, patient_of=None) -> "InhouseProbands | None":
        db = Path(db)
        if not db.is_file():
            return None
        by_variant: dict[str, set] = {}
        try:
            conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
            try:
                rows = conn.execute("SELECT variant_id, sample_id FROM observations "
                                    "WHERE status='1' AND genome_build IN ('hg38', 'GRCh38')").fetchall()
            finally:
                conn.close()
        except sqlite3.Error:
            return None
        for vid, sid in rows:
            by_variant.setdefault(str(vid), set()).add(str(sid))
        return cls(by_variant, patient_of, signature=_db_signature(db))

    def signature(self) -> str:
        return self._signature

    def patient(self, sample_id: str) -> str:
        if sample_id not in self._patients:
            from .recurrence import patient_key
            fn = self._patient_of or _sample_mrn
            self._patients[sample_id] = patient_key(sample_id, fn(sample_id))
        return self._patients[sample_id]

    def count(self, variant_key: str, sample_id: str = "") -> int:
        try:
            from ..manual_acmg import normalize_variant_id
            key = normalize_variant_id(variant_key)
        except Exception:  # noqa: BLE001
            key = variant_key
        samples = self.by_variant.get(key)
        if not samples:
            return 0
        own = self.patient(sample_id) if sample_id else None
        return len({self.patient(s) for s in samples if s != sample_id} - {own})


def _db_signature(db: Path) -> str:
    parts = []
    for path in (Path(db), Path(f"{db}-wal")):
        try:
            st = path.stat()
            parts.append(f"{st.st_size}:{st.st_mtime_ns}")
        except OSError:
            parts.append("-")
    return "|".join(parts)


def open_probands(db: Path | None = None):
    if db is None:
        try:
            from ... import config
            db = Path(config.MANUAL_ACMG_DB)
        except Exception:  # noqa: BLE001
            return None
    return InhouseProbands.from_db(db)


# ── per-case context ───────────────────────────────────────────────────────
class CaseRecurrence(References):
    """Delegates to shared references; WES recurrence excludes the case's own patient."""

    def __init__(self, refs: References, store=None, sample_id: str = "", patient: str = "",
                 probands: "InhouseProbands | None" = None):
        self._refs = refs
        self._store = store
        self._sample_id = sample_id
        self._patient = patient
        self._probands = probands

    @property
    def available(self):  # type: ignore[override]
        tables = set(self._refs.available)
        if self._store is not None:
            tables.add("recurrence")
        return frozenset(tables)

    @property
    def version(self):  # type: ignore[override]
        info = dict(getattr(self._refs, "version", {}) or {})
        if self._store is not None:
            info["wes_recurrence"] = self._store.signature()
        return info

    def __getattr__(self, name):
        return getattr(self._refs, name)

    def lof_mechanism(self, gene):
        return self._refs.lof_mechanism(gene)

    def truncating_downstream(self, gene, protein_pos):
        return self._refs.truncating_downstream(gene, protein_pos)

    def clinvar_missense(self, gene, protein_pos):
        return self._refs.clinvar_missense(gene, protein_pos)

    def clinvar_missense_window(self, gene, start, end):
        return self._refs.clinvar_missense_window(gene, start, end)

    def benign_missense_window(self, gene, start, end):
        return self._refs.benign_missense_window(gene, start, end)

    def missense_gene_counts(self, gene):
        return self._refs.missense_gene_counts(gene)

    def exon_skip(self, transcript, exon_number):
        return self._refs.exon_skip(transcript, exon_number)

    def ba1_exception(self, variant_key):
        return self._refs.ba1_exception(variant_key)

    def vcep(self, gene):
        return self._refs.vcep(gene)

    def collagen_gly(self, gene):
        return self._refs.collagen_gly(gene)

    def wes_recurrence(self, variant_key):
        if self._store is None:
            return None
        return self._store.lookup(variant_key, exclude_patient=self._patient or None)

    def inhouse_probands(self, variant_key):
        if self._probands is None:
            return self._refs.inhouse_probands(variant_key)
        return self._probands.count(variant_key, self._sample_id)

    def clinvar_affected(self, variation_id):
        return self._refs.clinvar_affected(variation_id)

    def clinvar_biallelic(self, variation_id):
        return self._refs.clinvar_biallelic(variation_id)


def omim_moi(gene: str) -> set | None:
    """Gene-level OMIM inheritance codes, or None when unknown."""
    if not gene:
        return None
    try:
        from .. import omim_store
        row = omim_store.lookup(gene=gene)
    except Exception:
        return None
    codes = P.moi_codes((row or {}).get("Inheritance") or "")
    return codes or None


def omim_status() -> str:
    """"loaded" when OMIM inheritance can be read, else "unavailable".

    omim_store fails silently (missing openpyxl, unreadable workbook), which
    would quietly drop every recessive-gene rule; callers check this first.
    """
    try:
        from .. import omim_store
        return "loaded" if omim_store.is_loaded() else "unavailable"
    except Exception:
        return "unavailable"


# ── evaluation ─────────────────────────────────────────────────────────────
def evaluate_rows(groups: dict[str, list[dict]], *, test_type: str, refs: References,
                  moi_fn=omim_moi, rulebook: dict | None = None):
    """Yield ``(ui_variant_id, result)`` for each grouped variant."""
    rb = rulebook or load_rulebook()
    moi_cache: dict[str, set | None] = {}
    for vid, rows in groups.items():
        gene = P.text(rows[0].get("GENE"))
        if gene not in moi_cache:
            moi_cache[gene] = moi_fn(gene) if moi_fn else None
        ctx = Context(test_type=test_type, moi=moi_cache[gene], refs=refs)
        yield vid, evaluate_variant(rows, ctx, rb)


def summary(result: dict) -> dict:
    """Small per-variant fields attached to UI payloads."""
    return {
        "classification": result.get("classification") or "",
        "points": result.get("points"),
        "engine_points": result.get("engine_points"),
        "engine_classification": result.get("engine_classification") or "",
        "criteria_text": result.get("criteria_text") or "",
        "clinvar_action": (result.get("clinvar_rule") or {}).get("action") or "none",
        "clinvar_source": (result.get("clinvar_rule") or {}).get("source") or "",
        "flags": [f.get("code") for f in result.get("flags") or [] if isinstance(f, dict)],
        "ps4_inhouse": result.get("ps4_inhouse"),
        "rank_points": result.get("rank_points"),
        "version": result.get("version") or "",
    }


def evidence_view(result: dict, status: str) -> dict:
    """Per-criterion reasons of one variant for the manual ACMG modal."""
    rule = result.get("clinvar_rule") or {}
    own = own_call(summary(result))
    return {
        "status": status,
        "version": result.get("version") or "",
        "classification": own[0], "points": own[1],
        "rank_classification": result.get("classification") or "", "rank_points": result.get("points"),
        "clinvar_action": rule.get("action") or "none", "clinvar_source": rule.get("source") or "",
        "clinvar_category": rule.get("category") or "", "clinvar_stars": rule.get("stars") or 0,
        "vcep": {k: (result.get("vcep") or {}).get(k) for k in ("vcep", "svi", "version", "url")} if result.get("vcep") else None,
        "evidence": [{k: ev.get(k) for k in ("code", "status", "strength", "points", "reason")}
                     for ev in result.get("evidence") or [] if isinstance(ev, dict) and ev.get("code")],
        "flags": [{"code": f.get("code"), "message": f.get("message")}
                  for f in result.get("flags") or [] if isinstance(f, dict)],
    }


def evidence_detail(sample_id: str, variant_id: str, *, rows_loader=None) -> dict | None:
    """Stored reasons for one variant (``current`` / ``stale``); variants outside the
    review list are evaluated on demand from ``rows_loader(variant_id)``."""
    case = case_inputs(sample_id)
    path = case["existing"]
    stored = lookup_many(path, [variant_id], detail=True).get(variant_id) if Path(path).is_file() else None
    refs = shared_references()
    if stored:
        status = "current" if is_current(path, case["review"], refs) else "stale"
        return evidence_view(stored, status)
    rows = rows_loader(variant_id) if rows_loader else None
    if not rows:
        return None
    groups = {vid: grp for vid, grp in group_variant_rows(rows).items() if vid == variant_id}
    if not groups:
        return None
    recurrence = shared_recurrence() if case["test_type"] == "WES" else None
    case_refs = CaseRecurrence(refs, recurrence, sample_id,
                               case_patient(recurrence, sample_id, case["mrn"]), shared_probands())
    for _, result in evaluate_rows(groups, test_type=case["test_type"], refs=case_refs):
        return evidence_view(result, "on_demand")
    return None


def expert_adopted(engine: dict, vcep_class: str = "") -> bool:
    """The stored class was adopted from a ClinGen expert panel (VCEP) call.

    Summaries written before ``clinvar_source`` was kept fall back to the
    variant's own ``CLINGEN_VCEP_CLASS`` (the engine prefers it over ClinVar).
    """
    engine = engine or {}
    if (engine.get("clinvar_action") or "none") != "adopt":
        return False
    source = engine.get("clinvar_source")
    if source:
        return source == "ClinGen VCEP"
    return bool(str(vcep_class or "").strip()) and "conflict" not in str(vcep_class).lower()


def own_call(engine: dict) -> tuple[str, object]:
    """(class, points) from the rule's own criteria, before the ClinVar adopt/floor.

    ``classification``/``points`` already carry the ClinVar ≥2★ (or ClinGen
    VCEP) call; the platform ranks with those but shows this one.  Results
    stored before ``engine_classification`` was kept are derived the same way
    as the engine: BA1 met → Benign, otherwise the Tavtigian bins.
    """
    engine = engine or {}
    points = engine.get("engine_points")
    if points is None:
        return engine.get("classification") or "", engine.get("points")
    name = engine.get("engine_classification")
    if not name:
        from .engine import _class_for
        criteria = {c.strip().split("_")[0] for c in str(engine.get("criteria_text") or "").split(",") if c.strip()}
        name = "Benign" if "BA1" in criteria else _class_for(int(points))
    return name, points


# ── storage ────────────────────────────────────────────────────────────────
def _review_signature(review_tsv: Path) -> dict:
    """Size and mtime only: the staged case is copied (mtime kept) to another path."""
    try:
        st = Path(review_tsv).stat()
    except OSError:
        return {"missing": True}
    return {"size": st.st_size, "mtime_ns": st.st_mtime_ns}


def core_signature(review_tsv: Path, refs: References) -> dict:
    version = dict(getattr(refs, "version", {}) or {})
    version.pop("wes_recurrence", None)
    return {
        "format": FORMAT_VERSION,
        "code": code_signature(),
        "review": _review_signature(review_tsv),
        "refs": version,
    }


_ZDICT_SAMPLE = 40          # results used to build the per-file compression dictionary
_ZDICT_BYTES = 32 * 1024


def _compressor(zdict: bytes):
    def compress(data: bytes) -> bytes:
        obj = zlib.compressobj(6, zdict=zdict) if zdict else zlib.compressobj(6)
        return obj.compress(data) + obj.flush()
    return compress


def _decompress(blob: bytes, zdict: bytes) -> bytes:
    obj = zlib.decompressobj(zdict=zdict) if zdict else zlib.decompressobj()
    return obj.decompress(blob) + obj.flush()


def write_results(path: Path, results, manifest: dict) -> int:
    """Atomically write results; details share one zlib dictionary per file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.unlink(missing_ok=True)
    conn = sqlite3.connect(tmp)
    count = 0
    try:
        conn.executescript(
            """
            PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;
            CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE zdict (id INTEGER PRIMARY KEY, data BLOB NOT NULL);
            CREATE TABLE results (
                variant_id TEXT PRIMARY KEY, gene TEXT, transcript TEXT,
                classification TEXT, points INTEGER, criteria_text TEXT,
                summary TEXT NOT NULL, detail BLOB NOT NULL) WITHOUT ROWID;
            """
        )
        iterator = iter(results)
        head = []
        for item in iterator:
            head.append(item)
            if len(head) >= _ZDICT_SAMPLE:
                break
        encoded = lambda result: json.dumps(result, ensure_ascii=False).encode("utf-8")  # noqa: E731
        zdict = b"".join(encoded(r) for _vid, r in head)[-_ZDICT_BYTES:]
        conn.execute("INSERT INTO zdict VALUES (1, ?)", (zdict,))
        compress = _compressor(zdict)

        def rows():
            for source in (head, iterator):
                for vid, result in source:
                    yield (vid, result.get("gene") or "", result.get("transcript") or "",
                           result.get("classification") or "", result.get("points"),
                           result.get("criteria_text") or "",
                           json.dumps(summary(result), ensure_ascii=False), compress(encoded(result)))

        batch = []
        for item in rows():
            batch.append(item)
            if len(batch) >= 5000:
                conn.executemany("INSERT OR REPLACE INTO results VALUES (?,?,?,?,?,?,?,?)", batch)
                count += len(batch)
                batch.clear()
        if batch:
            conn.executemany("INSERT OR REPLACE INTO results VALUES (?,?,?,?,?,?,?,?)", batch)
            count += len(batch)
        manifest = dict(manifest, variants=count, written_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
        conn.execute("INSERT INTO meta VALUES ('manifest', ?)", (json.dumps(manifest, ensure_ascii=False),))
        conn.commit()
    finally:
        conn.close()
    os.replace(tmp, path)
    return count


def read_manifest(path: Path) -> dict | None:
    path = Path(path)
    if not path.is_file():
        return None
    try:
        conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        try:
            row = conn.execute("SELECT value FROM meta WHERE key='manifest'").fetchone()
        finally:
            conn.close()
        return json.loads(row[0]) if row else None
    except (sqlite3.Error, ValueError):
        return None


def is_current(path: Path, review_tsv: Path, refs: References) -> bool:
    manifest = read_manifest(path)
    if not manifest or manifest.get("core") != core_signature(review_tsv, refs):
        return False
    # Built while OMIM could not be read: recessive rules were skipped.
    return manifest.get("omim") != "unavailable"


def lookup_many(path: Path, variant_ids, *, detail: bool = False) -> dict[str, dict]:
    """``{variant_id: summary}`` (or full result when ``detail``) for stored variants."""
    ids = [v for v in dict.fromkeys(variant_ids) if v]
    path = Path(path)
    if not ids or not path.is_file():
        return {}
    column = "detail" if detail else "summary"
    out: dict[str, dict] = {}
    try:
        conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        try:
            zdict = b""
            if detail:
                hit = conn.execute("SELECT data FROM zdict WHERE id=1").fetchone()
                zdict = bytes(hit[0]) if hit else b""
            for start in range(0, len(ids), 900):
                chunk = ids[start:start + 900]
                marks = ",".join("?" * len(chunk))
                for vid, value in conn.execute(
                        f"SELECT variant_id, {column} FROM results WHERE variant_id IN ({marks})", chunk):
                    text = _decompress(value, zdict).decode("utf-8") if detail else value
                    out[vid] = json.loads(text)
        finally:
            conn.close()
    except (sqlite3.Error, ValueError, zlib.error):
        return {}
    return out


def iter_results(path: Path, *, detail: bool = True):
    """Yield ``(variant_id, result)`` for every stored variant (full result when ``detail``)."""
    path = Path(path)
    if not path.is_file():
        return
    conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        zdict = b""
        if detail:
            hit = conn.execute("SELECT data FROM zdict WHERE id=1").fetchone()
            zdict = bytes(hit[0]) if hit else b""
        column = "detail" if detail else "summary"
        for vid, value in conn.execute(f"SELECT variant_id, {column} FROM results"):
            yield vid, json.loads(_decompress(value, zdict).decode("utf-8") if detail else value)
    finally:
        conn.close()


def _load_json(path: Path) -> dict:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def case_inputs(sample_id: str) -> dict:
    """Review TSV, results paths, effective test type (WES/WGS) and MRN of a published case."""
    from .. import sample_layout, test_types

    manifest = _load_json(sample_layout.review_manifest(sample_id))
    meta = _load_json(sample_layout.state_file(sample_id, "sample_metadata.json"))
    source = _load_json(sample_layout.state_file(sample_id, "pipeline_source.json"))
    identity = str(meta.get("lis_id") or meta.get("sample_id") or sample_id)
    test_type = test_types.normalize_test_type(
        meta.get("test_type") or manifest.get("test_type") or source.get("seq_type") or "",
        sample_id=identity, default="WES")
    return {
        "review": sample_layout.review_tsv(sample_id),
        "out": sample_layout.state_file(sample_id, RESULTS_NAME, for_write=True),
        "existing": sample_layout.state_file(sample_id, RESULTS_NAME),
        "test_type": "WGS" if test_types.is_wgs_type(test_type) else "WES",
        "mrn": str(meta.get("mrn") or "").strip(),
    }


# ── one case ───────────────────────────────────────────────────────────────
def open_references(refs_dir: Path | None = None) -> SqliteReferences:
    from .ref_store import REFS_DB_NAME, default_refs_dir
    refs_dir = Path(refs_dir) if refs_dir else default_refs_dir()
    return SqliteReferences(refs_dir / REFS_DB_NAME, ba1_exceptions=load_ba1_exceptions())


def open_recurrence(db: Path | None = None):
    from .recurrence import WesRecurrence
    if db is None:
        try:
            from ... import config
            db = Path(config.ACMG_WES_RECURRENCE_DB)
        except Exception:
            return None
    db = Path(db)
    if not db.is_file():
        return None
    return WesRecurrence(db, readonly=True)


def case_patient(store, sample_id: str, mrn: str = "") -> str:
    """The case's own patient key in the recurrence table (stored row wins)."""
    from .recurrence import patient_key
    if store is not None:
        stored = store.patient_of(sample_id)
        if stored:
            return stored
    return patient_key(sample_id, mrn)


def build_case(*, review_tsv: Path, out: Path, sample_id: str, test_type: str,
               refs: SqliteReferences | None = None, recurrence=None, mrn: str = "",
               moi_fn=omim_moi, probands=None) -> dict:
    """Evaluate every variant in ``review_tsv`` and write ``out``; returns stats."""
    started = time.time()
    refs = refs if refs is not None else open_references()
    test_type = (test_type or "WES").upper()
    probands = open_probands() if probands is None else probands
    case_refs = CaseRecurrence(refs, recurrence if test_type == "WES" else None,
                               sample_id, case_patient(recurrence, sample_id, mrn), probands)
    with open(review_tsv, encoding="utf-8", errors="replace", newline="") as handle:
        groups = group_variant_rows(csv.DictReader(handle, delimiter="\t"))
    classes: dict[str, int] = {}

    def _results():
        for vid, result in evaluate_rows(groups, test_type=test_type, refs=case_refs, moi_fn=moi_fn):
            classes[result["classification"]] = classes.get(result["classification"], 0) + 1
            yield vid, result

    omim = omim_status() if moi_fn is omim_moi else ("custom" if moi_fn else "none")
    manifest = {
        "sample_id": sample_id,
        "test_type": test_type,
        "core": core_signature(review_tsv, refs),
        "recurrence": recurrence.signature() if recurrence is not None and test_type == "WES" else None,
        "reference_tables": sorted(refs.available),
        "omim": omim,
        "probands": probands.signature() if probands is not None else None,
    }
    count = write_results(out, _results(), manifest)
    return {"sample_id": sample_id, "variants": count, "classes": classes, "omim": omim,
            "seconds": round(time.time() - started, 1), "reference_tables": sorted(refs.available)}


# ── web side: attach results to UI payloads ────────────────────────────────
_SHARED: dict = {}


def _shared(kind: str, path: Path | None, factory):
    """Cache one opened resource per file signature (reopened when the file changes)."""
    try:
        st = Path(path).stat() if path else None
        sig = (str(path), st.st_size, st.st_mtime_ns) if st else (str(path), None)
    except OSError:
        sig = (str(path), None)
    hit = _SHARED.get(kind)
    if hit and hit[0] == sig:
        return hit[1]
    value = factory()
    _SHARED[kind] = (sig, value)
    return value


def shared_references() -> SqliteReferences:
    from .ref_store import REFS_DB_NAME, default_refs_dir
    path = default_refs_dir() / REFS_DB_NAME
    return _shared("refs", path, open_references)


def shared_recurrence():
    try:
        from ... import config
        path = Path(config.ACMG_WES_RECURRENCE_DB)
    except Exception:
        return None
    return _shared("recurrence", path, lambda: open_recurrence(path))


def shared_probands():
    try:
        from ... import config
        db = Path(config.MANUAL_ACMG_DB)
    except Exception:  # noqa: BLE001
        return None
    sig = _db_signature(db)
    hit = _SHARED.get("probands")
    if hit and hit[0] == sig:
        return hit[1]
    value = InhouseProbands.from_db(db)
    _SHARED["probands"] = (sig, value)
    return value


def _review_rows(review_tsv: Path, wanted: set) -> list[dict]:
    """Rows of the review TSV whose UI variant id is in ``wanted``."""
    out = []
    try:
        with open(review_tsv, encoding="utf-8", errors="replace", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if ui_variant_id(row) in wanted:
                    out.append(row)
    except OSError:
        return []
    return out


def attach(variants: dict, *, results_path: Path, review_tsv: Path, rows=None,
           test_type: str = "WES", sample_id: str = "", mrn: str = "") -> dict:
    """Add ``variant["engine_acmg"]`` (summary + status) where missing.

    Stored results come from the case's results file (status ``current`` or
    ``stale``); with ``rows`` (gene search) the remaining variants are
    evaluated on demand (status ``on_demand``).  Returns counts per status.
    """
    counts: dict[str, int] = {}
    missing = [vid for vid, v in variants.items() if isinstance(v, dict) and "engine_acmg" not in v]
    if not missing:
        return counts
    refs = shared_references()
    probands = shared_probands()
    test_type = "WGS" if str(test_type).upper() in ("WGS", "TITAN-WGS") else "WES"
    _case_refs = []

    def case_refs():
        if not _case_refs:
            recurrence = shared_recurrence() if test_type == "WES" else None
            _case_refs.append(CaseRecurrence(refs, recurrence, sample_id,
                                             case_patient(recurrence, sample_id, mrn), probands))
        return _case_refs[0]

    stored = lookup_many(results_path, missing)
    if stored:
        status = "current" if is_current(results_path, review_tsv, refs) else "stale"
        for vid, value in stored.items():
            variants[vid]["engine_acmg"] = dict(value, status=status)
            counts[status] = counts.get(status, 0) + 1
        # PS4 A changes whenever a reviewer marks the variant Causative in another
        # case: re-evaluate just those variants instead of waiting for a rebuild.
        if probands is not None:
            changed = {vid for vid, value in stored.items()
                       if value.get("ps4_inhouse") is not None
                       and probands.count(vid, sample_id) != value.get("ps4_inhouse")}
            if changed:
                groups = group_variant_rows(_review_rows(review_tsv, changed))
                for vid, result in evaluate_rows(groups, test_type=test_type, refs=case_refs(),
                                                 moi_fn=lambda gene: omim_moi(gene)):
                    variants[vid]["engine_acmg"] = dict(summary(result), status=status, ps4_refreshed=True)
                    counts["ps4_refreshed"] = counts.get("ps4_refreshed", 0) + 1
    if rows:
        wanted = {vid for vid in missing if vid not in stored}
        groups = {vid: grp for vid, grp in group_variant_rows(rows).items() if vid in wanted}
        if groups:
            for vid, result in evaluate_rows(groups, test_type=test_type, refs=case_refs()):
                variants[vid]["engine_acmg"] = dict(summary(result), status="on_demand")
                counts["on_demand"] = counts.get("on_demand", 0) + 1
    return counts
