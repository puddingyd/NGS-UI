"""Build and read the rule engine's reference database (plan §5.2.11).

One SQLite file, ``acmg_refs.sqlite``, built atomically from:

* ClinGen dosage curation (haploinsufficiency score per gene);
* ClinVar ``variant_summary.txt.gz`` (GRCh38 P/LP truncating and missense
  variants, protein positions on the MANE transcript);
* MANE summary (RefSeq ↔ Ensembl MANE transcript per gene);
* MANE Ensembl GTF (per-exon CDS length for splice-site exon skipping).

Tables: ``lof_genes``, ``truncating``, ``missense`` (P/LP, with category),
``benign_missense`` (B/LB), ``missense_genes`` (per-gene ≥1★ missense counts),
``exons``, ``meta``.
:class:`SqliteReferences` implements the :class:`References` interface on top
of it and of the WES recurrence store.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import sqlite3
import time
from functools import lru_cache
from pathlib import Path

from . import parse as P
from .refs import References

SCHEMA_VERSION = 2
GENE_COUNT_MIN_STARS = 1      # missense_genes counts (PP2)
REFS_DB_NAME = "acmg_refs.sqlite"
BA1_EXCEPTIONS_PATH = Path(__file__).with_name("ba1_exceptions.tsv")
VCEP_SPECS_PATH = Path(__file__).with_name("vcep_specs.json")
COLLAGEN_PATH = Path(__file__).with_name("collagen_gly.json")

_NAME_RE = re.compile(r"^(?P<nm>N[MR]_\d+)(?:\.\d+)?\((?P<gene>[^)]+)\):(?P<c>c\.[^ ]+)(?:\s+\((?P<p>p\.[^)]+)\))?")
_P_RE = re.compile(r"^p\.(?P<ref>[A-Z][a-z]{2})(?P<pos>\d+)(?P<alt>[A-Z][a-z]{2}|fs|\*|=)")
_STARS = (
    ("practice guideline", 4),
    ("reviewed by expert panel", 3),
    ("criteria provided, multiple submitters, no conflicts", 2),
    ("criteria provided, conflicting classifications", 1),
    ("criteria provided, conflicting interpretations", 1),
    ("criteria provided, single submitter", 1),
)


def _open(path: Path):
    path = Path(path)
    if path.suffix in (".gz", ".bgz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, encoding="utf-8", errors="replace", newline="")


def file_signature(path: Path | None) -> dict:
    if not path:
        return {}
    path = Path(path)
    try:
        st = path.stat()
    except OSError:
        return {"path": str(path), "missing": True}
    return {"path": str(path), "size": st.st_size, "mtime": int(st.st_mtime)}


def review_stars(status: str) -> int:
    low = (status or "").strip().lower()
    for prefix, stars in _STARS:
        if low.startswith(prefix):
            return stars
    return 0


# ── source parsers ─────────────────────────────────────────────────────────
def parse_clingen_hi(path: Path) -> dict[str, str]:
    """Gene symbol -> haploinsufficiency score string ("3", "30", "40", ...)."""
    out: dict[str, str] = {}
    header: list[str] | None = None
    with _open(path) as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith("#"):
                cols = line.lstrip("#").split("\t")
                if len(cols) > 3 and any("Haploinsufficiency Score" in c for c in cols):
                    header = cols
                continue
            if not header or not line:
                continue
            fields = line.split("\t")
            row = dict(zip(header, fields))
            gene = (row.get("Gene Symbol") or "").strip()
            score = (row.get("Haploinsufficiency Score") or "").strip()
            if gene and score:
                out[gene] = score
    return out


def parse_mane_summary(path: Path) -> dict[str, dict]:
    """RefSeq NM (unversioned) -> {gene, ensembl, status} for MANE transcripts."""
    out: dict[str, dict] = {}
    with _open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row = {k.lstrip("#"): v for k, v in row.items() if k}
            nm = (row.get("RefSeq_nuc") or "").split(".")[0]
            if not nm:
                continue
            out[nm] = {
                "gene": row.get("symbol") or "",
                "ensembl": (row.get("Ensembl_nuc") or "").split(".")[0],
                "status": row.get("MANE_status") or "",
            }
    return out


def _variant_key(row: dict) -> str:
    chrom = (row.get("Chromosome") or "").strip()
    pos = (row.get("PositionVCF") or "").strip()
    ref = (row.get("ReferenceAlleleVCF") or "").strip()
    alt = (row.get("AlternateAlleleVCF") or "").strip()
    if not chrom or not pos.isdigit() or not ref or not alt or "na" in (ref.lower(), alt.lower()):
        return ""
    chrom = "M" if chrom.upper() == "MT" else chrom
    return f"chr{chrom}-{pos}-{ref}-{alt}"


def iter_clinvar_plp(path: Path):
    """Yield GRCh38 ClinVar P/LP records with parsed HGVS (non-conflicting)."""
    return iter_clinvar_records(path, ("P", "LP"))


def iter_clinvar_records(path: Path, categories=("P", "LP", "LB", "B")):
    """Yield GRCh38 ClinVar records of the given categories (non-conflicting)."""
    with _open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row = {k.lstrip("#"): v for k, v in row.items() if k}
            if (row.get("Assembly") or "") != "GRCh38":
                continue
            category = P.clinvar_category(row.get("ClinicalSignificance"))
            if category not in categories:
                continue
            m = _NAME_RE.match((row.get("Name") or "").strip())
            if not m:
                continue
            yield {
                "nm": m.group("nm"),
                "gene": m.group("gene").strip(),
                "hgvs_c": m.group("c"),
                "hgvs_p": m.group("p") or "",
                "stars": review_stars(row.get("ReviewStatus") or ""),
                "variant_key": _variant_key(row),
                "category": category,
            }


def classify_clinvar(record: dict) -> tuple[str, str, int | None]:
    """-> (kind, alt_aa, protein_pos) with kind in truncating/splice/missense/other."""
    pm = _P_RE.match(record["hgvs_p"]) if record["hgvs_p"] else None
    if pm:
        alt = pm.group("alt")
        pos = int(pm.group("pos"))
        if alt in ("Ter", "*", "fs"):
            return "truncating", alt, pos
        if alt != "=" and alt != pm.group("ref"):
            return "missense", alt, pos
        return "other", alt, pos
    if P.canonical_splice_hit(record["hgvs_c"]) and not re.search(r"(dup|ins|delins)", record["hgvs_c"]):
        return "splice", "", None
    return "other", "", None


def parse_mane_gtf(path: Path) -> dict[tuple[str, int], dict]:
    """(ENST unversioned, exon number) -> {cds_length, frameshift, fraction, last_exon}.

    Exons are numbered 5'→3' from their coordinates and strand (works whether
    or not the GTF quotes ``exon_number`` or has it at all); CDS segments are
    assigned to the exon that contains them.
    """
    exons: dict[str, list[tuple[int, int]]] = {}
    strands: dict[str, str] = {}
    cds: dict[str, list[tuple[int, int]]] = {}
    with _open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("exon", "CDS"):
                continue
            m = re.search(r'transcript_id\s+"?([^";\s]+)', fields[8])
            tx = m.group(1).split(".")[0] if m else ""
            if not tx.startswith("ENST"):
                continue
            span = (int(fields[3]), int(fields[4]))
            strands[tx] = fields[6]
            (exons if fields[2] == "exon" else cds).setdefault(tx, []).append(span)
    out: dict[tuple[str, int], dict] = {}
    for tx, spans in exons.items():
        ordered = sorted(set(spans), reverse=strands.get(tx) == "-")
        lengths = [0] * len(ordered)
        for c_start, c_end in cds.get(tx, ()):
            for i, (e_start, e_end) in enumerate(ordered):
                if e_start <= c_start and c_end <= e_end:
                    lengths[i] += c_end - c_start + 1
                    break
        total = sum(lengths)
        for i, length in enumerate(lengths):
            out[(tx, i + 1)] = {
                "cds_length": length,
                "frameshift": length % 3 != 0,
                "fraction": (length / total) if total else 0.0,
                "last_exon": i == len(ordered) - 1,
            }
    return out


# ── build ──────────────────────────────────────────────────────────────────
def build_reference_db(
    out_dir: Path,
    *,
    clinvar_summary: Path,
    mane_summary: Path,
    clingen_hi: Path,
    mane_gtf: Path | None = None,
    min_truncating: int = 2,
    clinvar_release: str = "",
) -> dict:
    """Build ``acmg_refs.sqlite`` in ``out_dir`` atomically; returns stats."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    final = out_dir / REFS_DB_NAME
    tmp = out_dir / f".{REFS_DB_NAME}.{os.getpid()}.tmp"
    if tmp.exists():
        tmp.unlink()
    hi = parse_clingen_hi(clingen_hi)
    mane = parse_mane_summary(mane_summary)
    conn = sqlite3.connect(tmp)
    try:
        conn.executescript(
            """
            PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;
            CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE lof_genes (gene TEXT PRIMARY KEY, decision TEXT NOT NULL, reason TEXT NOT NULL,
                                    hi_score TEXT, truncating_count INTEGER NOT NULL);
            CREATE TABLE truncating (gene TEXT NOT NULL, protein_pos INTEGER NOT NULL,
                                     variant_key TEXT, stars INTEGER NOT NULL);
            CREATE TABLE missense (gene TEXT NOT NULL, protein_pos INTEGER NOT NULL, ref_aa TEXT NOT NULL,
                                   alt_aa TEXT NOT NULL, variant_key TEXT, stars INTEGER NOT NULL,
                                   category TEXT NOT NULL);
            CREATE TABLE benign_missense (gene TEXT NOT NULL, protein_pos INTEGER NOT NULL,
                                          variant_key TEXT, stars INTEGER NOT NULL);
            CREATE TABLE missense_genes (gene TEXT PRIMARY KEY, plp INTEGER NOT NULL,
                                         blb INTEGER NOT NULL) WITHOUT ROWID;
            CREATE TABLE exons (transcript TEXT NOT NULL, exon_number INTEGER NOT NULL, cds_length INTEGER NOT NULL,
                                frameshift INTEGER NOT NULL, fraction REAL NOT NULL, last_exon INTEGER NOT NULL,
                                PRIMARY KEY (transcript, exon_number)) WITHOUT ROWID;
            """
        )
        truncating_keys: dict[str, set] = {}
        counts = {"clinvar_plp": 0, "truncating": 0, "missense": 0, "benign_missense": 0}
        trunc_rows, missense_rows, benign_rows = [], [], []
        gene_missense: dict[str, dict[str, set]] = {}
        for record in iter_clinvar_records(clinvar_summary):
            kind, alt, pos = classify_clinvar(record)
            gene = record["gene"]
            on_mane = record["nm"] in mane
            if record["category"] in ("LB", "B"):
                if on_mane and pos is not None and kind == "missense":
                    benign_rows.append((gene, pos, record["variant_key"], record["stars"]))
                    if record["stars"] >= GENE_COUNT_MIN_STARS:
                        gene_missense.setdefault(gene, {"plp": set(), "blb": set()})["blb"].add(
                            record["variant_key"] or record["hgvs_c"])
                continue
            counts["clinvar_plp"] += 1
            if kind in ("truncating", "splice") and record["stars"] >= 1:
                truncating_keys.setdefault(gene, set()).add(record["variant_key"] or record["hgvs_c"])
            if not on_mane or pos is None:
                continue
            if kind == "truncating" and record["stars"] >= 1:
                trunc_rows.append((gene, pos, record["variant_key"], record["stars"]))
            elif kind == "missense":
                ref_aa = _P_RE.match(record["hgvs_p"]).group("ref")
                missense_rows.append((gene, pos, ref_aa, alt, record["variant_key"], record["stars"],
                                      record["category"]))
                if record["stars"] >= GENE_COUNT_MIN_STARS:
                    gene_missense.setdefault(gene, {"plp": set(), "blb": set()})["plp"].add(
                        record["variant_key"] or record["hgvs_c"])
        conn.executemany("INSERT INTO truncating VALUES (?,?,?,?)", trunc_rows)
        conn.executemany("INSERT INTO missense VALUES (?,?,?,?,?,?,?)", missense_rows)
        conn.executemany("INSERT INTO benign_missense VALUES (?,?,?,?)", benign_rows)
        conn.executemany("INSERT INTO missense_genes VALUES (?,?,?)",
                         [(g, len(v["plp"]), len(v["blb"])) for g, v in sorted(gene_missense.items())])
        counts["truncating"], counts["missense"] = len(trunc_rows), len(missense_rows)
        counts["benign_missense"], counts["missense_genes"] = len(benign_rows), len(gene_missense)
        lof_rows = []
        for gene in sorted(set(hi) | set(truncating_keys)):
            score = hi.get(gene, "")
            n_trunc = len(truncating_keys.get(gene, ()))
            if score == "3":
                lof_rows.append((gene, "yes", "ClinGen HI=3", score, n_trunc))
            elif score == "40":
                lof_rows.append((gene, "no", "ClinGen HI=40（劑量不敏感）", score, n_trunc))
            elif n_trunc >= min_truncating:
                lof_rows.append((gene, "yes", f"ClinVar 有 {n_trunc} 個致病截斷變異（≥1★）", score, n_trunc))
        conn.executemany("INSERT INTO lof_genes VALUES (?,?,?,?,?)", lof_rows)
        counts["lof_genes_yes"] = sum(1 for r in lof_rows if r[1] == "yes")
        if mane_gtf:
            exon_rows = [(tx, n, v["cds_length"], int(v["frameshift"]), v["fraction"], int(v["last_exon"]))
                         for (tx, n), v in parse_mane_gtf(mane_gtf).items()]
            conn.executemany("INSERT INTO exons VALUES (?,?,?,?,?,?)", exon_rows)
            counts["exons"] = len(exon_rows)
        conn.executescript(
            """
            CREATE INDEX truncating_gene ON truncating(gene, protein_pos);
            CREATE INDEX missense_residue ON missense(gene, protein_pos);
            CREATE INDEX benign_missense_residue ON benign_missense(gene, protein_pos);
            """
        )
        meta = {
            "schema_version": str(SCHEMA_VERSION),
            "built_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "clinvar_release": clinvar_release,
            "sources": json.dumps({
                "clinvar_summary": file_signature(clinvar_summary),
                "mane_summary": file_signature(mane_summary),
                "clingen_hi": file_signature(clingen_hi),
                "mane_gtf": file_signature(mane_gtf),
            }, ensure_ascii=False),
            "counts": json.dumps(counts),
            "min_truncating": str(min_truncating),
            "gene_count_min_stars": str(GENE_COUNT_MIN_STARS),
        }
        conn.executemany("INSERT INTO meta VALUES (?,?)", sorted(meta.items()))
        conn.commit()
        check = conn.execute("PRAGMA quick_check").fetchone()[0]
        if check != "ok":
            raise RuntimeError(f"acmg_refs quick_check failed: {check}")
    finally:
        conn.close()
    os.replace(tmp, final)
    return counts


# ── read ───────────────────────────────────────────────────────────────────
def _missense_entry(alt_aa, variant_key, stars, category) -> dict:
    entry = {"alt_aa": alt_aa, "variant_key": variant_key, "stars": stars}
    if category:
        entry["category"] = category
    return entry


@lru_cache(maxsize=4)
def _load_gene_table(path: Path) -> dict:
    try:
        return json.loads(Path(path).read_text(encoding="utf-8")).get("genes") or {}
    except (OSError, ValueError):
        return {}


def load_vcep_specs(path: Path = VCEP_SPECS_PATH) -> dict:
    """gene -> VCEP spec (built by scripts/acmg_engine/update_vcep_specs.py)."""
    return _load_gene_table(Path(path))


def load_collagen(path: Path = COLLAGEN_PATH) -> dict:
    """gene -> {"accession", "length", "positions"} with positions as a frozenset."""
    out = {}
    for gene, entry in _load_gene_table(Path(path)).items():
        out[gene] = dict(entry, positions=frozenset(entry.get("positions") or ()))
    return out


def load_ba1_exceptions(path: Path = BA1_EXCEPTIONS_PATH) -> set[str]:
    out: set[str] = set()
    try:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                key = line.split("\t")[0].strip()
                if key and key != "variant_key":
                    out.add(key)
    except OSError:
        pass
    return out


class SqliteReferences(References):
    """References backed by ``acmg_refs.sqlite`` (+ optional WES recurrence store)."""

    def __init__(self, refs_db: Path | None, recurrence=None, ba1_exceptions: set | None = None,
                 cases_db: Path | None = None):
        self.refs_db = Path(refs_db) if refs_db else None
        self._conn = None
        if self.refs_db and self.refs_db.is_file():
            self._conn = sqlite3.connect(f"file:{self.refs_db}?mode=ro", uri=True, check_same_thread=False)
        # PS4 B: ClinVar affected cases (clinvar_cases.sqlite next to the refs DB, optional)
        if cases_db is None and self.refs_db is not None:
            from .clinvar_cases import CASES_DB_NAME
            cases_db = self.refs_db.with_name(CASES_DB_NAME)
        self.cases_db = Path(cases_db) if cases_db else None
        self._cases = None
        self._cases_pm3 = False
        if self.cases_db and self.cases_db.is_file():
            self._cases = sqlite3.connect(f"file:{self.cases_db}?mode=ro", uri=True, check_same_thread=False)
            self._cases_pm3 = "pm3_submitters" in {r[1] for r in self._cases.execute("PRAGMA table_info(cases)")}
        self.recurrence = recurrence
        self._ba1 = load_ba1_exceptions() if ba1_exceptions is None else set(ba1_exceptions)
        self._vcep = load_vcep_specs()
        self._collagen = load_collagen()
        self._tables = set()
        self._missense_category = False
        if self._conn is not None:
            self._tables = {"lof_genes", "truncating_positions", "missense_by_residue"}
            if self._conn.execute("SELECT 1 FROM exons LIMIT 1").fetchone():
                self._tables.add("exons")
            names = {r[0] for r in self._conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            if {"benign_missense", "missense_genes"} <= names:      # schema ≥2
                self._tables.update({"missense_window", "missense_genes"})
            self._missense_category = "category" in {
                r[1] for r in self._conn.execute("PRAGMA table_info(missense)")}
        if recurrence is not None:
            self._tables.add("recurrence")
        if self._cases is not None:
            self._tables.add("clinvar_cases")
        self.lof_mechanism = lru_cache(maxsize=None)(self._lof_mechanism)

    @property
    def available(self) -> frozenset:  # type: ignore[override]
        return frozenset(self._tables)

    @property
    def version(self) -> dict:  # type: ignore[override]
        info = {"refs_db": file_signature(self.refs_db)}
        if self._conn is not None:
            info.update(dict(self._conn.execute("SELECT key, value FROM meta WHERE key IN "
                                                "('schema_version','built_at','clinvar_release')")))
        if self._cases is not None:
            info["clinvar_cases"] = file_signature(self.cases_db)
            release = self._cases.execute("SELECT value FROM meta WHERE key='release'").fetchone()
            if release:
                info["clinvar_cases_release"] = release[0]
        if self.recurrence is not None:
            info["wes_recurrence"] = self.recurrence.signature()
        return info

    def clinvar_affected(self, variation_id):
        if self._cases is None:
            return None
        if not variation_id:
            return {}
        hit = self._cases.execute("SELECT submitters, individuals FROM cases WHERE variation_id=?",
                                  (str(variation_id),)).fetchone()
        return {"submitters": int(hit[0]), "individuals": int(hit[1])} if hit else {}

    def clinvar_biallelic(self, variation_id):
        if self._cases is None or not self._cases_pm3:
            return None
        if not variation_id:
            return {}
        hit = self._cases.execute("SELECT pm3_submitters, pm3_hom, pm3_comphet, pm3_detail FROM cases "
                                  "WHERE variation_id=?", (str(variation_id),)).fetchone()
        if not hit or not hit[0]:
            return {}
        try:
            detail = json.loads(hit[3] or "[]")
        except ValueError:
            detail = []
        return {"submitters": int(hit[0]), "hom": int(hit[1]), "comphet": int(hit[2]), "detail": detail}

    def _lof_mechanism(self, gene):
        if self._conn is None or not gene:
            return None
        hit = self._conn.execute("SELECT decision, reason FROM lof_genes WHERE gene=?", (gene,)).fetchone()
        return (hit[0], hit[1]) if hit else None

    def truncating_downstream(self, gene, protein_pos):
        if self._conn is None:
            return None
        hit = self._conn.execute("SELECT 1 FROM truncating WHERE gene=? AND protein_pos>? LIMIT 1",
                                 (gene, int(protein_pos))).fetchone()
        return bool(hit)

    def clinvar_missense(self, gene, protein_pos):
        if self._conn is None:
            return None
        category = "category" if self._missense_category else "NULL"
        rows = self._conn.execute(
            f"SELECT alt_aa, variant_key, stars, {category} FROM missense WHERE gene=? AND protein_pos=?",
            (gene, int(protein_pos))).fetchall()
        return [_missense_entry(a, k, s, c) for a, k, s, c in rows]

    def clinvar_missense_window(self, gene, start, end):
        if self._conn is None or "missense_window" not in self._tables:
            return None
        rows = self._conn.execute(
            "SELECT protein_pos, alt_aa, variant_key, stars, category FROM missense "
            "WHERE gene=? AND protein_pos BETWEEN ? AND ?", (gene, int(start), int(end))).fetchall()
        return [dict(_missense_entry(a, k, s, c), protein_pos=p) for p, a, k, s, c in rows]

    def benign_missense_window(self, gene, start, end):
        if self._conn is None or "missense_window" not in self._tables:
            return None
        rows = self._conn.execute(
            "SELECT protein_pos, variant_key, stars FROM benign_missense "
            "WHERE gene=? AND protein_pos BETWEEN ? AND ?", (gene, int(start), int(end))).fetchall()
        return [{"protein_pos": p, "variant_key": k, "stars": s} for p, k, s in rows]

    def missense_gene_counts(self, gene):
        if self._conn is None or "missense_genes" not in self._tables or not gene:
            return None
        hit = self._conn.execute("SELECT plp, blb FROM missense_genes WHERE gene=?", (gene,)).fetchone()
        return (int(hit[0]), int(hit[1])) if hit else (0, 0)

    def exon_skip(self, transcript, exon_number):
        if self._conn is None or "exons" not in self._tables or not exon_number:
            return None
        hit = self._conn.execute(
            "SELECT frameshift, fraction, last_exon FROM exons WHERE transcript=? AND exon_number=?",
            (str(transcript).split(".")[0], int(exon_number))).fetchone()
        if not hit:
            return None
        return {"frameshift": bool(hit[0]), "fraction": float(hit[1]), "last_exon": bool(hit[2])}

    def cds_exons(self, transcript):
        if self._conn is None or "exons" not in self._tables or not transcript:
            return None
        rows = self._conn.execute("SELECT cds_length FROM exons WHERE transcript=? ORDER BY exon_number",
                                  (str(transcript).split(".")[0],)).fetchall()
        return [int(r[0]) for r in rows] or None

    def ba1_exception(self, variant_key):
        return variant_key in self._ba1

    def vcep(self, gene):
        return self._vcep.get(gene)

    def collagen_gly(self, gene):
        return self._collagen.get(gene)

    def wes_recurrence(self, variant_key):
        if self.recurrence is None:
            return None
        return self.recurrence.lookup(variant_key)


def default_refs_dir() -> Path:
    try:
        from ... import config  # noqa: WPS433 - lazy to keep the engine importable standalone
        return Path(config.ACMG_REFS_DIR)
    except Exception:  # pragma: no cover - fallback outside the app
        home = Path(os.environ.get("NGS_UI_HOME") or Path.home() / "NGS_UI")
        return Path(os.environ.get("NGS_UI_ACMG_REFS_DIR") or home / "biotools" / "acmg_refs")


def load_references(refs_dir: Path | None = None, recurrence=None) -> SqliteReferences:
    refs_dir = Path(refs_dir) if refs_dir else default_refs_dir()
    return SqliteReferences(refs_dir / REFS_DB_NAME, recurrence=recurrence)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()
