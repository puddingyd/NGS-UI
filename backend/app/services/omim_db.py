"""Platform-owned OMIM database (SQLite) replacing OMIM.xlsx as the source.

Built once from the curator workbook (``scripts/omim_db.py init``) and then
updated by adopting OMIM page captures from the browser extension.
``omim_store`` reads this database whenever the file exists, otherwise it
keeps reading OMIM.xlsx, so nothing changes until the operator runs init.

Model
-----
* ``genes``          one row per workbook gene row (gene MIM, symbol,
                      gene-level Inheritance, Done date).
* ``gene_diseases``  the skeleton: one row per disease of a gene, seeded
                      from the workbook ``OMIM_disease`` lines, in order.
                      ``label`` is the workbook line ("Name (MIM)(INH)").
* ``legacy_cells``   the original Disease1..16 cell texts, verbatim.
                      ``gene_diseases.legacy_slot`` points at the cell the
                      migration matched to that disease.
* ``phenotype_texts`` adopted capture content per phenotype MIM
                      (Description + Clinical Synopsis + OMIM dates); shared
                      by every gene linked to that MIM.
* ``review_items``   things a curator must confirm (migration mismatches,
                      diseases OMIM no longer lists, …). Never auto-applied.
* ``change_log``     every adoption / removal, with before/after values.

Rendering (``load_rows``) produces the same row dicts ``omim_store`` always
served (OMIM_id, gene_symbol, OMIM_disease, Inheritance, Disease1..N, Done)
so sample loading, DOCX and the case list need no changes. Each Disease
cell is, in priority order: adopted phenotype text → the legacy workbook
cell → the label alone (shown with the "*" missing-detail marker).
``Disease_legacy_slots`` maps old workbook slot numbers to rendered slots so
legacy slot-number report ticks keep pointing at the same disease.
"""
from __future__ import annotations

import difflib
import io
import json
import re
import shutil
import sqlite3
import threading
from contextlib import contextmanager
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

from ..config import OMIM_DB
from . import disease_picks

SCHEMA_VERSION = 1
MAX_RENDERED_SLOTS = 16
EXCEL_CELL_LIMIT = 32767
_PLACEHOLDER_LABELS = {"nil", "na", "none", "-"}
_INHERITANCE_SPLIT_RE = re.compile(r"\s*[,;/]\s*")

_write_lock = threading.RLock()

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS genes (
    gene_mim     TEXT PRIMARY KEY,
    gene_symbol  TEXT NOT NULL DEFAULT '',
    alias_symbols TEXT NOT NULL DEFAULT '',
    inheritance  TEXT NOT NULL DEFAULT '',
    done         TEXT NOT NULL DEFAULT '',
    legacy_done  TEXT NOT NULL DEFAULT '',
    row_order    INTEGER NOT NULL DEFAULT 0,
    updated_at   TEXT
);
CREATE INDEX IF NOT EXISTS genes_symbol ON genes(gene_symbol);
CREATE TABLE IF NOT EXISTS gene_diseases (
    id             INTEGER PRIMARY KEY,
    gene_mim       TEXT NOT NULL,
    ordinal        INTEGER NOT NULL,
    label          TEXT NOT NULL,
    phenotype_mim  TEXT NOT NULL DEFAULT '',
    inheritance    TEXT NOT NULL DEFAULT '',
    legacy_slot    INTEGER,
    source         TEXT NOT NULL DEFAULT 'xlsx',
    active         INTEGER NOT NULL DEFAULT 1,
    updated_at     TEXT
);
CREATE INDEX IF NOT EXISTS gene_diseases_gene ON gene_diseases(gene_mim, ordinal);
CREATE INDEX IF NOT EXISTS gene_diseases_mim ON gene_diseases(phenotype_mim);
CREATE TABLE IF NOT EXISTS legacy_cells (
    gene_mim  TEXT NOT NULL,
    slot      INTEGER NOT NULL,
    text      TEXT NOT NULL,
    PRIMARY KEY (gene_mim, slot)
);
CREATE TABLE IF NOT EXISTS phenotype_texts (
    phenotype_mim          TEXT PRIMARY KEY,
    title                  TEXT NOT NULL DEFAULT '',
    description            TEXT NOT NULL DEFAULT '',
    description_source     TEXT NOT NULL DEFAULT '',
    synopsis               TEXT NOT NULL DEFAULT '',
    synopsis_status        TEXT NOT NULL DEFAULT '',
    omim_edit_date         TEXT NOT NULL DEFAULT '',
    capture_id             TEXT NOT NULL DEFAULT '',
    content_hash           TEXT NOT NULL DEFAULT '',
    adopted_at             TEXT NOT NULL,
    adopted_by             TEXT NOT NULL DEFAULT ''
);
CREATE TABLE IF NOT EXISTS review_items (
    id           INTEGER PRIMARY KEY,
    kind         TEXT NOT NULL,
    gene_mim     TEXT NOT NULL DEFAULT '',
    disease_id   INTEGER,
    detail       TEXT NOT NULL DEFAULT '{}',
    status       TEXT NOT NULL DEFAULT 'open',
    created_at   TEXT NOT NULL,
    resolved_at  TEXT,
    resolved_by  TEXT,
    resolution   TEXT
);
CREATE INDEX IF NOT EXISTS review_open ON review_items(status, kind);
CREATE TABLE IF NOT EXISTS adoptions (
    id          INTEGER PRIMARY KEY,
    capture_id  TEXT NOT NULL,
    mim         TEXT NOT NULL,
    entry_type  TEXT NOT NULL DEFAULT '',
    title       TEXT NOT NULL DEFAULT '',
    adopted_at  TEXT NOT NULL,
    adopted_by  TEXT NOT NULL DEFAULT '',
    summary     TEXT NOT NULL DEFAULT '{}'
);
CREATE INDEX IF NOT EXISTS adoptions_capture ON adoptions(capture_id);
CREATE TABLE IF NOT EXISTS change_log (
    id          INTEGER PRIMARY KEY,
    at          TEXT NOT NULL,
    user        TEXT NOT NULL DEFAULT '',
    action      TEXT NOT NULL,
    gene_mim    TEXT NOT NULL DEFAULT '',
    phenotype_mim TEXT NOT NULL DEFAULT '',
    capture_id  TEXT NOT NULL DEFAULT '',
    before      TEXT,
    after       TEXT
);
"""

REVIEW_KINDS = {
    "renamed_match": "OMIM_disease 與 Disease 欄名稱不同（已配對，請確認）",
    "slot_not_in_omim_disease": "Disease 欄有、OMIM_disease 沒有的疾病",
    "missing_mim": "OMIM_disease 這行沒有 phenotype MIM",
    "overflow": "疾病超過 16 個，卡片只顯示前 16 個",
    "not_in_omim": "OMIM 基因頁已不列出這個疾病",
    "unmatched_label": "基因頁擷取時無法自動對應的舊疾病名稱",
}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _today() -> str:
    return date.today().isoformat()


def db_path() -> Path:
    return OMIM_DB


def exists(path: Path | None = None) -> bool:
    return (path or db_path()).is_file()


_schema_checked: set[str] = set()


def _ensure_schema(path: Path) -> None:
    """Add tables introduced after a database was built (idempotent)."""
    key = str(path)
    if key in _schema_checked or not path.is_file():
        return
    con = sqlite3.connect(key, timeout=30)
    try:
        con.executescript(SCHEMA)
        gene_cols = {r[1] for r in con.execute("PRAGMA table_info(genes)")}
        if "legacy_done" not in gene_cols:
            # Workbook Done date, kept apart from Done (which adoption updates)
            # so stale curated text is judged by when it was really curated.
            con.execute("ALTER TABLE genes ADD COLUMN legacy_done TEXT NOT NULL DEFAULT ''")
            con.execute("UPDATE genes SET legacy_done = done")
        if not con.execute("SELECT 1 FROM adoptions LIMIT 1").fetchone():
            # Back-fill adoptions recorded only in change_log (pre-v10.08).
            con.execute(
                "INSERT INTO adoptions(capture_id, mim, entry_type, adopted_at, adopted_by, summary) "
                "SELECT capture_id, CASE WHEN phenotype_mim != '' THEN phenotype_mim ELSE gene_mim END, "
                "'', MIN(at), user, '{}' FROM change_log WHERE capture_id != '' GROUP BY capture_id"
            )
        con.commit()
    finally:
        con.close()
    _schema_checked.add(key)


@contextmanager
def connect(path: Path | None = None, *, write: bool = False):
    path = path or db_path()
    _ensure_schema(Path(path))
    con = sqlite3.connect(str(path), timeout=30)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA foreign_keys = ON")
    try:
        if write:
            con.execute("BEGIN IMMEDIATE")
        yield con
        if write:
            con.commit()
    except Exception:
        if write:
            con.rollback()
        raise
    finally:
        con.close()


def _norm_gene_mim(value) -> str:
    text = str(value or "").strip()
    if text.endswith(".0"):
        text = text[:-2]
    return text


def _label_parts(label: str) -> tuple[str, str, str]:
    """(name, phenotype MIM, inheritance) from a "Name (MIM)(INH)" line."""
    line = disease_picks.first_line(label)
    mim = disease_picks.phenotype_mim(line)
    name = line
    inheritance = ""
    groups = []
    previous = None
    while previous != name:
        previous = name
        match = re.search(r"\s*\(([^()]*)\)\s*$", name)
        if match:
            groups.insert(0, match.group(1).strip())
            name = name[:match.start()]
    for group in groups:
        if re.fullmatch(r"\d{6}", group):
            continue
        tokens = [t for t in _INHERITANCE_SPLIT_RE.split(group) if t]
        if tokens and all(len(t) <= 14 for t in tokens):
            inheritance = group
    return name.strip(), mim, inheritance


def format_label(name: str, phenotype_mim: str, inheritance: str) -> str:
    name = str(name or "").strip()
    out = name
    if phenotype_mim:
        out += f" ({phenotype_mim})"
    out += f"({inheritance})" if phenotype_mim else (f" ({inheritance})" if inheritance else "")
    return out


# --------------------------------------------------------------------------
# Migration from OMIM.xlsx
# --------------------------------------------------------------------------

def _read_workbook(xlsx_path: Path) -> list[dict]:
    import openpyxl

    wb = openpyxl.load_workbook(str(xlsx_path), read_only=True, data_only=True)
    ws = wb.worksheets[0]
    rows = ws.iter_rows(values_only=True)
    headers = [str(h or "").strip() for h in next(rows, ())]
    h2i = {h: i for i, h in enumerate(headers) if h}
    slot_fields = [f"Disease{i}" for i in range(1, 65) if f"Disease{i}" in h2i]

    def cell(raw, name):
        i = h2i.get(name)
        if i is None or i >= len(raw) or raw[i] is None:
            return None
        return raw[i]

    out = []
    for order, raw in enumerate(rows, start=1):
        if not raw:
            continue
        gene_mim = _norm_gene_mim(cell(raw, "OMIM_id"))
        symbol = str(cell(raw, "gene_symbol") or "").strip()
        if not gene_mim and not symbol:
            continue
        done = cell(raw, "Done")
        if hasattr(done, "date"):
            done = done.date().isoformat()
        out.append({
            "order": order,
            "gene_mim": gene_mim,
            "gene_symbol": symbol,
            "omim_disease": str(cell(raw, "OMIM_disease") or ""),
            "inheritance": str(cell(raw, "Inheritance") or "").strip(),
            "done": str(done or "").strip(),
            "slots": [
                (int(f[7:]), str(cell(raw, f)).strip())
                for f in slot_fields
                if cell(raw, f) is not None and str(cell(raw, f)).strip()
                and str(cell(raw, f)).strip() != "NA"
            ],
        })
    wb.close()
    return out


def match_rows(lines: list[str], slots: list[tuple[int, str]]) -> dict:
    """Pair OMIM_disease lines with legacy Disease cells.

    Returns {"pairs": [(line_idx, slot, kind)], "unmatched_lines": [...],
    "unmatched_slots": [...]}. ``kind`` is exact / mim / name / similar.
    Similarity pairing is only used when one side lacks a phenotype MIM,
    so two different MIMs are never merged (e.g. SCA1 vs SCA10).
    """
    line_keys = [(i, disease_picks.pick_key(l)) for i, l in enumerate(lines)]
    slot_keys = [(s, disease_picks.pick_key(t)) for s, t in slots]
    used_lines: set[int] = set()
    used_slots: set[int] = set()
    pairs: list[tuple[int, int, str]] = []

    def run(kind, predicate):
        for s, skey in slot_keys:
            if s in used_slots:
                continue
            hits = [i for i, lkey in line_keys if i not in used_lines and predicate(lkey, skey)]
            if len(hits) == 1:
                used_lines.add(hits[0])
                used_slots.add(s)
                pairs.append((hits[0], s, kind))

    run("exact", lambda a, b: a == b)
    run("mim", lambda a, b: a.split("|")[0] and a.split("|")[0] == b.split("|")[0])
    run("name", lambda a, b: a.split("|")[1] and a.split("|")[1] == b.split("|")[1])

    candidates = []
    for i, lkey in line_keys:
        if i in used_lines:
            continue
        for s, skey in slot_keys:
            if s in used_slots:
                continue
            lmim, lname = lkey.split("|", 1)
            smim, sname = skey.split("|", 1)
            if lmim and smim and lmim != smim:
                continue
            ratio = difflib.SequenceMatcher(None, lname, sname).ratio()
            if ratio >= 0.6:
                candidates.append((ratio, i, s))
    for _ratio, i, s in sorted(candidates, reverse=True):
        if i in used_lines or s in used_slots:
            continue
        used_lines.add(i)
        used_slots.add(s)
        pairs.append((i, s, "similar"))

    return {
        "pairs": sorted(pairs),
        "unmatched_lines": [i for i, _ in line_keys if i not in used_lines],
        "unmatched_slots": [s for s, _ in slot_keys if s not in used_slots],
    }


def build_from_xlsx(xlsx_path: Path, out_path: Path) -> dict:
    """Create a fresh database at ``out_path`` from the curator workbook."""
    rows = _read_workbook(xlsx_path)
    tmp = out_path.with_suffix(out_path.suffix + ".building")
    if tmp.exists():
        tmp.unlink()
    tmp.parent.mkdir(parents=True, exist_ok=True)
    stats = {
        "genes": 0, "diseases": 0, "legacy_cells": 0, "duplicate_gene_mim": 0,
        "pairs_exact": 0, "pairs_mim": 0, "pairs_name": 0, "pairs_similar": 0,
        "label_only": 0, "newly_visible": 0, "slot_only": 0, "missing_mim": 0,
        "overflow_genes": 0, "placeholder_lines": 0,
    }
    now = _now()
    con = sqlite3.connect(str(tmp))
    try:
        con.executescript(SCHEMA)
        con.execute("INSERT INTO meta(key, value) VALUES ('schema_version', ?)", (str(SCHEMA_VERSION),))
        con.execute("INSERT INTO meta(key, value) VALUES ('source_xlsx', ?)", (str(xlsx_path),))
        con.execute("INSERT INTO meta(key, value) VALUES ('created_at', ?)", (now,))
        merged: dict[str, dict] = {}
        for row in rows:
            key = row["gene_mim"] or f"symbol:{row['gene_symbol']}"
            first = merged.get(key)
            if first is None:
                merged[key] = dict(row, aliases=[])
                continue
            # Same gene MIM under an old and a new symbol (e.g. BLM / RECQL3):
            # keep one gene, remember the other symbol as an alias, and keep
            # whichever row carries more curated Disease cells.
            stats["duplicate_gene_mim"] += 1
            if row["gene_symbol"] and row["gene_symbol"] not in first["aliases"]:
                first["aliases"].append(row["gene_symbol"])
            if len(row["slots"]) > len(first["slots"]):
                first["slots"] = row["slots"]
                if row["omim_disease"].strip():
                    first["omim_disease"] = row["omim_disease"]
                first["done"] = row["done"] or first["done"]
            elif not first["omim_disease"].strip():
                first["omim_disease"] = row["omim_disease"]
        rows = list(merged.values())
        seen: set[str] = set()

        def review(kind, gene_mim, disease_id, detail):
            con.execute(
                "INSERT INTO review_items(kind, gene_mim, disease_id, detail, created_at) VALUES (?,?,?,?,?)",
                (kind, gene_mim, disease_id, json.dumps(detail, ensure_ascii=False), now),
            )

        for row in rows:
            gene_mim = row["gene_mim"] or f"symbol:{row['gene_symbol']}"
            seen.add(gene_mim)
            stats["genes"] += 1
            con.execute(
                "INSERT INTO genes(gene_mim, gene_symbol, alias_symbols, inheritance, done, legacy_done, "
                "row_order, updated_at) VALUES (?,?,?,?,?,?,?,?)",
                (gene_mim, row["gene_symbol"], ",".join(row["aliases"]), row["inheritance"],
                 row["done"], row["done"], row["order"], now),
            )
            for slot, text in row["slots"]:
                con.execute(
                    "INSERT INTO legacy_cells(gene_mim, slot, text) VALUES (?,?,?)",
                    (gene_mim, slot, text),
                )
                stats["legacy_cells"] += 1

            lines = []
            for line in str(row["omim_disease"]).splitlines():
                line = line.strip()
                if not line:
                    continue
                if line.lower() in _PLACEHOLDER_LABELS:
                    stats["placeholder_lines"] += 1
                    continue
                lines.append(line)
            matched = match_rows(lines, row["slots"])
            slot_for_line = {i: (s, kind) for i, s, kind in matched["pairs"]}
            slot_text = dict(row["slots"])
            ordinal = 0
            for i, line in enumerate(lines):
                ordinal += 1
                _name, mim, inheritance = _label_parts(line)
                slot, kind = slot_for_line.get(i, (None, ""))
                cur = con.execute(
                    "INSERT INTO gene_diseases(gene_mim, ordinal, label, phenotype_mim, inheritance, legacy_slot, source, updated_at) "
                    "VALUES (?,?,?,?,?,?, 'xlsx', ?)",
                    (gene_mim, ordinal, line, mim, inheritance, slot, now),
                )
                disease_id = cur.lastrowid
                stats["diseases"] += 1
                if kind:
                    stats[f"pairs_{kind}"] += 1
                else:
                    stats["label_only"] += 1
                    if row["slots"]:
                        stats["newly_visible"] += 1
                if kind == "similar":
                    review("renamed_match", gene_mim, disease_id, {
                        "gene": row["gene_symbol"], "omim_disease": line,
                        "disease_cell": disease_picks.first_line(slot_text[slot]), "slot": slot,
                    })
                if not mim:
                    stats["missing_mim"] += 1
                    review("missing_mim", gene_mim, disease_id, {
                        "gene": row["gene_symbol"], "omim_disease": line,
                    })
            for slot in matched["unmatched_slots"]:
                ordinal += 1
                label = disease_picks.first_line(slot_text[slot])
                _name, mim, inheritance = _label_parts(label)
                cur = con.execute(
                    "INSERT INTO gene_diseases(gene_mim, ordinal, label, phenotype_mim, inheritance, legacy_slot, source, updated_at) "
                    "VALUES (?,?,?,?,?,?, 'xlsx_slot', ?)",
                    (gene_mim, ordinal, label, mim, inheritance, slot, now),
                )
                stats["diseases"] += 1
                stats["slot_only"] += 1
                review("slot_not_in_omim_disease", gene_mim, cur.lastrowid, {
                    "gene": row["gene_symbol"], "disease_cell": label, "slot": slot,
                    "omim_disease": lines,
                })
            if ordinal > MAX_RENDERED_SLOTS:
                stats["overflow_genes"] += 1
                review("overflow", gene_mim, None, {
                    "gene": row["gene_symbol"], "count": ordinal,
                })
        con.commit()
    finally:
        con.close()
    stats["review_items"] = _count_open_reviews(tmp)
    return {"path": str(tmp), "stats": stats}


def _count_open_reviews(path: Path) -> dict:
    with connect(path) as con:
        rows = con.execute(
            "SELECT kind, COUNT(*) AS n FROM review_items WHERE status='open' GROUP BY kind"
        ).fetchall()
    return {r["kind"]: r["n"] for r in rows}


def install(built_path: Path, target: Path | None = None) -> Path | None:
    """Atomically move a built database into place, backing up any old one."""
    target = target or db_path()
    backup = None
    if target.exists():
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        backup = target.with_name(f"{target.stem}.{stamp}.bak{target.suffix}")
        shutil.copy2(target, backup)
    Path(built_path).replace(target)
    return backup


# --------------------------------------------------------------------------
# Rendering for omim_store
# --------------------------------------------------------------------------

def compose_cell(label: str, text: dict | None) -> str:
    parts = [label]
    if text:
        if text.get("description"):
            parts.append(text["description"])
        if text.get("synopsis"):
            parts.append(text["synopsis"])
    return "\n\n".join(p for p in parts if p)


def _render_gene(gene: dict, diseases: list[dict], legacy: dict[int, str],
                 texts: dict[str, dict], slot_count: int) -> dict:
    row = {
        "OMIM_id": gene["gene_mim"] if not gene["gene_mim"].startswith("symbol:") else "",
        "gene_symbol": gene["gene_symbol"],
        "Inheritance": gene["inheritance"],
        "Done": gene["done"],
        "Legacy_done": gene.get("legacy_done", gene["done"]),
        "OMIM_disease": "\n".join(d["label"] for d in diseases),
    }
    meta = []
    legacy_map = {}
    for idx in range(1, MAX_RENDERED_SLOTS + 1):
        row[f"Disease{idx}"] = ""
    for idx, d in enumerate(diseases, start=1):
        text = texts.get(d["phenotype_mim"]) if d["phenotype_mim"] else None
        if text:
            cell = compose_cell(d["label"], text)
            source = "capture"
        elif d["legacy_slot"] and legacy.get(d["legacy_slot"]):
            cell = legacy[d["legacy_slot"]]
            source = "xlsx"
        else:
            cell = d["label"]
            source = "label"
        if d["legacy_slot"]:
            legacy_map[str(d["legacy_slot"])] = idx
        meta.append({
            "slot": idx,
            "disease_id": d["id"],
            "phenotype_mim": d["phenotype_mim"],
            "source": source,
            "omim_edit_date": (text or {}).get("omim_edit_date", ""),
            "adopted_at": (text or {}).get("adopted_at", ""),
        })
        if idx <= slot_count:
            row[f"Disease{idx}"] = cell
    row["Disease_legacy_slots"] = legacy_map
    row["Disease_meta"] = meta
    return row


def load_rows(path: Path | None = None, slot_count: int = MAX_RENDERED_SLOTS) -> tuple[dict, dict]:
    """Render every gene into omim_store row dicts: (by_omim_id, by_gene)."""
    with connect(path) as con:
        genes = [dict(r) for r in con.execute("SELECT * FROM genes ORDER BY row_order")]
        diseases: dict[str, list[dict]] = {}
        for r in con.execute(
            "SELECT * FROM gene_diseases WHERE active = 1 ORDER BY gene_mim, ordinal, id"
        ):
            diseases.setdefault(r["gene_mim"], []).append(dict(r))
        legacy: dict[str, dict[int, str]] = {}
        for r in con.execute("SELECT gene_mim, slot, text FROM legacy_cells"):
            legacy.setdefault(r["gene_mim"], {})[r["slot"]] = r["text"]
        texts = {r["phenotype_mim"]: dict(r) for r in con.execute("SELECT * FROM phenotype_texts")}
    by_id: dict[int, dict] = {}
    by_gene: dict[str, dict] = {}
    for gene in genes:
        row = _render_gene(
            gene, diseases.get(gene["gene_mim"], []), legacy.get(gene["gene_mim"], {}),
            texts, slot_count,
        )
        try:
            by_id.setdefault(int(row["OMIM_id"]), row)
        except (TypeError, ValueError):
            pass
        if row["gene_symbol"]:
            by_gene.setdefault(row["gene_symbol"], row)
        for alias in filter(None, (gene.get("alias_symbols") or "").split(",")):
            by_gene.setdefault(alias, row)
    return by_id, by_gene


def signature(path: Path | None = None) -> tuple:
    path = path or db_path()
    try:
        st = path.stat()
    except OSError:
        return ("omim-db", None)
    wal = path.with_name(path.name + "-wal")
    wal_m = wal.stat().st_mtime_ns if wal.exists() else 0
    return ("omim-db", st.st_mtime_ns, st.st_size, wal_m)


# --------------------------------------------------------------------------
# Adopting captures
# --------------------------------------------------------------------------

def _log(con, user, action, gene_mim="", phenotype_mim="", capture_id="", before=None, after=None):
    con.execute(
        "INSERT INTO change_log(at, user, action, gene_mim, phenotype_mim, capture_id, before, after) "
        "VALUES (?,?,?,?,?,?,?,?)",
        (_now(), user, action, gene_mim, phenotype_mim, capture_id,
         json.dumps(before, ensure_ascii=False) if before is not None else None,
         json.dumps(after, ensure_ascii=False) if after is not None else None),
    )


def _gene_row(con, gene_mim: str, symbol: str, user: str, capture_id: str) -> dict:
    row = con.execute("SELECT * FROM genes WHERE gene_mim = ?", (gene_mim,)).fetchone()
    if row:
        return dict(row)
    if symbol:
        row = con.execute("SELECT * FROM genes WHERE gene_symbol = ?", (symbol,)).fetchone()
        if row:
            return dict(row)
    order = con.execute("SELECT COALESCE(MAX(row_order), 0) + 1 FROM genes").fetchone()[0]
    con.execute(
        "INSERT INTO genes(gene_mim, gene_symbol, inheritance, done, row_order, updated_at) VALUES (?,?,?,?,?,?)",
        (gene_mim, symbol, "", "", order, _now()),
    )
    _log(con, user, "add_gene", gene_mim, "", capture_id, None, {"gene_symbol": symbol})
    return dict(con.execute("SELECT * FROM genes WHERE gene_mim = ?", (gene_mim,)).fetchone())


def _merge_inheritance(values) -> str:
    seen = []
    for value in values:
        for token in _INHERITANCE_SPLIT_RE.split(str(value or "")):
            token = token.strip()
            if token and token not in seen:
                seen.append(token)
    return ", ".join(seen)


def _refresh_gene(con, gene_mim: str, *, done: bool) -> None:
    labels = con.execute(
        "SELECT inheritance FROM gene_diseases WHERE gene_mim = ? AND active = 1 ORDER BY ordinal",
        (gene_mim,),
    ).fetchall()
    inheritance = _merge_inheritance(r["inheritance"] for r in labels)
    if done:
        con.execute(
            "UPDATE genes SET inheritance = ?, done = ?, updated_at = ? WHERE gene_mim = ?",
            (inheritance, _today(), _now(), gene_mim),
        )
    else:
        con.execute(
            "UPDATE genes SET inheritance = ?, updated_at = ? WHERE gene_mim = ?",
            (inheritance, _now(), gene_mim),
        )


def _upsert_relationship(con, gene_mim: str, rel: dict, *, user: str, capture_id: str) -> dict:
    """Make sure the gene has a disease row for this relationship."""
    mim = rel.get("phenotype_mim") or ""
    name = rel.get("phenotype") or ""
    inheritance = rel.get("inheritance") or ""
    label = format_label(name, mim, inheritance)
    rows = [dict(r) for r in con.execute(
        "SELECT * FROM gene_diseases WHERE gene_mim = ? AND active = 1 ORDER BY ordinal", (gene_mim,)
    )]
    target = None
    if mim:
        same = [r for r in rows if r["phenotype_mim"] == mim]
        exact = [r for r in same if disease_picks.name_key(r["label"]) == disease_picks.name_key(label)]
        if len(exact) == 1:
            target = exact[0]
        elif len(same) == 1:
            target = same[0]
    if target is None:
        named = [r for r in rows if not r["phenotype_mim"]
                 and disease_picks.name_key(r["label"]) == disease_picks.name_key(label)]
        target = named[0] if len(named) == 1 else None
    if target is None:
        # A stale workbook line (old name, no MIM) whose curated Disease cell
        # already carries the current name / MIM, e.g. CERT1.
        via_cell = []
        for r in rows:
            if r["phenotype_mim"] or not r["legacy_slot"]:
                continue
            cell = con.execute(
                "SELECT text FROM legacy_cells WHERE gene_mim = ? AND slot = ?",
                (gene_mim, r["legacy_slot"]),
            ).fetchone()
            if not cell:
                continue
            first = disease_picks.first_line(cell["text"])
            if (mim and disease_picks.phenotype_mim(first) == mim) or (
                disease_picks.name_key(first) == disease_picks.name_key(label)
            ):
                via_cell.append(r)
        target = via_cell[0] if len(via_cell) == 1 else None
    if target is None:
        ordinal = (max((r["ordinal"] for r in rows), default=0)) + 1
        cur = con.execute(
            "INSERT INTO gene_diseases(gene_mim, ordinal, label, phenotype_mim, inheritance, source, updated_at) "
            "VALUES (?,?,?,?,?, 'capture', ?)",
            (gene_mim, ordinal, label, mim, inheritance, _now()),
        )
        _log(con, user, "add_disease", gene_mim, mim, capture_id, None, {"label": label})
        return {"action": "added", "label": label, "disease_id": cur.lastrowid}
    if target["label"] != label or target["phenotype_mim"] != mim or target["inheritance"] != inheritance:
        con.execute(
            "UPDATE gene_diseases SET label = ?, phenotype_mim = ?, inheritance = ?, updated_at = ? WHERE id = ?",
            (label, mim, inheritance, _now(), target["id"]),
        )
        _log(con, user, "update_label", gene_mim, mim, capture_id,
             {"label": target["label"]}, {"label": label})
        return {"action": "updated", "label": label, "previous": target["label"], "disease_id": target["id"]}
    return {"action": "unchanged", "label": label, "disease_id": target["id"]}


def adopt_capture(parsed: dict, *, user: str, capture_id: str, path: Path | None = None) -> dict:
    """Apply a parsed capture to the database. Returns a change summary."""
    mim = parsed.get("mim") or ""
    entry_type = parsed.get("entry_type") or ""
    relationships = parsed.get("relationships") or []
    summary: dict = {"mim": mim, "entry_type": entry_type, "genes": [], "text": None, "reviews": 0}
    if not mim or entry_type == "moved_removed":
        raise ValueError("this entry cannot be adopted (no MIM, or moved/removed entry)")
    with _write_lock, connect(path, write=True) as con:
        if entry_type in {"gene", "gene_phenotype"}:
            gene = _gene_row(con, mim, (parsed.get("symbol") or "").split(",")[0].strip(), user, capture_id)
            listed = []
            changes = []
            for rel in relationships:
                if rel.get("gene_mim") and rel["gene_mim"] != mim:
                    continue
                result = _upsert_relationship(con, gene["gene_mim"], rel, user=user, capture_id=capture_id)
                listed.append(result["disease_id"])
                changes.append(result)
            if relationships:
                for row in con.execute(
                    "SELECT * FROM gene_diseases WHERE gene_mim = ? AND active = 1", (gene["gene_mim"],)
                ).fetchall():
                    if row["id"] in listed:
                        continue
                    kind = "not_in_omim" if row["phenotype_mim"] else "unmatched_label"
                    # Best guess of which OMIM-listed disease an old label meant.
                    guesses = sorted(
                        ((difflib.SequenceMatcher(None, disease_picks.name_key(row["label"]),
                                                  disease_picks.name_key(r.get("phenotype") or "")).ratio(),
                          format_label(r.get("phenotype"), r.get("phenotype_mim"), r.get("inheritance")))
                         for r in relationships),
                        reverse=True,
                    )
                    suggestion = guesses[0][1] if guesses and guesses[0][0] >= 0.6 else ""
                    exists_open = con.execute(
                        "SELECT 1 FROM review_items WHERE disease_id = ? AND kind = ? AND status = 'open'",
                        (row["id"], kind),
                    ).fetchone()
                    if not exists_open:
                        con.execute(
                            "INSERT INTO review_items(kind, gene_mim, disease_id, detail, created_at) VALUES (?,?,?,?,?)",
                            (kind, gene["gene_mim"], row["id"], json.dumps({
                                "gene": gene["gene_symbol"], "label": row["label"], "capture_id": capture_id,
                                "suggestion": suggestion,
                            }, ensure_ascii=False), _now()),
                        )
                # The gene page is now the reference: migration-time questions
                # about this gene are settled (OMIM-listed rows are confirmed;
                # unlisted rows carry the not_in_omim / unmatched_label item).
                con.execute(
                    "UPDATE review_items SET status='resolved', resolved_at=?, resolved_by=?, "
                    "resolution='gene_page_synced' WHERE status='open' AND gene_mim = ? AND kind IN "
                    "('renamed_match','missing_mim','slot_not_in_omim_disease')",
                    (_now(), user, gene["gene_mim"]),
                )
                if listed:
                    marks = ",".join("?" * len(listed))
                    con.execute(
                        "UPDATE review_items SET status='resolved', resolved_at=?, resolved_by=?, "
                        f"resolution='listed_by_omim' WHERE status='open' AND disease_id IN ({marks}) "
                        "AND kind IN ('not_in_omim','unmatched_label')",
                        (_now(), user, *listed),
                    )
                summary["reviews"] = con.execute(
                    "SELECT COUNT(*) FROM review_items WHERE status='open' AND gene_mim = ?",
                    (gene["gene_mim"],),
                ).fetchone()[0]
            _refresh_gene(con, gene["gene_mim"], done=True)
            summary["genes"].append({"gene_mim": gene["gene_mim"], "gene": gene["gene_symbol"], "changes": changes})
        if entry_type != "gene":
            own = [r for r in relationships if r.get("phenotype_mim") == mim]
            if entry_type == "gene_phenotype":
                own = []  # skeleton already synced above
            for rel in own:
                gene_mim = rel.get("gene_mim") or ""
                if not gene_mim:
                    continue
                gene = _gene_row(con, gene_mim, rel.get("gene") or "", user, capture_id)
                result = _upsert_relationship(con, gene["gene_mim"], rel, user=user, capture_id=capture_id)
                _refresh_gene(con, gene["gene_mim"], done=True)
                summary["genes"].append({
                    "gene_mim": gene["gene_mim"], "gene": gene["gene_symbol"], "changes": [result],
                })
            if not own and entry_type != "gene_phenotype":
                # No relationship row (e.g. collapsed table): attach text to
                # every gene already linked to this phenotype MIM.
                for row in con.execute(
                    "SELECT DISTINCT d.gene_mim, g.gene_symbol FROM gene_diseases d JOIN genes g "
                    "ON g.gene_mim = d.gene_mim WHERE d.phenotype_mim = ? AND d.active = 1", (mim,),
                ).fetchall():
                    _refresh_gene(con, row["gene_mim"], done=True)
                    summary["genes"].append({"gene_mim": row["gene_mim"], "gene": row["gene_symbol"], "changes": []})
            before = con.execute("SELECT * FROM phenotype_texts WHERE phenotype_mim = ?", (mim,)).fetchone()
            synopsis = parsed.get("synopsis") or ""
            synopsis_status = parsed.get("synopsis_status") or ""
            if not synopsis and synopsis_status == "not_loaded" and before and before["synopsis"]:
                # The synopsis was not expanded in this capture: keep the one
                # already adopted instead of wiping it.
                synopsis, synopsis_status = before["synopsis"], before["synopsis_status"]
                summary["kept_previous_synopsis"] = True
            values = {
                "title": parsed.get("title") or "",
                "description": parsed.get("description") or "",
                "description_source": parsed.get("description_source") or "",
                "synopsis": synopsis,
                "synopsis_status": synopsis_status,
                "omim_edit_date": parsed.get("omim_edit_date") or "",
                "content_hash": parsed.get("content_hash") or "",
            }
            con.execute(
                "INSERT INTO phenotype_texts(phenotype_mim, title, description, description_source, synopsis, "
                "synopsis_status, omim_edit_date, capture_id, content_hash, adopted_at, adopted_by) "
                "VALUES (:mim, :title, :description, :description_source, :synopsis, :synopsis_status, "
                ":omim_edit_date, :capture_id, :content_hash, :adopted_at, :user) "
                "ON CONFLICT(phenotype_mim) DO UPDATE SET title=excluded.title, description=excluded.description, "
                "description_source=excluded.description_source, synopsis=excluded.synopsis, "
                "synopsis_status=excluded.synopsis_status, omim_edit_date=excluded.omim_edit_date, "
                "capture_id=excluded.capture_id, content_hash=excluded.content_hash, "
                "adopted_at=excluded.adopted_at, adopted_by=excluded.adopted_by",
                {**values, "mim": mim, "capture_id": capture_id, "adopted_at": _now(), "user": user},
            )
            _log(con, user, "adopt_text", "", mim, capture_id,
                 dict(before) if before else None, {k: values[k] for k in ("title", "omim_edit_date", "content_hash")})
            summary["text"] = "replaced" if before else "added"
            # Migration reviews about this MIM are settled by adopting OMIM's text.
            con.execute(
                "UPDATE review_items SET status='resolved', resolved_at=?, resolved_by=?, resolution='adopted_capture' "
                "WHERE status='open' AND kind IN ('renamed_match','missing_mim') AND disease_id IN "
                "(SELECT id FROM gene_diseases WHERE phenotype_mim = ?)",
                (_now(), user, mim),
            )
        con.execute(
            "INSERT INTO adoptions(capture_id, mim, entry_type, title, adopted_at, adopted_by, summary) "
            "VALUES (?,?,?,?,?,?,?)",
            (capture_id, mim, entry_type, parsed.get("title") or "", _now(), user,
             json.dumps(summary, ensure_ascii=False)),
        )
    return summary


def adopted_capture_ids(path: Path | None = None) -> set[str]:
    if not exists(path):
        return set()
    with connect(path) as con:
        return {r["capture_id"] for r in con.execute("SELECT DISTINCT capture_id FROM adoptions")}


def gene_sync_map(path: Path | None = None) -> dict[str, str]:
    """gene MIM -> last time its OMIM gene page was adopted."""
    if not exists(path):
        return {}
    with connect(path) as con:
        return {r["mim"]: r["at"] for r in con.execute(
            "SELECT mim, MAX(adopted_at) AS at FROM adoptions "
            "WHERE entry_type IN ('gene','gene_phenotype') GROUP BY mim"
        )}


def open_reviews_by_gene(path: Path | None = None) -> dict[str, list[dict]]:
    if not exists(path):
        return {}
    out: dict[str, list[dict]] = {}
    with connect(path) as con:
        for r in con.execute(
            "SELECT r.id, r.kind, r.gene_mim, r.disease_id, r.detail, d.label AS current_label, "
            "d.active AS disease_active FROM review_items r LEFT JOIN gene_diseases d ON d.id = r.disease_id "
            "WHERE r.status = 'open' ORDER BY r.id"
        ):
            item = dict(r)
            try:
                item["detail"] = json.loads(item.get("detail") or "{}")
            except json.JSONDecodeError:
                item["detail"] = {}
            item["kind_label"] = REVIEW_KINDS.get(item["kind"], item["kind"])
            out.setdefault(item["gene_mim"], []).append(item)
    return out


def recently_adopted_mims(days: int = 7, path: Path | None = None) -> set[str]:
    """MIMs adopted within ``days`` (kept out of the "next page" suggestion)."""
    if not exists(path):
        return set()
    since = (datetime.now(timezone.utc) - timedelta(days=days)).isoformat(timespec="seconds")
    with connect(path) as con:
        return {r["mim"] for r in con.execute("SELECT DISTINCT mim FROM adoptions WHERE adopted_at >= ?", (since,))}


def list_adoptions(limit: int = 200, path: Path | None = None) -> list[dict]:
    if not exists(path):
        return []
    with connect(path) as con:
        rows = con.execute(
            "SELECT * FROM adoptions ORDER BY adopted_at DESC, id DESC LIMIT ?", (int(limit),)
        ).fetchall()
    out = []
    for r in rows:
        item = dict(r)
        try:
            summary = json.loads(item.get("summary") or "{}")
        except json.JSONDecodeError:
            summary = {}
        item["genes"] = [g.get("gene") for g in summary.get("genes") or [] if g.get("gene")]
        item["text"] = summary.get("text")
        item["reviews"] = summary.get("reviews") or 0
        item.pop("summary", None)
        out.append(item)
    return out


def adoption_blockers(parsed: dict) -> list[str]:
    """Reasons a capture should not be adopted in bulk (single adopt still allowed)."""
    reasons = []
    entry_type = parsed.get("entry_type") or ""
    if not parsed.get("mim") or entry_type == "moved_removed":
        reasons.append("已移除或搬移的 entry")
    if any("找不到 entry" in w or "不一致" in w for w in parsed.get("warnings") or []):
        reasons.append("頁面與 MIM 對不上")
    if entry_type in {"gene", "gene_phenotype"}:
        if not parsed.get("relationships"):
            reasons.append("基因頁沒有解析到疾病清單")
    else:
        if parsed.get("synopsis_status") == "not_loaded":
            reasons.append("Clinical Synopsis 沒有擷取到")
        if not parsed.get("description"):
            reasons.append("沒有 Description／Clinical Features")
    return reasons


def preview_adoption(parsed: dict, path: Path | None = None) -> dict:
    """What adopting would change, without writing (for the panel / page)."""
    mim = parsed.get("mim") or ""
    out = {"mim": mim, "text_exists": False, "same_content": False, "genes": []}
    if not exists(path):
        out["db_missing"] = True
        return out
    with connect(path) as con:
        text = con.execute("SELECT * FROM phenotype_texts WHERE phenotype_mim = ?", (mim,)).fetchone()
        if text:
            out["text_exists"] = True
            out["same_content"] = text["content_hash"] == parsed.get("content_hash")
            out["adopted_at"] = text["adopted_at"]
        for rel in parsed.get("relationships") or []:
            if parsed.get("entry_type") not in {"gene", "gene_phenotype"} and rel.get("phenotype_mim") != mim:
                continue
            gene_mim = rel.get("gene_mim") or mim
            label = format_label(rel.get("phenotype"), rel.get("phenotype_mim"), rel.get("inheritance"))
            rows = con.execute(
                "SELECT label FROM gene_diseases WHERE gene_mim = ? AND active = 1 AND phenotype_mim = ?",
                (gene_mim, rel.get("phenotype_mim") or ""),
            ).fetchall() if rel.get("phenotype_mim") else []
            fills = ""
            if not rows:
                # An old MIM-less line with the same name gets this MIM.
                for r in con.execute(
                    "SELECT label FROM gene_diseases WHERE gene_mim = ? AND active = 1 AND phenotype_mim = ''",
                    (gene_mim,),
                ):
                    if disease_picks.name_key(r["label"]) == disease_picks.name_key(label):
                        fills = r["label"]
                        break
            out["genes"].append({
                "gene": rel.get("gene") or "", "gene_mim": gene_mim, "label": label,
                "current_labels": [r["label"] for r in rows], "fills_mim_for": fills,
            })
    return out


# --------------------------------------------------------------------------
# Review list
# --------------------------------------------------------------------------

def list_reviews(status: str = "open", kind: str = "", limit: int = 500, path: Path | None = None) -> dict:
    with connect(path) as con:
        params: list = [status]
        sql = ("SELECT r.*, g.gene_symbol, d.label AS current_label, d.active AS disease_active "
               "FROM review_items r LEFT JOIN genes g ON g.gene_mim = r.gene_mim "
               "LEFT JOIN gene_diseases d ON d.id = r.disease_id WHERE r.status = ?")
        if kind:
            sql += " AND r.kind = ?"
            params.append(kind)
        sql += " ORDER BY r.kind, g.gene_symbol, r.id LIMIT ?"
        params.append(int(limit))
        items = []
        for r in con.execute(sql, params):
            item = dict(r)
            item["detail"] = json.loads(item.get("detail") or "{}")
            item["kind_label"] = REVIEW_KINDS.get(item["kind"], item["kind"])
            items.append(item)
        counts = {r["kind"]: r["n"] for r in con.execute(
            "SELECT kind, COUNT(*) AS n FROM review_items WHERE status = 'open' GROUP BY kind"
        )}
    return {"items": items, "open_counts": counts}


def resolve_review(item_id: int, action: str, *, user: str, path: Path | None = None) -> dict:
    """``keep`` closes the item; ``remove`` also deactivates its disease row."""
    if action not in {"keep", "remove"}:
        raise ValueError("action must be keep or remove")
    with _write_lock, connect(path, write=True) as con:
        item = con.execute("SELECT * FROM review_items WHERE id = ?", (int(item_id),)).fetchone()
        if not item:
            raise KeyError(item_id)
        if item["status"] != "open":
            return {"id": item_id, "status": item["status"]}
        if action == "remove":
            if not item["disease_id"]:
                raise ValueError("this item has no disease row to remove")
            disease = con.execute("SELECT * FROM gene_diseases WHERE id = ?", (item["disease_id"],)).fetchone()
            con.execute("UPDATE gene_diseases SET active = 0, updated_at = ? WHERE id = ?", (_now(), item["disease_id"]))
            _log(con, user, "remove_disease", item["gene_mim"], disease["phenotype_mim"] if disease else "",
                 "", {"label": disease["label"] if disease else ""}, None)
            _refresh_gene(con, item["gene_mim"], done=False)
        con.execute(
            "UPDATE review_items SET status='resolved', resolved_at=?, resolved_by=?, resolution=? WHERE id = ?",
            (_now(), user, action, int(item_id)),
        )
    return {"id": item_id, "status": "resolved", "resolution": action}


def status(path: Path | None = None) -> dict:
    path = path or db_path()
    if not exists(path):
        return {"enabled": False, "path": str(path)}
    with connect(path) as con:
        meta = {r["key"]: r["value"] for r in con.execute("SELECT key, value FROM meta")}
        counts = {
            "genes": con.execute("SELECT COUNT(*) FROM genes").fetchone()[0],
            "diseases": con.execute("SELECT COUNT(*) FROM gene_diseases WHERE active = 1").fetchone()[0],
            "phenotype_texts": con.execute("SELECT COUNT(*) FROM phenotype_texts").fetchone()[0],
            "open_reviews": con.execute("SELECT COUNT(*) FROM review_items WHERE status='open'").fetchone()[0],
            "changes": con.execute("SELECT COUNT(*) FROM change_log").fetchone()[0],
        }
    return {"enabled": True, "path": str(path), "meta": meta, "counts": counts}


# --------------------------------------------------------------------------
# Export to the workbook layout
# --------------------------------------------------------------------------

def export_xlsx(path: Path | None = None) -> bytes:
    """OMIM.xlsx in the familiar A–U layout, plus an overflow sheet."""
    import openpyxl

    _by_id, _by_gene = load_rows(path, slot_count=10_000)
    with connect(path) as con:
        order = [r["gene_mim"] for r in con.execute("SELECT gene_mim FROM genes ORDER BY row_order")]
    rows = {}
    for row in list(_by_id.values()) + list(_by_gene.values()):
        key = row["OMIM_id"] or f"symbol:{row['gene_symbol']}"
        rows.setdefault(key, row)
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "工作表1"
    headers = ["OMIM_id", "gene_symbol", "OMIM_disease", "Inheritance"] + [
        f"Disease{i}" for i in range(1, MAX_RENDERED_SLOTS + 1)
    ] + ["Done"]
    ws.append(headers)
    overflow = wb.create_sheet("Overflow")
    overflow.append(["OMIM_id", "gene_symbol", "slot", "reason", "full_text"])
    for gene_mim in order:
        row = rows.get(gene_mim)
        if not row:
            continue
        values = [
            int(row["OMIM_id"]) if str(row["OMIM_id"]).isdigit() else row["OMIM_id"],
            row["gene_symbol"], row["OMIM_disease"], row["Inheritance"],
        ]
        n = 1
        while row.get(f"Disease{n}"):
            n += 1
        for i in range(1, MAX_RENDERED_SLOTS + 1):
            text = row.get(f"Disease{i}") or None
            if text and len(text) > EXCEL_CELL_LIMIT:
                overflow.append([values[0], row["gene_symbol"], i, "cell_over_32767", text[:EXCEL_CELL_LIMIT]])
                text = text[:EXCEL_CELL_LIMIT - 40] + "\n…（全文見 Overflow 工作表）"
            values.append(text)
        for i in range(MAX_RENDERED_SLOTS + 1, n):
            overflow.append([values[0], row["gene_symbol"], i, "over_16_diseases",
                             (row.get(f"Disease{i}") or "")[:EXCEL_CELL_LIMIT]])
        done = row.get("Done") or ""
        try:
            values.append(datetime.strptime(done[:10], "%Y-%m-%d").date() if done else None)
        except ValueError:
            values.append(done)
        ws.append(values)
        if values[-1] and not isinstance(values[-1], str):
            ws.cell(ws.max_row, len(headers)).number_format = "yyyy/mm/dd"
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()
