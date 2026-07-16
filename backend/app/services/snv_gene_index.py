"""SQLite gene index for complete SNV/Indel TSV search.

The main reviewer screen uses the compact snv_indel.review.tsv, while
gene search needs the full 03_acmg raw TSV. WGS raw TSVs can be
1-2 GB, so tertiary jobs prebuild this tiny lookup database once and the
web request fetches only rows for requested genes.
"""
from __future__ import annotations

import json
import os
import sqlite3
import time
from pathlib import Path

from . import panel_deadzone
from .snv_rows import is_reportable_raw_row


INDEX_NAME = "snv_gene_index.sqlite"
SCHEMA_VERSION = "3"


def index_path_for(raw_tsv: Path) -> Path:
    return raw_tsv.with_name(INDEX_NAME)


def _log_perf(event: str, started: float, **fields) -> None:
    parts = [f"[perf] {event}", f"elapsed={time.perf_counter() - started:.3f}s"]
    parts.extend(f"{key}={value}" for key, value in fields.items())
    print(" ".join(parts), flush=True)


def _source_meta(raw_tsv: Path) -> dict[str, str]:
    raw_tsv = Path(raw_tsv).resolve()
    st = raw_tsv.stat()
    return {
        "schema_version": SCHEMA_VERSION,
        "source_path": str(raw_tsv),
        "source_mtime_ns": str(st.st_mtime_ns),
        "source_size": str(st.st_size),
    }


def _read_meta(conn: sqlite3.Connection) -> dict[str, str]:
    try:
        rows = conn.execute("SELECT key, value FROM meta").fetchall()
    except sqlite3.Error:
        return {}
    return {str(k): str(v) for k, v in rows}


def is_current(raw_tsv: Path, index_path: Path | None = None) -> bool:
    index_path = index_path or index_path_for(raw_tsv)
    if not raw_tsv.is_file() or not index_path.is_file():
        return False
    expected = _source_meta(raw_tsv)
    try:
        with sqlite3.connect(index_path) as conn:
            meta = _read_meta(conn)
            return all(meta.get(key) == value for key, value in expected.items())
    except sqlite3.Error:
        return False


def build_index(raw_tsv: Path, index_path: Path | None = None) -> Path:
    """Build or replace the per-sample gene index for *raw_tsv*."""
    started = time.perf_counter()
    raw_tsv = Path(raw_tsv)
    if not raw_tsv.is_file():
        raise FileNotFoundError(raw_tsv)
    index_path = index_path or index_path_for(raw_tsv)
    tmp = index_path.with_suffix(index_path.suffix + ".tmp")
    try:
        tmp.unlink()
    except FileNotFoundError:
        pass
    tmp.parent.mkdir(parents=True, exist_ok=True)

    scanned = indexed = 0
    with sqlite3.connect(tmp) as conn:
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")
        conn.execute("PRAGMA temp_store=MEMORY")
        conn.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
        conn.execute(
            "CREATE TABLE variants ("
            "gene TEXT NOT NULL, "
            "variant_id TEXT NOT NULL, "
            "offset INTEGER NOT NULL, "
            "length INTEGER NOT NULL)"
        )
        source_meta = _source_meta(raw_tsv)
        batch: list[tuple[str, str, int, int]] = []
        with raw_tsv.open("rb") as fh:
            header_line = fh.readline()
            header = header_line.decode("utf-8").rstrip("\n\r").split("\t")
            source_meta["header_json"] = json.dumps(
                header, ensure_ascii=False, separators=(",", ":")
            )
            conn.executemany(
                "INSERT INTO meta(key, value) VALUES (?, ?)",
                sorted(source_meta.items()),
            )
            while True:
                offset = fh.tell()
                line = fh.readline()
                if not line:
                    break
                scanned += 1
                length = len(line)
                values = line.decode("utf-8").rstrip("\n\r").split("\t")
                row = dict(zip(header, values))
                if not is_reportable_raw_row(row):
                    continue
                gene, _hgnc_id = panel_deadzone.canonical_gene_symbol(
                    row.get("GENE", ""),
                    row.get("HGNC_ID", ""),
                )
                if not gene:
                    continue
                chrom = (row.get("CHROM") or "").strip()
                pos = (row.get("POS") or "").strip()
                ref = (row.get("REF") or "").strip()
                alt = (row.get("ALT") or "").strip()
                if not (chrom and pos and ref and alt):
                    continue
                variant_id = "-".join((chrom, pos, ref, alt))
                batch.append((
                    gene.upper(),
                    variant_id,
                    offset,
                    length,
                ))
                indexed += 1
                if len(batch) >= 10000:
                    conn.executemany(
                        "INSERT INTO variants(gene, variant_id, offset, length) "
                        "VALUES (?, ?, ?, ?)",
                        batch,
                    )
                    batch.clear()
        if batch:
            conn.executemany(
                "INSERT INTO variants(gene, variant_id, offset, length) "
                "VALUES (?, ?, ?, ?)",
                batch,
            )
        conn.execute("CREATE INDEX idx_variants_gene ON variants(gene)")
        conn.execute("CREATE INDEX idx_variants_id ON variants(variant_id)")
        conn.commit()
    os.replace(tmp, index_path)
    _log_perf(
        "snv_gene_index.build",
        started,
        tsv=raw_tsv.name,
        index=index_path.name,
        scanned=scanned,
        indexed=indexed,
    )
    return index_path


def query_rows(
    raw_tsv: Path,
    genes: list[str],
    index_path: Path | None = None,
) -> list[dict[str, str]] | None:
    """Return raw TSV rows for genes, or None when the index is unavailable."""
    started = time.perf_counter()
    raw_tsv = Path(raw_tsv)
    index_path = Path(index_path) if index_path else index_path_for(raw_tsv)
    canonical = {
        panel_deadzone.canonical_gene_symbol(g)[0].upper()
        for g in genes
        if g and g.strip()
    }
    canonical.discard("")
    if not canonical:
        return []
    if not is_current(raw_tsv, index_path):
        _log_perf(
            "snv_gene_index.query",
            started,
            status="missing_or_stale",
            index=index_path.name,
            genes=len(canonical),
        )
        return None

    placeholders = ",".join("?" for _ in canonical)
    with sqlite3.connect(index_path) as conn:
        meta = _read_meta(conn)
        header = json.loads(meta.get("header_json") or "[]")
        rows = conn.execute(
            f"SELECT offset, length FROM variants WHERE gene IN ({placeholders})",
            sorted(canonical),
        ).fetchall()
    out = []
    with raw_tsv.open("rb") as fh:
        for offset, length in rows:
            fh.seek(int(offset))
            values = fh.read(int(length)).decode("utf-8").rstrip("\n\r").split("\t")
            out.append(dict(zip(header, values)))
    _log_perf(
        "snv_gene_index.query",
        started,
        status="hit",
        index=index_path.name,
        genes=len(canonical),
        rows=len(out),
    )
    return out


def query_rows_by_ids(
    raw_tsv: Path,
    variant_ids: set[str] | list[str],
    index_path: Path | None = None,
) -> list[dict[str, str]] | None:
    """Return raw TSV rows for exact variant ids, or None when unavailable."""
    started = time.perf_counter()
    raw_tsv = Path(raw_tsv)
    index_path = Path(index_path) if index_path else index_path_for(raw_tsv)
    wanted = sorted({str(vid).strip() for vid in variant_ids if str(vid).strip()})
    if not wanted:
        return []
    if not is_current(raw_tsv, index_path):
        _log_perf(
            "snv_gene_index.query_ids",
            started,
            status="missing_or_stale",
            index=index_path.name,
            ids=len(wanted),
        )
        return None

    placeholders = ",".join("?" for _ in wanted)
    with sqlite3.connect(index_path) as conn:
        meta = _read_meta(conn)
        header = json.loads(meta.get("header_json") or "[]")
        rows = conn.execute(
            f"SELECT offset, length FROM variants WHERE variant_id IN ({placeholders})",
            wanted,
        ).fetchall()
    out = []
    with raw_tsv.open("rb") as fh:
        for offset, length in rows:
            fh.seek(int(offset))
            values = fh.read(int(length)).decode("utf-8").rstrip("\n\r").split("\t")
            out.append(dict(zip(header, values)))
    _log_perf(
        "snv_gene_index.query_ids",
        started,
        status="hit",
        index=index_path.name,
        ids=len(wanted),
        rows=len(out),
    )
    return out
