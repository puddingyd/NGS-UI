"""Sparse post-processing overlay for immutable pipeline SNV TSVs."""
from __future__ import annotations

import csv
import json
import os
import sqlite3
from pathlib import Path
from typing import Iterable


OVERLAY_NAME = "snv_annotations.sqlite"
SCHEMA_VERSION = "1"
KEY_FIELDS = (
    "CHROM", "POS", "REF", "ALT", "GENE", "TRANSCRIPT",
    "HGVS_C", "HGVS_P", "CONSEQUENCE",
)


def row_key(row: dict[str, str]) -> str:
    return json.dumps(
        [str(row.get(field) or "") for field in KEY_FIELDS],
        ensure_ascii=False,
        separators=(",", ":"),
    )


def source_signature(path: Path) -> dict[str, str]:
    path = Path(path).resolve()
    stat = path.stat()
    return {
        "source_path": str(path),
        "source_mtime_ns": str(stat.st_mtime_ns),
        "source_size": str(stat.st_size),
    }


def _read_meta(conn: sqlite3.Connection) -> dict[str, str]:
    try:
        rows = conn.execute("SELECT key, value FROM meta").fetchall()
    except sqlite3.Error:
        return {}
    return {str(key): str(value) for key, value in rows}


def is_current(raw_tsv: Path, overlay_path: Path) -> bool:
    if not Path(raw_tsv).is_file() or not Path(overlay_path).is_file():
        return False
    expected = {"schema_version": SCHEMA_VERSION, **source_signature(Path(raw_tsv))}
    try:
        with sqlite3.connect(overlay_path) as conn:
            meta = _read_meta(conn)
    except sqlite3.Error:
        return False
    return all(meta.get(key) == value for key, value in expected.items())


def build_overlay(raw_tsv: Path, annotated_tsv: Path, overlay_path: Path) -> Path:
    """Build an overlay while allowing *annotated_tsv* to omit raw rows.

    Post-processors preserve row order.  Legacy UI TSVs may have removed ``*``
    or alternate-contig rows, so the comparison advances through raw rows until
    each annotated key is found.  It never has to load a multi-GB TSV in memory.
    """
    raw_tsv = Path(raw_tsv)
    annotated_tsv = Path(annotated_tsv)
    overlay_path = Path(overlay_path)
    if not raw_tsv.is_file():
        raise FileNotFoundError(raw_tsv)
    if not annotated_tsv.is_file():
        raise FileNotFoundError(annotated_tsv)
    overlay_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = overlay_path.with_suffix(overlay_path.suffix + ".tmp")
    try:
        tmp.unlink()
    except FileNotFoundError:
        pass

    raw_count = annotated_count = overlay_count = skipped_raw = 0
    with raw_tsv.open("r", encoding="utf-8", newline="") as raw_fh, \
            annotated_tsv.open("r", encoding="utf-8", newline="") as annotated_fh, \
            sqlite3.connect(tmp) as conn:
        raw_reader = csv.DictReader(raw_fh, delimiter="\t")
        annotated_reader = csv.DictReader(annotated_fh, delimiter="\t")
        raw_fields = raw_reader.fieldnames or []
        annotated_fields = annotated_reader.fieldnames or []
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")
        conn.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
        conn.execute(
            "CREATE TABLE annotations ("
            "row_key TEXT PRIMARY KEY, payload_json TEXT NOT NULL)"
        )
        batch: list[tuple[str, str]] = []
        raw_iter = iter(raw_reader)
        raw_row = next(raw_iter, None)
        if raw_row is not None:
            raw_count += 1
        for annotated_row in annotated_reader:
            annotated_count += 1
            wanted_key = row_key(annotated_row)
            while raw_row is not None and row_key(raw_row) != wanted_key:
                skipped_raw += 1
                raw_row = next(raw_iter, None)
                if raw_row is not None:
                    raw_count += 1
            if raw_row is None:
                raise ValueError(
                    "annotated TSV row order/key does not match raw TSV at "
                    f"annotated row {annotated_count}: {wanted_key}"
                )
            payload = {
                field: str(annotated_row.get(field) or "")
                for field in annotated_fields
                if str(annotated_row.get(field) or "") != str(raw_row.get(field) or "")
            }
            if payload:
                batch.append((wanted_key, json.dumps(payload, ensure_ascii=False)))
                overlay_count += 1
            if len(batch) >= 10000:
                conn.executemany(
                    "INSERT OR REPLACE INTO annotations(row_key, payload_json) VALUES (?, ?)",
                    batch,
                )
                batch.clear()
            raw_row = next(raw_iter, None)
            if raw_row is not None:
                raw_count += 1
        if batch:
            conn.executemany(
                "INSERT OR REPLACE INTO annotations(row_key, payload_json) VALUES (?, ?)",
                batch,
            )
        while raw_row is not None:
            skipped_raw += 1
            raw_row = next(raw_iter, None)
            if raw_row is not None:
                raw_count += 1

        meta = {
            "schema_version": SCHEMA_VERSION,
            **source_signature(raw_tsv),
            "raw_fields_json": json.dumps(raw_fields, ensure_ascii=False),
            "annotated_fields_json": json.dumps(annotated_fields, ensure_ascii=False),
            "overlay_fields_json": json.dumps(
                sorted(set(annotated_fields) - set(raw_fields)), ensure_ascii=False
            ),
            "raw_rows": str(raw_count),
            "annotated_rows": str(annotated_count),
            "overlay_rows": str(overlay_count),
            "skipped_raw_rows": str(skipped_raw),
        }
        conn.executemany(
            "INSERT INTO meta(key, value) VALUES (?, ?)", sorted(meta.items())
        )
        conn.commit()
    os.replace(tmp, overlay_path)
    return overlay_path


class OverlayReader:
    """Reusable read connection that merges sparse fields onto raw rows."""

    def __init__(self, raw_tsv: Path, overlay_path: Path | None):
        self.raw_tsv = Path(raw_tsv)
        self.path = Path(overlay_path) if overlay_path else None
        self.conn: sqlite3.Connection | None = None
        self.fields: list[str] = []

    def __enter__(self) -> "OverlayReader":
        if self.path and is_current(self.raw_tsv, self.path):
            self.conn = sqlite3.connect(self.path)
            meta = _read_meta(self.conn)
            try:
                self.fields = list(json.loads(meta.get("annotated_fields_json") or "[]"))
            except (TypeError, json.JSONDecodeError):
                self.fields = []
        return self

    def __exit__(self, *_args) -> None:
        if self.conn is not None:
            self.conn.close()
            self.conn = None

    @property
    def active(self) -> bool:
        return self.conn is not None

    def apply(self, row: dict[str, str]) -> dict[str, str]:
        if self.conn is None:
            return row
        result = self.conn.execute(
            "SELECT payload_json FROM annotations WHERE row_key = ?",
            (row_key(row),),
        ).fetchone()
        if not result:
            return row
        merged = dict(row)
        try:
            payload = json.loads(result[0])
        except (TypeError, json.JSONDecodeError):
            return row
        if isinstance(payload, dict):
            merged.update({str(k): str(v) for k, v in payload.items()})
        return merged

    def apply_many(self, rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
        values = list(rows)
        if self.conn is None or not values:
            return values
        payloads: dict[str, dict] = {}
        keys = list(dict.fromkeys(row_key(row) for row in values))
        for start in range(0, len(keys), 800):
            chunk = keys[start:start + 800]
            placeholders = ",".join("?" for _ in chunk)
            for key, raw_payload in self.conn.execute(
                f"SELECT row_key, payload_json FROM annotations "
                f"WHERE row_key IN ({placeholders})",
                chunk,
            ):
                try:
                    payload = json.loads(raw_payload)
                except (TypeError, json.JSONDecodeError):
                    continue
                if isinstance(payload, dict):
                    payloads[str(key)] = {
                        str(field): str(value) for field, value in payload.items()
                    }
        out: list[dict[str, str]] = []
        for row in values:
            payload = payloads.get(row_key(row))
            if not payload:
                out.append(row)
                continue
            merged = dict(row)
            merged.update(payload)
            out.append(merged)
        return out


def overlay_signature(raw_tsv: Path, overlay_path: Path | None) -> dict[str, object]:
    if not overlay_path or not Path(overlay_path).is_file():
        return {"exists": False}
    path = Path(overlay_path)
    stat = path.stat()
    return {
        "exists": True,
        "path": str(path),
        "mtime_ns": stat.st_mtime_ns,
        "size": stat.st_size,
        "current": is_current(Path(raw_tsv), path),
    }
