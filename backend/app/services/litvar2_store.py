"""Local LitVar2 bulk-data index and variant lookup helpers.

The official LitVar2 dataset is a large gzip-compressed JSON export.  This
module streams either a JSON array or newline-delimited JSON without loading
the file into memory, and builds a compact SQLite index containing only the
fields needed by the reviewer card:

* exact LitVar2 variant ID / rsID
* gene + HGVS aliases for rsID-less fallback
* PMID count and the first five PMIDs in source order
* the bulk file's HTTP Last-Modified date

The live database is replaced atomically only after the new index passes
``PRAGMA quick_check``.  Existing tertiary readers therefore keep using a
complete old inode while an update is being downloaded or built.
"""
from __future__ import annotations

import email.utils
import gzip
import hashlib
import json
import os
import re
import shutil
import sqlite3
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Iterator, TextIO

from ..config import (
    LITVAR2_BULK_PATH,
    LITVAR2_BULK_URL,
    LITVAR2_DB,
    LITVAR2_DIR,
    LITVAR2_MANIFEST_PATH,
)


SCHEMA_VERSION = "1"
CHUNK_SIZE = 1024 * 1024
DOWNLOAD_RETRIES = 4
DOWNLOAD_TIMEOUT_SECONDS = 120
_RSID_RE = re.compile(r"(?i)\brs(\d+)\b")
_HGVS_RE = re.compile(
    r"(?i)(?:^|:)(?:c|g|m|n|p|r)\.[A-Za-z0-9_*?+>\-=()\[\];:.]+$"
)
Progress = Callable[[str, dict[str, object]], None]


class LitVar2Error(RuntimeError):
    """Raised when a bulk download or index build is not publishable."""


def _notify(progress: Progress | None, step: str, **fields: object) -> None:
    if progress:
        progress(step, fields)


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    os.replace(tmp, path)


def _read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _flatten_strings(value) -> Iterator[str]:
    if value is None:
        return
    if isinstance(value, str):
        text = value.strip()
        if text:
            yield text
        return
    if isinstance(value, dict):
        for nested in value.values():
            yield from _flatten_strings(nested)
        return
    if isinstance(value, (list, tuple, set)):
        for nested in value:
            yield from _flatten_strings(nested)
        return
    text = str(value).strip()
    if text:
        yield text


def _read_more(handle: TextIO, buffer: str, cursor: int) -> tuple[str, int, bool]:
    if cursor:
        buffer = buffer[cursor:]
        cursor = 0
    chunk = handle.read(CHUNK_SIZE)
    return buffer + chunk, cursor, not bool(chunk)


def iter_json_values(handle: TextIO) -> Iterator[dict]:
    """Yield records from a JSON array or sequential/JSONL objects.

    ``json.JSONDecoder.raw_decode`` is retried as chunks arrive.  Memory use is
    bounded by the largest individual LitVar2 record rather than the complete
    decompressed dataset.
    """
    decoder = json.JSONDecoder()
    buffer = ""
    cursor = 0
    eof = False
    mode: str | None = None

    while True:
        while True:
            while cursor < len(buffer) and buffer[cursor].isspace():
                cursor += 1
            if cursor < len(buffer) or eof:
                break
            buffer, cursor, eof = _read_more(handle, buffer, cursor)

        if mode is None:
            if cursor >= len(buffer):
                return
            if buffer[cursor] == "[":
                mode = "array"
                cursor += 1
            else:
                mode = "values"

        while True:
            while cursor < len(buffer) and (
                buffer[cursor].isspace()
                or (mode == "array" and buffer[cursor] == ",")
            ):
                cursor += 1
            if cursor < len(buffer) or eof:
                break
            buffer, cursor, eof = _read_more(handle, buffer, cursor)

        if mode == "array" and cursor < len(buffer) and buffer[cursor] == "]":
            return
        if cursor >= len(buffer) and eof:
            return

        try:
            value, end = decoder.raw_decode(buffer, cursor)
        except json.JSONDecodeError as exc:
            if eof:
                raise LitVar2Error(
                    f"invalid/truncated LitVar2 JSON near character {exc.pos}: {exc.msg}"
                ) from exc
            buffer, cursor, eof = _read_more(handle, buffer, cursor)
            continue

        cursor = end
        if isinstance(value, dict):
            yield value

        if cursor > CHUNK_SIZE * 2:
            buffer = buffer[cursor:]
            cursor = 0


def iter_bulk_records(path: Path) -> Iterator[dict]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from iter_json_values(handle)


def normalize_rsid(value: str) -> str:
    match = _RSID_RE.search(str(value or ""))
    return f"rs{match.group(1)}" if match else ""


def normalize_gene(value: str) -> str:
    return str(value or "").strip().upper()


def normalize_hgvs(value: str) -> str:
    text = urllib.parse.unquote(str(value or "")).strip()
    if ":" in text:
        prefix, suffix = text.rsplit(":", 1)
        if re.search(r"(?i)(?:^|\.)(?:c|g|m|n|p|r)$", suffix.split(".", 1)[0]):
            text = suffix
        elif re.match(r"(?i)^(?:c|g|m|n|p|r)\.", suffix):
            text = suffix
    text = re.sub(r"\s+", "", text)
    return text.upper()


def _looks_hgvs(value: str) -> bool:
    text = urllib.parse.unquote(str(value or "")).strip()
    return bool(_HGVS_RE.search(text.replace(" ", "")))


def _human_record(record: dict) -> bool:
    tax_ids = set(
        _flatten_strings(record.get("data_tax_id"))
    ) | set(_flatten_strings(record.get("data_species_id")))
    return not tax_ids or "9606" in tax_ids


def _record_rsids(record: dict) -> list[str]:
    values: list[str] = []
    for raw in (
        record.get("rsid"),
        record.get("data_snp_id"),
        record.get("_id"),
    ):
        for text in _flatten_strings(raw):
            rsid = normalize_rsid(text)
            if not rsid and text.isdigit():
                rsid = f"rs{text}"
            if rsid and rsid not in values:
                values.append(rsid)
    return values


def _record_genes(record: dict) -> list[str]:
    values: list[str] = []
    for text in _flatten_strings(record.get("gene")):
        gene = normalize_gene(text)
        if gene and gene not in values:
            values.append(gene)
    return values


def _record_hgvs(record: dict) -> list[str]:
    values: list[str] = []
    for field in ("hgvs", "hgvs_prot", "name", "all_hgvs", "synonyms"):
        for text in _flatten_strings(record.get(field)):
            if not _looks_hgvs(text):
                continue
            hgvs = normalize_hgvs(text)
            if hgvs and hgvs not in values:
                values.append(hgvs)
    return values


def _record_pmids(record: dict) -> list[str]:
    values: list[str] = []
    for raw in _flatten_strings(record.get("pmids")):
        for token in re.findall(r"\d+", raw):
            if token and token not in values:
                values.append(token)
    return values


def _to_int(value, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _database_meta(conn: sqlite3.Connection) -> dict[str, str]:
    try:
        rows = conn.execute("SELECT key, value FROM meta").fetchall()
    except sqlite3.Error:
        return {}
    return {str(key): str(value) for key, value in rows}


def database_metadata(
    path: Path = LITVAR2_DB,
    *,
    verify: bool = False,
) -> dict[str, object]:
    path = Path(path)
    if not path.is_file():
        return {"available": False, "path": str(path)}
    try:
        uri = f"file:{urllib.parse.quote(str(path.resolve()))}?mode=ro"
        with sqlite3.connect(uri, uri=True) as conn:
            meta = _database_meta(conn)
            tables = {
                str(row[0]) for row in conn.execute(
                    "SELECT name FROM sqlite_master WHERE type = 'table'"
                )
            }
            quick = conn.execute("PRAGMA quick_check").fetchone() if verify else ("ok",)
    except sqlite3.Error as exc:
        return {
            "available": False,
            "path": str(path),
            "error": str(exc),
        }
    valid = (
        meta.get("schema_version") == SCHEMA_VERSION
        and {"meta", "variants", "hgvs_aliases"} <= tables
        and bool(quick)
        and quick[0] == "ok"
    )
    return {
        "available": valid,
        "path": str(path),
        "dataset_date": meta.get("dataset_date", ""),
        "record_count": _to_int(meta.get("record_count")),
        "alias_count": _to_int(meta.get("alias_count")),
        "source_last_modified": meta.get("source_last_modified", ""),
        "source_sha256": meta.get("source_sha256", ""),
        **({"error": "schema/database validation failed"} if not valid else {}),
    }


def build_database(
    bulk_path: Path,
    out_path: Path,
    *,
    dataset_date: str,
    source_url: str = "",
    source_last_modified: str = "",
    source_sha256: str = "",
    progress: Progress | None = None,
) -> dict[str, object]:
    """Stream *bulk_path* and build a new slim SQLite database at *out_path*."""
    bulk_path = Path(bulk_path)
    out_path = Path(out_path)
    if not bulk_path.is_file():
        raise FileNotFoundError(bulk_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.unlink(missing_ok=True)

    processed = indexed = alias_count = 0
    started = time.monotonic()
    _notify(progress, "building-index", processed=0, indexed=0)

    with sqlite3.connect(out_path) as conn:
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")
        conn.execute("PRAGMA temp_store=FILE")
        conn.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
        conn.execute(
            "CREATE TABLE variants ("
            "id INTEGER PRIMARY KEY,"
            "litvar_id TEXT NOT NULL,"
            "rsid TEXT NOT NULL,"
            "primary_gene TEXT NOT NULL,"
            "preferred_hgvs TEXT NOT NULL,"
            "pmids_count INTEGER NOT NULL,"
            "top_pmids_json TEXT NOT NULL)"
        )
        conn.execute(
            "CREATE TABLE hgvs_aliases ("
            "gene TEXT NOT NULL,"
            "hgvs TEXT NOT NULL,"
            "variant_id INTEGER NOT NULL)"
        )

        variant_batch: list[tuple] = []
        alias_batch: list[tuple[str, str, int]] = []
        next_id = 1
        for record in iter_bulk_records(bulk_path):
            processed += 1
            if not _human_record(record):
                continue
            litvar_id = str(record.get("_id") or "").strip()
            if not litvar_id:
                continue
            rsids = _record_rsids(record)
            genes = _record_genes(record)
            hgvs_values = _record_hgvs(record)
            pmids = _record_pmids(record)
            pmids_count = max(_to_int(record.get("pmids_count")), len(pmids))
            preferred_hgvs = normalize_hgvs(record.get("hgvs") or "")
            if not preferred_hgvs and hgvs_values:
                preferred_hgvs = hgvs_values[0]
            variant_id = next_id
            next_id += 1
            variant_batch.append((
                variant_id,
                litvar_id,
                rsids[0] if rsids else "",
                genes[0] if genes else "",
                preferred_hgvs,
                pmids_count,
                json.dumps(pmids[:5], separators=(",", ":")),
            ))
            indexed += 1

            seen_aliases: set[tuple[str, str]] = set()
            for gene in genes or [""]:
                for hgvs in hgvs_values:
                    key = (gene, hgvs)
                    if key in seen_aliases:
                        continue
                    seen_aliases.add(key)
                    alias_batch.append((gene, hgvs, variant_id))
                    alias_count += 1

            if len(variant_batch) >= 10000:
                conn.executemany(
                    "INSERT INTO variants VALUES (?, ?, ?, ?, ?, ?, ?)",
                    variant_batch,
                )
                variant_batch.clear()
            if len(alias_batch) >= 50000:
                conn.executemany(
                    "INSERT INTO hgvs_aliases VALUES (?, ?, ?)",
                    alias_batch,
                )
                alias_batch.clear()
            if processed % 100000 == 0:
                conn.commit()
                _notify(
                    progress,
                    "building-index",
                    processed=processed,
                    indexed=indexed,
                    aliases=alias_count,
                    elapsed_seconds=round(time.monotonic() - started, 1),
                )

        if variant_batch:
            conn.executemany(
                "INSERT INTO variants VALUES (?, ?, ?, ?, ?, ?, ?)",
                variant_batch,
            )
        if alias_batch:
            conn.executemany(
                "INSERT INTO hgvs_aliases VALUES (?, ?, ?)",
                alias_batch,
            )
        if not indexed:
            raise LitVar2Error("LitVar2 bulk file produced zero human variant records")

        _notify(
            progress,
            "building-indexes",
            processed=processed,
            indexed=indexed,
            aliases=alias_count,
        )
        conn.execute("CREATE UNIQUE INDEX variants_litvar_id_idx ON variants(litvar_id)")
        conn.execute("CREATE INDEX variants_rsid_idx ON variants(rsid)")
        conn.execute(
            "CREATE INDEX hgvs_alias_lookup_idx ON hgvs_aliases(gene, hgvs)"
        )
        conn.execute("CREATE INDEX hgvs_alias_variant_idx ON hgvs_aliases(variant_id)")
        meta = {
            "schema_version": SCHEMA_VERSION,
            "dataset_date": dataset_date,
            "source_url": source_url,
            "source_last_modified": source_last_modified,
            "source_sha256": source_sha256,
            "record_count": str(indexed),
            "source_record_count": str(processed),
            "alias_count": str(alias_count),
            "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
        conn.executemany(
            "INSERT INTO meta(key, value) VALUES (?, ?)",
            sorted(meta.items()),
        )
        conn.commit()
        conn.execute("ANALYZE")
        conn.commit()

    _notify(progress, "validating-index", indexed=indexed, aliases=alias_count)
    check = database_metadata(out_path, verify=True)
    if not check.get("available"):
        raise LitVar2Error(f"new LitVar2 SQLite failed validation: {check}")
    return check


def _remote_metadata(url: str) -> dict[str, object]:
    request = urllib.request.Request(
        url,
        method="HEAD",
        headers={"User-Agent": "NCKUH-NGS-UI-LitVar2/1.0"},
    )
    with urllib.request.urlopen(
        request,
        timeout=DOWNLOAD_TIMEOUT_SECONDS,
    ) as response:
        headers = response.headers
        return {
            "source_url": url,
            "source_last_modified": headers.get("Last-Modified", ""),
            "source_etag": headers.get("ETag", ""),
            "source_content_length": _to_int(headers.get("Content-Length")),
        }


def _dataset_date(last_modified: str) -> str:
    if last_modified:
        try:
            dt = email.utils.parsedate_to_datetime(last_modified)
            return dt.date().isoformat()
        except (TypeError, ValueError, OverflowError):
            pass
    return datetime.now(timezone.utc).date().isoformat()


def _download(
    url: str,
    destination: Path,
    *,
    expected_size: int = 0,
    progress: Progress | None = None,
) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(DOWNLOAD_RETRIES):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "NCKUH-NGS-UI-LitVar2/1.0"},
            )
            digest = hashlib.sha256()
            downloaded = 0
            with urllib.request.urlopen(
                request,
                timeout=DOWNLOAD_TIMEOUT_SECONDS,
            ) as response, destination.open("wb") as output:
                total = _to_int(response.headers.get("Content-Length")) or expected_size
                while True:
                    chunk = response.read(CHUNK_SIZE)
                    if not chunk:
                        break
                    output.write(chunk)
                    digest.update(chunk)
                    downloaded += len(chunk)
                    if downloaded == len(chunk) or downloaded % (64 * CHUNK_SIZE) < CHUNK_SIZE:
                        _notify(
                            progress,
                            "downloading",
                            downloaded_bytes=downloaded,
                            total_bytes=total,
                            percent=round(downloaded * 100 / total, 1) if total else None,
                        )
            if expected_size and downloaded != expected_size:
                raise LitVar2Error(
                    f"download size mismatch: expected {expected_size}, got {downloaded}"
                )
            return digest.hexdigest()
        except (OSError, urllib.error.URLError, LitVar2Error) as exc:
            last_error = exc
            destination.unlink(missing_ok=True)
            if attempt + 1 < DOWNLOAD_RETRIES:
                time.sleep(min(2 ** (attempt + 1), 16))
    raise LitVar2Error(f"LitVar2 bulk download failed: {last_error}")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(CHUNK_SIZE)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def update_database(
    *,
    directory: Path = LITVAR2_DIR,
    bulk_path: Path = LITVAR2_BULK_PATH,
    db_path: Path = LITVAR2_DB,
    manifest_path: Path = LITVAR2_MANIFEST_PATH,
    source_url: str = LITVAR2_BULK_URL,
    local_bulk: Path | None = None,
    dataset_date: str = "",
    force: bool = False,
    progress: Progress | None = None,
) -> dict[str, object]:
    """Download/check the official bulk export and atomically publish SQLite."""
    directory = Path(directory)
    bulk_path = Path(bulk_path)
    db_path = Path(db_path)
    manifest_path = Path(manifest_path)
    directory.mkdir(parents=True, exist_ok=True)

    remote: dict[str, object]
    candidate_bulk = bulk_path
    candidate_is_temporary = False
    if local_bulk is not None:
        local_bulk = Path(local_bulk)
        if not local_bulk.is_file():
            raise FileNotFoundError(local_bulk)
        remote = {
            "source_url": source_url or str(local_bulk),
            "source_last_modified": "",
            "source_etag": "",
            "source_content_length": local_bulk.stat().st_size,
        }
        source_sha256 = _sha256_file(local_bulk)
        if local_bulk.resolve() != bulk_path.resolve():
            candidate_bulk = bulk_path.with_suffix(bulk_path.suffix + ".building")
            candidate_bulk.unlink(missing_ok=True)
            shutil.copyfile(local_bulk, candidate_bulk)
            candidate_is_temporary = True
        else:
            candidate_bulk = local_bulk
    else:
        _notify(progress, "checking-source", source_url=source_url)
        remote = _remote_metadata(source_url)
        current = _read_json(manifest_path)
        unchanged = (
            not force
            and bulk_path.is_file()
            and db_path.is_file()
            and current.get("source_last_modified")
            and current.get("source_last_modified") == remote.get("source_last_modified")
            and (
                not remote.get("source_content_length")
                or current.get("source_content_length")
                == remote.get("source_content_length")
            )
            and database_metadata(db_path).get("available")
        )
        if unchanged:
            result = database_metadata(db_path)
            result.update(action="already-current")
            _notify(progress, "already-current", **result)
            return result

        download_tmp = bulk_path.with_suffix(bulk_path.suffix + ".download")
        source_sha256 = _download(
            source_url,
            download_tmp,
            expected_size=_to_int(remote.get("source_content_length")),
            progress=progress,
        )
        candidate_bulk = download_tmp
        candidate_is_temporary = True

    version_date = dataset_date or _dataset_date(
        str(remote.get("source_last_modified") or "")
    )
    db_tmp = db_path.with_suffix(db_path.suffix + ".building")
    try:
        result = build_database(
            candidate_bulk,
            db_tmp,
            dataset_date=version_date,
            source_url=str(remote.get("source_url") or source_url),
            source_last_modified=str(remote.get("source_last_modified") or ""),
            source_sha256=source_sha256,
            progress=progress,
        )
        _notify(progress, "promoting", dataset_date=version_date)
        if candidate_bulk.resolve() != bulk_path.resolve():
            os.replace(candidate_bulk, bulk_path)
            candidate_is_temporary = False
        os.replace(db_tmp, db_path)
    finally:
        db_tmp.unlink(missing_ok=True)
        if candidate_is_temporary:
            candidate_bulk.unlink(missing_ok=True)

    manifest = {
        **remote,
        "dataset_date": version_date,
        "source_sha256": source_sha256,
        "bulk_path": str(bulk_path),
        "database_path": str(db_path),
        "record_count": result.get("record_count", 0),
        "alias_count": result.get("alias_count", 0),
        "updated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    _atomic_json(manifest_path, manifest)
    return {**result, "action": "updated"}


def open_readonly(path: Path = LITVAR2_DB) -> sqlite3.Connection:
    path = Path(path)
    meta = database_metadata(path)
    if not meta.get("available"):
        raise LitVar2Error(f"LitVar2 database unavailable: {meta}")
    uri = f"file:{urllib.parse.quote(str(path.resolve()))}?mode=ro"
    conn = sqlite3.connect(uri, uri=True)
    conn.row_factory = sqlite3.Row
    return conn


def _rows_to_hits(rows: Iterable[sqlite3.Row]) -> dict[int, sqlite3.Row]:
    return {int(row["id"]): row for row in rows}


def _result_from_row(
    row: sqlite3.Row,
    *,
    dataset_date: str,
    match_method: str,
    query: str,
) -> dict[str, object]:
    try:
        pmids = [
            str(value) for value in json.loads(row["top_pmids_json"])
            if str(value).isdigit()
        ][:5]
    except (TypeError, json.JSONDecodeError):
        pmids = []
    litvar_id = str(row["litvar_id"])
    query_text = query or str(row["rsid"] or row["preferred_hgvs"] or "")
    params = urllib.parse.urlencode({
        "variant": litvar_id,
        "query": query_text,
    })
    return {
        "status": "hit",
        "litvar_id": litvar_id,
        "rsid": str(row["rsid"] or ""),
        "gene": str(row["primary_gene"] or ""),
        "hgvs": str(row["preferred_hgvs"] or ""),
        "pmids_count": _to_int(row["pmids_count"]),
        "pmids": pmids,
        "dataset_date": dataset_date,
        "match_method": match_method,
        "url": f"https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?{params}",
    }


def _ambiguous_result(
    hits: dict[int, sqlite3.Row],
    *,
    dataset_date: str,
    match_method: str,
) -> dict[str, object]:
    candidates = []
    for _row_id, row in sorted(hits.items()):
        query = (
            str(row["rsid"] or "")
            if match_method == "rsid"
            else " ".join(filter(None, (
                str(row["primary_gene"] or ""),
                str(row["preferred_hgvs"] or ""),
            )))
        )
        candidates.append(_result_from_row(
            row,
            dataset_date=dataset_date,
            match_method=match_method,
            query=query,
        ))
    return {
        "status": "ambiguous",
        "dataset_date": dataset_date,
        "match_method": match_method,
        "candidates": candidates,
    }


def lookup_variant(
    conn: sqlite3.Connection,
    *,
    rsids: Iterable[str] = (),
    genes: Iterable[str] = (),
    hgvs_values: Iterable[str] = (),
) -> dict[str, object]:
    """Resolve one genomic variant conservatively against the local index."""
    meta = _database_meta(conn)
    dataset_date = meta.get("dataset_date", "")
    normalized_rsids = list(dict.fromkeys(
        value for value in (normalize_rsid(raw) for raw in rsids) if value
    ))
    rsid_hits: dict[int, sqlite3.Row] = {}
    for rsid in normalized_rsids:
        rows = conn.execute(
            "SELECT * FROM variants WHERE rsid = ? ORDER BY id",
            (rsid,),
        ).fetchall()
        rsid_hits.update(_rows_to_hits(rows))
    if len(rsid_hits) == 1:
        return _result_from_row(
            next(iter(rsid_hits.values())),
            dataset_date=dataset_date,
            match_method="rsid",
            query=normalized_rsids[0],
        )
    if len(rsid_hits) > 1:
        return _ambiguous_result(
            rsid_hits,
            dataset_date=dataset_date,
            match_method="rsid",
        )

    normalized_genes = list(dict.fromkeys(
        value for value in (normalize_gene(raw) for raw in genes) if value
    ))
    normalized_hgvs = list(dict.fromkeys(
        value for value in (normalize_hgvs(raw) for raw in hgvs_values) if value
    ))
    hgvs_hits: dict[int, sqlite3.Row] = {}
    for gene in normalized_genes:
        for hgvs in normalized_hgvs:
            rows = conn.execute(
                "SELECT v.* FROM hgvs_aliases a "
                "JOIN variants v ON v.id = a.variant_id "
                "WHERE a.gene = ? AND a.hgvs = ? ORDER BY v.id",
                (gene, hgvs),
            ).fetchall()
            hgvs_hits.update(_rows_to_hits(rows))
    if len(hgvs_hits) == 1:
        query = " ".join(filter(None, (
            normalized_genes[0] if normalized_genes else "",
            normalized_hgvs[0] if normalized_hgvs else "",
        )))
        return _result_from_row(
            next(iter(hgvs_hits.values())),
            dataset_date=dataset_date,
            match_method="gene_hgvs",
            query=query,
        )
    if len(hgvs_hits) > 1:
        return _ambiguous_result(
            hgvs_hits,
            dataset_date=dataset_date,
            match_method="gene_hgvs",
        )
    return {
        "status": "no_match",
        "dataset_date": dataset_date,
        "match_method": "gene_hgvs" if normalized_hgvs else "",
    }
