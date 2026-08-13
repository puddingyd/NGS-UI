"""MRN-level patient document storage.

Original uploads live outside pipeline output and phenotype sidecars.  The
SQLite catalog keeps reviewer-visible names and audit events while physical
files use UUID names, so renaming never rewrites a large medical record.
"""
from __future__ import annotations

import hashlib
import io
import json
import os
import re
import shutil
import sqlite3
import threading
import unicodedata
import uuid
import warnings
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator, Protocol

from PIL import Image, ImageOps, UnidentifiedImageError

from ..config import (
    PATIENT_DOCUMENTS_DIR,
    PATIENT_DOCUMENTS_MIN_FREE_GB,
    PHENOTYPE_DIR,
)


_TOKEN_RE = re.compile(r"^[A-Za-z0-9_-]{1,32}$")
_CONTROL_RE = re.compile(r"[\x00-\x1f\x7f]+")
_WRITE_LOCK = threading.RLock()
_CHUNK_SIZE = 1024 * 1024
_PREVIEW_MAX_SIDE = 2400

_FORMATS = {
    "PDF": {
        "extension": ".pdf",
        "extensions": {".pdf"},
        "content_type": "application/pdf",
        "image": False,
    },
    "JPEG": {
        "extension": ".jpg",
        "extensions": {".jpg", ".jpeg"},
        "content_type": "image/jpeg",
        "image": True,
    },
    "PNG": {
        "extension": ".png",
        "extensions": {".png"},
        "content_type": "image/png",
        "image": True,
    },
    "TIFF": {
        "extension": ".tiff",
        "extensions": {".tif", ".tiff"},
        "content_type": "image/tiff",
        "image": True,
    },
}


class UploadLike(Protocol):
    filename: str | None

    async def read(self, size: int = -1) -> bytes: ...

    async def close(self) -> None: ...


class DocumentError(RuntimeError):
    pass


class InvalidDocument(DocumentError):
    pass


class DocumentConflict(DocumentError):
    pass


class DocumentNotFound(DocumentError):
    pass


class DocumentStorageFull(DocumentError):
    pass


class PreviewUnavailable(DocumentError):
    pass


class _ZipStreamBuffer:
    """Minimal unseekable sink used by zipfile for incremental output."""

    def __init__(self) -> None:
        self._buffer = bytearray()
        self._position = 0

    def write(self, data: bytes) -> int:
        self._buffer.extend(data)
        self._position += len(data)
        return len(data)

    def tell(self) -> int:
        return self._position

    def flush(self) -> None:
        return None

    def seekable(self) -> bool:
        return False

    def drain(self) -> bytes:
        data = bytes(self._buffer)
        self._buffer.clear()
        return data


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def validate_mrn(value: str) -> str:
    mrn = str(value or "").strip()
    if not mrn:
        raise InvalidDocument("請先填寫病歷號")
    if not _TOKEN_RE.fullmatch(mrn):
        raise InvalidDocument("病歷號只能是英數 / - / _（最多 32 字）")
    return mrn


def _root() -> Path:
    return Path(PATIENT_DOCUMENTS_DIR)


def _db_path() -> Path:
    return _root() / "documents.sqlite"


def _files_root() -> Path:
    return _root() / "files"


def _trash_root() -> Path:
    return _root() / "trash"


def _incoming_root() -> Path:
    return _root() / ".incoming"


def _mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    try:
        path.chmod(0o700)
    except OSError:
        pass


def _conn() -> sqlite3.Connection:
    _mkdir(_root())
    conn = sqlite3.connect(_db_path(), timeout=30)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=WAL")
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS documents (
            id                  TEXT PRIMARY KEY,
            mrn                 TEXT NOT NULL,
            source_sample_id    TEXT NOT NULL DEFAULT '',
            display_name        TEXT NOT NULL,
            original_name       TEXT NOT NULL,
            stored_name         TEXT NOT NULL UNIQUE,
            content_type        TEXT NOT NULL,
            file_format         TEXT NOT NULL,
            size_bytes          INTEGER NOT NULL,
            sha256              TEXT NOT NULL,
            image_width         INTEGER,
            image_height        INTEGER,
            image_pages         INTEGER,
            created_at          TEXT NOT NULL,
            created_by_user_id  INTEGER,
            created_by_username TEXT NOT NULL DEFAULT '',
            updated_at          TEXT NOT NULL,
            updated_by_user_id  INTEGER,
            updated_by_username TEXT NOT NULL DEFAULT '',
            deleted_at          TEXT,
            deleted_by_user_id  INTEGER,
            deleted_by_username TEXT NOT NULL DEFAULT ''
        );
        CREATE UNIQUE INDEX IF NOT EXISTS uq_documents_active_name
            ON documents(mrn, display_name COLLATE NOCASE)
            WHERE deleted_at IS NULL;
        CREATE UNIQUE INDEX IF NOT EXISTS uq_documents_active_content
            ON documents(mrn, sha256)
            WHERE deleted_at IS NULL;
        CREATE INDEX IF NOT EXISTS ix_documents_mrn_created
            ON documents(mrn, created_at DESC);

        CREATE TABLE IF NOT EXISTS document_events (
            id          INTEGER PRIMARY KEY AUTOINCREMENT,
            document_id TEXT NOT NULL,
            action      TEXT NOT NULL,
            username    TEXT NOT NULL DEFAULT '',
            user_id     INTEGER,
            occurred_at TEXT NOT NULL,
            details     TEXT NOT NULL DEFAULT '{}'
        );
        CREATE INDEX IF NOT EXISTS ix_document_events_document
            ON document_events(document_id, occurred_at DESC);
        """
    )
    return conn


def _event(
    conn: sqlite3.Connection,
    document_id: str,
    action: str,
    user: dict | None,
    details: dict | None = None,
) -> None:
    conn.execute(
        """INSERT INTO document_events
           (document_id, action, username, user_id, occurred_at, details)
           VALUES (?, ?, ?, ?, ?, ?)""",
        (
            document_id,
            action,
            str((user or {}).get("username") or ""),
            (user or {}).get("id"),
            _now(),
            json.dumps(details or {}, ensure_ascii=False, separators=(",", ":")),
        ),
    )


def _public(row: sqlite3.Row | dict) -> dict:
    data = dict(row)
    return {
        "id": data["id"],
        "mrn": data["mrn"],
        "source_sample_id": data.get("source_sample_id") or "",
        "display_name": data["display_name"],
        "original_name": data["original_name"],
        "content_type": data["content_type"],
        "file_format": data["file_format"],
        "size_bytes": int(data["size_bytes"] or 0),
        "sha256": data["sha256"],
        "previewable": bool(_FORMATS.get(data["file_format"], {}).get("image")),
        "image_width": data.get("image_width"),
        "image_height": data.get("image_height"),
        "image_pages": int(data.get("image_pages") or 1),
        "created_at": data["created_at"],
        "created_by_username": data.get("created_by_username") or "",
        "updated_at": data["updated_at"],
        "updated_by_username": data.get("updated_by_username") or "",
    }


def _active_row(conn: sqlite3.Connection, document_id: str) -> sqlite3.Row:
    row = conn.execute(
        "SELECT * FROM documents WHERE id=? AND deleted_at IS NULL",
        (str(document_id or ""),),
    ).fetchone()
    if row is None:
        raise DocumentNotFound("找不到文件")
    return row


def _stored_path(row: sqlite3.Row | dict) -> Path:
    data = dict(row)
    mrn = validate_mrn(data["mrn"])
    stored_name = str(data["stored_name"] or "")
    if Path(stored_name).name != stored_name:
        raise DocumentNotFound("文件路徑不合法")
    path = _files_root() / mrn / stored_name
    expected_parent = (_files_root() / mrn).resolve()
    if path.resolve().parent != expected_parent:
        raise DocumentNotFound("文件路徑不合法")
    return path


def _detect(path: Path) -> dict:
    with path.open("rb") as fh:
        prefix = fh.read(8)
    if prefix.startswith(b"%PDF-"):
        return {
            "file_format": "PDF",
            "content_type": _FORMATS["PDF"]["content_type"],
            "extension": ".pdf",
            "image_width": None,
            "image_height": None,
            "image_pages": None,
        }
    magic_format = ""
    if prefix.startswith(b"\x89PNG\r\n\x1a\n"):
        magic_format = "PNG"
    elif prefix.startswith(b"\xff\xd8\xff"):
        magic_format = "JPEG"
    elif prefix[:4] in {b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+"}:
        magic_format = "TIFF"
    try:
        with Image.open(path) as image:
            file_format = str(image.format or "").upper()
            if file_format not in {"JPEG", "PNG", "TIFF"}:
                raise InvalidDocument("只支援 PDF、JPG、PNG、TIF、TIFF")
            width, height = image.size
            pages = int(getattr(image, "n_frames", 1) or 1)
            image.verify()
    except InvalidDocument:
        raise
    except Image.DecompressionBombError:
        # Preserve the original even when its declared pixel dimensions are
        # too large for safe decoding. The preview endpoint will reject it,
        # while download remains available and upload still has no byte cap.
        if not magic_format:
            raise InvalidDocument("檔案內容不是有效的 PDF、JPG、PNG 或 TIFF")
        file_format = magic_format
        width = height = None
        pages = 1
    except (UnidentifiedImageError, OSError, SyntaxError, ValueError) as exc:
        raise InvalidDocument("檔案內容不是有效的 PDF、JPG、PNG 或 TIFF") from exc
    spec = _FORMATS[file_format]
    return {
        "file_format": file_format,
        "content_type": spec["content_type"],
        "extension": spec["extension"],
        "image_width": int(width) if width is not None else None,
        "image_height": int(height) if height is not None else None,
        "image_pages": pages,
    }


def _clean_leaf_name(value: str) -> str:
    value = unicodedata.normalize("NFC", str(value or ""))
    value = value.replace("\\", "/").split("/")[-1]
    value = _CONTROL_RE.sub(" ", value).strip().strip(".")
    value = re.sub(r"\s+", " ", value)
    return value


def normalize_display_name(
    requested: str,
    *,
    original_name: str,
    file_format: str,
) -> str:
    spec = _FORMATS.get(file_format)
    if not spec:
        raise InvalidDocument("不支援的文件格式")
    name = _clean_leaf_name(requested) or _clean_leaf_name(original_name)
    if not name:
        name = "Document"
    suffix = Path(name).suffix.lower()
    if suffix:
        if suffix not in spec["extensions"]:
            allowed = " / ".join(sorted(spec["extensions"]))
            raise InvalidDocument(f"副檔名必須是 {allowed}")
    else:
        original_suffix = Path(_clean_leaf_name(original_name)).suffix.lower()
        name += original_suffix if original_suffix in spec["extensions"] else spec["extension"]
    if len(name) > 200:
        stem = Path(name).stem[: 199 - len(Path(name).suffix)].rstrip()
        name = stem + Path(name).suffix
    if not Path(name).stem.strip():
        raise InvalidDocument("檔名不可為空")
    return name


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _check_disk_reserve(next_chunk_bytes: int = 0) -> None:
    _mkdir(_root())
    reserve = int(PATIENT_DOCUMENTS_MIN_FREE_GB) * 1024**3
    if reserve <= 0:
        return
    free = shutil.disk_usage(_root()).free
    if free - int(next_chunk_bytes or 0) < reserve:
        raise DocumentStorageFull(
            f"文件儲存空間不足；系統需保留至少 {PATIENT_DOCUMENTS_MIN_FREE_GB} GB"
        )


async def save_upload(
    upload: UploadLike,
    *,
    mrn: str,
    display_name: str = "",
    source_sample_id: str = "",
    user: dict | None = None,
) -> dict:
    mrn = validate_mrn(mrn)
    _mkdir(_incoming_root())
    incoming = _incoming_root() / f"{uuid.uuid4().hex}.upload"
    size = 0
    try:
        with incoming.open("xb") as out:
            try:
                incoming.chmod(0o600)
            except OSError:
                pass
            while True:
                chunk = await upload.read(_CHUNK_SIZE)
                if not chunk:
                    break
                _check_disk_reserve(len(chunk))
                out.write(chunk)
                size += len(chunk)
        if size <= 0:
            raise InvalidDocument("空白檔案")
        detected = _detect(incoming)
        original_name = _clean_leaf_name(upload.filename or "") or (
            "Document" + detected["extension"]
        )
        final_display_name = normalize_display_name(
            display_name,
            original_name=original_name,
            file_format=detected["file_format"],
        )
        digest = _sha256(incoming)
        document_id = str(uuid.uuid4())
        stored_name = document_id + detected["extension"]
        patient_dir = _files_root() / mrn
        _mkdir(patient_dir)
        final_path = patient_dir / stored_name
        now = _now()
        row = {
            "id": document_id,
            "mrn": mrn,
            "source_sample_id": str(source_sample_id or "").strip(),
            "display_name": final_display_name,
            "original_name": original_name,
            "stored_name": stored_name,
            "content_type": detected["content_type"],
            "file_format": detected["file_format"],
            "size_bytes": size,
            "sha256": digest,
            "image_width": detected["image_width"],
            "image_height": detected["image_height"],
            "image_pages": detected["image_pages"],
            "created_at": now,
            "created_by_user_id": (user or {}).get("id"),
            "created_by_username": str((user or {}).get("username") or ""),
            "updated_at": now,
            "updated_by_user_id": (user or {}).get("id"),
            "updated_by_username": str((user or {}).get("username") or ""),
        }
        with _WRITE_LOCK:
            with _conn() as conn:
                duplicate = conn.execute(
                    "SELECT display_name FROM documents WHERE mrn=? AND sha256=? AND deleted_at IS NULL",
                    (mrn, digest),
                ).fetchone()
                if duplicate:
                    raise DocumentConflict(
                        f"相同檔案已存在：{duplicate['display_name']}"
                    )
                os.replace(incoming, final_path)
                try:
                    final_path.chmod(0o600)
                except OSError:
                    pass
                try:
                    conn.execute(
                        """INSERT INTO documents (
                            id, mrn, source_sample_id, display_name, original_name,
                            stored_name, content_type, file_format, size_bytes,
                            sha256, image_width, image_height, image_pages,
                            created_at, created_by_user_id, created_by_username,
                            updated_at, updated_by_user_id, updated_by_username
                        ) VALUES (
                            :id, :mrn, :source_sample_id, :display_name, :original_name,
                            :stored_name, :content_type, :file_format, :size_bytes,
                            :sha256, :image_width, :image_height, :image_pages,
                            :created_at, :created_by_user_id, :created_by_username,
                            :updated_at, :updated_by_user_id, :updated_by_username
                        )""",
                        row,
                    )
                    _event(conn, document_id, "upload", user, {
                        "display_name": final_display_name,
                        "size_bytes": size,
                        "sha256": digest,
                    })
                except Exception as exc:
                    if final_path.exists():
                        os.replace(final_path, incoming)
                    if isinstance(exc, sqlite3.IntegrityError):
                        raise DocumentConflict("同名或相同內容的文件已存在") from exc
                    raise
        return _public(row)
    finally:
        try:
            incoming.unlink(missing_ok=True)
        except OSError:
            pass
        try:
            await upload.close()
        except Exception:
            pass


def list_documents(mrn: str) -> list[dict]:
    mrn = validate_mrn(mrn)
    with _conn() as conn:
        rows = conn.execute(
            """SELECT * FROM documents
               WHERE mrn=? AND deleted_at IS NULL
               ORDER BY created_at DESC, id DESC""",
            (mrn,),
        ).fetchall()
    return [_public(row) for row in rows]


def stream_archive(
    mrn: str,
    *,
    user: dict | None = None,
) -> tuple[Iterator[bytes], int]:
    """Return a streaming ZIP of every active document for one MRN.

    The ZIP is never staged as a second file on disk. Images and PDFs are
    already compressed, so ZIP_STORED avoids wasting CPU while still giving
    the reviewer one portable archive with the current display names.
    """
    mrn = validate_mrn(mrn)
    with _conn() as conn:
        rows = conn.execute(
            """SELECT * FROM documents
               WHERE mrn=? AND deleted_at IS NULL
               ORDER BY created_at ASC, rowid ASC""",
            (mrn,),
        ).fetchall()
        if not rows:
            raise DocumentNotFound("目前沒有文件可下載")
        members: list[tuple[Path, str]] = []
        for row in rows:
            path = _stored_path(row)
            if not path.is_file():
                raise DocumentNotFound(f"文件內容不存在：{row['display_name']}")
            members.append((path, row["display_name"]))
            _event(conn, row["id"], "archive_download", user, {
                "display_name": row["display_name"],
                "mrn": mrn,
            })

    def _generate() -> Iterator[bytes]:
        sink = _ZipStreamBuffer()
        with zipfile.ZipFile(
            sink,
            mode="w",
            compression=zipfile.ZIP_STORED,
            allowZip64=True,
            strict_timestamps=False,
        ) as archive:
            for path, display_name in members:
                with path.open("rb") as source, archive.open(
                    display_name,
                    mode="w",
                    force_zip64=True,
                ) as target:
                    while True:
                        chunk = source.read(_CHUNK_SIZE)
                        if not chunk:
                            break
                        target.write(chunk)
                        output = sink.drain()
                        if output:
                            yield output
                output = sink.drain()
                if output:
                    yield output
        output = sink.drain()
        if output:
            yield output

    return _generate(), len(members)


def rename_document(
    document_id: str,
    display_name: str,
    *,
    user: dict | None = None,
) -> dict:
    with _WRITE_LOCK:
        with _conn() as conn:
            row = _active_row(conn, document_id)
            new_name = normalize_display_name(
                display_name,
                original_name=row["display_name"],
                file_format=row["file_format"],
            )
            old_name = row["display_name"]
            if new_name != old_name:
                now = _now()
                try:
                    conn.execute(
                        """UPDATE documents
                           SET display_name=?, updated_at=?, updated_by_user_id=?,
                               updated_by_username=?
                           WHERE id=?""",
                        (
                            new_name,
                            now,
                            (user or {}).get("id"),
                            str((user or {}).get("username") or ""),
                            row["id"],
                        ),
                    )
                except sqlite3.IntegrityError as exc:
                    raise DocumentConflict("同名文件已存在") from exc
                _event(conn, row["id"], "rename", user, {
                    "old_name": old_name,
                    "new_name": new_name,
                })
            updated = _active_row(conn, document_id)
            return _public(updated)


def delete_document(document_id: str, *, user: dict | None = None) -> dict:
    with _WRITE_LOCK:
        with _conn() as conn:
            row = _active_row(conn, document_id)
            source = _stored_path(row)
            if not source.is_file():
                raise DocumentNotFound("文件內容不存在")
            _mkdir(_trash_root())
            trash = _trash_root() / row["stored_name"]
            if trash.exists():
                raise DocumentConflict("文件回收區已有同名內容")
            os.replace(source, trash)
            try:
                now = _now()
                conn.execute(
                    """UPDATE documents
                       SET deleted_at=?, deleted_by_user_id=?, deleted_by_username=?,
                           updated_at=?, updated_by_user_id=?, updated_by_username=?
                       WHERE id=?""",
                    (
                        now,
                        (user or {}).get("id"),
                        str((user or {}).get("username") or ""),
                        now,
                        (user or {}).get("id"),
                        str((user or {}).get("username") or ""),
                        row["id"],
                    ),
                )
                _event(conn, row["id"], "delete", user, {
                    "display_name": row["display_name"],
                })
            except Exception:
                if trash.exists() and not source.exists():
                    os.replace(trash, source)
                raise
            try:
                source.parent.rmdir()
            except OSError:
                pass
            return {"id": row["id"], "deleted": True}


def document_file(
    document_id: str,
    *,
    user: dict | None = None,
    action: str = "download",
) -> tuple[Path, dict]:
    with _conn() as conn:
        row = _active_row(conn, document_id)
        path = _stored_path(row)
        if not path.is_file():
            raise DocumentNotFound("文件內容不存在")
        _event(conn, row["id"], action, user, {
            "display_name": row["display_name"],
        })
        return path, _public(row)


def render_preview(
    document_id: str,
    *,
    page: int = 0,
    user: dict | None = None,
) -> tuple[bytes, dict]:
    if page < 0:
        raise InvalidDocument("預覽頁碼不合法")
    path, info = document_file(document_id, user=user, action="preview")
    if not info["previewable"]:
        raise PreviewUnavailable("這個檔案格式不支援圖片預覽")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", Image.DecompressionBombWarning)
            with Image.open(path) as image:
                pages = int(getattr(image, "n_frames", 1) or 1)
                if page >= pages:
                    raise InvalidDocument("預覽頁碼超出範圍")
                image.seek(page)
                frame = ImageOps.exif_transpose(image.copy())
                frame.thumbnail(
                    (_PREVIEW_MAX_SIDE, _PREVIEW_MAX_SIDE),
                    Image.Resampling.LANCZOS,
                )
                if frame.mode in {"RGBA", "LA"}:
                    frame = frame.convert("RGBA")
                else:
                    frame = frame.convert("RGB")
                out = io.BytesIO()
                frame.save(out, format="PNG", optimize=True)
                info["image_pages"] = pages
                return out.getvalue(), info
    except InvalidDocument:
        raise
    except (Image.DecompressionBombWarning, Image.DecompressionBombError) as exc:
        raise PreviewUnavailable("圖片尺寸過大，請下載原檔查看") from exc
    except (UnidentifiedImageError, OSError, SyntaxError, ValueError) as exc:
        raise PreviewUnavailable("無法產生圖片預覽，請下載原檔查看") from exc


def _sidecar_moves(old_mrn: str, new_mrn: str) -> list[tuple[Path, Path]]:
    root = Path(PHENOTYPE_DIR)
    if not root.is_dir():
        return []
    candidates = [
        root / f"{old_mrn}_clinical_presentation.txt",
        root / f"{old_mrn}_phenotype.txt",
        *root.glob(f"*_{old_mrn}_clinical_presentation.txt"),
        *root.glob(f"*_{old_mrn}_phenotype.txt"),
    ]
    moves: list[tuple[Path, Path]] = []
    seen: set[Path] = set()
    for source in candidates:
        if source in seen or not source.is_file():
            continue
        seen.add(source)
        name = source.name
        if name == f"{old_mrn}_clinical_presentation.txt":
            target_name = f"{new_mrn}_clinical_presentation.txt"
        elif name == f"{old_mrn}_phenotype.txt":
            target_name = f"{new_mrn}_phenotype.txt"
        else:
            clinical_suffix = f"_{old_mrn}_clinical_presentation.txt"
            phenotype_suffix = f"_{old_mrn}_phenotype.txt"
            if name.endswith(clinical_suffix):
                target_name = name[:-len(clinical_suffix)] + f"_{new_mrn}_clinical_presentation.txt"
            elif name.endswith(phenotype_suffix):
                target_name = name[:-len(phenotype_suffix)] + f"_{new_mrn}_phenotype.txt"
            else:
                continue
        moves.append((source, root / target_name))
    return moves


def has_patient_data(mrn: str) -> bool:
    mrn = validate_mrn(mrn)
    with _conn() as conn:
        count = conn.execute(
            "SELECT COUNT(*) FROM documents WHERE mrn=? AND deleted_at IS NULL",
            (mrn,),
        ).fetchone()[0]
    return bool(count or _sidecar_moves(mrn, "__probe__"))


def move_mrn(
    old_mrn: str,
    new_mrn: str,
    *,
    user: dict | None = None,
) -> dict:
    old_mrn = validate_mrn(old_mrn)
    new_mrn = validate_mrn(new_mrn)
    if old_mrn == new_mrn:
        return {"old_mrn": old_mrn, "new_mrn": new_mrn, "documents": 0, "sidecars": []}

    sidecar_moves = _sidecar_moves(old_mrn, new_mrn)
    target_sidecars = _sidecar_moves(new_mrn, "__probe__")
    source_dir = _files_root() / old_mrn
    target_dir = _files_root() / new_mrn
    moved_paths: list[tuple[Path, Path]] = []
    with _WRITE_LOCK:
        with _conn() as conn:
            source_rows = conn.execute(
                "SELECT * FROM documents WHERE mrn=?",
                (old_mrn,),
            ).fetchall()
            target_active = conn.execute(
                "SELECT COUNT(*) FROM documents WHERE mrn=? AND deleted_at IS NULL",
                (new_mrn,),
            ).fetchone()[0]
            source_active = any(row["deleted_at"] is None for row in source_rows)
            source_has_data = bool(source_active or sidecar_moves or source_dir.exists())
            target_has_data = bool(target_active or target_sidecars or target_dir.exists())
            if source_has_data and target_has_data:
                raise DocumentConflict("新病歷號已經有文件，不能自動合併病人資料")
            for source, target in sidecar_moves:
                if target.exists():
                    raise DocumentConflict(
                        f"新病歷號已經有 {target.name}，不能自動合併病人資料"
                    )
            if source_dir.exists() and target_dir.exists():
                raise DocumentConflict("新病歷號已經有文件目錄，不能自動合併病人資料")

            try:
                if source_dir.exists():
                    _mkdir(target_dir.parent)
                    os.replace(source_dir, target_dir)
                    moved_paths.append((source_dir, target_dir))
                for source, target in sidecar_moves:
                    os.replace(source, target)
                    moved_paths.append((source, target))
                conn.execute(
                    "UPDATE documents SET mrn=?, updated_at=? WHERE mrn=?",
                    (new_mrn, _now(), old_mrn),
                )
                for row in source_rows:
                    _event(conn, row["id"], "mrn_move", user, {
                        "old_mrn": old_mrn,
                        "new_mrn": new_mrn,
                    })
            except Exception:
                for source, target in reversed(moved_paths):
                    try:
                        if target.exists() and not source.exists():
                            os.replace(target, source)
                    except OSError:
                        pass
                raise
    return {
        "old_mrn": old_mrn,
        "new_mrn": new_mrn,
        "documents": len(source_rows),
        "sidecars": [target.name for _, target in sidecar_moves],
    }
