import asyncio
import io
import sqlite3
import sys
import zipfile
from pathlib import Path

import pytest
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from app.services import patient_documents as store


class MemoryUpload:
    def __init__(self, filename: str, content: bytes):
        self.filename = filename
        self._body = io.BytesIO(content)

    async def read(self, size: int = -1) -> bytes:
        return self._body.read(size)

    async def close(self) -> None:
        self._body.close()


@pytest.fixture()
def document_store(tmp_path, monkeypatch):
    documents = tmp_path / "patient_documents"
    phenotype = tmp_path / "patient_phenotype"
    monkeypatch.setattr(store, "PATIENT_DOCUMENTS_DIR", documents)
    monkeypatch.setattr(store, "PHENOTYPE_DIR", phenotype)
    monkeypatch.setattr(store, "PATIENT_DOCUMENTS_MIN_FREE_GB", 0)
    return documents, phenotype


def _image_bytes(fmt="PNG", *, pages=1):
    out = io.BytesIO()
    images = [Image.new("RGB", (40 + i, 30 + i), (120, 20 + i, 40)) for i in range(pages)]
    if pages > 1:
        images[0].save(out, format=fmt, save_all=True, append_images=images[1:])
    else:
        images[0].save(out, format=fmt)
    return out.getvalue()


def _upload(filename, content, *, mrn="12345678", display_name="", user=None):
    return asyncio.run(store.save_upload(
        MemoryUpload(filename, content),
        mrn=mrn,
        display_name=display_name,
        source_sample_id="26WE0001",
        user=user or {"id": 7, "username": "reviewer"},
    ))


def test_upload_rename_preview_download_and_soft_delete(document_store):
    documents, _ = document_store
    saved = _upload(
        "scan.png",
        _image_bytes(),
        display_name="門診截圖.png",
    )

    assert saved["display_name"] == "門診截圖.png"
    assert saved["source_sample_id"] == "26WE0001"
    assert saved["previewable"] is True
    assert saved["image_width"] == 40
    assert store.list_documents("12345678") == [saved]

    preview, preview_info = store.render_preview(saved["id"], page=0)
    assert preview.startswith(b"\x89PNG")
    assert preview_info["image_pages"] == 1

    renamed = store.rename_document(
        saved["id"],
        "追蹤影像.png",
        user={"id": 8, "username": "editor"},
    )
    assert renamed["display_name"] == "追蹤影像.png"
    path, download_info = store.document_file(saved["id"], action="download")
    assert path.is_file()
    assert path.name.startswith(saved["id"])
    assert download_info["display_name"] == "追蹤影像.png"

    result = store.delete_document(
        saved["id"],
        user={"id": 8, "username": "editor"},
    )
    assert result == {"id": saved["id"], "deleted": True}
    assert store.list_documents("12345678") == []
    assert (documents / "trash" / path.name).is_file()

    with sqlite3.connect(documents / "documents.sqlite") as conn:
        actions = [row[0] for row in conn.execute(
            "SELECT action FROM document_events WHERE document_id=? ORDER BY id",
            (saved["id"],),
        )]
    assert actions == ["upload", "preview", "rename", "download", "delete"]


def test_tiff_preview_supports_multiple_pages(document_store):
    _documents, _phenotype = document_store
    saved = _upload(
        "multipage.tiff",
        _image_bytes("TIFF", pages=2),
        display_name="掃描病歷.tiff",
    )

    assert saved["file_format"] == "TIFF"
    assert saved["image_pages"] == 2
    page_two, info = store.render_preview(saved["id"], page=1)
    assert page_two.startswith(b"\x89PNG")
    assert info["image_pages"] == 2
    with pytest.raises(store.InvalidDocument, match="頁碼"):
        store.render_preview(saved["id"], page=2)


def test_download_all_streams_zip_with_display_names(document_store):
    documents, _phenotype = document_store
    first = _upload(
        "scan.png",
        _image_bytes(),
        display_name="門診截圖.png",
    )
    second = _upload(
        "record.pdf",
        b"%PDF-1.4\narchive fixture\n%%EOF\n",
        display_name="病歷摘要.pdf",
    )

    chunks, count = store.stream_archive(
        "12345678",
        user={"id": 8, "username": "downloader"},
    )
    body = b"".join(chunks)

    assert count == 2
    with zipfile.ZipFile(io.BytesIO(body)) as archive:
        assert archive.namelist() == ["門診截圖.png", "病歷摘要.pdf"]
        assert archive.read("門診截圖.png").startswith(b"\x89PNG")
        assert archive.read("病歷摘要.pdf").startswith(b"%PDF-")

    with sqlite3.connect(documents / "documents.sqlite") as conn:
        archived_ids = [row[0] for row in conn.execute(
            """SELECT document_id FROM document_events
               WHERE action='archive_download' ORDER BY id"""
        )]
    assert archived_ids == [first["id"], second["id"]]


def test_pdf_is_downloadable_but_not_previewable(document_store):
    _documents, _phenotype = document_store
    saved = _upload(
        "record.pdf",
        b"%PDF-1.4\nminimal test fixture\n%%EOF\n",
        display_name="病歷摘要.pdf",
    )

    assert saved["previewable"] is False
    with pytest.raises(store.PreviewUnavailable):
        store.render_preview(saved["id"])


def test_duplicate_content_and_mismatched_extension_are_rejected(document_store):
    _documents, _phenotype = document_store
    content = _image_bytes()
    _upload("first.png", content, display_name="第一張.png")
    with pytest.raises(store.DocumentConflict, match="相同檔案"):
        _upload("second.png", content, display_name="第二張.png")
    with pytest.raises(store.InvalidDocument, match="副檔名"):
        _upload("third.png", _image_bytes(), display_name="錯誤.pdf", mrn="87654321")


def test_mrn_move_moves_documents_and_patient_sidecars(document_store):
    documents, phenotype = document_store
    saved = _upload("scan.png", _image_bytes(), display_name="影像.png", mrn="OLD123")
    phenotype.mkdir(parents=True)
    (phenotype / "OLD123_clinical_presentation.txt").write_text("clinical\n", encoding="utf-8")
    (phenotype / "26WE0001_OLD123_phenotype.txt").write_text("phenotype\n", encoding="utf-8")

    result = store.move_mrn(
        "OLD123",
        "NEW456",
        user={"id": 7, "username": "reviewer"},
    )

    assert result["documents"] == 1
    assert sorted(result["sidecars"]) == [
        "26WE0001_NEW456_phenotype.txt",
        "NEW456_clinical_presentation.txt",
    ]
    assert store.list_documents("OLD123") == []
    assert store.list_documents("NEW456")[0]["id"] == saved["id"]
    assert not (documents / "files" / "OLD123").exists()
    assert (documents / "files" / "NEW456").is_dir()
    assert (phenotype / "NEW456_clinical_presentation.txt").read_text(encoding="utf-8") == "clinical\n"
    assert (phenotype / "26WE0001_NEW456_phenotype.txt").is_file()


def test_mrn_move_refuses_to_merge_existing_patients(document_store):
    _documents, phenotype = document_store
    _upload("old.png", _image_bytes(), display_name="舊病人.png", mrn="OLD123")
    _upload("new.png", _image_bytes("PNG", pages=1) + b"different", display_name="新病人.png", mrn="NEW456")
    phenotype.mkdir(parents=True, exist_ok=True)

    with pytest.raises(store.DocumentConflict, match="不能自動合併"):
        store.move_mrn("OLD123", "NEW456")


def test_mrn_move_refuses_to_merge_sidecar_with_target_documents(document_store):
    _documents, phenotype = document_store
    phenotype.mkdir(parents=True)
    (phenotype / "OLD123_clinical_presentation.txt").write_text("old\n", encoding="utf-8")
    _upload("new.png", _image_bytes(), display_name="新病人.png", mrn="NEW456")

    with pytest.raises(store.DocumentConflict, match="不能自動合併"):
        store.move_mrn("OLD123", "NEW456")
    assert (phenotype / "OLD123_clinical_presentation.txt").is_file()


def test_invalid_and_unsupported_uploads_are_rejected(document_store):
    _documents, _phenotype = document_store
    with pytest.raises(store.InvalidDocument, match="病歷號"):
        _upload("scan.png", _image_bytes(), mrn="bad/mrn")
    with pytest.raises(store.InvalidDocument, match="有效"):
        _upload("notes.txt", b"not an image", mrn="12345678")
