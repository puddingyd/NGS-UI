from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_documents_buttons_and_shared_assets_are_loaded_on_both_pages():
    main = (ROOT / "frontend" / "index.html").read_text(encoding="utf-8")
    phenotype = (ROOT / "frontend" / "phenotype" / "index.html").read_text(encoding="utf-8")

    assert 'id="btn-patient-documents"' in main
    assert ">Documents</button>" in main
    assert 'src="./documents.js"' in main
    assert 'href="./documents.css"' in main
    assert 'id="btn-patient-documents"' in phenotype
    assert 'src="../documents.js"' in phenotype
    assert 'href="../documents.css"' in phenotype


def test_documents_frontend_supports_paste_rename_delete_and_tiff_preview():
    script = (ROOT / "frontend" / "documents.js").read_text(encoding="utf-8")

    assert 'addEventListener("paste"' in script
    assert "Screenshot_" in script
    assert ".tif,.tiff" in script
    assert 'data-pdoc-preview' in script
    assert 'method: "PATCH"' in script
    assert 'method: "DELETE"' in script
    assert 'method: "POST"' in script
    assert 'credentials: "same-origin"' in script


def test_documents_api_requires_authentication_and_streams_downloads():
    router = (ROOT / "backend" / "app" / "routers" / "documents.py").read_text(encoding="utf-8")
    service = (ROOT / "backend" / "app" / "services" / "patient_documents.py").read_text(encoding="utf-8")

    assert "dependencies=[Depends(current_user)]" in router
    assert "FileResponse(" in router
    assert "while True:" in service
    assert "await upload.read(_CHUNK_SIZE)" in service
    assert "PATIENT_DOCUMENTS_MIN_FREE_GB" in service
    assert "_PREVIEW_MAX_SIDE" in service
