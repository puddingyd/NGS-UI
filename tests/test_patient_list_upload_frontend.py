from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
INDEX_HTML = (REPO_ROOT / "frontend" / "index.html").read_text(encoding="utf-8")
STYLE_CSS = (REPO_ROOT / "frontend" / "style.css").read_text(encoding="utf-8")


def test_patient_list_picker_accepts_multiple_files():
    assert 'id="upload-list-file" type="file" accept=".xlsx,.xlsm" multiple' in INDEX_HTML
    assert "const files = Array.from(file.files || []);" in APP_JS
    assert "for (let index = 0; index < files.length; index += 1)" in APP_JS
    assert "files.slice(index + 1)" in APP_JS


def test_current_batch_results_are_separate_from_collapsed_history():
    assert '<section id="patient-list-current-upload" hidden>' in INDEX_HTML
    assert '<details id="patient-list-history-details"' in INDEX_HTML
    assert "_renderPatientListCurrentUpload(results);" in APP_JS
    assert 'historyDetails?.addEventListener("toggle"' in APP_JS
    assert "if (historyDetails.open) _renderPatientListHistory();" in APP_JS
    assert ".patient-list-history-details[open] > summary::before" in STYLE_CSS
