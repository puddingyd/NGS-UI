from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
STYLE_CSS = (REPO_ROOT / "frontend" / "style.css").read_text(encoding="utf-8")
INDEX_HTML = (REPO_ROOT / "frontend" / "index.html").read_text(encoding="utf-8")


def test_ploidy_modal_is_reviewer_first_and_keeps_raw_values_collapsed():
    assert "Aneuploidy signal detected" in APP_JS
    assert "No aneuploidy signal" in APP_JS
    assert "需要複核的染色體" in INDEX_HTML
    assert "顯示全部染色體" in INDEX_HTML
    assert "顯示原始測量值" in INDEX_HTML
    assert "<th>ALT</th>" in INDEX_HTML
    assert "<th>QUAL</th>" in INDEX_HTML
    assert "Raw ratio" in INDEX_HTML


def test_ploidy_modal_explains_source_specific_confidence_and_derived_ratio():
    assert 'confidence === "low"' in APP_JS
    assert 'confidence === "suspect"' in APP_JS
    assert "Derived ratio = DC / autosomeDepthOfCoverage" in APP_JS
    assert "非 DRAGEN 原始欄位" in APP_JS
    assert ".ploidy-quality-low" in STYLE_CSS
    assert ".ploidy-quality-suspect" in STYLE_CSS

