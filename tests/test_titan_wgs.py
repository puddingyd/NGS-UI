from pathlib import Path

from docx import Document

from app.services import docx_export, patient_list_store, sample_loader, test_types


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
INDEX_HTML = (REPO_ROOT / "frontend" / "index.html").read_text(encoding="utf-8")
STYLE_CSS = (REPO_ROOT / "frontend" / "style.css").read_text(encoding="utf-8")


def test_year_t_sample_ids_are_titan_wgs():
    for sample_id in ("25T00001", "26T00039-dragen", "27T99999"):
        assert test_types.is_titan_sample_id(sample_id)
        assert test_types.normalize_test_type("WES", sample_id=sample_id) == "TITAN-WGS"
        assert sample_loader._effective_test_type(
            {"lis_id": sample_id, "test_type": "WGS"}, sample_id,
        ) == "TITAN-WGS"

    assert test_types.normalize_test_type("WGS", sample_id="26WE0001") == "WGS"
    assert patient_list_store._test_type_from_name("WES", "26T00039") == "TITAN-WGS"


def test_titan_wgs_uses_wgs_report_wording():
    assert test_types.is_wgs_type("TITAN-WGS")
    doc = Document()
    docx_export._section_test_info(doc, "TITAN-WGS", health=True)
    docx_export._section_methods(doc, "TITAN-WGS", health=True)
    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    assert "全基因體定序檢測" in text
    assert "Illumina NovaSeq X Plus" in text
    assert "平均定序深度 ≧ 27X" in text


def test_frontend_has_three_test_filters_and_only_buttons():
    assert 'const SAMPLE_TEST_TYPES = ["WES", "WGS", "TITAN-WGS"];' in APP_JS
    assert 'return "TITAN-WGS";' in APP_JS
    assert 'sampleTestFilters.clear();' in APP_JS
    assert '_caseListTestFilters.clear();' in APP_JS
    assert INDEX_HTML.count('value="TITAN-WGS"') == 4
    for test_type in ("WES", "WGS", "TITAN-WGS"):
        assert INDEX_HTML.count(f'data-test-type="{test_type}"') == 2


def test_titan_wgs_has_diagnostic_analysis_visibility_toggle():
    assert 'id="diagnostic-analysis-toggle-row"' in INDEX_HTML
    assert 'id="btn-toggle-diagnostic-analysis"' in INDEX_HTML
    assert INDEX_HTML.count("diagnostic-analysis-sidebar-link") == 5
    assert "function applyDiagnosticAnalysisVisibility()" in APP_JS
    assert "function setupDiagnosticAnalysisToggle()" in APP_JS
    assert '"titan-diagnostic-analysis-hidden"' in APP_JS
    assert 'btn.textContent = visible ? "▾ 隱藏診斷分析" : "▸ 顯示診斷分析";' in APP_JS

    for selector in (
        "#counseling-card",
        "#clinical-card",
        "#phenotype-card",
        "#dead-zone-card",
        ".btn-export-clinical",
        ".btn-print-report",
        "#sec-causative",
        "#sec-other",
        "#sec-candidate",
        "#card-snv",
        "#card-cnv-sv",
        "#card-mito",
        "#card-str",
        "#card-roh",
    ):
        assert f"body.titan-diagnostic-analysis-hidden {selector}" in STYLE_CSS

    # Health-screening content must not be hidden by the diagnostic-only rule.
    for selector in ("#sec-acmg-sf", "#sec-stroke", "#sec-carrier", "#sec-pharmcat"):
        assert f"body.titan-diagnostic-analysis-hidden {selector}" not in STYLE_CSS


def test_titan_wgs_opens_all_secondary_findings_by_default():
    assert "function applyTitanSecondaryFindingsDefault()" in APP_JS
    assert 'currentSampleTestType() === "TITAN-WGS"' in APP_JS
    assert '(isPanel && currentSampleTestType() === "TITAN-WGS")' in APP_JS
    assert 'currentSampleTestType() === "TITAN-WGS" || hostId === "cat-pharmcat-c"' in APP_JS


def test_analysis_heading_matches_report_heading_size():
    assert ".section-header h2 { margin: 0; font-size: 1.5em; }" in STYLE_CSS
