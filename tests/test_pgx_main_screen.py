from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
STYLE_CSS = (REPO_ROOT / "frontend" / "style.css").read_text(encoding="utf-8")


def test_pgx_main_screen_uses_shared_report_projection():
    assert "function renderPharmcatReportBody(pc)" in APP_JS
    assert "function renderPharmcatAnalysisBody(pc)" in APP_JS
    assert "const reportView = pc.report_view || {};" in APP_JS
    assert "用藥建議概覽" in APP_JS
    assert "藥物建議摘要" in APP_JS
    assert "基因型與表現型" in APP_JS
    assert "完整用藥建議" in APP_JS
    assert "21 個 CPIC Level A 基因詳細結果" in APP_JS
    assert "function renderPharmcatBody(pc)" not in APP_JS


def test_pgx_main_screen_warns_when_json_is_missing_and_limits_tsv_fallback():
    assert "除 MT-RNR1 可由 TSV 補值外" in APP_JS
    assert ".pgx-source-warning" in STYLE_CSS
    assert ".pgx-report-table" in STYLE_CSS
