from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
STYLE_CSS = (REPO_ROOT / "frontend" / "style.css").read_text(encoding="utf-8")
INDEX_HTML = (REPO_ROOT / "frontend" / "index.html").read_text(encoding="utf-8")


def test_pgx_main_screen_uses_shared_report_projection():
    assert "function renderPharmcatReportBody(pc)" in APP_JS
    assert "function renderPharmcatAnalysisBody(pc)" in APP_JS
    assert "const reportView = pc.report_view || {};" in APP_JS
    assert "用藥建議概覽" in APP_JS
    assert "藥物建議摘要" in APP_JS
    assert "基因型與表現型" in APP_JS
    assert "完整用藥建議" in APP_JS
    assert "CPIC Level A 基因詳細結果" in APP_JS
    assert "21 個 CPIC Level A 基因詳細結果" not in APP_JS
    assert "analysis_action_categories" in APP_JS
    assert "可對應之 star alleles" in APP_JS
    assert "不代表病人同時具有全部 alleles" in APP_JS
    assert "pgx-gene-genotype" in APP_JS
    assert '["藥物", "CPIC/FDA 建議"]' in APP_JS
    assert 'class="pgx-star-alleles-note"' not in APP_JS
    assert "function renderPharmcatBody(pc)" not in APP_JS


def test_pgx_main_screen_warns_when_json_is_missing_and_limits_tsv_fallback():
    assert "除 MT-RNR1 可由 TSV 補值外" in APP_JS
    assert ".pgx-source-warning" in STYLE_CSS
    assert ".pgx-report-table" in STYLE_CSS
    assert ".pgx-gene-head-attention" in STYLE_CSS
    assert ".pgx-gene-genotype" in STYLE_CSS
    assert ".pgx-star-alleles-header" in STYLE_CSS


def test_secondary_cleanup_is_a_dgx_command_not_a_ui_delete():
    assert "顯示 DGX2 清理指令" in INDEX_HTML
    assert "secondary-clean-result-panel" in INDEX_HTML
    assert 'apiFetch("/secondary/nf-work/cleanup-command")' in APP_JS
    assert 'apiFetch("/secondary/nf-work/cleanup", { method: "POST" })' not in APP_JS
    assert "function _secondaryWriteClipboard(text)" in APP_JS
    assert 'document.execCommand("copy")' in APP_JS
    assert 'prompt("請複製以下 DGX2 清理指令' not in APP_JS
