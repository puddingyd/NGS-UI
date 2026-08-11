from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
STYLE_CSS = (REPO_ROOT / "frontend" / "style.css").read_text(encoding="utf-8")
FLOW_MD = (REPO_ROOT / "frontend" / "ANALYSIS_FLOW.md").read_text(encoding="utf-8")


def test_snv_prioritization_splits_editable_filter_from_vertical_tiers():
    assert "- review_filter | AF < 0.01 · VAF > 0.2 · ClinVar rescue" in FLOW_MD
    assert 'analysisFlowValue(flow, name, "review_filter")' in APP_JS
    assert '<div class="analysis-flow-review-filter-title">Filter</div>' in APP_JS
    assert ">Review TSV filter</div>" not in APP_JS
    assert '<div class="analysis-flow-tier-strip-title">Category</div>' in APP_JS
    assert '"1B": "1B — Loss-of-function"' in APP_JS
    assert '"Other": "2 - Other"' in APP_JS
    assert 'class="analysis-flow-snv-priority-arrow"' in APP_JS
    assert ".analysis-flow-snv-priority-flow" in STYLE_CSS
    assert "grid-template-columns: minmax(0, 0.85fr) auto minmax(0, 1.55fr);" in STYLE_CSS
    assert ".analysis-flow-tier-strip {\n  display: grid;\n  gap: 4px;" in STYLE_CSS
    assert "border-top: 1px solid var(--border);" in STYLE_CSS
