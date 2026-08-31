from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
APP = (ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
HTML = (ROOT / "frontend" / "index.html").read_text(encoding="utf-8")
CSS = (ROOT / "frontend" / "style.css").read_text(encoding="utf-8")


def test_acmg_card_is_modal_only_and_shows_source_next_to_result():
    card_block = APP[APP.index("function renderVariantCard"):APP.index(
        "// ---- helpers used by renderVariantCard"
    )]
    assert 'class="v acmg-summary-btn js-acmg-open"' in card_block
    assert 'class="acmg-summary-source"' in card_block
    assert '<textarea class="acmg-crit"' not in card_block
    assert '<select class="acmg-class' not in card_block
    assert '<input class="acmg-score"' not in card_block
    assert "editAcmgDisplayClass" in card_block
    summary_css = CSS[CSS.index(".acmg-summary-btn"):CSS.index(".acmg-modal-card")]
    assert "border: 0;" in summary_css
    assert "white-space: nowrap;" in summary_css
    assert "line-height: 1.5;" in summary_css
    assert "min-height: 30px" not in summary_css
    assert "border-left" not in summary_css


def test_vus_subtiers_keep_formal_class_and_have_directional_colors():
    assert 'return "VUS-low"' in APP
    assert 'return "VUS-mid"' in APP
    assert 'return "VUS-high"' in APP
    assert "VUS-low（0–1）" in HTML
    assert "正式分類仍為 Uncertain significance" in HTML
    low = CSS[CSS.index(".sig-vus-low"):CSS.index(".sig-vus-mid")]
    high = CSS[CSS.index(".sig-vus-high"):CSS.index(".sig-lb")]
    assert "linear-gradient(90deg" in low
    assert "#b9dfbd" in low
    assert "linear-gradient(90deg" in high
    assert "#ef9a9a" in high


def test_acmg_save_restores_the_same_variant_after_backend_reordering():
    assert "captureAcmgCardOrigin(trigger)" in APP
    assert "card.dataset.variantId = id;" in APP
    assert "activeTierTab = updatedTier;" in APP
    assert 'blockBody.classList.add("open")' in APP
    assert "restoreAcmgSavedCard(editedId, origin, updatedTier);" in APP
    assert "acmg-just-saved" not in APP
    assert ".variant-card.acmg-just-saved" not in CSS


def test_acmg_modal_has_three_apply_sources_and_all_criteria_editor():
    acmg_modal = HTML[HTML.index('id="acmg-modal"'):HTML.index('id="observed-modal"')]
    assert 'id="acmg-modal"' in acmg_modal
    assert "Manual ACMG/AMP variant classification" in acmg_modal
    assert 'id="acmg-source-summaries"' in HTML
    assert 'id="acmg-criteria-list"' in HTML
    assert acmg_modal.count('class="btn btn-primary js-acmg-save"') == 2
    assert acmg_modal.count('data-close="acmg-modal"') == 2
    assert "全部關閉" not in acmg_modal
    assert acmg_modal.index('id="acmg-criteria-list"') < acmg_modal.index(
        "分類採 Tavtigian natural-scale points"
    )
    assert 'renderColumn("pathogenic", "Pathogenic criteria")' in APP
    assert 'renderColumn("benign", "Benign criteria")' in APP
    assert "data-evidence-group=" in APP
    assert ".acmg-criteria-column.pathogenic" in CSS
    assert ".acmg-criteria-column.benign" in CSS
    assert '["manual", "Manual"]' in APP
    assert '["genebe", "GeneBe"]' in APP
    assert '["inhouse", "In-house"]' in APP
    assert "ClinGen SVI recommends retiring" in (
        ROOT / "backend" / "app" / "services" / "manual_acmg.py"
    ).read_text(encoding="utf-8")


def test_acmg_source_summaries_always_show_evidence_strength():
    helper_start = APP.index("function formatAcmgCriteriaSummary(source, catalog)")
    helper_end = APP.index("function calculateAcmgPreview", helper_start)
    helper = APP[helper_start:helper_end]
    summaries_start = APP.index("function renderAcmgSourceSummaries()")
    summaries_end = APP.index("function renderAcmgCriteriaEditor()", summaries_start)
    summaries = APP[summaries_start:summaries_end]

    assert "parseAcmgCriteriaText" in helper
    assert "strengthLabels" in helper
    assert "`${code}_${label}`" in helper
    assert '.replace(/[\\s-]+/g, "_")' in helper
    assert "formatAcmgCriteriaSummary(source, _acmgCatalog)" in summaries
    assert "reusable_criteria_text" in summaries
    assert "formatAcmgCriteriaSummary({" in summaries


def test_observed_is_a_separate_badge_and_modal():
    assert 'Observed (${Number(v.observed_count)})' in APP
    assert 'id="observed-modal"' in HTML
    acmg_modal = HTML[HTML.index('id="acmg-modal"'):HTML.index('id="observed-modal"')]
    assert "Observed in other cases" not in acmg_modal
    assert ".badge-observed" in CSS
