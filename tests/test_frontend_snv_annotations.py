from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
STYLE_CSS = (REPO_ROOT / "frontend" / "style.css").read_text(encoding="utf-8")
INDEX_HTML = (REPO_ROOT / "frontend" / "index.html").read_text(encoding="utf-8")


def test_core_snv_fields_have_annotation_hints():
    assert 'const scoreHint = _annotationHint("Score"' in APP_JS
    assert 'ClinVar (2026-07-20)${clinvarExternalLink}${clinvarChange}' in APP_JS
    assert 'const clinvarHint = _annotationHint("ClinVar"' not in APP_JS
    assert 'class="v acmg-summary-btn js-acmg-open"' in APP_JS
    assert '["manual", "Manual"]' in APP_JS
    assert '["erepo", "ERepo"]' in APP_JS
    assert '["genebe", "GeneBe"]' in APP_JS
    assert '["inhouse", "In-house"]' in APP_JS
    assert "ClinGen VCEP experts 評估" in APP_JS


def test_all_displayed_in_silico_tools_have_calibration_metadata():
    expected = {
        "pknn", "alphamissense", "pangolin", "esm1b", "varity", "bayesdel",
        "revel", "spliceai", "cadd", "mutpred2", "vest4", "dann",
        "phactboost", "phylop", "gerp", "sift", "loftool",
    }
    for key in expected:
        assert f"  {key}: {{" in APP_JS
        assert f'case "{key}"' in APP_JS
    assert '"PP3_Strong: ≥ 0.990"' in APP_JS
    assert '"PP3_Strong: ≥ 0.932"' in APP_JS
    assert '"BP4_Strong: ≤ 0.036"' in APP_JS
    assert '"PP3: ≥ 0.20"' in APP_JS
    assert "Reference: <a href=" in APP_JS
    assert "PMID:" in APP_JS


def test_reference_tooltips_are_clickable_and_keyboard_accessible():
    assert 'data-tip-html=' in APP_JS
    assert '<button type="button" class="${className}"' in APP_JS
    assert 'target="_blank" rel="noopener"' in APP_JS
    assert 'tipEl.addEventListener("mouseenter", cancelHide)' in APP_JS
    assert "pointer-events: auto" in STYLE_CSS


def test_clinvar_card_keeps_fixed_baseline_label_without_info_icon():
    assert '<span class="k">ClinVar (2026-07-20)${clinvarExternalLink}${clinvarChange}</span>' in APP_JS
    assert '_annotationHint("ClinVar"' not in APP_JS


def test_clinvar_external_link_prefers_baseline_then_latest_and_hides_without_id():
    helper = APP_JS[
        APP_JS.index("function _normalizeClinvarVariationId"):
        APP_JS.index("function _litvar2TitleHtml")
    ]
    assert 'Object.prototype.hasOwnProperty.call(' in helper
    assert '"clinvar_variation_id_old"' in helper
    assert 'v?.clinvar_variation_id' in helper
    assert 'v?.clinvar_latest_variation_id' in helper
    assert 'const variationId = baselineId || latestId;' in helper
    assert 'if (!variationId) return "";' in helper
    assert 'baselineId ? "在 ClinVar 開啟" : "在最新版 ClinVar 開啟"' in helper
    assert 'https://www.ncbi.nlm.nih.gov/clinvar/variation/${encodeURIComponent(variationId)}/' in helper
    assert 'class="clinvar-external-link"' in helper
    assert 'target="_blank" rel="noopener"' in helper
    assert 'EXTERNAL_LINK_ICON_SVG' in helper
    assert '.clinvar-external-link' in STYLE_CSS


def test_reviewer_requested_predictor_annotation_and_coloring_rules():
    assert '"Uncertain: > -1 to < 1"' in APP_JS
    assert '"PP3_Supporting: 1 to < 2"' in APP_JS
    assert '"PP3_Moderate: 2 to < 4"' in APP_JS
    assert '"PP3_Strong: ≥ 4"' in APP_JS
    assert '"BP4_Supporting: > -2 to ≤ -1"' in APP_JS
    assert '"BP4_Moderate: > -4 to ≤ -2"' in APP_JS
    assert '"BP4_Strong: ≤ -4"' in APP_JS
    assert 'case "pknn": return _pknnEvidence(v.PKNN_evidence);' in APP_JS
    assert 'return _evidence("sig-vus", value || "Uncertain");' in APP_JS
    assert 'case "dann": return _evidence("", "No calibrated evidence");' in APP_JS

    loftee_start = APP_JS.index('extras.push({ key: "LOFTEE"')
    loftee_hint = APP_JS.index('hint: _annotationHint("LOFTEE"', loftee_start)
    assert "cls:" not in APP_JS[loftee_start:loftee_hint]


def test_reviewer_requested_annotation_copy_is_streamlined():
    assert "Variant score：ACMG score 轉換成 0–100" in APP_JS
    assert "PKNN_LLR 僅顯示數值" not in APP_JS
    assert "目前文獻為 preprint" not in APP_JS
    assert "因此所有數值以黃色標示為 contextual evidence" not in APP_JS
    assert "這是 gene-level intolerance，不是 variant-level PP3/BP4 evidence" not in APP_JS
    assert "Total score 可超過 100；手動修改 ACMG 不會即時重算排序 score" not in APP_JS


def test_predictor_display_order_and_primary_count_match_review_workflow():
    ordered = [
        'key: "pknn"', 'key: "alphamissense"', 'key: "pangolin"',
        'key: "revel"', 'key: "spliceai"', 'key: "esm1b"',
        'key: "varity"', 'key: "bayesdel"', 'key: "cadd"',
        'key: "dann"', 'key: "mutpred2"', 'key: "vest4"',
        'key: "phactboost"', 'key: "phylop"',
        'key: "gerp"', 'key: "sift"', 'key: "loftool"',
    ]
    positions = [APP_JS.index(token, APP_JS.index("const IN_SILICO_TOOLS")) for token in ordered]
    assert positions == sorted(positions)
    assert "const IN_SILICO_PRIMARY_COUNT = 3;" in APP_JS


def test_clinvar_change_arrows_only_use_clinically_meaningful_boundary():
    assert 'direction === "UP_TO_PLP"' in APP_JS
    assert 'direction === "DOWN_FROM_PLP"' in APP_JS
    assert '${isUp ? "↑" : "↓"}' in APP_JS
    assert '最新版分類：${latestClassification}' in APP_JS
    assert 'Review status：${reviewStatus}' in APP_JS
    assert '最新版日期：${latestClinvarDate || "未提供"}' in APP_JS
    assert "v.clinvar_latest_review_status" in APP_JS
    assert '.clinvar-change.upgrade' in STYLE_CSS
    assert '.clinvar-change.downgrade' in STYLE_CSS


def test_research_only_checkbox_controls_academic_dbnsfp_and_spliceai():
    assert '<span>Research-only</span>' in INDEX_HTML
    assert 'with_research_only: !!extra?.checked' in APP_JS
    assert 'academic dbNSFP 5.3a' in INDEX_HTML
    assert 'post-processing 另補 SpliceAI' in INDEX_HTML


def test_predicted_suspect_tier_labels_without_trigger_badges():
    assert 'const TIER_ORDER = ["1A", "1B", "1C", "2"];' in APP_JS
    assert '"1C": "1C — Predicted suspect"' in APP_JS
    assert '"2":  "2 — Other"' in APP_JS
    assert "ClinVar P/LP 0★ or CONF" not in APP_JS
    assert "badge-suspect" not in APP_JS
    assert "Extra-VEP rescue" not in APP_JS
