from pathlib import Path

from app.services import sample_loader


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
INDEX_HTML = (REPO_ROOT / "frontend" / "index.html").read_text(encoding="utf-8")


def test_secondary_findings_only_keep_requested_snv_panels():
    assert sample_loader.SECONDARY_SNV_PANELS == {
        "acmg_sf": "ACMG_SF_v3.3",
        "stroke": "WGS__神經科__Stroke",
        "carrier": "carrier_mackenzie_1300+",
    }

    defs = APP_JS.split("const SECONDARY_PANEL_DEFS = [", 1)[1].split("];", 1)[0]
    assert 'key: "acmg_sf"' in defs
    assert 'key: "stroke"' in defs
    assert 'key: "carrier"' in defs
    for removed in ("lipid_fh", "hereditary_cancer", "proactive"):
        assert removed not in defs

    for removed_id in (
        "cat-lipid-fh-c",
        "cat-hereditary-cancer-c",
        "cat-proactive-c",
        "sec-lipid-fh",
        "sec-hereditary-cancer",
        "sec-proactive",
    ):
        assert removed_id not in INDEX_HTML


def test_health_export_picker_has_four_requested_options_and_defaults():
    picker = APP_JS.split("function _pickHealthReportSections()", 1)[1].split(
        "return new Promise", 1,
    )[0]
    assert picker.count("key:") == 4
    assert '{ key: "acmg_sf", title: "ACMG 疾病風險基因", checked: true }' in picker
    assert '{ key: "stroke", title: "中風相關基因", checked: false }' in picker
    assert '{ key: "carrier", title: "帶因者篩查", checked: false }' in picker
    assert '{ key: "pgx", title: "藥物基因體學", checked: true }' in picker


def test_registration_status_and_analysis_queue_require_hpo_in_frontend():
    assert "分析已排入" not in APP_JS
    assert 'fd.set("run_analysis"' not in APP_JS
    assert 'const hasHpo  = Array.isArray(phenoEdit.hpo) && phenoEdit.hpo.length > 0;' in APP_JS
    assert 'if (!hasHpo)' in APP_JS
    assert 'if (!Array.isArray(phenoEdit.hpo) || !phenoEdit.hpo.length) return;' in APP_JS
    assert 'if (sampleInput) sampleInput.value = LIS_ID || "";' in APP_JS


def test_overlapping_secondary_variants_show_in_each_panel_with_global_status():
    assert "function _secondaryPanelsForVariant(id)" in APP_JS
    assert "_secondaryCanonicalPanel" not in APP_JS
    assert "Explicit dismissal wins" in APP_JS
    assert "panels.forEach(key =>" in APP_JS
    assert 'return ids.filter(id => _isSecondaryEligible(id));' in APP_JS
    assert 'function _syncVariantCheckboxes(selector, id, idx, checked, source = null)' in APP_JS


def test_secondary_candidates_use_main_snv_retrieval_tiers_but_default_clinvar_only():
    base = {"alt_af": 0.35, "zygosity": "Heterozygous", "CLNSIG": "Benign"}
    assert sample_loader._is_secondary_snv_candidate({**base, "tier": "1A"}) is True
    assert sample_loader._is_secondary_snv_candidate({**base, "tier": "1B"}) is True
    assert sample_loader._is_secondary_snv_candidate({**base, "tier": "1C"}) is True
    assert sample_loader._is_secondary_snv_candidate({**base, "tier": "2"}) is False
    assert sample_loader._is_secondary_snv_candidate({
        **base,
        "tier": "2",
        "CLNSIG": "Likely_pathogenic",
    }) is True
    assert sample_loader._is_secondary_snv_candidate({
        **base,
        "tier": "1C",
        "alt_af": 0.1,
    }) is False

    eligible = APP_JS.split("function _isSecondaryEligible(id)", 1)[1].split(
        "function _secondarySection", 1,
    )[0]
    selected = APP_JS.split("function isSecondarySelected(id, panel)", 1)[1].split(
        "function getPanelStatus", 1,
    )[0]
    assert '["1A", "1B", "1C"].includes' in eligible
    assert "return _isClinvarPlp(v);" in selected


def test_tertiary_log_height_is_about_122_percent_of_original():
    style = (REPO_ROOT / "frontend" / "style.css").read_text(encoding="utf-8")
    assert "#dragen-job-log {\n  max-height: 342px;\n}" in style
