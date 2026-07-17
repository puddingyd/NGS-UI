from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
STYLE_CSS = (REPO_ROOT / "frontend" / "style.css").read_text(encoding="utf-8")


def test_core_snv_fields_have_annotation_hints():
    assert 'const scoreHint = _annotationHint("Score"' in APP_JS
    assert 'const clinvarHint = _annotationHint("ClinVar"' in APP_JS
    assert 'source = "manual"' in APP_JS
    assert 'source = "GeneBe"' in APP_JS
    assert 'source = "in-house"' in APP_JS
    assert "in-house: ${_shortAcmgClass(inHouse)}" in APP_JS
    assert "GeneBe: ${_shortAcmgClass(geneBe)}" in APP_JS


def test_all_displayed_in_silico_tools_have_calibration_metadata():
    expected = {
        "pknn", "alphamissense", "pangolin", "esm1b", "varity", "bayesdel",
        "revel", "spliceai", "metarnn", "dann", "phactboost", "phylop", "gerp",
        "sift", "loftool",
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


def test_missing_clinvar_version_is_explicit():
    assert "version date: 三級輸出未提供" in APP_JS
    assert "annotation_versions.json sidecar" in APP_JS


def test_predictor_display_order_and_primary_count_match_review_workflow():
    ordered = [
        'key: "pknn"', 'key: "alphamissense"', 'key: "pangolin"',
        'key: "revel"', 'key: "spliceai"', 'key: "esm1b"',
        'key: "varity"', 'key: "bayesdel"', 'key: "metarnn"',
        'key: "dann"', 'key: "phactboost"', 'key: "phylop"',
        'key: "gerp"', 'key: "sift"', 'key: "loftool"',
    ]
    positions = [APP_JS.index(token, APP_JS.index("const IN_SILICO_TOOLS")) for token in ordered]
    assert positions == sorted(positions)
    assert "const IN_SILICO_PRIMARY_COUNT = 3;" in APP_JS


def test_predicted_suspect_tier_labels_without_trigger_badges():
    assert 'const TIER_ORDER = ["1A", "1B", "1C", "2"];' in APP_JS
    assert '"1C": "1C — Predicted suspect"' in APP_JS
    assert '"2":  "2 — Other"' in APP_JS
    assert "ClinVar P/LP 0★ or CONF" not in APP_JS
    assert "badge-suspect" not in APP_JS
    assert "Extra-VEP rescue" not in APP_JS
