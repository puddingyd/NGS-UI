from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
INDEX_HTML = (REPO_ROOT / "frontend" / "index.html").read_text(encoding="utf-8")
STYLE_CSS = (REPO_ROOT / "frontend" / "style.css").read_text(encoding="utf-8")


def _function_body(name: str) -> str:
    marker = f"function {name}("
    start = APP_JS.index(marker)
    # Some functions destructure an options object in their parameter list;
    # use the signature's closing `) {` rather than the first `{`.
    brace = APP_JS.index(") {", start) + 2
    depth = 0
    for index in range(brace, len(APP_JS)):
        char = APP_JS[index]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return APP_JS[brace + 1:index]
    raise AssertionError(f"unterminated JavaScript function: {name}")


def test_main_snv_filter_does_not_apply_gnomad_af_again():
    body = _function_body("_passesMainSnvDisplayFilters")
    setup_body = _function_body("setupSnvDisplayFilters")

    assert 'id="filter-gnomad-af"' not in INDEX_HTML
    assert "filter-gnomad-af" not in body
    assert "filter-gnomad-af" not in setup_body
    # The independent raw-TSV gene-search filter remains available.
    assert 'id="gene-search-filter-gnomad-af"' in INDEX_HTML


def test_modifier_filter_rescues_clinvar_plp():
    body = " ".join(_function_body("_passesMainSnvDisplayFilters").split())

    assert 'String(v.impact || "").toUpperCase() === "MODIFIER"' in body
    assert "&& !_isClinvarPlp(v)) return false;" in body


def test_nckuh_common_filter_defaults_hidden_for_snv_and_mito():
    checkbox = '<input id="filter-nckuh-common" type="checkbox" />'
    assert checkbox in INDEX_HTML
    assert INDEX_HTML.index('id="filter-in-panel-only"') < INDEX_HTML.index(
        'id="filter-nckuh-common"'
    ) < INDEX_HTML.index('id="filter-vaf"')
    assert "AF_NCKUH ≥ 0.05 &amp; AC ≥ 50" in INDEX_HTML
    assert "const NCKUH_COMMON_AF_THRESHOLD = 0.05;" in APP_JS
    assert "const NCKUH_COMMON_AC_THRESHOLD = 50;" in APP_JS

    snv_body = _function_body("_passesMainSnvDisplayFilters")
    mito_body = _function_body("_mitoIdsForTier")
    assert '_passesNckuhCommonFilter(v, id, "snv")' in snv_body
    assert '"mito"' in mito_body


def test_nckuh_common_rescues_known_pathogenic_and_primary_review_status_only():
    body = _function_body("_isNckuhCommonRescued")

    assert "_hasPrimaryReviewerStatus(id)" in body
    assert "_isClinvarPlp(v)" in body
    assert "mitomap_pathogenic" in body
    assert "mitomap_reported" in body
    assert "effective_acmg_class" in body
    assert "zygosity" not in body.lower()
    assert "recessive" not in body.lower()


def test_scope_labels_counts_and_filter_heights_match_visible_cards():
    assert "<span>疾病相關</span>" in INDEX_HTML
    assert "<span>臨床相關</span>" in INDEX_HTML
    count_body = _function_body("updateInPanelCount")
    assert "ignoreInPanelOnly: true" in count_body
    assert "ignoreDiseaseAssociated: true" in count_body
    assert "min-height: 27px;" in STYLE_CSS
