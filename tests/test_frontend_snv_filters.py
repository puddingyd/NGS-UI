from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
INDEX_HTML = (REPO_ROOT / "frontend" / "index.html").read_text(encoding="utf-8")


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
