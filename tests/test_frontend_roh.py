from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
APP_JS = (ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
INDEX_HTML = (ROOT / "frontend" / "index.html").read_text(encoding="utf-8")
STYLE_CSS = (ROOT / "frontend" / "style.css").read_text(encoding="utf-8")


def test_roh_is_a_lazy_side_channel_with_source_specific_rendering():
    assert "`/samples/${encodeURIComponent(sid)}/roh`" in APP_JS
    assert "function renderRohCard()" in APP_JS
    assert "DRAGEN large ROH SNVs" in APP_JS
    assert "開啟 AutoMap PDF" in APP_JS
    assert 'id="roh-content"' in INDEX_HTML


def test_roh_card_has_filters_ideogram_and_table_styles():
    assert 'data-roh-filter="${key}"' in APP_JS
    assert "ROH_GRCH38_CHROM_LENGTHS" in APP_JS
    assert ".roh-ideogram" in STYLE_CSS
    assert ".roh-table" in STYLE_CSS

