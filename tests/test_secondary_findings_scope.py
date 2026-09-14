from pathlib import Path

import pytest

from app.adapters.snv_tsv import classify_tier
from app.services import docx_export, phenotype_scorer, sample_loader


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
INDEX_HTML = (REPO_ROOT / "frontend" / "index.html").read_text(encoding="utf-8")


def test_secondary_findings_only_keep_requested_snv_panels():
    assert sample_loader.SECONDARY_SNV_PANELS == {
        "acmg_sf": "ACMG_SF_v3.3",
        "hereditary_cancer": "WES-I__腫瘤醫學__遺傳癌症 v2.0",
        "stroke": "WGS__神經科__Stroke",
        "carrier": "carrier_mackenzie_1300+",
    }

    defs = APP_JS.split("const SECONDARY_PANEL_DEFS = [", 1)[1].split("];", 1)[0]
    assert 'key: "acmg_sf"' in defs
    assert 'key: "hereditary_cancer"' in defs
    assert 'key: "stroke"' in defs
    assert 'key: "carrier"' in defs
    assert defs.index('key: "acmg_sf"') < defs.index('key: "hereditary_cancer"') < defs.index('key: "stroke"')
    for removed in ("lipid_fh", "proactive"):
        assert removed not in defs

    for removed_id in (
        "cat-lipid-fh-c",
        "cat-proactive-c",
        "sec-lipid-fh",
        "sec-proactive",
    ):
        assert removed_id not in INDEX_HTML
    for prefix, suffix in (("cat-", "-c"), ("sec-", "")):
        ids = [f'id="{prefix}{name}{suffix}"' for name in ("acmg-sf", "hereditary-cancer", "stroke", "carrier")]
        positions = [INDEX_HTML.index(element_id) for element_id in ids]
        assert positions == sorted(positions)
    assert 'data-target="cat-hereditary-cancer-c">遺傳癌症 v2.0</button>' in INDEX_HTML
    assert 'data-target="cat-carrier-c">Carrier screening</button>' in INDEX_HTML


def test_health_export_picker_has_five_requested_options_and_defaults():
    picker = APP_JS.split("function _pickHealthReportSections()", 1)[1].split(
        "return new Promise", 1,
    )[0]
    assert picker.count("key:") == 5
    assert '{ key: "acmg_sf", title: "ACMG 疾病風險基因（ACMG SF）", checked: true }' in picker
    assert '{ key: "hereditary_cancer", title: "遺傳癌症 v2.0", checked: false }' in picker
    assert '{ key: "stroke", title: "中風相關基因", checked: false }' in picker
    assert '{ key: "carrier", title: "帶因者篩查", checked: false }' in picker
    assert '{ key: "pgx", title: "藥物基因體學（PGx）", checked: true }' in picker
    assert picker.index('key: "acmg_sf"') < picker.index('key: "hereditary_cancer"') < picker.index('key: "stroke"')


def test_cancer_panel_uses_full_fixed_v2_gene_list_in_health_report(monkeypatch):
    panel_name = "WES-I__腫瘤醫學__遺傳癌症 v2.0"
    panels = phenotype_scorer._load_panels_from_dir(REPO_ROOT / "phenotype_data/gene_panels")[0]
    genes = panels[panel_name]
    assert len(genes) == 209
    assert {"BRCA1", "BRCA2", "CHEK2", "ATM"} <= genes
    requested_panels = []

    def lookup(name, kind):
        requested_panels.append((name, kind))
        return {"genes": sorted(panels[name])}

    monkeypatch.setattr(phenotype_scorer, "genes_for_key", lookup)
    sections = docx_export._health_panel_gene_sections({"hereditary_cancer"})
    assert requested_panels == [(panel_name, "panel")]
    assert sections == [("遺傳癌症 v2.0", sorted(genes))]


@pytest.mark.parametrize("changes,expected", [
    ({"CLNSIG": "Pathogenic"}, True),
    ({"CLNSIG": "Likely_pathogenic", "alt_af": 0.2}, True),
    ({"tier": "1A"}, True),
    ({"tier": "1B"}, True),
    ({"tier": "1C"}, True),
    ({"tier": "2"}, False),
    ({"tier": "1C", "alt_af": 0.19}, False),
    ({"CLNSIG": "Pathogenic", "alt_af": None}, False),
    ({"CLNSIG": "Pathogenic", "zygosity": "0/0"}, False),
    ({"CLNSIG": "Pathogenic", "gene_symbol": "NOT_IN_PANEL"}, False),
])
def test_all_four_secondary_panels_share_acmg_candidate_rules(monkeypatch, changes, expected):
    # Equal gene membership isolates the candidate rule from panel contents.
    monkeypatch.setattr(phenotype_scorer, "genes_for_key", lambda *args, **kwargs: {"genes": ["BRCA1"]})
    variant = {"gene_symbol": "BRCA1", "CLNSIG": "Benign", "tier": "2",
               "alt_af": 0.35, "zygosity": "Heterozygous", **changes}
    categories = sample_loader._build_secondary_snv_categories({"shared": variant})
    assert list(categories) == ["acmg_sf", "hereditary_cancer", "stroke", "carrier"]
    assert all(ids == (["shared"] if expected else []) for ids in categories.values())


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

    pknn_tier = classify_tier({"PKNN_LLR": "1"})
    assert pknn_tier == "1C"
    assert sample_loader._is_secondary_snv_candidate({
        **base,
        "tier": pknn_tier,
    }) is True

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
