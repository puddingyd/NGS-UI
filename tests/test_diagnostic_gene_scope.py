"""WES exclusions affect diagnostic scope lists, not panel definitions."""
import json

import pytest
from docx import Document

from app.services import docx_export, phenotype_scorer, sample_loader


MITO_GENES = """
MT-ATP6 MT-ATP8 MT-CO1 MT-CO2 MT-CO3 MT-CYB
MT-ND1 MT-ND2 MT-ND3 MT-ND4 MT-ND4L MT-ND5 MT-ND6
MT-RNR1 MT-RNR2 MT-TA MT-TC MT-TD MT-TE MT-TF MT-TG MT-TH MT-TI MT-TK
MT-TL1 MT-TL2 MT-TM MT-TN MT-TP MT-TQ MT-TR MT-TS1 MT-TS2 MT-TT MT-TV MT-TW MT-TY
""".split()
EXCLUDED = MITO_GENES + [
    "SMN1", "SMN2", "STRC", "CYP21A2", "NCF1", "IKBKG", "OPN1LW", "OPN1MW", "MUC1",
]
RETAINED = ["MTOR", "MTFMT", "POLG", "TWNK", "PMS2", "PKD1", "RPGR", "FMR1", "GBA1"]
HPO = [{"phenotype": "HP:0001250", "label": "Seizure"}]
PANELS = [{"name": "WES-I__科別__Test panel"}]


@pytest.fixture
def scope_genes(monkeypatch):
    genes = EXCLUDED + RETAINED
    monkeypatch.setattr(phenotype_scorer, "load", lambda: None)
    monkeypatch.setattr(phenotype_scorer, "genes_for_key", lambda key, **kw: {"genes": genes})
    monkeypatch.setattr(phenotype_scorer, "panel_output_name", lambda name: "Test panel")
    monkeypatch.setattr(docx_export.panel_deadzone, "is_disease_associated_gene", lambda gene: True)
    return genes


@pytest.mark.parametrize("test_type", ["WES", " wes ", "WGS", "TITAN-WGS"])
@pytest.mark.parametrize("mode", ["grouped", "merged"])
def test_word_and_pdf_scope_lists_agree_and_preserve_other_text(scope_genes, test_type, mode):
    expected = sorted(RETAINED if test_type.strip().upper() == "WES" else scope_genes)
    sample = {"meta": {"Test": test_type}, "patient_phenotype": HPO, "selected_panels": PANELS}
    original = json.dumps(sample, sort_keys=True)

    payload = sample_loader._report_gene_list(HPO, PANELS, test_type=test_type)
    assert payload == {"grouped": [
        {"name": "Seizure", "genes": expected},
        {"name": "Test panel", "genes": expected},
    ], "merged": expected}

    doc = Document()
    docx_export._section_annotations(doc, sample, mode)
    paragraphs = [p.text for p in doc.paragraphs]
    baseline = Document()
    docx_export._section_annotations(baseline, {"meta": {"Test": "WGS"}}, mode)
    # All annotations before the gene list remain byte-for-byte identical.
    assert paragraphs[:9] == [p.text for p in baseline.paragraphs][:9]
    genes_text = "    " + ", ".join(expected)
    assert paragraphs[9:] == (["    Seizure, Test panel:", genes_text] if mode == "merged" else [
        "    Seizure:", genes_text, "", "    Test panel:", genes_text,
    ])
    assert json.dumps(sample, sort_keys=True) == original
    assert phenotype_scorer.genes_for_key("HP:0001250")["genes"] == scope_genes


@pytest.mark.parametrize("sample_id,meta,excluded", [
    ("SYNTHETIC-WES", {"test_type": "WES"}, True),
    ("SYNTHETIC-WGS", {"test_type": "WGS"}, False),
    ("26T-SYNTHETIC", {"test_type": "WES"}, False),
    ("SYNTHETIC-DEFAULT", {}, True),
])
def test_pdf_export_uses_sample_test_type(scope_genes, tmp_path, monkeypatch, sample_id, meta, excluded):
    metadata = tmp_path / "sample_metadata.json"
    metadata.write_text(json.dumps(meta))
    monkeypatch.setattr(sample_loader.sample_layout, "state_file", lambda *a: metadata)
    monkeypatch.setattr(sample_loader, "_load_pheno_context", lambda *a: (tmp_path, tmp_path, HPO, PANELS, {}))
    result = sample_loader.load_report_gene_list(sample_id)
    assert result["merged"] == sorted(RETAINED if excluded else scope_genes)


def test_filtered_empty_panel_keeps_existing_empty_list_wording(monkeypatch):
    monkeypatch.setattr(docx_export, "_genes_for_term_or_panel", lambda key: ["MT-ATP6", "SMN1"])
    monkeypatch.setattr(phenotype_scorer, "panel_output_name", lambda name: "Test panel")
    for mode in ["grouped", "merged"]:
        doc = Document()
        docx_export._render_gene_list(doc, {"meta": {"Test": "WES"}, "selected_panels": PANELS}, mode)
        assert [p.text for p in doc.paragraphs] == ["    Test panel:", "    （無對應基因）"]
