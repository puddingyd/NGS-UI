from docx import Document

from app.services import docx_export, phenotype_scorer, sample_loader


def test_fixed_panel_report_name_is_short_and_old_key_is_compatible():
    phenotype_scorer.reload_db()

    new_name = "WES-I__腫瘤醫學__遺傳癌症 v2.0"
    old_name = "WES-I__腫瘤醫學__遺傳癌症"
    assert phenotype_scorer.panel_output_name(new_name) == "遺傳癌症 v2.0"
    assert phenotype_scorer.panel_output_name(old_name) == "遺傳癌症 v2.0"
    assert phenotype_scorer.genes_for_key(old_name, kind="panel")["genes"]
    assert phenotype_scorer.normalize_panel_entries(
        [{"name": old_name, "weight": 2}],
    ) == [{"name": new_name, "weight": 2}]


def test_custom_panel_metadata_controls_report_name_and_source():
    phenotype_scorer.reload_db()

    panel = phenotype_scorer.genes_for_key("ARVC_panelapp", kind="panel")
    assert panel["source"] == "PanelApp"
    assert panel["output_name"] == "ARVC"
    assert phenotype_scorer.panel_output_name("ARVC_panelapp") == "ARVC"


def test_pdf_gene_list_payload_uses_panel_output_name(monkeypatch):
    monkeypatch.setattr(phenotype_scorer, "load", lambda: (0, 0))
    monkeypatch.setattr(
        phenotype_scorer,
        "genes_for_key",
        lambda key, kind="": {"genes": ["BRCA1"]},
    )
    monkeypatch.setattr(
        phenotype_scorer,
        "panel_output_name",
        lambda name: "遺傳癌症 v2.0",
    )
    monkeypatch.setattr(
        sample_loader.panel_deadzone,
        "canonical_panel_gene_symbol",
        lambda gene: gene,
    )
    monkeypatch.setattr(
        sample_loader.panel_deadzone,
        "is_disease_associated_gene",
        lambda gene: True,
    )

    payload = sample_loader._report_gene_list(
        [], [{"name": "WES-I__腫瘤醫學__遺傳癌症 v2.0"}],
    )
    assert payload["grouped"] == [
        {"name": "遺傳癌症 v2.0", "genes": ["BRCA1"]},
    ]


def test_diagnosis_docx_gene_list_uses_panel_output_name(monkeypatch):
    monkeypatch.setattr(
        docx_export,
        "_genes_for_term_or_panel",
        lambda key: ["BRCA1"],
    )
    monkeypatch.setattr(
        phenotype_scorer,
        "panel_output_name",
        lambda name: "遺傳癌症 v2.0",
    )
    doc = Document()
    docx_export._render_gene_list(
        doc,
        {
            "patient_phenotype": [],
            "selected_panels": [
                {"name": "WES-I__腫瘤醫學__遺傳癌症 v2.0"},
            ],
            "meta": {"Test": "WES"},
        },
        "grouped",
    )

    assert [paragraph.text for paragraph in doc.paragraphs] == [
        "遺傳癌症 v2.0:",
        "BRCA1",
    ]
