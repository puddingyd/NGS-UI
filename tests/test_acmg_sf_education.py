import io
from pathlib import Path
import re

from docx import Document
from docx.oxml.ns import qn
import pytest

from app.services import acmg_sf_education, docx_export


def test_catalogue_covers_the_report_panel_and_cites_every_disease():
    catalogue = acmg_sf_education.load_catalogue()
    panel = Path(__file__).resolve().parents[1] / "phenotype_data/custom_panels/ACMG_SF_v3.3.txt"
    panel_genes = {line.strip() for line in panel.read_text().splitlines()
                   if line.strip() and not line.startswith("#")}
    conditions = catalogue["conditions"]
    genes = {gene for condition in conditions for gene in condition["genes"]}
    assert len(conditions) == 38
    assert len(genes) == 84
    assert genes == panel_genes
    assert genes == {gene for group in docx_export._ACMG_SF_GROUPS for gene in group["genes"]}
    assert len({c["id"] for c in conditions}) == len(conditions)
    category_ids = [category["id"] for category in catalogue["categories"]]
    assert list(dict.fromkeys(c["category"] for c in conditions)) == category_ids
    for condition in conditions:
        assert re.fullmatch(r"[a-z][a-z0-9_]*", condition["id"])
        assert condition["genes"] and len(set(condition["genes"])) == len(condition["genes"])
        for field in ("title", "english", "inheritance", "clinical_course", "management", "notes"):
            assert condition[field].strip()
        assert condition["references"]
        for key in condition["references"]:
            source = catalogue["sources"][key]
            assert source["title"] and source["publisher"]
            assert source["url"].startswith("https://")


def test_education_index_has_working_links_and_no_patient_result_labels():
    doc = Document()
    normal_xml = doc.styles["Normal"].element.xml
    acmg_sf_education.render_acmg_sf_education(doc)
    # Round-trip through a DOCX package, including its hyperlink relationships.
    buffer = io.BytesIO()
    doc.save(buffer)
    doc = Document(io.BytesIO(buffer.getvalue()))
    assert doc.styles["Normal"].element.xml == normal_xml
    catalogue = acmg_sf_education.load_catalogue()
    assert len(doc.tables) == len(catalogue["categories"])
    assert sum(len(table.rows) - 1 for table in doc.tables) == 38
    assert len(doc.element.xpath(".//w:tblHeader")) == 6
    for table in doc.tables:
        assert [cell.text for cell in table.rows[0].cells] == ["編號", "相關疾病", "基因", "遺傳模式"]
        assert all(row._tr.find(qn("w:trPr")).find(qn("w:cantSplit")) is not None for row in table.rows)
    bookmarks = doc.element.xpath(".//w:bookmarkStart/@w:name")
    assert len(bookmarks) == len(set(bookmarks))
    assert set(doc.element.xpath(".//w:hyperlink/@w:anchor")) <= set(bookmarks)
    assert {f"acmgsf_{c['id']}" for c in catalogue["conditions"]} <= set(bookmarks)
    external_urls = {rel.target_ref for rel in doc.part.rels.values() if rel.is_external}
    assert {source["url"] for source in catalogue["sources"].values()} == external_urls
    text = "\n".join(doc.element.xpath(".//w:t/text()"))
    for status in ("檢出", "未檢出", "Positive", "Negative", "No Result"):
        assert status not in text
    for disease in catalogue["conditions"]:
        assert disease["clinical_course"] in text
        assert disease["management"] in text
    assert "p.Cys282Tyr" in text
    assert "少數僅有一份 CASQ2" in text
    assert "半顯性" in text


@pytest.mark.parametrize("sections,with_finding", [
    (["acmg_sf"], False),
    (["acmg_sf", "pgx"], False),
    (["acmg_sf", "pgx"], True),
    (["pgx"], False),
    (["stroke", "pgx"], True),
    (["carrier"], False),
])
def test_health_export_selects_and_orders_education(monkeypatch, tmp_path, sections, with_finding):
    from app.services import sample_layout

    variant = {"id": "test-variant", "gene_symbol": "LDLR", "HGVS_C": "c.1A>G",
               "ACMG_classification": "Pathogenic", "Zygosity": "Heterozygous"}
    variants = {variant["id"]: variant} if with_finding else {}
    categories = {key: list(variants) for key in docx_export._HEALTH_DISEASE_SECTIONS}
    monkeypatch.setattr(docx_export.sample_loader, "load_sample", lambda *args, **kwargs: {"meta": {"Test": "WGS"}})
    monkeypatch.setattr(docx_export.sample_loader, "load_sample_secondary_snv",
                        lambda *args, **kwargs: {"variants": variants, "categories": categories})
    monkeypatch.setattr(docx_export.sample_loader, "load_sample_pgx", lambda *args: {})
    monkeypatch.setattr(docx_export.report_store, "load", lambda *args: {
        "edits": {}, "secondary_findings": {key: {"selected": list(variants)} for key in categories},
    })
    monkeypatch.setattr(sample_layout, "state_dir", lambda *args: tmp_path)
    monkeypatch.setattr(docx_export.phenotype_scorer, "genes_for_key", lambda *args, **kwargs: {"genes": ["LDLR"]})
    groups = [{"drug": "Clopidogrel", "genes": {"CYP2C19": {"phenotype": "Poor Metabolizer"}},
               "recommendations": [{"source": "CPIC", "level": "Strong", "recommendation": "Use an alternative antiplatelet agent."}]}]
    # Only the already-tested PGx presenter is substituted; appendix assembly,
    # variant references, education and full recommendation rendering are real.
    monkeypatch.setattr(docx_export, "_render_health_pgx_section", lambda *args: groups)
    payload = docx_export.build_health_docx("education-test", sections=sections)
    doc = Document(io.BytesIO(payload))
    text = "\n".join(doc.element.xpath(".//w:t/text()"))
    title = acmg_sf_education.load_catalogue()["title"]
    assert (title in text) == ("acmg_sf" in sections)
    if "acmg_sf" in sections:
        assert len(doc.tables) == 6
        assert "38 個疾病群組" in text
        if with_finding:
            assert text.index("變異位點參考資料") < text.index(title)
        else:
            assert "變異位點參考資料" not in text
        if "pgx" in sections:
            assert text.index(title) < text.index("ACMG SF 疾病簡介參考資料") < text.index("完整用藥建議")
            assert "Use an alternative antiplatelet agent." in text
    if "pgx" in sections and with_finding and "acmg_sf" not in sections:
        assert text.index("變異位點參考資料") < text.index("完整用藥建議")
    if not with_finding and "acmg_sf" not in sections and "pgx" not in sections:
        assert "附錄" not in text
