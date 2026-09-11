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
    assert len(conditions) == 36
    assert len(genes) == 84
    assert genes == panel_genes
    assert genes == {gene for group in docx_export._ACMG_SF_GROUPS for gene in group["genes"]}
    assert len({c["id"] for c in conditions}) == len(conditions)
    category_ids = [category["id"] for category in catalogue["categories"]]
    assert list(dict.fromkeys(c["category"] for c in conditions)) == category_ids
    for condition in conditions:
        assert re.fullmatch(r"[a-z][a-z0-9_]*", condition["id"])
        assert condition["genes"] and len(set(condition["genes"])) == len(condition["genes"])
        for field in ("title", "english", "inheritance", "index_inheritance", "clinical_course", "management"):
            assert condition[field].strip()
        assert condition["references"]
        for key in condition["references"]:
            source = catalogue["sources"][key]
            assert source["title"] and source["publisher"]
            assert source["url"].startswith("https://")
    assert len(catalogue["source_order"]) == len(set(catalogue["source_order"])) == 45
    assert set(catalogue["source_order"]) == set(catalogue["sources"])


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
    assert sum(len(table.rows) - 1 for table in doc.tables) == 36
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


def test_reviewed_index_and_anchored_comments_are_applied():
    doc = Document()
    acmg_sf_education.render_acmg_sf_education(doc)
    catalogue = acmg_sf_education.load_catalogue()
    by_id = {c["id"]: c for c in catalogue["conditions"]}
    paragraphs = [p.text for p in doc.paragraphs]
    text = "\n".join(doc.element.xpath(".//w:t/text()"))
    assert paragraphs[0] == "ACMG疾病風險基因與相關疾病簡介"
    assert paragraphs[paragraphs.index("疾病索引") + 1] == catalogue["reading_note"]
    for removed in ("如何閱讀遺傳模式", "內容更新：", "疾病短文", "資料來源：",
                    "文中編號對應", "息肉", "半顯性", "p.Cys282Tyr", "少數僅有一份 CASQ2"):
        assert removed not in text
    assert paragraphs.count("其他疾病") == 2
    assert "ACMG SF 疾病簡介參考資料" in paragraphs
    for title in ("疾病介紹", "ACMG SF 疾病簡介參考資料"):
        assert next(p for p in doc.paragraphs if p.text == title).paragraph_format.page_break_before
    assert not doc.element.xpath(".//w:br[@w:type='page']")
    assert {c["id"] for c in catalogue["conditions"] if c["notes"]} == {"fh", "pgl"}
    assert sum(p.startswith("補充說明：") for p in paragraphs) == 2
    assert "瘜肉相關問題與癌症風險需分別考慮" in by_id["pjs"]["management"]
    assert not by_id["pjs"]["notes"]

    merged = by_id["polyposis"]
    assert set(merged["genes"]) == {"APC", "MUTYH", "BMPR1A", "SMAD4"}
    assert not {"apc", "mutyh", "jps"} & set(by_id)
    for title in ("家族性腺瘤性瘜肉症", "MUTYH 相關瘜肉症", "幼年型瘜肉症候群"):
        assert title in merged["clinical_course"]
    assert "APC、BMPR1A、SMAD4：體染色體顯性遺傳" in merged["inheritance"]
    assert "MUTYH：體染色體隱性遺傳" in merged["inheritance"]
    assert "出現深褐色的小斑點" in by_id["pjs"]["clinical_course"]
    assert "女性也可能在成年後出現雙腿僵硬、走路困難，或難以控制排尿、排便等症狀" in by_id["ald"]["clinical_course"]
    assert by_id["fh"]["index_inheritance"] == "體染色體顯性、體染色體隱性"
    assert by_id["pgl"]["index_inheritance"] == "體染色體顯性"
    for cid in ("ald", "fabry", "otc"):
        assert by_id[cid]["index_inheritance"] == by_id[cid]["inheritance"] == "X 染色體性聯遺傳"

    index_rows = [row for table in doc.tables for row in table.rows[1:]]
    for number, (condition, row) in enumerate(zip(catalogue["conditions"], index_rows), 1):
        assert row.cells[0].text == f"{number:02}"
        assert "".join(row.cells[1]._tc.xpath(".//w:t/text()")) == condition["title"]
        assert row.cells[2].text == "、".join(condition["genes"])
        assert row.cells[3].text == condition["index_inheritance"]
        # Chinese and English share one heading and one 12-point bold run.
        heading = next(p for p in doc.paragraphs
                       if p._p.xpath("./w:bookmarkStart/@w:name") == [f"acmgsf_{condition['id']}"])
        assert heading.text == f"{number:02}　{condition['title']} {condition['english']}"
        assert all(run.font.size.pt == 12 and run.bold for run in heading.runs)
        assert not heading._p.xpath(".//w:br")
        assert condition["english"] not in paragraphs
    # Finished exports contain neither unresolved comments nor revision markup.
    assert not doc.element.xpath(".//w:ins|.//w:del|.//w:rPrChange|.//w:pPrChange|.//w:commentRangeStart")


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
        assert "此處列出清單中之84個基因" in text
        assert "36　遺傳性轉甲狀腺素蛋白類澱粉沉積症" in text
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
