from docx import Document

from app.services import docx_export


def _gjb2_variant(
    variant_id: str,
    *,
    rs_id: str,
    hgvs_c: str,
    hgvs_p: str,
    acmg: str,
) -> dict:
    return {
        "id": variant_id,
        "gene_symbol": "GJB2",
        "ensembl_transcript": "ENST00000382848",
        "refseq_transcript": "NM_004004.6",
        "rs_id": rs_id,
        "exon": "2/2",
        "HGVS_C": hgvs_c,
        "HGVS_P": hgvs_p,
        "zygosity": "Heterozygous",
        "CLNSIG": "Pathogenic",
        "ACMG_classification": acmg,
        "Disease1": "Deafness, autosomal recessive 1A (220290)(AR)",
    }


def test_diagnosis_groups_same_gene_snvs_and_combines_acmg_wording():
    doc = Document()
    first = _gjb2_variant(
        "v1",
        rs_id="rs80338943",
        hgvs_c="c.235del",
        hgvs_p="p.Leu79CysfsTer3",
        acmg="Pathogenic",
    )
    second = _gjb2_variant(
        "v2",
        rs_id="rs72474224",
        hgvs_c="c.109G>A",
        hgvs_p="p.Val37Ile",
        acmg="Likely pathogenic",
    )

    docx_export._section_results(
        doc,
        {"variants": {"v1": first, "v2": second}},
        {"status": {"v1": "1", "v2": "1"}, "edits": {}},
        "WES",
    )

    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    assert text.count("GJB2 (ENST00000382848; NM_004004.6)") == 1
    assert text.count("ACMG&AMP指引") == 1
    assert "c.235del" in text
    assert "c.109G>A" in text
    assert text.count("GJB2為Deafness, autosomal recessive 1A的致病基因之一") == 1
    assert "此為致病性及疑似致病性之變異位點，與臨床症狀相關。" in text

    paragraphs = [paragraph.text for paragraph in doc.paragraphs]
    patho_index = paragraphs.index(
        "    2. 此為致病性及疑似致病性之變異位點，與臨床症狀相關。"
    )
    assert paragraphs[patho_index + 1:patho_index + 5] == [
        "",
        "    第二類：其他變異位點",
        "    未找到其他變異位點。",
        "",
    ]
    assert paragraphs[patho_index + 5].startswith("    建議比對臨床表徵")


def test_ploidy_finding_report_keeps_standard_headings_and_patient_wording():
    doc = Document()
    finding = {
        "id": "PLOIDY-chr21-GAIN-test",
        "CHROM": "chr21",
        "dosage_call": "gain",
        "interpretation": "possible trisomy 21",
        "NDC": 1.346,
        "filter": "SUSPECT",
        "pipeline_source": "NCKUH_PLOIDY_MOSDEPTH",
    }
    docx_export._section_results(
        doc,
        {"ploidy_findings": {finding["id"]: finding}},
        {"status": {finding["id"]: "1"}, "edits": {finding["id"]: {
            "disease": "唐氏症", "ACMG_class_sv": "5",
        }}},
        "WES",
    )
    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    assert "第一類：與臨床症狀相關基因之已知致病性變異位點" in text
    assert "第二類：其他變異位點" in text
    assert "與臨床症狀相關之變異或染色體劑量訊號" not in text
    assert "[GRCh38] chr21 trisomy" in text
    assert "拷貝數" in text and "3（疑似）" not in text
    assert "     1    21" in text and " 3 " in text
    assert "    1. 此為第 21 號染色體三體，與唐氏症相關。" in text
    assert "ACMG 分類" not in text
    assert "Ploidy VCF" not in text
    assert "NDC" not in text and "FILTER" not in text
    assert "參考資料:" not in text
    assert "chr21:1-46709983" not in text


def test_second_class_snv_and_cnv_descriptions_do_not_claim_full_clinical_match():
    doc = Document()
    snv = _gjb2_variant(
        "snv1", rs_id="rs80338943", hgvs_c="c.235del",
        hgvs_p="p.Leu79CysfsTer3", acmg="Pathogenic",
    )
    cnv = {
        "id": "cnv1", "source": "cnv", "CHROM": "17", "POS": 100,
        "END": 200, "sv_type": "DEL", "copy_number": 1,
        "zygosity": "het", "acmg_class": 5, "genes": [],
    }
    docx_export._section_results(
        doc,
        {"variants": {"snv1": snv}, "cnv_variants": {"cnv1": cnv}},
        {"status": {"snv1": "2", "cnv1": "2"}, "edits": {}},
        "WGS",
    )
    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    tail = "無法完全解釋受檢者全部之臨床症狀，其臨床意義須由醫師配合其他相關資料進行最佳綜合判斷。"
    assert f"此為致病性之變異位點，{tail}" in text
    assert f"此為致病性之變異，{tail}" in text
    assert "此為致病性之變異位點，與臨床症狀相關。" not in text
    assert "此為致病性之變異，與臨床症狀相關。" not in text


def test_ploidy_loss_without_manual_edits_stays_unclassified():
    doc = Document()
    docx_export._ploidy_variant_block(
        doc,
        {"CHROM": "chr18", "dosage_call": "loss", "NDC": 0.7, "filter": "SUSPECT"},
        tier="1", edits={},
    )
    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    assert "[GRCh38] chr18 monosomy" in text
    assert "1（疑似）" not in text
    assert "    1. 此為第 18 號染色體單體。" in text
    assert "ACMG 分類" not in text
    assert "Ploidy VCF" not in text


def test_second_class_ploidy_keeps_original_second_class_heading():
    doc = Document()
    finding = {"id": "PLOIDY-chr21-GAIN-test", "CHROM": "chr21", "dosage_call": "gain"}
    docx_export._section_results(
        doc,
        {"ploidy_findings": {finding["id"]: finding}},
        {"status": {finding["id"]: "2"}, "edits": {finding["id"]: {"disease": "唐氏症"}}},
        "WES",
    )
    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    assert "第二類：其他變異位點" in text
    assert "第二類：其他變異或染色體劑量訊號" not in text
    assert "    1. 此為第 21 號染色體三體，與唐氏症相關。" in text


def test_diagnosis_joins_all_checked_snv_diseases_and_mim_numbers():
    doc = Document()
    variant = _gjb2_variant(
        "v1",
        rs_id="rs80338943",
        hgvs_c="c.235del",
        hgvs_p="p.Leu79CysfsTer3",
        acmg="Pathogenic",
    )
    variant["Disease1"] = "Disease A (600001)(AD)"
    variant["Disease2"] = "Disease B (600001)(AR)"
    variant["Disease3"] = "Disease C (600003)(XLR)"
    edits = {"report_diseases": {"1": True, "2": True}}

    assert docx_export._picked_diseases_for_snv(variant, edits) == [
        "Disease A (600001)(AD)",
        "Disease B (600001)(AR)",
    ]

    docx_export._section_results(
        doc,
        {"variants": {"v1": variant}},
        {"status": {"v1": "1"}, "edits": {"v1": edits}},
        "WES",
    )

    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    assert (
        "GJB2為Disease A、Disease B的致病基因之一，"
        "其遺傳模式屬於體染色體顯性遺傳、體染色體隱性遺傳 "
        "(Phenotype MIM number: 600001、600001)。"
    ) in text
    assert "Disease C" not in text


def test_diagnosis_without_snv_disease_ticks_keeps_first_slot_fallback():
    variant = _gjb2_variant(
        "v1", rs_id="", hgvs_c="c.1A>G", hgvs_p="p.Met1Val", acmg="VUS"
    )
    variant["Disease1"] = "Disease A (600001)(AD)"
    variant["Disease2"] = "Disease B (600002)(AR)"

    assert docx_export._picked_diseases_for_snv(variant, {}) == [
        "Disease A (600001)(AD)"
    ]
    assert "Disease B" not in docx_export._omim_block_for_snv(variant, {})


def test_diagnosis_without_first_or_second_category_uses_compact_negative_summary():
    doc = Document()

    docx_export._section_results(
        doc,
        {
            "patient_phenotype": [
                {"phenotype": "HP:0001263", "label": "Global developmental delay"},
                {"phenotype": "HP:0002119", "label": "Ventriculomegaly"},
            ],
            "selected_panels": [{
                "name": "WES-I__兒科__先天神經肌肉疾病",
            }],
        },
        {"status": {}, "edits": {}},
        "WES",
    )

    paragraphs = [paragraph.text for paragraph in doc.paragraphs]
    assert paragraphs == [
        "三、檢測結果",
        "  檢體說明:",
        "    檢體類別：血液",
        "  綜合說明:",
        "    在非特定 (Global developmental delay, Ventriculomegaly, 先天神經肌肉疾病) 檢驗套組中未找到已知致病性位點。",
        "    建議持續追蹤。",
        "  參考資料:",
        "    依據疾病資料庫中目前記載，本次檢測套組所涵蓋的基因，未檢測到具有足夠疾病關連性的致病變異。",
        "    此報告僅供參考，臨床判斷仍應以病患的實際狀況為主。",
        "",
    ]


def test_em_dash_uses_full_width_for_ascii_table_padding():
    assert docx_export._str_width("—") == 2
    assert docx_export._pad_right("—", 13) == "—" + (" " * 11)


def test_cnv_report_disease_is_selected_union_plus_free_text():
    edits = {
        "report_disease_items": {
            "omim:GENE:1:1": {
                "label": "Disease A", "source": "omim", "gene": "GENE",
                "phenotype_mim": "600001", "inheritance": "AD",
            },
            "overlap:p_loss:Disease B": {
                "label": "Disease B", "source": "overlap", "overlap_type": "p_loss",
            },
            "duplicate": {"label": "disease a", "source": "overlap"},
        },
        "disease": "Disease A、Manual disease",
    }

    assert docx_export._cnv_report_disease(edits) == (
        "Disease A、Disease B、Manual disease"
    )


def test_single_gene_cnv_keeps_one_point_and_uses_combined_disease_text():
    doc = Document()
    variant = {
        "id": "cnv1", "source": "cnv", "CHROM": "17", "POS": 100,
        "END": 200, "sv_type": "DEL", "copy_number": 1, "zygosity": "het",
        "acmg_class": 5,
        "genes": [{
            "gene": "COL1A1", "location": "exon1-exon2", "omim_id": "120150",
            "omim_phenotype": "Legacy disease (600000)(AD)",
            "omim_inheritance": "AD",
        }],
    }
    edits = {
        "report_disease_items": {
            "omim:COL1A1:120150:1": {"label": "Disease A", "source": "omim"},
            "overlap:p_loss:Disease B": {"label": "Disease B", "source": "overlap"},
        },
        "disease": "Manual disease",
    }

    docx_export._cnv_variant_block(
        doc, variant, tier="1", is_wgs=True, edits=edits
    )

    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    assert "COL1A1為Disease A、Disease B、Manual disease的致病基因之一" in text
    assert "Phenotype MIM number" not in text


def test_wes_single_gene_cnv_uses_exon_range_and_cnv_wording():
    doc = Document()
    variant = {
        "id": "cnv1", "source": "cnv", "CHROM": "17", "POS": 50180000,
        "END": 50200000, "sv_type": "DEL", "copy_number": 1,
        "zygosity": "het", "acmg_class": 5,
        "genes": [{
            "gene": "COL1A1", "location": "txStart-intron28",
            "exon_count": 51, "omim_id": "120150",
            "omim_phenotype": "Osteogenesis imperfecta (166200)(AD)",
            "omim_inheritance": "AD",
        }],
    }

    docx_export._cnv_variant_block(
        doc, variant, tier="1", is_wgs=False, edits={}
    )

    paragraphs = [paragraph.text for paragraph in doc.paragraphs]
    assert (
        "    1. 此片段位於第 17 號染色體上 COL1A1 基因，"
        "涵蓋 Exon 1 至 Exon 28 區域。"
    ) in paragraphs
    assert "    3. 此為致病性之變異，與臨床症狀相關。" in paragraphs
    assert any("內含子(Intron)，則無法" in paragraph for paragraph in paragraphs)
    assert all("內含子(Intron) ，" not in paragraph for paragraph in paragraphs)

    reference = docx_export._cnv_reference_text(
        variant, {}, docx_export._omim_genes(variant), "缺失", False
    )
    assert "COL1A1 基因之 Exon 1 至 Exon 28 區域" in reference
    assert "評測此變異為「Pathogenic」" in reference
    assert "評測此變異位點" not in reference


def test_wes_single_gene_cnv_marks_breakpoint_exon_as_partial():
    assert docx_export._wes_exon_span({
        "location": "txStart-exon4", "exon_count": 10,
    }) == "Exon 1 至部分 Exon 4 區域"
    assert docx_export._wes_exon_span({
        "location": "intron2-intron4", "exon_count": 10,
    }) == "Exon 3 至 Exon 4 區域"
    assert docx_export._wes_exon_span({
        "location": "exon4-exon4", "exon_count": 10,
    }) == "部分 Exon 4 區域"


def test_wgs_single_gene_cnv_retains_location_wording():
    doc = Document()
    variant = {
        "id": "cnv1", "source": "cnv", "CHROM": "17", "POS": 100,
        "END": 200, "sv_type": "DEL", "copy_number": 1,
        "zygosity": "het", "acmg_class": 5,
        "genes": [{
            "gene": "COL1A1", "location": "txStart-intron28",
            "exon_count": 51, "omim_id": "120150",
            "omim_phenotype": "Osteogenesis imperfecta (166200)(AD)",
            "omim_inheritance": "AD",
        }],
    }

    docx_export._cnv_variant_block(
        doc, variant, tier="1", is_wgs=True, edits={}
    )
    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    assert "COL1A1 基因之基因起始至 Intron 28 區域" in text
    assert "涵蓋 Exon 1 至 Exon 28 區域" not in text


def test_wes_multi_gene_cnv_retains_location_wording_with_one_omim_gene():
    doc = Document()
    variant = {
        "id": "cnv1", "source": "cnv", "CHROM": "17", "POS": 100,
        "END": 200, "sv_type": "DEL", "copy_number": 1,
        "zygosity": "het", "acmg_class": 5, "gene_count": 2,
        "genes": [{
            "gene": "COL1A1", "location": "txStart-intron28",
            "exon_count": 51, "omim_id": "120150",
            "omim_phenotype": "Osteogenesis imperfecta (166200)(AD)",
            "omim_inheritance": "AD",
        }, {
            "gene": "OTHER", "location": "txStart-txEnd", "omim_id": "",
        }],
    }

    docx_export._cnv_variant_block(
        doc, variant, tier="1", is_wgs=False, edits={}
    )
    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    assert "COL1A1 基因之基因起始至 Intron 28 區域" in text
    assert "涵蓋 Exon 1 至 Exon 28 區域" not in text


def test_wes_single_gene_sv_retains_existing_location_and_site_wording():
    doc = Document()
    variant = {
        "id": "sv1", "source": "sv", "CHROM": "17", "POS": 100,
        "END": 200, "sv_type": "DEL", "copy_number": 1,
        "zygosity": "het", "acmg_class": 5,
        "genes": [{
            "gene": "COL1A1", "location": "txStart-intron28",
            "exon_count": 51, "omim_id": "120150",
            "omim_phenotype": "Osteogenesis imperfecta (166200)(AD)",
            "omim_inheritance": "AD",
        }],
    }

    docx_export._cnv_variant_block(
        doc, variant, tier="1", is_wgs=False, edits={}
    )
    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    assert "COL1A1 基因之基因起始至 Intron 28 區域" in text
    assert "此為致病性之變異位點，與臨床症狀相關。" in text
    reference = docx_export._cnv_reference_text(
        variant, {}, docx_export._omim_genes(variant), "缺失", False
    )
    assert "評測此變異位點為「Pathogenic」" in reference


def test_single_overlap_disease_does_not_inherit_gene_omim_metadata():
    doc = Document()
    variant = {
        "id": "cnv1", "source": "cnv", "CHROM": "17", "POS": 100,
        "END": 200, "sv_type": "DEL", "copy_number": 1, "zygosity": "het",
        "acmg_class": 5,
        "genes": [{
            "gene": "COL1A1", "location": "exon1-exon2", "omim_id": "120150",
            "omim_phenotype": "Legacy disease (600000)(AD)",
            "omim_inheritance": "AD",
        }],
    }
    edits = {
        "report_disease_items": {
            "overlap:p_loss:Overlap disease": {
                "label": "Overlap disease", "source": "overlap",
            },
        },
    }

    docx_export._cnv_variant_block(
        doc, variant, tier="1", is_wgs=True, edits=edits
    )

    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    assert "COL1A1為Overlap disease的致病基因之一" in text
    assert "Phenotype MIM number" not in text
    assert "顯性遺傳" not in text
