import gzip
import json

from docx import Document

from app.adapters import pgx_tsv, snv_tsv
from app.services import docx_export


def test_acmg_ad_narrative_is_risk_not_diagnosis():
    lines = docx_export._health_acmg_narrative(
        "LDLR",
        disease_text="Familial hypercholesterolemia",
        inheritance_codes=["AD"],
        acmg="Pathogenic",
        zygosity="Heterozygous",
        same_gene_count=1,
        sex_karyotype="",
    )

    assert "具有較高的相關疾病風險" in lines[1]
    assert "遺傳諮詢或門診相關專科" in lines[2]
    assert "已罹患" not in "".join(lines)


def test_acmg_ar_single_variant_uses_carrier_wording():
    lines = docx_export._health_acmg_narrative(
        "ATP7B",
        disease_text="Wilson disease",
        inheritance_codes=["AR"],
        acmg="Likely pathogenic",
        zygosity="Heterozygous",
        same_gene_count=1,
        sex_karyotype="",
    )

    assert "本次僅檢出一個符合報告條件之變異" in lines[0]
    assert "符合報告條件的變異" not in lines[0]
    assert "檢測結果符合帶因者狀態" in lines[0]
    assert "帶因者" in lines[1]
    assert "多數不會出現典型疾病症狀" not in lines[1]
    assert "然而" not in lines[1]
    assert "故仍可能" in lines[1]
    assert "資料庫尚未收錄" in lines[1]


def test_structure_label_uses_intron_when_exon_is_placeholder():
    assert docx_export._structure_label({
        "exon": ".",
        "intron": "10/10",
    }) == "intron10"
    assert snv_tsv._clean_vep_rank(".") == ""


def test_acmg_ar_multiple_variants_does_not_assume_phase():
    lines = docx_export._health_acmg_narrative(
        "ATP7B",
        disease_text="Wilson disease",
        inheritance_codes=["AR"],
        acmg="Pathogenic",
        zygosity="Heterozygous",
        same_gene_count=2,
        sex_karyotype="",
    )

    assert "體染色體隱性遺傳" in lines[0]
    assert "此為致病性之變異位點。由於兩變異之相位尚未確認" in lines[1]
    assert "此為致病性之變異位點，由於兩變異之相位尚未確認" not in lines[1]
    assert "由於兩變異之相位尚未確認" in lines[1]
    assert "建議進行家族成員檢測，以釐清" in lines[1]
    assert "兩變異是否位於不同等位基因" in lines[1]
    assert "體染色體隱性疾病之雙等位基因致病型態" in lines[1]
    assert "本結果表示受檢者可能具有較高的相關疾病風險" not in lines[1]
    assert "相位分析" in lines[2]
    assert "遺傳諮詢或門診相關專科" in lines[2]


def test_acmg_x_linked_uses_karyotype():
    xy = docx_export._health_acmg_narrative(
        "ABCD1",
        disease_text="X-linked adrenoleukodystrophy",
        inheritance_codes=["XL"],
        acmg="Pathogenic",
        zygosity="Hemizygous",
        same_gene_count=1,
        sex_karyotype="XY",
    )
    xx = docx_export._health_acmg_narrative(
        "ABCD1",
        disease_text="X-linked adrenoleukodystrophy",
        inheritance_codes=["XL"],
        acmg="Pathogenic",
        zygosity="Heterozygous",
        same_gene_count=1,
        sex_karyotype="XX",
    )

    assert "性聯遺傳" in xy[0] and "可能具有較高" in xy[1]
    assert "依本次檢測所判定" not in xy[1]
    assert "性聯遺傳" in xx[0] and "異型合子女性" in xx[1]
    assert "可能無症狀" in xx[1] and "臨床表現" in xx[1]


def test_health_karyotype_prefers_ploidy_sidecar(tmp_path, monkeypatch):
    sample_dir = tmp_path / "S1"
    sample_dir.mkdir()
    with gzip.open(sample_dir / "ploidy.vcf.gz", "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n##estimatedSexKaryotype=XY\n#CHROM\tPOS\n")
    monkeypatch.setattr(docx_export, "TERTIARY_OUTPUT_ROOT", tmp_path)

    assert docx_export._health_sex_karyotype("S1", {"Sex": "F"}) == "XY"


def test_health_karyotype_does_not_read_uncopied_source_sibling(tmp_path, monkeypatch):
    sample_dir = tmp_path / "S1"
    source_dir = tmp_path / "source"
    sample_dir.mkdir()
    source_dir.mkdir()
    source_vcf = source_dir / "VAL-3.hard-filtered.vcf.gz"
    source_vcf.touch()
    with gzip.open(source_dir / "VAL-3.ploidy.vcf.gz", "wt", encoding="utf-8") as handle:
        handle.write("##estimatedSexKaryotype=XX\n#CHROM\tPOS\n")
    (sample_dir / "pipeline_source.json").write_text(
        json.dumps({"source_vcf_path": str(source_vcf)}),
        encoding="utf-8",
    )
    monkeypatch.setattr(docx_export, "TERTIARY_OUTPUT_ROOT", tmp_path)

    assert docx_export._health_sex_karyotype("S1", {"Sex": "M"}) == "XY"


def test_health_karyotype_preserves_nonstandard_call(tmp_path, monkeypatch):
    sample_dir = tmp_path / "S1"
    sample_dir.mkdir()
    with gzip.open(sample_dir / "ploidy.vcf.gz", "wt", encoding="utf-8") as handle:
        handle.write("##estimatedSexKaryotype=XXY\n#CHROM\tPOS\n")
    monkeypatch.setattr(docx_export, "TERTIARY_OUTPUT_ROOT", tmp_path)

    assert docx_export._health_sex_karyotype("S1", {"Sex": "M"}) == "XXY"


def test_pgx_summary_prefers_cpic_for_same_gene_and_drug():
    pgx = {
        "guideline_annotations": [
            {
                "section": "CPIC Guideline Annotation",
                "classification": "Strong",
                "alternate_drug_available": True,
                "genes": ["CYP2C19"],
                "drug": "clopidogrel",
                "recommendation": "Use prasugrel or ticagrelor at standard dose.",
            },
            {
                "section": "FDA PGx Association",
                "fda_category": "therapeutic_management",
                "genes": ["CYP2C19"],
                "drug": "clopidogrel",
                "recommendation": "Consider an alternative.",
            },
        ],
    }

    alerts = docx_export._pgx_summary_alerts(pgx)

    assert [row["source"] for row in alerts] == ["CPIC"]
    assert alerts[0]["recommendation"] == "Use prasugrel or ticagrelor at standard dose."


def test_pgx_report_genes_are_fixed_cpic_level_a_scope():
    genes = docx_export._pgx_report_genes({})

    assert len(genes) == 22
    assert "VKORC1" in genes
    assert "NAT2" in genes
    assert "CYP3A4" not in genes


def test_pgx_full_groups_include_cpic_and_fda_for_json_only_gene():
    pgx = {
        "genes": {
            "VKORC1": {
                "additional": True,
                "details": {
                    "label": "rs9923231 variant (T)/rs9923231 variant (T)",
                    "phenotypes": ["-1639 AA"],
                    "allele1_function": "Higher coumarin sensitivity",
                    "allele2_function": "Higher coumarin sensitivity",
                },
            },
        },
        "guideline_annotations": [
            {
                "section": "CPIC Guideline Annotation",
                "classification": "Strong",
                "dosing_information": True,
                "genes": ["VKORC1"],
                "drug": "warfarin",
                "recommendation": "Use genotype-guided dosing.",
            },
            {
                "section": "FDA PGx Association",
                "fda_category": "therapeutic_management",
                "classification": "Unspecified",
                "genes": ["VKORC1"],
                "drug": "warfarin",
                "recommendation": "Select initial dosage using clinical and genetic factors.",
            },
        ],
    }

    assert docx_export._pgx_full_groups(pgx) == [{
        "gene": "VKORC1",
        "diplotype": "rs9923231 T/T（-1639 A/A）",
        "phenotype": "Higher coumarin sensitivity",
        "recommendations": [
            {
                "drug": "warfarin",
                "source": "CPIC",
                "recommendation": "Use genotype-guided dosing.",
                "level": "Strong",
            },
            {
                "drug": "warfarin",
                "source": "FDA PGx Association",
                "recommendation": "Select initial dosage using clinical and genetic factors.",
                "level": "Therapeutic Management",
            },
        ],
    }]


def test_pgx_summary_excludes_standard_cpic_action():
    pgx = {
        "guideline_annotations": [{
            "section": "CPIC Guideline Annotation",
            "classification": "Strong",
            "dosing_information": False,
            "alternate_drug_available": False,
            "other_prescribing_guidance": False,
            "genes": ["DPYD"],
            "drug": "fluorouracil",
            "recommendation": "Use label-recommended dosage and administration.",
        }],
    }

    assert docx_export._pgx_summary_alerts(pgx) == []


def test_pgx_not_recommended_is_actionable_alternative():
    assert docx_export._pgx_action_zh("Ivacaftor is not recommended") == "考慮替代藥物"


def test_pgx_full_groups_exclude_standard_and_uncertain_results():
    pgx = {
        "genes": {
            "DPYD": {"details": {"label": "Reference/Reference", "phenotypes": ["Normal Metabolizer"]}},
            "CACNA1S": {"details": {"label": "Reference/Reference", "phenotypes": ["Uncertain Susceptibility"]}},
        },
        "guideline_annotations": [
            {
                "section": "CPIC Guideline Annotation",
                "classification": "Strong",
                "genes": ["DPYD"],
                "drug": "fluorouracil",
                "recommendation": "Use label-recommended dosage and administration.",
            },
            {
                "section": "CPIC Guideline Annotation",
                "classification": "Strong",
                "other_prescribing_guidance": True,
                "genes": ["CACNA1S"],
                "drug": "sevoflurane",
                "recommendation": "Clinical findings should guide use.",
            },
        ],
    }

    assert docx_export._pgx_full_groups(pgx) == []


def test_pgx_full_groups_deduplicate_same_recommendation_at_strongest_level():
    base = {
        "section": "CPIC Guideline Annotation",
        "alternate_drug_available": True,
        "genes": ["CYP2C19"],
        "drug": "clopidogrel",
        "recommendation": "Avoid clopidogrel if possible.",
    }
    pgx = {
        "genes": {
            "CYP2C19": {"details": {"label": "*2/*2", "phenotypes": ["Poor Metabolizer"]}},
        },
        "guideline_annotations": [
            {**base, "classification": "Moderate"},
            {**base, "classification": "Strong"},
        ],
    }

    groups = docx_export._pgx_full_groups(pgx)

    assert len(groups[0]["recommendations"]) == 1
    assert groups[0]["recommendations"][0]["level"] == "Strong"


def test_pgx_genotype_rows_use_fixed_cpic_level_a_scope():
    pgx = {
        "genes": {
            "CYP2C19": {"details": {"label": "*2/*2", "phenotypes": ["Poor Metabolizer"]}},
            "CYP3A4": {"details": {"label": "*1/*1", "phenotypes": ["Normal Metabolizer"]}},
        },
    }

    rows = docx_export._pgx_genotype_rows(pgx)

    assert len(rows) == 22
    assert ["CYP2C19", "*2/*2", "Poor Metabolizer"] in rows
    assert not any(row[0] == "CYP3A4" for row in rows)
    assert any(row[0] == "ABCG2" and row[1:] == ["—", "No phenotype assigned"] for row in rows)


def test_pgx_full_recommendation_table_width_matches_genotype_table():
    columns = docx_export._PGX_FULL_RECOMMENDATION_COLUMNS

    assert columns == [
        ("藥物", 17, "buffered"),
        ("基因與表型", 18, "word-buffered"),
        ("CPIC/FDA 建議", 49, "word-buffered"),
    ]
    assert sum(col[1] for col in columns) + 2 == 86


def test_pgx_gene_result_falls_back_to_tsv_for_mt_rnr1():
    pgx = {
        "genes": {
            "MT-RNR1": {
                "diplotype": "m.1555A>G positive",
                "phenotype": "HIGH risk",
                "details": {
                    "label": "Unknown",
                    "phenotypes": ["No Result"],
                },
            },
        },
    }

    assert docx_export._pgx_gene_result(pgx, "MT-RNR1") == ("m.1555A>G positive", "HIGH risk")


def test_pgx_genotype_rows_treat_dash_phenotype_as_not_assigned():
    pgx = {
        "genes": {
            "VKORC1": {
                "details": {
                    "label": "rs9923231 T/T",
                    "phenotypes": ["—"],
                },
            },
        },
    }

    rows = docx_export._pgx_genotype_rows(pgx)

    assert ["VKORC1", "rs9923231 T/T", "No phenotype assigned"] in rows


def test_pgx_gene_result_falls_back_to_allele_function_for_ifnl3():
    pgx = {
        "genes": {
            "IFNL3": {
                "details": {
                    "label": "rs12979860 reference (C)/rs12979860 reference (C)",
                    "phenotypes": ["n/a"],
                    "allele1_function": "Favorable response allele",
                    "allele2_function": "Favorable response allele",
                },
            },
        },
    }

    assert docx_export._pgx_gene_result(pgx, "IFNL3") == (
        "rs12979860 reference (C)/rs12979860 reference (C)",
        "Favorable response allele",
    )


def test_pgx_gene_result_keeps_vkorc1_sensitivity_as_phenotype():
    pgx = {
        "genes": {
            "VKORC1": {
                "details": {
                    "label": "rs9923231 variant (T)/rs9923231 variant (T)",
                    "phenotypes": ["-1639 AA"],
                    "allele1_function": "Higher coumarin sensitivity",
                    "allele2_function": "Higher coumarin sensitivity",
                },
            },
        },
    }

    assert docx_export._pgx_gene_result(pgx, "VKORC1") == (
        "rs9923231 T/T（-1639 A/A）",
        "Higher coumarin sensitivity",
    )


def test_pgx_drug_groups_are_drug_centric_and_keep_strongest_cpic():
    pgx = {
        "genes": {
            "CYP2C19": {"details": {"label": "*2/*2", "phenotypes": ["Poor Metabolizer"]}},
        },
        "guideline_annotations": [
            {
                "section": "CPIC Guideline Annotation",
                "classification": "Moderate",
                "alternate_drug_available": True,
                "genes": ["CYP2C19"],
                "drug": "clopidogrel",
                "recommendation": "Avoid clopidogrel if possible. Consider an alternative P2Y12 inhibitor.",
            },
            {
                "section": "CPIC Guideline Annotation",
                "classification": "Strong",
                "alternate_drug_available": True,
                "genes": ["CYP2C19"],
                "drug": "clopidogrel",
                "recommendation": "Avoid clopidogrel if possible. Use prasugrel or ticagrelor.",
            },
            {
                "section": "FDA Label Annotation",
                "alternate_drug_available": True,
                "genes": ["CYP2C19"],
                "drug": "clopidogrel",
                "recommendation": "Consider use of another platelet P2Y12 inhibitor.",
            },
            {
                "section": "FDA PGx Association",
                "fda_category": "therapeutic_management",
                "genes": ["CYP2C19"],
                "drug": "clopidogrel",
                "recommendation": "Consider use of another platelet P2Y12 inhibitor.",
            },
        ],
    }

    groups = docx_export._pgx_drug_groups(pgx)

    assert len(groups) == 1
    assert groups[0]["drug"] == "Clopidogrel"
    assert groups[0]["action"] == "考慮替代藥物"
    assert groups[0]["source_level"] == "CPIC Strong"
    assert docx_export._pgx_gene_phenotype_text(groups[0]) == "CYP2C19 Poor Metabolizer"
    assert [
        (rec["source"], rec["level"])
        for rec in groups[0]["recommendations"]
    ] == [
        ("CPIC", "Strong"),
        ("FDA PGx Association", "Therapeutic Management"),
        ("FDA Label", "Unspecified"),
    ]


def test_pgx_recommendation_text_combines_fda_sources_after_content():
    text = docx_export._pgx_recommendation_text({
        "recommendations": [
            {
                "source": "CPIC",
                "level": "Strong",
                "recommendation": '"Use alternative therapy."',
            },
            {
                "source": "FDA Label",
                "level": "Unspecified",
                "recommendation": '"Use label dosing."',
            },
            {
                "source": "FDA PGx Association",
                "level": "Therapeutic Management",
                "recommendation": '"Adjust dose."',
            },
        ],
    })

    assert text == (
        "Use alternative therapy. (CPIC Strong)；"
        "Adjust dose. Use label dosing. (FDA Therapeutic Management / FDA Label)"
    )
    assert '"' not in text

    cleaned = docx_export._pgx_clean_recommendation_text(
        '&quot;[...]increase monitoring for adverse reactions.&quot; '
        'In CYP2C19 poor metabolizers ... the starting dose should be reduced. '
        'No dosage adjustment is needed...For known poor metabolizers, '
        '...thioridazine is contraindicated...in patients.'
    )
    assert cleaned == (
        "increase monitoring for adverse reactions. "
        "In CYP2C19 poor metabolizers the starting dose should be reduced. "
        "No dosage adjustment is needed. For known poor metabolizers, "
        "thioridazine is contraindicated in patients."
    )
    assert "..." not in cleaned and "[...]" not in cleaned


def test_pgx_summary_rows_exclude_fda_label_only_drugs():
    rows = docx_export._pgx_summary_rows_from_drug_groups([
        {
            "drug": "Belzutifan",
            "genes": {"CYP2C19": {}},
            "action": "調整劑量並監測",
            "source_level": "FDA Label",
            "recommendations": [
                {
                    "gene": "CYP2C19",
                    "source": "FDA Label",
                    "level": "Unspecified",
                    "action": "調整劑量並監測",
                    "source_priority": 3,
                },
            ],
        },
        {
            "drug": "Deutetrabenazine",
            "genes": {"CYP2D6": {}},
            "action": "調整劑量並監測",
            "source_level": "FDA Therapeutic Management",
            "recommendations": [
                {
                    "gene": "CYP2D6",
                    "source": "FDA PGx Association",
                    "level": "Therapeutic Management",
                    "action": "調整劑量並監測",
                    "source_priority": 2,
                },
            ],
        },
    ])

    assert rows == [{
        "drug": "Deutetrabenazine",
        "gene": "CYP2D6",
        "action": "調整劑量並監測",
        "source_level": "FDA Therapeutic Management",
    }]


def test_pgx_adapter_keeps_json_when_tsv_is_absent(tmp_path, monkeypatch):
    monkeypatch.setattr(pgx_tsv, "_compact_report_json", lambda _path: {
        "pharmcat_version": "3.2.0",
        "data_version": "2026-02-09",
        "timestamp": "2026-06-30",
        "messages": [],
        "genes": {
            "VKORC1": {
                "label": "rs9923231 variant (T)/rs9923231 variant (T)",
                "phenotypes": ["-1639 AA"],
            },
        },
        "guideline_annotations": [{
            "section": "FDA PGx Association",
            "genes": ["VKORC1"],
        }],
    })

    result = pgx_tsv.load_pgx(
        tmp_path / "missing.pgx.tsv",
        tmp_path / "pharmcat.report.json",
    )

    assert result["pharmcat_version"] == "3.2.0"
    assert result["additional_genes"] == ["VKORC1"]
    assert result["genes"]["VKORC1"]["details"]["phenotypes"] == ["-1639 AA"]


def test_pgx_only_annotation_lists_all_cpic_level_a_genes():
    doc = Document()

    docx_export._section_health_annotations(doc, {"pgx"}, {})

    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    assert "CPIC) Level A gene" in text
    assert "ABCG2" in text and "VKORC1" in text
    assert "CYP3A4" not in text
    assert "FDA" not in text


def test_health_pathogenicities_translate_and_combine():
    assert docx_export._health_pathogenicities_zh([
        "Pathogenic",
        "Uncertain significance",
    ]) == "致病性及不確定意義"


def test_vkorc1_display_moves_position_phenotype_into_genotype():
    assert docx_export._pgx_display_genotype(
        "VKORC1",
        "rs9923231 variant (T)/rs9923231 variant (T)",
        "-1639 AA",
    ) == ("rs9923231 T/T（-1639 A/A）", "—")


def test_pgx_summary_rows_translate_and_sort_by_drug():
    alerts = [
        {
            "gene": "G6PD",
            "drug": "rasburicase",
            "source": "CPIC",
            "level": "Strong",
            "recommendation": "Avoid rasburicase.",
            "priority": True,
        },
        {
            "gene": "G6PD",
            "drug": "dapsone",
            "source": "CPIC",
            "level": "Strong",
            "recommendation": "Avoid dapsone.",
            "priority": True,
        },
    ]

    assert docx_export._pgx_summary_rows(alerts) == [
        {
            "drug": "Dapsone",
            "gene": "G6PD",
            "action": "考慮替代藥物",
            "source_level": "CPIC Strong",
        },
        {
            "drug": "Rasburicase",
            "gene": "G6PD",
            "action": "考慮替代藥物",
            "source_level": "CPIC Strong",
        },
    ]


def test_health_bundle_name_follows_selected_sections():
    assert docx_export._health_test_bundle_name({"acmg_sf", "pgx"}) == (
        "ACMG疾病風險基因及藥物基因體學基因篩檢"
    )
    assert docx_export._health_test_bundle_name({"pgx"}) == "藥物基因體學基因篩檢"
    assert docx_export._health_test_bundle_name({"acmg_sf"}) == "ACMG疾病風險基因篩檢"


def test_health_acmg_section_titles_use_current_wording():
    assert dict(docx_export._HEALTH_SECTION_ORDER)["acmg_sf"] == "第一類：ACMG疾病風險基因"
    assert docx_export._HEALTH_ACMG_GENE_LIST_TITLE.startswith("第一類：ACMG疾病風險基因")
    assert "可採取醫療處置之疾病風險基因" not in dict(docx_export._HEALTH_SECTION_ORDER)["acmg_sf"]


def test_health_header_and_acmg_caution_follow_current_template():
    doc = Document()

    docx_export._section_health_patient_header(doc, {"LIS_ID": "SF1-dragen"})

    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    assert "<<基因醫學部基因檢測檢驗分析研究報告>>" in text
    assert "研究報告）" not in text
    assert "ACMG SF 3.3" in docx_export._HEALTH_ACMG_CAUTION
    assert "分析與 ACMG SF 3.3 所列遺傳性疾病相關的風險基因" in docx_export._HEALTH_ACMG_CAUTION
    assert "可採取醫療處置之遺傳性疾病" not in docx_export._HEALTH_ACMG_CAUTION


def test_health_pgx_action_categories_remove_below_from_residual_drug_note():
    doc = Document()
    pgx = {
        "genes": {
            "CYP2C19": {"details": {"label": "*2/*2", "phenotypes": ["Poor Metabolizer"]}},
        },
        "guideline_annotations": [{
            "section": "CPIC Guideline Annotation",
            "classification": "Strong",
            "alternate_drug_available": True,
            "genes": ["CYP2C19"],
            "drug": "clopidogrel",
            "recommendation": "Use an alternative antiplatelet therapy.",
        }],
    }

    docx_export._render_health_pgx_section(doc, "藥物基因體學", pgx)

    text = "\n".join(paragraph.text for paragraph in doc.paragraphs)
    assert "其餘未列於下方之藥物" not in text
    assert "其餘未列之藥物" in text
    assert "完整用藥建議" in text
