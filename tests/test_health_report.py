import gzip
import json

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

    assert "帶因者" in lines[1]
    assert "資料庫尚未收錄" in lines[1]


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
    assert "相位分析" in lines[1]
    assert "相位分析" in lines[2]


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


def test_pgx_summary_prefers_cpic_and_marks_priority():
    pgx = {
        "guideline_annotations": [{
            "fda_category": "therapeutic_management",
            "genes": ["CYP2C19"],
            "drug": "clopidogrel",
            "recommendation": "Consider an alternative.",
        }],
    }
    groups = [{
        "gene": "CYP2C19",
        "recommendations": [{
            "drug": "clopidogrel",
            "level": "Strong",
            "recommendation": "Use an alternative antiplatelet agent.",
        }],
    }]

    alerts = docx_export._pgx_summary_alerts(pgx, groups)

    assert len(alerts) == 1
    assert alerts[0]["source"] == "CPIC"
    assert alerts[0]["priority"] is True


def test_pgx_report_genes_union_tsv_calls_and_fda_annotations():
    pgx = {
        "gene_order": ["CYP2C19", "CACNA1S"],
        "genes": {
            "CYP2C19": {"additional": False},
            "CACNA1S": {"additional": False},
            "VKORC1": {"additional": True},
        },
        "guideline_annotations": [
            {"section": "FDA PGx Association", "genes": ["CYP2C19", "VKORC1"]},
            {"section": "DPWG Guideline Annotation", "genes": ["DPYD"]},
        ],
    }

    assert docx_export._pgx_report_genes(pgx) == ["CACNA1S", "CYP2C19", "VKORC1"]


def test_pgx_fda_groups_include_json_only_gene():
    pgx = {
        "genes": {
            "VKORC1": {
                "additional": True,
                "details": {
                    "label": "rs9923231 variant (T)/rs9923231 variant (T)",
                    "phenotypes": ["-1639 AA"],
                },
            },
        },
        "guideline_annotations": [{
            "section": "FDA PGx Association",
            "fda_category": "therapeutic_management",
            "genes": ["VKORC1"],
            "drug": "warfarin",
            "recommendation": "Select initial dosage using clinical and genetic factors.",
        }],
    }

    assert docx_export._pgx_fda_groups(pgx) == [{
        "gene": "VKORC1",
        "diplotype": "rs9923231 variant (T)/rs9923231 variant (T)",
        "phenotype": "-1639 AA",
        "recommendations": [{
            "drug": "warfarin",
            "recommendation": "Select initial dosage using clinical and genetic factors.",
            "level": "Therapeutic management",
        }],
    }]


def test_pgx_summary_rows_translate_and_group_recommendations():
    alerts = [
        {
            "gene": "G6PD",
            "drug": "dapsone",
            "source": "CPIC",
            "level": "Strong",
            "recommendation": "Avoid dapsone.",
            "priority": True,
        },
        {
            "gene": "G6PD",
            "drug": "rasburicase",
            "source": "CPIC",
            "level": "Strong",
            "recommendation": "Avoid rasburicase.",
            "priority": True,
        },
    ]

    assert docx_export._pgx_summary_rows(alerts) == [{
        "gene": "G6PD",
        "action": "建議考慮替代藥物。",
        "source": "CPIC",
        "level": "Strong",
        "drugs": ["Dapsone", "Rasburicase"],
    }]


def test_health_bundle_name_follows_selected_sections():
    assert docx_export._health_test_bundle_name({"acmg_sf", "pgx"}) == (
        "可採取醫療處置之疾病風險基因及藥物基因體學基因篩檢"
    )
    assert docx_export._health_test_bundle_name({"pgx"}) == "藥物基因體學基因篩檢"
    assert docx_export._health_test_bundle_name({"acmg_sf"}) == "可採取醫療處置之疾病風險基因篩檢"
