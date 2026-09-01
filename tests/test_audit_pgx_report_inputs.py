import csv
import json

from scripts import audit_pgx_report_inputs as pgx_audit


def _source(label, allele1, allele2, phenotypes=None):
    return {
        "label": label,
        "allele1": {"name": allele1},
        "allele2": {"name": allele2},
        "phenotypes": phenotypes or [],
    }


def test_audit_projects_json_alleles_and_mt_rnr1_tsv_fallback(tmp_path):
    sample_dir = tmp_path / "26T00001-dragen" / "07_pgx"
    sample_dir.mkdir(parents=True)
    genes = {
        gene: {"sourceDiplotypes": [_source("Reference/Reference", "Reference", "Reference")]}
        for gene in pgx_audit.REPORT_GENES
    }
    genes["ABCG2"] = {"sourceDiplotypes": [_source(
        "rs2231142 reference (G)/rs2231142 variant (T)",
        "rs2231142 reference (G)",
        "rs2231142 variant (T)",
    )]}
    genes["CFTR"] = {"sourceDiplotypes": [_source(
        "Reference/Reference",
        "ivacaftor non-responsive CFTR sequence",
        "ivacaftor non-responsive CFTR sequence",
    )]}
    genes["DPYD"] = {
        "sourceDiplotypes": [
            _source("c.85T>C (*9A)/Reference", "c.85T>C (*9A)", "Reference"),
            _source("c.1627A>G (*5)", "c.1627A>G (*5)", ""),
        ],
        "variants": [
            {
                "dbSnpId": "rs1801265", "chromosome": "chr1", "position": 1,
                "call": "A/G", "referenceAllele": "A",
                "alleles": ["c.85T>C (*9A)"], "phased": False,
            },
            {
                "dbSnpId": "rs3918290", "chromosome": "chr1", "position": 2,
                "call": "C/T", "referenceAllele": "C",
                    "alleles": ["c.1627A>G (*5)"], "phased": False,
            },
        ],
    }
    genes["HLA-A"] = {"sourceDiplotypes": [_source(
        "*02:07/*11:01", "*02:07", "*11:01", ["*31:01 negative"],
    )]}
    genes["HLA-B"] = {"sourceDiplotypes": [_source(
        "*58:01/*58:01", "*58:01", "*58:01",
        ["*15:02 negative", "*57:01 negative", "*58:01 positive"],
    )]}
    genes["MT-RNR1"] = {"sourceDiplotypes": [_source(
        "Unknown", "Unknown", "", ["No Result"],
    )]}
    genes["VKORC1"] = {"sourceDiplotypes": [_source(
        "rs9923231 reference (C)/rs9923231 variant (T)",
        "rs9923231 reference (C)",
        "rs9923231 variant (T)",
    )]}
    (sample_dir / "26T00001.pharmcat.report.json").write_text(
        json.dumps({"genes": genes}), encoding="utf-8",
    )
    with (sample_dir / "26T00001.pgx.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("GENE", "DIPLOTYPE", "PHENOTYPE"),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow({"GENE": "MT-RNR1", "DIPLOTYPE": "Reference", "PHENOTYPE": "Normal"})

    report = pgx_audit.audit(tmp_path, "26T")
    rows = {row["gene"]: row for row in report["cases"]}

    assert report["sample_count"] == 1
    assert report["complete_file_pairs"] == 1
    assert (rows["ABCG2"]["allele1"], rows["ABCG2"]["allele2"]) == (
        "C (Reference)", "A (Variant)",
    )
    assert (rows["CFTR"]["allele1"], rows["CFTR"]["allele2"]) == (
        "Reference", "Reference",
    )
    assert rows["DPYD"]["test"] == "DPYD（相位未定）"
    assert rows["DPYD"]["allele1"] == rows["DPYD"]["allele2"] == ""
    assert rows["DPYD"]["result_span"] == (
        "檢出變異：c.85T>C (*9A) (het)；c.1627A>G (*5) (het)"
    )
    assert rows["HLA-B*58:01"]["allele1"] == "Positive"
    assert rows["HLA-B*58:01"]["allele2"] == "Positive"
    assert rows["MT-RNR1"]["source_rule"] == "TSV fallback"
    assert (rows["MT-RNR1"]["allele1"], rows["MT-RNR1"]["allele2"]) == (
        "Reference", "N/A",
    )
    assert (rows["VKORC1"]["allele1"], rows["VKORC1"]["allele2"]) == (
        "G (Reference)", "A (Variant)",
    )
    assert [row["label"] for row in report["source_diplotypes"] if row["gene"] == "DPYD"] == [
        "c.85T>C (*9A)/Reference", "c.1627A>G (*5)",
    ]
    assert report["nonreference_variants"] == [
        {
            "sample": "26T00001-dragen", "gene": "DPYD",
            "rsid": "rs1801265", "chromosome": "chr1", "position": "1",
            "call": "A/G", "reference": "A", "alleles": "c.85T>C (*9A)",
            "phased": "false", "phase_set": "",
        },
        {
            "sample": "26T00001-dragen", "gene": "DPYD",
            "rsid": "rs3918290", "chromosome": "chr1", "position": "2",
            "call": "C/T", "reference": "C", "alleles": "c.1627A>G (*5)",
            "phased": "false", "phase_set": "",
        },
    ]

    output_dir = tmp_path / "audit"
    pgx_audit.write_outputs(report, output_dir)
    assert (output_dir / "pgx_audit_report.json").is_file()
    assert (output_dir / "pgx_audit_observed.tsv").is_file()
    assert (output_dir / "pgx_audit_source_diplotypes.tsv").is_file()
    assert (output_dir / "pgx_audit_nonreference_variants.tsv").is_file()
