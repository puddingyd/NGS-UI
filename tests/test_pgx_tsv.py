import json

from app.adapters.pgx_tsv import load_pgx
from app.services import sample_loader


def test_load_pgx_preserves_compact_fda_association(tmp_path):
    tsv = tmp_path / "sample.pgx.tsv"
    tsv.write_text(
        "GENE\tDRUG\tPHENOTYPE\tGUIDELINE_SOURCE\tRECOMMENDATION\tCPIC_LEVEL\n"
        "CYP2C19\tclopidogrel\tPoor Metabolizer\tCPIC\tUse an alternative\tStrong\n",
        encoding="utf-8",
    )
    report = {
        "drugs": {
            "FDA PGx Association": {
                "clopidogrel": {
                    "name": "clopidogrel",
                    "source": "FDA_ASSOCIATION",
                    "guidelines": [{
                        "name": "Pharmacogenetic Associations for which the Data Support Therapeutic Management Recommendations",
                        "source": "FDA_ASSOCIATION",
                        "url": "https://example.test/fda",
                        "annotations": [{
                            "drugRecommendation": "<p>Consider an alternative.</p>",
                            "classification": "Unspecified",
                            "lookupKey": [{"CYP2C19": "Poor Metabolizer"}],
                            "dosingInformation": True,
                        }],
                    }],
                },
            },
        },
    }
    json_path = tmp_path / "sample.pharmcat.report.json"
    json_path.write_text(json.dumps(report), encoding="utf-8")

    result = load_pgx(tsv, json_path)

    assert result["pharmcat_available"] is True
    assert result["genes"]["CYP2C19"]["drugs"][0]["drug"] == "clopidogrel"
    assert result["guideline_annotations"] == [{
        "section": "FDA PGx Association",
        "source": "FDA_ASSOCIATION",
        "drug": "clopidogrel",
        "guideline": "Pharmacogenetic Associations for which the Data Support Therapeutic Management Recommendations",
        "url": "https://example.test/fda",
        "classification": "Unspecified",
        "recommendation": "Consider an alternative.",
        "implications": [],
        "genes": ["CYP2C19"],
        "fda_category": "therapeutic_management",
        "dosing_information": True,
        "alternate_drug_available": False,
        "other_prescribing_guidance": False,
    }]


def test_load_pgx_without_tsv_still_returns_json_annotations(tmp_path):
    report = {
        "drugs": {
            "FDA PGx Association": {
                "drug-a": {
                    "guidelines": [{
                        "name": "Potential Impact on Safety or Response",
                        "annotations": [{
                            "drugRecommendation": "Monitor response.",
                            "lookupKey": [{"GENE1": "result"}],
                        }],
                    }],
                },
            },
        },
    }
    json_path = tmp_path / "sample.pharmcat.report.json"
    json_path.write_text(json.dumps(report), encoding="utf-8")

    result = load_pgx(tmp_path / "missing.tsv", json_path)

    assert result["pharmcat_available"] is True
    assert result["guideline_annotations"][0]["fda_category"] == "potential_impact"


def test_load_pgx_marks_missing_json_explicitly(tmp_path):
    result = load_pgx(tmp_path / "missing.tsv", tmp_path / "missing.json")

    assert result["pharmcat_available"] is False


def test_load_pgx_omits_unrenderable_object_messages(tmp_path):
    report = {
        "messages": [{"unexpected": "object"}, {"message": "Readable global message"}],
        "genes": {
            "CYP2B6": {
                "sourceDiplotypes": [{"label": "*1/*4"}],
                "messages": [
                    {"unexpected": "object"},
                    {"text": "Readable gene message"},
                ],
            },
        },
    }
    json_path = tmp_path / "sample.pharmcat.report.json"
    json_path.write_text(json.dumps(report), encoding="utf-8")

    result = load_pgx(tmp_path / "missing.tsv", json_path)

    assert result["messages"] == ["Readable global message"]
    assert result["genes"]["CYP2B6"]["details"]["messages"] == ["Readable gene message"]


def test_load_pgx_preserves_all_source_diplotypes_and_variant_phase(tmp_path):
    report = {
        "genes": {
            "DPYD": {
                "phased": False,
                "effectivelyPhased": False,
                "sourceDiplotypes": [
                    {
                        "label": "c.1627A>G(*5)",
                        "allele1": {"name": "c.1627A>G(*5)"},
                        "allele2": None,
                        "phenotypes": ["Indeterminate"],
                    },
                    {
                        "label": "c.1896T>C",
                        "allele1": {"name": "c.1896T>C"},
                        "allele2": None,
                        "phenotypes": ["Indeterminate"],
                    },
                ],
                "variants": [
                    {
                        "dbSnpId": "rs1801159",
                        "chromosome": "chr1",
                        "position": 1,
                        "call": "T/C",
                        "referenceAllele": "T",
                        "alleles": ["c.1627A>G(*5)"],
                        "phased": False,
                        "phaseSet": None,
                    },
                    {
                        "dbSnpId": "rs17376848",
                        "chromosome": "chr1",
                        "position": 2,
                        "call": "A/G",
                        "referenceAllele": "A",
                        "alleles": ["c.1896T>C"],
                        "phased": False,
                        "phaseSet": None,
                    },
                ],
            },
        },
    }
    json_path = tmp_path / "sample.pharmcat.report.json"
    json_path.write_text(json.dumps(report), encoding="utf-8")

    result = load_pgx(tmp_path / "missing.tsv", json_path)
    details = result["genes"]["DPYD"]["details"]

    assert [row["allele1_name"] for row in details["source_diplotypes"]] == [
        "c.1627A>G(*5)", "c.1896T>C",
    ]
    assert details["effectively_phased"] is False
    assert details["variants"] == [
        {
            "rsid": "rs1801159", "chr": "chr1", "pos": 1,
            "call": "T/C", "alleles": "c.1627A>G(*5)",
            "allele_names": ["c.1627A>G(*5)"], "reference": "T",
            "phased": False, "phase_set": "",
        },
        {
            "rsid": "rs17376848", "chr": "chr1", "pos": 2,
            "call": "A/G", "alleles": "c.1896T>C",
            "allele_names": ["c.1896T>C"], "reference": "A",
            "phased": False, "phase_set": "",
        },
    ]


def test_staged_pgx_loader_attaches_shared_report_view(tmp_path, monkeypatch):
    sample_dir = tmp_path / "S1"
    sample_dir.mkdir()
    json_path = sample_dir / "sample.pharmcat.report.json"
    json_path.write_text(json.dumps({
        "pharmcatVersion": "3.2.0",
        "genes": {
            "CYP2C19": {
                "sourceDiplotypes": [{
                    "label": "*2/*2",
                    "phenotypes": ["Poor Metabolizer"],
                }],
            },
        },
    }), encoding="utf-8")
    monkeypatch.setattr(sample_loader.sample_layout, "state_dir", lambda _sample: sample_dir)
    monkeypatch.setattr(
        sample_loader.sample_layout,
        "pgx_tsv",
        lambda _sample: sample_dir / "missing.pgx.tsv",
    )
    monkeypatch.setattr(sample_loader.sample_layout, "pharmcat_json", lambda _sample: json_path)

    payload = sample_loader.load_sample_pgx("S1")

    assert payload is not None
    assert payload["pgx"]["pharmcat_available"] is True
    assert len(payload["pgx"]["report_view"]["genotype_rows"]) == 21
    assert payload["pgx"]["report_view"]["genotype_rows"][4] == {
        "gene": "CYP2C19",
        "genotype": "*2/*2",
        "phenotype": "Poor Metabolizer",
    }
