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
