import json

from app.adapters.pgx_tsv import load_pgx


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

    assert result["guideline_annotations"][0]["fda_category"] == "potential_impact"

