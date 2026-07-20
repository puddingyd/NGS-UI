import csv

import pytest

from app.adapters.snv_tsv import (
    TIERS,
    classify_tier,
    load_snv_tsv,
    predicted_suspect_evidence,
)


def _row(**overrides):
    row = {
        "CLINVAR_SIG": "",
        "CLINVAR_STARS": "0",
        "LOFTEE": "",
        "ACMG_SCORE": "0",
        "PKNN_LLR": "",
        "ALPHAMISSENSE": "",
        "BAYESDEL_NOAF": "",
        "PANGOLIN_SCORE": "",
        "REVEL": "",
        "SPLICEAI_MAX": "",
    }
    row.update(overrides)
    return row


def test_snv_tiers_drop_legacy_clinvar_bucket_and_rename_other():
    assert TIERS == ["1A", "1B", "1C", "2"]
    assert classify_tier(_row(
        CLINVAR_SIG="Pathogenic", CLINVAR_STARS="0",
    )) == "2"
    assert classify_tier(_row(
        CLINVAR_SIG="Conflicting_classifications_of_pathogenicity",
        CLINVAR_SIGCONF="Pathogenic(1)|Uncertain_significance(2)",
    )) == "2"


def test_1a_and_1b_keep_priority_over_predicted_suspect():
    assert classify_tier(_row(
        CLINVAR_SIG="Pathogenic", CLINVAR_STARS="1", REVEL="0.95",
    )) == "1A"
    assert classify_tier(_row(LOFTEE="HC", SPLICEAI_MAX="0.8")) == "1B"


@pytest.mark.parametrize(
    "scores",
    [
        {"ACMG_SCORE": "4"},
        {"PKNN_LLR": "1"},
        {"ALPHAMISSENSE": "0.906"},
        {"BAYESDEL_NOAF": "0.27"},
        {"ALPHAMISSENSE": "0.792", "BAYESDEL_NOAF": "0.13"},
        {"PANGOLIN_SCORE": "-0.20"},
        {"REVEL": "0.773"},
        {"REVEL": "0.644", "ALPHAMISSENSE": "0.792"},
        {"REVEL": "0.644", "BAYESDEL_NOAF": "0.13"},
        {"SPLICEAI_MAX": "0.20"},
    ],
)
def test_each_planned_trigger_enters_predicted_suspect(scores):
    assert classify_tier(_row(**scores)) == "1C"


@pytest.mark.parametrize(
    "scores",
    [
        {},
        {"ACMG_SCORE": "3"},
        {"PKNN_LLR": "0.9999"},
        {"ALPHAMISSENSE": "0.905"},
        {"BAYESDEL_NOAF": "0.269"},
        {"ALPHAMISSENSE": "0.792", "BAYESDEL_NOAF": "0.129"},
        {"PANGOLIN_SCORE": "0.199"},
        {"REVEL": "0.772"},
        {"REVEL": "0.644", "ALPHAMISSENSE": "0.791"},
        {"SPLICEAI_MAX": "0.199"},
    ],
)
def test_below_threshold_or_missing_scores_remain_other(scores):
    assert classify_tier(_row(**scores)) == "2"


def test_evidence_keeps_core_and_extra_reasons_separate():
    evidence = predicted_suspect_evidence(_row(
        ACMG_SCORE="5",
        PKNN_LLR="1.25",
        PANGOLIN_SCORE="-0.31",
        REVEL="0.81",
        SPLICEAI_MAX="0.42",
    ))

    assert evidence["acmg_trigger"] is True
    assert evidence["core_trigger"] is True
    assert evidence["extra_trigger"] is True
    assert evidence["core_reasons"] == [
        "P-KNN LLR 1.25 ≥ 1",
        "|Pangolin -0.31| ≥ 0.2",
    ]
    assert evidence["extra_reasons"] == [
        "REVEL 0.81 ≥ 0.773",
        "SpliceAI 0.42 ≥ 0.2",
    ]
    assert evidence["reasons"][0] == "ACMG points 5 ≥ 4"


def test_pknn_includes_boundary_and_uses_max_multi_value():
    assert classify_tier(_row(PKNN_LLR="0.9999")) == "2"
    assert classify_tier(_row(PKNN_LLR="1")) == "1C"
    assert classify_tier(_row(PKNN_LLR=".&0.4&1.01")) == "1C"


def test_loader_returns_new_categories_and_does_not_rank_by_trigger_source(tmp_path):
    path = tmp_path / "sample.snv_indel.review.tsv"
    rows = [
        {
            "CHROM": "chr1", "POS": "10", "REF": "A", "ALT": "G",
            "GENE": "GENE1", "ACMG_CRITERIA": ".", "ACMG_SCORE": "3",
            "DP": "30", "REVEL": "0.80",
        },
        {
            "CHROM": "chr1", "POS": "20", "REF": "C", "ALT": "T",
            "GENE": "GENE2", "ACMG_CRITERIA": ".", "ACMG_SCORE": "0",
            "DP": "30", "PANGOLIN_SCORE": "0.30",
        },
        {
            "CHROM": "chr1", "POS": "30", "REF": "G", "ALT": "A",
            "GENE": "GENE3", "ACMG_CRITERIA": ".", "ACMG_SCORE": "0",
            "DP": "30",
        },
        {
            "CHROM": "chr1", "POS": "40", "REF": "T", "ALT": "C",
            "GENE": "GENE4", "ACMG_CRITERIA": ".", "ACMG_SCORE": "0",
            "DP": "30", "CLINVAR_SIG": "Pathogenic", "CLINVAR_STARS": "0",
        },
    ]
    fieldnames = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    _, categories = load_snv_tsv(path, test_type="WES")

    assert list(categories) == ["1A", "1B", "1C", "2"]
    # The Extra-only row stays ahead because its score is higher; predictor
    # source must not create a hidden Core-before-Extra ordering rule.
    assert categories["1C"] == ["chr1-10-A-G", "chr1-20-C-T"]
    assert categories["2"] == ["chr1-30-G-A", "chr1-40-T-C"]
