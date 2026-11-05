"""Phase-2 report scenarios reconstructed with synthetic data, not patient files."""
import copy
import json

import pytest

from app.adapters.pgx_tsv import load_pgx
from app.services import docx_export, pgx_concordance as pc


def gene(label, **details):
    return {"details": {"label": label, "source_diplotypes": [{"label": label}], **details}}


@pytest.mark.parametrize("platform,dragen", [
    ("B (reference)/B (reference)", "B(wildtype)"),
    ("B (reference)", "B(wildtype)/B(wildtype)"),
    ("Radlowo/Radlowo", "Radlowo"),
])
def test_g6pd_identical_homozygote_and_haploid_are_silent(platform, dragen):
    result = pc.compare_gene("G6PD", gene(platform), {"genotype": dragen})
    assert result["status"] == "same" and result["badges"] == []


def test_g6pd_true_heterozygous_call_is_not_collapsed_to_reference():
    result = pc.compare_gene("G6PD", gene("B (reference)/Radlowo"), {"genotype": "B(wildtype)"})
    assert result["status"] == "different"


@pytest.mark.parametrize("platform,dragen", [
    ("c.1896T>C (heterozygous)", "c.1896T>C/Reference"),
    ("c.1627A>G (*5) (het)", "Reference/c.1627A>G(*5)"),
    ("c.1896T>C (hom)", "c.1896T>C/c.1896T>C"),
    ("c.1627A>G (*5)/c.1896T>C", "c.1896T>C/c.1627A>G"),
    ("c.1905+1G>A (het)", "Reference/c.1905+1G>A"),
])
def test_explicit_hgvs_zygosity_and_star_annotation_are_equivalent(platform, dragen):
    result = pc.compare_gene("DPYD", gene(platform), {"genotype": dragen})
    assert result["status"] == "same" and result["badges"] == []


def test_explicit_homozygote_does_not_match_heterozygote():
    row = pc.compare_gene("DPYD", gene("c.1896T>C (hom)"), {"genotype": "c.1896T>C/Reference"})
    assert row["status"] == "different"


def component_case(gene="DPYD"):
    details = {
        "label": "c.1627A>G (*5)", "phased": False, "effectively_phased": False,
        "source_diplotypes": [
            {"label": "c.1627A>G (*5)", "allele1_name": "c.1627A>G (*5)"},
            {"label": "c.1896T>C", "allele1_name": "c.1896T>C"},
        ],
        "comparison_variants": [
            {"allele_names": ["c.1627A>G (*5)"], "call": "T/C", "reference": "T"},
            {"allele_names": ["c.1896T>C"], "call": "A/G", "reference": "A"},
        ],
    }
    return {"gene": gene, "details": details}


@pytest.mark.parametrize("gene_name", ["DPYD", "RYR1"])
def test_complete_component_sources_are_compared_and_presented_without_inventing_phase(gene_name):
    payload = component_case(gene_name)
    row = pc.compare_gene(gene_name, payload, {"genotype": "c.1627A>G/c.1896T>C"})
    assert row["status"] == "phase_compatible" and row["badges"] == ["相位未定"]
    pgx = {"genes": {gene_name: payload}}
    display, _ = docx_export._pgx_gene_result(pgx, gene_name)
    assert "c.1627A>G" in display and "c.1896T>C" in display and "相位未定" in display
    result = next(r for r in docx_export._pgx_health_genotype_rows(pgx) if r["test"].startswith(gene_name))
    assert "c.1627A>G" in result["result_span"] and "c.1896T>C" in result["result_span"]
    assert pc.compare_gene(gene_name, payload, {"genotype": "c.1627A>G/Reference"})["status"] == "different"


def test_observation_compatibility_preserves_dosage_and_does_not_assign_cis_trans():
    payload = component_case()
    row = pc.compare_gene("DPYD", payload, {"genotype": "[c.1627A>G+c.1896T>C]/Reference"})
    assert row["badges"] == ["相位未定"]
    payload["details"]["comparison_variants"][0]["call"] = "C/C"
    assert pc.compare_gene("DPYD", payload, {"genotype": "c.1627A>G/c.1896T>C"})["status"] == "different"


def test_known_phase_cannot_be_verified_from_a_bag_of_source_components():
    payload = component_case()
    payload["details"]["phased"] = True
    assert pc.compare_gene("DPYD", payload, {"genotype": "c.1627A>G/c.1896T>C"})["status"] == "uncomparable"


def test_alternative_complete_diplotypes_are_not_treated_as_components():
    payload = gene("*1/*28", source_diplotypes=[{"label": "*1/*28"}, {"label": "*1/*80+*28"}])
    row = pc.compare_gene("UGT1A1", payload, {"genotype": "*1/*28"})
    assert row["status"] == "same" and row["badges"] == ["相位未定"]
    assert docx_export._pgx_unphased_variants({"genes": {"DPYD": payload}}, "DPYD") == []


def nudt15():
    payload = gene("*1/*3", definition_source="CLINPGX", definition_version="2026-02-09-10-28",
                   phased=False, effectively_phased=False,
                   comparison_variants=[
                       {"chr": "chr13", "pos": 48037782, "rsid": "rs746071566", "reference": "AGGAGTC", "call": "AGGAGTC/AGGAGTCGGAGTC"},
                       {"chr": "chr13", "pos": 48045719, "rsid": "rs116855232", "reference": "C", "call": "C/T"},
                   ])
    dragen = {"genotype": "*2/*1;*3/*6", "starAlleleDatabaseSource": "PharmGKB Database Version: Snapshot-2022.01.01"}
    return payload, dragen


def test_nudt15_validated_two_heterozygous_unphased_sites_get_only_phase_warning():
    payload, dragen = nudt15()
    original = copy.deepcopy(payload)
    row = pc.compare_gene("NUDT15", payload, dragen)
    assert row["status"] == "phase_compatible" and row["badges"] == ["相位未定"]
    assert payload == original  # no phenotype/function equivalence is inferred


@pytest.mark.parametrize("changed", ["missing", "homozygous", "phased", "phase_set", "warnings", "different_version", "different_definition", "different_candidate", "wrong_alleles"])
def test_nudt15_does_not_hide_differences_without_the_validated_evidence(changed):
    payload, dragen = nudt15(); details = payload["details"]
    v = details["comparison_variants"][0]
    if changed == "missing": details["comparison_variants"] = []
    elif changed == "homozygous": v["call"] = "AGGAGTC/AGGAGTC"
    elif changed == "phased": v.update(phased=True, call="AGGAGTC|AGGAGTCGGAGTC")
    elif changed == "phase_set": v["phase_set"] = "synthetic-PS"
    elif changed == "warnings": v["warnings"] = ["unresolved"]
    elif changed == "different_version": details["definition_version"] = "2099-01-01"
    elif changed == "different_definition": dragen["starAlleleDatabaseSource"] = "unknown"
    elif changed == "different_candidate": dragen["genotype"] = "*2/*1;*4/*4"
    elif changed == "wrong_alleles": v["call"] = "AGGAGTC/A"
    assert pc.compare_gene("NUDT15", payload, dragen)["status"] == "different"


def test_adapter_keeps_all_comparison_evidence_even_when_display_variants_are_limited(tmp_path):
    payload, dragen = nudt15()
    raw = []
    for i in range(100):
        raw.append({"chromosome": "chr13", "position": i + 1, "call": "C/C", "referenceAllele": "C", "alleles": []})
    for v in payload["details"]["comparison_variants"]:
        raw.append({"chromosome": v["chr"], "position": v["pos"], "dbSnpId": v["rsid"], "call": v["call"], "referenceAllele": v["reference"]})
    report = {"genes": {"NUDT15": {"alleleDefinitionSource": "CLINPGX", "alleleDefinitionVersion": "2026-02-09-10-28",
                                    "sourceDiplotypes": [{"label": "*1/*3"}], "variants": raw}}}
    path = tmp_path / "synthetic.report.json"; path.write_text(json.dumps(report))
    adapted = load_pgx(tmp_path / "missing.tsv", path)["genes"]["NUDT15"]
    assert len(adapted["details"]["comparison_variants"]) == 102
    assert len(adapted["details"]["variants"]) == 2
    assert pc.compare_gene("NUDT15", adapted, dragen)["badges"] == ["相位未定"]


def test_true_cyp2d6_copy_number_or_distribution_differences_are_still_visible():
    for platform, target in [("*10/*36", "*10+*36/*5"),
                             ("*36+*10/*36+*10", "*10/*10x2+*36x2"),
                             ("*1/*36+*10", "*1/*10+*36x2"), ("*2/*2", "*2/*75")]:
        assert pc.compare_gene("CYP2D6", gene(platform), {"genotype": target})["status"] == "different"


def test_old_report_does_not_create_global_reference_exemptions():
    assert pc.compare_gene("UGT1A1", gene("*1/*6"), {"genotype": "*1/*1"})["status"] == "different"
    assert pc.compare_gene("RYR1", gene("c.11266C>G (het)"), {"genotype": "Reference/Reference"})["status"] == "different"
