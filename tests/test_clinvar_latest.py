import csv
import gzip
import importlib.util
import json
from pathlib import Path

from app.adapters.snv_tsv import _row_to_variant
from app.services import clinvar_latest_store, sample_loader


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "annotate_clinvar_latest.py"
SPEC = importlib.util.spec_from_file_location("annotate_clinvar_latest", SCRIPT)
annotator = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(annotator)

UPDATE_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "update_clinvar_latest_db.py"
UPDATE_SPEC = importlib.util.spec_from_file_location("update_clinvar_latest_db", UPDATE_SCRIPT)
updater = importlib.util.module_from_spec(UPDATE_SPEC)
assert UPDATE_SPEC and UPDATE_SPEC.loader
UPDATE_SPEC.loader.exec_module(updater)


def _write_vcf(path: Path) -> None:
    rows = [
        ("1", 100, "1001", "A", "G", "Pathogenic", "reviewed_by_expert_panel"),
        ("1", 200, "1002", "C", "T", "Uncertain_significance", "criteria_provided,_single_submitter"),
        ("1", 300, "1003", "G", "A", "Pathogenic", "reviewed_by_expert_panel"),
        ("1", 400, "1004", "T", "C", "Uncertain_significance", "criteria_provided,_single_submitter"),
        ("1", 600, "1006", "A", "C", "Likely_pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts"),
    ]
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n##fileDate=2026-08-03\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, identifier, ref, alt, significance, review in rows:
            handle.write(
                f"{chrom}\t{pos}\t{identifier}\t{ref}\t{alt}\t.\tPASS\t"
                f"CLNSIG={significance};CLNREVSTAT={review};CLNDN=Disease_{pos}\n"
            )


def _write_tsv(path: Path) -> None:
    fields = [
        "CHROM", "POS", "REF", "ALT", "CLINVAR_SIG", "CLINVAR_STARS",
        "CLINVAR_DN", "CLINVAR_SIGCONF", "CLINVAR_VARIATION_ID",
    ]
    rows = [
        ["chr1", "100", "A", "G", "Uncertain significance", "1", "old", "", "1001"],
        ["chr1", "200", "C", "T", "Pathogenic", "3", "old", "", "1002"],
        ["chr1", "300", "G", "A", "Likely pathogenic", "2", "old", "", "1003"],
        ["chr1", "400", "T", "C", "Conflicting classifications", "1", "old", "x", "1004"],
        ["chr1", "500", "G", "T", "Pathogenic", "2", "old", "", "9999"],
        ["chr1", "600", "A", "C", "", "", "", "", ""],
        ["chr1", "700", "C", "G", "", "", "", "", ""],
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def test_weekly_clinvar_comparison_and_arrow_policy(tmp_path):
    vcf = tmp_path / "clinvar.vcf.gz"
    db = tmp_path / "clinvar.sqlite"
    tsv = tmp_path / "sample.tsv"
    marker = tmp_path / "sample.clinvar_comparison.json"
    _write_vcf(vcf)
    _write_tsv(tsv)
    built = clinvar_latest_store.build_database(vcf, db)

    result = annotator.annotate_tsv(
        tsv, db, marker, baseline_release="2026-07-20"
    )

    assert built["release_date"] == "2026-08-03"
    with tsv.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["CLINVAR_CHANGE"] == "UP_TO_PLP"
    assert rows[1]["CLINVAR_CHANGE"] == "DOWN_FROM_PLP"
    assert rows[2]["CLINVAR_CHANGE"] == ""  # LP -> P
    assert rows[3]["CLINVAR_CHANGE"] == ""  # conflicting -> VUS
    assert rows[4]["CLINVAR_CHANGE"] == "DOWN_FROM_PLP"  # verified no record
    assert rows[5]["CLINVAR_CHANGE"] == "UP_TO_PLP"  # new record by allele
    assert rows[6]["CLINVAR_LATEST_APPLIED"] == ""  # unknown miss is preserved
    assert rows[0]["CLINVAR_SIG"] == "Uncertain significance"
    assert rows[0]["CLINVAR_LATEST_SIG"] == "Pathogenic"
    assert rows[0]["CLINVAR_LATEST_REVIEW_STATUS"] == "reviewed_by_expert_panel"
    assert rows[1]["CLINVAR_SIG"] == "Pathogenic"
    assert rows[1]["CLINVAR_LATEST_SIG"] == "Uncertain_significance"
    assert rows[4]["CLINVAR_SIG"] == "Pathogenic"
    assert rows[4]["CLINVAR_LATEST_SIG"] == ""
    assert rows[5]["CLINVAR_SIG"] == ""
    assert rows[5]["CLINVAR_LATEST_SIG"] == "Likely_pathogenic"
    assert result["latest_release"] == "2026-08-03"
    assert result["stats"]["up_to_plp"] == 2
    assert result["stats"]["down_from_plp"] == 2


def test_updater_atomically_publishes_database_vcf_and_manifest(tmp_path):
    source = tmp_path / "source.vcf.gz"
    _write_vcf(source)
    db = tmp_path / "published" / "clinvar.sqlite"
    vcf = tmp_path / "published" / "clinvar.vcf.gz"
    manifest_path = tmp_path / "published" / "manifest.json"

    result = updater.update(
        vcf_file=source,
        url="https://example.invalid/clinvar.vcf.gz",
        db_path=db,
        vcf_path=vcf,
        manifest_path=manifest_path,
    )

    assert db.is_file() and vcf.is_file() and manifest_path.is_file()
    assert result["release_date"] == "2026-08-03"
    assert json.loads(manifest_path.read_text(encoding="utf-8"))["record_count"] == 5
    by_key, _by_id = clinvar_latest_store.lookup_records(db, ["1:100:A:G"])
    assert by_key["1:100:A:G"]["review_status"] == "reviewed_by_expert_panel"
    assert clinvar_latest_store.metadata(db)["schema_version"] == "2"


def test_adapter_exposes_latest_and_report_export_has_fixed_baseline_restore():
    variant = _row_to_variant({
        "CHROM": "1", "POS": "100", "REF": "A", "ALT": "G",
        "CLINVAR_SIG": "Uncertain significance", "CLINVAR_STARS": "1",
        "CLINVAR_DN": "old disease", "CLINVAR_VARIATION_ID": "9001",
        "CLINVAR_BASE_SIG": "Uncertain significance",
        "CLINVAR_BASE_STARS": "1", "CLINVAR_BASE_DN": "old disease",
        "CLINVAR_BASE_VARIATION_ID": "9001", "CLINVAR_LATEST_APPLIED": "1",
        "CLINVAR_LATEST_SIG": "Pathogenic", "CLINVAR_LATEST_STARS": "3",
        "CLINVAR_LATEST_DN": "new disease",
        "CLINVAR_LATEST_VARIATION_ID": "1001",
        "CLINVAR_LATEST_REVIEW_STATUS": "reviewed_by_expert_panel",
        "CLINVAR_CHANGE": "UP_TO_PLP",
    })

    assert variant["clinvar_change"] == "UP_TO_PLP"
    assert variant["CLNSIG"] == "Uncertain significance"
    assert variant["CLNSIG_old"] == "Uncertain significance"
    assert variant["clinvar_stars_old"] == 1
    assert variant["clinvar_dn_old"] == "old disease"
    assert variant["clinvar_variation_id_old"] == "9001"
    assert variant["clinvar_latest_sig"] == "Pathogenic"
    assert variant["clinvar_latest_stars"] == 3
    assert variant["clinvar_latest_review_status"] == "reviewed_by_expert_panel"
    assert variant["tier"] == "2"  # weekly upgrade must not replace baseline tier
    export_source = (
        Path(__file__).resolve().parents[1] / "backend/app/services/docx_export.py"
    ).read_text(encoding="utf-8")
    assert "clinvar_latest_store.restore_pipeline_variants(variants)" in export_source
    assert "clinvar_baseline=True" in export_source
    assert 'fixed["clinvar_date"] = "2026-07-20"' in export_source


def test_health_report_secondary_candidates_use_pipeline_clinvar_baseline():
    variants = {
        "weekly-upgrade-only": {
            "CLNSIG": "Pathogenic", "CLNSIG_old": "Uncertain significance",
            "clinvar_latest_applied": True, "tier": "1A",
        },
        "weekly-downgrade": {
            "CLNSIG": "Uncertain significance", "CLNSIG_old": "Pathogenic",
            "clinvar_latest_applied": True, "tier": "2",
        },
        "weekly-upgrade-with-loftee": {
            "CLNSIG": "Pathogenic", "CLNSIG_old": "Uncertain significance",
            "clinvar_latest_applied": True, "tier": "1A", "loftee_hc": True,
        },
    }

    restored = sample_loader._pipeline_clinvar_secondary_variants(variants)

    assert restored["weekly-upgrade-only"]["CLNSIG"] == "Uncertain significance"
    assert restored["weekly-upgrade-only"]["tier"] == "2"
    assert restored["weekly-downgrade"]["CLNSIG"] == "Pathogenic"
    assert restored["weekly-upgrade-with-loftee"]["tier"] == "1B"
    assert variants["weekly-upgrade-only"]["CLNSIG"] == "Pathogenic"


def test_erepo_dot_placeholders_are_treated_as_absent():
    variant = _row_to_variant({
        "CHROM": "1", "POS": "101", "REF": "A", "ALT": "T",
        "CLINGEN_VCEP_CLASS": ".", "CLINGEN_VCEP_CRITERIA": ".",
        "CLINGEN_VCEP_PANEL": ".", "CLINGEN_AGREEMENT": ".",
    })

    assert variant["clingen_vcep_class"] == ""
    assert variant["clingen_vcep_criteria"] == ""
    assert variant["clingen_vcep_panel"] == ""


def test_erepo_exposes_expert_class_and_derived_criteria_score():
    variant = _row_to_variant({
        "CHROM": "1", "POS": "102", "REF": "G", "ALT": "A",
        "CLINGEN_VCEP_CLASS": "Likely_pathogenic",
        "CLINGEN_VCEP_CRITERIA": "PVS1,PM2_Supporting",
        "CLINGEN_VCEP_PANEL": "Hearing Loss VCEP",
    })

    assert variant["clingen_vcep_class"] == "Likely pathogenic"
    assert variant["clingen_vcep_score"] == 9
    assert variant["clingen_vcep_criteria"] == "PVS1,PM2_Supporting"
