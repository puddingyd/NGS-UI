import json

from app.adapters.snv_tsv import _min_multi
from app.services import annotation_versions


def test_reads_preferred_sample_sidecar(tmp_path):
    raw = tmp_path / "VAL-1.snv_indel.acmg.tsv"
    raw.write_text("CHROM\tPOS\n", encoding="utf-8")
    (tmp_path / "VAL-1.annotation_versions.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "databases": {
                    "clinvar": {
                        "release_date": "2026-05-10",
                        "source": "clinvar_20260510.vcf.gz",
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    result = annotation_versions.load_annotation_versions(raw)

    assert result["clinvar_date"] == "2026-05-10"
    assert result["clinvar"]["release_date"] == "2026-05-10"
    assert result["metadata_path"].endswith("VAL-1.annotation_versions.json")


def test_accepts_pipeline_source_migration_fallback(tmp_path):
    raw = tmp_path / "S1.snv_indel.acmg.tsv"
    raw.touch()
    state = tmp_path / "08_postprocessing"
    state.mkdir()
    (state / "pipeline_source.json").write_text(
        json.dumps(
            {
                "annotation_versions": {
                    "clinvar": {"version": "clinvar_20260426.vcf.gz"}
                }
            }
        ),
        encoding="utf-8",
    )

    result = annotation_versions.load_annotation_versions(raw, state)

    assert result["clinvar_date"] == "2026-04-26"


def test_missing_or_invalid_metadata_does_not_invent_date(tmp_path):
    raw = tmp_path / "S1.snv_indel.acmg.tsv"
    raw.touch()
    (tmp_path / "S1.annotation_versions.json").write_text(
        json.dumps({"databases": {"clinvar": {"release_date": "latest"}}}),
        encoding="utf-8",
    )

    assert annotation_versions.load_annotation_versions(raw) == {}


def test_lower_is_more_pathogenic_predictors_take_worst_transcript():
    # ESM1b and SIFT both use lower scores for greater predicted impact.
    assert _min_multi(".&0.08&0.001&0.32") == 0.001
