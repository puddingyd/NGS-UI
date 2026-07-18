import json

from app import config
from app.workers import exomiser_lirical


def test_worker_skips_before_vcf_preparation_without_hpo(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", tmp_path)
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", tmp_path / "pipeline")
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", tmp_path / "old_pipeline")

    sample_dir = tmp_path / "S1"
    analysis_dir = sample_dir / "analyses" / "default"
    analysis_dir.mkdir(parents=True)
    (sample_dir / "sample_metadata.json").write_text(
        json.dumps({"active_analysis": "default", "genome_build": "hg38"}),
        encoding="utf-8",
    )
    (analysis_dir / "analysis.json").write_text(
        json.dumps({"hpo": [], "selected_panels": [{"name": "panel"}]}),
        encoding="utf-8",
    )

    updates = []
    monkeypatch.setattr(
        exomiser_lirical.job_store,
        "update",
        lambda job_id, patch: updates.append((job_id, patch)) or patch,
    )
    monkeypatch.setattr(
        exomiser_lirical.vcf_writer,
        "needs_rebuild",
        lambda sample_id: (_ for _ in ()).throw(AssertionError("VCF preparation must be skipped")),
    )

    result = exomiser_lirical.run_exomiser_lirical("job-1", "S1")

    assert result["status"] == "succeeded"
    assert result["step"] == "skipped:no-hpo"
    assert result["skipped"] is True
    assert updates == [("job-1", result)]
