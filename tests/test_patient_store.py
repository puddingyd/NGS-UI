from app.services import patient_store, sample_loader


def test_delete_unregisters_sample_without_removing_pipeline_outputs(tmp_path, monkeypatch):
    monkeypatch.setattr(patient_store, "TERTIARY_OUTPUT_ROOT", tmp_path)
    monkeypatch.setattr(patient_store, "PIPELINE_OUT_ROOT", tmp_path / "pipeline")
    monkeypatch.setattr(sample_loader, "TERTIARY_OUTPUT_ROOT", tmp_path)

    from app.services import dragen_jobs, patient_list_store

    monkeypatch.setattr(dragen_jobs, "active_sample_ids", lambda: set())
    monkeypatch.setattr(patient_list_store, "load_roster", lambda: {})
    monkeypatch.setattr(
        patient_list_store,
        "lookup_with_key",
        lambda lis_id, source_sample_id, roster=None: ({}, lis_id),
    )
    monkeypatch.setattr(
        patient_list_store,
        "lookup_candidates",
        lambda lis_id, roster_lis_id="", source_sample_id="": [lis_id],
    )

    sample_dir = tmp_path / "S1"
    sample_dir.mkdir()
    (sample_dir / "snv_indel.annotated.tsv").write_text(
        "CHROM\tPOS\tREF\tALT\n",
        encoding="utf-8",
    )
    (sample_dir / "pgx.tsv").write_text("gene\tphenotype\n", encoding="utf-8")
    (sample_dir / "sample_metadata.json").write_text(
        '{"sample_id":"S1","lis_id":"S1","name":"Patient","mrn":"MRN"}',
        encoding="utf-8",
    )
    (sample_dir / "case_summary.json").write_text("{}", encoding="utf-8")
    analysis_dir = sample_dir / "analyses" / "default"
    analysis_dir.mkdir(parents=True)
    (analysis_dir / "analysis.json").write_text("{}", encoding="utf-8")

    result = patient_store.delete("S1")

    assert result["unregistered"] is True
    assert sample_dir.is_dir()
    assert (sample_dir / "snv_indel.annotated.tsv").is_file()
    assert (sample_dir / "pgx.tsv").is_file()
    assert not (sample_dir / "sample_metadata.json").exists()
    assert not (sample_dir / "case_summary.json").exists()
    assert not (sample_dir / "analyses").exists()
    assert [row["lis_id"] for row in sample_loader.list_unregistered()] == ["S1"]
