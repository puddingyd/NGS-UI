from app import config
from app.services import dragen_jobs, sample_layout, sample_loader


def test_delete_unified_pipeline_output_invalidates_original_state_dir(tmp_path, monkeypatch):
    unified_root = tmp_path / "unified"
    legacy_pipeline_root = tmp_path / "legacy_pipeline"
    legacy_ui_root = tmp_path / "legacy_ui"
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", unified_root)
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", legacy_pipeline_root)
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", legacy_ui_root)
    monkeypatch.setattr(dragen_jobs, "PIPELINE_OUT_ROOT", unified_root)
    monkeypatch.setattr(dragen_jobs, "LEGACY_PIPELINE_OUT_ROOT", legacy_pipeline_root)
    monkeypatch.setattr(dragen_jobs, "list_jobs", lambda limit=1000: [])

    raw_tsv = unified_root / "S1" / "03_acmg" / "S1.snv_indel.acmg.tsv"
    raw_tsv.parent.mkdir(parents=True)
    raw_tsv.write_text("CHROM\tPOS\tREF\tALT\n", encoding="utf-8")
    marker = sample_layout.write_layout_marker("S1", source_id="S1", raw_tsv=raw_tsv)
    state_dir = marker.parent

    invalidated = []
    removed_rows = []
    monkeypatch.setattr(sample_loader, "invalidate_sample_cache", invalidated.append)
    monkeypatch.setattr(sample_loader, "remove_case_table_row", removed_rows.append)

    result = dragen_jobs.delete_pipeline_output("S1")

    assert not (unified_root / "S1").exists()
    assert result == {
        "sample_id": "S1",
        "source_sample_id": "S1",
        "deleted": [str(unified_root / "S1")],
    }
    assert invalidated == [state_dir]
    assert removed_rows == ["S1"]
