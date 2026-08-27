import json

from app import config
from app.services import (
    patient_list_store,
    patient_phenotype_store,
    patient_store,
    phenotype_io,
    sample_loader,
)


def _legacy_sample(tmp_path, monkeypatch, sample_id="S1"):
    old_ui = tmp_path / "old_ui"
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", old_ui)
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", tmp_path / "pipeline")
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", tmp_path / "old_pipeline")
    sample = old_ui / sample_id
    sample.mkdir(parents=True)
    (sample / "snv_indel.annotated.tsv").write_text(
        "CHROM\tPOS\tREF\tALT\n", encoding="utf-8",
    )
    return sample


def test_mrn_only_snapshot_wins_over_legacy_lis_file(tmp_path, monkeypatch):
    monkeypatch.setattr(patient_phenotype_store, "PHENOTYPE_DIR", tmp_path)
    (tmp_path / "S1_MRN1_phenotype.txt").write_text(
        "phenotype\thpo_name\tweight\nHP:0000001\tLegacy\t1\n",
        encoding="utf-8",
    )
    saved = patient_phenotype_store.save(
        mrn="MRN1",
        code="S2",
        hpo=[{"phenotype": "HP:0001250", "label": "Seizure", "weight": 2}],
        panels=[{"name": "WES-I__Neuro", "weight": 3}],
    )

    assert saved["filename"] == "MRN1_phenotype.txt"
    assert not (tmp_path / "S2_MRN1_phenotype.txt").exists()
    loaded = patient_phenotype_store.load(
        mrn="MRN1", code="S1", code_candidates=["S2"],
    )
    assert loaded["filename"] == "MRN1_phenotype.txt"
    assert loaded["hpo"][0]["phenotype"] == "HP:0001250"
    assert loaded["panels"][0] == {"name": "WES-I__Neuro", "weight": 3.0}


def test_empty_snapshot_is_a_header_only_authoritative_marker(tmp_path, monkeypatch):
    monkeypatch.setattr(patient_phenotype_store, "PHENOTYPE_DIR", tmp_path)
    (tmp_path / "S1_MRN1_phenotype.txt").write_text(
        "phenotype\thpo_name\tweight\nHP:0001250\tStale\t1\n",
        encoding="utf-8",
    )

    patient_phenotype_store.save(mrn="MRN1", code="S1", hpo=[], panels=[])
    loaded = patient_phenotype_store.load(mrn="MRN1", code="S1")

    assert loaded["content"] == "phenotype\thpo_name\tweight\n"
    assert phenotype_io.parse(loaded["content"]) == ([], [])


def test_main_phenotype_save_survives_sample_unregister(tmp_path, monkeypatch):
    sample = _legacy_sample(tmp_path, monkeypatch)
    phenotype_root = tmp_path / "patient_phenotype"
    monkeypatch.setattr(patient_phenotype_store, "PHENOTYPE_DIR", phenotype_root)
    (sample / "sample_metadata.json").write_text(
        json.dumps({"sample_id": "S1", "lis_id": "S1", "mrn": "MRN1"}),
        encoding="utf-8",
    )
    monkeypatch.setattr(sample_loader, "update_case_table_row", lambda _sample_id: None)
    monkeypatch.setattr(sample_loader, "remove_case_table_row", lambda _sample_id: None)
    monkeypatch.setattr(sample_loader, "invalidate_sample_cache", lambda _path: None)

    patient_phenotype_store.save(
        mrn="MRN1",
        code="S1",
        hpo=[{"phenotype": "HP:0001250", "label": "Seizure", "weight": 2}],
        panels=[{"name": "Neuro", "weight": 1}],
    )
    analysis_dir = sample / "analyses" / "default"
    analysis_dir.mkdir(parents=True)
    (analysis_dir / "analysis.json").write_text(
        json.dumps({"hpo": [{"phenotype": "HP:0001250"}]}),
        encoding="utf-8",
    )

    assert (phenotype_root / "MRN1_phenotype.txt").is_file()

    patient_store.delete("S1")
    assert (phenotype_root / "MRN1_phenotype.txt").is_file()
    assert not (sample / "analyses").exists()


def test_unregistered_new_lis_reuses_roster_mrn_snapshot(tmp_path, monkeypatch):
    _legacy_sample(tmp_path, monkeypatch, sample_id="NEWLIS")
    phenotype_root = tmp_path / "patient_phenotype"
    monkeypatch.setattr(patient_phenotype_store, "PHENOTYPE_DIR", phenotype_root)
    patient_phenotype_store.save(
        mrn="MRN1",
        hpo=[{"phenotype": "HP:0001250", "label": "Seizure", "weight": 1}],
        panels=[{"name": "Neuro", "weight": 1}],
    )
    monkeypatch.setattr(patient_list_store, "load_roster", lambda: {})
    monkeypatch.setattr(
        patient_list_store,
        "lookup_with_key",
        lambda lis_id, source_sample_id, roster=None: (
            {"mrn": "MRN1", "name": "Patient"}, lis_id,
        ),
    )
    monkeypatch.setattr(
        patient_list_store,
        "lookup_candidates",
        lambda lis_id, roster_lis_id="", source_sample_id="": [lis_id],
    )
    from app.services import dragen_jobs
    monkeypatch.setattr(dragen_jobs, "active_sample_ids", lambda: set())

    rows = sample_loader.list_unregistered()

    assert len(rows) == 1
    assert rows[0]["phenotype"]["path"].endswith("MRN1_phenotype.txt")
    assert rows[0]["phenotype"]["hpo"][0]["phenotype"] == "HP:0001250"
    assert rows[0]["phenotype"]["panels"][0]["name"] == "Neuro"
