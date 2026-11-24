import json

import pytest

from app import config
from app.services import (
    patient_list_store,
    patient_phenotype_store,
    patient_store,
    phenotype_io,
    sample_loader,
    sample_layout,
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


@pytest.fixture
def recovery_tree(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "DATA_ROOT", tmp_path / "data")
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", tmp_path / "old_ui")
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", tmp_path / "pipeline")
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", tmp_path / "old_pipeline")
    monkeypatch.setattr(patient_phenotype_store, "PHENOTYPE_DIR", tmp_path / "phenotype")

    def registered(sid, mrn="MRN1", *, hpo=None, panels=None, updated="2026-11-24T00:00:00Z", legacy=False, v3=False):
        sample = (tmp_path / "pipeline" / sid / "08_postprocessing" if v3
                  else tmp_path / "old_ui" / sid)
        sample.mkdir(parents=True)
        prefix = f"{sid}." if v3 else ""
        meta = {"sample_id": sid, "lis_id": sid, "mrn": mrn, "active_analysis": "experiment"}
        snapshot = {"hpo": hpo or [], "selected_panels": panels or [], "updated_at": updated}
        if legacy:
            meta.update(snapshot)
        else:
            default = sample / "analyses" / "default"
            default.mkdir(parents=True)
            (default / f"{prefix}analysis.json").write_text(json.dumps(snapshot), encoding="utf-8")
        if v3:
            (sample / f"{sid}.layout.json").write_text(json.dumps({"layout_version": 3}), encoding="utf-8")
        (sample / f"{prefix}sample_metadata.json").write_text(json.dumps(meta), encoding="utf-8")
        return sample

    return registered


@pytest.mark.parametrize("v3", [False, True])
def test_missing_mrn_file_recovers_saved_default_via_public_load(recovery_tree, v3):
    from fastapi import FastAPI
    from fastapi.testclient import TestClient
    from app.routers.phenotype_tool import router

    sample = recovery_tree("S1", v3=v3,
        hpo=[{"phenotype": "HP:0001250", "label": "Seizure", "weight": 3}],
        panels=[{"name": "Neuro", "weight": 2}, "Custom_panel"])
    experiment = sample / "analyses" / "experiment"
    experiment.mkdir()
    (experiment / "analysis.json").write_text(json.dumps({"hpo": [{"phenotype": "HP:OTHER"}]}))
    before = {path: path.read_bytes() for path in sample.rglob("*.json")}
    app = FastAPI()
    app.include_router(router)  # The standalone phenotype load is intentionally public.
    with TestClient(app) as client:
        response = client.get("/api/phenotype-tool/load", params={"mrn": "MRN1"})
        assert response.status_code == 200
        loaded = response.json()
        assert loaded["filename"] == "MRN1_phenotype.txt"
        assert loaded["mrn"] == "MRN1"
        assert loaded["hpo"] == [{"phenotype": "HP:0001250", "label": "Seizure", "weight": 3.0}]
        assert loaded["panels"] == [{"name": "Neuro", "weight": 2.0}, {"name": "Custom_panel", "weight": 1.0}]
        assert "HP:OTHER" not in loaded["content"]
        assert client.get("/api/phenotype-tool/load", params={"mrn": "UNKNOWN"}).status_code == 404
    assert {path: path.read_bytes() for path in sample.rglob("*.json")} == before
    assert not list(sample.rglob("pheno_score.tsv"))
    canonical = patient_phenotype_store.PHENOTYPE_DIR / "MRN1_phenotype.txt"
    assert canonical.read_text(encoding="utf-8") == loaded["content"]
    assert patient_phenotype_store.backfill_registered() == {"created": 0, "errors": 0}


def test_batch_recovery_handles_legacy_metadata_and_keeps_patients_separate(recovery_tree):
    recovery_tree("S1", "MRN1", legacy=True, hpo=["HP:0001250"], panels=["PanelA"])
    recovery_tree("S2", "MRN2", hpo=[{"hpo_id": "HP:0001263", "hpo_name": "Delay"}], panels=["PanelB"])
    result = patient_phenotype_store.backfill_registered()
    assert result == {"created": 2, "errors": 0}
    one = patient_phenotype_store.load(mrn="MRN1")
    two = patient_phenotype_store.load(mrn="MRN2")
    assert one["hpo"][0]["phenotype"] == "HP:0001250"
    assert one["panels"][0]["name"] == "PanelA"
    assert two["hpo"][0]["phenotype"] == "HP:0001263"
    assert two["panels"][0]["name"] == "PanelB"


def test_recovery_uses_latest_default_including_explicit_empty(recovery_tree):
    recovery_tree("OLD", hpo=[{"phenotype": "HP:OLD"}], updated="2026-11-23T12:00:00+08:00")
    recovery_tree("NEW", panels=["NewPanel"], updated="2026-11-23T05:00:00Z")
    assert patient_phenotype_store.load(mrn="MRN1")["panels"][0]["name"] == "NewPanel"
    (patient_phenotype_store.PHENOTYPE_DIR / "MRN1_phenotype.txt").unlink()
    recovery_tree("LATEST", hpo=[], panels=[], updated="2026-11-24T00:00:00Z")
    assert patient_phenotype_store.load(mrn="MRN1")["content"] == "phenotype\thpo_name\tweight\n"


@pytest.mark.parametrize("legacy_file", [False, True])
def test_recovery_never_overwrites_existing_patient_or_empty_legacy_file(recovery_tree, legacy_file):
    recovery_tree("S1", panels=["SavedDefault"])
    root = patient_phenotype_store.PHENOTYPE_DIR
    root.mkdir()
    path = root / ("S1_MRN1_phenotype.txt" if legacy_file else "MRN1_phenotype.txt")
    body = "phenotype\thpo_name\tweight\n"
    path.write_text(body)
    assert patient_phenotype_store.backfill_registered()["created"] == 0
    assert patient_phenotype_store.load(mrn="MRN1")["content"] == body
    assert path.read_text() == body


def test_recovery_cannot_replace_a_concurrent_manual_save(recovery_tree, monkeypatch):
    recovery_tree("S1", panels=["OldDefault"])
    original = patient_phenotype_store._atomic_write

    def manual_save_first(hpo, panels, out, *, overwrite=True):
        assert overwrite is False
        patient_phenotype_store.save(mrn="MRN1", panels=["ManualPanel"])
        return original(hpo, panels, out, overwrite=overwrite)

    # Inject the competing save just before the create-only publication.
    def inject(hpo, panels, out, *, overwrite=True):
        if overwrite:
            return original(hpo, panels, out)
        return manual_save_first(hpo, panels, out, overwrite=False)

    monkeypatch.setattr(patient_phenotype_store, "_atomic_write", inject)
    loaded = patient_phenotype_store.load(mrn="MRN1")
    assert loaded["panels"][0]["name"] == "ManualPanel"


def test_invalid_default_blocks_older_specimen_without_changing_analysis(recovery_tree):
    recovery_tree("OLD", panels=["OldPanel"])
    new = recovery_tree("NEW", panels=["NewPanel"])
    (new / "analyses" / "default" / "analysis.json").write_text("{invalid")
    result = patient_phenotype_store.backfill_registered()
    assert result == {"created": 0, "errors": 1}
    assert not (patient_phenotype_store.PHENOTYPE_DIR / "MRN1_phenotype.txt").exists()


def test_no_saved_default_does_not_import_an_experimental_version(recovery_tree):
    sample = recovery_tree("S1", panels=["ExperimentOnly"])
    (sample / "analyses" / "default").rename(sample / "analyses" / "experiment")
    assert patient_phenotype_store.load(mrn="MRN1") == {}
    assert not (patient_phenotype_store.PHENOTYPE_DIR / "MRN1_phenotype.txt").exists()


def test_recovery_write_failure_returns_error_and_cleans_temporary_file(recovery_tree, monkeypatch):
    from fastapi import FastAPI
    from fastapi.testclient import TestClient
    from app.routers.phenotype_tool import router

    recovery_tree("S1", panels=["SavedPanel"])
    def fail_link(*args):
        raise OSError("synthetic storage failure")
    monkeypatch.setattr(patient_phenotype_store.os, "link", fail_link)
    app = FastAPI()
    app.include_router(router)
    with TestClient(app) as client:
        response = client.get("/api/phenotype-tool/load", params={"mrn": "MRN1"})
    assert response.status_code == 500
    assert "補建或讀取失敗" in response.json()["detail"]
    assert not list(patient_phenotype_store.PHENOTYPE_DIR.iterdir())
