import json
import threading
import time
from pathlib import Path

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from app import config
from app.auth import current_user
from app.routers import samples
from app.services import dragen_jobs, patient_list_store, patient_phenotype_store, sample_layout, unregistered_samples as store


@pytest.fixture
def roots(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "DATA_ROOT", tmp_path / "data")
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", tmp_path / "new")
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", tmp_path / "old_ui")
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", tmp_path / "old_pipeline")
    monkeypatch.setattr(patient_list_store, "PATIENT_LIST_DIR", tmp_path / "roster")
    monkeypatch.setattr(patient_phenotype_store, "PHENOTYPE_DIR", tmp_path / "phenotype")
    monkeypatch.setattr(dragen_jobs, "active_sample_ids", lambda: set())
    monkeypatch.setattr(store, "_last_error_at", 0)
    return tmp_path


def sample(sample_id="S1", *, unified=True):
    if unified:
        root = sample_layout.unified_sample_dir(sample_id)
        raw = root / "03_acmg" / f"{sample_id}.snv_indel.acmg.tsv"
        raw.parent.mkdir(parents=True)
        raw.write_text("CHROM\tPOS\tREF\tALT\n")
        sample_layout.write_layout_marker(sample_id, source_id=sample_id, raw_tsv=raw)
    else:
        root = sample_layout.legacy_ui_root() / sample_id
        root.mkdir(parents=True)
        (root / "snv_indel.annotated.tsv").write_text("CHROM\tPOS\tREF\tALT\n")
    return sample_layout.state_dir(sample_id)


def test_warm_list_never_scans_or_reads_phenotype_or_stats_source_vcf(roots, monkeypatch):
    post = sample()
    sample("LEGACY", unified=False)
    (post / "S1.pipeline_source.json").write_text(json.dumps({"source_vcf_path": "/must/not/stat/source.vcf.gz"}))
    patient_list_store.link_patient(code="S1", mrn="MRN1")
    monkeypatch.setattr(patient_phenotype_store, "load", lambda **_: pytest.fail("list must not read phenotype"))
    original_stat = Path.stat

    def stat(path, *args, **kwargs):
        assert str(path) != "/must/not/stat/source.vcf.gz"
        return original_stat(path, *args, **kwargs)

    monkeypatch.setattr(Path, "stat", stat)
    first = store.listing(force=True)
    assert {r["lis_id"] for r in first["items"]} == {"S1", "LEGACY"}
    assert all("phenotype" not in r and "source_vcf_size" not in r for r in first["items"])
    monkeypatch.setattr(sample_layout, "iter_sample_ids", lambda: pytest.fail("warm list must not scan"))
    monkeypatch.setattr(store, "describe", lambda *_: pytest.fail("warm list must not read sample files"))
    second = store.listing()
    assert second["items"] == first["items"]
    assert second["refreshing"] is False
    # A new roster link updates the response without invalidating/scanning pipeline.
    patient_list_store.link_patient(code="LEGACY", mrn="MRN2")
    assert next(r for r in store.listing()["items"] if r["lis_id"] == "LEGACY")["roster"]["mrn"] == "MRN2"
    assert "mrn" not in store._path().read_text()  # demographics aren't persisted twice


def test_active_jobs_are_excluded_even_when_the_index_is_fresh(roots, monkeypatch):
    sample()
    store.refresh(force=True)
    monkeypatch.setattr(dragen_jobs, "active_sample_ids", lambda: {"S1"})
    assert store.listing()["items"] == []
    with pytest.raises(RuntimeError):
        store.detail("S1")


def test_incremental_updates_and_force_refresh_reconcile_file_changes(roots, monkeypatch):
    post = sample()
    first = store.refresh(force=True)
    meta = post / "S1.sample_metadata.json"
    meta.write_text("{}")
    store.refresh_samples(["S1"])
    assert store.listing()["items"] == []
    assert store._snapshot()["scanned_at"] == first["scanned_at"]
    meta.unlink()
    store.refresh_samples(["S1"])
    assert [r["lis_id"] for r in store.listing()["items"]] == ["S1"]
    sample("S2")
    assert len(store.listing()["items"]) == 1
    assert len(store.listing(force=True)["items"]) == 2
    # Files on disk, not the cache, decide whether a stale selection is allowed.
    meta.write_text("{}")
    with pytest.raises(FileNotFoundError):
        store.detail("S1")


def test_prefixed_layout_precedence_and_incomplete_samples(roots):
    post = sample()
    (post / "sample_metadata.json").write_text("{}")
    sample("OLD", unified=False)
    raw = sample_layout.unified_sample_dir("INCOMPLETE") / "03_acmg" / "INCOMPLETE.snv_indel.acmg.tsv"
    raw.parent.mkdir(parents=True)
    raw.touch()
    assert [r["lis_id"] for r in store.listing(force=True)["items"]] == ["OLD"]
    (post / "sample_metadata.json").unlink()
    (post / "S1.pipeline_source.json").write_text(json.dumps({"pipeline_type": "dragen"}))
    (post / "pipeline_source.json").write_text(json.dumps({"pipeline_type": "old"}))
    assert store.describe("S1")["pipeline_type"] == "dragen"


def test_detail_loads_current_patient_and_preserves_explicit_empty_snapshot(roots):
    post = sample("S1-dragen")
    vcf = roots / "source.vcf.gz"
    vcf.write_bytes(b"vcf")
    (post / "S1-dragen.pipeline_source.json").write_text(json.dumps({"source_vcf_path": str(vcf)}))
    patient_list_store.link_patient(code="S1", mrn="MRN1")
    patient_phenotype_store.save(mrn="MRN1", hpo=[{"phenotype": "HP:0001250", "label": "Seizure", "weight": 3}])
    first = store.detail("S1-dragen")
    assert first["phenotype"]["hpo"][0]["weight"] == 3
    assert first["source_vcf_size"] == 3
    patient_phenotype_store.save(mrn="MRN1", hpo=[])
    assert store.detail("S1-dragen")["phenotype"]["hpo"] == []
    patient_phenotype_store.save(mrn="MRN2", panels=[{"name": "Neuro", "weight": 2}])
    assert store.detail("S1-dragen", mrn="MRN2")["phenotype"]["panels"][0]["name"] == "Neuro"
    # Never fall back to this sample's MRN1 data when MRN3 was explicitly entered.
    assert store.detail("S1-dragen", mrn="MRN3")["phenotype"] is None


def test_legacy_filename_mrn_still_loads_on_selection(roots):
    sample("S1", unified=False)
    root = roots / "phenotype"
    root.mkdir()
    (root / "S1_MRN1_phenotype.txt").write_text("phenotype\thpo_name\tweight\nHP:0001250\tSeizure\t1\n")
    assert store.detail("S1")["phenotype"]["mrn"] == "MRN1"


def test_stale_snapshot_is_served_while_background_refresh_is_scheduled(roots, monkeypatch):
    sample()
    store.refresh(force=True)
    data = store._snapshot()
    data["scanned_at"] = time.time() - store.MAX_AGE_SECONDS - 1
    store._path().write_text(json.dumps(data))
    scheduled = []
    monkeypatch.setattr(store, "start_refresh", lambda: scheduled.append(True) or True)
    result = store.listing()
    assert result["refreshing"] is True
    assert result["items"][0]["lis_id"] == "S1"
    assert scheduled == [True]
    monkeypatch.setattr(store, "start_refresh", lambda: False)
    failed = store.listing()
    assert failed["items"] == result["items"]
    assert failed["error"]


def test_cold_cache_and_background_single_flight(roots, monkeypatch):
    entered = threading.Event()
    finish = threading.Event()
    complete = threading.Event()
    original = store.refresh
    calls = []

    def slow_refresh(**kwargs):
        calls.append(True)
        entered.set()
        assert finish.wait(5)
        try:
            return original(**kwargs)
        finally:
            complete.set()

    monkeypatch.setattr(store, "refresh", slow_refresh)
    result = store.listing()
    assert entered.wait(2)
    assert result["items"] == [] and result["refreshing"]
    assert store.start_refresh() is True
    assert calls == [True]
    finish.set()
    assert complete.wait(5)
    # Wait for the actual thread's finally rather than racing the next test.
    with store._REFRESH_LOCK:
        pass
    assert store._snapshot()["items"] == []


def test_authenticated_routes_return_summary_then_one_detail(roots):
    sample()
    app = FastAPI()
    app.include_router(samples.router)
    app.dependency_overrides[current_user] = lambda: {"id": 1}
    client = TestClient(app)
    result = client.get("/api/samples/unregistered?refresh=true")
    assert result.status_code == 200
    assert "phenotype" not in result.json()["items"][0]
    detail = client.get("/api/samples/unregistered/S1")
    assert detail.status_code == 200 and "phenotype" in detail.json()
    assert client.get("/api/samples/unregistered/UNKNOWN").status_code == 404


def test_register_rechecks_active_jobs_instead_of_trusting_cached_list(roots, monkeypatch):
    from app.services import patient_store
    sample()
    assert store.listing(force=True)["items"]
    monkeypatch.setattr(dragen_jobs, "active_sample_ids", lambda: {"S1"})
    with pytest.raises(RuntimeError, match="三級分析"):
        patient_store.register(lis_id="S1", name="Test", mrn="MRN1", hpo=[], panels=[])
    assert not sample_layout.state_file("S1", "sample_metadata.json").exists()


def test_unregister_updates_persisted_index_without_full_rescan(roots, monkeypatch):
    from app.services import patient_store, sample_loader
    post = sample()
    (post / "S1.sample_metadata.json").write_text("{}")
    assert store.listing(force=True)["items"] == []
    monkeypatch.setattr(sample_layout, "iter_sample_ids", lambda: pytest.fail("must update one row"))
    monkeypatch.setattr(sample_loader, "remove_case_table_row", lambda *_: None)
    patient_store.delete("S1")
    assert store.listing()["items"][0]["lis_id"] == "S1"
