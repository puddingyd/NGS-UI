import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from app import config
from app.routers import phenotype_tool
from app.services import (
    clinical_presentation_store,
    dragen_jobs,
    patient_list_store,
    patient_phenotype_store,
    sample_loader,
)


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setattr(patient_list_store, "PATIENT_LIST_DIR", tmp_path / "patient_list")
    monkeypatch.setattr(patient_phenotype_store, "PHENOTYPE_DIR", tmp_path / "phenotype")
    monkeypatch.setattr(clinical_presentation_store, "PHENOTYPE_DIR", tmp_path / "phenotype")
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", tmp_path / "samples")
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", tmp_path / "pipeline")
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", tmp_path / "old_pipeline")
    monkeypatch.setattr(dragen_jobs, "active_sample_ids", lambda: set())
    sample = tmp_path / "samples" / "26WE0092-dragen"
    sample.mkdir(parents=True)
    (sample / "snv_indel.annotated.tsv").write_text("CHROM\tPOS\tREF\tALT\n")
    app = FastAPI()
    app.include_router(phenotype_tool.router)
    return TestClient(app)


@pytest.mark.parametrize("mode", ["phenotype", "clinical", "ids_only"])
def test_explicit_save_link_feeds_new_case_roster_and_patient_snapshot(client, mode):
    ids = {"code": "8BB126WE0092", "mrn": "00123456"}
    if mode == "phenotype":
        response = client.post("/api/phenotype-tool/save", json={
            **ids, "content": "phenotype\thpo_name\tweight\nHP:0001250\tSeizure\t1\n",
        })
        assert response.status_code == 200
    elif mode == "clinical":
        response = client.post("/api/phenotype-tool/clinical-presentation/save", json={
            **ids, "content": "Clinical test text",
        })
        assert response.status_code == 200
    # Autosave alone never links a partially typed identity.
    assert patient_list_store.load_roster() == {}

    response = client.post("/api/phenotype-tool/patient-link", json=ids)
    assert response.status_code == 200
    rows = sample_loader.list_unregistered()
    assert len(rows) == 1
    assert rows[0]["lis_id"] == "26WE0092-dragen"
    assert rows[0]["roster"]["mrn"] == "00123456"
    if mode == "phenotype":
        assert rows[0]["phenotype"]["hpo"][0]["phenotype"] == "HP:0001250"
    elif mode == "clinical":
        clinical = clinical_presentation_store.load(mrn=rows[0]["roster"]["mrn"])
        assert clinical["content"] == "Clinical test text\n"


def test_link_endpoint_reports_invalid_and_conflicting_ids(client):
    assert client.post("/api/phenotype-tool/patient-link", json={"code": "S1"}).status_code == 400
    assert client.post("/api/phenotype-tool/patient-link", json={"code": "S1", "mrn": "MRN1"}).status_code == 200
    assert client.post("/api/phenotype-tool/patient-link", json={"code": "S1", "mrn": "MRN2"}).status_code == 409
    assert patient_list_store.lookup("S1")["mrn"] == "MRN1"
