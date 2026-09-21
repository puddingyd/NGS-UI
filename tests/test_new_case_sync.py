import json

from fastapi import FastAPI
from fastapi.testclient import TestClient

from app import config
from app.auth import current_user
from app.routers import phenotype_tool, samples
from app.services import analyses_store, clinical_presentation_store, patient_list_store, patient_phenotype_store, patient_store


def test_mrn_sync_reads_saved_content_and_registration_keeps_clinical_text(tmp_path, monkeypatch):
    monkeypatch.setattr(patient_phenotype_store, "PHENOTYPE_DIR", tmp_path / "phenotype")
    monkeypatch.setattr(clinical_presentation_store, "PHENOTYPE_DIR", tmp_path / "phenotype")
    monkeypatch.setattr(patient_list_store, "PATIENT_LIST_DIR", tmp_path / "roster")
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", tmp_path / "pipeline")
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", tmp_path / "old_pipeline")
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", tmp_path / "samples")
    captured = {}

    def register(**kwargs):
        captured.update(kwargs)
        return {"sample_id": kwargs["lis_id"]}

    monkeypatch.setattr(patient_store, "register", register)
    monkeypatch.setattr(analyses_store, "read_version", lambda *_: {})
    app = FastAPI()
    app.include_router(phenotype_tool.router)
    app.include_router(samples.router)
    app.dependency_overrides[current_user] = lambda: {"id": 1, "username": "test"}
    client = TestClient(app)
    mrn = "00123456"
    patient_phenotype_store.save(mrn=mrn, code="OLD-SAMPLE", hpo=[
        {"phenotype": "HP:0001250", "label": "Seizure", "weight": 3},
    ], panels=[{"name": "Neuro", "weight": 2}])
    clinical_presentation_store.save(mrn=mrn, code="OLD-SAMPLE", content="Saved clinical text")

    saved = client.get(f"/api/phenotype-tool/load?mrn={mrn}")
    assert saved.status_code == 200
    phenotype = saved.json()
    assert phenotype["hpo"][0]["weight"] == 3
    assert phenotype["panels"][0]["name"] == "Neuro"
    clinical = client.get(f"/api/phenotype-tool/clinical-presentation/load?mrn={mrn}")
    assert clinical.status_code == 200
    assert clinical.json()["content"].strip() == "Saved clinical text"

    response = client.post("/api/samples", data={
        "lis_id": "NEW-SAMPLE", "name": "Test", "mrn": mrn,
        "hpo_json": json.dumps(phenotype["hpo"]), "panels_json": json.dumps(phenotype["panels"]),
        "phenotype_explicit": "true",
    })
    assert response.status_code == 200
    assert captured["hpo"] == phenotype["hpo"]
    assert captured["panels"] == phenotype["panels"]
    assert captured["clinical_description"] == "Saved clinical text"
    assert client.get("/api/phenotype-tool/load?mrn=OTHER").status_code == 404
    assert client.get("/api/phenotype-tool/clinical-presentation/load?mrn=OTHER").status_code == 404
