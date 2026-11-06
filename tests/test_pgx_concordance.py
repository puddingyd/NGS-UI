"""Synthetic calls only: comparison, reviewer persistence and report consistency."""
import io
import copy
import json
import re
import sys
import threading
import time

import pytest
from docx import Document
from fastapi import FastAPI
from fastapi.testclient import TestClient
from starlette.middleware.sessions import SessionMiddleware

from app.auth import current_user
from app.adapters.pgx_tsv import load_pgx
from app.routers import samples
from app.services import docx_export, pgx_concordance as pc, sample_layout, sample_loader, specialized
from app.services import pgx_reannotation as annotations


def gene(label, **extra):
    return {"details": {"label": label, "source_diplotypes": [{"label": label}], **extra}}


def test_comparison_scope_matches_shared_presenter():
    assert pc.GENES == docx_export._PGX_CPIC_LEVEL_A_GENES


@pytest.mark.parametrize("left,right", [
    ("*36+*10/*36+*10", "*10+*36/*10+*36"),
    ("*1/*2", "*2/*1"), ("*1x2/*2", "*2/*1+*1"),
    ("CYP2D6*36+CYP2D6*10/CYP2D6*1", "*1/*10+*36"),
])
def test_order_is_ignored_with_dosage_preserved(left, right):
    assert pc.canonical("CYP2D6", left) == pc.canonical("CYP2D6", right)


def test_dosage_phase_grouping_and_combination_are_not_flattened():
    for left, right in [("*1/*2", "*1x2/*2"), ("*1+*2/*3", "*1/*2+*3"),
                        ("[*1+*2]/*3", "*1+*2/*3"), ("Reference", "Reference/Reference")]:
        assert pc.canonical("CYP2D6", left) != pc.canonical("CYP2D6", right)


def test_hgvs_intronic_plus_is_not_treated_as_a_duplication():
    assert pc.canonical("DPYD", "c.1905+1G>A/Reference") == pc.canonical("DPYD", "Reference/c.1905+1G>A")


def test_same_results_are_silent_and_quality_phase_are_independent():
    row = pc.compare_gene("CYP2B6", gene("*1/*6"), {"genotype": "*1/*6;*4/*9"})
    assert row["status"] == "same" and row["badges"] == ["相位未定"]
    row = pc.compare_gene("CFTR", gene("Reference/Reference"), {
        "genotype": "ivacaftornon-responsiveCFTRsequence/ivacaftornon-responsiveCFTRsequence",
        "missingVariantSites": [{"id": "synthetic:LowGQ"}],
    })
    assert row["badges"] == ["LowGQ"]
    assert pc.compare_gene("G6PD", gene("B (reference)/B (reference)"), {"genotype": "B(wildtype)/B(wildtype)"})["badges"] == []


def test_nat2_versioned_legacy_mapping_is_silent_without_generic_suffix_stripping():
    dr = {"genotype": "*7/*13;*7A/*13A;*7B/*4",
          "starAlleleDatabaseSource": "PharmGKB Database Version: Snapshot-2022.01.01"}
    row = pc.compare_gene("NAT2", gene("*4/*7", definition_version="2026-02-09"), dr)
    assert row["status"] == "same" and row["badges"] == []
    assert pc.canonical("NAT2", "*7C/*13A", legacy=True) != pc.canonical("NAT2", "*4/*7")


def test_mito_uses_existing_tsv_fallback_without_source_badge():
    payload = gene("Unknown")
    payload["diplotype"] = "Reference"
    assert pc.compare_gene("MT-RNR1", payload, {"genotype": "Reference"})["badges"] == []


def test_unphased_dpyd_variants_are_not_mistaken_for_alternative_diplotypes():
    payload = gene("c.1A>G", source_diplotypes=[{"label": "c.1A>G"}, {"label": "c.2C>T"}])
    assert pc.compare_gene("DPYD", payload, {"genotype": "c.1A>G/Reference"})["status"] == "uncomparable"


@pytest.fixture
def case(tmp_path, monkeypatch):
    from app import config
    monkeypatch.setattr(config, "DATA_ROOT", tmp_path / "data")
    monkeypatch.setattr(annotations, "_launch", lambda jid: None)
    monkeypatch.setattr(annotations, "runtime", lambda: {"jar": "synthetic-pinned-engine", "image": ""})
    def synthetic_engine(raw, choices, directory, run):
        out = copy.deepcopy(raw)
        for g, choice in choices.items():
            a = choice["genotype"].split("/")
            out["genes"][g]["sourceDiplotypes"] = [{"label": choice["genotype"],
                "allele1": {"name": a[0]}, "allele2": {"name": a[1]},
                "phenotypes": ["Normal Metabolizer"], "activityScore": "2.0"}]
        for drugs in out["drugs"].values():
            for payload in drugs.values():
                for guideline in payload["guidelines"]:
                    for row in guideline["annotations"]:
                        row["drugRecommendation"] = "Use the recommended starting dose."
                        row["lookupKey"][0]["CYP2C9"] = "Normal Metabolizer"
        return out
    monkeypatch.setattr(annotations, "reannotate", synthetic_engine)
    post = tmp_path / "post"
    post.mkdir()
    monkeypatch.setattr(sample_layout, "state_dir", lambda sid, **kw: post)
    monkeypatch.setattr(sample_layout, "state_file", lambda sid, name, **kw: post / f"{sid}.{name}")
    monkeypatch.setattr(sample_layout, "source_sample_id", lambda sid: "SYNTHETIC_SOURCE")
    vcf = tmp_path / "run/vcf.gz/SYNTHETIC_SOURCE.hard-filtered.vcf.gz"
    vcf.parent.mkdir(parents=True)
    vcf.touch()
    source = {"pipeline_type": "dragen", "source_sample_id": "SYNTHETIC_SOURCE", "source_vcf_path": str(vcf)}
    monkeypatch.setattr(sample_layout, "_pipeline_source", lambda sid: source)
    target = tmp_path / "run/other/SYNTHETIC_SOURCE/germline_seq/SYNTHETIC_SOURCE.targeted.json"
    target.parent.mkdir(parents=True)
    data = {"sampleId": "SYNTHETIC_SOURCE", "genomeBuild": "hg38", "locusAnnotations": [{
        "gene": "CYP2C9", "genotype": "*1/*1", "genotypeQuality": 0,
        "phenotypeDatabaseAnnotation": "Normal Metabolizer",
        "missingVariantSites": [{"id": "synthetic:LowGQ", "alleleIds": "*33"}],
    }]}
    target.write_text(json.dumps(data))
    report = {"pharmcatVersion": "synthetic", "dataVersion": "2026-01-01",
              "matcherMetadata": {"sampleId": "SYNTHETIC_SOURCE", "genomeBuild": "GRCh38.p14"},
              "genes": {"CYP2C9": {"sourceDiplotypes": [{"label": "*1/*33", "allele1": {"name": "*1"},
                           "allele2": {"name": "*33"}, "phenotypes": ["Intermediate Metabolizer"]}]}},
              "drugs": {"CPIC Guideline Annotation": {"synthetic-drug": {"guidelines": [{
                  "source": "CPIC", "annotations": [{"drugRecommendation": "Use a lower dose.",
                    "classification": "Strong", "dosingInformation": True,
                    "lookupKey": [{"CYP2C9": "Intermediate Metabolizer", "VKORC1": "test"}]}],
              }]}}}}
    report_path = post / "report.json"
    report_path.write_text(json.dumps(report))
    monkeypatch.setattr(sample_layout, "pharmcat_json", lambda sid: report_path)
    monkeypatch.setattr(sample_layout, "pgx_tsv", lambda sid: post / "missing.tsv")
    (post / "SYNTHETIC.sample_metadata.json").write_text("{}")
    from app.workers import dragen_run
    monkeypatch.setattr(dragen_run, "TERTIARY_JOBS_DIR", tmp_path / "jobs")
    return post, report_path, target, data, source


def loaded():
    return sample_loader.load_sample_pgx("SYNTHETIC")["pgx"]


def choose(pgx=None, source="dragen", genotype="*1/*1", complete=True):
    pgx = pgx or loaded()
    result = pc.save_choice("SYNTHETIC", "CYP2C9", source, genotype,
                           pgx["concordance"]["genes"]["CYP2C9"]["fingerprint"],
                           {"id": 7, "username": "synthetic-reviewer"})["pgx"]
    if complete and result.get("reannotation", {}).get("busy"):
        from app.workers.pgx_reannotate import Worker
        Worker(pc._review("SYNTHETIC")["annotation_job"]).run()
        return loaded()
    return result


def test_select_restore_and_docx_share_results_without_rewriting_inputs(case):
    post, report, target, _, _ = case
    originals = report.read_bytes(), target.read_bytes()
    baseline = loaded()
    assert baseline["concordance"]["genes"]["CYP2C9"]["badges"] == ["結果不一致", "LowGQ"]
    assert baseline["report_view"]["drug_groups"]
    selected = choose(baseline)
    assert docx_export._pgx_gene_result(selected, "CYP2C9") == ("*1/*1", "Normal Metabolizer")
    assert selected["report_view"]["drug_groups"]
    assert not selected["recommendations_pending"]
    assert selected["genes"]["CYP2C9"]["activity_score"] == "2.0"
    row = next(x for x in selected["report_view"]["health_genotype_rows"] if x["test"] == "CYP2C9")
    assert (row["allele1"], row["allele2"]) == ("*1", "*1")
    doc = Document()
    groups = docx_export._render_health_pgx_section(doc, "PGx", selected)
    docx_export._render_health_pgx_appendix(doc, "完整用藥建議", groups)
    stream = io.BytesIO(); doc.save(stream)
    roundtrip = Document(io.BytesIO(stream.getvalue()))
    text = "\n".join(roundtrip.element.xpath(".//w:t/text()"))
    assert "用藥建議待重新註解" not in text
    assert "recommended starting dose" in text
    assert "Use a lower dose" not in text
    restored = choose(source="platform")
    assert docx_export._pgx_gene_result(restored, "CYP2C9")[0] == "*1/*33"
    assert restored["report_view"]["drug_groups"]
    review = json.loads((post / "SYNTHETIC.pgx_review.json").read_text())
    assert review["choices"] == {} and len(review["audit"]) == 3
    assert review["audit"][0]["reviewer_username"] == "synthetic-reviewer"
    assert originals == (report.read_bytes(), target.read_bytes())


@pytest.mark.parametrize("change", ["target", "report", "pipeline"])
def test_input_changes_invalidate_choice_and_block_export_until_reconfirmed(case, change):
    post, report, target, data, source = case
    choose()
    if change == "target":
        data["softwareVersion"] = "changed"; target.write_text(json.dumps(data))
    elif change == "report":
        data = json.loads(report.read_text()); data["dataVersion"] = "changed"; report.write_text(json.dumps(data))
    else:
        source["job_id"] = "new-pipeline-run"
    stale = loaded()
    assert stale["review_blocked"] and stale["concordance"]["genes"]["CYP2C9"]["stale"]
    with pytest.raises(ValueError, match="重新確認"):
        pc.assert_exportable(stale)
    assert not choose(source="platform")["review_blocked"]


def test_missing_target_after_selection_does_not_silently_restore_export(case):
    choose()
    case[2].unlink()
    stale = loaded()
    assert stale["review_blocked"]
    assert not choose(source="platform")["review_blocked"]


def test_null_dragen_phenotype_is_recomputed_by_pharmcat(case):
    _, _, target, data, _ = case
    data["locusAnnotations"][0]["phenotypeDatabaseAnnotation"] = None
    target.write_text(json.dumps(data))
    pgx = choose()
    assert docx_export._pgx_gene_result(pgx, "CYP2C9") == ("*1/*1", "Normal Metabolizer")


def test_corrupt_review_record_blocks_export_instead_of_being_treated_as_no_override(case):
    (case[0] / "SYNTHETIC.pgx_review.json").write_text("not json")
    assert loaded()["review_blocked"]


def test_rule_version_change_requires_reconfirmation_of_previous_choice(case, monkeypatch):
    choose()
    monkeypatch.setattr(pc, "RULE_VERSION", "synthetic-next-version")
    pgx = loaded()
    assert pgx["review_blocked"] and pgx["concordance"]["genes"]["CYP2C9"]["stale"]
    assert not choose(source="platform")["review_blocked"]


def test_local_validated_snapshot_is_preferred_and_corruption_never_falls_back(case):
    post, _, target, data, _ = case
    specialized.write_bundle(result={"loci": {}}, payload=data, evidence=[], sources={}, build="hg38",
                             sample_id="SYNTHETIC", source_sample_id="SYNTHETIC_SOURCE", post_dir=post, status="complete")
    data["sampleId"] = "WRONG"; target.write_text(json.dumps(data))
    assert loaded()["concordance"]["genes"]["CYP2C9"]["status"] == "different"
    (post / "SYNTHETIC.specialized.source.json").write_text("{}")
    pgx = loaded()
    assert "checksum" in pgx["concordance"]["message"]
    assert pgx["concordance"]["genes"] == {}


def test_source_identity_build_and_duplicate_conflicts_are_not_compared(case):
    _, _, target, data, _ = case
    data["sampleId"] = "WRONG"; target.write_text(json.dumps(data))
    assert "樣本不符" in loaded()["concordance"]["message"]
    data["sampleId"] = "SYNTHETIC_SOURCE"; data["genomeBuild"] = "hg19"; target.write_text(json.dumps(data))
    assert "build 不符" in loaded()["concordance"]["message"]
    other = target.parents[3] / "vcf.gz/SYNTHETIC_SOURCE.targeted.json"
    other.write_text("{}")
    assert "不同 targeted.json" in loaded()["concordance"]["message"]


def test_active_pipeline_and_arbitrary_genotype_cannot_be_overridden(case):
    with pytest.raises(ValueError):
        choose(genotype="*999/*999")
    from app.workers import dragen_run
    locks = dragen_run._acquire_sample_locks(["SYNTHETIC_SOURCE"])
    try:
        with pytest.raises(RuntimeError, match="already running"):
            choose()
    finally:
        dragen_run._release_sample_locks(locks)


def test_choice_endpoint_requires_login_confirmation_and_fresh_server_candidate(case):
    app = FastAPI(); app.add_middleware(SessionMiddleware, secret_key="synthetic-only"); app.include_router(samples.router)
    client = TestClient(app)
    endpoint = "/api/samples/SYNTHETIC/pgx/choice"
    assert client.put(endpoint, json={}).status_code == 401
    app.dependency_overrides[current_user] = lambda: {"id": 7, "username": "test"}
    assert client.put(endpoint, json={}).status_code == 400
    body = {"gene": "CYP2C9", "source": "dragen", "genotype": "*1/*1", "confirmed": True,
            "fingerprint": loaded()["concordance"]["genes"]["CYP2C9"]["fingerprint"]}
    assert client.put(endpoint, json={**body, "fingerprint": "outdated"}).status_code == 409
    response = client.put(endpoint, json=body)
    assert response.status_code == 200
    assert response.json()["pgx"]["reannotation"]["busy"]
    assert response.json()["pgx"]["report_view"]["genotype_rows"][5]["genotype"] == "*1/*33"


def test_reannotation_endpoint_requires_login_and_only_accepts_bounded_actions(case):
    app = FastAPI(); app.add_middleware(SessionMiddleware, secret_key="synthetic-only"); app.include_router(samples.router)
    client = TestClient(app); endpoint = "/api/samples/SYNTHETIC/pgx/reannotation"
    assert client.post(endpoint, json={"action": "restore"}).status_code == 401
    app.dependency_overrides[current_user] = lambda: {"id": 7, "username": "test"}
    assert client.post(endpoint, json={"action": "execute", "command": "anything"}).status_code == 409
    assert client.post(endpoint, json={"action": "retry"}).status_code == 409
    assert client.post(endpoint, json={"action": "restore"}).status_code == 200


def test_health_download_uses_saved_alleles_and_rejects_change_before_archive(case, monkeypatch):
    from app import config
    choose()
    monkeypatch.setattr(sample_loader, "load_sample", lambda *a, **k: {"meta": {"LIS_ID": "SYNTHETIC", "Test": "WGS"}, "variants": {}})
    monkeypatch.setattr(docx_export.report_store, "load", lambda sid: {})
    archive = case[0] / "archive"; archive.mkdir()
    monkeypatch.setattr(config, "REPORT_OUTPUT_DIR", archive)
    app = FastAPI(); app.include_router(samples.router)
    app.dependency_overrides[current_user] = lambda: {"id": 7, "username": "test"}
    client = TestClient(app)
    url = "/api/samples/SYNTHETIC/health-report.docx?sections=pgx"
    response = client.get(url)
    assert response.status_code == 200
    doc = Document(io.BytesIO(response.content))
    lines = [re.sub(r"\s+", " ", p.text.strip()) for p in doc.paragraphs]
    assert "CYP2C9 *1 *1" in lines and "CYP2C9 *1 *33" not in lines
    count = len(list(archive.iterdir()))
    original = docx_export._section_methods
    def change_during_render(*args, **kwargs):
        original(*args, **kwargs)
        data = json.loads(case[1].read_text()); data["dataVersion"] = "changed-during-export"
        case[1].write_text(json.dumps(data))
    monkeypatch.setattr(docx_export, "_section_methods", change_during_render)
    response = client.get(url)
    assert response.status_code == 400 and "重新確認" in response.json()["detail"]
    assert len(list(archive.iterdir())) == count


def test_pending_keeps_complete_baseline_and_blocks_export(case):
    before = loaded()
    waiting = choose(complete=False)
    assert waiting["reannotation"]["busy"]
    assert waiting["report_view"]["drug_groups"] == before["report_view"]["drug_groups"]
    assert docx_export._pgx_gene_result(waiting, "CYP2C9") == docx_export._pgx_gene_result(before, "CYP2C9")
    assert waiting["concordance"]["genes"]["CYP2C9"]["requested_source"] == "dragen"
    assert "SYNTHETIC_SOURCE" in annotations.active_ids()
    with pytest.raises(ValueError, match="更新尚未完成"):
        pc.assert_exportable(waiting)
    with pytest.raises(ValueError, match="更新中"):
        choose(source="platform")


def test_failed_update_retains_baseline_then_retry_publishes(case, monkeypatch):
    from app.workers.pgx_reannotate import Worker
    engine = annotations.reannotate
    def fail(*args):
        raise ValueError("synthetic failure")
    monkeypatch.setattr(annotations, "reannotate", fail)
    pgx = choose()
    assert pgx["reannotation"]["status"] == "failed" and pgx["review_blocked"]
    assert docx_export._pgx_gene_result(pgx, "CYP2C9")[0] == "*1/*33"
    assert pgx["report_view"]["drug_groups"]
    monkeypatch.setattr(annotations, "reannotate", engine)
    annotations.action("SYNTHETIC", "retry", {"username": "retry-reviewer"})
    Worker(pc._review("SYNTHETIC")["annotation_job"]).run()
    assert not loaded()["review_blocked"]
    assert docx_export._pgx_gene_result(loaded(), "CYP2C9")[0] == "*1/*1"


def test_restore_cancels_pending_work_without_publishing(case):
    from app.workers.pgx_reannotate import Worker
    choose(complete=False)
    jid = pc._review("SYNTHETIC")["annotation_job"]
    annotations.action("SYNTHETIC", "restore", {"username": "restore-reviewer"})
    Worker(jid).run()
    assert annotations.job_status(jid)["status"] == "cancelled"
    assert not loaded()["review_blocked"]
    assert not pc._review("SYNTHETIC")["choices"]


def test_source_changed_during_work_cannot_publish(case, monkeypatch):
    engine = annotations.reannotate
    def change(*args):
        result = engine(*args)
        raw = json.loads(case[1].read_text()); raw["timestamp"] = "new-source"
        case[1].write_text(json.dumps(raw))
        return result
    monkeypatch.setattr(annotations, "reannotate", change)
    pgx = choose()
    assert pgx["review_blocked"] and pgx["reannotation"]["status"] == "failed"
    assert not pc._review("SYNTHETIC").get("published")


def test_cached_revision_reused_and_corruption_rejected(case):
    first = choose()
    pointer = pc._review("SYNTHETIC")["published"]
    choose(source="platform")
    second = choose(complete=False)
    assert second["reannotation"]["status"] == "complete"
    assert pc._review("SYNTHETIC")["published"] == pointer
    assert second["guideline_annotations"] == first["guideline_annotations"]
    (annotations.revision_dir("SYNTHETIC", pointer["job_id"]) / "report.json").write_text("{}")
    assert loaded()["review_blocked"]
    assert not annotations.action("SYNTHETIC", "restore", {})["pgx"]["review_blocked"]


def test_hospital_build_is_recorded_and_never_reuses_other_build_cache(case, monkeypatch):
    choose()
    official = pc._review("SYNTHETIC")["published"]
    choose(source="platform")
    monkeypatch.setattr(annotations, "runtime", lambda: {"jar": "synthetic-hospital-jar", "image": "synthetic.sif",
                                                        "sha256": annotations.HOSPITAL_JAR_SHA256})
    result = choose()
    hospital = pc._review("SYNTHETIC")["published"]
    assert hospital != official and not result["review_blocked"]
    manifest = annotations.read(annotations.revision_dir("SYNTHETIC", hospital["job_id"]) / "manifest.json")
    assert manifest["engine_sha256"] == annotations.HOSPITAL_JAR_SHA256
    choose(source="platform")
    assert choose(complete=False)["reannotation"]["status"] == "complete"
    assert pc._review("SYNTHETIC")["published"] == hospital


def test_same_choice_with_new_build_stays_blocked_until_new_revision_complete(case, monkeypatch):
    from app.workers.pgx_reannotate import Worker
    before = choose()
    old = pc._review("SYNTHETIC")["published"]
    monkeypatch.setattr(annotations, "runtime", lambda: {"jar": "synthetic-hospital-jar", "image": "synthetic.sif",
                                                        "sha256": annotations.HOSPITAL_JAR_SHA256})
    waiting = choose(complete=False)
    assert waiting["review_blocked"] and waiting["reannotation"]["busy"]
    assert waiting["guideline_annotations"] == before["guideline_annotations"]
    assert pc._review("SYNTHETIC")["published"] == old
    Worker(pc._review("SYNTHETIC")["annotation_job"]).run()
    assert not loaded()["review_blocked"]
    assert pc._review("SYNTHETIC")["published"] != old


def test_unavailable_runtime_keeps_request_and_complete_previous_result_then_recovers(case, monkeypatch):
    before = choose()
    pointer = pc._review("SYNTHETIC")["published"]
    runtime = annotations.runtime
    def unavailable():
        raise ValueError("PharmCAT 程式校驗不符")
    def must_not_launch(jid):
        pytest.fail("unregistered engine must not launch")
    monkeypatch.setattr(annotations, "runtime", unavailable)
    monkeypatch.setattr(annotations, "_launch", must_not_launch)
    failed = choose(complete=False)
    assert failed["review_blocked"] and failed["reannotation"]["status"] == "failed"
    assert pc._review("SYNTHETIC")["choices"]["CYP2C9"]["genotype"] == "*1/*1"
    assert pc._review("SYNTHETIC")["published"] == pointer
    assert failed["guideline_annotations"] == before["guideline_annotations"]
    monkeypatch.setattr(annotations, "runtime", runtime)
    recovered = annotations.action("SYNTHETIC", "retry", {})["pgx"]
    assert not recovered["review_blocked"]
    assert pc._review("SYNTHETIC")["published"] == pointer


def test_build_changed_after_queue_cannot_silently_publish(case, monkeypatch):
    from app.workers.pgx_reannotate import Worker
    choose(complete=False)
    jid = pc._review("SYNTHETIC")["annotation_job"]
    monkeypatch.setattr(annotations, "runtime", lambda: {"jar": "synthetic-hospital-jar", "image": "synthetic.sif",
                                                        "sha256": annotations.HOSPITAL_JAR_SHA256})
    Worker(jid).run()
    result = loaded()
    assert result["review_blocked"] and result["reannotation"]["status"] == "failed"
    assert not pc._review("SYNTHETIC").get("published")


def test_old_manual_choices_require_update_without_mixed_advice(case):
    choose()
    path = case[0] / "SYNTHETIC.pgx_review.json"
    review = json.loads(path.read_text()); review.pop("published"); review.pop("annotation_job")
    path.write_text(json.dumps(review))
    pgx = loaded()
    assert pgx["reannotation"]["status"] == "needs_update" and pgx["review_blocked"]
    assert docx_export._pgx_gene_result(pgx, "CYP2C9")[0] == "*1/*33"


def test_second_gene_failure_keeps_previous_published_set_then_reuses_it(case, monkeypatch):
    from app.workers.pgx_reannotate import Worker
    raw = json.loads(case[1].read_text())
    raw["genes"]["HLA-B"] = {"sourceDiplotypes": [{"label": "*40:01/*51:01",
        "allele1": {"name": "*40:01"}, "allele2": {"name": "*51:01"}, "phenotypes": ["*15:02 negative"]}]}
    case[1].write_text(json.dumps(raw))
    target = json.loads(case[2].read_text()); target["locusAnnotations"].append({"gene": "HLA-B", "genotype": "*15:02/*40:01"})
    case[2].write_text(json.dumps(target))
    published = choose()
    row = published["concordance"]["genes"]["HLA-B"]
    pending = pc.save_choice("SYNTHETIC", "HLA-B", "dragen", "*15:02/*40:01", row["fingerprint"], {})["pgx"]
    assert pending["reannotation"]["busy"] and pending["review_blocked"]
    for g in ("CYP2C9", "HLA-B"):
        assert docx_export._pgx_gene_result(pending, g) == docx_export._pgx_gene_result(published, g)
    assert pending["guideline_annotations"] == published["guideline_annotations"]
    def fail(*args):
        raise ValueError("synthetic second revision failure")
    monkeypatch.setattr(annotations, "reannotate", fail)
    Worker(pc._review("SYNTHETIC")["annotation_job"]).run()
    assert loaded()["reannotation"]["status"] == "failed"
    back = pc.save_choice("SYNTHETIC", "HLA-B", "platform", "", row["fingerprint"], {})["pgx"]
    assert not back["review_blocked"]
    assert docx_export._pgx_gene_result(back, "CYP2C9")[0] == "*1/*1"


def test_queued_pgx_blocks_both_unregister_and_output_deletion(case):
    from app.services import patient_store, dragen_jobs
    choose(complete=False)
    for remove in (patient_store.delete, dragen_jobs.delete_pipeline_output):
        with pytest.raises(RuntimeError, match="PGx 更新"):
            remove("SYNTHETIC")
    assert (case[0] / "SYNTHETIC.sample_metadata.json").exists()


def test_dead_worker_is_not_treated_as_forever_running(case):
    choose(complete=False)
    jid = pc._review("SYNTHETIC")["annotation_job"]
    path = annotations.job_dir(jid) / "state.json"
    job = annotations.read(path); job.update(status="running", created=0); annotations.atomic(path, job)
    assert not annotations.active_ids()
    pgx = loaded()
    assert pgx["review_blocked"] and not pgx["reannotation"]["busy"]
    assert pgx["reannotation"]["status"] == "failed"


def test_restore_stops_owned_tool_process_group(case, monkeypatch):
    from app.workers.pgx_reannotate import Worker
    marker = case[0] / "owned-child.ready"
    stopped = case[0] / "owned-child.stopped"
    child = ("import signal,sys,time; from pathlib import Path; "
             f"signal.signal(signal.SIGTERM,lambda *args:(Path({str(stopped)!r}).write_text('stopped'),sys.exit(0))); "
             f"Path({str(marker)!r}).write_text('ready'); time.sleep(60)")
    program = ("import subprocess,sys,time; from pathlib import Path; "
               f"p=subprocess.Popen([sys.executable,'-c',{child!r}]); time.sleep(60)")
    monkeypatch.setattr(annotations, "command", lambda *args: [sys.executable, "-c", program])
    monkeypatch.setattr(annotations, "reannotate", lambda raw, choices, directory, run: run(["synthetic"]))
    choose(complete=False)
    jid = pc._review("SYNTHETIC")["annotation_job"]
    errors = []
    def restore():
        try:
            for _ in range(100):
                if marker.exists():
                    annotations.action("SYNTHETIC", "restore", {"username": "cancel-reviewer"})
                    return
                time.sleep(.05)
            raise AssertionError("tool never started")
        except Exception as exc:
            errors.append(exc)
    thread = threading.Thread(target=restore); thread.start()
    Worker(jid).run(); thread.join(timeout=10)
    assert not errors and not thread.is_alive()
    assert annotations.job_status(jid)["status"] == "cancelled"
    assert not loaded()["review_blocked"]
    assert stopped.read_text() == "stopped"
