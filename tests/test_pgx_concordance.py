"""Synthetic calls only: comparison, reviewer persistence and report consistency."""
import io
import json
import re

import pytest
from docx import Document
from fastapi import FastAPI
from fastapi.testclient import TestClient
from starlette.middleware.sessions import SessionMiddleware

from app.auth import current_user
from app.adapters.pgx_tsv import load_pgx
from app.routers import samples
from app.services import docx_export, pgx_concordance as pc, sample_layout, sample_loader, specialized


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


def choose(pgx=None, source="dragen", genotype="*1/*1"):
    pgx = pgx or loaded()
    return pc.save_choice("SYNTHETIC", "CYP2C9", source, genotype,
                          pgx["concordance"]["genes"]["CYP2C9"]["fingerprint"],
                          {"id": 7, "username": "synthetic-reviewer"})["pgx"]


def test_select_restore_and_docx_share_results_without_rewriting_inputs(case):
    post, report, target, _, _ = case
    originals = report.read_bytes(), target.read_bytes()
    baseline = loaded()
    assert baseline["concordance"]["genes"]["CYP2C9"]["badges"] == ["結果不一致", "LowGQ"]
    assert baseline["report_view"]["drug_groups"]
    selected = choose(baseline)
    assert docx_export._pgx_gene_result(selected, "CYP2C9") == ("*1/*1", "Normal Metabolizer")
    assert not selected["report_view"]["drug_groups"]
    assert set(selected["recommendations_pending"]) == {"CYP2C9", "VKORC1"}
    row = next(x for x in selected["report_view"]["health_genotype_rows"] if x["test"] == "CYP2C9")
    assert (row["allele1"], row["allele2"]) == ("*1", "*1")
    doc = Document()
    docx_export._render_health_pgx_section(doc, "PGx", selected)
    stream = io.BytesIO(); doc.save(stream)
    roundtrip = Document(io.BytesIO(stream.getvalue()))
    text = "\n".join(p.text for p in roundtrip.paragraphs)
    assert "用藥建議待重新註解" in text and "未發現符合" not in text
    assert "Use a lower dose" not in text
    restored = choose(source="platform")
    assert docx_export._pgx_gene_result(restored, "CYP2C9")[0] == "*1/*33"
    assert restored["report_view"]["drug_groups"]
    review = json.loads((post / "SYNTHETIC.pgx_review.json").read_text())
    assert review["choices"] == {} and len(review["audit"]) == 2
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


def test_null_dragen_phenotype_never_inherits_original_phenotype(case):
    _, _, target, data, _ = case
    data["locusAnnotations"][0]["phenotypeDatabaseAnnotation"] = None
    target.write_text(json.dumps(data))
    pgx = choose()
    assert docx_export._pgx_gene_result(pgx, "CYP2C9") == ("*1/*1", "")


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
    assert response.json()["pgx"]["report_view"]["genotype_rows"][5]["genotype"] == "*1/*1"


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
