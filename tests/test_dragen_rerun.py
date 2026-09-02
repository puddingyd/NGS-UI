import json
import sqlite3
from pathlib import Path

import pytest

from app.services import dragen_jobs
from app.workers import dragen_run


def _make_pipeline_tree(root, sample_id, *, with_pgx=True):
    sample = root / sample_id
    for name in dragen_run.REQUIRED_PIPELINE_STAGE_DIRS:
        (sample / name).mkdir(parents=True)
    (sample / "03_acmg" / f"{sample_id}.snv_indel.acmg.tsv").write_text(
        "CHROM\tPOS\tREF\tALT\n",
        encoding="utf-8",
    )
    (sample / "04_mito" / f"{sample_id}.mito.tsv").touch()
    (sample / "05_str" / f"{sample_id}.str.tsv").touch()
    (sample / "06_cnv_sv" / f"{sample_id}.cnv.annotated.tsv").touch()
    (sample / "06_cnv_sv" / f"{sample_id}.sv.annotated.tsv").touch()
    if with_pgx:
        (sample / "07_pgx").mkdir()
        (sample / "07_pgx" / f"{sample_id}.pgx.tsv").touch()
    return sample


def _sqlite_meta(path, source_path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(path) as conn:
        conn.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
        conn.execute(
            "INSERT INTO meta(key, value) VALUES ('source_path', ?)",
            (str(source_path),),
        )


def test_nextflow_context_is_shared_across_batches_per_mode(tmp_path, monkeypatch):
    monkeypatch.setattr(dragen_run, "TERTIARY_NF_WORK_ROOT", tmp_path)
    first_batch = [{
        "mode": "dragen",
        "sample_id": "S1-dragen",
        "source_sample_id": "S1",
    }, {
        "mode": "dragen",
        "sample_id": "S2-dragen",
        "source_sample_id": "S2",
    }]
    overlapping_batch = [dict(first_batch[0])]
    unrelated_batch = [{
        "mode": "dragen",
        "sample_id": "S3-dragen",
        "source_sample_id": "S3",
    }]
    inhouse_batch = [{
        "mode": "inhouse",
        "sample_id": "S1-nckuh",
        "source_sample_id": "S1",
    }]

    first = dragen_run._nextflow_context(first_batch)

    assert first == dragen_run._nextflow_context(overlapping_batch)
    assert first == dragen_run._nextflow_context(unrelated_batch)
    assert first != dragen_run._nextflow_context(inhouse_batch)
    assert first == (
        tmp_path / "contexts" / "shared-dragen" / "launch",
        tmp_path / "contexts" / "shared-dragen" / "work",
    )


def test_shared_context_adopts_prior_session_with_most_sample_overlap(
    tmp_path,
    monkeypatch,
):
    repo = tmp_path / "repo"
    work_root = tmp_path / "nf_work"
    launch = work_root / "contexts" / "shared-dragen" / "launch"
    best_session = "11111111-2222-3333-4444-555555555555"
    other_session = "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
    best_nf = work_root / "contexts" / "old-batch" / "launch" / ".nextflow"
    other_nf = work_root / "contexts" / "other-batch" / "launch" / ".nextflow"
    best_work = work_root / "contexts" / "old-batch" / "work"
    other_work = work_root / "contexts" / "other-batch" / "work"
    for nf_root, session_id in (
        (best_nf, best_session),
        (other_nf, other_session),
    ):
        (nf_root / "cache" / session_id).mkdir(parents=True)
        (nf_root / "cache" / session_id / "index.run").touch()
    best_work.mkdir(parents=True)
    other_work.mkdir(parents=True)
    best_sheet = tmp_path / "best.csv"
    best_sheet.write_text(
        "sample_id,pipeline_type,input_dir,seq_type,hpo\n"
        "S1,dragen,/input,WGS,\n"
        "S2,dragen,/input,WGS,\n",
        encoding="utf-8",
    )
    other_sheet = tmp_path / "other.csv"
    other_sheet.write_text(
        "sample_id,pipeline_type,input_dir,seq_type,hpo\n"
        "S1,dragen,/input,WGS,\n",
        encoding="utf-8",
    )
    (best_nf / "history").write_text(
        "2026-07-01 01:00:00\t1m\tbest\tERR\thash\t"
        f"{best_session}\tnextflow run main.nf --pipeline_type dragen "
        f"--samplesheet {best_sheet} -work-dir {best_work} -resume\n",
        encoding="utf-8",
    )
    (other_nf / "history").write_text(
        "2026-07-02 01:00:00\t1m\tother\tOK\thash\t"
        f"{other_session}\tnextflow run main.nf --pipeline_type dragen "
        f"--samplesheet {other_sheet} -work-dir {other_work} -resume\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(dragen_run, "REPO_ROOT", repo)
    monkeypatch.setattr(dragen_run, "TERTIARY_NF_WORK_ROOT", work_root)

    imported = dragen_run._seed_shared_nextflow_context(
        [
            {
                "mode": "dragen",
                "sample_id": "S1-dragen",
                "source_sample_id": "S1",
            },
            {
                "mode": "dragen",
                "sample_id": "S2-dragen",
                "source_sample_id": "S2",
            },
        ],
        launch_dir=launch,
    )

    assert imported == (best_session, best_work.resolve())
    assert (
        launch / ".nextflow" / "cache" / best_session / "index.run"
    ).is_file()
    assert best_session in (launch / ".nextflow" / "history").read_text()
    dragen_run._write_shared_work_pointer(launch, imported[1])
    assert dragen_run._nextflow_context([{
        "mode": "dragen",
        "sample_id": "S1-dragen",
        "source_sample_id": "S1",
    }])[1] == best_work.resolve()
    assert dragen_run._seed_shared_nextflow_context(
        [{
            "mode": "dragen",
            "sample_id": "S1-dragen",
            "source_sample_id": "S1",
        }],
        launch_dir=launch,
    ) is None


def test_sample_lock_rejects_concurrent_rerun_and_releases(tmp_path, monkeypatch):
    monkeypatch.setattr(dragen_run, "TERTIARY_JOBS_DIR", tmp_path / "jobs")
    first = dragen_run._acquire_sample_locks(["S1-dragen", "S1"])
    try:
        with pytest.raises(RuntimeError, match="already running"):
            dragen_run._acquire_sample_locks(["S1"])
    finally:
        dragen_run._release_sample_locks(first)

    second = dragen_run._acquire_sample_locks(["S1"])
    dragen_run._release_sample_locks(second)


def test_nextflow_cache_lock_is_separate_per_mode(tmp_path, monkeypatch):
    monkeypatch.setattr(dragen_run, "TERTIARY_JOBS_DIR", tmp_path / "jobs")
    dragen_lock = dragen_run._acquire_nextflow_cache_lock(
        "dragen",
        blocking=False,
    )
    try:
        with pytest.raises(BlockingIOError):
            dragen_run._acquire_nextflow_cache_lock(
                "dragen",
                blocking=False,
            )
        nckuh_lock = dragen_run._acquire_nextflow_cache_lock(
            "inhouse",
            blocking=False,
        )
        dragen_run._release_nextflow_cache_lock(nckuh_lock)
    finally:
        dragen_run._release_nextflow_cache_lock(dragen_lock)

    again = dragen_run._acquire_nextflow_cache_lock(
        "dragen",
        blocking=False,
    )
    dragen_run._release_nextflow_cache_lock(again)


def test_nextflow_batch_progress_uses_known_sample_count_while_total_grows():
    tokens = ("SNV_ANNOTATE:VEP_ANNOTATE",)

    first = dragen_run._parse_nextflow_progress_line(
        "[aa/000001] SNV_ANNOTATE:VEP_ANNOTATE (1) | 1 of 1",
        tokens,
        expected_total=40,
    )
    halfway = dragen_run._parse_nextflow_progress_line(
        "[bb/000002] SNV_ANNOTATE:VEP_ANNOTATE (20) | 20 of 20",
        tokens,
        expected_total=40,
    )
    finished = dragen_run._parse_nextflow_progress_line(
        "[cc/000003] SNV_ANNOTATE:VEP_ANNOTATE (40) | 40 of 40 ✔",
        tokens,
        expected_total=40,
    )

    assert first is not None
    assert first["event"] == "start"
    assert first["progress_total"] == 40
    assert first["fraction"] == pytest.approx(1 / 40)
    assert halfway is not None
    assert halfway["event"] == "start"
    assert halfway["fraction"] == pytest.approx(0.5)
    assert finished is not None
    assert finished["event"] == "done"
    assert finished["fraction"] == 1.0


def test_nextflow_complete_marker_finishes_optional_shorter_process():
    parsed = dragen_run._parse_nextflow_progress_line(
        "\x1b[32m[dd/000004] PGX_ANNOTATE:PGX_PARSE (38) | 38 of 38 ✔\x1b[0m",
        ("PGX_ANNOTATE:PGX_PARSE",),
        expected_total=40,
    )

    assert parsed is not None
    assert parsed["event"] == "done"
    assert parsed["fraction"] == 1.0
    assert parsed["reported_total"] == 38
    assert parsed["progress_total"] == 40


def test_start_job_rejects_active_sample_before_spawning(tmp_path, monkeypatch):
    vcf = tmp_path / "S1.hard-filtered.vcf.gz"
    vcf.touch()
    monkeypatch.setattr(
        dragen_jobs,
        "active_sample_ids",
        lambda: {"S1", "S1-dragen"},
    )

    with pytest.raises(RuntimeError, match="already running"):
        dragen_jobs.start_job(
            str(vcf),
            "S1",
            source_sample_id="S1",
            mode="dragen",
        )


def test_staging_validation_requires_str_and_requested_pgx(tmp_path):
    _make_pipeline_tree(tmp_path, "S1", with_pgx=True)

    outputs = dragen_run._validate_staged_sample(
        tmp_path,
        "S1",
        require_pgx=True,
    )
    assert set(outputs) == {
        "snv_indel.acmg.tsv",
        "mito.tsv",
        "str.tsv",
        "cnv.annotated.tsv",
        "sv.annotated.tsv",
        "PGx/PharmCAT",
    }

    (tmp_path / "S1" / "05_str" / "S1.str.tsv").unlink()
    with pytest.raises(RuntimeError, match="str.tsv"):
        dragen_run._validate_staged_sample(
            tmp_path,
            "S1",
            require_pgx=True,
        )


def test_v36_validation_checks_research_dbnsfp_branch(tmp_path):
    required = {
        "CHROM", "POS", "REF", "ALT", "GENE", "TRANSCRIPT",
        "TRANSCRIPT_TYPE", "HGVS_C", "HGVS_P", "CONSEQUENCE", "IMPACT",
        "HGNC_ID", "ACMG_CRITERIA", "ACMG_SCORE", "ACMG_CLASS", "ACMG_NOTES",
        "STRAND_BIAS", "CLINGEN_VCEP_CLASS", "CLINGEN_VCEP_CRITERIA",
        "CLINGEN_VCEP_PANEL", "REVEL", "MUTPRED2", "MUTPRED2_PRED", "VEST4",
        "CADD_PHRED", "DBNSFP_VERSION", "CLINGEN_AGREEMENT", "PVS1_STRENGTH",
        "PVS1_REASON",
    }
    fields = sorted(required)
    fields.extend(f"DUMMY_{index}" for index in range(81 - len(fields)))
    values = ["."] * len(fields)
    values[fields.index("DBNSFP_VERSION")] = "5.3a"
    path = tmp_path / "sample.snv_indel.acmg.tsv"
    path.write_text(
        "\t".join(fields) + "\n" + "\t".join(values) + "\n",
        encoding="utf-8",
    )

    dragen_run._validate_acmg_tsv(
        path, strict_v31=True, expect_academic_dbnsfp=True
    )
    with pytest.raises(RuntimeError, match="expected 4.9c"):
        dragen_run._validate_acmg_tsv(
            path, strict_v31=True, expect_academic_dbnsfp=False
        )

    sidecar = tmp_path / "sample.annotation_versions.json"
    sidecar.write_text(
        json.dumps({"databases": {"clinvar": {"release_date": "2026-05-10"}}}),
        encoding="utf-8",
    )
    dragen_run._ensure_v36_annotation_versions(path)
    payload = json.loads(sidecar.read_text(encoding="utf-8"))
    assert payload["databases"]["clinvar"]["release_date"] == "2026-07-20"


def test_production_postprocessing_is_spliceai_only():
    script = (
        Path(__file__).resolve().parents[1] / "scripts/run_stopgaps.sh"
    ).read_text(encoding="utf-8")
    active = script[script.index("# 3. SpliceAI-only"):script.index("# 3b. GIAB")]
    assert "annotate_spliceai.py" in active
    assert "annotate_extra_vep.py" not in active
    assert "--dbnsfp" not in active
    assert '--academic_dbnsfp", "true"' in Path(dragen_run.__file__).read_text(encoding="utf-8")

def test_rebase_staged_derived_paths_uses_future_live_paths(tmp_path):
    stage_post = tmp_path / "stage" / "08_postprocessing"
    final_post = tmp_path / "live" / "08_postprocessing"
    stage_raw = tmp_path / "stage" / "03_acmg" / "SRC.snv_indel.acmg.tsv"
    final_raw = tmp_path / "live" / "03_acmg" / "SRC.snv_indel.acmg.tsv"
    stage_raw.parent.mkdir(parents=True)
    stage_raw.touch()
    overlay = stage_post / "S1.snv_annotations.sqlite"
    index = stage_post / "S1.snv_gene_index.sqlite"
    _sqlite_meta(overlay, stage_raw)
    _sqlite_meta(index, stage_raw)
    manifest = stage_post / "S1.snv_indel.review.tsv.source.json"
    manifest.write_text(
        json.dumps({
            "overlay": {
                "exists": True,
                "path": str(overlay),
                "current": True,
            },
        }),
        encoding="utf-8",
    )

    dragen_run._rebase_staged_derived_paths(
        sample_id="S1",
        stage_post_dir=stage_post,
        final_raw_tsv=final_raw,
        final_post_dir=final_post,
    )

    for path in (overlay, index):
        with sqlite3.connect(path) as conn:
            value = conn.execute(
                "SELECT value FROM meta WHERE key = 'source_path'"
            ).fetchone()[0]
        assert value == str(final_raw.absolute())
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert payload["overlay"]["path"] == str(
        (final_post / "S1.snv_annotations.sqlite").absolute()
    )
    assert payload["overlay"]["mtime_ns"] == overlay.stat().st_mtime_ns
    assert payload["overlay"]["size"] == overlay.stat().st_size


def test_generation_promotion_preserves_reviewer_state_and_can_roll_back(tmp_path):
    live = tmp_path / "live" / "S1-dragen"
    staged = _make_pipeline_tree(tmp_path / "stage", "S1", with_pgx=True)
    rollback = tmp_path / "rollback" / "S1-dragen"
    for name in dragen_run.PIPELINE_STAGE_DIRS:
        source = staged / name
        if source.is_dir():
            (source / "generation.txt").write_text("new", encoding="utf-8")
        target = live / name
        target.mkdir(parents=True)
        (target / "generation.txt").write_text("old", encoding="utf-8")

    live_post = live / "08_postprocessing"
    live_post.mkdir()
    metadata = live_post / "S1-dragen.sample_metadata.json"
    metadata.write_text('{"comment":"keep me"}', encoding="utf-8")
    stale_ploidy = live_post / "S1-dragen.ploidy.vcf.gz"
    stale_ploidy.write_text("old ploidy", encoding="utf-8")

    staged_post = staged / "08_postprocessing"
    staged_post.mkdir()
    (staged_post / "S1-dragen.pipeline_source.json").write_text(
        '{"source":"new"}',
        encoding="utf-8",
    )
    (staged_post / "S1-dragen.layout.json").write_text(
        '{"layout_version":3}',
        encoding="utf-8",
    )

    operations = dragen_run._promote_staged_generation(
        sample_id="S1-dragen",
        staged_sample_dir=staged,
        live_sample_dir=live,
        rollback_sample_dir=rollback,
    )

    assert (live / "03_acmg" / "generation.txt").read_text() == "new"
    assert json.loads(metadata.read_text())["comment"] == "keep me"
    assert not stale_ploidy.exists()
    assert json.loads(
        (live_post / "S1-dragen.pipeline_source.json").read_text()
    )["source"] == "new"

    dragen_run._rollback_promotion_operations(operations)

    assert (live / "03_acmg" / "generation.txt").read_text() == "old"
    assert stale_ploidy.read_text() == "old ploidy"
    assert json.loads(metadata.read_text())["comment"] == "keep me"


def test_generation_promotion_removes_old_pgx_when_rerun_disables_it(tmp_path):
    live = tmp_path / "live" / "S1-dragen"
    staged = _make_pipeline_tree(tmp_path / "stage", "S1", with_pgx=False)
    rollback = tmp_path / "rollback" / "S1-dragen"
    old_pgx = live / "07_pgx"
    old_pgx.mkdir(parents=True)
    (old_pgx / "old.pgx.tsv").touch()
    staged_post = staged / "08_postprocessing"
    staged_post.mkdir()
    (staged_post / "S1-dragen.layout.json").write_text(
        '{"layout_version":3}',
        encoding="utf-8",
    )

    operations = dragen_run._promote_staged_generation(
        sample_id="S1-dragen",
        staged_sample_dir=staged,
        live_sample_dir=live,
        rollback_sample_dir=rollback,
    )

    assert not old_pgx.exists()
    dragen_run._rollback_promotion_operations(operations)
    assert (old_pgx / "old.pgx.tsv").is_file()
