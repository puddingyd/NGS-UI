from pathlib import Path

import pytest

from app.services import dragen_jobs


def _configure_job_roots(tmp_path: Path, monkeypatch) -> Path:
    jobs = tmp_path / "jobs"
    monkeypatch.setattr(dragen_jobs, "TERTIARY_JOBS_DIR", jobs)
    monkeypatch.setattr(
        dragen_jobs,
        "LEGACY_DRAGEN_JOBS_DIR",
        tmp_path / "legacy-jobs",
    )
    return jobs


def _write_job(jobs: Path, *, job_id: str = "job-1", state: str = "done") -> Path:
    job_dir = jobs / job_id
    job_dir.mkdir(parents=True)
    dragen_jobs.save_state(job_id, {
        "job_id": job_id,
        "sample_id": "S1-dragen",
        "source_sample_id": "S1",
        "samples": [{
            "sample_id": "S1-dragen",
            "source_sample_id": "S1",
        }],
        "sample_count": 1,
        "state": state,
        "step": "done" if state == "done" else "nextflow",
        "created_at": "2026-09-22 10:00:00",
        "finished_at": "2026-09-22 11:00:00" if state == "done" else None,
    })
    return job_dir


def test_pipeline_log_view_returns_complete_job_log(tmp_path, monkeypatch):
    jobs = _configure_job_roots(tmp_path, monkeypatch)
    job_dir = _write_job(jobs)
    lines = [f"line-{index}" for index in range(2505)]
    (job_dir / "log.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    result = dragen_jobs.get_pipeline_output_log("S1-dragen")

    assert result["job_id"] == "job-1"
    assert result["log"].startswith("line-0\n")
    assert "line-2000\n" in result["log"]
    assert result["log"].endswith("line-2504\n")


def test_nextflow_log_snapshot_stays_bound_to_its_job(tmp_path, monkeypatch):
    jobs = _configure_job_roots(tmp_path, monkeypatch)
    _write_job(jobs)
    launch = tmp_path / "shared-launch"
    launch.mkdir()
    previous_signature = dragen_jobs.nextflow_log_signature(launch)
    live_log = launch / ".nextflow.log"
    live_log.write_text("first job full nextflow log\n", encoding="utf-8")

    snapshot = dragen_jobs.snapshot_nextflow_log(
        "job-1",
        launch,
        previous_signature,
    )
    assert snapshot == jobs / "job-1" / ".nextflow.log"

    live_log.write_text("later shared run\n", encoding="utf-8")
    path, filename = dragen_jobs.get_pipeline_nextflow_log("S1-dragen")

    assert path.read_text(encoding="utf-8") == "first job full nextflow log\n"
    assert filename == "S1-dragen.job-1.nextflow.log"


def test_unchanged_shared_nextflow_log_is_not_misattributed(
    tmp_path,
    monkeypatch,
):
    jobs = _configure_job_roots(tmp_path, monkeypatch)
    job_dir = _write_job(jobs)
    launch = tmp_path / "shared-launch"
    launch.mkdir()
    (launch / ".nextflow.log").write_text("previous run\n", encoding="utf-8")
    previous_signature = dragen_jobs.nextflow_log_signature(launch)

    snapshot = dragen_jobs.snapshot_nextflow_log(
        "job-1",
        launch,
        previous_signature,
    )

    assert snapshot is None
    assert not (job_dir / ".nextflow.log").exists()


def test_completed_job_never_falls_back_to_shared_nextflow_log(
    tmp_path,
    monkeypatch,
):
    jobs = _configure_job_roots(tmp_path, monkeypatch)
    job_dir = _write_job(jobs)
    launch = tmp_path / "shared-launch"
    launch.mkdir()
    (launch / ".nextflow.log").write_text("another job\n", encoding="utf-8")
    state = dragen_jobs.load_state("job-1") or {}
    state["nextflow_launch_dir"] = str(launch)
    dragen_jobs.save_state("job-1", state)
    assert not (job_dir / ".nextflow.log").exists()

    with pytest.raises(FileNotFoundError, match="尚未保存"):
        dragen_jobs.get_pipeline_nextflow_log("S1-dragen")
