from pathlib import Path

from app.services import dragen_jobs
from app.workers import tertiary_index_update


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_scheduled_worker_refreshes_the_tertiary_index(monkeypatch, capsys):
    monkeypatch.setattr(
        tertiary_index_update.dragen_jobs,
        "refresh_index",
        lambda: {
            "updated_at": "2026-07-31 02:00:00",
            "scan_duration_sec": 12.5,
            "dragen": [{}, {}],
            "inhouse": [{}],
        },
    )

    assert tertiary_index_update.main() == 0
    output = capsys.readouterr().out
    assert "dragen=2" in output
    assert "inhouse=1" in output
    assert "updated_at=2026-07-31 02:00:00" in output


def test_refresh_index_persists_results_under_a_process_lock(tmp_path, monkeypatch):
    index_path = tmp_path / "pipeline_vcf_index.json"
    monkeypatch.setattr(dragen_jobs, "PIPELINE_VCF_INDEX_PATH", index_path)
    monkeypatch.setattr(dragen_jobs, "list_dragen_vcfs", lambda: [{"sample_id": "D1"}])
    monkeypatch.setattr(dragen_jobs, "list_inhouse_vcfs", lambda: [{"sample_id": "N1"}])

    idx = dragen_jobs.refresh_index()

    assert idx["dragen"] == [{"sample_id": "D1"}]
    assert idx["inhouse"] == [{"sample_id": "N1"}]
    assert dragen_jobs.load_index() == idx
    assert (tmp_path / "pipeline_vcf_index.refresh.lock").is_file()


def test_systemd_timer_runs_daily_at_0200_taipei():
    timer = (
        REPO_ROOT / "deploy/systemd/ngs-ui-tertiary-index-update.timer"
    ).read_text(encoding="utf-8")
    service = (
        REPO_ROOT / "deploy/systemd/ngs-ui-tertiary-index-update.service"
    ).read_text(encoding="utf-8")

    assert "OnCalendar=*-*-* 02:00:00 Asia/Taipei" in timer
    assert "Persistent=true" in timer
    assert "Unit=ngs-ui-tertiary-index-update.service" in timer
    assert (
        "ExecStart=/usr/bin/env python3 -m app.workers.tertiary_index_update"
        in service
    )
