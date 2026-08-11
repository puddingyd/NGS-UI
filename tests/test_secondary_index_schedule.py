from pathlib import Path

from app.services import secondary_analysis
from app.workers import secondary_index_update


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_scheduled_worker_refreshes_the_secondary_index(monkeypatch, capsys):
    monkeypatch.setattr(
        secondary_index_update.secondary_analysis,
        "refresh_index",
        lambda: {
            "updated_at": "2026-08-11 02:00:00",
            "scan_duration_sec": 8.5,
            "wes": [{}, {}],
            "wgs": [{"lane_count": 3}, {}],
        },
    )

    assert secondary_index_update.main() == 0
    output = capsys.readouterr().out
    assert "wes=2" in output
    assert "wgs=2" in output
    assert "wgs_lanes=4" in output
    assert "updated_at=2026-08-11 02:00:00" in output


def test_refresh_index_persists_results_under_a_process_lock(tmp_path, monkeypatch):
    index_path = tmp_path / "secondary_fastq_index.json"
    monkeypatch.setattr(secondary_analysis, "SECONDARY_FASTQ_INDEX_PATH", index_path)
    monkeypatch.setattr(secondary_analysis, "list_wes_fastqs", lambda: [{"sample_id": "W1"}])
    monkeypatch.setattr(secondary_analysis, "list_wgs_fastqs", lambda: [{"sample_id": "G1"}])

    idx = secondary_analysis.refresh_index()

    assert idx["wes"] == [{"sample_id": "W1"}]
    assert idx["wgs"] == [{"sample_id": "G1"}]
    assert secondary_analysis.load_index() == idx
    assert (tmp_path / "secondary_fastq_index.refresh.lock").is_file()


def test_systemd_timer_runs_daily_at_0200_taipei():
    timer = (
        REPO_ROOT / "deploy/systemd/ngs-ui-secondary-index-update.timer"
    ).read_text(encoding="utf-8")
    service = (
        REPO_ROOT / "deploy/systemd/ngs-ui-secondary-index-update.service"
    ).read_text(encoding="utf-8")

    assert "OnCalendar=*-*-* 02:00:00 Asia/Taipei" in timer
    assert "Persistent=true" in timer
    assert "Unit=ngs-ui-secondary-index-update.service" in timer
    assert (
        "ExecStart=/usr/bin/env python3 -m app.workers.secondary_index_update"
        in service
    )
