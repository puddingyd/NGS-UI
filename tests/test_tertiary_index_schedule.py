from pathlib import Path

from app.services import dragen_jobs, file_index
from app.workers import tertiary_index_update

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_legacy_worker_uses_shared_daily_update(monkeypatch, capsys):
    calls = []
    def refresh(**kwargs):
        calls.append(kwargs)
        return {"meta":{"sample_count":3,"generation":"shared","updated_at":"2026-11-06 02:00:00"}}
    monkeypatch.setattr(file_index, "refresh", refresh)
    assert tertiary_index_update.main() == 0
    assert calls[0]["force"] is True
    assert "shared file index" in capsys.readouterr().out


def test_legacy_refresh_returns_shared_view_without_writing_old_cache(tmp_path, monkeypatch):
    old = tmp_path / "old-index.json"
    monkeypatch.setattr(dragen_jobs, "PIPELINE_VCF_INDEX_PATH", old)
    expected = {'dragen': [{'sample_id': 'D1'}], 'inhouse': [{'sample_id': 'N1'}]}
    monkeypatch.setattr(file_index, "refresh", lambda:{"tertiary":expected})
    assert dragen_jobs.refresh_index() == expected
    assert not old.exists()


def test_installer_migrates_to_single_taipei_timer():
    installer=(REPO_ROOT/"scripts/install_tertiary_index_timer.sh").read_text()
    assert 'install_file_index_timer.sh' in installer
    timer=(REPO_ROOT/"deploy/systemd/ngs-ui-file-index-update.timer").read_text()
    assert "OnCalendar=*-*-* 02:00:00 Asia/Taipei" in timer
    assert "Persistent=true" in timer
    assert "Unit=ngs-ui-file-index-update.service" in timer
