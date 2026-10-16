from __future__ import annotations

import csv
from pathlib import Path

import pytest

from scripts import prepare_wes_prefix_reruns as reruns


def _pair(folder: Path, sample: str) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    for read in (1, 2):
        (folder / f"{sample}_S1_R{read}_001.fastq.gz").touch()


def _configure(monkeypatch, tmp_path: Path) -> Path:
    raw = tmp_path / "raw" / "20261016_run" / "Analysis" / "1" / "Data" / "fastq"
    secondary = reruns.secondary
    monkeypatch.setattr(secondary, "SECONDARY_WES_FASTQ_ROOTS", [tmp_path / "raw"])
    monkeypatch.setattr(secondary, "SECONDARY_WGS_FASTQ_ROOTS", [])
    monkeypatch.setattr(secondary, "SECONDARY_OUTPUT_ROOT", tmp_path / "output")
    monkeypatch.setattr(secondary, "SECONDARY_SAMPLESHEET_STAGING_ROOT", tmp_path / "staging")
    monkeypatch.setattr(secondary, "SECONDARY_DGX_SAMPLESHEET_STAGING_ROOT", Path("/dgx/staging"))
    return raw


def test_prepares_three_sheets_and_current_wes_launch(monkeypatch, tmp_path):
    raw = _configure(monkeypatch, tmp_path)
    for sample in ("23WE0001", "24WE0002", "25WE0003", "22WE0004"):
        _pair(raw, sample)

    results = reruns.prepare(rescan=True)

    assert [result["batch_name"] for result in results] == list(reruns.BATCHES.values())
    for result in results:
        with Path(result["samplesheet_path"]).open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == 1
        assert list(rows[0]) == ["sample", "fastq_1", "fastq_2", "sex"]
        assert rows[0]["sex"] == "unknown"
        assert rows[0]["sample"].startswith(result["batch_name"].split("_")[1])
        assert result["dgx_launch_file"] == f"/dgx/staging/{result['batch_name']}/run_on_dgx2.sh"
        command = Path(result["launch_file"]).read_text(encoding="utf-8")
        assert "-profile dgx" in command
        assert "--run_gcnv true" in command
        assert "--run_manta" in command
        assert "--run_expansionhunter" in command
        assert "--run_automap" in command
        assert "secondary_qc_report.py" in command


def test_duplicate_sample_from_different_runs_stops_all_batches(monkeypatch, tmp_path):
    raw = _configure(monkeypatch, tmp_path)
    for sample in ("23WE0001", "24WE0002", "25WE0003"):
        _pair(raw, sample)
    _pair(tmp_path / "raw" / "another_run", "23WE0001")

    with pytest.raises(ValueError, match="同一 Sample ID"):
        reruns.prepare(rescan=True)
    assert not (tmp_path / "staging").exists()


def test_missing_prefix_or_existing_batch_stops_before_writes(monkeypatch, tmp_path):
    raw = _configure(monkeypatch, tmp_path)
    _pair(raw, "23WE0001")
    _pair(raw, "24WE0002")
    with pytest.raises(ValueError, match="25WE"):
        reruns.prepare(rescan=True)
    assert not (tmp_path / "staging").exists()

    _pair(raw, "25WE0003")
    (tmp_path / "output" / "WES_24WE_re").mkdir(parents=True)
    with pytest.raises(FileExistsError, match="WES_24WE_re"):
        reruns.prepare(rescan=True)
    assert not (tmp_path / "staging").exists()


def test_incomplete_pair_stops_all_batches(monkeypatch, tmp_path):
    raw = _configure(monkeypatch, tmp_path)
    for sample in ("23WE0001", "24WE0002", "25WE0003"):
        _pair(raw, sample)
    (raw / "23WE0999_S2_R1_001.fastq.gz").touch()

    with pytest.raises(ValueError, match="23WE0999"):
        reruns.prepare(rescan=True)
    assert not (tmp_path / "staging").exists()


def test_fresh_index_avoids_full_scan(monkeypatch, tmp_path):
    raw = _configure(monkeypatch, tmp_path)
    for sample in ("23WE0001", "24WE0002", "25WE0003"):
        _pair(raw, sample)
    secondary = reruns.secondary
    index_path = tmp_path / "index.json"
    monkeypatch.setattr(secondary, "SECONDARY_FASTQ_INDEX_PATH", index_path)
    secondary.save_index({
        "schema_version": secondary._FASTQ_INDEX_SCHEMA_VERSION,
        "updated_at": secondary._now(),
        "wes": secondary.list_wes_fastqs(),
        "wgs": [],
    })
    monkeypatch.setattr(secondary, "list_wes_fastqs", lambda **kwargs: pytest.fail("should use index"))

    results = reruns.prepare(dry_run=True)

    assert [result["sample_count"] for result in results] == [1, 1, 1]
    assert not (tmp_path / "staging").exists()


def test_stale_index_requires_explicit_rescan(monkeypatch, tmp_path):
    _configure(monkeypatch, tmp_path).mkdir(parents=True)
    secondary = reruns.secondary
    index_path = tmp_path / "index.json"
    monkeypatch.setattr(secondary, "SECONDARY_FASTQ_INDEX_PATH", index_path)
    secondary.save_index({
        "schema_version": secondary._FASTQ_INDEX_SCHEMA_VERSION,
        "updated_at": "2020-01-01 00:00:00",
        "wes": [],
        "wgs": [],
    })

    with pytest.raises(ValueError, match="已過期"):
        reruns.prepare(dry_run=True)
