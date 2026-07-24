from __future__ import annotations

import csv
from pathlib import Path

from backend.app.services import secondary_analysis as secondary


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"")
    return path


def _lane_pair(folder: Path, sample: str, sample_number: int, lane: int) -> tuple[Path, Path]:
    stem = f"{sample}_S{sample_number}_L{lane:03d}"
    return (
        _touch(folder / f"{stem}_R1_001.fastq.gz"),
        _touch(folder / f"{stem}_R2_001.fastq.gz"),
    )


def test_wgs_index_groups_lanes_and_ignores_merged(monkeypatch, tmp_path):
    run = tmp_path / "20260611_LH00873_0018_A23NJJYLT4" / "fastq.gz"
    for lane in (1, 2, 3):
        _lane_pair(run, "26G00114", 1, lane)
    _touch(run / "26G00114_R1_merged.fastq.gz")
    _touch(run / "26G00114_R2_merged.fastq.gz")
    monkeypatch.setattr(secondary, "SECONDARY_WGS_FASTQ_ROOTS", [tmp_path])

    rows = secondary.list_wgs_fastqs()

    assert len(rows) == 1
    assert rows[0]["sample_id"] == "26G00114"
    assert rows[0]["lane_count"] == 3
    assert rows[0]["fastq_file_count"] == 6
    assert [lane["lane"] for lane in rows[0]["lanes"]] == ["L001", "L002", "L003"]
    assert all("merged" not in lane["fastq_1"] for lane in rows[0]["lanes"])


def test_legacy_per_lane_index_is_invalidated(monkeypatch, tmp_path):
    index_path = tmp_path / "secondary_fastq_index.json"
    index_path.write_text('{"updated_at":"2026-07-12 09:00:00","wes":[],"wgs":[{"lane":"L001"}]}')
    monkeypatch.setattr(secondary, "SECONDARY_FASTQ_INDEX_PATH", index_path)

    assert secondary.load_index() is None


def test_wgs_samplesheet_expands_group_to_lane_rows(monkeypatch, tmp_path):
    raw_root = tmp_path / "raw"
    folder = raw_root / "20260611_run" / "fastq.gz"
    pairs = [_lane_pair(folder, "26G00114", 1, lane) for lane in (1, 2)]
    staging = tmp_path / "staging"
    monkeypatch.setattr(secondary, "SECONDARY_WES_FASTQ_ROOTS", [])
    monkeypatch.setattr(secondary, "SECONDARY_WGS_FASTQ_ROOTS", [raw_root])
    monkeypatch.setattr(secondary, "SECONDARY_OUTPUT_ROOT", tmp_path / "output")
    monkeypatch.setattr(secondary, "SECONDARY_SAMPLESHEET_STAGING_ROOT", staging)
    monkeypatch.setattr(secondary, "SECONDARY_DGX_SAMPLESHEET_STAGING_ROOT", Path("/dgx/staging"))

    sample = {
        "sample_id": "26G00114",
        "source_sample_id": "26G00114",
        "run": "20260611_run",
        "input_dir": str(folder),
        "lane_count": 2,
        "lanes": [
            {"lane": f"L{idx:03d}", "fastq_1": str(pair[0]), "fastq_2": str(pair[1])}
            for idx, pair in zip((1, 2), pairs)
        ],
    }

    result = secondary.create_samplesheet("WGS", [sample], batch_name="260611_WGS")

    with Path(result["samplesheet_path"]).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert result["sample_count"] == 1
    assert result["samplesheet_row_count"] == 2
    assert [row["lane"] for row in rows] == ["L001", "L002"]
    assert {row["sample"] for row in rows} == {"26G00114"}
    assert all("merged" not in row["fastq_1"] for row in rows)


def test_wes_samplesheet_remains_single_row(monkeypatch, tmp_path):
    raw_root = tmp_path / "wes"
    f1 = _touch(raw_root / "SAMPLE001_S1_R1_001.fastq.gz")
    f2 = _touch(raw_root / "SAMPLE001_S1_R2_001.fastq.gz")
    monkeypatch.setattr(secondary, "SECONDARY_WES_FASTQ_ROOTS", [raw_root])
    monkeypatch.setattr(secondary, "SECONDARY_WGS_FASTQ_ROOTS", [])
    monkeypatch.setattr(secondary, "SECONDARY_OUTPUT_ROOT", tmp_path / "output")
    monkeypatch.setattr(secondary, "SECONDARY_SAMPLESHEET_STAGING_ROOT", tmp_path / "staging")
    monkeypatch.setattr(secondary, "SECONDARY_DGX_SAMPLESHEET_STAGING_ROOT", Path("/dgx/staging"))

    result = secondary.create_samplesheet("WES", [{
        "sample_id": "SAMPLE001",
        "fastq_1": str(f1),
        "fastq_2": str(f2),
    }], batch_name="260611_WES")

    with Path(result["samplesheet_path"]).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert result["sample_count"] == 1
    assert result["samplesheet_row_count"] == 1
    assert list(rows[0]) == ["sample", "fastq_1", "fastq_2", "sex"]


def test_wgs_launch_command_runs_extended_analysis_by_default():
    command = secondary._launch_command("260719_WGS", "WGS", has_one_sample=True)

    assert "-profile dgx_single" in command
    assert (
        '--seq_type WGS \\\n'
        '    --out_dir "${OUT_DIR}" \\\n'
        "    --run_manta \\\n"
        "    --run_expansionhunter \\\n"
        "    --run_automap \\\n"
        '    -w "${WORK_DIR}" \\\n'
        "    -resume"
    ) in command
    assert "--run_gcnv" not in command


def test_wes_launch_command_runs_gcnv_and_extended_analysis_by_default():
    command = secondary._launch_command("260719_WES", "WES", has_one_sample=True)

    assert "--run_gcnv true" in command
    assert "--run_manta" in command
    assert "--run_expansionhunter" in command
    assert "--run_automap" in command


def test_cleanup_secondary_nextflow_work_returns_guarded_dgx_command(monkeypatch):
    work_root = Path("/raid/DGM/work")
    monkeypatch.setattr(secondary, "SECONDARY_DGX_WORK_ROOT", work_root)

    result = secondary.cleanup_nf_work_command()

    assert result["path"] == "/raid/DGM/work"
    assert 'if [ ! -d "${SECONDARY_WORK_ROOT}" ]' in result["command"]
    assert "pgrep -af '[n]extflow'" in result["command"]
    assert 'find "${SECONDARY_WORK_ROOT}" -mindepth 1 -maxdepth 1 -print' in result["command"]
    assert 'read -r -p "確定刪除以上二級分析 Nextflow 暫存？[y/N] "' in result["command"]
    assert '-exec rm -rf -- {} +' in result["command"]
