from __future__ import annotations

import csv
import os
import subprocess
from pathlib import Path

import pytest

from backend.app.services import secondary_analysis as secondary


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"")
    return path


def _fake_command(path: Path, body: str) -> None:
    path.write_text(f"#!/bin/sh\n{body}\n", encoding="utf-8")
    path.chmod(0o755)


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


def test_server_path_to_dgx_maps_home_raw_root_without_changing_direct_raw_paths(monkeypatch):
    monkeypatch.setattr(secondary, "SECONDARY_DGX_RAW_ROOT", Path("/datalake_Raw/datalake_Raw"))

    assert secondary._server_path_to_dgx(
        "/home/datalake_Raw/NextSeq2000/run/Analysis/1/Data/fastq/sample_R1_001.fastq.gz"
    ) == (
        "/datalake_Raw/datalake_Raw/NextSeq2000/run/Analysis/1/Data/fastq/"
        "sample_R1_001.fastq.gz"
    )
    assert secondary._server_path_to_dgx(
        "/datalake_Raw/Other/Reanalysis/sample.R1.clean.fastq.gz"
    ) == "/datalake_Raw/Other/Reanalysis/sample.R1.clean.fastq.gz"


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
    assert result["qc_report_path"].endswith("/260611_WES/pipeline_info/report_summary.csv")


def test_wgs_launch_command_runs_extended_analysis_by_default():
    command = secondary._launch_command("260719_WGS", "WGS")

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
    assert "secondary_qc_report.py" not in command


def test_wes_launch_command_runs_gcnv_and_extended_analysis_by_default():
    command = secondary._launch_command("260719_WES", "WES")

    assert "-profile dgx \\" in command
    assert "-profile dgx_single" not in command
    assert "--run_gcnv true" in command
    assert "--run_manta" in command
    assert "--run_expansionhunter" in command
    assert "--run_automap" in command


def test_multi_sample_launch_also_uses_dgx_single_profile():
    command = secondary._launch_command("260719_WGS_MULTI", "WGS")

    assert "-profile dgx_single" in command
    assert "-profile dgx \\" not in command


def test_launch_command_uses_group_writable_umask_after_environment():
    command = secondary._launch_command("260831_WGS", "WGS")

    source_position = command.index("source ")
    umask_position = command.index("umask 0002")
    mkdir_position = command.index('mkdir -p "${OUT_DIR}"')

    assert source_position < umask_position < mkdir_position
    assert 'ORIGINAL_UMASK="$(umask)"' in command
    assert 'umask "${ORIGINAL_UMASK}"\n    exec bash -i' in command


@pytest.mark.parametrize("failure,expected,stage", [
    ("", ["config", "preflight", "nextflow", "qc"], None),
    ("missing_script", [], "QC preflight"),
    ("config", ["config"], "QC preflight"),
    ("preflight", ["config", "preflight"], "QC preflight"),
    ("nextflow", ["config", "preflight", "nextflow"], "Nextflow"),
    ("qc", ["config", "preflight", "nextflow", "qc"], "QC report"),
])
def test_generated_wes_runner_stage_order_and_errors(monkeypatch, tmp_path, failure, expected, stage):
    """Execute the actual generated Bash, replacing only external programs."""
    batch = "BATCH_WES"
    for constant, folder in [
        ("SECONDARY_DGX_OUTPUT_ROOT", "output"),
        ("SECONDARY_DGX_LAUNCH_ROOT", "launch"),
        ("SECONDARY_DGX_WORK_ROOT", "work"),
        ("SECONDARY_DGX_SAMPLESHEET_STAGING_ROOT", "staging"),
    ]:
        monkeypatch.setattr(secondary, constant, tmp_path / folder)
    _touch(tmp_path / "staging" / batch / "samplesheet.csv")
    code_dir = tmp_path / "pipeline code"
    qc_script = _touch(code_dir / "scripts/secondary_qc_report.py")
    if failure == "missing_script":
        qc_script.unlink()
    env_script = tmp_path / "env.sh"
    env_script.write_text(f'export PIPELINE_CODE="{code_dir}"\nexport PIPELINE_CONFIG="{code_dir}/nextflow_main.config"\n')
    monkeypatch.setattr(secondary, "SECONDARY_DGX_ENV_SCRIPT", env_script)
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    events = tmp_path / "events"
    _fake_command(bin_dir / "nextflow", '''
kind=nextflow
for arg in "$@"; do
    if [ "$arg" = config ]; then kind=config; fi
done
echo "$kind" >> "$EVENTS"
if [ "$kind" = nextflow ]; then echo run-log > .nextflow.log; fi
if [ "$kind" = "$FAIL_STAGE" ]; then exit 7; fi
echo "params.wes_targets = '/ref/targets.bed'"
''')
    _fake_command(bin_dir / "python3", '''
kind=qc
for arg in "$@"; do
    if [ "$arg" = --check-only ]; then kind=preflight; fi
done
echo "$kind" >> "$EVENTS"
if [ "$kind" = "$FAIL_STAGE" ]; then exit 9; fi
echo "QC: FAIL=1 ERROR=0"
''')
    # The production trap deliberately keeps tmux open. End that final shell in the test.
    _fake_command(bin_dir / "bash", "exit 0")
    launch_dir = tmp_path / "launch" / batch
    launch_dir.mkdir(parents=True)
    (launch_dir / ".nextflow.log").write_text("stale-log")
    generated = secondary._launch_command(batch, "WES")
    runner = generated.split("<<'NGS2_EOF'\n", 1)[1].split("\nNGS2_EOF", 1)[0]
    result = subprocess.run(["/bin/bash", "-c", runner], capture_output=True, text=True,
                            env={**os.environ, "PATH": f"{bin_dir}:{os.environ['PATH']}",
                                 "EVENTS": str(events), "FAIL_STAGE": failure}, timeout=20)
    assert result.returncode == 0  # Only the test's final interactive-shell replacement.
    assert (events.read_text().splitlines() if events.exists() else []) == expected
    if stage:
        assert f"FAILED: {stage};" in result.stdout
        assert "[NGS2] DONE:" not in result.stdout
    else:
        assert "[NGS2] DONE: Nextflow and QC report completed." in result.stdout
        assert "QC: FAIL=1 ERROR=0" in result.stdout
    copied_log = tmp_path / "output" / batch / "nextflow.log"
    if "nextflow" in expected:
        assert copied_log.read_text().strip() == "run-log"
    else:
        assert not copied_log.exists()  # Preflight failures must not copy an older run's log.


def test_cleanup_secondary_nextflow_work_returns_guarded_dgx_command(monkeypatch):
    work_root = Path("/raid/DGM/work")
    monkeypatch.setattr(secondary, "SECONDARY_DGX_WORK_ROOT", work_root)

    result = secondary.cleanup_nf_work_command()

    assert result["path"] == "/raid/DGM/work"
    assert 'if [ ! -d "${SECONDARY_WORK_ROOT}" ]' in result["command"]
    assert "pgrep -af '[n]extflow'" in result["command"]
    assert 'find "${SECONDARY_WORK_ROOT}" -mindepth 1 -maxdepth 1 -print' in result["command"]
    assert 'read -r -p "確定刪除以上二級分析 Nextflow 暫存？[y/N] "' in result["command"]
    assert 'if ! find "${SECONDARY_WORK_ROOT}"' in result["command"]
    assert '-exec rm -rf -- {} +; then' in result["command"]
    assert '-maxdepth 1 -print -quit)' in result["command"]
    assert "部分檔案無法刪除" in result["command"]
    assert "目錄仍有殘留項目" in result["command"]


def test_cleanup_secondary_nextflow_work_reports_success_only_after_empty(
    monkeypatch, tmp_path
):
    work_root = tmp_path / "work"
    _touch(work_root / "batch" / "task" / ".command.out")
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    _fake_command(fake_bin / "pgrep", "exit 1")
    monkeypatch.setattr(secondary, "SECONDARY_DGX_WORK_ROOT", work_root)
    command = secondary.cleanup_nf_work_command()["command"]
    env = os.environ.copy()
    env["PATH"] = f"{fake_bin}{os.pathsep}{env['PATH']}"

    completed = subprocess.run(
        ["bash", "-c", command],
        input="y\n",
        text=True,
        capture_output=True,
        env=env,
        check=False,
    )

    assert completed.returncode == 0
    assert list(work_root.iterdir()) == []
    assert f"已清理：{work_root}" in completed.stdout


def test_cleanup_secondary_nextflow_work_does_not_report_false_success(
    monkeypatch, tmp_path
):
    work_root = tmp_path / "work"
    leftover = _touch(work_root / "batch" / "task" / ".command.out")
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    _fake_command(fake_bin / "pgrep", "exit 1")
    _fake_command(fake_bin / "rm", "exit 1")
    monkeypatch.setattr(secondary, "SECONDARY_DGX_WORK_ROOT", work_root)
    command = secondary.cleanup_nf_work_command()["command"]
    env = os.environ.copy()
    env["PATH"] = f"{fake_bin}{os.pathsep}{env['PATH']}"

    completed = subprocess.run(
        ["bash", "-c", command],
        input="y\n",
        text=True,
        capture_output=True,
        env=env,
        check=False,
    )

    assert completed.returncode == 1
    assert leftover.exists()
    assert "部分檔案無法刪除" in completed.stdout
    assert "已清理：" not in completed.stdout
