from __future__ import annotations

import csv
import os
import shlex
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
    assert 'umask "${ORIGINAL_UMASK}"\n    exit "${status}"' in command
    assert "exec bash -i" in command


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
    assert result.returncode == (0 if not failure else (7 if failure in {"config", "nextflow"} else (2 if failure == "missing_script" else 9)))
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
    monkeypatch.setattr(secondary, "SECONDARY_DGX_WORK_ROOT", Path("/raid/DGM/work"))
    result = secondary.cleanup_nf_work_command()
    argv = shlex.split(result["command"])
    assert result["path"] == "/raid/DGM/work"
    assert argv[:2] == ["python3", "-c"]
    assert argv[3:5] == ["cleanup", "/raid/DGM/work"]
    assert "def acquire_lease" in argv[2]


def test_generated_cleanup_executes_embedded_guard_and_preserves_busy_batch(monkeypatch, tmp_path):
    root = tmp_path / "work"
    busy = _touch(root / "busy" / "task" / ".command.out")
    idle = _touch(root / "idle" / "task" / ".command.out")
    monkeypatch.setattr(secondary, "SECONDARY_DGX_WORK_ROOT", root)
    monkeypatch.setattr(secondary, "SECONDARY_DGX_OUTPUT_ROOT", tmp_path / "output")
    argv = shlex.split(secondary.cleanup_nf_work_command()["command"])
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    _fake_command(fake_bin / "ps", "exit 0")
    from backend.app.services import nextflow_work
    with nextflow_work.acquire_lease(root, "busy"):
        completed = subprocess.run(argv, input="y\n", text=True, capture_output=True, timeout=20,
                                   env={**os.environ, "PATH": f"{fake_bin}:{os.environ['PATH']}"})
    assert completed.returncode == 0, completed.stderr
    assert busy.exists()
    assert not idle.exists()
    assert "已清理 1；保留 1；失敗 0" in completed.stdout


def test_generated_secondary_launch_holds_lease_until_runner_finishes(monkeypatch, tmp_path):
    from backend.app.services import nextflow_work
    monkeypatch.setattr(secondary, "SECONDARY_DGX_WORK_ROOT", tmp_path / "work")
    # Parse the actual tmux command, keeping the embedded Python unchanged.
    command = secondary._launch_command("BATCH", "WGS")
    # Source includes newlines, so use shell parsing on everything after the runner heredoc.
    tmux_args = shlex.split(command.split("\nNGS2_EOF\n", 1)[1])
    pane = tmux_args[5]
    assert len(pane) < 4096  # Do not put the embedded module into a tmux argument.
    argv = shlex.split(pane.removesuffix("; exec bash -i"))
    helper = tmp_path / 'guard.py'
    helper.write_text(command.split("<<'NGS2_GUARD_EOF'\n", 1)[1].split("\nNGS2_GUARD_EOF", 1)[0])
    argv[1] = str(helper)
    runner = tmp_path / "runner.sh"
    runner.write_text('echo ready; read answer')
    argv[-1] = str(runner)
    proc = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    try:
        assert proc.stdout.readline().strip() == "ready"
        with pytest.raises(BlockingIOError):
            nextflow_work.open_lock(nextflow_work._lock_path(tmp_path / "work", "BATCH"), blocking=False)
    finally:
        proc.communicate("done\n", timeout=10)
    handle = nextflow_work.open_lock(nextflow_work._lock_path(tmp_path / "work", "BATCH"), blocking=False)
    handle.close()
