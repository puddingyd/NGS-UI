from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

from scripts import secondary_qc_report as qc


def write_bed(tmp_path, text="chr1\t10\t20\nchr1\t15\t25\nchr1\t40\t45\nchr2\t0\t5\n"):
    path = tmp_path / "targets.bed"
    path.write_text(text)
    return path


def metrics(path):
    path.write_text(
        "## METRICS CLASS\tpicard.sam.DuplicationMetrics\n"
        "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tPERCENT_DUPLICATION\n"
        "lib1\t10\t5\t2\t1\t0.2\n"
        "lib2\t0\t10\t0\t4\t0.4\n\n"
        "## HISTOGRAM\nBIN\tVALUE\n1\t100\n"
    )


def test_bed_union_boundaries_and_zero_coverage_contig(tmp_path):
    targets = qc.Targets(write_bed(tmp_path))
    assert targets.length == 25
    assert targets.regions["chr1"] == [(10, 25), (40, 45)]
    assert not targets.overlap("chr1", 5, 10)
    assert not targets.overlap("chr1", 25, 40)
    assert targets.overlap("chr1", 24, 41)
    hist = qc.depth_histogram(["chr1\t11\t2\n", "chr1\t12\t3\n"], targets)
    assert hist == Counter({0: 23, 2: 1, 3: 1})
    with pytest.raises(ValueError, match="duplicate"):
        qc.depth_histogram(["chr1\t11\t2\n"] * 2, targets)


def test_on_target_cigar_requires_real_aligned_bases_and_counts_once(tmp_path):
    targets = qc.Targets(write_bed(tmp_path))
    def sam(pos, cigar):
        return f"read\t0\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t*\t*\n"
    assert qc.count_target_reads([
        sam(6, "5M15D5M"), sam(6, "5M15N5M"),  # Target lies entirely in D/N.
        sam(16, "2S30M3S"),  # Spans two target intervals, still one read.
        sam(41, "2=1X2="),
    ], targets) == 2


def test_duplication_weights_multiple_libraries(tmp_path):
    path = tmp_path / "dup.txt"
    metrics(path)
    assert qc.duplication_fraction(path) == (12, 40)


def test_exact_thresholds_and_rounding_do_not_change_pass_fail():
    # Mean 50, exactly 90% >=10X, target exactly 40%, mapping exactly 95%.
    hist = Counter({10: 8, 420: 1, 0: 1})
    result = qc.summarize("S1", 30000000, 28500000, 11400000, (12, 40), hist, 10)
    assert result["row"]["QC"] == "PASS"
    assert result["row"]["Uniformity"] == "90.00%"
    assert result["counts"]["uniformity_depth_cutoff"] == 10
    failed = qc.summarize("S1", 30000000, 28499999, 11400000, (12, 40), hist, 10)
    assert failed["row"]["Mapping rate"] == "95.00%"
    assert failed["row"]["QC"] == "FAIL"
    assert failed["failed_checks"] == ["Mapping rate >= 95%"]
    # Unrounded mean is 50.001, hence bases at 10X fail the uniformity cutoff.
    near = qc.summarize("S1", 30000000, 30000000, 12000000,
                        (0, 1), Counter({10: 999, 40011: 1}), 1000)
    assert near["row"]["Mean depth"] == "50.00"
    assert near["counts"]["uniformity_depth_cutoff"] == 11
    assert near["row"]["Uniformity"] == "0.10%"


def test_all_zero_coverage_is_not_reported_as_uniform():
    result = qc.summarize("S1", 0, 0, 0, (0, 0), Counter({0: 10}), 10)
    assert result["row"]["Uniformity"] == "NA"
    assert result["row"]["QC"] == "FAIL"
    assert len(result["failed_checks"]) == 5


@pytest.mark.parametrize("target,hist,territory,field,display,failed_check", [
    (11999999, Counter({50: 1000}), 1000, "On target rate", "40.00%", "On target rate >= 40%"),
    (12000000, Counter({50: 999, 49: 1}), 1000, "Mean depth", "50.00", "Mean depth >= 50X"),
    (12000000, Counter({100: 899999, 0: 100001}), 1000000, "Uniformity", "90.00%", "Uniformity >= 90%"),
])
def test_values_rounded_to_threshold_still_fail(target, hist, territory, field, display, failed_check):
    result = qc.summarize("S1", 30000000, 30000000, target, (0, 1), hist, territory)
    assert result["row"][field] == display
    assert result["failed_checks"] == [failed_check]


def test_stream_nonzero_exit_rejects_partial_stdout(tmp_path):
    program = tmp_path / "partial.py"
    program.write_text("import sys\nprint('partial data')\nprint('read error', file=sys.stderr)\nsys.exit(7)\n")
    tool = qc.Samtools.__new__(qc.Samtools)
    tool.prefix = [sys.executable, str(program)]
    with pytest.raises(RuntimeError, match=r"failed \(7\): read error"):
        with tool.stream(["depth"]) as lines:
            assert list(lines) == ["partial data\n"]


def test_samplesheet_deduplicates_lanes_and_rejects_traversal(tmp_path):
    path = tmp_path / "samplesheet.csv"
    path.write_text("sample,lane\nS2,L001\nS1,L001\nS2,L002\n")
    assert qc.samples_from_sheet(path) == ["S2", "S1"]
    path.write_text("sample\n../S1\n")
    with pytest.raises(ValueError, match="Invalid"):
        qc.samples_from_sheet(path)


def test_resolved_config_literal_paths_and_override(tmp_path):
    path = tmp_path / "config.flat"
    path.write_text(
        "params.wes_targets = '/ref/kit targets.bed'\n"
        "process.withName:SAMTOOLS.*.container = '/sif/samtools_1.23.1.sif'\n"
        "untrusted = __import__('os').system('false')\n"
    )
    assert qc.resolved_settings(path) == ("/ref/kit targets.bed", "/sif/samtools_1.23.1.sif")
    path.write_text(path.read_text() + "process.withName:SAMTOOLS_OTHER.container = '/other.sif'\n")
    with pytest.raises(ValueError, match="Multiple"):
        qc.resolved_settings(path)
    assert qc.resolved_settings(path, runtime_override=True) == ("/ref/kit targets.bed", None)


@pytest.fixture
def bam_batch(tmp_path):
    samtools = shutil.which("samtools")
    if not samtools:
        pytest.skip("Real Samtools is required for BAM integration tests")
    out = tmp_path / "batch"
    alignment = out / "S1" / "02_alignment"
    alignment.mkdir(parents=True)
    sam = tmp_path / "fixture.sam"
    lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:100", "@SQ\tSN:chr2\tLN:100",
             "@RG\tID:rg1\tSM:S1"]
    def read(name, flag=0, pos=11, cigar="10M", length=10, mq=60, quality="I", mate="*", mate_pos=0, tlen=0):
        chrom = "*" if flag & 4 else "chr1"
        lines.append(f"{name}\t{flag}\t{chrom}\t{pos}\t{mq}\t{cigar}\t{mate}\t{mate_pos}\t{tlen}\t"
                     + "A" * length + "\t" + quality * length + "\tRG:Z:rg1")
    read("normal")
    read("pair", 99, mate="=", mate_pos=16, tlen=15)
    read("pair", 147, pos=16, mate="=", mate_pos=11, tlen=-15)
    read("duplicate", 1024)
    read("secondary", 256)
    read("supplementary", 2048)
    read("qcfail", 512)
    read("lowmq", mq=19)
    read("lowbq", quality="4")  # Phred 19.
    read("off", pos=71)
    read("deletion_only", pos=6, cigar="5M15D5M")
    read("skip_only", pos=6, cigar="5M15N5M")
    read("two_intervals", pos=16, cigar="30M", length=30)
    read("equals_mismatch", pos=41, cigar="2=1X2=", length=5)
    read("unmapped", flag=4, pos=0, cigar="*", mq=0)
    sam.write_text("\n".join(lines) + "\n")
    bam = alignment / "S1.aligned.sorted.bam"
    subprocess.run([samtools, "sort", "-o", str(bam), str(sam)], check=True, capture_output=True)
    subprocess.run([samtools, "index", str(bam)], check=True, capture_output=True)
    metrics(alignment / "S1.duplicate_metrics.txt")
    bed = write_bed(tmp_path)
    sheet = out / "samplesheet.csv"
    sheet.write_text("sample,lane\nS1,L001\nS1,L002\nMISSING,L001\n")
    args = ["--out-dir", str(out), "--samplesheet", str(sheet), "--target-bed", str(bed), "--samtools", samtools]
    return out, bed, args


def test_real_bam_counts_quality_filters_overlap_and_partial_failure(bam_batch):
    out, _, args = bam_batch
    assert qc.main(args + ["--check-only"]) == 0
    assert not (out / "pipeline_info/report_summary.csv").exists()
    assert qc.main(args) == 2
    with (out / "pipeline_info/report_summary.csv").open() as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == qc.FIELDS
        rows = list(reader)
    assert rows == [
        {"Sample ID": "S1", "Total reads": "13", "Duplicated rate": "30.00%", "Mapping rate": "92.31%",
         "On target rate": "75.00%", "Mean depth": "1.80", "Uniformity": "80.00%", "QC": "FAIL"},
        dict(zip(qc.FIELDS, ["MISSING", "NA", "NA", "NA", "NA", "NA", "NA", "ERROR"])),
    ]
    cached = json.loads((out / "S1/03_alignment_qc/S1.report_qc.json").read_text())
    assert cached["counts"]["depth_sum"] == 45  # Normal 10 + pair union 15 + multi 15 + =/X 5.
    assert cached["counts"]["target_bases"] == 25
    assert cached["counts"]["zero_depth_bases"] == 5  # Entire chr2 target is uncovered.
    assert json.loads((out / "pipeline_info/report_summary.details.json").read_text())["state"] == "ERROR"


def test_cache_reuse_invalidation_force_and_fail_exit_success(bam_batch, capsys):
    out, bed, args = bam_batch
    (out / "samplesheet.csv").write_text("sample\nS1\n")
    assert qc.main(args) == 0  # Measured FAIL is successful report generation.
    cache = out / "S1/03_alignment_qc/S1.report_qc.json"
    original = cache.read_bytes()
    capsys.readouterr()
    assert qc.main(args) == 0
    assert "(cached)" in capsys.readouterr().out
    assert cache.read_bytes() == original
    bed.write_text(bed.read_text() + "chr2\t5\t10\n")
    assert qc.main(args) == 0
    assert "(cached)" not in capsys.readouterr().out
    assert json.loads(cache.read_text())["counts"]["target_bases"] == 30
    assert qc.main(args + ["--force"]) == 0
    assert "(cached)" not in capsys.readouterr().out


def test_target_reference_mismatch_is_error(bam_batch):
    out, bed, args = bam_batch
    bed.write_text("1\t0\t20\n")
    assert qc.main(args) == 2
    detail = json.loads((out / "pipeline_info/report_summary.details.json").read_text())
    assert "BED contig/coordinates" in detail["samples"][0]["error"]


def test_container_runtime_and_resolved_config_run_same_calculation(bam_batch, monkeypatch, tmp_path):
    out, bed, args = bam_batch
    native_samtools = args[-1]
    engine_dir = tmp_path / "runtime"
    engine_dir.mkdir()
    engine = engine_dir / "apptainer"
    calls = tmp_path / "container_calls.jsonl"
    # Exercise the real container command construction while using local Samtools.
    engine.write_text(f"#!{sys.executable}\n" +
                      "import json, os, sys\n" +
                      f"with open({str(calls)!r}, 'a') as f: f.write(json.dumps(sys.argv[1:]) + '\\n')\n" +
                      "i = sys.argv.index('samtools')\n" +
                      f"os.execv({native_samtools!r}, ['samtools'] + sys.argv[i + 1:])\n")
    engine.chmod(0o755)
    monkeypatch.setenv("PATH", str(engine_dir))
    sif = tmp_path / "samtools image.sif"
    sif.touch()
    config = tmp_path / "config.flat"
    config.write_text(f"params.wes_targets = {str(bed)!r}\nprocess.withName:SAMTOOLS.*.container = {str(sif)!r}\n")
    (out / "samplesheet.csv").write_text("sample\nS1\n")
    assert qc.main(args[:4] + ["--nextflow-config", str(config)]) == 0
    detail = json.loads((out / "pipeline_info/report_summary.details.json").read_text())
    assert detail["samples"][0]["counts"]["depth_sum"] == 45
    invoked = [json.loads(line) for line in calls.read_text().splitlines()]
    assert all(call[:2] == ["exec", "--cleanenv"] for call in invoked)
    assert all(f"{out}:{out}" in call and str(sif) in call for call in invoked)
    assert all("--nv" not in call for call in invoked)
