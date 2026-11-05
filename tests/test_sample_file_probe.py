"""Synthetic path acceptance checks for the standalone server probe."""
import importlib.util
import json
from pathlib import Path


spec = importlib.util.spec_from_file_location(
    "sample_file_probe", Path(__file__).resolve().parents[1] / "scripts/probe_sample_files.py"
)
probe = importlib.util.module_from_spec(spec)
spec.loader.exec_module(probe)


def touch(root, relative, content="x"):
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


def test_four_lanes_fastq_only_incomplete_unknown_and_orphan(tmp_path):
    for lane in range(1, 5):
        for read in (1, 2):
            touch(tmp_path, f"run/fastq.gz/26T01345_S2_L{lane:03d}_R{read}_001.fastq.gz")
    touch(tmp_path, "run/fastq.gz/ONLY_S3_L001_R1_001.fastq.gz")
    touch(tmp_path, "run/fastq.gz/unrecognized.fq.gz")
    touch(tmp_path, "run/vcf.gz/ORPHAN.hard-filtered.vcf.gz")
    touch(tmp_path, "run/bam/26T01345.bam")
    touch(tmp_path, "run/bam/26T01345.bai")
    touch(tmp_path, "run2/bam/26T01345.bam")
    report = probe.scan({"all": [str(tmp_path)]})
    assert report["total_fastq_samples"] == 2
    full, only = report["samples"]
    assert full["fastq_count"] == 8
    assert full["incomplete_pairs"] == []
    assert [f["kind"] for f in full["files"]].count("BAM") == 2
    assert any(f["kind"] == "BAI" for f in full["files"])
    assert only["files"] == []
    assert only["incomplete_pairs"][0]["reads"] == ["1"]
    assert len(report["unrecognized_fastqs"]) == 1
    assert report["unmatched_files"][0]["source_id"] == "ORPHAN"
    assert not report["scan_errors"]


def test_reanalysis_and_ui_marker_preserve_exact_source(tmp_path):
    touch(tmp_path, "Reanalysis/RENAMED/ORIGINAL.R1.clean.fastq.gz")
    touch(tmp_path, "Reanalysis/RENAMED/ORIGINAL.R2.clean.fastq.gz")
    touch(tmp_path, "tertiary/RENAMED-nckuh/08_postprocessing/RENAMED-nckuh.pipeline_source.json",
          json.dumps({"source_sample_id": "ORIGINAL"}))
    touch(tmp_path, "tertiary/RENAMED-nckuh/08_postprocessing/pipeline_source.json",
          json.dumps({"source_sample_id": "WRONG_OLD_SOURCE"}))
    touch(tmp_path, "tertiary/RENAMED-nckuh/08_postprocessing/RENAMED-nckuh.cnv.review.tsv")
    report = probe.scan({"all": [str(tmp_path)]}, "ORIGINAL")
    sample = report["samples"][0]
    assert sample["sample_id"] == "RENAMED"
    assert sample["source_ids"] == ["ORIGINAL"]
    assert sample["files"][0]["source_id"] == "ORIGINAL"
    assert sample["files"][0]["ui_sample_id"] == "RENAMED-nckuh"


def test_custom_prs_path_external_to_roots_and_overlapping_roots(tmp_path):
    raw = tmp_path / "raw"
    touch(raw, "run/fastq.gz/A_S1_R1_001.fastq.gz")
    touch(raw, "run/fastq.gz/A_S1_R2_001.fastq.gz")
    touch(raw, "batch/A/04_snv_indel/A.ensemble.fixed.vcf.gz")
    touch(tmp_path, "gvcfs/A.g.vcf.gz")
    touch(tmp_path, "gvcfs/A.g.vcf.gz.tbi")
    cfg = touch(tmp_path, "prs.json", json.dumps({"gvcf_templates": {
        "inhouse": [str(tmp_path / "gvcfs/{source}.g.vcf.gz")]
    }}))
    report = probe.scan({"overlap": [str(raw), str(raw / "run")]}, prs_config=cfg)
    sample = report["samples"][0]
    assert sample["fastq_count"] == 2
    assert {f["kind"] for f in sample["files"]} == {"SNV_VCF", "gVCF", "VCF_INDEX"}
    assert not report["scan_errors"]


def test_missing_root_broken_link_and_directory_link_are_visible(tmp_path):
    (tmp_path / "broken.bam").symlink_to(tmp_path / "absent.bam")
    (tmp_path / "dirlink").symlink_to(tmp_path, target_is_directory=True)
    report = probe.scan({"all": [str(tmp_path), str(tmp_path / "missing")]})
    assert len(report["scan_errors"]) == 3
    assert report["samples"] == []
