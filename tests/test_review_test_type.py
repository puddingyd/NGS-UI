"""build_snv_review_tsv.py records the test type the backend loader asks for."""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_snv_review_tsv.py"
_spec = importlib.util.spec_from_file_location("build_snv_review_tsv", SCRIPT)
brt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(brt)

RAW = (
    "CHROM\tPOS\tREF\tALT\tGENE\tTRANSCRIPT\tHGVS_C\tHGVS_P\tCONSEQUENCE\tIMPACT\tDP_DV\n"
    "chr1\t1000\tA\tG\tAAA\tNM_1\tc.1\tp.1\tmis\tMODERATE\t30\n"
    "chr1\t1100\tC\tT\tAAA\tNM_1\tc.2\tp.2\tmis\tMODERATE\t12\n"
)


@pytest.mark.parametrize("requested, sample, expected", [
    ("WGS", "26T00028-dragen", "TITAN-WGS"),   # worker says WGS, loader wants TITAN-WGS
    ("WES", "26T00028-dragen", "TITAN-WGS"),
    ("WGS", "VAL-24-WGS", "WGS"),
    ("WES", "25WE0001-nckuh", "WES"),
    ("", "25WE0001-nckuh", "WES"),
])
def test_resolve_test_type(tmp_path, requested, sample, expected):
    post = tmp_path / "stage" / "08_postprocessing"    # staging: parent != sample ID
    post.mkdir(parents=True)
    raw = tmp_path / "raw.tsv"
    raw.write_text(RAW, encoding="utf-8")
    assert brt.resolve_test_type(raw, post, requested, sample) == expected


def test_sample_inferred_from_postprocessing_dir_and_metadata(tmp_path):
    post = tmp_path / "26T00029-dragen" / "08_postprocessing"
    post.mkdir(parents=True)
    raw = tmp_path / "raw.tsv"
    raw.write_text(RAW, encoding="utf-8")
    assert brt.resolve_test_type(raw, post) == "TITAN-WGS"
    other = tmp_path / "S1" / "08_postprocessing"
    other.mkdir(parents=True)
    (other / "S1.sample_metadata.json").write_text(
        json.dumps({"test_type": "wgs"}), encoding="utf-8")
    assert brt.resolve_test_type(raw, other) == "WGS"


def test_manifest_records_titan_wgs_and_keeps_low_dp_rows(tmp_path, monkeypatch):
    bed = tmp_path / "cds.bed"
    bed.write_text("chr1\t900\t2000\n", encoding="utf-8")
    monkeypatch.setenv("NGS_UI_CDS_CANDIDATE_BED", str(bed))
    post = tmp_path / "stage" / "08_postprocessing"
    post.mkdir(parents=True)
    raw = tmp_path / "raw.tsv"
    raw.write_text(RAW, encoding="utf-8")
    review = post / "26T00028-dragen.snv_indel.review.tsv"
    manifest = post / "26T00028-dragen.snv_indel.review.tsv.source.json"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--tsv", str(raw), "--output-dir", str(post),
         "--output-path", str(review), "--manifest-path", str(manifest),
         "--test-type", "WGS", "--sample", "26T00028-dragen",
         "--gpn-msa-db", str(tmp_path / "missing.bgz")],
        check=True, capture_output=True, env={**__import__("os").environ,
                                               "NGS_UI_CDS_CANDIDATE_BED": str(bed)},
    )
    assert json.loads(manifest.read_text(encoding="utf-8"))["test_type"] == "TITAN-WGS"
    rows = review.read_text(encoding="utf-8").splitlines()
    assert len(rows) == 3          # DP 12 row kept: no WES DP>=20 floor


def test_run_stopgaps_passes_sample_to_review_builder():
    script = (SCRIPT.parent / "run_stopgaps.sh").read_text(encoding="utf-8")
    call = script[script.index('"$SCRIPT_DIR/build_snv_review_tsv.py"'):]
    call = call[:call.index("step_done")]
    assert '--sample "$SID"' in call
