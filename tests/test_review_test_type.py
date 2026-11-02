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
         "--gpn-star-dir", str(tmp_path / "missing.bgz")],
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


@pytest.mark.parametrize('built_type,loaded_type', [
    ('WGS','TITAN-WGS'), ('TITAN-WGS','WGS'),
])
def test_equivalent_wgs_labels_reuse_review_without_raw_scan(tmp_path,monkeypatch,built_type,loaded_type):
    from app.services import snv_review
    bed=tmp_path/'cds.bed';bed.write_text('chr1\t900\t2000\n')
    monkeypatch.setenv('NGS_UI_CDS_CANDIDATE_BED',str(bed))
    raw=tmp_path/'raw.tsv';raw.write_text(RAW)
    options=dict(output_dir=tmp_path/'post',gpn_msa_db=tmp_path/'no-msa',gpn_star_dir=tmp_path/'no-star')
    review=snv_review.ensure_review_tsv(raw,test_type=built_type,**options)
    before=review.read_bytes(),review.stat().st_mtime_ns
    def no_scan(*args,**kwargs):
        raise AssertionError('Equivalent WGS labels must reuse the existing review TSV')
    monkeypatch.setattr(snv_review,'_load_bed',no_scan)
    monkeypatch.setattr(snv_review.gpn_msa,'annotate_review_tsv',no_scan)
    assert snv_review.ensure_review_tsv(raw,test_type=loaded_type,**options)==review
    assert (review.read_bytes(),review.stat().st_mtime_ns)==before
    assert len(review.read_text().splitlines())==3  # WGS keeps DP12.


def test_wes_depth_policy_still_rebuilds_wgs_review(tmp_path,monkeypatch,capsys):
    from app.services import snv_review
    bed=tmp_path/'cds.bed';bed.write_text('chr1\t900\t2000\n')
    monkeypatch.setenv('NGS_UI_CDS_CANDIDATE_BED',str(bed))
    raw=tmp_path/'raw.tsv';raw.write_text(RAW)
    options=dict(output_dir=tmp_path/'post',gpn_msa_db=tmp_path/'no-msa',gpn_star_dir=tmp_path/'no-star')
    review=snv_review.ensure_review_tsv(raw,test_type='WGS',**options)
    capsys.readouterr()
    snv_review.ensure_review_tsv(raw,test_type='WES',**options)
    assert len(review.read_text().splitlines())==2
    log=capsys.readouterr().out
    assert 'snv_review.rebuild_start' in log and 'test_type' in log


@pytest.mark.parametrize('changed', ['raw','bed'])
def test_equivalent_labels_do_not_hide_real_input_changes(tmp_path,monkeypatch,changed):
    from app.services import snv_review
    bed=tmp_path/'cds.bed';bed.write_text('chr1\t900\t2000\n')
    monkeypatch.setenv('NGS_UI_CDS_CANDIDATE_BED',str(bed))
    raw=tmp_path/'raw.tsv';raw.write_text(RAW)
    options=dict(output_dir=tmp_path/'post',gpn_msa_db=tmp_path/'no-msa',gpn_star_dir=tmp_path/'no-star')
    review=snv_review.ensure_review_tsv(raw,test_type='WGS',**options)
    if changed=='raw':
        raw.write_text(RAW+'chr1\t1200\tA\tT\tBBB\tNM_2\tc.3\tp.3\tmis\tMODERATE\t30\n')
        expected_lines=4
    else:
        bed.write_text('chr1\t900\t1050\n')
        expected_lines=2
    snv_review.ensure_review_tsv(raw,test_type='TITAN-WGS',**options)
    assert len(review.read_text().splitlines())==expected_lines
