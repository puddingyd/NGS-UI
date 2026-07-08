import gzip

from app.services.ploidy import load_sample_ploidy, read_karyotype
from app.workers.dragen_run import _copy_dragen_ploidy_vcf, _find_dragen_ploidy_vcf


def _write_ploidy(path, call):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(f"##estimatedSexKaryotype={call}\n#CHROM\tPOS\n")


def test_find_and_copy_dragen_ploidy_by_source_vcf_basename(tmp_path):
    source_vcf = tmp_path / "VAL-3.hard-filtered.vcf.gz"
    source_vcf.touch()
    expected = tmp_path / "VAL-3.ploidy.vcf.gz"
    _write_ploidy(expected, "XY")
    destination_dir = tmp_path / "tertiary_output" / "LIS-001"
    destination_dir.mkdir(parents=True)

    assert _find_dragen_ploidy_vcf(source_vcf, "VAL-3") == expected
    copied = _copy_dragen_ploidy_vcf(source_vcf, "VAL-3", destination_dir)

    assert copied == (expected, destination_dir / "ploidy.vcf.gz")
    assert read_karyotype(destination_dir / "ploidy.vcf.gz") == "XY"


def test_find_dragen_ploidy_never_uses_another_sample(tmp_path):
    source_vcf = tmp_path / "VAL-3.hard-filtered.vcf.gz"
    source_vcf.touch()
    _write_ploidy(tmp_path / "VAL-4.ploidy.vcf.gz", "XX")

    assert _find_dragen_ploidy_vcf(source_vcf, "VAL-3") is None


def test_missing_dragen_ploidy_removes_stale_ui_sidecar(tmp_path):
    source_vcf = tmp_path / "VAL-3.hard-filtered.vcf.gz"
    source_vcf.touch()
    destination_dir = tmp_path / "tertiary_output" / "LIS-001"
    destination_dir.mkdir(parents=True)
    stale = destination_dir / "ploidy.vcf.gz"
    _write_ploidy(stale, "XX")

    assert _copy_dragen_ploidy_vcf(source_vcf, "VAL-3", destination_dir) is None
    assert not stale.exists()


def test_load_sample_ploidy_preserves_x_and_xxy(tmp_path):
    _write_ploidy(tmp_path / "ploidy.vcf.gz", "XXY")
    assert load_sample_ploidy(tmp_path) == {
        "karyotype": "XXY",
        "source": "ploidy.vcf.gz",
    }

    _write_ploidy(tmp_path / "ploidy.vcf.gz", "X")
    assert read_karyotype(tmp_path / "ploidy.vcf.gz") == "X"
