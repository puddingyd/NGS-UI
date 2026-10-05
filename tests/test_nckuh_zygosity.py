import gzip
import json

import pytest

from app.services import sample_layout, snv_zygosity


def _ploidy(path, karyotype, source=""):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(f"##estimatedSexKaryotype={karyotype}\n")
        if source:
            handle.write(f"##source={source}\n")
        handle.write("#CHROM\tPOS\n")


def _setup_sample(tmp_path, monkeypatch, *, karyotype="XY", mode="inhouse"):
    monkeypatch.setattr(
        sample_layout, "state_file",
        lambda _sample_id, name, **_kwargs: tmp_path / name,
    )
    _ploidy(tmp_path / "ploidy.vcf.gz", karyotype)
    (tmp_path / "pipeline_source.json").write_text(
        json.dumps({"pipeline_type": mode}), encoding="utf-8"
    )


def test_existing_nckuh_variants_are_corrected_on_read(tmp_path, monkeypatch):
    _setup_sample(tmp_path, monkeypatch)
    variants = {
        "x": {"CHROM": "chrX", "POS": 60_000_000, "REF": "A", "zygosity": "hom"},
        "x_alias": {"CHROM": "X", "POS": 60_000_001, "REF": "A", "zygosity": "Homozygous"},
        "par1": {"CHROM": "chrX", "POS": 2_781_479, "REF": "A", "zygosity": "hom"},
        "nonpar": {"CHROM": "chrX", "POS": 2_781_480, "REF": "A", "zygosity": "hom"},
        "par2": {"CHROM": "chrX", "POS": 155_701_383, "REF": "A", "zygosity": "hom"},
        "crosses_par2": {"CHROM": "chrX", "POS": 155_701_382, "REF": "AT", "zygosity": "hom"},
        "flag": {"CHROM": "chrX", "POS": 60_000_002, "REF": "A", "zygosity": "hom", "haploid_het": True},
        "somatic": {"CHROM": "chrX", "POS": 60_000_003, "REF": "A", "zygosity": "hom", "somatic": True},
        "het": {"CHROM": "chrX", "POS": 60_000_004, "REF": "A", "zygosity": "het"},
        "autosome": {"CHROM": "chr1", "POS": 60_000_000, "REF": "A", "zygosity": "hom"},
    }

    assert snv_zygosity.normalize_loaded_variants(variants, "S1-nckuh") == 3
    assert {key for key, variant in variants.items() if variant["zygosity"] == "hemi"} == {
        "x", "x_alias", "nonpar"
    }
    assert snv_zygosity.normalize_loaded_variants(variants, "S1-nckuh") == 0


@pytest.mark.parametrize("karyotype", ["XX", "XXY", "X", ""])
def test_non_xy_ploidy_keeps_hom(tmp_path, monkeypatch, karyotype):
    _setup_sample(tmp_path, monkeypatch, karyotype=karyotype)
    variants = {"x": {"CHROM": "chrX", "POS": 60_000_000, "REF": "A", "zygosity": "hom"}}

    assert snv_zygosity.normalize_loaded_variants(variants, "S1-nckuh") == 0
    assert variants["x"]["zygosity"] == "hom"


def test_source_must_be_nckuh_and_saved_ploidy_must_exist(tmp_path, monkeypatch):
    _setup_sample(tmp_path, monkeypatch, mode="dragen")
    variants = {"x": {"CHROM": "chrX", "POS": 60_000_000, "REF": "A", "zygosity": "hom"}}
    assert snv_zygosity.normalize_loaded_variants(variants, "S1-nckuh") == 0

    (tmp_path / "pipeline_source.json").write_text(
        json.dumps({"pipeline_type": "inhouse"}), encoding="utf-8"
    )
    (tmp_path / "ploidy.vcf.gz").unlink()
    assert snv_zygosity.normalize_loaded_variants(variants, "S1-nckuh") == 0


def test_nckuh_ploidy_source_supports_legacy_sample_without_sidecar(tmp_path, monkeypatch):
    _setup_sample(tmp_path, monkeypatch)
    (tmp_path / "pipeline_source.json").unlink()
    _ploidy(tmp_path / "ploidy.vcf.gz", "XY", "NCKUH_PLOIDY_MOSDEPTH")
    variants = {"x": {"CHROM": "chrX", "POS": 60_000_000, "REF": "A", "zygosity": "hom"}}

    assert snv_zygosity.normalize_loaded_variants(variants, "S1") == 1
    assert variants["x"]["zygosity"] == "hemi"
