import gzip

import pytest

from app.services.snv_overlay import OverlayReader, build_overlay
from app.services.snv_zygosity import copy_nckuh_work_tsv


HEADER = "CHROM\tPOS\tREF\tALT\tGENE\tTRANSCRIPT\tHGVS_C\tHGVS_P\tCONSEQUENCE\tZYGOSITY\tGT_DV\tGT_HC\tHAPLOID_HET\n"


def _row(chrom, pos, zygosity="hom", ref="A", flag="."):
    return f"{chrom}\t{pos}\t{ref}\tG\tTEST\tENST1\tc.1A>G\tp.Lys1Arg\tmissense_variant\t{zygosity}\t1/1\t1/1\t{flag}\n"


def _ploidy(path, karyotype):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(f"##estimatedSexKaryotype={karyotype}\n#CHROM\tPOS\n")


def test_xy_nckuh_chr_x_nonpar_hom_is_hemi_in_overlay(tmp_path):
    raw = tmp_path / "raw.tsv"
    work = tmp_path / "work.tsv"
    overlay_path = tmp_path / "overlay.sqlite"
    ploidy = tmp_path / "ploidy.vcf.gz"
    _ploidy(ploidy, "XY")
    raw.write_text(
        HEADER
        + _row("chrX", 2_781_479)  # PAR1 endpoint
        + _row("chrX", 2_781_480)
        + _row("X", 60_000_000)
        + _row("chrX", 155_701_383)  # PAR2 start
        + _row("chrX", 155_701_382, ref="AT")  # REF overlaps PAR2
        + _row("chr1", 60_000_000)
        + _row("chrX", 60_000_001, zygosity="het")
        + _row("chrX", 60_000_002, flag="DV"),
        encoding="utf-8",
    )

    assert copy_nckuh_work_tsv(raw, work, ploidy) == 2
    assert raw.read_text(encoding="utf-8").count("\themi\t") == 0
    rows = [line.split("\t") for line in work.read_text(encoding="utf-8").splitlines()[1:]]
    assert [row[9] for row in rows] == ["hom", "hemi", "hemi", "hom", "hom", "hom", "het", "hom"]
    assert rows[1][10:12] == ["1/1", "1/1"]
    build_overlay(raw, work, overlay_path)
    raw_rows = [
        dict(zip(HEADER.strip().split("\t"), line.split("\t")))
        for line in raw.read_text(encoding="utf-8").splitlines()[1:]
    ]
    with OverlayReader(raw, overlay_path) as overlay:
        corrected = overlay.apply_many(raw_rows)
    assert corrected[1]["ZYGOSITY"] == "hemi"


@pytest.mark.parametrize("karyotype", ["XX", "XXY", "X", ""])
def test_non_xy_ploidy_does_not_reclassify(tmp_path, karyotype):
    raw = tmp_path / "raw.tsv"
    work = tmp_path / "work.tsv"
    ploidy = tmp_path / "ploidy.vcf.gz"
    raw.write_text(HEADER + _row("chrX", 60_000_000), encoding="utf-8")
    _ploidy(ploidy, karyotype)

    assert copy_nckuh_work_tsv(raw, work, ploidy) == 0
    assert work.read_bytes() == raw.read_bytes()
    assert copy_nckuh_work_tsv(raw, work, None) == 0
    assert work.read_bytes() == raw.read_bytes()
