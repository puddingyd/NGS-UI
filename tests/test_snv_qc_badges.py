import csv

from app.adapters.snv_tsv import _row_to_variant, load_snv_tsv, merge_snv_variant_row


def _row(**updates):
    row = {
        "CHROM": "chr1",
        "POS": "100",
        "REF": "A",
        "ALT": "G",
        "GENE": "TEST1",
        "TRANSCRIPT": "ENST1",
        "TRANSCRIPT_TYPE": "MANE_SELECT",
        "HGVS_C": "c.1A>G",
        "HGVS_P": "p.Lys1Arg",
        "CONSEQUENCE": "missense_variant",
        "ACMG_CRITERIA": "",
        "ACMG_SCORE": "0",
        "ACMG_CLASS": "Uncertain significance",
        "DP_DV": "31",
        "AD_DV": "0,8",
        "VAF_DV": "0.258",
        "STRAND_BIAS": "WARN(FS=61.2,SOR=3.4)",
    }
    row.update(updates)
    return row


def test_strand_bias_and_alt_support_are_normalized():
    variant = _row_to_variant(_row())

    assert variant["strand_bias_status"] == "warn"
    assert variant["strand_bias_fs"] == 61.2
    assert variant["strand_bias_sor"] == 3.4
    assert variant["strand_bias_threshold"] == "SNV: FS>60 or SOR>3.0"
    assert variant["depth"] == 31
    assert variant["alt_depth"] == 8
    assert variant["low_alt_support"] is True


def test_haploid_het_is_exposed_as_manual_review_flag():
    variant = _row_to_variant(_row(HAPLOID_HET="DV,HC"))

    assert variant["haploid_het"] is True
    assert variant["haploid_het_callers"] == "DV,HC"
    assert _row_to_variant(_row(HAPLOID_HET="."))["haploid_het"] is False


def test_missing_alt_ad_stays_missing_instead_of_becoming_zero():
    variant = _row_to_variant(_row(DP_DV="10", AD_DV="10,.", VAF_DV="."))

    assert variant["AD"] == "10,."
    assert variant["alt_depth"] is None
    assert variant["alt_af"] is None
    assert variant["low_alt_support"] is False


def test_pipeline_zygosity_is_used_without_chr_x_reinterpretation():
    variant = _row_to_variant(_row(CHROM="chrX", ZYGOSITY="hom", GT_DV="1/1"))

    assert variant["zygosity"] == "hom"


def test_manual_and_legacy_strand_bias_are_distinct():
    assert _row_to_variant(_row(STRAND_BIAS="."))["strand_bias_status"] == "manual"
    legacy = _row()
    legacy.pop("STRAND_BIAS")
    assert _row_to_variant(legacy)["strand_bias_status"] == ""


def test_zero_depth_wgs_is_low_but_wes_is_filtered(tmp_path):
    row = _row(DP_DV="0", AD_DV="0,0", VAF_DV="0")
    path = tmp_path / "sample.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    wgs, _ = load_snv_tsv(path, test_type="WGS")
    assert wgs["chr1-100-A-G"]["depth"] == 0
    assert wgs["chr1-100-A-G"]["low_depth"] is True
    wes, _ = load_snv_tsv(path, test_type="WES")
    assert wes == {}


def test_callers_none_is_not_reportable(tmp_path):
    row = _row(CALLERS="NONE")
    path = tmp_path / "sample.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    variants, _ = load_snv_tsv(path, test_type="WGS")
    assert variants == {}


def test_dragen_reference_genotype_row_is_not_reportable(tmp_path):
    row = _row(CALLERS="DRAGEN", ZYGOSITY="ref", GT_DV="0/0")
    path = tmp_path / "sample.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    variants, _ = load_snv_tsv(path, test_type="WGS")
    assert variants == {}


def test_phase_unknown_overlapping_rows_remain_separate(tmp_path):
    deletion = _row(
        POS="100", REF="AT", ALT="A", HGVS_C="c.100_101del",
        HGVS_P="p.Test34del", ZYGOSITY="het", GT_DV="0/1",
    )
    overlapping_snv = _row(
        POS="101", REF="T", ALT="C", HGVS_C="c.101T>C",
        HGVS_P="p.Test34Ala", ZYGOSITY="het", GT_DV="0/1",
    )
    path = tmp_path / "sample.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(deletion), delimiter="\t")
        writer.writeheader()
        writer.writerow(deletion)
        writer.writerow(overlapping_snv)

    variants, _ = load_snv_tsv(path, test_type="WGS")
    assert set(variants) == {"chr1-100-AT-A", "chr1-101-T-C"}


def test_each_transcript_option_keeps_its_own_disease_associated_scope():
    variants = {}
    gnpda2 = _row_to_variant(_row(
        GENE="GNPDA2",
        TRANSCRIPT="ENST00000609092",
        TRANSCRIPT_TYPE="BEST_CONSEQUENCE",
        HGVS_C="c.353G>C",
        HGVS_P="p.Ter118SerextTer3",
        CONSEQUENCE="stop_lost",
    ))
    guf1 = _row_to_variant(_row(
        GENE="GUF1",
        TRANSCRIPT="ENST00000281543",
        TRANSCRIPT_TYPE="MANE_SELECT",
        HGVS_C="c.514C>G",
        HGVS_P="p.Gln172Glu",
        CONSEQUENCE="missense_variant",
    ))

    merge_snv_variant_row(variants, gnpda2)
    merged = merge_snv_variant_row(variants, guf1)
    options = {option["gene_symbol"]: option for option in merged["transcript_options"]}

    assert options["GNPDA2"]["disease_associated"] is False
    assert options["GUF1"]["disease_associated"] is True
    assert merged["gene_symbol"] == "GNPDA2"
    assert merged["disease_associated"] is False
