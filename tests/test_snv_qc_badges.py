import csv

from app.adapters.snv_tsv import _row_to_variant, load_snv_tsv


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
