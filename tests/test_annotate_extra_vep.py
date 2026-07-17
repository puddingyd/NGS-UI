import csv
import gzip
import importlib.util
from pathlib import Path

from app.adapters.snv_tsv import _row_to_variant


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "annotate_extra_vep.py"
SPEC = importlib.util.spec_from_file_location("annotate_extra_vep", SCRIPT)
extra_vep = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(extra_vep)


def _write_dbnsfp(path: Path, *, with_revel: bool) -> None:
    fields = ["#chr", "pos(1-based)", "MetaRNN_score", "MetaRNN_pred"]
    if with_revel:
        fields.append("REVEL_score")
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\t".join(fields) + "\n")


def test_select_dbnsfp_fields_enables_revel_when_header_has_it(tmp_path):
    dbnsfp = tmp_path / "dbNSFP.tsv.gz"
    _write_dbnsfp(dbnsfp, with_revel=True)

    selected, missing = extra_vep.select_dbnsfp_fields(dbnsfp)

    assert selected == ["MetaRNN_score", "MetaRNN_pred", "REVEL_score"]
    assert missing == []


def test_select_dbnsfp_fields_safely_skips_revel_for_legacy_branch(tmp_path):
    dbnsfp = tmp_path / "dbNSFP4.9c.tsv.gz"
    _write_dbnsfp(dbnsfp, with_revel=False)

    selected, missing = extra_vep.select_dbnsfp_fields(dbnsfp)

    assert selected == ["MetaRNN_score", "MetaRNN_pred"]
    assert missing == ["REVEL_score"]


def test_parse_and_merge_revel_from_picked_vep_transcript(tmp_path):
    vcf = tmp_path / "vep.vcf"
    csq_fields = [
        "Allele", "Consequence", "PICK", "MetaRNN_score", "MetaRNN_pred",
        "REVEL_score", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
        "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL",
    ]
    non_picked = "G|missense_variant|0|0.20|T|0.10|0.01|0.02|0.03|0.04"
    picked = "G|missense_variant|1|0.70&0.80|D|0.65&0.73|0.01|0.02|0.30|0.04"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: {"|".join(csq_fields)}">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr1\t100\t.\tA\tG\t.\tPASS\tCSQ={non_picked},{picked}\n",
        encoding="utf-8",
    )

    annotations = extra_vep.parse_vep_vcf(vcf)

    assert annotations[("chr1", "100", "A", "G")] == {
        "METARNN": "0.8",
        "METARNN_PRED": "D",
        "REVEL": "0.73",
        "SPLICEAI_MAX": "0.3000",
    }

    source = tmp_path / "source.tsv"
    output = tmp_path / "output.tsv"
    source.write_text(
        "CHROM\tPOS\tREF\tALT\nchr1\t100\tA\tG\n",
        encoding="utf-8",
    )
    assert extra_vep.merge_into_tsv(source, output, annotations) == (1, 1)
    with output.open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["REVEL"] == "0.73"
    assert row["METARNN"] == "0.8"
    assert row["SPLICEAI_MAX"] == "0.3000"


def test_snv_payload_exposes_worst_revel_score():
    variant = _row_to_variant({
        "CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "G",
        "REVEL": "0.65&0.73",
    })

    assert variant["REVEL_score"] == 0.73


def test_default_dbnsfp_is_the_revel_enabled_531a_database():
    assert Path(extra_vep.DEFAULT_DBNSFP).name == "dbNSFP5.3.1a_grch38.gz"
