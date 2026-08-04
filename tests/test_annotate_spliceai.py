import csv
import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "annotate_spliceai.py"
SPEC = importlib.util.spec_from_file_location("annotate_spliceai", SCRIPT)
spliceai = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(spliceai)


def test_spliceai_parser_and_merge_write_no_dbnsfp_fields(tmp_path):
    vcf = tmp_path / "vep.vcf"
    fields = [
        "Allele", "PICK", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
        "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL",
    ]
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: {"|".join(fields)}">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\t.\tA\tG\t.\tPASS\tCSQ=G|1|0.01|0.70|0.20|0.30\n",
        encoding="utf-8",
    )
    annotations = spliceai.parse_vep_vcf(vcf)
    assert annotations[("1", "100", "A", "G")] == "0.7000"

    tsv = tmp_path / "sample.tsv"
    tsv.write_text("CHROM\tPOS\tREF\tALT\n1\t100\tA\tG\n", encoding="utf-8")
    assert spliceai.merge_into_tsv(tsv, tsv, annotations) == (1, 1)
    with tsv.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        row = next(reader)
        assert reader.fieldnames == ["CHROM", "POS", "REF", "ALT", "SPLICEAI_MAX"]
    assert row["SPLICEAI_MAX"] == "0.7000"
