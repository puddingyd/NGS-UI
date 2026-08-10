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


def test_spliceai_sites_use_final_review_predicate_including_clinvar_change(tmp_path):
    tsv = tmp_path / "working.tsv"
    fields = [
        "CHROM", "POS", "REF", "ALT", "GNOMAD_G_AF", "DP",
        "CLINVAR_SIG", "CLINVAR_CHANGE",
    ]
    rows = [
        ["chr1", "10", "A", "G", "0.001", "30", "", ""],
        ["chr1", "20", "C", "T", "0.2", "30", "", ""],
        ["chr1", "30", "G", "A", "0.2", "30", "Pathogenic", ""],
        ["chr1", "40", "T", "C", "0.2", "30", "", "UP_TO_PLP"],
        ["chr1", "50", "A", "C", "0.001", "10", "", "UP_TO_PLP"],
        ["chr1", "60", "A", "T", "0.01", "30", "", ""],
    ]
    with tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)
    sites = tmp_path / "sites.vcf"

    count, dropped = spliceai.tsv_to_sites(
        tsv,
        sites,
        test_type="WES",
        candidate_bed=None,
    )

    assert count == 3
    assert dropped == 3
    records = [
        line for line in sites.read_text(encoding="utf-8").splitlines()
        if not line.startswith("#")
    ]
    assert [line.split("\t")[1] for line in records] == ["10", "30", "40"]
