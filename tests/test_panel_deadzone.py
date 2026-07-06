from app.services import panel_deadzone
from app.services import snv_gene_index


def test_variant_gene_canonicalization_keeps_ercc6_identity():
    assert panel_deadzone.canonical_gene_symbol("ERCC6", "HGNC:3438") == (
        "ERCC6",
        "HGNC:3438",
    )
    assert panel_deadzone.canonical_gene_symbol("ERCC6") == (
        "ERCC6",
        "HGNC:3438",
    )


def test_variant_gene_canonicalization_uses_hgnc_aliases_not_positional_aliases():
    assert panel_deadzone.canonical_gene_symbol("RAD26") == (
        "ERCC6",
        "HGNC:3438",
    )
    assert panel_deadzone.canonical_gene_symbol("NDUFA4") == (
        "COXFA4",
        "HGNC:7687",
    )
    assert panel_deadzone.canonical_gene_symbol("PGBD3", "HGNC:19400") == (
        "PGBD3",
        "HGNC:19400",
    )


def test_snv_gene_index_uses_variant_identity_for_ercc6(tmp_path):
    raw_tsv = tmp_path / "snv_indel.annotated.tsv"
    raw_tsv.write_text(
        "CHROM\tPOS\tREF\tALT\tGENE\tHGNC_ID\n"
        "chr10\t49458903\tA\tG\tERCC6\tHGNC:3438\n",
        encoding="utf-8",
    )

    snv_gene_index.build_index(raw_tsv)

    rows = snv_gene_index.query_rows(raw_tsv, ["ERCC6"])
    assert rows and rows[0]["GENE"] == "ERCC6"
    assert snv_gene_index.query_rows(raw_tsv, ["PGBD3"]) == []
