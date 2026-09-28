from app.services import panel_deadzone
from app.services import snv_gene_index


def test_panel_canonicalization_matches_legacy_rules_for_reference_symbols():
    current = set(panel_deadzone.hgnc_id_to_symbol().values())
    by_upper = panel_deadzone._current_symbol_by_upper()
    aliases = panel_deadzone._panel_gene_alias()
    inputs = current | set(aliases) | {"", "unlisted_gene", " c7ORF50 "}
    for value in inputs:
        sym = value.strip()
        if not sym or sym in current:
            expected = sym
        elif sym.upper() in by_upper:
            expected = by_upper[sym.upper()]
        else:
            mapped = aliases.get(sym) or aliases.get(sym.upper())
            expected = by_upper.get(mapped.upper(), mapped) if mapped else sym
        assert panel_deadzone.canonical_panel_gene_symbol(value) == expected


def test_panel_canonicalization_does_not_scan_hgnc_values(monkeypatch):
    class NoValueScan(dict):
        def values(self):
            raise AssertionError("per-gene lookup must not scan the reference")

    monkeypatch.setattr(panel_deadzone, "hgnc_id_to_symbol", lambda: NoValueScan())
    monkeypatch.setattr(panel_deadzone, "symbol_to_hgnc_id", lambda: {"C7orf50": "HGNC:1"})
    monkeypatch.setattr(panel_deadzone, "_current_symbol_by_upper", lambda: {"C7ORF50": "C7orf50"})
    monkeypatch.setattr(panel_deadzone, "_panel_gene_alias", lambda: {"OLD": "C7orf50"})
    for value in ("C7orf50", "C7ORF50", "old"):
        assert panel_deadzone.canonical_panel_gene_symbol(value) == "C7orf50"
    assert panel_deadzone.canonical_panel_gene_symbol("unlisted") == "unlisted"


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
