import csv

from app.services import manual_acmg, sample_layout, sample_loader, snv_gene_index


def _row(*, gene: str, transcript: str, transcript_type: str, hgvs_c: str, hgvs_p: str):
    return {
        "CHROM": "chr16",
        "POS": "2097465",
        "REF": "AGAA",
        "ALT": "A",
        "GENE": gene,
        "TRANSCRIPT": transcript,
        "REFSEQ_TRANSCRIPT": transcript,
        "TRANSCRIPT_TYPE": transcript_type,
        "HGVS_C": hgvs_c,
        "HGVS_P": hgvs_p,
        "CONSEQUENCE": "frameshift_variant",
        "ACMG_SCORE": "10",
        "ACMG_CLASS": "Pathogenic",
        "ZYGOSITY": "het",
        "DP_DV": "80",
        "AD_DV": "45,35",
        "VAF_DV": "0.4375",
    }


def test_case_summary_keeps_all_transcripts_and_honors_saved_selection(
    tmp_path, monkeypatch
):
    sample_dir = tmp_path / "S1"
    sample_dir.mkdir()
    raw_tsv = sample_dir / "snv.tsv"
    raw_tsv.write_text("placeholder\n", encoding="utf-8")
    rows = [
        _row(
            gene="TSC2",
            transcript="NM_000548.5",
            transcript_type="MANE_SELECT",
            hgvs_c="c.4258_4261del",
            hgvs_p="p.Ser1420GlyfsTer55",
        ),
        _row(
            gene="PKD1",
            transcript="NM_001009944.3",
            transcript_type="BEST_CONSEQUENCE",
            hgvs_c="c.1_4del",
            hgvs_p="p.Test1fs",
        ),
    ]

    monkeypatch.setattr(sample_layout, "snv_raw_tsv", lambda sample_id: raw_tsv)
    monkeypatch.setattr(
        sample_layout, "snv_gene_index_path", lambda sample_id: sample_dir / "index.sqlite"
    )
    monkeypatch.setattr(
        sample_layout, "snv_overlay_path", lambda sample_id: sample_dir / "overlay.sqlite"
    )
    monkeypatch.setattr(
        sample_layout,
        "state_file",
        lambda sample_id, name: sample_dir / name,
    )
    monkeypatch.setattr(snv_gene_index, "query_rows_by_ids", lambda *args: rows)
    monkeypatch.setattr(sample_loader, "_enrich_snv_variants", lambda *args: {})
    monkeypatch.setattr(manual_acmg, "bulk_current", lambda *args: {})
    monkeypatch.setattr(
        sample_loader.omim_store,
        "lookup_cached",
        lambda **kwargs: {
            "Disease1": (
                "Tuberous sclerosis-2 (613254)(AD)"
                if kwargs.get("gene") == "TSC2"
                else "Polycystic kidney disease 1 (173900)(AD)"
            )
        },
    )

    variant_id = "chr16-2097465-AGAA-A"
    variants = sample_loader._case_snv_variants_by_id(sample_dir, {variant_id})
    variant = variants[variant_id]
    options = {option["gene_symbol"]: option for option in variant["transcript_options"]}

    assert set(options) == {"TSC2", "PKD1"}
    edits = {"selected_transcript_key": options["TSC2"]["key"]}
    assert sample_loader._case_variant_label(variant, edits).startswith(
        "TSC2:NM_000548.5:c.4258_4261del:p.Ser1420GlyfsTer55, P, het"
    )
    edits["report_diseases"] = {"1": True}
    assert sample_loader._case_selected_diseases(variant, edits) == [
        "Tuberous sclerosis-2 (613254)(AD)"
    ]


def test_case_summary_signature_includes_schema_version(tmp_path, monkeypatch):
    sample_dir = tmp_path / "S1"
    sample_dir.mkdir()
    monkeypatch.setattr(sample_loader.omim_store, "cache_signature", lambda: ())
    monkeypatch.setattr(sample_loader.gene_disease_store, "cache_signature", lambda: ())

    signature = sample_loader._case_summary_signature(sample_dir)

    assert signature[0] == ["case_summary_version", sample_loader.CASE_SUMMARY_VERSION]


def test_case_summary_review_fallback_reads_every_transcript(tmp_path):
    review_tsv = tmp_path / "review.tsv"
    rows = [
        _row(
            gene="TSC2",
            transcript="NM_000548.5",
            transcript_type="MANE_SELECT",
            hgvs_c="c.4258_4261del",
            hgvs_p="p.Ser1420GlyfsTer55",
        ),
        _row(
            gene="PKD1",
            transcript="NM_001009944.3",
            transcript_type="BEST_CONSEQUENCE",
            hgvs_c="c.1_4del",
            hgvs_p="p.Test1fs",
        ),
    ]
    with review_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    matched = sample_loader._scan_snv_review_rows_by_ids(
        review_tsv, {"chr16-2097465-AGAA-A"}
    )

    assert [row["GENE"] for row in matched] == ["TSC2", "PKD1"]
