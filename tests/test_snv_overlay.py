import csv
import sqlite3

from app.services import snv_gene_index, snv_overlay, snv_review


FIELDS = [
    "CHROM", "POS", "REF", "ALT", "GENE", "TRANSCRIPT", "HGVS_C",
    "HGVS_P", "CONSEQUENCE", "CLINVAR_SIG", "GNOMAD_G_AF", "DP",
]


def _write(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_sparse_overlay_merges_annotations_and_allows_filtered_legacy_rows(tmp_path, monkeypatch):
    raw = tmp_path / "raw.tsv"
    annotated = tmp_path / "annotated.tsv"
    overlay_path = tmp_path / "snv_annotations.sqlite"
    rows = [
        dict(zip(FIELDS, ["chr1", "10", "A", "*", "SKIP", "T0", "", "", "intron_variant", "", "", "30"])),
        dict(zip(FIELDS, ["chr1", "20", "G", "A", "MUTYH", "T1", "c.1G>A", "", "intron_variant", "Pathogenic", "0.35", "30"])),
    ]
    _write(raw, FIELDS, rows)
    annotated_fields = FIELDS + ["GENEBE_CLASSIFICATION"]
    enriched = dict(rows[1])
    enriched["GENEBE_CLASSIFICATION"] = "Pathogenic"
    _write(annotated, annotated_fields, [enriched])

    snv_overlay.build_overlay(raw, annotated, overlay_path)
    assert snv_overlay.is_current(raw, overlay_path)
    with snv_overlay.OverlayReader(raw, overlay_path) as overlay:
        merged = overlay.apply(rows[1])
    assert merged["GENEBE_CLASSIFICATION"] == "Pathogenic"

    monkeypatch.setattr(snv_review, "_candidate_bed_path", lambda: tmp_path / "missing.bed")
    review = snv_review.ensure_review_tsv(
        raw,
        test_type="WGS",
        output_dir=tmp_path / "08_postprocessing",
        overlay_path=overlay_path,
    )
    with review.open("r", encoding="utf-8", newline="") as handle:
        review_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(review_rows) == 1
    assert review_rows[0]["GENE"] == "MUTYH"
    assert review_rows[0]["GENEBE_CLASSIFICATION"] == "Pathogenic"


def test_gene_index_uses_explicit_postprocessing_path_and_skips_star(tmp_path):
    raw = tmp_path / "03_acmg" / "S1.snv_indel.acmg.tsv"
    raw.parent.mkdir()
    rows = [
        dict(zip(FIELDS, ["chr1", "10", "A", "*", "MUTYH", "T0", "", "", "intron_variant", "", "", "30"])),
        dict(zip(FIELDS, ["chr1", "20", "G", "A", "MUTYH", "T1", "", "", "intron_variant", "", "", "30"])),
    ]
    _write(raw, FIELDS, rows)
    index = tmp_path / "08_postprocessing" / "snv_gene_index.sqlite"
    snv_gene_index.build_index(raw, index)
    found = snv_gene_index.query_rows(raw, ["MUTYH"], index)
    assert found is not None
    assert [row["ALT"] for row in found] == ["A"]
    with sqlite3.connect(index) as conn:
        assert conn.execute("SELECT count(*) FROM variants").fetchone()[0] == 1
