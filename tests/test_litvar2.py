import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from app.adapters.snv_tsv import _row_to_variant
from app.services import litvar2_store


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "annotate_litvar2.py"
REPO_ROOT = Path(__file__).resolve().parents[1]
APP_JS = (REPO_ROOT / "frontend" / "app.js").read_text(encoding="utf-8")
INDEX_HTML = (REPO_ROOT / "frontend" / "index.html").read_text(encoding="utf-8")
SPEC = importlib.util.spec_from_file_location("annotate_litvar2", SCRIPT)
annotate_litvar2 = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(annotate_litvar2)


def _records():
    return [
        {
            "_id": "litvar-123",
            "rsid": "rs123",
            "gene": ["GENE1"],
            "all_hgvs": ["NM_000001.1:c.1A>G", "NP_000001.1:p.Lys1Arg"],
            "pmids": ["9005", "9004", "9003", "9002", "9001", "9000"],
            "pmids_count": 8,
            "data_tax_id": [9606],
        },
        {
            "_id": "litvar-hgvs",
            "gene": ["GENE2"],
            "all_hgvs": ["NM_000002.1:c.2C>T"],
            "pmids": ["8100"],
            "pmids_count": 1,
            "data_tax_id": "9606",
        },
        {
            "_id": "litvar-ambiguous-a",
            "gene": ["GENE3"],
            "all_hgvs": ["NM_000003.1:c.3G>A"],
            "pmids": [],
            "data_tax_id": 9606,
        },
        {
            "_id": "litvar-ambiguous-b",
            "gene": ["GENE3"],
            "all_hgvs": ["NM_000003.2:c.3G>A"],
            "pmids": [],
            "data_tax_id": 9606,
        },
        {
            "_id": "mouse-only",
            "rsid": "rs999",
            "gene": ["MouseGene"],
            "pmids": ["1"],
            "data_tax_id": 10090,
        },
    ]


def _write_bulk(path: Path, records=None, *, json_lines=False) -> Path:
    records = _records() if records is None else records
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        if json_lines:
            for record in records:
                handle.write(json.dumps(record) + "\n")
        else:
            json.dump(records, handle)
    return path


def _build_db(tmp_path: Path) -> tuple[Path, Path]:
    bulk = _write_bulk(tmp_path / "source.json.gz")
    db = tmp_path / "litvar2.sqlite"
    litvar2_store.build_database(
        bulk,
        db,
        dataset_date="2026-07-01",
        source_url="https://example.test/litvar.json.gz",
    )
    return bulk, db


def test_bulk_build_and_lookup_preserves_source_pmid_order(tmp_path):
    _bulk, db = _build_db(tmp_path)

    with litvar2_store.open_readonly(db) as conn:
        hit = litvar2_store.lookup_variant(conn, rsids=["RS123"])
        fallback = litvar2_store.lookup_variant(
            conn,
            genes=["gene2"],
            hgvs_values=["NM_000002.1:c.2C>T"],
        )
        ambiguous = litvar2_store.lookup_variant(
            conn,
            genes=["GENE3"],
            hgvs_values=["c.3G>A"],
        )
        mouse = litvar2_store.lookup_variant(conn, rsids=["rs999"])

    assert hit["status"] == "hit"
    assert hit["pmids"] == ["9005", "9004", "9003", "9002", "9001"]
    assert hit["pmids_count"] == 8
    assert hit["dataset_date"] == "2026-07-01"
    assert hit["url"].startswith(
        "https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?"
    )
    assert fallback["status"] == "hit"
    assert fallback["match_method"] == "gene_hgvs"
    assert ambiguous["status"] == "ambiguous"
    assert mouse["status"] == "no_match"


def test_stream_parser_accepts_json_lines(tmp_path):
    bulk = _write_bulk(tmp_path / "lines.json.gz", json_lines=True)
    assert [row["_id"] for row in litvar2_store.iter_bulk_records(bulk)] == [
        row["_id"] for row in _records()
    ]


def test_failed_rebuild_preserves_live_bulk_and_database(tmp_path):
    source = _write_bulk(tmp_path / "good-source.json.gz")
    live_bulk = tmp_path / "litvar2_variants.json.gz"
    live_db = tmp_path / "litvar2.sqlite"
    manifest = tmp_path / "litvar2_manifest.json"
    litvar2_store.update_database(
        directory=tmp_path,
        bulk_path=live_bulk,
        db_path=live_db,
        manifest_path=manifest,
        local_bulk=source,
        dataset_date="2026-07-01",
    )
    old_bulk_hash = hashlib.sha256(live_bulk.read_bytes()).hexdigest()

    broken = tmp_path / "broken.json.gz"
    with gzip.open(broken, "wt", encoding="utf-8") as handle:
        handle.write("[{")
    with pytest.raises(litvar2_store.LitVar2Error):
        litvar2_store.update_database(
            directory=tmp_path,
            bulk_path=live_bulk,
            db_path=live_db,
            manifest_path=manifest,
            local_bulk=broken,
            dataset_date="2026-08-01",
            force=True,
        )

    assert hashlib.sha256(live_bulk.read_bytes()).hexdigest() == old_bulk_hash
    assert litvar2_store.database_metadata(live_db)["dataset_date"] == "2026-07-01"


def test_postprocessing_only_annotates_review_candidates(tmp_path, monkeypatch):
    _bulk, db = _build_db(tmp_path)
    tsv = tmp_path / "working.tsv"
    fields = [
        "CHROM", "POS", "REF", "ALT", "GENE", "TRANSCRIPT", "HGVS_C",
        "HGVS_P", "CONSEQUENCE", "RS_ID", "DP_DV", "GNOMAD_G_AF",
        "CLINVAR_SIG", "ACMG_CRITERIA",
    ]
    rows = [
        {
            "CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "G",
            "GENE": "GENE1", "TRANSCRIPT": "ENST1", "HGVS_C": "c.1A>G",
            "HGVS_P": "p.Lys1Arg", "CONSEQUENCE": "missense_variant",
            "RS_ID": "rs123", "DP_DV": "30", "GNOMAD_G_AF": "0.0001",
            "CLINVAR_SIG": "", "ACMG_CRITERIA": "",
        },
        {
            "CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "G",
            "GENE": "GENE1", "TRANSCRIPT": "ENST2", "HGVS_C": "c.1A>G",
            "HGVS_P": "p.Lys1Arg", "CONSEQUENCE": "missense_variant",
            "RS_ID": "rs123", "DP_DV": "30", "GNOMAD_G_AF": "0.0001",
            "CLINVAR_SIG": "", "ACMG_CRITERIA": "",
        },
        {
            "CHROM": "chr2", "POS": "200", "REF": "C", "ALT": "T",
            "GENE": "GENE2", "TRANSCRIPT": "ENST3", "HGVS_C": "c.2C>T",
            "HGVS_P": "", "CONSEQUENCE": "missense_variant",
            "RS_ID": "", "DP_DV": "10", "GNOMAD_G_AF": "0.0001",
            "CLINVAR_SIG": "", "ACMG_CRITERIA": "",
        },
    ]
    with tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    monkeypatch.setattr(annotate_litvar2.snv_review, "load_candidate_bed", lambda: None)

    stats = annotate_litvar2.annotate_tsv(tsv, db, test_type="WES")
    with tsv.open(encoding="utf-8", newline="") as handle:
        annotated = list(csv.DictReader(handle, delimiter="\t"))

    assert stats["candidate_variants"] == 1
    assert stats["hit"] == 1
    assert annotated[0]["LITVAR2_STATUS"] == "hit"
    assert annotated[1]["LITVAR2_ID"] == annotated[0]["LITVAR2_ID"]
    assert annotated[2]["LITVAR2_STATUS"] == ""
    payload = _row_to_variant(annotated[0])["litvar2"]
    assert payload["pmids"] == ["9005", "9004", "9003", "9002", "9001"]
    assert payload["url"].startswith("https://www.ncbi.nlm.nih.gov/research/litvar2/")


def test_adapter_rejects_untrusted_litvar_url():
    payload = _row_to_variant({
        "CHROM": "chr1",
        "POS": "1",
        "REF": "A",
        "ALT": "G",
        "LITVAR2_URL": "https://example.com/not-litvar",
        "LITVAR2_PMIDS_TOP5": "123,456",
        "LITVAR2_PMID_COUNT": "2",
    })["litvar2"]
    assert payload["url"] == ""
    assert payload["pmids"] == ["123", "456"]


def test_frontend_places_litvar_below_acmg_with_required_links():
    card = APP_JS[
        APP_JS.index("function renderVariantCard"):
        APP_JS.index("// ---- helpers used by renderVariantCard")
    ]
    assert card.index('class="k">ACMG') < card.index("${renderLitvar2(v)}")
    renderer = APP_JS[
        APP_JS.index("function renderLitvar2"):
        APP_JS.index("// ---------- Render: sample header")
    ]
    assert "https://pubmed.ncbi.nlm.nih.gov/" in renderer
    assert "and ${remaining} others" in renderer
    assert "litvar2-title-link" not in renderer
    assert "litvar2-external-link" in renderer
    assert "LITVAR2_EXTERNAL_ICON_SVG" in renderer
    assert 'title="在 LitVar2 開啟"' in renderer
    assert 'let value = "NA (請重跑三級)"' in renderer
    assert 'data.status === "no_match"' in renderer
    assert 'value = "No reference"' in renderer


def test_manual_update_button_is_between_pgx_and_index_refresh():
    pgx = INDEX_HTML.index('id="dragen-pgx"')
    litvar = INDEX_HTML.index('id="dragen-litvar2-update-btn"')
    refresh = INDEX_HTML.index('id="dragen-refresh-btn"')
    assert pgx < litvar < refresh


def test_update_status_is_in_tertiary_modal_header():
    head = INDEX_HTML[
        INDEX_HTML.index('<div class="dragen-modal-head">'):
        INDEX_HTML.index('<div class="dragen-form">')
    ]
    assert 'id="dragen-litvar2-status"' in head
