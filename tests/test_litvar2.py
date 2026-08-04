import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from app.adapters.snv_tsv import _row_to_variant
from app.services import litvar2_on_demand, litvar2_store, snv_gene_index


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
    assert [candidate["litvar_id"] for candidate in ambiguous["candidates"]] == [
        "litvar-ambiguous-a",
        "litvar-ambiguous-b",
    ]
    assert all(candidate["gene"] == "GENE3" for candidate in ambiguous["candidates"])
    assert all(candidate["hgvs"].upper().endswith("C.3G>A") for candidate in ambiguous["candidates"])
    assert all(
        candidate["url"].startswith("https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?")
        for candidate in ambiguous["candidates"]
    )
    browser_ambiguous = litvar2_on_demand._browser_payload(ambiguous)
    assert [candidate["id"] for candidate in browser_ambiguous["candidates"]] == [
        "litvar-ambiguous-a",
        "litvar-ambiguous-b",
    ]
    assert mouse["status"] == "no_match"


def test_stream_parser_accepts_json_lines(tmp_path):
    bulk = _write_bulk(tmp_path / "lines.json.gz", json_lines=True)
    assert [row["_id"] for row in litvar2_store.iter_bulk_records(bulk)] == [
        row["_id"] for row in _records()
    ]


def test_legacy_ambiguous_cache_is_refreshed_for_candidate_details(tmp_path):
    cache = tmp_path / "litvar2_on_demand.sqlite"
    db_fingerprint = "db-v1"
    litvar2_on_demand._write_results(
        cache,
        [
            (
                "chr1-1-A-G",
                db_fingerprint,
                "identifiers-a",
                {"status": "ambiguous", "dataset_date": "2026-07-01"},
                "manual",
            ),
            (
                "chr1-2-C-T",
                db_fingerprint,
                "identifiers-b",
                {"status": "no_match", "dataset_date": "2026-07-01"},
                "manual",
            ),
        ],
    )

    cached = litvar2_on_demand._cached_rows(
        cache,
        ["chr1-1-A-G", "chr1-2-C-T"],
        db_fingerprint,
    )

    assert "chr1-1-A-G" not in cached
    assert cached["chr1-2-C-T"]["status"] == "no_match"


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

    marker = tmp_path / "sample.litvar2_annotation.json"
    stats = annotate_litvar2.annotate_tsv(
        tsv,
        db,
        test_type="WES",
        marker_path=marker,
    )
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
    marker_payload = json.loads(marker.read_text(encoding="utf-8"))
    assert marker_payload["status"] == "complete"
    assert marker_payload["scope"] == "review_candidates"
    assert marker_payload["dataset_date"] == "2026-07-01"
    assert marker_payload["candidate_variants"] == 1


def test_postprocessing_round_trips_all_ambiguous_candidates(tmp_path, monkeypatch):
    _bulk, db = _build_db(tmp_path)
    tsv = tmp_path / "ambiguous.tsv"
    fields = [
        "CHROM", "POS", "REF", "ALT", "GENE", "TRANSCRIPT", "HGVS_C",
        "HGVS_P", "CONSEQUENCE", "RS_ID", "DP_DV", "GNOMAD_G_AF",
        "CLINVAR_SIG", "ACMG_CRITERIA",
    ]
    with tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({
            "CHROM": "chr3", "POS": "300", "REF": "G", "ALT": "A",
            "GENE": "GENE3", "TRANSCRIPT": "ENST3", "HGVS_C": "c.3G>A",
            "HGVS_P": "", "CONSEQUENCE": "missense_variant", "RS_ID": "",
            "DP_DV": "30", "GNOMAD_G_AF": "0.0001", "CLINVAR_SIG": "",
            "ACMG_CRITERIA": "",
        })
    monkeypatch.setattr(annotate_litvar2.snv_review, "load_candidate_bed", lambda: None)

    stats = annotate_litvar2.annotate_tsv(tsv, db, test_type="WES")
    with tsv.open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))

    assert stats["ambiguous"] == 1
    assert row["LITVAR2_STATUS"] == "ambiguous"
    raw_candidates = json.loads(row["LITVAR2_CANDIDATES_JSON"])
    assert len(raw_candidates) == 2
    payload = _row_to_variant(row)["litvar2"]
    assert [candidate["id"] for candidate in payload["candidates"]] == [
        "litvar-ambiguous-a",
        "litvar-ambiguous-b",
    ]


def test_on_demand_lookup_uses_marker_and_versioned_sample_cache(tmp_path, monkeypatch):
    bulk, db = _build_db(tmp_path)
    raw_tsv = tmp_path / "sample.snv_indel.acmg.tsv"
    index_path = tmp_path / "sample.snv_gene_index.sqlite"
    marker_path = tmp_path / "sample.litvar2_annotation.json"
    cache_path = tmp_path / "sample.litvar2_on_demand.sqlite"
    fields = ["CHROM", "POS", "REF", "ALT", "GENE", "HGNC_ID", "RS_ID", "HGVS_C", "HGVS_P"]
    with raw_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({
            "CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "G",
            "GENE": "GENE1", "HGNC_ID": "", "RS_ID": "rs123",
            "HGVS_C": "NM_000001.1:c.1A>G", "HGVS_P": "NP_000001.1:p.Lys1Arg",
        })
    snv_gene_index.build_index(raw_tsv, index_path)
    marker_path.write_text(json.dumps({
        "schema_version": 1,
        "status": "complete",
        "scope": "review_candidates",
        "dataset_date": "2026-07-01",
    }), encoding="utf-8")

    monkeypatch.setattr(litvar2_store, "LITVAR2_DB", db)
    monkeypatch.setattr(litvar2_on_demand.sample_layout, "snv_raw_tsv", lambda _sid: raw_tsv)
    monkeypatch.setattr(litvar2_on_demand.sample_layout, "review_tsv", lambda _sid: raw_tsv)
    monkeypatch.setattr(
        litvar2_on_demand.sample_layout,
        "snv_gene_index_path",
        lambda _sid: index_path,
    )
    monkeypatch.setattr(
        litvar2_on_demand.sample_layout,
        "litvar2_marker_path",
        lambda _sid, for_write=False: marker_path,
    )
    monkeypatch.setattr(
        litvar2_on_demand.sample_layout,
        "litvar2_on_demand_path",
        lambda _sid, for_write=False: cache_path,
    )

    first = litvar2_on_demand.lookup_variants(
        "sample",
        ["chr1-100-A-G"],
        trigger="gene_search",
    )
    assert first["eligible"] is True
    assert first["queried"] == 1
    assert first["results"]["chr1-100-A-G"]["status"] == "hit"
    assert first["results"]["chr1-100-A-G"]["dataset_date"] == "2026-07-01"
    assert cache_path.is_file()

    def fail_if_queried(*_args, **_kwargs):
        raise AssertionError("current-version cache should avoid a second DB lookup")

    real_lookup = litvar2_store.lookup_variant
    monkeypatch.setattr(litvar2_store, "lookup_variant", fail_if_queried)
    second = litvar2_on_demand.lookup_variants(
        "sample",
        ["chr1-100-A-G"],
        trigger="gene_search",
    )
    assert second["cached"] == 1
    variants = {"chr1-100-A-G": {"litvar2": {}}}
    assert litvar2_on_demand.apply_cached(variants, "sample") == 1
    assert variants["chr1-100-A-G"]["litvar2"]["pmid_count"] == 8

    next_db = tmp_path / "litvar2-next.sqlite"
    litvar2_store.build_database(
        bulk,
        next_db,
        dataset_date="2026-08-01",
        source_url="https://example.test/litvar.json.gz",
    )
    monkeypatch.setattr(litvar2_store, "LITVAR2_DB", next_db)
    monkeypatch.setattr(litvar2_store, "lookup_variant", real_lookup)
    assert litvar2_on_demand.apply_cached(variants, "sample") == 0
    refreshed = litvar2_on_demand.lookup_variants(
        "sample", ["chr1-100-A-G"], trigger="gene_search",
    )
    assert refreshed["cached"] == 0
    assert refreshed["queried"] == 1
    assert refreshed["results"]["chr1-100-A-G"]["dataset_date"] == "2026-08-01"


def test_manual_on_demand_lookup_is_allowed_without_postprocessing_marker(tmp_path, monkeypatch):
    _bulk, db = _build_db(tmp_path)
    raw_tsv = tmp_path / "sample.tsv"
    raw_tsv.write_text(
        "CHROM\tPOS\tREF\tALT\tGENE\tHGNC_ID\tRS_ID\tHGVS_C\tHGVS_P\n"
        "chr1\t100\tA\tG\tGENE1\t\trs123\tNM_000001.1:c.1A>G\tNP_000001.1:p.Lys1Arg\n",
        encoding="utf-8",
    )
    index_path = snv_gene_index.build_index(raw_tsv, tmp_path / "index.sqlite")
    missing_marker = tmp_path / "missing-marker.json"
    cache_path = tmp_path / "cache.sqlite"
    monkeypatch.setattr(litvar2_store, "LITVAR2_DB", db)
    monkeypatch.setattr(litvar2_on_demand.sample_layout, "snv_raw_tsv", lambda _sid: raw_tsv)
    monkeypatch.setattr(litvar2_on_demand.sample_layout, "review_tsv", lambda _sid: raw_tsv)
    monkeypatch.setattr(litvar2_on_demand.sample_layout, "snv_gene_index_path", lambda _sid: index_path)
    monkeypatch.setattr(
        litvar2_on_demand.sample_layout,
        "litvar2_marker_path",
        lambda _sid, for_write=False: missing_marker,
    )
    monkeypatch.setattr(
        litvar2_on_demand.sample_layout,
        "litvar2_on_demand_path",
        lambda _sid, for_write=False: cache_path,
    )

    automatic = litvar2_on_demand.lookup_variants(
        "sample", ["chr1-100-A-G"], trigger="gene_search",
    )
    assert automatic["eligible"] is False
    assert not cache_path.exists()
    manual = litvar2_on_demand.lookup_variants(
        "sample", ["chr1-100-A-G"], trigger="manual", force=True,
    )
    assert manual["results"]["chr1-100-A-G"]["status"] == "hit"
    assert cache_path.is_file()


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


def test_adapter_preserves_safe_ambiguous_candidates():
    candidates = [
        {
            "litvar_id": "candidate-a",
            "rsid": "rs300",
            "gene": "GENE3",
            "hgvs": "c.3G>A",
            "pmids_count": 7,
            "pmids": ["101", "102", "103", "104", "105"],
            "url": "https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?variant=candidate-a",
        },
        {
            "litvar_id": "candidate-b",
            "rsid": "",
            "gene": "GENE3",
            "hgvs": "c.3G>A",
            "pmids_count": 1,
            "pmids": ["201"],
            "url": "https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?variant=candidate-b",
        },
    ]
    payload = _row_to_variant({
        "CHROM": "chr3",
        "POS": "3",
        "REF": "G",
        "ALT": "A",
        "LITVAR2_STATUS": "ambiguous",
        "LITVAR2_DATASET_DATE": "2026-08-01",
        "LITVAR2_CANDIDATES_JSON": json.dumps(candidates),
    })["litvar2"]
    assert payload["status"] == "ambiguous"
    assert [candidate["id"] for candidate in payload["candidates"]] == [
        "candidate-a",
        "candidate-b",
    ]
    assert payload["candidates"][0]["pmids"] == ["101", "102", "103", "104", "105"]
    assert payload["candidates"][0]["pmid_count"] == 7


def test_frontend_places_litvar_below_acmg_with_required_links():
    card = APP_JS[
        APP_JS.index("function renderVariantCard"):
        APP_JS.index("// ---- helpers used by renderVariantCard")
    ]
    assert card.index('class="k">ACMG') < card.index("${renderLitvar2(v, id)}")
    renderer = APP_JS[
        APP_JS.index("function _litvar2TitleHtml"):
        APP_JS.index("// ---------- Render: sample header")
    ]
    assert "https://pubmed.ncbi.nlm.nih.gov/" in renderer
    assert "and ${remaining} others" in renderer
    assert "litvar2-title-link" not in renderer
    assert "litvar2-external-link" in renderer
    assert "EXTERNAL_LINK_ICON_SVG" in renderer
    assert "LITVAR2_REFRESH_ICON_SVG" in renderer
    assert "litvar2-refresh-btn" in renderer
    assert "litvar2-ambiguous-details" in renderer
    assert "Ambiguous match (${candidates.length} records)" in renderer
    assert "Ambiguous match（請按重新整理取得候選明細）" in renderer
    assert "litvar2-candidate-link" in renderer
    assert "_litvar2PmidsHtml(candidate, url)" in renderer
    assert 'const externalUrl = url || LITVAR2_HOME_URL' in renderer
    assert '"在 LitVar2 開啟此 variant" : "開啟 LitVar2 首頁"' in renderer
    assert 'LITVAR2_HOME_URL = "https://www.ncbi.nlm.nih.gov/research/litvar2/"' in APP_JS
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
