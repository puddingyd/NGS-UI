from app import config
from app.services import sample_layout
import json


def _roots(tmp_path, monkeypatch):
    new = tmp_path / "new"
    old_ui = tmp_path / "old_ui"
    old_pipeline = tmp_path / "old_pipeline"
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", new)
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", old_ui)
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", old_pipeline)
    return new, old_ui, old_pipeline


def test_legacy_state_wins_until_layout_marker_exists(tmp_path, monkeypatch):
    new, old_ui, _old_pipeline = _roots(tmp_path, monkeypatch)
    legacy = old_ui / "S1"
    legacy.mkdir(parents=True)
    (legacy / "snv_indel.annotated.tsv").write_text(
        "CHROM\tPOS\tREF\tALT\n", encoding="utf-8"
    )

    assert sample_layout.state_dir("S1") == legacy
    assert sample_layout.snv_raw_tsv("S1") == legacy / "snv_indel.annotated.tsv"

    raw = new / "S1" / "03_acmg" / "S1.snv_indel.acmg.tsv"
    raw.parent.mkdir(parents=True)
    raw.write_text("CHROM\tPOS\tREF\tALT\n", encoding="utf-8")
    sample_layout.write_layout_marker("S1", source_id="S1", raw_tsv=raw)

    assert sample_layout.state_dir("S1") == new / "S1" / "08_postprocessing"
    assert sample_layout.snv_raw_tsv("S1") == raw


def test_pipeline_aux_files_are_read_directly(tmp_path, monkeypatch):
    new, _old_ui, _old_pipeline = _roots(tmp_path, monkeypatch)
    sample = new / "S2"
    raw = sample / "03_acmg" / "S2.snv_indel.acmg.tsv"
    cnv = sample / "06_cnv_sv" / "S2.cnv.annotated.tsv"
    pgx = sample / "07_pgx" / "S2.pgx.tsv"
    for path in (raw, cnv, pgx):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x\n", encoding="utf-8")
    sample_layout.write_layout_marker("S2", source_id="S2", raw_tsv=raw)

    assert sample_layout.cnv_tsv("S2") == cnv
    assert sample_layout.pgx_tsv("S2") == pgx


def test_v2_unprefixed_state_remains_readable_and_writable(tmp_path, monkeypatch):
    new, _old_ui, _old_pipeline = _roots(tmp_path, monkeypatch)
    post = new / "S3" / "08_postprocessing"
    post.mkdir(parents=True)
    (post / "layout.json").write_text(
        json.dumps({"layout_version": 2, "sample_id": "S3"}),
        encoding="utf-8",
    )
    legacy_meta = post / "sample_metadata.json"
    legacy_meta.write_text("{}", encoding="utf-8")

    assert sample_layout.uses_unified_layout("S3")
    assert sample_layout.state_file("S3", "sample_metadata.json") == legacy_meta
    assert sample_layout.state_file(
        "S3", "sample_metadata.json", for_write=True
    ) == legacy_meta


def test_v3_prefixed_state_wins_when_both_exist(tmp_path, monkeypatch):
    new, _old_ui, _old_pipeline = _roots(tmp_path, monkeypatch)
    post = new / "S4" / "08_postprocessing"
    post.mkdir(parents=True)
    (post / "layout.json").write_text(
        json.dumps({"layout_version": 2, "sample_id": "S4"}),
        encoding="utf-8",
    )
    old_meta = post / "sample_metadata.json"
    old_meta.write_text('{"old":true}', encoding="utf-8")
    new_meta = post / "S4.sample_metadata.json"
    new_meta.write_text('{"new":true}', encoding="utf-8")
    (post / "S4.layout.json").write_text(
        json.dumps({"layout_version": 3, "sample_id": "S4"}),
        encoding="utf-8",
    )

    assert sample_layout.layout_marker_path("S4") == post / "S4.layout.json"
    assert sample_layout.state_file("S4", "sample_metadata.json") == new_meta


def test_full_reprocess_copies_v2_state_without_deleting_old_files(tmp_path, monkeypatch):
    new, _old_ui, _old_pipeline = _roots(tmp_path, monkeypatch)
    post = new / "S5" / "08_postprocessing"
    analysis = post / "analyses" / "default"
    analysis.mkdir(parents=True)
    old_meta = post / "sample_metadata.json"
    old_analysis = analysis / "analysis.json"
    old_marker = post / "layout.json"
    old_meta.write_text('{"sample_id":"S5"}', encoding="utf-8")
    old_analysis.write_text('{"hpo":[]}', encoding="utf-8")
    old_marker.write_text('{"layout_version":2}', encoding="utf-8")

    copied = sample_layout.promote_state_tree_to_v3("S5")

    assert len(copied) == 2
    assert (post / "S5.sample_metadata.json").is_file()
    assert (analysis / "S5.analysis.json").is_file()
    assert not (post / "S5.layout.json").exists()
    assert old_meta.is_file()
    assert old_analysis.is_file()
    assert old_marker.is_file()


def test_staged_layout_marker_keeps_raw_path_relative(tmp_path):
    sample = tmp_path / "stage" / "SRC"
    raw = sample / "03_acmg" / "SRC.snv_indel.acmg.tsv"
    raw.parent.mkdir(parents=True)
    raw.touch()

    marker = sample_layout.write_layout_marker_in_sample_dir(
        sample,
        "S1-dragen",
        source_id="SRC",
        raw_tsv=raw,
    )
    payload = json.loads(marker.read_text(encoding="utf-8"))

    assert marker == (
        sample / "08_postprocessing" / "S1-dragen.layout.json"
    )
    assert payload["sample_id"] == "S1-dragen"
    assert payload["source_sample_id"] == "SRC"
    assert payload["raw_tsv"] == "03_acmg/SRC.snv_indel.acmg.tsv"
