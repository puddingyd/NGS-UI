from app import config
from app.services import sample_layout


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
