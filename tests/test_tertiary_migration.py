import csv
import importlib.util
from pathlib import Path

from app import config
from app.services import sample_layout


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "migrate_tertiary_output_layout.py"
SPEC = importlib.util.spec_from_file_location("migrate_tertiary_output_layout", SCRIPT)
migration = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(migration)


def _write_tsv(path, *, enriched=False):
    fields = [
        "CHROM", "POS", "REF", "ALT", "GENE", "TRANSCRIPT", "HGVS_C",
        "HGVS_P", "CONSEQUENCE", "CLINVAR_SIG", "GNOMAD_G_AF", "DP",
    ]
    if enriched:
        fields.append("GENEBE_CLASSIFICATION")
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "CHROM": "chr1", "POS": "20", "REF": "G", "ALT": "A",
        "GENE": "MUTYH", "TRANSCRIPT": "T1", "HGVS_C": "c.1G>A",
        "HGVS_P": "", "CONSEQUENCE": "intron_variant",
        "CLINVAR_SIG": "Pathogenic", "GNOMAD_G_AF": "0.35", "DP": "30",
        "GENEBE_CLASSIFICATION": "Pathogenic",
    }
    if not enriched:
        row.pop("GENEBE_CLASSIFICATION")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def test_canary_migration_activates_marker_last_and_keeps_legacy(tmp_path, monkeypatch):
    old_pipeline = tmp_path / "old_pipeline"
    old_ui = tmp_path / "old_ui"
    target = tmp_path / "target"
    raw = old_pipeline / "S1" / "03_acmg" / "S1.snv_indel.acmg.tsv"
    old_full = old_ui / "S1" / "snv_indel.annotated.tsv"
    _write_tsv(raw)
    _write_tsv(old_full, enriched=True)
    pipeline_str = old_pipeline / "S1" / "05_str" / "S1.str.tsv"
    pipeline_str.parent.mkdir(parents=True)
    pipeline_str.write_text("same\n", encoding="utf-8")
    (old_ui / "S1" / "str.tsv").write_text("same\n", encoding="utf-8")
    pipeline_cnv = old_pipeline / "S1" / "06_cnv_sv" / "S1.cnv.annotated.tsv"
    pipeline_cnv.parent.mkdir(parents=True)
    pipeline_cnv.write_text("pipeline\n", encoding="utf-8")
    (old_ui / "S1" / "cnv.annotated.tsv").write_text("locally enriched\n", encoding="utf-8")
    (old_ui / "S1" / "sample_metadata.json").write_text(
        '{"sample_id":"S1","lis_id":"S1","test_type":"WGS"}',
        encoding="utf-8",
    )
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", target)
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", old_ui)
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", old_pipeline)

    assert migration.migrate_one(
        "S1",
        source_pipeline_root=old_pipeline,
        source_ui_root=old_ui,
        target_root=target,
        apply=True,
    )

    post = target / "S1" / "08_postprocessing"
    assert (post / "S1.layout.json").is_file()
    assert (post / "S1.snv_annotations.sqlite").is_file()
    assert (post / "S1.snv_indel.review.tsv").is_file()
    assert (post / "S1.snv_gene_index.sqlite").is_file()
    assert not (post / "snv_indel.annotated.tsv").exists()
    assert not (post / "str.tsv").exists()
    assert (post / "S1.cnv.annotated.tsv").read_text(encoding="utf-8") == "locally enriched\n"
    assert old_full.is_file()
    assert sample_layout.state_dir("S1") == post
