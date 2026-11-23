"""Before/after simulation of replacing the pipeline ACMG with the rule engine (step 6) — synthetic only."""
import csv

from app.services.acmg_engine import replacement_diff as rd
from test_acmg_engine_validation import LOF, OTHER, publish, refs, row

# missense with only a pipeline ACMG (no GeneBe): falls through to the pipeline today
PIPE_ONLY = row(POS="4000", REF="T", ALT="C", CONSEQUENCE="missense_variant", IMPACT="MODERATE",
                HGVS_P="ENSP00000000001.1:p.Leu90Pro", LOFTEE=".", ACMG_SCORE="5",
                ACMG_CLASS="Uncertain_significance")


def test_helpers():
    assert rd.marked_statuses({"status": {"a": "C,0", "b": "1", "c": "0", "d": "2"}}) == {"a": "C", "b": "1", "d": "2"}
    assert rd.tier_with("2", 4) == "1C" and rd.tier_with("2", 3) == "2" and rd.tier_with("1B", -5) == "1B"
    assert rd.tier_with("1C", None) == "1C"
    manual = rd.manual_choice(None, {"reusable_criteria": {"PVS1": {}}, "reusable_classification": "Pathogenic",
                                     "reusable_score": 10, "criteria": {"PVS1": {}}})
    assert manual == {"source": "manual_global", "class": "P", "points": 10.0}
    # case-specific criteria only (nothing reusable) → falls through
    assert rd.manual_choice(None, {"reusable_criteria": {}, "criteria": {"PS2": {}}}) is None
    assert rd.effective(None, None, {"class": None, "points": None}, "pipeline")["source"] == "none"
    assert rd._direction("VUS", "LP") == "up" and rd._direction("LP", "VUS") == "down"
    assert rd._direction(None, "B") == "gained"


def test_run(tmp_path, monkeypatch):
    from app import config

    root = tmp_path / "root"
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", root)
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", tmp_path / "legacy_ui")
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", tmp_path / "legacy_pipe")
    shared = refs()
    sid = "24WE0001-nckuh"
    publish(root, sid, [LOF, OTHER, PIPE_ONLY], shared)
    publish(root, "24WE0002-nckuh", [LOF], shared)

    metas = {sid: {"status": {"chr1-4000-T-C": "1", "chr1-1000-C-T": "C", "chr1-9999-A-G": "2"}, "edits": {}}}
    out = tmp_path / "out"
    report = rd.run(out_dir=out, refs=shared, meta_fn=lambda s: metas.get(s, {}), manual_fn=lambda ids: {},
                    pheno_fn=lambda s: {"GENE1": 10.0})

    assert report["cases"] == {"WES": 2}
    wes = report["by_test_type"]["WES"]
    assert wes["occurrences"] == 4
    # OTHER has a GeneBe score but no GeneBe class: like the platform, it falls through to the pipeline
    assert wes["source_before"]["genebe"] == 2 and wes["source_before"]["pipeline"] == 2
    assert wes["source_after"]["engine"] == 2 and wes["source_after"].get("pipeline", 0) == 0
    assert sum(wes["direction"].values()) == 2
    assert wes["tier_transition"].get("1C→2") == 1           # pipeline 5 points → engine below 4

    marked = report["marked"]["WES"]
    assert marked["n"] == 2                                   # the unknown variant is ignored
    assert marked["by_status"] == {"1": 1, "C": 1}
    assert marked["ranks"]["before"]["n"] == 2
    assert sum(marked["rank_move"]["A"].values()) == 2
    assert report["top20_overlap"]["WES"]["cases"] == 1

    with open(out / "marked_variants.tsv", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {r["variant"] for r in rows} == {"chr1-4000-T-C", "chr1-1000-C-T"}
    pipe = next(r for r in rows if r["variant"] == "chr1-4000-T-C")
    assert pipe["source_before"] == "pipeline" and pipe["source_after"] == "engine"
    assert pipe["tier_before"] == "1C" and pipe["tier_after"] == "2"
    lof = next(r for r in rows if r["variant"] == "chr1-1000-C-T")
    assert lof["source_after"] == "genebe"                    # A keeps GeneBe; B ranks by the engine
    with open(out / "tier_changes.tsv", encoding="utf-8") as handle:
        tiers = list(csv.DictReader(handle, delimiter="\t"))
    assert len(tiers) == 1 and "sample" not in tiers[0]
    text = (out / "report.md").read_text(encoding="utf-8")
    assert "24WE0001" not in text and "會出現在診斷報告" in text


def test_stale_results_are_skipped(tmp_path, monkeypatch):
    from app import config
    from app.services.acmg_engine import StaticReferences

    root = tmp_path / "root"
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", root)
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", tmp_path / "legacy_ui")
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", tmp_path / "legacy_pipe")
    publish(root, "24WE0001-nckuh", [LOF], refs())
    newer = StaticReferences(lof_genes={"GENE1": ("yes", "x")}, version={"refs_db": {"test": 2}})
    report = rd.run(out_dir=tmp_path / "out", refs=newer, meta_fn=lambda s: {}, manual_fn=lambda ids: {},
                    pheno_fn=lambda s: {})
    assert report["cases"] == {"skipped_stale": 1}
