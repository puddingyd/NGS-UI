"""Offline validation of the ACMG rule engine (plan §5.2.10) — synthetic cases only."""
import gzip
import json
import sqlite3
from pathlib import Path

from app.services.acmg_engine import StaticReferences, runner, validation
from app.services.acmg_engine.recurrence import WesRecurrence

COLUMNS = ["CHROM", "POS", "REF", "ALT", "GENE", "TRANSCRIPT", "TRANSCRIPT_TYPE", "HGVS_C", "HGVS_P",
           "CONSEQUENCE", "IMPACT", "EXON", "INTRON", "PROTEIN_POSITION", "NMD", "LOFTEE", "LOFTEE_FILTER",
           "GNOMAD_E_AF", "GNOMAD_G_AF", "PANGOLIN_SCORE", "PHYLOP100", "CLINVAR_SIG", "CLINVAR_STARS",
           "GIAB_STRATA", "ACMG_SCORE", "ACMG_CLASS", "ACMG_CRITERIA", "GENEBE_ACMG_SCORE",
           "GENEBE_ACMG_CLASS", "GENEBE_ACMG_CRITERIA", "CLINGEN_VCEP_CRITERIA", "CLINGEN_VCEP_CLASS"]


def row(**values):
    base = dict.fromkeys(COLUMNS, ".")
    base.update(CHROM="chr1", POS="1000", REF="C", ALT="T", GENE="GENE1", TRANSCRIPT="ENST00000000001.1",
                TRANSCRIPT_TYPE="MANE_SELECT", HGVS_C="ENST00000000001.1:c.100C>T",
                HGVS_P="ENSP00000000001.1:p.Arg34Ter", CONSEQUENCE="stop_gained", IMPACT="HIGH",
                EXON="2/10", PROTEIN_POSITION="34/500", NMD="", LOFTEE="HC", PANGOLIN_SCORE="0.01",
                PHYLOP100="1.0", GIAB_STRATA="")
    base.update(values)
    return base


LOF = row(CLINVAR_SIG="Pathogenic", CLINVAR_STARS="3", ACMG_SCORE="10", ACMG_CLASS="Pathogenic",
          ACMG_CRITERIA="PVS1,PM2_Supporting", GENEBE_ACMG_SCORE="9", GENEBE_ACMG_CLASS="Likely_pathogenic",
          GENEBE_ACMG_CRITERIA="PVS1,PM2", CLINGEN_VCEP_CRITERIA="PVS1,PM2_Supporting",
          CLINGEN_VCEP_CLASS="Pathogenic")
BENIGN = row(POS="3000", REF="A", ALT="G", CONSEQUENCE="synonymous_variant", IMPACT="LOW",
             HGVS_P="ENSP00000000001.1:p.Gly80=", LOFTEE=".", GNOMAD_E_AF="0.2",
             CLINVAR_SIG="Benign", CLINVAR_STARS="2", ACMG_SCORE="-8", ACMG_CLASS="Benign",
             GENEBE_ACMG_SCORE="-5", GENEBE_ACMG_CLASS="Likely_benign")
OTHER = row(POS="2000", REF="G", ALT="A", CONSEQUENCE="missense_variant", IMPACT="MODERATE",
            HGVS_P="ENSP00000000001.1:p.Gly60Asp", LOFTEE=".", ACMG_SCORE="1", GENEBE_ACMG_SCORE="2")


def write_review(path: Path, rows):
    lines = ["\t".join(COLUMNS)] + ["\t".join(r.get(c, ".") for c in COLUMNS) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def refs():
    return StaticReferences(lof_genes={"GENE1": ("yes", "ClinGen HI=3")}, version={"refs_db": {"test": 1}})


def publish(root: Path, sid: str, rows, shared):
    post = root / sid / "08_postprocessing"
    post.mkdir(parents=True)
    (post / f"{sid}.layout.json").write_text(json.dumps({"layout_version": 3}))
    review = write_review(post / f"{sid}.snv_indel.review.tsv", rows)
    (post / f"{sid}.snv_indel.review.tsv.source.json").write_text(json.dumps({"test_type": "WES"}))
    runner.build_case(review_tsv=review, out=post / f"{sid}.snv_acmg_engine.sqlite", sample_id=sid,
                      test_type="WES", refs=shared, recurrence=None, moi_fn=lambda gene: {"AR"})


def manual_db(path: Path, sid: str):
    conn = sqlite3.connect(path)
    conn.executescript(
        """
        CREATE TABLE manual_revisions (id INTEGER PRIMARY KEY, genome_build TEXT, variant_id TEXT,
          criteria_json TEXT, criteria_text TEXT, acmg_score REAL, acmg_classification TEXT,
          reusable_criteria_json TEXT, reusable_criteria_text TEXT, reusable_acmg_score REAL,
          reusable_acmg_classification TEXT, reviewer_user_id INTEGER, reviewer_username TEXT,
          source_sample_id TEXT, created_at TEXT);
        CREATE TABLE manual_current (genome_build TEXT, variant_id TEXT, revision_id INTEGER);
        CREATE TABLE observations (genome_build TEXT, variant_id TEXT, sample_id TEXT, status TEXT,
          reviewer_user_id INTEGER, reviewer_username TEXT, updated_at TEXT);
        """)
    crit = json.dumps({"PVS1": {"enabled": True, "strength": "very_strong"},
                       "PS3": {"enabled": True, "strength": "strong"}})
    conn.execute("INSERT INTO manual_revisions VALUES (1,'hg38','chr1-1000-C-T',?, 'PVS1,PS3', 12, "
                 "'Pathogenic', ?, 'PVS1,PS3', 12, 'Pathogenic', 1, 'reviewer', ?, '2026-11-01')",
                 (crit, crit, sid))
    conn.execute("INSERT INTO manual_current VALUES ('hg38','chr1-1000-C-T',1)")
    conn.execute("INSERT INTO observations VALUES ('hg38','chr1-1000-C-T',?, '1', 1, 'reviewer', 'x')", (sid,))
    conn.commit()
    conn.close()


def test_parse_criteria_and_auc():
    assert validation.parse_criteria("PVS1,PM2_Supporting;PP3 (Strong)") == {
        "PVS1": "very_strong", "PM2": "supporting", "PP3": "strong"}
    assert validation.parse_criteria("PS3_Moderate,BA1") == {"PS3": "moderate", "BA1": "stand_alone"}
    assert validation.auc([3, 4], [1, 2]) == 1.0
    assert validation.auc([1, 2], [1, 2]) == 0.5
    assert validation.auc([], [1]) is None
    assert validation.geno_score(None) == 0 and validation.geno_score(10) == 100 and validation.geno_score(-20) == 0


def test_run_reports_every_section(tmp_path, monkeypatch):
    from app import config

    root = tmp_path / "root"
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", root)
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", tmp_path / "legacy_ui")
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", tmp_path / "legacy_pipe")
    shared = refs()
    publish(root, "24WE0001-nckuh", [LOF, OTHER, BENIGN], shared)
    publish(root, "24WE0002-nckuh", [LOF, OTHER], shared)
    db = tmp_path / "manual.sqlite"
    manual_db(db, "24WE0001-nckuh")
    rec_db = tmp_path / "rec.sqlite"
    store = WesRecurrence(rec_db)
    store.upsert_sample("24WE0001-nckuh", "mrn:A", {}, ["chr1-2000-G-A"])
    store.upsert_sample("24WE0002-nckuh", "mrn:B", {}, ["chr1-2000-G-A"])
    store.close()
    old = tmp_path / "variant_summary_t0.txt.gz"
    with gzip.open(old, "wt") as handle:
        handle.write("#AlleleID\tAssembly\tClinicalSignificance\tReviewStatus\tChromosome\tPositionVCF\t"
                     "ReferenceAlleleVCF\tAlternateAlleleVCF\n")
        handle.write("1\tGRCh38\tBenign\tcriteria provided, multiple submitters, no conflicts\t1\t3000\tA\tG\n")

    out = tmp_path / "out"
    report = validation.run(out_dir=out, refs=shared, recurrence=None, recurrence_db=rec_db, manual_db=db,
                            t0_clinvar=old, moi_fn=lambda gene: {"AR"},
                            pheno_fn=lambda sid: {"GENE1": 50.0},
                            rule_variants={"pp2_loose": {"pp2": {"min_plp": 0, "max_benign_fraction": 1.0}}})

    rv = report["rule_variants"]
    assert set(rv["subsets"]["all"]) == {"current", "pp2_loose"}
    cur = rv["subsets"]["all"]["current"]
    assert cur["plp"] == 1 and cur["blb"] == 1 and set(cur["criteria"]) == {"PM1", "PP2"}
    assert cur["criteria"]["PP2"]["plp_n"] == 1
    assert "## 7." in validation.render_markdown(report)
    assert report["cases"] == {"WES": 2}
    labelled = report["clinvar"]["all"]["all"]
    assert labelled["plp"] == 1 and labelled["blb"] == 1
    assert labelled["sources"]["engine"]["auc"] == 1.0
    new = report["clinvar"]["new_since_t0"]["all"]
    assert new["plp"] == 1 and new["blb"] == 0                     # benign already labelled at T0
    assert report["vcep"]["n"] == 1 and report["vcep"]["criteria"]["PVS1"]["tp"] == 1
    assert report["manual"]["counts"]["matched"] == 1
    assert report["manual"]["criteria"]["PVS1"]["tp"] == 1
    assert report["manual"]["reference_only_criteria"] == {"PS3": 1}
    ranking = report["ranking"]["by_test_type"]["WES"]["with_pheno"]["engine"]
    assert ranking["n"] == 1 and ranking["top1"] == 1.0
    lof = report["coverage"]["WES"]["lof_hc_pvs1"]["recessive"]
    assert lof["n"] == 2 and lof["engine"] == 1.0 and lof["pipeline"] == 1.0
    breakdown = report["ranking_breakdown"]["WES"]
    assert breakdown["n"] == 1 and breakdown["engine_better"] + breakdown["same"] + breakdown["engine_worse"] == 1
    assert report["criteria_usage"]["genebe"]["PLP"]["criteria"]["PVS1"] == 1.0
    assert report["vcep"]["frequency_crosstab"] == {}
    presence = report["coverage"]["WES"]["predictor_presence"]["truncating"]
    assert presence["n"] == 2 and presence["PANGOLIN_SCORE"] == 1.0
    assert report["coverage"]["WES"]["columns_absent"]["SPLICEAI_MAX"] > 0
    assert report["recurrence"]["wes_patients"] == 2
    assert report["recurrence"]["wes_thresholds"]["5%"]["variants"] == 1

    for name in ("report.md", "summary.json", "clinvar_labelled.tsv", "recurrent_flagged.tsv",
                 "wes_recurrence_candidates.tsv"):
        text = (out / name).read_text(encoding="utf-8")
        assert "24WE000" not in text and "mrn:" not in text and "reviewer" not in text, name
    assert "## 4." in (out / "report.md").read_text(encoding="utf-8")
