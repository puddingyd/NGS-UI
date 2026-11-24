"""ACMG rule engine wiring (plan §5.2, step 4) — synthetic review TSVs only."""
import json
import os
import subprocess
import sys
from pathlib import Path

from app.services.acmg_engine import StaticReferences, runner
from app.services.acmg_engine import parse as P
from app.services.acmg_engine.recurrence import WesRecurrence, update_locked

REPO = Path(__file__).resolve().parents[1]
COLUMNS = ["CHROM", "POS", "REF", "ALT", "GENE", "TRANSCRIPT", "TRANSCRIPT_TYPE", "HGVS_C", "HGVS_P",
           "CONSEQUENCE", "IMPACT", "EXON", "INTRON", "PROTEIN_POSITION", "NMD", "LOFTEE", "LOFTEE_FILTER",
           "GNOMAD_E_AF", "GNOMAD_G_AF", "PANGOLIN_SCORE", "PHYLOP100", "CLINVAR_SIG", "CLINVAR_STARS", "GIAB_STRATA"]


def row(**values):
    base = dict.fromkeys(COLUMNS, ".")
    base.update(CHROM="chr1", POS="1000", REF="C", ALT="T", GENE="GENE1", TRANSCRIPT="ENST00000000001.1",
                TRANSCRIPT_TYPE="MANE_SELECT", HGVS_C="ENST00000000001.1:c.100C>T",
                HGVS_P="ENSP00000000001.1:p.Arg34Ter", CONSEQUENCE="stop_gained", IMPACT="HIGH",
                EXON="2/10", PROTEIN_POSITION="34/500", NMD="", LOFTEE="HC", PANGOLIN_SCORE="0.01",
                PHYLOP100="1.0", GIAB_STRATA="")
    base.update(values)
    return base


LOF = row()
LOF_OTHER_TX = row(TRANSCRIPT="ENST00000000009.1", TRANSCRIPT_TYPE=".", CONSEQUENCE="intron_variant",
                   IMPACT="MODIFIER", HGVS_P=".", LOFTEE=".")
MISSENSE = row(POS="2000", REF="G", ALT="A", CONSEQUENCE="missense_variant", IMPACT="MODERATE",
               HGVS_P="ENSP00000000001.1:p.Gly60Asp", LOFTEE=".")
MITO = row(CHROM="chrM", POS="73", REF="A", ALT="G")


def write_review(path: Path, rows):
    lines = ["\t".join(COLUMNS)] + ["\t".join(r.get(c, ".") for c in COLUMNS) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def refs():
    return StaticReferences(lof_genes={"GENE1": ("yes", "ClinGen HI=3")}, version={"refs_db": {"test": 1}})


def test_group_rows_and_ui_ids():
    groups = runner.group_variant_rows([LOF, LOF_OTHER_TX, MISSENSE, MITO, row(POS="1000.0")])
    assert list(groups) == ["chr1-1000-C-T", "chr1-2000-G-A"]
    assert len(groups["chr1-1000-C-T"]) == 3     # two transcripts + "1000.0" normalised to the same id


def test_build_case_and_lookup(tmp_path):
    review = write_review(tmp_path / "S1.snv_indel.review.tsv", [LOF, LOF_OTHER_TX, MISSENSE, MITO])
    out = tmp_path / "S1.snv_acmg_engine.sqlite"
    stats = runner.build_case(review_tsv=review, out=out, sample_id="S1", test_type="WES",
                              refs=refs(), recurrence=None, moi_fn=lambda gene: {"AR"})
    assert stats["variants"] == 2
    summaries = runner.lookup_many(out, ["chr1-1000-C-T", "chr1-2000-G-A", "chrX-1-A-G"])
    assert set(summaries) == {"chr1-1000-C-T", "chr1-2000-G-A"}
    lof = summaries["chr1-1000-C-T"]
    assert lof["criteria_text"].startswith("PVS1") and lof["classification"] == "Likely pathogenic"
    detail = runner.lookup_many(out, ["chr1-1000-C-T"], detail=True)["chr1-1000-C-T"]
    assert detail["transcript"].startswith("ENST00000000001")          # MANE row chosen
    assert any(e["code"] == "PVS1" and e["status"] == "met" for e in detail["evidence"])
    manifest = runner.read_manifest(out)
    assert manifest["variants"] == 2 and manifest["test_type"] == "WES"
    assert runner.is_current(out, review, refs())
    os.utime(review, ns=(1, 1_000_000_000))
    assert not runner.is_current(out, review, refs())                   # review rebuilt → stale


def test_review_signature_survives_copy_to_live_path(tmp_path):
    (tmp_path / "stage").mkdir()
    review = write_review(tmp_path / "stage" / "S1.snv_indel.review.tsv", [LOF])
    out = tmp_path / "stage" / "S1.snv_acmg_engine.sqlite"
    runner.build_case(review_tsv=review, out=out, sample_id="S1", test_type="WES", refs=refs(), recurrence=None,
                      moi_fn=None)
    import shutil
    live = tmp_path / "live"
    shutil.copytree(tmp_path / "stage", live)                           # copy2 keeps mtime
    assert runner.is_current(live / out.name, live / review.name, refs())


def test_recurrence_excludes_own_patient(tmp_path):
    store = WesRecurrence(tmp_path / "rec.sqlite")
    key = P.variant_key(LOF)
    store.upsert_sample("S1", "mrn:A", {}, [key])
    store.upsert_sample("S2", "mrn:B", {}, [key])
    store.upsert_sample("S3", "mrn:C", {}, [])
    assert store.lookup(key) == (2, 3)
    assert store.lookup(key, exclude_patient="mrn:A") == (1, 2)
    assert store.lookup(key, exclude_patient="mrn:Z") == (2, 3)       # patient not in the table
    case = runner.CaseRecurrence(refs(), store, "S1", runner.case_patient(store, "S1"))
    assert case.wes_recurrence(key) == (1, 2)
    assert "recurrence" in case.available and case.lof_mechanism("GENE1")[0] == "yes"
    store.close()


def test_attach_stored_and_on_demand(tmp_path, monkeypatch):
    review = write_review(tmp_path / "S1.snv_indel.review.tsv", [LOF])
    out = tmp_path / "S1.snv_acmg_engine.sqlite"
    shared = refs()
    monkeypatch.setattr(runner, "shared_references", lambda: shared)
    monkeypatch.setattr(runner, "shared_recurrence", lambda: None)
    runner.build_case(review_tsv=review, out=out, sample_id="S1", test_type="WES", refs=shared,
                      recurrence=None, moi_fn=None)
    variants = {"chr1-1000-C-T": {}, "chr1-2000-G-A": {}, "chr1-3000-A-C": {}}
    counts = runner.attach(variants, results_path=out, review_tsv=review, rows=[MISSENSE], sample_id="S1")
    assert counts == {"current": 1, "on_demand": 1}
    assert variants["chr1-1000-C-T"]["engine_acmg"]["status"] == "current"
    assert variants["chr1-2000-G-A"]["engine_acmg"]["status"] == "on_demand"
    assert "engine_acmg" not in variants["chr1-3000-A-C"]               # no rows, not stored
    again = runner.attach(variants, results_path=out, review_tsv=review, sample_id="S1")
    assert again == {}                                                   # already attached → untouched
    missing = {"chr1-1000-C-T": {}}
    runner.attach(missing, results_path=tmp_path / "absent.sqlite", review_tsv=review, sample_id="S1")
    assert missing == {"chr1-1000-C-T": {}}


def _env(tmp_path):
    return dict(os.environ, PYTHONPATH=str(REPO / "backend"), NGS_UI_DATA_ROOT=str(tmp_path / "data"),
                NGS_UI_TERTIARY_ROOT=str(tmp_path / "root"), NGS_UI_ACMG_REFS_DIR=str(tmp_path / "refs"),
                NGS_UI_ACMG_WES_RECURRENCE_DB=str(tmp_path / "data" / "rec.sqlite"),
                NGS_UI_PATIENT_LIST_DIR=str(tmp_path / "no_roster"))


def _cli(tmp_path, *args):
    # A fastapi that cannot be imported: the backfill must not need the web stack.
    blocker = tmp_path / "no_web"
    (blocker / "fastapi").mkdir(parents=True, exist_ok=True)
    (blocker / "fastapi" / "__init__.py").write_text("raise ImportError('fastapi blocked in test')\n")
    env = _env(tmp_path)
    env["PYTHONPATH"] = os.pathsep.join([str(blocker), env["PYTHONPATH"]])
    return subprocess.run([sys.executable, str(REPO / "scripts/acmg_engine/build_engine_results.py"), *map(str, args)],
                          capture_output=True, text=True, env=env)


def test_cli_staged_case(tmp_path):
    review = write_review(tmp_path / "S1.snv_indel.review.tsv", [LOF, MISSENSE])
    proc = _cli(tmp_path, "--review", review, "--out", tmp_path / "S1.snv_acmg_engine.sqlite",
                "--sample", "S1", "--test-type", "WES")
    assert proc.returncode == 0, proc.stderr
    assert "[acmg-engine] sample=S1 variants=2" in proc.stdout


def _published_case(root: Path, sid: str, rows, test_type="WES"):
    post = root / sid / "08_postprocessing"
    post.mkdir(parents=True)
    (post / f"{sid}.layout.json").write_text(json.dumps({"layout_version": 3}))
    write_review(post / f"{sid}.snv_indel.review.tsv", rows)
    (post / f"{sid}.snv_indel.review.tsv.source.json").write_text(json.dumps({"test_type": test_type}))
    return post


def test_cli_backfill_only_rebuilds_stale(tmp_path):
    root = tmp_path / "root"
    post = _published_case(root, "24WE0001-nckuh", [LOF])
    (root / "24WE0002-nckuh" / "08_postprocessing").mkdir(parents=True)    # no review TSV
    refused = _cli(tmp_path, "--all")                                     # no OMIM in this sandbox
    assert refused.returncode == 2 and "OMIM" in refused.stderr
    assert not (post / "24WE0001-nckuh.snv_acmg_engine.sqlite").exists()
    first = json.loads(_cli(tmp_path, "--all", "--allow-no-omim").stdout)
    assert first["built"] == 1 and first["skipped"] == {"沒有 review TSV": 1}
    assert first["omim"] == "unavailable"
    assert (post / "24WE0001-nckuh.snv_acmg_engine.sqlite").is_file()
    second = json.loads(_cli(tmp_path, "--all", "--allow-no-omim").stdout)
    assert second["built"] == 1                    # built without OMIM → never counted as current
    assert (tmp_path / "data" / "jobs" / "tertiary" / ".sample_locks" / "24WE0001-nckuh.lock").is_file()
    forced = json.loads(_cli(tmp_path, "--sample", "24WE0001-nckuh", "--force", "--allow-no-omim").stdout)
    assert forced["built"] == 1


def test_results_current_only_when_omim_was_read(tmp_path):
    review = write_review(tmp_path / "S1.snv_indel.review.tsv", [LOF])
    out = tmp_path / "S1.snv_acmg_engine.sqlite"
    runner.build_case(review_tsv=review, out=out, sample_id="S1", test_type="WES", refs=refs(),
                      recurrence=None, moi_fn=lambda gene: {"AD"})
    assert runner.read_manifest(out)["omim"] == "custom" and runner.is_current(out, review, refs())
    runner.write_results(out, iter(()), dict(runner.read_manifest(out), omim="unavailable"))
    assert not runner.is_current(out, review, refs())


def test_update_locked_adds_and_removes(tmp_path, monkeypatch):
    root = tmp_path / "root"
    _published_case(root, "24WE0001-nckuh", [LOF])
    db = tmp_path / "rec.sqlite"
    added = update_locked(["24WE0001-nckuh"], db=db, root=root, roster={})
    assert added["added"] == ["24WE0001-nckuh"] and added["patients"] == 1
    removed = update_locked([], remove=["24WE0001-nckuh"], db=db, root=root)
    assert removed["removed"] == ["24WE0001-nckuh"] and removed["patients"] == 0


def test_sample_loader_attaches_published_results(tmp_path, monkeypatch):
    from app import config
    from app.services import sample_loader

    root = tmp_path / "root"
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", root)
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", tmp_path / "legacy_ui")
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", tmp_path / "legacy_pipe")
    post = _published_case(root, "24WE0001-nckuh", [LOF])
    shared = refs()
    monkeypatch.setattr(runner, "shared_references", lambda: shared)
    monkeypatch.setattr(runner, "shared_recurrence", lambda: None)
    runner.build_case(review_tsv=post / "24WE0001-nckuh.snv_indel.review.tsv",
                      out=post / "24WE0001-nckuh.snv_acmg_engine.sqlite", sample_id="24WE0001-nckuh",
                      test_type="WES", refs=shared, recurrence=None, moi_fn=None)
    variants = {"chr1-1000-C-T": {}, "chr1-2000-G-A": {}}
    sample_loader._attach_engine_acmg(variants, "24WE0001-nckuh", rows=[MISSENSE], test_type="WES")
    assert variants["chr1-1000-C-T"]["engine_acmg"]["status"] == "current"
    assert variants["chr1-2000-G-A"]["engine_acmg"]["status"] == "on_demand"
    # step 7: the engine is the "In-house" ACMG slot
    lof = variants["chr1-1000-C-T"]
    assert lof["inhouse_acmg_source"] == "engine"
    assert lof["ACMG_score"] == lof["engine_acmg"]["points"]
    assert lof["ACMG_classification"] in ("Pathogenic", "Likely pathogenic", "Uncertain significance",
                                          "Likely benign", "Benign")
    assert lof["geno_score"] is not None
    broken = {"chr1-1000-C-T": {}}
    monkeypatch.setattr(runner, "shared_references", lambda: (_ for _ in ()).throw(OSError("boom")))
    sample_loader._attach_engine_acmg(broken, "24WE0001-nckuh")      # never raises into loading
    assert broken == {"chr1-1000-C-T": {}}


def test_code_signature_ignores_tooling_files():
    names = {p.name for p in runner._PACKAGE_DIR.iterdir() if p.is_file()}
    assert {"engine.py", "rulebook_v1.json"} <= names
    hashed = {n for n in names if n in runner._RESULT_CODE_FILES or n.endswith((".json", ".tsv"))}
    assert "validation.py" not in hashed and "runner.py" not in hashed and "recurrence.py" not in hashed


def test_effective_acmg_prefers_inhouse_engine_over_genebe(monkeypatch):
    """Step 7: reviewer → In-house (rule engine) → GeneBe; ranking uses the same points."""
    from app.services import manual_acmg, sample_loader

    monkeypatch.setattr(manual_acmg, "bulk_current", lambda *a, **k: {})
    monkeypatch.setattr(manual_acmg, "bulk_observed_counts", lambda *a, **k: {})
    variants = {
        "chr1-1-A-G": {"ACMG_classification": "Uncertain significance", "ACMG_score": 2, "tier": "2",
                       "genebe_acmg_class": "Likely pathogenic", "genebe_acmg_score": 7, "pheno_score": 10,
                       "engine_acmg": {"classification": "Likely pathogenic", "points": 8,
                                       "criteria_text": "PVS1_Strong,PM2_Supporting,PP3_Supporting"}},
        "chr1-2-A-G": {"ACMG_classification": "", "ACMG_score": None, "tier": "2",
                       "genebe_acmg_class": "Likely benign", "genebe_acmg_score": -2},
    }
    sample_loader._promote_engine_to_inhouse(variants)
    first = variants["chr1-1-A-G"]
    assert first["pipeline_acmg_class"] == "Uncertain significance" and first["ACMG_score"] == 8
    assert first["tier"] == "1C"
    categories = {"1A": [], "1B": [], "1C": [], "2": list(variants)}
    sample_loader._apply_effective_acmg(variants, categories, sample_id="S", genome_build="hg38", meta={})
    assert variants["chr1-1-A-G"]["effective_acmg_source"] == "in-house"
    assert variants["chr1-1-A-G"]["effective_acmg_score"] == 8
    assert variants["chr1-1-A-G"]["total_score"] == manual_acmg.acmg_to_variant_score(8) + 10
    assert variants["chr1-2-A-G"]["effective_acmg_source"] == "GeneBe"     # no in-house value → GeneBe
    sample_loader._promote_engine_to_inhouse(variants)                      # idempotent
    assert variants["chr1-1-A-G"]["pipeline_acmg_class"] == "Uncertain significance"
