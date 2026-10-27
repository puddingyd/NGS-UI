import csv
import fcntl
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

from app.adapters.snv_tsv import _row_to_variant, classify_tier
from app.services import gpn_msa, gpn_star, snv_review

REPO = Path(__file__).resolve().parents[1]
FIELDS = ["CHROM", "POS", "REF", "ALT", "GENE", "CLINVAR_SIG", "GNOMAD_G_AF", "DP", "GPN_MSA_SCORE"]


def row(pos=10, ref="A", alt="G", **values):
    return {"CHROM": "chr1", "POS": str(pos), "REF": ref, "ALT": alt,
            "GENE": "KEEP", "CLINVAR_SIG": "", "GNOMAD_G_AF": "0.001",
            "DP": "30", "GPN_MSA_SCORE": "-9", **values}


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def make_bundle(root, model=gpn_star.MODEL):
    entries = []
    for chrom in sorted(gpn_star.PRIMARY_CONTIGS):
        path = gpn_star.shard_path(root, chrom, model)
        path.parent.mkdir(parents=True, exist_ok=True)
        pl.DataFrame({"chrom": [chrom] * 5, "pos": [10, 10, 11, 12, 13],
                      "ref": ["A"] * 5, "alt": ["G", "C", "G", "G", "G"],
                      "llr_calibrated": [-1.2345678 if model == gpn_star.MODEL else -7.0, 2.5, 0.0, float("nan"), None]},
                     schema={"chrom": pl.String, "pos": pl.Int64, "ref": pl.String,
                             "alt": pl.String, "llr_calibrated": pl.Float32}).write_parquet(path, row_group_size=2)
        entries.append({"chrom": chrom, "path": path.relative_to(root).as_posix(),
                        "size": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "rows": 5, "assembly": "hg38", "score_set": model, "score_type": "llr"})
    gpn_star.manifest_path(root, model).write_text(json.dumps({"schema_version": 1,
        "model": model, "assembly": "hg38", "score_field": "llr_calibrated",
        "revision": "a" * 40, "files": entries}))
    return root


@pytest.fixture
def bundle(tmp_path):
    return make_bundle(tmp_path / "gpn_star")


@pytest.fixture
def dual_bundle(bundle):
    make_bundle(bundle, gpn_star.V_MODEL)
    return bundle


def test_real_parquet_exact_alleles_precision_zero_missing_and_unsupported(tmp_path, bundle):
    review = tmp_path / "review.tsv"
    rows = [row(), row(GENE="OTHER_TRANSCRIPT"), row(alt="C"), row(ref="T"),
            row(pos=11), row(pos=12), row(pos=13), row(pos=30, alt="AT"),
            row(ref="AA", alt="GG"), row(alt="*"), row(CHROM="chrM"),
            row(CHROM="chrX"), row(CHROM="chrY"), row(alt="A")]
    write_tsv(review, rows)
    stats = gpn_star.annotate_review_tsv(review, bundle, required=True)
    actual = read_tsv(review)
    assert stats["status"] == "complete"
    assert stats["annotated_rows"] == 6
    assert stats["unsupported_rows"] == 5
    assert stats["unmatched_variants"] == 3
    assert float(actual[0]["GPN_STAR_SCORE"]) == pytest.approx(-1.2345678)
    assert actual[0]["GPN_STAR_SCORE"] == actual[1]["GPN_STAR_SCORE"]
    assert actual[2]["GPN_STAR_SCORE"] == "2.5"
    assert actual[3]["GPN_STAR_SCORE"] == ""
    assert actual[4]["GPN_STAR_SCORE"] == "0"
    assert all(r["GPN_STAR_SCORE"] == "" for r in actual[5:11])
    assert [{k: r[k] for k in FIELDS} for r in actual] == rows
    assert _row_to_variant(actual[0])["GPN_STAR_score"] == pytest.approx(-1.2345678)
    assert _row_to_variant(actual[4])["GPN_STAR_score"] == 0
    assert classify_tier({"GPN_STAR_SCORE": "-99"}) == "2"


@pytest.mark.parametrize("failure,expected", [("db", "skipped_missing_db"),
    ("shard", "skipped_missing_shard"), ("schema", "failed"),
    ("model", "failed"), ("polars", "skipped_missing_polars"), ("query", "failed")])
def test_failure_preserves_review(tmp_path, bundle, monkeypatch, failure, expected):
    review = tmp_path / "review.tsv"
    write_tsv(review, [row()])
    before = review.read_bytes()
    if failure == "db":
        bundle = tmp_path / "absent"
    elif failure == "shard":
        gpn_star.shard_path(bundle, "Y").unlink()
    elif failure == "model":
        marker = bundle / "manifest.json"
        payload = json.loads(marker.read_text()); payload["model"] = "wrong-model"
        marker.write_text(json.dumps(payload))
    elif failure == "schema":
        path = gpn_star.shard_path(bundle, "1")
        path.write_bytes(b"broken parquet")
    else:
        def fail(*args):
            raise ImportError("polars missing") if failure == "polars" else RuntimeError("query failed")
        monkeypatch.setattr(gpn_star, "query_scores", fail)
    stats = gpn_star.annotate_review_tsv(review, bundle)
    assert stats["status"] == expected
    assert review.read_bytes() == before
    assert not list(tmp_path.glob(".review.tsv.*.tmp"))
    with pytest.raises(Exception):
        gpn_star.annotate_review_tsv(review, bundle, required=True)


def test_only_retained_snvs_queried_and_msa_manifest_migrates_without_raw_scan(tmp_path, bundle, monkeypatch):
    raw = tmp_path / "raw.tsv"
    write_tsv(raw, [row(), row(pos=20, GNOMAD_G_AF="0.2"), row(pos=30, alt="AT", CLINVAR_SIG="Pathogenic")])
    before = raw.read_bytes()
    monkeypatch.setattr(snv_review, "_candidate_bed_path", lambda: tmp_path / "no.bed")
    captured = []
    query = gpn_star.query_scores
    def spy(root, keys, model):
        captured.append(keys)
        return query(root, keys, model)
    monkeypatch.setattr(gpn_star, "query_scores", spy)
    review = snv_review.ensure_review_tsv(raw, test_type="WGS", output_dir=tmp_path / "post", gpn_star_dir=tmp_path / "absent", research_only=True)
    marker = review.with_suffix(review.suffix + ".source.json")
    old = json.loads(marker.read_text())
    old.pop("gpn_star"); old.pop("gpn_star_annotation")
    marker.write_text(json.dumps(old))
    def no_scan(*args):
        raise AssertionError("raw/filter rebuild unexpected")
    monkeypatch.setattr(snv_review, "_load_bed", no_scan)
    snv_review.ensure_review_tsv(raw, test_type="WGS", output_dir=review.parent, gpn_star_dir=bundle)
    assert captured == [{("1", 10, "A", "G")}]
    rows = read_tsv(review)
    assert [r["POS"] for r in rows] == ["10", "30"]
    assert rows[0]["GPN_STAR_SCORE"] and rows[1]["GPN_STAR_SCORE"] == ""
    assert raw.read_bytes() == before
    payload = json.loads(marker.read_text())
    assert payload["gpn_star_annotation"]["status"] == "complete"
    assert payload["gpn_msa"] == old["gpn_msa"]
    # The unchanged source is a cache hit, not another scan/query.
    snv_review.ensure_review_tsv(raw, test_type="WGS", output_dir=review.parent, gpn_star_dir=bundle)
    assert len(captured) == 1
    with (bundle / "manifest.json").open() as handle:
        updated = json.load(handle)
    updated["revision"] = "b" * 40
    (bundle / "manifest.json").write_text(json.dumps(updated))
    snv_review.ensure_review_tsv(raw, test_type="WGS", output_dir=review.parent, gpn_star_dir=bundle)
    assert len(captured) == 2
    assert json.loads(marker.read_text())["gpn_star_annotation"]["revision"] == "b" * 40
    previous_review = review.read_bytes()
    updated["revision"] = "c" * 40
    updated["model"] = "wrong-model"
    (bundle / "manifest.json").write_text(json.dumps(updated))
    snv_review.ensure_review_tsv(raw, test_type="WGS", output_dir=review.parent, gpn_star_dir=bundle)
    failed = json.loads(marker.read_text())
    assert failed["gpn_star_annotation"]["status"] == "failed"
    assert failed["gpn_star_values"]["revision"] == "b" * 40
    assert review.read_bytes() == previous_review


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def msa_db(tmp_path, monkeypatch):
    path = tmp_path / "msa.bgz"
    path.touch()
    path.with_suffix(".bgz.tbi").touch()
    monkeypatch.setattr(gpn_msa, "validate_database", lambda _db: (path, "tabix"))
    calls = []
    def query(_db, keys, **kwargs):
        calls.append(keys)
        return {key: "-8.25" for key in keys}
    monkeypatch.setattr(gpn_msa, "query_scores", query)
    return path, calls


def test_standard_review_runs_msa_without_touching_star_resources(tmp_path, msa_db, monkeypatch):
    raw = tmp_path / "raw.tsv"
    write_tsv(raw, [row(), row(pos=20, GNOMAD_G_AF="0.2")])
    monkeypatch.setattr(snv_review, "_candidate_bed_path", lambda: tmp_path / "no.bed")
    def forbidden(*args, **kwargs):
        raise AssertionError("standard case must not inspect or query Star")
    monkeypatch.setattr(gpn_star, "database_signature", forbidden)
    monkeypatch.setattr(gpn_star, "annotate_review_tsv", forbidden)
    args = dict(test_type="WGS", output_dir=tmp_path / "post", gpn_msa_db=msa_db[0])
    review = snv_review.ensure_review_tsv(raw, **args)
    assert msa_db[1] == [{("1", 10, "A", "G")}]
    assert read_tsv(review)[0]["GPN_MSA_SCORE"] == "-8.25"
    assert "GPN_STAR_SCORE" not in read_tsv(review)[0]
    payload = json.loads(review.with_suffix(review.suffix + ".source.json").read_text())
    assert payload["gpn_star_research_only"] is False
    assert payload["gpn_star_annotation"]["status"] == "skipped_not_research_only"
    snv_review.ensure_review_tsv(raw, **args)
    assert len(msa_db[1]) == 1
    with pytest.raises(ValueError, match="Research-only"):
        snv_review.ensure_review_tsv(raw, **args, require_gpn_star=True)


def test_research_both_scores_and_independent_cache_refresh(tmp_path, bundle, msa_db, monkeypatch):
    raw = tmp_path / "raw.tsv"
    write_tsv(raw, [row()])
    monkeypatch.setattr(snv_review, "_candidate_bed_path", lambda: tmp_path / "no.bed")
    star_calls = []
    query = gpn_star.query_scores
    def spy(*args):
        star_calls.append(1)
        return query(*args)
    monkeypatch.setattr(gpn_star, "query_scores", spy)
    args = dict(test_type="WGS", output_dir=tmp_path / "post", gpn_msa_db=msa_db[0], gpn_star_dir=bundle)
    review = snv_review.ensure_review_tsv(raw, **args, research_only=True)
    assert read_tsv(review)[0]["GPN_MSA_SCORE"] == "-8.25"
    assert float(read_tsv(review)[0]["GPN_STAR_SCORE"]) == pytest.approx(-1.2345678)
    # Loader calls inherit the explicit mode; neither provider reruns on a hit.
    snv_review.ensure_review_tsv(raw, **args)
    assert len(msa_db[1]) == len(star_calls) == 1
    msa_db[0].write_text("new msa")
    snv_review.ensure_review_tsv(raw, **args)
    assert len(msa_db[1]) == 2 and len(star_calls) == 1
    manifest = bundle / "manifest.json"
    payload = json.loads(manifest.read_text()); payload["revision"] = "b" * 40
    manifest.write_text(json.dumps(payload))
    snv_review.ensure_review_tsv(raw, **args)
    assert len(msa_db[1]) == len(star_calls) == 2
    # Explicitly disabling future annotation preserves existing displayed values.
    before = review.read_bytes()
    snv_review.ensure_review_tsv(raw, **args, research_only=False)
    snv_review.ensure_review_tsv(raw, **args)
    assert review.read_bytes() == before
    assert len(star_calls) == 2
    marker = review.with_suffix(review.suffix + ".source.json")
    saved = json.loads(marker.read_text())
    assert saved["gpn_star_research_only"] is False
    assert saved["gpn_star_values"]["revision"] == "b" * 40
    # A pre-gating Star manifest must not accidentally opt into new research.
    saved.pop("gpn_star_research_only")
    marker.write_text(json.dumps(saved))
    snv_review.ensure_review_tsv(raw, **args)
    assert len(star_calls) == 2


def test_missing_star_does_not_remove_msa(tmp_path, msa_db, monkeypatch):
    raw = tmp_path / "raw.tsv"
    write_tsv(raw, [row()])
    monkeypatch.setattr(snv_review, "_candidate_bed_path", lambda: tmp_path / "no.bed")
    review = snv_review.ensure_review_tsv(raw, test_type="WGS", output_dir=tmp_path / "post",
        gpn_msa_db=msa_db[0], gpn_star_dir=tmp_path / "absent", research_only=True)
    assert read_tsv(review)[0]["GPN_MSA_SCORE"] == "-8.25"
    payload = json.loads(review.with_suffix(review.suffix + ".source.json").read_text())
    assert payload["gpn_msa_annotation"]["status"] == "complete"
    assert payload["gpn_star_annotation"]["status"] == "skipped_missing_db"


def test_review_cli_requires_research_flag_for_star_and_restores_msa_option(tmp_path, bundle):
    raw = tmp_path / "raw.tsv"
    write_tsv(raw, [row()])
    command = [sys.executable, "scripts/build_snv_review_tsv.py", "--tsv", str(raw),
               "--output-dir", str(tmp_path / "post"), "--test-type", "WGS",
               "--gpn-msa-db", str(tmp_path / "missing.bgz"), "--gpn-star-dir", str(bundle)]
    env = {**os.environ, "NGS_UI_CDS_CANDIDATE_BED": str(tmp_path / "no.bed")}
    ordinary = subprocess.run(command, cwd=REPO, env=env, capture_output=True, text=True)
    assert ordinary.returncode == 0, ordinary.stderr
    assert "skipped_not_research_only" in ordinary.stdout
    review = tmp_path / "post" / "snv_indel.review.tsv"
    assert "GPN_STAR_SCORE" not in read_tsv(review)[0]
    research = subprocess.run(command + ["--research-only"], cwd=REPO, env=env, capture_output=True, text=True)
    assert research.returncode == 0, research.stderr
    assert read_tsv(review)[0]["GPN_STAR_SCORE"]


@pytest.mark.parametrize("research", [False, True])
def test_stopgaps_passes_research_mode_to_actual_review_builder(tmp_path, dual_bundle, research):
    bundle = dual_bundle
    raw = tmp_path / "raw.tsv"
    work = tmp_path / "work.tsv"
    write_tsv(raw, [row()]); work.write_bytes(raw.read_bytes())
    original = raw.read_bytes()
    msa = tmp_path / "msa.bgz"
    msa.touch(); msa.with_suffix(".bgz.tbi").touch()
    tabix = tmp_path / "tabix"
    tabix.write_text("#!/bin/sh\nprintf '1\\t10\\tA\\tG\\t-8.25\\n'\n")
    tabix.chmod(0o755)
    post = tmp_path / "08_postprocessing"
    command = ["bash", "scripts/run_stopgaps.sh", "--work-tsv", str(work),
               "--raw-tsv", str(raw), "--post-dir", str(post), "--sample", "CASE",
               "--seq-type", "WGS", "--gpn-msa-db", str(msa), "--gpn-star-dir", str(bundle),
               "--skip-spliceai", "--skip-clinvar-latest", "--skip-genebe",
               "--skip-giab", "--skip-inhouse-af", "--skip-mane-refseq", "--skip-litvar2",
               "--skip-cnv", "--skip-candidate-bed"]
    if research:
        command.append("--research-only")
    env = {**os.environ, "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"],
           "TABIX_BIN": str(tabix), "NGS_UI_DATA_ROOT": str(tmp_path / "data")}
    proc = subprocess.run(command, cwd=REPO, env=env, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    actual = read_tsv(post / "CASE.snv_indel.review.tsv")[0]
    assert actual["GPN_MSA_SCORE"] == "-8.25"
    assert bool(actual.get("GPN_STAR_SCORE")) is research
    assert bool(actual.get("GPN_STAR_V_SCORE")) is research
    if research:
        assert actual["GPN_STAR_V_SCORE"] == "-7"
    payload = json.loads((post / "CASE.snv_indel.review.tsv.source.json").read_text())
    assert payload["gpn_star_research_only"] is research
    assert raw.read_bytes() == original


def test_backfill_dry_run_apply_preserves_state_and_stale_signatures(tmp_path, dual_bundle):
    bundle = dual_bundle
    root = tmp_path / "tertiary"
    post = root / "CASE-nckuh" / "08_postprocessing"
    review = post / "CASE-nckuh.snv_indel.review.tsv"
    write_tsv(review, [row(), row(pos=30, alt="AT")])
    original = review.read_bytes()
    marker = review.with_suffix(review.suffix + ".source.json")
    old = {"raw_mtime_ns": 123, "overlay": {"old": True}, "candidate_bed": {"stale": True}, "gpn_msa": {"old": True}}
    marker.write_text(json.dumps(old))
    state = post / "CASE-nckuh.sample_metadata.json"
    state.write_text('{"reviewer": "preserve"}')
    command = [sys.executable, "scripts/backfill_gpn_star.py", "--root", str(root),
               "--gpn-star-dir", str(bundle), "CASE", "--research-only"]
    env = {**os.environ, "NGS_UI_DATA_ROOT": str(tmp_path / "data")}
    rejected = subprocess.run(command[:-1] + ["--apply"], cwd=REPO, env=env, capture_output=True, text=True)
    assert rejected.returncode == 2
    assert "requires --research-only" in rejected.stderr
    assert review.read_bytes() == original
    dry = subprocess.run(command, cwd=REPO, env=env, capture_output=True, text=True)
    assert dry.returncode == 0, dry.stderr
    assert '"status": "dry_run"' in dry.stdout
    assert review.read_bytes() == original
    assert not list((tmp_path / "data").rglob("*.lock"))
    lock_path = tmp_path / "data" / "jobs" / "tertiary" / ".sample_locks" / "CASE-nckuh.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        blocked = subprocess.run(command + ["--apply"], cwd=REPO, env=env, capture_output=True, text=True)
        assert blocked.returncode == 1
        assert "tertiary analysis is already running" in blocked.stdout
        assert review.read_bytes() == original
    applied = subprocess.run(command + ["--apply"], cwd=REPO, env=env, capture_output=True, text=True)
    assert applied.returncode == 0, applied.stderr
    assert read_tsv(review)[0]["GPN_STAR_SCORE"]
    assert read_tsv(review)[0]["GPN_STAR_V_SCORE"] == "-7"
    assert json.loads(marker.read_text())["gpn_star_research_only"] is True
    assert all(json.loads(marker.read_text())[k] == v for k, v in old.items())
    assert state.read_text() == '{"reviewer": "preserve"}'
    # Partial backfill reports failure explicitly but keeps successful old V values.
    gpn_star.manifest_path(bundle, gpn_star.V_MODEL).unlink()
    before = read_tsv(review)
    partial = subprocess.run(command + ["--apply"], cwd=REPO, env=env, capture_output=True, text=True)
    assert partial.returncode == 1
    result = json.loads(partial.stdout)
    assert result["status"] == "incomplete"
    assert result["models"]["gpn_star"]["status"] == "complete"
    assert result["models"]["gpn_star_v"]["status"] == "skipped_missing_db"
    assert read_tsv(review) == before
    assert json.loads(marker.read_text())["gpn_star_v_values"]["revision"] == "a" * 40


def test_backfill_rolls_back_when_manifest_publish_fails(tmp_path, monkeypatch):
    mod = load_script("backfill_gpn_star")
    review, staged, marker = [tmp_path / n for n in ("review.tsv", "new.tsv", "manifest.json")]
    review.write_text("old"); staged.write_text("new"); marker.write_text("old marker")
    replace = os.replace
    def fail(source, destination):
        if Path(destination) == marker:
            raise OSError("simulated disk failure")
        replace(source, destination)
    monkeypatch.setattr(mod.os, "replace", fail)
    with pytest.raises(OSError):
        mod.publish_pair(review, marker, staged, {})
    assert review.read_text() == "old"
    assert marker.read_text() == "old marker"


def test_download_verification_and_activation(tmp_path, bundle, monkeypatch):
    mod = load_script("download_gpn_star_scores")
    payload = json.loads((bundle / "manifest.json").read_text())
    files = [{**entry, "path": "data/" + entry["path"]} for entry in payload["files"]]
    release = {"repository": {"id": mod.REPOSITORY, "license": "apache-2.0"}, "parquet": {"files": files}}
    assert len(mod.select_files(release)) == 24
    dest = tmp_path / "install"
    def fake_curl(command, **kwargs):
        target = Path(command[command.index("--output") + 1])
        source = bundle / command[-1].split("/data/", 1)[1]
        target.write_bytes(source.read_bytes())
    monkeypatch.setattr(mod.subprocess, "run", fake_curl)
    mod.install(dest, release, "a" * 40, "https://example.invalid")
    assert gpn_star.validate_database(dest)["revision"] == "a" * 40
    mod.install(dest, release, "a" * 40, "https://example.invalid")  # restart/idempotence
    gpn_star.shard_path(dest, "1").write_bytes(b"corrupt")
    with pytest.raises(ValueError):
        mod.verify_shard(gpn_star.shard_path(dest, "1"), files[0])


def test_interrupted_download_resumes_before_activation(tmp_path, bundle, monkeypatch):
    mod = load_script("download_gpn_star_scores")
    payload = json.loads((bundle / "manifest.json").read_text())
    files = [{**entry, "path": "data/" + entry["path"]} for entry in payload["files"]]
    release = {"repository": {"id": mod.REPOSITORY}, "parquet": {"files": files}}
    dest = tmp_path / "incomplete"
    interrupted = False
    def fake_curl(command, **kwargs):
        nonlocal interrupted
        target = Path(command[command.index("--output") + 1])
        source = bundle / command[-1].split("/data/", 1)[1]
        if not interrupted:
            target.write_bytes(source.read_bytes()[:100])
            interrupted = True
            raise subprocess.CalledProcessError(18, command)
        assert "--continue-at" in command
        if target.exists():
            assert target.stat().st_size == 100
        target.write_bytes(source.read_bytes())
    monkeypatch.setattr(mod.subprocess, "run", fake_curl)
    with pytest.raises(subprocess.CalledProcessError):
        mod.install(dest, release, "a" * 40, "https://example.invalid")
    assert not (dest / "manifest.json").exists()
    mod.install(dest, release, "a" * 40, "https://example.invalid")
    assert gpn_star.validate_database(dest)["revision"] == "a" * 40


def test_duplicate_alleles_fail_without_replacing_review(tmp_path, bundle):
    path = gpn_star.shard_path(bundle, "1")
    data = pl.read_parquet(path)
    pl.concat([data, data.head(1)]).write_parquet(path)
    marker = bundle / "manifest.json"
    payload = json.loads(marker.read_text())
    entry = next(f for f in payload["files"] if f["chrom"] == "1")
    entry["size"] = path.stat().st_size
    entry["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    marker.write_text(json.dumps(payload))
    review = tmp_path / "review.tsv"
    write_tsv(review, [row()])
    before = review.read_bytes()
    stats = gpn_star.annotate_review_tsv(review, bundle)
    assert stats["status"] == "failed"
    assert "duplicate" in stats["message"]
    assert review.read_bytes() == before


def test_v_upgrade_queries_only_v_and_preserves_msa_m_and_raw(tmp_path, bundle, msa_db, monkeypatch):
    raw = tmp_path / 'raw.tsv'
    write_tsv(raw, [row(), row(pos=11), row(ref='T'), row(pos=30, alt='AT')])
    original = raw.read_bytes()
    monkeypatch.setattr(snv_review, '_candidate_bed_path', lambda: tmp_path / 'no.bed')
    args = dict(test_type='WGS', output_dir=tmp_path / 'post', gpn_star_dir=bundle, gpn_msa_db=msa_db[0])
    review = snv_review.ensure_review_tsv(raw, **args, research_only=True)
    before = read_tsv(review)
    marker = review.with_suffix(review.suffix + '.source.json')
    old = json.loads(marker.read_text())
    # Simulate a pre-V manifest, then install V without changing M.
    for key in list(old):
        if key in {'gpn_star_v', 'gpn_star_v_annotation', 'gpn_star_v_values'}:
            old.pop(key)
    marker.write_text(json.dumps(old))
    make_bundle(bundle, gpn_star.V_MODEL)
    calls = []
    query = gpn_star.query_scores
    def spy(root, keys, model):
        calls.append(model)
        assert model == gpn_star.V_MODEL
        return query(root, keys, model)
    monkeypatch.setattr(gpn_star, 'query_scores', spy)
    def forbidden(*args):
        raise AssertionError('must not rebuild from raw')
    monkeypatch.setattr(snv_review, '_load_bed', forbidden)
    snv_review.ensure_review_tsv(raw, **args)
    snv_review.ensure_review_tsv(raw, **args)
    assert calls == [gpn_star.V_MODEL]
    assert len(msa_db[1]) == 1
    actual = read_tsv(review)
    assert [r['GPN_STAR_V_SCORE'] for r in actual] == ['-7', '0', '', '']
    assert [{k: v for k, v in r.items() if k != 'GPN_STAR_V_SCORE'} for r in actual] == before
    assert raw.read_bytes() == original
    assert _row_to_variant(actual[0])['GPN_STAR_V_score'] == -7
    assert _row_to_variant(actual[1])['GPN_STAR_V_score'] == 0
    assert classify_tier({'GPN_STAR_V_SCORE': '-99', 'GPN_STAR_SCORE': '-99'}) == '2'
    saved = json.loads(marker.read_text())
    assert saved['gpn_star_values'] == old['gpn_star_values']
    assert saved['gpn_star_v_values']['model'] == gpn_star.V_MODEL
    # V failure preserves both the old V value and its successful provenance.
    manifest = gpn_star.manifest_path(bundle, gpn_star.V_MODEL)
    data = json.loads(manifest.read_text()); data['model'] = 'wrong'
    manifest.write_text(json.dumps(data))
    snv_review.ensure_review_tsv(raw, **args)
    failed = json.loads(marker.read_text())
    assert failed['gpn_star_v_annotation']['status'] == 'failed'
    assert failed['gpn_star_v_values'] == saved['gpn_star_v_values']
    assert read_tsv(review) == actual


def test_v_can_succeed_when_m_is_unavailable_and_strict_requires_both(tmp_path, dual_bundle, monkeypatch):
    raw = tmp_path / 'raw.tsv'
    write_tsv(raw, [row()])
    (dual_bundle / 'manifest.json').unlink()
    monkeypatch.setattr(snv_review, '_candidate_bed_path', lambda: tmp_path / 'no.bed')
    args = dict(test_type='WGS', output_dir=tmp_path / 'post', gpn_star_dir=dual_bundle, research_only=True)
    review = snv_review.ensure_review_tsv(raw, **args)
    assert read_tsv(review)[0]['GPN_STAR_V_SCORE'] == '-7'
    assert 'GPN_STAR_SCORE' not in read_tsv(review)[0]
    saved = json.loads(review.with_suffix(review.suffix + '.source.json').read_text())
    assert saved['gpn_star_annotation']['status'] == 'skipped_missing_db'
    assert saved['gpn_star_v_annotation']['status'] == 'complete'
    with pytest.raises(FileNotFoundError):
        snv_review.ensure_review_tsv(raw, **args, require_gpn_star=True)


def test_download_v_addition_keeps_active_m_and_partial_files_separate(tmp_path, dual_bundle, monkeypatch):
    mod = load_script('download_gpn_star_scores')
    files = []
    for model in gpn_star.MODELS:
        payload = json.loads(gpn_star.manifest_path(dual_bundle, model).read_text())
        files.extend({**f, 'path': 'data/' + f['path']} for f in payload['files'])
    release = {'repository': {'id': mod.REPOSITORY}, 'parquet': {'files': files}}
    dest = tmp_path / 'install'
    interrupted = False
    transfers = []
    def curl(command, **kwargs):
        nonlocal interrupted
        target = Path(command[command.index('--output') + 1])
        source = dual_bundle / command[-1].split('/data/', 1)[1]
        transfers.append(source)
        if gpn_star.V_MODEL in source.parts and not interrupted:
            target.write_bytes(source.read_bytes()[:100])
            interrupted = True
            raise subprocess.CalledProcessError(18, command)
        target.write_bytes(source.read_bytes())
    monkeypatch.setattr(mod.subprocess, 'run', curl)
    mod.install(dest, release, 'a' * 40, 'https://example.invalid')
    old_marker = (dest / 'manifest.json').read_bytes()
    old_shard = gpn_star.shard_path(dest, '1').read_bytes()
    with pytest.raises(subprocess.CalledProcessError):
        mod.install(dest, release, 'a' * 40, 'https://example.invalid', gpn_star.V_MODEL)
    assert not gpn_star.manifest_path(dest, gpn_star.V_MODEL).exists()
    assert (dest / 'manifest.json').read_bytes() == old_marker
    assert gpn_star.shard_path(dest, '1').read_bytes() == old_shard
    mod.install(dest, release, 'a' * 40, 'https://example.invalid', gpn_star.V_MODEL)
    assert gpn_star.validate_database(dest, model=gpn_star.V_MODEL)['model'] == gpn_star.V_MODEL
    assert len([p for p in transfers if gpn_star.MODEL in p.parts]) == 24
    assert (dest / 'manifest.json').read_bytes() == old_marker
    assert gpn_star.query_scores(dest, {('1', 10, 'A', 'G')}, gpn_star.V_MODEL)[('1', 10, 'A', 'G')] == '-7'
