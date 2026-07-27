import csv
import gzip
import importlib.util
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "annotate_acmg_genebe.py"
SPEC = importlib.util.spec_from_file_location("annotate_acmg_genebe", SCRIPT)
genebe = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(genebe)


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "CHROM", "POS", "REF", "ALT", "CLINVAR_SIG",
        "GNOMAD_G_AF", "DP", "ACMG_CLASS",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_genebe_db(path: Path) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(
            "#chr\tpos\tref\talt\tacmg_classification\t"
            "acmg_score\tacmg_criteria\n"
        )
        handle.write("chr1\t100\tA\tG\tLikely_benign\t-4\tBP4_Strong\n")


def test_db_classification_matches_official_seven_column_spelling():
    assert genebe.db_classification("Likely benign") == "Likely_benign"
    assert genebe.db_classification("Uncertain significance") == "VUS"
    assert genebe.db_classification(score=7) == "Likely_pathogenic"
    assert genebe._score_for_export("-4.0") == "-4"


def test_api_cache_success_and_negative_ttl(tmp_path):
    cache = tmp_path / "api.sqlite"
    hit_key = "chr1:100:A:G"
    miss_key = "chr1:200:C:T"
    now = datetime(2026, 7, 27, tzinfo=timezone.utc)

    genebe.cache_api_outcomes(
        cache,
        {hit_key: ("7", "PM2,PP3", "Likely pathogenic")},
        {miss_key},
        negative_ttl_days=30,
        now=now,
    )

    hits, negative = genebe.api_cache_lookup(
        cache,
        {hit_key, miss_key},
        now_epoch=int(now.timestamp()),
    )
    assert hits == {
        hit_key: ("7", "PM2,PP3", "Likely_pathogenic"),
    }
    assert negative == {miss_key}

    _, expired = genebe.api_cache_lookup(
        cache,
        {miss_key},
        now_epoch=int(now.timestamp()) + 31 * 86400,
    )
    assert expired == set()


def test_pending_tsv_is_exact_schema_sorted_and_deduplicated(tmp_path):
    db = tmp_path / "genebe_hg38.tsv.gz"
    _write_genebe_db(db)
    pending_dir = tmp_path / "pending"
    hits = {
        "chr2:20:G:A": ("-4.0", "BP4_Strong", "Likely benign"),
        "chr1:10:A:T": ("7", "PM2,PP3", "Likely pathogenic"),
    }

    output = genebe.write_pending_tsv(
        pending_dir,
        hits,
        source_db=db,
        queried_count=3,
        no_result_count=1,
        failed_count=0,
        sif=tmp_path / "genebe.sif",
    )

    assert output is not None
    assert output.read_text(encoding="utf-8").splitlines() == [
        "#chr\tpos\tref\talt\tacmg_classification\tacmg_score\tacmg_criteria",
        "chr1\t10\tA\tT\tLikely_pathogenic\t7\tPM2,PP3",
        "chr2\t20\tG\tA\tLikely_benign\t-4\tBP4_Strong",
    ]
    sidecar = json.loads(output.with_suffix(".json").read_text(encoding="utf-8"))
    assert sidecar["rows"] == 2
    assert sidecar["columns"] == list(genebe.DB_EXPORT_FIELDS)


def test_success_cache_can_be_consolidated_for_later_import(tmp_path):
    cache = tmp_path / "api.sqlite"
    genebe.cache_api_outcomes(
        cache,
        {
            "chr2:20:G:A": ("-4", "BP4_Strong", "Likely_benign"),
            "chr1:10:A:T": ("7", "PM2,PP3", "Likely_pathogenic"),
        },
        set(),
        negative_ttl_days=30,
    )
    output = tmp_path / "all_api_rows.tsv"
    exporter = SCRIPT.with_name("export_genebe_api_cache.py")

    subprocess.run(
        [sys.executable, str(exporter), "--cache", str(cache), "--out", str(output)],
        check=True,
        text=True,
        capture_output=True,
    )

    assert output.read_text(encoding="utf-8").splitlines() == [
        "#chr\tpos\tref\talt\tacmg_classification\tacmg_score\tacmg_criteria",
        "chr1\t10\tA\tT\tLikely_pathogenic\t7\tPM2,PP3",
        "chr2\t20\tG\tA\tLikely_benign\t-4\tBP4_Strong",
    ]


def test_api_candidates_share_review_rescue_and_wes_depth_rules(
    tmp_path,
    monkeypatch,
):
    from app.services import snv_review

    tsv = tmp_path / "working.tsv"
    _write_tsv(tsv, [
        {
            "CHROM": "chr1", "POS": "10", "REF": "A", "ALT": "G",
            "CLINVAR_SIG": "Pathogenic", "GNOMAD_G_AF": "0.5", "DP": "30",
            "ACMG_CLASS": "VUS",
        },
        {
            "CHROM": "chr1", "POS": "20", "REF": "C", "ALT": "T",
            "CLINVAR_SIG": "", "GNOMAD_G_AF": "0.001", "DP": "10",
            "ACMG_CLASS": "VUS",
        },
    ])
    unresolved = {"chr1:10:A:G", "chr1:20:C:T"}
    monkeypatch.setattr(snv_review, "load_candidate_bed", lambda: None)

    assert genebe.collect_api_candidates(
        tsv, unresolved, test_type="WES",
    ) == {"chr1:10:A:G"}
    assert genebe.collect_api_candidates(
        tsv, unresolved, test_type="WGS",
    ) == unresolved


def test_api_batch_keeps_credentials_out_of_host_command_arguments(
    tmp_path,
    monkeypatch,
):
    sif = tmp_path / "genebe.sif"
    sif.touch()
    captured = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["env"] = kwargs["env"]
        output = Path(command[-1])
        output.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t200\t.\tC\tT\t.\t.\tacmg_score=7;acmg_criteria=PM2,PP3\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(genebe.subprocess, "run", fake_run)
    hits = genebe.run_api_batch(
        {"chr1:200:C:T"},
        sif=sif,
        username="secret-user",
        api_key="secret-key",
        timeout_seconds=5,
        retries=1,
    )

    assert "secret-user" not in " ".join(captured["command"])
    assert "secret-key" not in " ".join(captured["command"])
    assert captured["env"]["APPTAINERENV_GENEBE_API_KEY"] == "secret-key"
    assert hits == {
        "chr1:200:C:T": ("7", "PM2,PP3", "Likely_pathogenic"),
    }


def test_main_uses_db_then_review_filtered_api_and_reuses_cache(
    tmp_path,
    monkeypatch,
):
    from app.services import snv_review

    db = tmp_path / "genebe_hg38.tsv.gz"
    _write_genebe_db(db)
    tsv = tmp_path / "working.tsv"
    _write_tsv(tsv, [
        {
            "CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "G",
            "CLINVAR_SIG": "", "GNOMAD_G_AF": "0.2", "DP": "30",
            "ACMG_CLASS": "VUS",
        },
        {
            "CHROM": "chr1", "POS": "200", "REF": "C", "ALT": "T",
            "CLINVAR_SIG": "", "GNOMAD_G_AF": "0.001", "DP": "30",
            "ACMG_CLASS": "VUS",
        },
        {
            "CHROM": "chr1", "POS": "300", "REF": "G", "ALT": "A",
            "CLINVAR_SIG": "", "GNOMAD_G_AF": "0.2", "DP": "30",
            "ACMG_CLASS": "VUS",
        },
    ])
    sif = tmp_path / "genebe.sif"
    sif.touch()
    cache = tmp_path / "api.sqlite"
    pending = tmp_path / "pending"
    calls: list[set[str]] = []

    def fake_api(candidates, **_kwargs):
        calls.append(set(candidates))
        return (
            {"chr1:200:C:T": ("7", "PM2,PP3", "Likely_pathogenic")},
            set(),
            0,
        )

    monkeypatch.setattr(genebe, "run_live_api", fake_api)
    monkeypatch.setattr(snv_review, "load_candidate_bed", lambda: None)
    monkeypatch.setenv("GENEBE_USER", "test-user")
    monkeypatch.setenv("GENEBE_API_KEY", "test-key")
    argv = [
        str(SCRIPT), "--tsv", str(tsv), "--genebe-db", str(db),
        "--no-sqlite", "--test-type", "WES", "--api-sif", str(sif),
        "--api-cache", str(cache), "--api-pending-dir", str(pending),
    ]
    monkeypatch.setattr(sys, "argv", argv)

    assert genebe.main() == 0
    assert calls == [{"chr1:200:C:T"}]
    with tsv.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["GENEBE_ACMG_SCORE"] == "-4"
    assert rows[1]["GENEBE_ACMG_CLASS"] == "Likely_pathogenic"
    assert rows[2]["GENEBE_ACMG_CLASS"] == ""
    assert len(list(pending.glob("*.tsv"))) == 1

    def should_not_call(*_args, **_kwargs):
        raise AssertionError("live API should be bypassed by its cache")

    monkeypatch.setattr(genebe, "run_live_api", should_not_call)
    monkeypatch.setattr(sys, "argv", argv)
    assert genebe.main() == 0
    assert len(list(pending.glob("*.tsv"))) == 1
