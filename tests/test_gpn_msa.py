import csv
import json
from pathlib import Path

import pytest

from app.services import gpn_msa, snv_review
from app.adapters.snv_tsv import _row_to_variant, classify_tier


FIELDS = [
    "CHROM", "POS", "REF", "ALT", "GENE", "CLINVAR_SIG",
    "GNOMAD_G_AF", "DP",
]


def _write(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def test_gpn_msa_is_joined_only_to_retained_review_snvs(tmp_path, monkeypatch):
    raw = tmp_path / "raw.tsv"
    _write(raw, [
        {
            "CHROM": "chr1", "POS": "10", "REF": "A", "ALT": "G",
            "GENE": "KEEP", "CLINVAR_SIG": "", "GNOMAD_G_AF": "0.001", "DP": "30",
        },
        {
            "CHROM": "chr1", "POS": "20", "REF": "C", "ALT": "T",
            "GENE": "DROP", "CLINVAR_SIG": "", "GNOMAD_G_AF": "0.2", "DP": "30",
        },
        {
            "CHROM": "chr1", "POS": "30", "REF": "A", "ALT": "AT",
            "GENE": "INDEL", "CLINVAR_SIG": "Pathogenic", "GNOMAD_G_AF": "", "DP": "30",
        },
    ])
    db = tmp_path / "scores.tsv.bgz"
    db.touch()
    Path(f"{db}.tbi").touch()
    monkeypatch.setattr(snv_review, "_candidate_bed_path", lambda: tmp_path / "missing.bed")
    monkeypatch.setattr(gpn_msa, "validate_database", lambda _db: (db, "tabix"))
    captured = {}

    def fake_query(_db, keys, *, tabix_bin):
        captured["keys"] = keys
        captured["tabix"] = tabix_bin
        return {("1", 10, "A", "G"): "-8.25"}

    monkeypatch.setattr(gpn_msa, "query_scores", fake_query)
    review = snv_review.ensure_review_tsv(
        raw,
        test_type="WGS",
        output_dir=tmp_path / "08_postprocessing",
        gpn_msa_db=db,
        require_gpn_msa=True,
    )

    with review.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["GENE"] for row in rows] == ["KEEP", "INDEL"]
    assert rows[0]["GPN_MSA_SCORE"] == "-8.25"
    assert rows[1]["GPN_MSA_SCORE"] == ""
    assert _row_to_variant(rows[0])["GPN_MSA_score"] == -8.25
    assert captured == {
        "keys": {("1", 10, "A", "G")},
        "tabix": "tabix",
    }
    assert "GPN_MSA_SCORE" not in raw.read_text(encoding="utf-8").splitlines()[0]
    manifest = json.loads(
        review.with_suffix(review.suffix + ".source.json").read_text(encoding="utf-8")
    )
    assert manifest["gpn_msa"]["database"]["path"] == str(db)


def test_gpn_msa_score_is_display_only_and_does_not_trigger_tier_1c():
    assert classify_tier({"GPN_MSA_SCORE": "-12"}) == "2"


def test_required_gpn_msa_rejects_missing_deployment_data(tmp_path):
    review = tmp_path / "review.tsv"
    _write(review, [])
    with pytest.raises(FileNotFoundError, match="required GPN-MSA input"):
        gpn_msa.annotate_review_tsv(
            review,
            tmp_path / "missing.tsv.bgz",
            required=True,
        )


def test_query_scores_matches_exact_alt_and_normalises_chr(monkeypatch, tmp_path):
    class Result:
        returncode = 0
        stderr = ""
        stdout = (
            "1\t10\tA\tC\t-1.2\n"
            "1\t10\tA\tG\t-8.25\n"
            "1\t10\tA\tT\tnan-value\n"
        )

    captured = {}

    def fake_run(command, **_kwargs):
        captured["command"] = command
        return Result()

    monkeypatch.setattr(gpn_msa.subprocess, "run", fake_run)
    db = tmp_path / "scores.tsv.bgz"
    found = gpn_msa.query_scores(
        db,
        {("1", 10, "A", "G")},
        tabix_bin="/usr/bin/tabix",
    )
    assert found == {("1", 10, "A", "G"): "-8.25"}
    assert captured["command"] == [
        "/usr/bin/tabix", str(db), "1:10-10",
    ]
