import json

from app.services import manual_acmg
from app.services import report_store
from app.services import sample_loader


def test_catalog_contains_all_28_criteria_and_keeps_pp5_bp6_usable():
    payload = manual_acmg.catalog()
    by_code = {item["code"]: item for item in payload["criteria"]}

    assert len(by_code) == 28
    assert set(by_code) == set(manual_acmg.CRITERIA_ORDER)
    assert by_code["PP5"]["deprecated_warning"]
    assert by_code["BP6"]["deprecated_warning"]
    assert by_code["PS2"]["scope"] == "case"
    assert by_code["PP4"]["scope"] == "case"
    assert by_code["PVS1"]["scope"] == "global"
    assert "disabled" not in by_code["PP5"]
    assert "disabled" not in by_code["BP6"]
    assert payload["vus_subclasses"] == [
        {"value": "VUS-low", "min_points": 0, "max_points": 1},
        {"value": "VUS-mid", "min_points": 2, "max_points": 3},
        {"value": "VUS-high", "min_points": 4, "max_points": 5},
    ]
    assert payload["vus_subclass_reference"]["url"].startswith("https://doi.org/")


def test_points_classification_strengths_and_ba1():
    result = manual_acmg.calculate({
        "PVS1": {"enabled": True, "strength": "very_strong"},
        "PM2": {"enabled": True, "strength": "supporting"},
        "PP3": {"enabled": True, "strength": "supporting"},
    })
    assert result["score"] == 10
    assert result["classification"] == "Pathogenic"
    assert result["criteria_text"] == "PVS1,PM2_Supporting,PP3"

    benign = manual_acmg.calculate({
        "BA1": {"enabled": True, "strength": "stand_alone"},
    })
    assert benign["score"] == -8
    assert benign["classification"] == "Benign"
    assert manual_acmg.acmg_to_variant_score(benign["score"]) == 10
    parsed, unknown = manual_acmg.parse_criteria_text(
        "PVS1_Strong PM2_Supporting|PP3"
    )
    assert set(parsed) == {"PVS1", "PM2", "PP3"}
    assert unknown == []


def test_vus_subclasses_follow_published_point_bins_without_changing_formal_class():
    criteria_by_points = {
        0: {},
        1: {"PP2": {"enabled": True, "strength": "supporting"}},
        2: {"PM1": {"enabled": True, "strength": "moderate"}},
        3: {
            code: {"enabled": True, "strength": "supporting"}
            for code in ("PP2", "PP3", "PP5")
        },
        4: {"PS1": {"enabled": True, "strength": "strong"}},
        5: {
            "PS1": {"enabled": True, "strength": "strong"},
            "PP2": {"enabled": True, "strength": "supporting"},
        },
    }
    expected = {
        0: "VUS-low",
        1: "VUS-low",
        2: "VUS-mid",
        3: "VUS-mid",
        4: "VUS-high",
        5: "VUS-high",
    }
    for points, subclass in expected.items():
        result = manual_acmg.calculate(criteria_by_points[points])
        assert result["score"] == points
        assert result["classification"] == "Uncertain significance"
        assert result["vus_subclass"] == subclass
    assert manual_acmg.vus_subclass("Likely pathogenic", 5) == ""
    assert manual_acmg.vus_subclass("VUS", 2.5) == ""


def test_manual_revision_is_append_only_and_current_moves(tmp_path):
    db = tmp_path / "manual.sqlite"
    first = manual_acmg.save_assertion(
        "GRCh38",
        "1-100-a-g",
        {"PP3": {"enabled": True, "strength": "supporting"}},
        reviewer_user_id=7,
        reviewer_username="alice",
        source_sample_id="CASE-A",
        path=db,
    )
    second = manual_acmg.save_assertion(
        "hg38",
        "chr1-100-A-G",
        {"PVS1": {"enabled": True, "strength": "very_strong"}},
        reviewer_user_id=8,
        reviewer_username="bob",
        source_sample_id="CASE-B",
        path=db,
    )

    assert first["classification"] == "Uncertain significance"
    assert first["vus_subclass"] == "VUS-low"
    assert second["revision_id"] > first["revision_id"]
    current = manual_acmg.current_assertion("hg38", "chr1-100-A-G", path=db)
    assert current["revision_id"] == second["revision_id"]
    assert current["reviewer_username"] == "bob"
    assert current["source_sample_id"] == "CASE-B"
    assert current["reusable_criteria_text"] == "PVS1"


def test_case_scoped_criteria_are_saved_but_not_globally_applied(tmp_path):
    db = tmp_path / "manual.sqlite"
    saved = manual_acmg.save_assertion(
        "hg38",
        "chr1-100-A-G",
        {
            "PVS1": {"enabled": True, "strength": "very_strong"},
            "PS2": {"enabled": True, "strength": "strong"},
            "PP4": {"enabled": True, "strength": "supporting"},
        },
        reviewer_user_id=7,
        reviewer_username="alice",
        source_sample_id="CASE-A",
        path=db,
    )

    assert saved["criteria_text"] == "PVS1,PS2,PP4"
    assert saved["score"] == 13
    assert saved["classification"] == "Pathogenic"
    assert saved["reusable_criteria_text"] == "PVS1"
    assert saved["reusable_score"] == 8
    assert saved["reusable_classification"] == "Likely pathogenic"


def test_observed_registry_only_keeps_active_causative_and_other(tmp_path):
    db = tmp_path / "manual.sqlite"
    manual_acmg.sync_observations(
        "hg38",
        "CASE-A",
        {
            "chr1-100-A-G": "1",
            "chr2-200-C-T": "2",
            "chr3-300-G-A": "C",
            "chr4-400-T-C": "0",
        },
        reviewer_user_id=7,
        reviewer_username="alice",
        path=db,
    )
    manual_acmg.sync_observations(
        "hg38",
        "CASE-B",
        {"chr1-100-A-G": "2"},
        reviewer_user_id=8,
        reviewer_username="bob",
        path=db,
    )

    counts = manual_acmg.bulk_observed_counts(
        "hg38",
        ["chr1-100-A-G", "chr2-200-C-T", "chr3-300-G-A"],
        exclude_sample_id="CASE-A",
        path=db,
    )
    assert counts == {"chr1-100-A-G": 1}
    cases = manual_acmg.observed_cases(
        "hg38", "chr1-100-A-G", exclude_sample_id="CASE-A", path=db
    )
    assert [(item["sample_id"], item["status_label"]) for item in cases] == [
        ("CASE-B", "Other")
    ]

    # Cancelling the status physically removes it from current observations.
    manual_acmg.sync_observations(
        "hg38",
        "CASE-B",
        {"chr1-100-A-G": "0"},
        reviewer_user_id=8,
        reviewer_username="bob",
        path=db,
    )
    assert manual_acmg.observed_cases(
        "hg38", "chr1-100-A-G", exclude_sample_id="CASE-A", path=db
    ) == []


def test_effective_overlay_prefers_sample_snapshot_and_reorders(monkeypatch):
    variants = {
        "chr1-100-A-G": {
            "tier": "2",
            "predicted_suspect_non_acmg": False,
            "genebe_acmg_class": "Likely benign",
            "genebe_acmg_score": -2,
            "genebe_acmg_criteria": "BP4",
            "ACMG_classification": "Uncertain significance",
            "ACMG_score": 1,
            "ACMG_criteria": "PP3",
            "pheno_score": 20,
        },
        "chr2-200-C-T": {
            "tier": "1A",
            "predicted_suspect_non_acmg": False,
            "ACMG_classification": "Uncertain significance",
            "ACMG_score": 0,
            "ACMG_criteria": "",
            "pheno_score": 0,
        },
    }
    categories = {"1A": ["chr2-200-C-T"], "1B": [], "1C": [], "2": ["chr1-100-A-G"]}
    monkeypatch.setattr(manual_acmg, "bulk_current", lambda *_a, **_k: {})
    monkeypatch.setattr(
        manual_acmg,
        "bulk_observed_counts",
        lambda *_a, **_k: {"chr1-100-A-G": 3},
    )
    meta = {
        "edits": {
            "chr1-100-A-G": {
                "manual_acmg": {
                    "classification": "Pathogenic",
                    "score": 10,
                    "criteria_text": "PVS1,PM2_Supporting,PP3",
                    "criteria": {},
                }
            }
        }
    }

    sample_loader._apply_effective_acmg(
        variants,
        categories,
        sample_id="CASE-A",
        genome_build="hg38",
        meta=meta,
    )

    assert variants["chr1-100-A-G"]["effective_acmg_source"] == "manual"
    assert variants["chr1-100-A-G"]["effective_acmg_scope"] == "sample"
    assert variants["chr1-100-A-G"]["effective_acmg_vus_subclass"] == ""
    assert variants["chr1-100-A-G"]["geno_score"] == 100
    assert variants["chr1-100-A-G"]["total_score"] == 120
    assert variants["chr1-100-A-G"]["observed_count"] == 3
    assert variants["chr1-100-A-G"]["tier"] == "1C"
    assert categories["1C"] == ["chr1-100-A-G"]
    assert variants["chr2-200-C-T"]["tier"] == "1A"
    assert variants["chr2-200-C-T"]["effective_acmg_vus_subclass"] == "VUS-low"


def test_report_save_preserves_newer_manual_snapshot_and_audits_status(
    tmp_path, monkeypatch
):
    meta_path = tmp_path / "sample_metadata.json"
    meta_path.write_text(json.dumps({
        "sample_id": "CASE-A",
        "genome_build": "hg38",
        "status": {},
        "edits": {
            "chr1-100-A-G": {
                "manual_acmg": {
                    "revision_id": 2,
                    "classification": "Pathogenic",
                    "score": 10,
                    "criteria_text": "PVS1,PM2_Supporting,PP3",
                },
            },
        },
    }), encoding="utf-8")
    monkeypatch.setattr(report_store, "_meta_path", lambda _sample: meta_path)
    synced = {}
    monkeypatch.setattr(
        manual_acmg,
        "sync_observations",
        lambda *args, **kwargs: synced.update({
            "statuses": args[2],
            "audit": kwargs["status_audit"],
        }),
    )

    report_store.save(
        "CASE-A",
        {
            "status": {"chr1-100-A-G": "1"},
            "edits": {
                "chr1-100-A-G": {
                    "manual_acmg": {
                        "revision_id": 1,
                        "classification": "Likely pathogenic",
                        "score": 6,
                        "criteria_text": "PVS1",
                    },
                },
            },
        },
        user={"id": 9, "username": "reviewer"},
    )

    stored = json.loads(meta_path.read_text(encoding="utf-8"))
    edit = stored["edits"]["chr1-100-A-G"]
    assert edit["manual_acmg"]["revision_id"] == 2
    assert edit["ACMG_classification"] == "Pathogenic"
    assert synced["statuses"] == {"chr1-100-A-G": "1"}
    assert synced["audit"]["chr1-100-A-G"]["reviewer_user_id"] == 9
    assert synced["audit"]["chr1-100-A-G"]["reviewer_username"] == "reviewer"
