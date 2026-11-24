"""院內常見 marking and ranking penalty (WES recurrence + in-house WGS AF)."""
from pathlib import Path

import pytest

from app.services import inhouse_common, manual_acmg, sample_loader
from app.services.acmg_engine import runner
from app.services.acmg_engine.recurrence import WesRecurrence


@pytest.fixture
def store(tmp_path, monkeypatch):
    db = WesRecurrence(tmp_path / "rec.sqlite")
    common = "chr1-100-A-G"   # 12 of 120 other patients → 10%
    medium = "chr1-200-C-T"   # 7 of 120 → 5.8%
    for i in range(121):
        keys = ["chr2-5-G-A"]
        if i < 12:
            keys.append(common)
        if i < 7:
            keys.append(medium)
        db.upsert_sample(f"23WE{i:04d}", f"mrn:{i}", {}, keys)
    db.upsert_sample("24WE9999", "mrn:0", {}, [common])          # same patient as #0: counted once
    monkeypatch.setattr(runner, "shared_recurrence", lambda: db)
    return db


def test_lookup_many_counts_distinct_patients_and_excludes_case(store):
    counts = store.lookup_many({"chr1-100-A-G", "chr1-200-C-T", "chr9-1-A-T"}, exclude_patient="mrn:0")
    assert counts == {"chr1-100-A-G": 11, "chr1-200-C-T": 6}
    assert store.lookup_many({"chr1-100-A-G"}) == {"chr1-100-A-G": 12}


def test_annotate_wes_levels_and_hiding(store):
    variants = {"chr1-100-A-G": {}, "chr1-200-C-T": {}, "chr2-5-G-A": {}, "chr3-1-A-T": {}}
    inhouse_common.annotate(variants, sample_id="25WE0001", test_type="WES", mrn="mrn:new")
    assert variants["chr1-100-A-G"]["wes_recurrence"] == {"count": 12, "patients": 121, "fraction": 0.0992}
    assert variants["chr1-100-A-G"]["inhouse_common"]["level"] == 1       # 9.9% < 10%
    assert variants["chr2-5-G-A"]["inhouse_common"]["level"] == 2
    assert variants["chr2-5-G-A"]["inhouse_common"]["hide_wes"] is True
    assert variants["chr1-200-C-T"]["inhouse_common"]["target"] == 4
    assert "inhouse_common" not in variants["chr3-1-A-T"]
    assert variants["chr3-1-A-T"]["wes_recurrence"]["count"] == 0
    # WGS cases ignore the WES table but use the in-house WGS AF (AN ≥ 2000).
    wgs = {"chr1-100-A-G": {"inhouse_af": 0.03, "inhouse_an": 2800},
           "chr1-200-C-T": {"inhouse_af": 0.06, "inhouse_an": 1500},
           "chr2-5-G-A": {"inhouse_af": 0.051, "inhouse_an": 2800}}
    inhouse_common.annotate(wgs, sample_id="25G00001", test_type="WGS")
    assert "wes_recurrence" not in wgs["chr1-100-A-G"]
    assert wgs["chr1-100-A-G"]["inhouse_common"]["level"] == 1
    assert "inhouse_common" not in wgs["chr1-200-C-T"]                     # AN too small
    assert wgs["chr2-5-G-A"]["inhouse_common"]["target"] == 8
    assert wgs["chr2-5-G-A"]["inhouse_common"]["hide_wes"] is False


def test_small_wes_table_is_not_used(tmp_path, monkeypatch):
    db = WesRecurrence(tmp_path / "rec.sqlite")
    for i in range(20):
        db.upsert_sample(f"23WE{i:04d}", f"mrn:{i}", {}, ["chr1-1-A-G"])
    monkeypatch.setattr(runner, "shared_recurrence", lambda: db)
    variants = {"chr1-1-A-G": {}}
    inhouse_common.annotate(variants, sample_id="25WE0001", test_type="WES")
    assert variants["chr1-1-A-G"] == {}


def test_penalty_counts_existing_bs1_ba1_and_exemptions():
    base = {"inhouse_common": {"level": 2, "target": 8, "basis": "x"}}
    v = {**base, "inhouse_common": dict(base["inhouse_common"])}
    assert inhouse_common.apply_penalty(v, scope="in_house", score=6, criteria_text="PM2_Supporting,PP3") == -2
    assert v["inhouse_common"]["penalty"] == 8
    v = {"inhouse_common": {"level": 2, "target": 8}}
    assert inhouse_common.apply_penalty(v, scope="in_house", score=-2, criteria_text="BS1,PP3_Moderate") == -6
    assert v["inhouse_common"]["penalty"] == 4
    v = {"inhouse_common": {"level": 1, "target": 4}}
    assert inhouse_common.apply_penalty(v, scope="genebe", score=-8, criteria_text="BA1") == -8
    assert v["inhouse_common"]["penalty"] == 0
    v = {"inhouse_common": {"level": 1, "target": 4}}
    assert inhouse_common.apply_penalty(v, scope="sample", score=7, criteria_text="") == 7
    assert v["inhouse_common"]["exempt"]
    v = {"inhouse_common": {"level": 1, "target": 4}, "CLNSIG": "Pathogenic", "clinvar_stars": 2}
    assert inhouse_common.apply_penalty(v, scope="in_house", score=7, criteria_text="") == 7
    v = {"inhouse_common": {"level": 1, "target": 4}, "CLNSIG": "Pathogenic", "clinvar_stars": 1}
    assert inhouse_common.apply_penalty(v, scope="in_house", score=7, criteria_text="") == 3


def test_effective_acmg_ranks_with_penalty_but_keeps_class_and_tier(store, monkeypatch):
    monkeypatch.setattr(manual_acmg, "bulk_current", lambda *a, **k: {})
    monkeypatch.setattr(manual_acmg, "bulk_observed_counts", lambda *a, **k: {})
    variants = {
        "chr2-5-G-A": {"ACMG_classification": "Likely pathogenic", "ACMG_score": 6, "ACMG_criteria": "PM2,PP3",
                       "tier": "1C", "pheno_score": 10},
        "chr3-1-A-T": {"ACMG_classification": "Likely pathogenic", "ACMG_score": 6, "ACMG_criteria": "PM2,PP3",
                       "tier": "1C", "pheno_score": 10},
    }
    categories = {"1A": [], "1B": [], "1C": list(variants), "2": []}
    sample_loader._apply_effective_acmg(variants, categories, sample_id="25WE0001", genome_build="hg38",
                                        meta={"test_type": "WES", "mrn": "new"})
    common, other = variants["chr2-5-G-A"], variants["chr3-1-A-T"]
    assert common["effective_acmg_class"] == "Likely pathogenic" and common["tier"] == "1C"
    assert common["inhouse_common"]["penalty"] == 8
    assert common["total_score"] == manual_acmg.acmg_to_variant_score(-2) + 10
    assert other["total_score"] == manual_acmg.acmg_to_variant_score(6) + 10
    assert categories["1C"] == ["chr3-1-A-T", "chr2-5-G-A"]
