"""Phase 0: cleaned HPO → gene table (plan §4.2.3 R1–R7, §10.7)."""
import csv
import json

import pytest

from app.services import hpo_gene_table as H
from app.services import pheno_table_diff, phenotype_scorer

OBO = """format-version: 1.2
data-version: hp/releases/2026-09-01

[Term]
id: HP:0000001
name: All

[Term]
id: HP:0000118
name: Phenotypic abnormality
is_a: HP:0000001 ! All

[Term]
id: HP:0000707
name: Abnormality of the nervous system
is_a: HP:0000118

[Term]
id: HP:0001250
name: Seizure
alt_id: HP:0009999
is_a: HP:0000707

[Term]
id: HP:0001249
name: Intellectual disability
is_a: HP:0000707

[Term]
id: HP:0000478
name: Abnormality of the eye
is_a: HP:0000118

[Term]
id: HP:0000005
name: Mode of inheritance
is_a: HP:0000001

[Term]
id: HP:0000006
name: Autosomal dominant inheritance
is_a: HP:0000005

[Term]
id: HP:0001111
name: Old seizure term
is_obsolete: true
replaced_by: HP:0001250

[Term]
id: HP:0002222
name: Gone term
is_obsolete: true

[Term]
id: HP:0003333
name: Orphan term
is_a: HP:0004444

[Typedef]
id: part_of
name: part of
"""

HEADER = ["database_id", "disease_name", "qualifier", "hpo_id", "reference", "evidence", "onset",
          "frequency", "sex", "modifier", "aspect", "biocuration"]


def _hpoa(rows):
    lines = ["#description: test", "#version: 2026-09-01", "\t".join(HEADER)]
    for r in rows:
        full = dict.fromkeys(HEADER, "")
        full.update(r)
        lines.append("\t".join(full[h] for h in HEADER))
    return "\n".join(lines) + "\n"


ROWS = [
    # D1 (GENEA): seizure 2/3, ID only as NOT, AD inheritance (aspect I)
    {"database_id": "OMIM:1", "disease_name": "D1", "hpo_id": "HP:0001250", "frequency": "2/3", "aspect": "P"},
    {"database_id": "OMIM:1", "disease_name": "D1", "hpo_id": "HP:0001249", "qualifier": "NOT", "aspect": "P"},
    {"database_id": "OMIM:1", "disease_name": "D1", "hpo_id": "HP:0000006", "aspect": "I"},
    # D2 (GENEB): ID 0/4 only → negative; eye Excluded → negative; seizure via alt_id
    {"database_id": "OMIM:2", "disease_name": "D2", "hpo_id": "HP:0001249", "frequency": "0/4", "aspect": "P"},
    {"database_id": "OMIM:2", "disease_name": "D2", "hpo_id": "HP:0000478", "frequency": "HP:0040285", "aspect": "P"},
    {"database_id": "OMIM:2", "disease_name": "D2", "hpo_id": "HP:0009999", "frequency": "HP:0040283", "aspect": "P"},
    # D3 (GENEC): ID 0/2 + 3/5 merged → kept 3/7, conflict (0% vs 60%); obsolete replaced; gone dropped; orphan
    {"database_id": "ORPHA:3", "disease_name": "D3", "hpo_id": "HP:0001249", "frequency": "0/2", "aspect": "P"},
    {"database_id": "ORPHA:3", "disease_name": "D3", "hpo_id": "HP:0001249", "frequency": "3/5", "aspect": "P"},
    {"database_id": "ORPHA:3", "disease_name": "D3", "hpo_id": "HP:0001111", "aspect": "P"},
    {"database_id": "ORPHA:3", "disease_name": "D3", "hpo_id": "HP:0002222", "aspect": "P"},
    {"database_id": "ORPHA:3", "disease_name": "D3", "hpo_id": "HP:0003333", "aspect": "P"},
    # D4 has no gene
    {"database_id": "OMIM:4", "disease_name": "D4", "hpo_id": "HP:0000478", "aspect": "P"},
]

G2D = """ncbi_gene_id\tgene_symbol\tassociation_type\tdisease_id\tsource
NCBIGene:1\tGENEA\tMENDELIAN\tOMIM:1\tx
NCBIGene:2\tGENEB\tMENDELIAN\tOMIM:2\tx
NCBIGene:3\tGENEC\tUNKNOWN\tORPHA:3\tx
"""


@pytest.fixture
def sources(tmp_path):
    (tmp_path / "hp.obo").write_text(OBO, encoding="utf-8")
    (tmp_path / "phenotype.hpoa").write_text(_hpoa(ROWS), encoding="utf-8")
    (tmp_path / "genes_to_disease.txt").write_text(G2D, encoding="utf-8")
    return tmp_path


def _build(src):
    return H.build(src / "hp.obo", src / "phenotype.hpoa", src / "genes_to_disease.txt")


def test_rules_and_propagation(sources):
    r = _build(sources)
    qc, genes = r["qc"], r["term_genes"]
    assert qc["R1_non_P_aspect"] == 1 and qc["R2_not"] == 1 and qc["R3_excluded"] == 1
    assert qc["R4_zero_pairs_to_negative"] == 1 and qc["R4_zero_rows_kept_merged"] == 1
    assert qc["R5_alt_id"] == 1 and qc["R5_replaced_by"] == 1 and qc["R5_dropped_obsolete"] == 1
    assert qc["R6_broken_chain"] == 1 and qc["R7_duplicate_pairs"] == 1 and qc["R7_frequency_conflict_pairs"] == 1
    # ID: GENEA only via NOT, GENEB only via 0/4 → only GENEC
    assert set(genes["HP:0001249"]) == {"GENEC"}
    # Seizure: GENEA (direct), GENEB (alt_id), GENEC (obsolete replaced_by)
    assert set(genes["HP:0001250"]) == {"GENEA", "GENEB", "GENEC"}
    # propagated to ancestors inside the phenotype subtree, not to All / inheritance
    assert set(genes["HP:0000118"]) == {"GENEA", "GENEB", "GENEC"}
    assert "HP:0000001" not in genes and "HP:0000006" not in genes
    assert "HP:0000478" not in genes          # Excluded (GENEB) and gene-less D4 only
    merged = r["merged"][("ORPHA:3", "HP:0001249")]
    assert merged["frequency"] == "3/7" and merged["conflict"]
    assert {k for k in r["negative"]} == {("OMIM:1", "HP:0001249"), ("OMIM:2", "HP:0001249"), ("OMIM:2", "HP:0000478")}


def test_snapshot_activate_and_scorer(sources, tmp_path, monkeypatch):
    kb = tmp_path / "kb"
    r = _build(sources)
    snap = H.write_snapshot(r, kb, {"hp_obo": sources / "hp.obo"})
    assert snap.name == "2026-09-01_2026-09-01"
    for name in ("hpo_gene.tsv", "term_status.tsv", "disease_terms.tsv", "negative_terms.tsv", "manifest.json",
                 "build_qc.md"):
        assert (snap / name).is_file()
    with open(snap / "negative_terms.tsv", encoding="utf-8") as fh:
        reasons = {row["reason"] for row in csv.DictReader(fh, delimiter="\t")}
    assert reasons == {"NOT", "zero_frequency", "excluded"}

    legacy = tmp_path / "phenotype_to_genes.txt"
    legacy.write_text("hpo_id\thpo_name\tncbi_gene_id\tgene_symbol\tdisease_id\n"
                      "HP:0001249\tID\t1\tGENEA\tOMIM:1\nHP:0000006\tAD\t1\tGENEA\tOMIM:1\n", encoding="utf-8")
    monkeypatch.setattr(H, "KB_ROOT", kb)
    monkeypatch.setattr(phenotype_scorer, "PHENO_TO_GENES_PATH", legacy)
    monkeypatch.setattr(phenotype_scorer, "_load_panels", lambda *a, **k: {})
    assert H.active_table().table_id == H.LEGACY_TABLE_ID
    phenotype_scorer.reload_db()
    # legacy: inheritance term carries genes and counts in the total
    assert phenotype_scorer.compute_pheno_score([("HP:0001249", 1), ("HP:0000006", 1)]) == {"GENEA": 100.0}

    H.activate(snap)
    phenotype_scorer.reload_db()
    try:
        assert phenotype_scorer.table_id().endswith("2026-09-01_2026-09-01")
        assert phenotype_scorer.excluded_from_scoring("HP:0000006")
        assert not phenotype_scorer.excluded_from_scoring("HP:0001250")
        # inheritance left out of the total; alt_id normalised; ID no longer gives GENEA
        scores = phenotype_scorer.compute_pheno_score([("HP:0009999", 1), ("HP:0000006", 1), ("HP:0001249", 1)])
        assert scores == {"GENEA": 50.0, "GENEB": 50.0, "GENEC": 100.0}
        assert phenotype_scorer.genes_for_key("HP:0009999")["gene_count"] == 3
        assert phenotype_scorer.gene_count("HP:0009999") == 3

        # pheno_score.tsv stamp
        out = phenotype_scorer.write_pheno_table("S1", {"GENEA": 50.0}, target_dir=tmp_path / "v")
        assert phenotype_scorer.pheno_table_is_current(out)
        out.with_name(phenotype_scorer.PHENO_TABLE_STAMP).unlink()
        assert not phenotype_scorer.pheno_table_is_current(out)     # unstamped = legacy table
    finally:
        assert H.deactivate()
        phenotype_scorer.reload_db()
    assert phenotype_scorer.pheno_table_is_current(out)


def test_frequency_helpers():
    assert H.parse_frequency("3/5") == ("fraction", 3, 5)
    assert H.parse_frequency("12.5%") == ("percent", 12.5)
    assert H.parse_frequency("HP:0040283") == ("code", "HP:0040283")
    assert H.frequency_category(("fraction", 0, 4)) == H.EXCLUDED
    assert H.frequency_category(("percent", 85.0)) == "HP:0040281"
    assert H.frequency_category(("blank",)) is None


def test_diff_report(sources, tmp_path, monkeypatch):
    snap = H.write_snapshot(_build(sources), tmp_path / "kb", {"hp_obo": sources / "hp.obo"})
    legacy = tmp_path / "phenotype_to_genes.txt"
    legacy.write_text("hpo_id\thpo_name\tncbi_gene_id\tgene_symbol\tdisease_id\n"
                      "HP:0001249\tID\t1\tGENEA\tOMIM:1\nHP:0001249\tID\t3\tGENEC\tORPHA:3\n"
                      "HP:0000006\tAD\t1\tGENEA\tOMIM:1\nHP:0001250\tS\t1\tGENEA\tOMIM:1\n", encoding="utf-8")
    monkeypatch.setattr(phenotype_scorer, "_load_panels", lambda *a, **k: {})
    hpo = [{"phenotype": "HP:0001249", "weight": 1}, {"phenotype": "HP:0000006", "weight": 1}]
    cases = {"S1": ("WES", {}, hpo, [])}
    variants = {"S1": ({"v1": ("GENEA", 2.0, "c.1A>G"), "v2": ("GENEC", 0.0, "c.2A>G"),
                        "v3": ("GENEX", 1.0, "c.3A>G")}, {"v1": "1"})}
    report = pheno_table_diff.run(out_dir=tmp_path / "out", old_table=legacy, new_snapshot=snap, sample_ids=["S1"],
                                  case_phenotype=cases.get, case_variants=lambda sid, meta: variants[sid])
    assert report["table"]["pairs_removed"] >= 2
    assert report["excluded_term_use"][0]["hpo_id"] == "HP:0000006"
    m = report["marked"]["WES|1"]
    assert m["in_panel_old"] == 1 and m["in_panel_new"] == 0 and m["lost_panel"] == 1
    assert m["rank_old"]["median"] == 1 and m["rank_new"]["median"] == 2     # GENEC now outranks GENEA
    lost = report["lost_panel_variants"][0]
    assert lost["gene"] == "GENEA" and [c.split()[0] for c in lost["causes"]] == ["HP:0000006", "HP:0001249"]
    text = (tmp_path / "out" / "report.md").read_text(encoding="utf-8")
    assert "S1" not in text and "GENEA" in text
    assert "S1" in (tmp_path / "out" / "marked_variants.tsv").read_text(encoding="utf-8")
    assert json.loads((tmp_path / "out" / "summary.json").read_text(encoding="utf-8"))["cases"]["WES"] == 1


def test_loader_recomputes_pheno_table_after_switch(tmp_path, monkeypatch):
    from app.services import sample_loader
    vdir = tmp_path / "v"
    vdir.mkdir()
    (vdir / "analysis.json").write_text("{}", encoding="utf-8")
    monkeypatch.setattr(sample_loader, "_sidecar_file", lambda sid, d, name, **k: d / name)
    monkeypatch.setattr(phenotype_scorer, "_LOADED", True)
    monkeypatch.setattr(phenotype_scorer, "_HPO_TO_GENES", {"HP:1": {"GENEA"}})
    monkeypatch.setattr(phenotype_scorer, "_NON_PHENOTYPE", set())
    table = {"id": "legacy:phenotype_to_genes.txt"}
    monkeypatch.setattr(phenotype_scorer, "table_id", lambda: table["id"])
    monkeypatch.setattr(phenotype_scorer, "write_pheno_table",
                        lambda sid, scores, target_dir=None, **k: calls.append(dict(scores)) or (target_dir / "x"))
    calls = []
    pheno = vdir / "pheno_score.tsv"
    pheno.write_text("gene_symbol\tpheno_score\nGENEA\t100.0\n", encoding="utf-8")
    import os
    os.utime(pheno, (2e9, 2e9))                       # newer than analysis.json
    sample_loader._ensure_pheno_table("S1", vdir, [("HP:1", 1)], [])
    assert calls == []                                # unstamped file, legacy table → current
    table["id"] = "hpo-gene-clean-v1:2026-09-01_2026-09-01"
    sample_loader._ensure_pheno_table("S1", vdir, [("HP:1", 1)], [])
    assert calls == [{"GENEA": 100.0}]                # table switched → recomputed
