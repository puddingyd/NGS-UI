"""Stage 6: literature validation helpers (plan §11) on a small synthetic knowledge base."""
import json

import numpy as np
import pytest

from app.services import hpo_update as U
from app.services.pheno_engine import benchmark as B
from app.services.pheno_engine import engine as E

from test_pheno_engine import HEADER, MONDO, OBO, OMIM_ROWS

ANN = {   # disease: (name, [(term, reference)])
    "OMIM:100001": ("Focal epilepsy A", [("HP:0007359", "PMID:111"), ("HP:0001249", "PMID:222")]),
    "OMIM:100002": ("Generalized epilepsy B", [("HP:0002197", "PMID:333")]),
    "OMIM:100003": ("Broad syndrome C", [("HP:0001250", "PMID:444"), ("HP:0001249", "PMID:444"),
                                         ("HP:0000316", "PMID:444")]),
    "OMIM:100004": ("Eye disease D", [("HP:0000316", "PMID:555")]),
}


@pytest.fixture(scope="module")
def kb(tmp_path_factory):
    src = tmp_path_factory.mktemp("bench")
    lines = ["#description: test", "#version: 2026-09-01", "\t".join(HEADER)]
    for disease, (name, terms) in ANN.items():
        for t, ref in terms:
            row = dict.fromkeys(HEADER, "")
            row.update(database_id=disease, disease_name=name, hpo_id=t, aspect="P", reference=ref, evidence="PCS")
            lines.append("\t".join(row[h] for h in HEADER))
    (src / "phenotype.hpoa").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (src / "hp.obo").write_text(OBO, encoding="utf-8")
    (src / "genes_to_disease.txt").write_text("ncbi_gene_id\tgene_symbol\tassociation_type\tdisease_id\tsource\n",
                                              encoding="utf-8")
    (src / "mondo.obo").write_text(MONDO, encoding="utf-8")
    folder, _, _ = U.build_snapshot(src / "hp.obo", src / "phenotype.hpoa", src / "genes_to_disease.txt",
                                    kb_root=src / "kb", mondo=src / "mondo.obo", allow_missing={"orphanet"},
                                    omim=(OMIM_ROWS, {"source": "test"}), progress=lambda *a: None)
    cases = src / "cases"
    packets = [
        ("GENEA", "case_a1", ["HP:0011153", "HP:0001249", "HP:0000316"], "OMIM:100001", "GENEA", "PMID:111"),
        ("GENEA", "case_a2", ["HP:0007359", "HP:0001249", "HP:0001250"], "OMIM:100001", "GENEA", "PMID:999"),
        ("GENEB", "case_b1", ["HP:0002197", "HP:0001250", "HP:0001249"], "OMIM:100002", "GENEB", "PMID:888"),
        ("GENEX", "case_x1", ["HP:0002197", "HP:0001250", "HP:0001249"], "OMIM:999999", "GENEX", "PMID:777"),
    ]
    for cohort, cid, terms, dis, gene, pmid in packets:
        (cases / cohort).mkdir(parents=True, exist_ok=True)
        feats = [{"type": {"id": t}} for t in terms] + [{"type": {"id": "HP:0000478"}, "excluded": True}]
        data = {"id": cid, "phenotypicFeatures": feats, "diseases": [{"term": {"id": dis}}],
                "interpretations": [{"diagnosis": {"disease": {"id": dis}, "genomicInterpretations": [
                    {"variantInterpretation": {"variationDescriptor": {"geneContext": {"symbol": gene}}}}]}}],
                "metaData": {"externalReferences": [{"id": pmid}]}}
        (cases / cohort / f"{cid}.json").write_text(json.dumps(data), encoding="utf-8")
    return folder, src


def test_load_cases_and_timesplit(kb):
    folder, src = kb
    eng = E.Engine(folder)
    cases = {c.id: c for c in B.load_cases(src / "cases")}
    a1 = cases["case_a1"]
    assert a1.cohort == "GENEA" and a1.genes == ["GENEA"] and a1.pmids == frozenset({"111"})
    assert "HP:0000478" not in a1.obs and a1.dis == ["OMIM:100001"]
    cited = B.hpoa_citations(src / "phenotype.hpoa")
    assert cited["OMIM:100001"] == {"111", "222"}
    elig = {c.id for c in B.timesplit_eligible(cases.values(), eng, cited)}
    # a1 is cited by its disease's annotations; x1's disease is not in the knowledge base
    assert elig == {"case_a2", "case_b1"}


def test_select_is_deterministic_and_capped(kb):
    _, src = kb
    cases = B.load_cases(src / "cases")
    one = [c.id for c in B.select(cases, per_cohort=1, seed=5)]
    assert one == [c.id for c in B.select(cases, per_cohort=1, seed=5)]
    assert len(one) == 3 and sum(i.startswith("case_a") for i in one) == 1
    assert "case_b1" not in [c.id for c in B.select(cases, exclude_cohorts={"GENEB"})]


def test_evaluate_methods_coverage_and_lopo(kb):
    folder, src = kb
    eng = E.Engine(folder)
    cases = B.load_cases(src / "cases")
    rows, cov = B.evaluate(eng, cases)
    assert cov["cases"] == 4 and cov["disease_missing"] == 1 and len(rows) == 3
    a2 = next(r for r in rows if r["case"] == "case_a2")
    assert a2["engine|d"] == 1 and a2["engine|g"] == 1
    for m in ("nobaseline", "resnik", "platform"):
        assert a2[f"{m}|g"] >= 1
    # LOPO: case_a1 cites PMID:111, the source of OMIM:100001's focal-seizure annotation
    specs = B.lopo_specs(eng, cases, B.disease_annotation_refs(folder))
    spec = specs["case_a1"]
    assert spec.removed == 1 and spec.total == 2
    B.lopo_baselines(eng, specs)
    assert np.isfinite(spec.mu) and spec.mu != eng.mu[spec.disease]
    full = next(r for r in rows if r["case"] == "case_a1")
    lopo, _ = B.evaluate(eng, [c for c in cases if c.id == "case_a1"], lopo=specs)
    assert lopo[0]["engine|d"] >= full["engine|d"] and lopo[0]["lopo_removed"] == 1
    # the unmodified engine result is unaffected by the LOPO run
    again, _ = B.evaluate(eng, [c for c in cases if c.id == "case_a1"])
    assert again[0]["engine|d"] == full["engine|d"]


def test_statistics():
    assert B.mcnemar_exact(0, 5) == pytest.approx(0.0625)
    assert B.mcnemar_exact(3, 3) == 1.0 and B.mcnemar_exact(0, 0) == 1.0
    m = B.metrics([1, 2, 20, None])
    assert m["n"] == 3 and m["top1"] == pytest.approx(1 / 3) and m["mrr"] == pytest.approx((1 + .5 + .05) / 3)
    rows = [{"cohort": f"c{i}", "a": 1.0, "b": 20.0} for i in range(6)] + [{"cohort": "c9", "a": 30.0, "b": 2.0}]
    p = B.paired(rows, "a", "b", boot=200)
    assert p["n"] == 7 and p["better"] == 6 and p["worse"] == 1 and p["delta_mrr"] > 0
    assert p["ci"][0] <= p["delta_mrr"] <= p["ci"][1]


def test_calibrate_points():
    ranks = [1] * 50 + [3] * 30 + [20] * 15 + [2000] * 5
    rows = B.calibrate(ranks, 4888, cohorts=[f"c{i % 10}" for i in range(len(ranks))], boot=100)
    first = rows[0]
    assert first["bin"] == "1" and first["k"] == 50
    assert first["points"] == pytest.approx(round(np.log(0.5 / (1 / 4888)) / np.log(2.08), 1))
    assert rows[-1]["bin"] == ">1000" and rows[-1]["k"] == 5
    assert rows[1]["points"] > rows[3]["points"]
    assert rows[4]["points"] is None and rows[4]["k"] == 0          # empty bin


def test_external_ranks(tmp_path):
    (tmp_path / "case_a2.tsv").write_text("rank\tgene_symbol\tdisease_id\n1\tGENEB\tOMIM:100002\n2\tGENEA\tOMIM:100001\n",
                                          encoding="utf-8")
    ext = B.external_ranks(tmp_path)
    assert ext["case_a2"]["genes"] == ["GENEB", "GENEA"]
    assert B._list_rank(ext["case_a2"]["genes"], "GENEA", 10) == 2
    assert B._list_rank(ext["case_a2"]["genes"], "GENEZ", 10) == 2 + (8 + 1) / 2


def test_duplicate_gene_map_rows(kb):
    folder, _ = kb
    links = [{"disease_id": "OMIM:100001", "gene_symbol": "GENEA", "source": "g2d"},
             {"disease_id": "OMIM:100001", "gene_symbol": "GENEA", "source": "omim"},
             {"disease_id": "OMIM:100001", "gene_symbol": "GENEA", "source": "g2d"}]
    eng = E.Engine(folder, links=links)
    assert eng.genes == ["GENEA"] and eng.link_info[("OMIM:100001", "GENEA")]["source"] == "g2d,omim"
