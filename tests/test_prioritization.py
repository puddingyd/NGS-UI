"""Phase 3: In-house prioritization = V′ + P(g, d) + C (plan §5.3, §6, §7)."""
import pytest

from app.services import omim_store
from app.services import prioritization as PR
from app.services.pheno_engine import engine as E

from test_pheno_engine import OMIM_ROWS, snap  # noqa: F401  (fixture)


@pytest.fixture
def omim(monkeypatch):
    rows = {r["gene_symbol"]: r for r in OMIM_ROWS}
    rows["GENEX"] = {"gene_symbol": "GENEX", "OMIM_disease": "X disease (300001)(XLR)"}
    rows["GENEM"] = {"gene_symbol": "GENEM", "OMIM_disease": "Mixed (100010)(AD, AR)"}
    rows["GENER"] = {"gene_symbol": "GENER", "OMIM_disease": "Recessive (100011)(AR)"}
    monkeypatch.setattr(omim_store, "lookup", lambda gene="", omim_id=None: rows.get(gene))
    monkeypatch.setattr(omim_store, "cache_signature", lambda: ("test", len(rows)))
    PR._OMIM_CACHE.update(sig=None, genes={})
    yield rows
    PR._OMIM_CACHE.update(sig=None, genes={})


@pytest.fixture
def eng(snap):  # noqa: F811
    return E.Engine(snap)


def test_pair_points_match_plan_table():
    assert PR.pair_points(10, 10) == pytest.approx(9.0, abs=0.1)
    assert PR.pair_points(7, 7) == pytest.approx(6.0, abs=0.1)
    assert PR.pair_points(10, 4) == pytest.approx(4.0, abs=0.1)
    assert PR.pair_points(4, 4) == pytest.approx(2.8, abs=0.1)
    assert PR.pair_points(3, 8) == PR.pair_points(8, 3)


def _one(var, moi, *, sex="XX", test="WES", P=5.0, others=()):
    dis = [PR.Disease(id="OMIM:1", P=P, moi=set(moi))]
    return PR.evaluate_gene([var, *others], dis, sex=sex, test_type=test)[var.id]


def test_inheritance_compatibility_points():
    het = PR.Var("a", 10, "het")
    assert _one(het, {"dominant"})["c"] == 0 and _one(het, {"dominant"})["total"] == 15
    assert _one(het, {"ar"})["c"] == -3 and _one(het, {"ar"})["mode"] == "single"
    assert _one(het, {"ar"}, test="WGS")["c"] == -4
    assert _one(het, {"ar"}, test="TITAN-WGS")["c"] == -4
    hom = PR.Var("b", 10, "hom")
    assert _one(hom, {"ar"})["c"] == 0 and _one(hom, {"ar"})["mode"] == "hom"
    # X-linked recessive: male 0, female het −3, female hom 0
    assert _one(het, {"xlr"}, sex="XY")["c"] == 0
    assert _one(het, {"xlr"}, sex="XX")["c"] == -3
    assert _one(het, {"xlr"}, sex="")["c"] == -3
    assert _one(hom, {"xlr"}, sex="XX")["c"] == 0
    assert _one(het, {"xld"})["c"] == 0 and _one(het, set())["mode"] == "unknown"
    # AD/AR: most compatible
    assert _one(het, {"dominant", "ar"})["mode"] == "dominant"
    # somatic points never take part in inheritance logic or pairs
    som = PR.Var("s", 10, "het", somatic=True)
    assert _one(som, {"ar"})["c"] == 0


def test_recessive_pair_and_best_partner():
    a, b, c = PR.Var("a", 10, "het"), PR.Var("b", 4, "het"), PR.Var("c", 1, "het")
    res = PR.evaluate_gene([a, b, c], [PR.Disease(id="OMIM:1", P=5, moi={"ar"})], sex="XX", test_type="WGS")
    # a: pair with b (4 + −1) beats single (10 − 4)? pair 4.0 − 1 = 3.0 < 6 → single
    assert res["a"]["mode"] == "single" and res["a"]["total"] == pytest.approx(11)
    # b pairs with the strongest other het (a): 4.0 − 1 > 4 − 4
    assert res["b"]["mode"] == "pair" and res["b"]["partner"] == "a"
    assert res["b"]["v"] == pytest.approx(PR.pair_points(4, 10))
    assert res["c"]["partner"] == "a"
    two_lp = PR.evaluate_gene([PR.Var("x", 7, "het"), PR.Var("y", 7, "het")],
                              [PR.Disease(id="d", P=0, moi={"ar"})], sex="XX", test_type="WGS")
    assert two_lp["x"]["mode"] == "pair" and two_lp["x"]["c"] == -1
    assert two_lp["x"]["total"] == pytest.approx(PR.pair_points(7, 7) - 1)


def test_best_disease_by_total(omim, eng):
    variants = {"v1": {"gene_symbol": "GENEA", "zygosity": "heterozygous", "prio_v_points": 6}}
    info = PR.annotate(variants, hpo=[{"phenotype": "HP:0007359", "weight": 1}], sex="F", test_type="WES",
                       engine=eng)
    assert info["mode"] == "prio" and info["terms"] == 1
    prio = variants["v1"]["prio"]
    # GENEA: focal epilepsy A (AD) and broad syndrome C (AR); patient has focal seizures → AD disease wins
    assert prio["disease_id"] == "OMIM:100001" and prio["mode"] == "dominant" and prio["c"] == 0
    rows = {r["disease_id"]: r for r in PR.pheno_result([{"phenotype": "HP:0007359"}], engine=eng)[1]
            .gene_disease_points("GENEA")}
    assert prio["p"] == pytest.approx(rows["OMIM:100001"]["P"])
    assert prio["total"] == pytest.approx(6 + prio["p"])
    assert prio["rank"] == 1 and variants["v1"]["old_rank"] is None


def test_no_disease_gene_gets_floor_and_ranks(omim, eng):
    variants = {
        "a": {"gene_symbol": "GENEB", "zygosity": "het", "prio_v_points": 2, "total_score": 60},
        "b": {"gene_symbol": "NODISEASE", "zygosity": "het", "prio_v_points": 10, "total_score": 80},
        "c": {"gene_symbol": "GENED", "zygosity": "het", "prio_v_points": 3, "total_score": 80},
    }
    PR.annotate(variants, hpo=["HP:0002197"], sex="M", test_type="WGS", engine=eng)
    assert variants["b"]["prio"]["p"] == PR.P_FLOOR and variants["b"]["prio"]["no_disease"]
    assert variants["c"]["prio"]["mode"] == "single" and variants["c"]["prio"]["c"] == -4
    order = sorted(variants, key=lambda i: PR.sort_key(variants[i], i, "prio"))
    assert order[0] == "a"                                  # matching phenotype beats a stronger variant
    assert variants["a"]["prio"]["rank"] == 1
    # old ranks: competition ranking on total_score
    assert variants["b"]["old_rank"] == 1 and variants["c"]["old_rank"] == 1 and variants["a"]["old_rank"] == 3


def test_panel_only_falls_back_to_old(omim, eng):
    variants = {"a": {"gene_symbol": "GENEA", "zygosity": "het", "prio_v_points": 4, "total_score": 70}}
    info = PR.annotate(variants, hpo=[], sex="", test_type="WES", engine=eng)
    assert info["mode"] == "old" and "panel" in info["reason"]
    prio = variants["a"]["prio"]
    assert prio["p"] is None and prio["pheno_counted"] is False
    assert prio["total"] == pytest.approx(4)               # AD disease available → most compatible, C = 0
    info = PR.annotate(variants, hpo=["HP:0000006"], engine=eng)      # inheritance term only → removed
    assert info["mode"] == "old"


def test_ranking_points_fallback():
    assert PR.ranking_points({"prio_v_points": 3.5}) == (3.5, False)
    assert PR.ranking_points({"effective_rank_score": "2"}) == (2.0, False)
    assert PR.ranking_points({}) == (0.0, True)


def test_resort_and_determinism(omim, eng):
    def make():
        return {
            "a": {"gene_symbol": "GENEA", "zygosity": "het", "prio_v_points": 1, "total_score": 90},
            "b": {"gene_symbol": "GENEB", "zygosity": "het", "prio_v_points": 1, "total_score": 10},
        }
    v1, v2 = make(), make()
    PR.annotate(v1, hpo=["HP:0002197"], engine=eng)
    PR.annotate(v2, hpo=["HP:0002197"], engine=eng)
    assert v1 == v2
    cats = {"2": ["a", "b"]}
    PR.resort(cats, v1, "prio")
    assert cats["2"] == ["b", "a"]
    PR.resort(cats, v1, "old")
    assert cats["2"] == ["a", "b"]


def test_norm_helpers():
    assert PR.norm_sex("XY", "F") == "XY"
    assert PR.norm_sex("", "女") == "XX"
    assert PR.norm_sex("XXY", "") == ""
    assert PR.norm_zygosity("Hemizygous") == "hemi"
    assert PR.moi_classes({"XL"}) == {"xlr", "xld"}
    assert PR.moi_classes({"SMu"}) == set()


def test_loader_applies_prioritization_and_sorts(omim, eng, monkeypatch):
    from app.services import ploidy, sample_loader
    monkeypatch.setattr(ploidy, "load_sample_ploidy", lambda sid: {"karyotype": "XX"})
    monkeypatch.setattr(E, "load_engine", lambda folder=None: eng)
    PR._PHENO_CACHE.clear()
    variants = {
        "a": {"gene_symbol": "GENEA", "zygosity": "het", "prio_v_points": 1, "total_score": 90},
        "b": {"gene_symbol": "GENEB", "zygosity": "het", "prio_v_points": 1, "total_score": 10},
    }
    cats = {"2": ["a", "b"], "1A": []}
    info = sample_loader._apply_prioritization(variants, cats, sample_id="S1", meta={"sex": "F"},
                                               hpo_list=[{"phenotype": "HP:0002197", "weight": 1}],
                                               test_type="WES")
    assert info["mode"] == "prio" and info["sex"] == "XX"
    assert cats["2"] == ["b", "a"]
    # panel only → old order, P not counted
    info = sample_loader._apply_prioritization(variants, cats, sample_id="S1", meta={},
                                               hpo_list=[], test_type="WES")
    assert info["mode"] == "old" and cats["2"] == ["a", "b"] and variants["a"]["prio"]["p"] is None
    # a broken engine never breaks the load
    monkeypatch.setattr(PR, "annotate", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("x")))
    info = sample_loader._apply_prioritization(variants, cats, sample_id="S1", meta={}, hpo_list=["HP:1"],
                                               test_type="WES")
    assert info["mode"] == "old" and "失敗" in info["reason"]


def test_gene_candidates_one_row_per_gene_and_pairs_once():
    a, b, r = PR.Var("a", 7, "het"), PR.Var("b", 7, "het"), PR.Var("r", 10, "ref")
    dis = [PR.Disease(id="OMIM:1", name="AR disease", P=8, moi={"ar"}),
           PR.Disease(id="OMIM:2", name="AD disease", P=2, moi={"dominant"})]
    cand = PR.gene_candidates([a, b, r], dis, sex="XX", test_type="WES")
    assert cand["mode"] == "pair" and cand["variants"] in (["a", "b"], ["b", "a"])
    assert cand["total"] == pytest.approx(PR.pair_points(7, 7) - 1 + 8)
    assert set(cand["v_each"]) == {"a", "b"}
    # the pair appears once; ref calls never form a candidate
    idents = [(o["disease_id"], o["mode"], tuple(sorted(o["variants"]))) for o in [cand, *cand["others"]]]
    assert len(idents) == len(set(idents))
    assert all("r" not in o["variants"] for o in [cand, *cand["others"]])
    assert PR.gene_candidates([r], dis, sex="XX", test_type="WES") is None
    assert PR.gene_candidates([a, b, r], dis, sex="XX", test_type="WES") == cand      # deterministic


def test_case_candidate_list(omim, eng):
    variants = {
        "a": {"gene_symbol": "GENEB", "zygosity": "het", "prio_v_points": 2},
        "b": {"gene_symbol": "GENEA", "zygosity": "het", "prio_v_points": 4},
        "c": {"gene_symbol": "GENEA", "zygosity": "het", "prio_v_points": 3},
        "d": {"gene_symbol": "NODISEASE", "zygosity": "hom", "prio_v_points": 10},
    }
    info = PR.annotate(variants, hpo=["HP:0002197"], sex="F", test_type="WES", engine=eng)
    cands = info["candidates"]
    assert info["candidate_count"] == 3 and [c["gene"] for c in cands].count("GENEA") == 1
    assert cands[0]["rank"] == 1 and cands == sorted(cands, key=lambda c: -c["total"])
    assert all(c["rank"] <= c2["rank"] for c, c2 in zip(cands, cands[1:]))
    nod = next(c for c in cands if c["gene"] == "NODISEASE")
    assert nod["p"] == PR.P_FLOOR and nod["no_disease"]
    # panel only: no list
    info = PR.annotate(variants, hpo=[], engine=eng)
    assert "candidates" not in info
    # search / secondary (rank=False): no list
    assert "candidates" not in PR.annotate(variants, hpo=["HP:0002197"], engine=eng, rank=False)


def test_candidate_limit():
    genes = [{"gene": f"G{i}", "total": float(i % 7), "p": 1.0, "v": 1.0} for i in range(30)]
    shown, total = PR.rank_candidates(genes, limit=10)
    assert total == 30 and len(shown) == 10 and shown[0]["rank"] == 1
    assert [g["rank"] for g in shown[:5]] == [1, 1, 1, 1, 5]


def test_floor_matches_engine():
    assert PR.P_FLOOR == E.P_MIN == E.gene_points(5000)
