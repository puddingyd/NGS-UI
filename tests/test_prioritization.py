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
    # a: pair with b (≈4.0, C 0) vs single (10 − 4 = 6) → single
    assert res["a"]["mode"] == "single" and res["a"]["total"] == pytest.approx(11)
    # b pairs with the strongest other het (a): 4.0 + 0 > 4 − 4
    assert res["b"]["mode"] == "pair" and res["b"]["partner"] == "a"
    assert res["b"]["v"] == pytest.approx(PR.pair_points(4, 10))
    assert res["c"]["partner"] == "a"
    two_lp = PR.evaluate_gene([PR.Var("x", 7, "het"), PR.Var("y", 7, "het")],
                              [PR.Disease(id="d", P=0, moi={"ar"})], sex="XX", test_type="WGS")
    assert two_lp["x"]["mode"] == "pair" and two_lp["x"]["c"] == 0
    assert two_lp["x"]["total"] == pytest.approx(PR.pair_points(7, 7))


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
    assert cand["total"] == pytest.approx(PR.pair_points(7, 7) + 8)
    assert set(cand["v_each"]) == {"a", "b"}
    # the pair appears once; ref calls never form a candidate
    idents = [(o["disease_id"], o["mode"], tuple(sorted(o["variants"]))) for o in [cand, *cand["others"]]]
    assert len(idents) == len(set(idents))
    assert all("r" not in o["variants"] for o in [cand, *cand["others"]])
    assert PR.gene_candidates([r], dis, sex="XX", test_type="WES") is None
    assert PR.gene_candidates([a, b, r], dis, sex="XX", test_type="WES") == cand      # deterministic
    assert cand["same_gene"] == []


def test_gene_candidates_list_other_hets_of_recessive_gene():
    # LP + VUS: the single het (8 − 3) beats the pair, but the VUS is still listed
    lp, vus, lb, ref = PR.Var("lp", 8, "het"), PR.Var("vus", 3.2, "het"), PR.Var("lb", -4, "het"), PR.Var("r", 9, "ref")
    ar = [PR.Disease(id="OMIM:1", name="AR disease", P=-0.3, moi={"ar"})]
    cand = PR.gene_candidates([lp, vus, lb, ref], ar, sex="XX", test_type="WES")
    assert cand["mode"] == "single" and cand["variants"] == ["lp"]
    assert cand["same_gene"] == ["vus", "lb"]            # strongest first, benign kept, ref never
    # a pair row lists the third variant only
    two = PR.gene_candidates([PR.Var("x", 8, "het"), PR.Var("y", 8, "het"), lb], ar, sex="XX", test_type="WES")
    assert two["mode"] == "pair" and two["same_gene"] == ["lb"]
    # dominant and homozygous rows do not list other variants
    ad = [PR.Disease(id="OMIM:2", name="AD disease", P=1, moi={"dominant"})]
    assert PR.gene_candidates([lp, vus], ad, sex="XX", test_type="WES")["same_gene"] == []
    hom = PR.gene_candidates([PR.Var("h", 8, "hom"), vus], ar, sex="XX", test_type="WES")
    assert hom["mode"] == "hom" and hom["same_gene"] == []


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


def _tx(csq, hgvs, ttype="protein_coding"):
    return {"Consequence": csq, "HGVS_C": hgvs, "transcript_type": ttype}


def test_type_prior_classes():
    utr = {"transcript_options": [_tx("5_prime_UTR_variant", "NM_1:c.-22G>T")], "SpliceAI_score": 0.02}
    tp = PR.type_prior(utr)
    assert tp["class"] == "noncoding" and tp["points"] == PR.TYPE_PRIOR["noncoding"] < 0
    assert PR.type_prior(dict(utr, SpliceAI_score=0.3)) is None            # splicing predicted: no prior
    assert PR.type_prior({"transcript_options": utr["transcript_options"]}) is None   # no prediction
    assert PR.type_prior(dict(utr, SpliceAI_score=None, Pangolin_score=-0.05))["class"] == "noncoding"
    assert PR.type_prior(dict(utr, Pangolin_score=-0.4)) is None           # |Pangolin| counts
    miss = {"transcript_options": [_tx("missense_variant", "NM_1:c.23G>T")], "SpliceAI_score": 0}
    assert PR.type_prior(miss) is None
    # protein changed in any coding transcript → nothing; non-coding transcripts are ignored
    mixed = {"transcript_options": [_tx("5_prime_UTR_variant", "NM_1:c.-5A>G"), _tx("missense_variant", "NM_2:c.4A>G")],
             "Pangolin_score": 0}
    assert PR.type_prior(mixed) is None
    lnc = {"transcript_options": [_tx("3_prime_UTR_variant", "NM_1:c.*40A>G"),
                                  _tx("non_coding_transcript_exon_variant", "NR_1:n.40A>G", "lncRNA")],
           "Pangolin_score": 0}
    assert PR.type_prior(lnc)["class"] == "noncoding"
    syn = {"Consequence": "splice_region_variant&synonymous_variant", "HGVS_C": "c.300G>A", "SpliceAI_score": 0.01}
    assert PR.type_prior(syn)["class"] == "synonymous"
    near = {"Consequence": "splice_region_variant&intron_variant", "HGVS_C": "c.123-5T>C", "SpliceAI_score": 0.0}
    deep = {"Consequence": "intron_variant", "HGVS_C": "c.123+40T>C", "SpliceAI_score": 0.0}
    assert PR.type_prior(near)["class"] == "intron_near" and PR.type_prior(deep)["class"] == "noncoding"
    assert PR.type_prior({"Consequence": "intron_variant", "HGVS_C": "c.123+2T>C", "SpliceAI_score": 0}) is None
    # several transcripts: the smallest deduction wins
    two = {"transcript_options": [_tx("intron_variant", "NM_1:c.200-8A>G"), _tx("5_prime_UTR_variant", "NM_2:c.-12A>G")],
           "SpliceAI_score": 0}
    assert PR.type_prior(two)["class"] == "intron_near"
    assert set(PR.TYPE_PRIOR) == {"synonymous", "intron_near", "noncoding"}       # UTR + deep intron pooled


def test_type_prior_breaks_the_kcnq2_tie(omim, eng):
    """A missense (V′ 1.984 after the P-KNN LLR) must beat a 5′UTR variant with the same ACMG points."""
    variants = {
        "mis": {"gene_symbol": "GENEA", "zygosity": "het", "prio_v_points": 1.984, "SpliceAI_score": 0.02,
                "transcript_options": [_tx("missense_variant", "NM_1:c.23G>T")]},
        "utr": {"gene_symbol": "GENEA", "zygosity": "het", "prio_v_points": 2, "SpliceAI_score": 0.02,
                "transcript_options": [_tx("5_prime_UTR_variant", "NM_1:c.-22G>T")]},
    }
    info = PR.annotate(variants, hpo=["HP:0007359"], sex="F", test_type="WES", engine=eng)
    assert variants["utr"]["prio"]["type_prior"]["class"] == "noncoding"
    assert variants["utr"]["prio"]["v"] == pytest.approx(2 + PR.TYPE_PRIOR["noncoding"])
    assert variants["mis"]["prio"]["type_prior"] is None
    assert variants["mis"]["prio"]["v"] == pytest.approx(1.984)
    cand = next(c for c in info["candidates"] if c["gene"] == "GENEA")
    assert cand["variants"] == ["mis"]


def test_candidate_limit():
    genes = [{"gene": f"G{i}", "total": float(i % 7), "p": 1.0, "v": 1.0} for i in range(30)]
    shown, total = PR.rank_candidates(genes, limit=10)
    assert total == 30 and len(shown) == 10 and shown[0]["rank"] == 1
    assert [g["rank"] for g in shown[:5]] == [1, 1, 1, 1, 5]


def test_floor_matches_engine():
    assert PR.P_FLOOR == E.P_MIN == E.gene_points(5000)


# ── CNV / SV and mtDNA (phase 5b / 5c) ─────────────────────────────────────
def test_annotsv_points_piecewise():
    assert PR.annotsv_points(1.0) == 10 and PR.annotsv_points(0.99) == 10
    assert PR.annotsv_points(0.90) == 6 and PR.annotsv_points(0.98) == 9
    assert PR.annotsv_points(0.0) == 0 and PR.annotsv_points(0.89) == 5
    assert PR.annotsv_points(-0.90) == -1 and PR.annotsv_points(-0.98) == -6
    assert PR.annotsv_points(-1.2) == -7 and PR.annotsv_points(None) is None
    assert PR.annotsv_points(0.45) == pytest.approx(0.45 / 0.89 * 5, abs=1e-3)


def test_event_zygosity():
    assert PR.event_zygosity({"sv_type": "DEL", "copy_number": 0, "CHROM": "chr1"}, "XX") == "hom"
    assert PR.event_zygosity({"sv_type": "DEL", "copy_number": 1, "CHROM": "chr1"}, "XY") == "het"
    assert PR.event_zygosity({"sv_type": "DUP", "CHROM": "chrX", "POS": 5_000_000, "END": 6_000_000}, "XY") == "hemi"
    assert PR.event_zygosity({"sv_type": "DEL", "CHROM": "chrX", "POS": 20_000, "END": 30_000}, "XY") == "het"


def test_cnv_events_scored_ranked_and_ar_genes(omim, eng):
    hpo = [{"phenotype": "HP:0007359", "weight": 1}]
    events = {
        "e1": {"sv_type": "DEL", "copy_number": 1, "CHROM": "chr1", "POS": 1, "END": 9, "ranking_score": 0.99,
               "_prio_genes": [("GENEA", "functional"), ("GENED", "functional")]},
        "e2": {"sv_type": "DUP", "copy_number": 3, "CHROM": "chr2", "POS": 1, "END": 9, "ranking_score": 0.0,
               "_prio_genes": []},
        "e3": {"sv_type": "DEL", "copy_number": 0, "CHROM": "chr3", "POS": 1, "END": 9, "ranking_score": 0.45,
               "_prio_genes": [("GENED", "unknown")]},
    }
    cats = {"CNV-1A": ["e2", "e3", "e1"]}
    info = PR.annotate_events(events, hpo=hpo, sex="XX", test_type="WES", source="cnv", engine=eng)
    PR.resort_events(cats, events, info["mode"])
    assert info["mode"] == "prio" and all("_prio_genes" not in e for e in events.values())
    e1 = events["e1"]["prio"]
    # GENEA's AD disease matches the phenotype: V′ 10 + P, C 0
    assert e1["gene"] == "GENEA" and e1["mode"] == "dominant" and e1["v"] == 10 and e1["c"] == 0
    assert e1["rank"] == 1 and cats["CNV-1A"][0] == "e1"
    # heterozygous deletion: recessive genes (GENEA's AR disease, GENED) kept for SNV pairing
    assert {g["gene"] for g in e1["ar_genes"]} == {"GENEA", "GENED"}
    assert all(d["moi_class"] == "ar" for g in e1["ar_genes"] for d in g["diseases"])
    # no affected gene → lowest phenotype points
    assert events["e2"]["prio"]["p"] == PR.P_FLOOR and events["e2"]["prio"]["mode"] == "no_gene"
    # homozygous deletion of an AR gene: C 0, not paired with SNVs
    e3 = events["e3"]["prio"]
    assert e3["mode"] == "hom" and e3["c"] == 0 and "ar_genes" not in e3
    # SV source never offers SNV pairing; no HPO → old mode, P not counted
    sv = {"s": {"sv_type": "DEL", "copy_number": 1, "CHROM": "chr1", "POS": 1, "END": 2, "ranking_score": 0.9,
                "_prio_genes": [("GENED", "functional")]}}
    info = PR.annotate_events(sv, hpo=[], source="sv", engine=eng)
    assert info["mode"] == "old" and sv["s"]["prio"]["p"] is None and "ar_genes" not in sv["s"]["prio"]


def test_mito_points_and_scores(omim, eng):
    assert PR.mito_points({"clinvar_pathogenic": True, "mitomap_pathogenic": True})[0] == 10
    assert PR.mito_points({"clinvar_pathogenic": True})[0] == 7
    assert PR.mito_points({"mitomap_pathogenic": True})[0] == 7
    assert PR.mito_points({"clinvar_benign": True, "mitomap_pathogenic": True})[0] == -7
    assert PR.mito_points({"gnomad_mito_af_max": 0.02, "mitomap_reported": True})[0] == -4
    assert PR.mito_points({"CLNSIG": "Uncertain_significance"})[0] == 2
    assert PR.mito_points({"gnomad_mito_af_max": 0.001})[0] == 1
    variants = {"m1": {"gene_symbol": "MT-ND4", "clinvar_pathogenic": True},
                "m2": {"gene_symbol": "MT-ND5", "gnomad_mito_af_max": 0.3}}
    cats = {"MITO-3": ["m2", "m1"]}
    info = PR.annotate_mito(variants, hpo=["HP:0007359"], engine=eng)
    PR.resort_events(cats, variants, info["mode"])
    assert info["mode"] == "prio" and variants["m1"]["prio"]["c"] == 0
    assert variants["m1"]["prio"]["p"] == PR.P_FLOOR          # gene without a ranked disease
    assert variants["m1"]["prio"]["rank"] == 1 and cats["MITO-3"] == ["m1", "m2"]


def test_explain_diseases(omim, eng):
    out = PR.explain_diseases(["HP:0011153", "HP:0000316"], ["OMIM:100001", "OMIM:100004", "OMIM:404", "OMIM:100001"],
                              engine=eng)
    assert sorted(t["hpo_id"] for t in out["terms"]) == ["HP:0000316", "HP:0011153"]
    assert set(out["diseases"]) == {"OMIM:100001", "OMIM:100004"}       # unknown and duplicate IDs dropped
    focal = {r["hpo_id"]: r for r in out["diseases"]["OMIM:100001"]["rows"]}
    # focal motor seizure under the disease's focal-onset seizure annotation: explained
    assert focal["HP:0011153"]["relation"] == "patient_specific" and focal["HP:0011153"]["match"] == "match"
    assert focal["HP:0011153"]["annotation"] == "HP:0007359"
    eye = {r["hpo_id"]: r for r in out["diseases"]["OMIM:100004"]["rows"]}
    assert eye["HP:0000316"]["relation"] == "same" and eye["HP:0000316"]["match"] == "match"
    assert eye["HP:0011153"]["match"] != "match"
    assert out["diseases"]["OMIM:100004"]["matched"] == 1
    assert PR.match_class("collateral", 0.4) == "partial" and PR.match_class("collateral", -0.1) == "none"
    assert PR.explain_diseases([], ["OMIM:100001"], engine=eng)["diseases"] == {}


def test_lineage_steps_and_graph(omim, eng):
    from app.services.pheno_engine import lineage as L
    lin = L.lineage(eng)
    t = eng.tid
    assert lin.parents(t["HP:0011153"]) == (t["HP:0007359"],)
    assert lin.path_up(t["HP:0011153"], t["HP:0001250"]) == [t["HP:0011153"], t["HP:0007359"], t["HP:0001250"]]
    assert lin.path_to_organ(t["HP:0011153"])[-1] == t["HP:0000707"] and lin.is_organ(t["HP:0000478"])
    out = PR.explain_diseases(["HP:0011153", "HP:0002197", "HP:0000316"], ["OMIM:100001"], engine=eng)
    rows = {r["hpo_id"]: r for r in out["diseases"]["OMIM:100001"]["rows"]}
    focal = rows["HP:0011153"]                       # disease annotated with the parent term: up 1
    assert (focal["relation"], focal["steps_patient"], focal["steps_annotation"]) == ("patient_specific", 1, 0)
    gen = rows["HP:0002197"]                         # collateral: both sides climb to the meeting term
    assert gen["relation"] == "collateral" and gen["mica"]
    assert gen["steps_patient"] == len(lin.path_up(t["HP:0002197"], t[gen["mica"]])) - 1 >= 1
    assert gen["steps_annotation"] == len(lin.path_up(t[gen["annotation"]], t[gen["mica"]])) - 1 >= 1
    eye = rows["HP:0000316"]
    assert eye["relation"] == "none" and eye["steps_patient"] is None and eye["mica"] == ""
    broad = PR.explain_diseases(["HP:0001250"], ["OMIM:100001"], engine=eng)["diseases"]["OMIM:100001"]["rows"][0]
    assert (broad["relation"], broad["steps_patient"], broad["steps_annotation"]) == ("patient_broader", 0, 1)

    g = PR.disease_graph(["HP:0011153", "HP:0002197", "HP:0000316"], "OMIM:100001", engine=eng)
    nodes = {n["id"]: n for n in g["nodes"]}
    assert "patient" in nodes["HP:0011153"]["roles"] and "annotation" in nodes["HP:0007359"]["roles"]
    assert "organ" in nodes["HP:0000707"]["roles"]
    assert "organ" in nodes["HP:0000478"]["roles"]                   # unmatched term still drawn up to its organ
    assert ["HP:0011153", "HP:0007359"] in g["edges"] and ["HP:0001250", "HP:0000707"] in g["edges"]
    assert all(e[0] in nodes and e[1] in nodes for e in g["edges"])
    pair = next(p for p in g["pairs"] if p["hpo_id"] == "HP:0002197")
    assert pair["path_patient"][0] == "HP:0002197" and pair["path_patient"][-1] == pair["mica"] == gen["mica"]
    assert pair["path_annotation"][0] == gen["annotation"] and pair["path_annotation"][-1] == pair["mica"]
    assert "meet" in nodes[pair["mica"]]["roles"]
    assert "meet" not in nodes["HP:0007359"]["roles"]             # direct parent/child pair: no meeting star
    assert pair["relation_zh"] and pair["match"] in ("match", "partial", "none")
    g2 = PR.disease_graph(["HP:0011153"], "OMIM:100001", engine=eng)
    assert "HP:0001249" not in {n["id"] for n in g2["nodes"]} and g2["others_total"] == 1   # hidden by default
    g2 = PR.disease_graph(["HP:0011153"], "OMIM:100001", all_annotations=True, engine=eng)
    assert "other_annotation" in {r for n in g2["nodes"] for r in n["roles"]} and g2["others_shown"] == 1
    assert PR.disease_graph(["HP:0011153"], "OMIM:999", engine=eng) is None
