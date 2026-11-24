"""Phase 2: phenotype engine (plan §4.3–4.8) on a small synthetic knowledge base."""
import csv
import json
import math

import numpy as np
import pytest

from app.services import hpo_update as U
from app.services.pheno_engine import engine as E

OBO = """format-version: 1.2
data-version: hp/releases/2026-09-01

[Term]
id: HP:0000001
name: All

[Term]
id: HP:0000118
name: Phenotypic abnormality
is_a: HP:0000001

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
id: HP:0007359
name: Focal-onset seizure
is_a: HP:0001250

[Term]
id: HP:0002197
name: Generalized-onset seizure
is_a: HP:0001250

[Term]
id: HP:0011153
name: Focal motor seizure
is_a: HP:0007359

[Term]
id: HP:0001249
name: Intellectual disability
is_a: HP:0000707

[Term]
id: HP:0000478
name: Abnormality of the eye
is_a: HP:0000118

[Term]
id: HP:0000316
name: Hypertelorism
is_a: HP:0000478

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
"""

HEADER = ["database_id", "disease_name", "qualifier", "hpo_id", "reference", "evidence", "onset",
          "frequency", "sex", "modifier", "aspect", "biocuration"]
ANN = {
    "OMIM:100001": ("Focal epilepsy A", ["HP:0007359", "HP:0001249"]),
    "OMIM:100002": ("Generalized epilepsy B", ["HP:0002197"]),
    "OMIM:100003": ("Broad syndrome C", ["HP:0001250", "HP:0001249", "HP:0000316"]),
    "OMIM:100004": ("Eye disease D", ["HP:0000316"]),
    "ORPHA:5": ("Focal epilepsy A", ["HP:0007359"]),
}
OMIM_ROWS = [
    {"gene_symbol": "GENEA", "OMIM_disease": "Focal epilepsy A (100001)(AD)\nBroad syndrome C (100003)(AR)"},
    {"gene_symbol": "GENEB", "OMIM_disease": "Generalized epilepsy B (100002)(AD)"},
    {"gene_symbol": "GENED", "OMIM_disease": "Eye disease D (100004)(AR)"},
    {"gene_symbol": "GENEZ", "OMIM_disease": "Not annotated (100099)(AD)"},
]
MONDO = """data-version: releases/2026-11-03

[Term]
id: MONDO:0000001
name: focal epilepsy a
xref: OMIM:100001 {source="MONDO:equivalentTo"}
xref: Orphanet:5 {source="MONDO:equivalentTo"}
"""


@pytest.fixture(scope="module")
def snap(tmp_path_factory):
    src = tmp_path_factory.mktemp("src")
    lines = ["#description: test", "#version: 2026-09-01", "\t".join(HEADER)]
    for disease, (name, terms) in ANN.items():
        for t in terms:
            row = dict.fromkeys(HEADER, "")
            row.update(database_id=disease, disease_name=name, hpo_id=t, aspect="P", reference=disease,
                       evidence="TAS", frequency="3/4" if t == "HP:0007359" else "")
            lines.append("\t".join(row[h] for h in HEADER))
    (src / "phenotype.hpoa").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (src / "hp.obo").write_text(OBO, encoding="utf-8")
    (src / "genes_to_disease.txt").write_text("ncbi_gene_id\tgene_symbol\tassociation_type\tdisease_id\tsource\n",
                                              encoding="utf-8")
    (src / "mondo.obo").write_text(MONDO, encoding="utf-8")
    folder, _, _ = U.build_snapshot(src / "hp.obo", src / "phenotype.hpoa", src / "genes_to_disease.txt",
                                    kb_root=src / "kb", mondo=src / "mondo.obo", allow_missing={"orphanet"},
                                    omim=(OMIM_ROWS, {"source": "test"}), progress=lambda *a: None)
    return folder


@pytest.fixture(scope="module")
def eng(snap):
    return E.Engine(snap)


def _ref_e(eng, q, d):
    """e(q, d) straight from the plan's definition (independent of the vectorised kernel)."""
    kb = eng.kb
    anc = [set(kb.anc(i).tolist()) for i in range(kb.T)]
    closure = [set(kb.cl_idx[kb.cl_ptr[i]:kb.cl_ptr[i + 1]].tolist()) for i in range(kb.N)]
    ic = lambda t: float(kb.ic[t]) if kb.n[t] > 0 else 0.0                     # noqa: E731
    if q in anc[d]:
        return float(kb.ic_star[q])
    if d in anc[q]:
        return ic(d)
    common = [a for a in anc[q] & anc[d] if kb.n[a] > 0]
    m = max(common, key=lambda a: (ic(a), -a))
    qa = q if kb.n[q] > 0 else max((a for a in anc[q] if kb.n[a] > 0), key=lambda a: (ic(a), -a))
    n_qd = sum(1 for c in closure if qa in c and d in c)
    lift = (n_qd + .5) * (kb.n[m] + .5) / ((kb.n[qa] + .5) * (kb.n[d] + .5))
    return 0.5 * min(1.0, lift) * ic(m)


def test_preprocess_normalises_and_removes(eng):
    pre = eng.preprocess([
        {"phenotype": "HP:0009999", "weight": 1},       # alt_id of Seizure, also ancestor of focal → redundant
        {"phenotype": "HP:0007359", "weight": 1},
        {"phenotype": "HP:0001111", "weight": 3},       # obsolete → Seizure (duplicate), weight 3 carried down
        {"phenotype": "HP:0002222"},                    # obsolete, no replacement
        {"phenotype": "HP:0000006"},                    # inheritance: not a phenotype term
        {"phenotype": "HP:0000118"},                    # root
        {"phenotype": "HP:1234567"},                    # unknown
        {"phenotype": "HP:0000316", "weight": "x"},     # bad weight → 1
        {"name": "some panel", "weight": 1},            # not an HPO entry: ignored
    ])
    kept = {t["hpo_id"]: t for t in pre["terms"]}
    assert set(kept) == {"HP:0007359", "HP:0000316"}
    assert kept["HP:0007359"]["weight"] == 3 and kept["HP:0007359"]["input_weight"] == 1
    assert kept["HP:0000316"]["weight"] == 1
    reasons = {(r["hpo_id"], r["reason"]) for r in pre["removed"]}
    assert ("HP:0002222", "obsolete") in reasons and ("HP:0000006", "non_phenotype") in reasons
    assert ("HP:0000118", "root") in reasons and ("HP:1234567", "unknown") in reasons
    assert ("HP:0001111", "duplicate") in reasons
    assert any(r["reason"] == "redundant_ancestor" and r["current_id"] == "HP:0001250" and r["kept"] == "HP:0007359"
               for r in pre["removed"])
    siblings = eng.preprocess(["HP:0007359", "HP:0002197"])                    # siblings both kept
    assert [t["hpo_id"] for t in siblings["terms"]] == ["HP:0002197", "HP:0007359"]


def test_scores_match_definition(eng):
    kb = eng.kb
    tid = eng.tid
    patient = ["HP:0011153", "HP:0001249", "HP:0000316"]       # focal motor (no annotation → IC*), ID, hypertelorism
    res = eng.score([(t, 1) for t in patient])
    for k, t in enumerate(patient_sorted := sorted(patient, key=lambda x: tid[x])):
        q = tid[t]
        for D in range(kb.N):
            ann = kb.dir_idx[kb.dir_ptr[D]:kb.dir_ptr[D + 1]]
            ref = max(_ref_e(eng, q, int(d)) for d in ann)
            assert res.s[k, D] == pytest.approx(ref)
    assert len(patient_sorted) == len(res.q)
    assert np.allclose(res.S_star, res.S - eng.mu * 3)
    assert (res.X >= 0).all() and (res.X <= 1 + 1e-12).all()
    # focal motor seizure is unannotated: its IC* is that of focal-onset seizure
    assert kb.ic_star[tid["HP:0011153"]] == pytest.approx(kb.ic[tid["HP:0007359"]])


def test_pair_detail_and_relations(eng):
    tid = eng.tid
    focal, gen, seiz = tid["HP:0007359"], tid["HP:0002197"], tid["HP:0001250"]
    assert eng.pair(focal, focal)["relation"] == "same"
    assert eng.pair(seiz, focal)["relation"] == "patient_broader"
    p = eng.pair(focal, seiz)
    assert p["relation"] == "patient_specific" and p["e"] == pytest.approx(float(eng.kb.ic[seiz]))
    c = eng.pair(focal, gen)
    assert c["relation"] == "collateral" and 0 < c["kappa"] <= 0.5
    # Seizure and "nervous system" reach the same diseases here (equal IC): the tie goes to the smaller ID
    assert c["ic_mica"] == pytest.approx(float(eng.kb.ic[seiz])) and eng.kb.terms[c["mica"]] == "HP:0000707"
    assert c["e"] == pytest.approx(eng.kernel.evec(focal)[gen])
    none = eng.pair(focal, tid["HP:0000316"])
    assert none["relation"] == "none" and none["e"] == 0.0


def test_rankings_points_and_explanation(eng):
    res = eng.score(["HP:0007359", "HP:0001249"])
    top = res.top_diseases(5)
    assert top[0]["disease_id"] == "OMIM:100001"                  # exact matches on both terms
    assert [d["rank"] for d in top] == list(range(1, len(top) + 1))
    genes = res.top_genes(5)
    assert genes[0]["gene"] == "GENEA" and genes[0]["best_disease"] == "OMIM:100001"
    assert genes[0]["P"] == 10.3 and genes[0]["rank"] == 1
    assert {g["gene"] for g in genes} == {"GENEA", "GENEB", "GENED"}
    assert res.gene("GENEZ")["status"] == "na" and res.gene("NOPE") is None
    lead = res.lead_nats()
    assert genes[0]["lead"] == (lead >= E.LEAD_NATS)
    # P(g, d): best disease keeps P(g); the gene's other disease is reduced but not below −3.1
    pts = {x["disease_id"]: x for x in res.gene_disease_points("GENEA")}
    assert pts["OMIM:100001"]["P"] == 10.3 and pts["OMIM:100001"]["best"]
    assert E.P_MIN <= pts["OMIM:100003"]["P"] <= 10.3
    # OMIM:100001 and ORPHA:5 are one MONDO group: same group rank
    i1, i5 = eng.did["OMIM:100001"], eng.did["ORPHA:5"]
    assert res.d_group_rank[i1] == res.d_group_rank[i5] == 1.0
    rows = res.explain(i1)
    assert [r["relation"] for r in rows] == ["same", "same"]
    assert sum(r["net"] for r in rows) == pytest.approx(res.S_star[i1])


def test_points_tables():
    assert E.gene_points(1) == 10.3 and E.gene_points(1.5) == 7.5 and E.gene_points(10) == 5.7
    assert E.gene_points(11) == 4.2 and E.gene_points(5000) == -3.1
    assert E.disease_points(3) == 8.9 and E.disease_points(2000) == -3.1


def test_weights_scale_contributions(eng):
    a = eng.score([("HP:0007359", 1), ("HP:0000316", 1)])
    b = eng.score([("HP:0007359", 2), ("HP:0000316", 1)])
    tid = eng.tid
    ka = [eng.kb.terms[q] for q in a.q].index("HP:0007359")
    assert np.allclose(b.S - a.S, a.s[ka])
    assert np.allclose(b.S_star, b.S - eng.mu * 3)
    assert tid  # silence


def test_empty_input(eng, tmp_path):
    res = eng.score(["HP:0000006"])
    assert res.empty and res.top_genes(3) == [] and res.gene("GENEA")["status"] == "no_terms"
    paths = E.write_outputs(res, tmp_path)
    assert [p.name for p in paths] == list(E.FILES)
    run = json.loads((tmp_path / "pheno_engine_run.json").read_text(encoding="utf-8"))
    assert run["removed"][0]["reason"] == "non_phenotype" and run["terms"] == []


def test_outputs_are_reproducible(eng, snap, tmp_path):
    items = [("HP:0011153", 2), ("HP:0001249", 1), ("HP:0000316", 1)]
    a, b = tmp_path / "a", tmp_path / "b"
    E.write_outputs(eng.score(items), a, prefix="S1.")
    E.write_outputs(E.Engine(snap).score(items), b, prefix="S1.")
    for name in E.FILES:
        assert (a / f"S1.{name}").read_bytes() == (b / f"S1.{name}").read_bytes()
    with open(a / "S1.pheno_engine_genes.tsv", encoding="utf-8") as handle:
        genes = list(csv.DictReader(handle, delimiter="\t"))
    assert genes[-1]["rank"] == "NA" and genes[-1]["gene"] == "GENEZ"
    with open(a / "S1.pheno_engine_explain.tsv", encoding="utf-8") as handle:
        explain = list(csv.DictReader(handle, delimiter="\t"))
    focal = [r for r in explain if r["disease_id"] == "OMIM:100001" and r["hpo_id"] == "HP:0011153"][0]
    assert focal["relation"] == "patient_specific" and focal["annotation_frequency"] == "3/4"
    run = json.loads((a / "S1.pheno_engine_run.json").read_text(encoding="utf-8"))
    assert run["engine_version"] == E.ENGINE_VERSION and run["snapshot"] == snap.name
    assert {t["hpo_id"]: t["weight"] for t in run["terms"]}["HP:0011153"] == 2 and math.isfinite(run["counts"]["information_total"])


def test_load_engine_caches(snap):
    a = E.load_engine(snap)
    assert a is not None and E.load_engine(snap) is a
    assert E.load_engine(snap.parent / "missing") is None


def test_run_for_analysis_writes_prefixed_files(eng, tmp_path, monkeypatch):
    from app.services import analyses_store, sample_layout
    monkeypatch.setattr(analyses_store, "active_version", lambda sid: "default")
    monkeypatch.setattr(analyses_store, "read_version",
                        lambda sid, ver: {"hpo": [{"phenotype": "HP:0007359", "weight": 1}], "selected_panels": []})
    monkeypatch.setattr(sample_layout, "analysis_file",
                        lambda sid, ver, name, for_write=False: tmp_path / ver / f"{sid}.{name}")
    out = E.run_for_analysis("S9", engine=eng)
    assert out["version"] == "default"
    assert sorted(p.name for p in (tmp_path / "default").iterdir()) == sorted(f"S9.{n}" for n in E.FILES)
    run = json.loads((tmp_path / "default" / "S9.pheno_engine_run.json").read_text(encoding="utf-8"))
    assert run["sample_id"] == "S9" and run["analysis_version"] == "default"
