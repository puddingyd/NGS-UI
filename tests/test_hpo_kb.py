"""Phase 1: phenotype knowledge base (plan §2.3, §4.2, §4.4, §4.5.4) — synthetic data."""
import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pytest

from app.services import hpo_gene_table as H
from app.services.pheno_engine import disease_gene as DG, kb as K, mondo as M, npzio

from test_hpo_gene_table import G2D, OBO, ROWS, _hpoa

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def src(tmp_path):
    (tmp_path / "hp.obo").write_text(OBO, encoding="utf-8")
    (tmp_path / "phenotype.hpoa").write_text(_hpoa(ROWS), encoding="utf-8")
    (tmp_path / "genes_to_disease.txt").write_text(G2D, encoding="utf-8")
    return tmp_path


def _kb(src):
    result = H.build(src / "hp.obo", src / "phenotype.hpoa", src / "genes_to_disease.txt")
    return result, K.build_kb(result)


def test_ic_kernel_and_baseline(src):
    result, kb = _kb(src)
    tid = {t: i for i, t in enumerate(kb.terms)}
    assert kb.diseases == ["OMIM:1", "OMIM:2", "OMIM:4", "ORPHA:3"] and kb.N == 4
    ln43 = math.log(4 / 3)
    assert kb.ic[tid["HP:0000118"]] == 0.0 and math.copysign(1, kb.ic[tid["HP:0000118"]]) == 1
    assert kb.ic[tid["HP:0001250"]] == pytest.approx(ln43) and kb.ic[tid["HP:0001249"]] == pytest.approx(math.log(4))
    assert "HP:0003333" not in tid                          # broken chain is not a phenotype term
    k = K.Kernel(kb)
    seizure = k.evec(tid["HP:0001250"])
    assert seizure[tid["HP:0000707"]] == pytest.approx(ln43)          # disease term is an ancestor
    assert seizure[tid["HP:0001249"]] == pytest.approx(0.5 * ln43)    # sibling, lift 1 → κ 0.5
    assert seizure[tid["HP:0000478"]] == 0.0                          # only the root in common
    assert k.evec(tid["HP:0000707"])[tid["HP:0001250"]] == pytest.approx(ln43)   # patient broader
    mu, pool = K.baseline(kb, k)
    assert pool == 3
    assert mu[kb.diseases.index("OMIM:1")] == pytest.approx(0.5 * ln43)
    brute = np.mean([k.best(int(q)) for q in np.unique(kb.dir_idx)], axis=0)
    assert np.allclose(mu, brute)
    qc = K.qc(kb, result["ontology"], mu)
    assert qc["passed"] and qc["N"] == 4 and qc["terms_p0"] == 0


def test_unannotated_patient_term_uses_ic_star(src):
    obo = OBO.replace("[Typedef]", "[Term]\nid: HP:0009001\nname: Focal thing\nis_a: HP:0001250\n\n[Typedef]")
    (src / "hp.obo").write_text(obo, encoding="utf-8")
    result, kb = _kb(src)
    tid = {t: i for i, t in enumerate(kb.terms)}
    q = tid["HP:0009001"]
    assert kb.n[q] == 0 and kb.ic_star[q] == pytest.approx(kb.ic[tid["HP:0001250"]])
    e = K.Kernel(kb).evec(q)
    assert e[tid["HP:0001250"]] == pytest.approx(kb.ic[tid["HP:0001250"]])   # disease term is an ancestor
    assert np.isfinite(e).all()


def test_omim_and_orphanet_links(tmp_path):
    rows = [{"gene_symbol": "GENEA", "OMIM_disease": "Disease A (100001)(AD)\n?Maybe B (100002)(AR)\n"
                                                    "{Risk C} (100003)\n[Blood group D] (100004)\n"
                                                    "Tumour E, somatic (100005)(SMu)"}]
    links = {r["disease_id"]: r for r in DG.omim_links(rows, canonical=lambda g: g)}
    assert links["OMIM:100001"]["use"] == 1 and links["OMIM:100001"]["flags"] == ""
    assert links["OMIM:100002"]["flags"] == "provisional" and links["OMIM:100002"]["use"] == 1
    assert links["OMIM:100003"]["flags"] == "susceptibility"
    assert links["OMIM:100004"]["use"] == 0 and links["OMIM:100005"]["use"] == 0
    xml = tmp_path / "p6.xml"
    xml.write_text("""<?xml version="1.0"?><JDBOR date="2026-07-01" version="1.3.36"><DisorderList>
<Disorder><OrphaCode>3</OrphaCode><DisorderGeneAssociationList>
 <DisorderGeneAssociation><Gene><Symbol>GENEC</Symbol></Gene>
  <DisorderGeneAssociationType><Name lang="en">Disease-causing germline mutation(s) (loss of function) in</Name></DisorderGeneAssociationType>
  <DisorderGeneAssociationStatus><Name lang="en">Assessed</Name></DisorderGeneAssociationStatus></DisorderGeneAssociation>
 <DisorderGeneAssociation><Gene><Symbol>GENES</Symbol></Gene>
  <DisorderGeneAssociationType><Name lang="en">Disease-causing somatic mutation(s) in</Name></DisorderGeneAssociationType></DisorderGeneAssociation>
 <DisorderGeneAssociation><Gene><Symbol>GENER</Symbol></Gene>
  <DisorderGeneAssociationType><Name lang="en">Major susceptibility factor in</Name></DisorderGeneAssociationType></DisorderGeneAssociation>
</DisorderGeneAssociationList></Disorder></DisorderList></JDBOR>""", encoding="utf-8")
    orpha, info = DG.orphanet_links(xml, canonical=lambda g: g)
    by_gene = {r["gene_symbol"]: r for r in orpha}
    assert info["version"] == "1.3.36" and by_gene["GENEC"]["use"] == 1 and by_gene["GENEC"]["disease_id"] == "ORPHA:3"
    assert by_gene["GENES"]["flags"] == "excluded_type" and by_gene["GENER"]["flags"] == "review"
    rows, stats = DG.cross_check(orpha + DG.omim_links(rows, canonical=lambda g: g),
                                 {"ORPHA:3": {"GENEC": "3", "GENEX": "9"}, "OMIM:100001": {"GENEA": "1"}},
                                 {"ORPHA:3", "OMIM:100001"}, canonical=lambda g: g)
    assert stats["ORPHA"]["both"] == 1 and ("ORPHA:3", "GENEX", 0, 1) in rows and stats["OMIM"]["both"] == 1


def test_mondo_groups(tmp_path):
    obo = tmp_path / "mondo.obo"
    obo.write_text("""data-version: releases/2026-11-03

[Term]
id: MONDO:0000001
name: thing one
xref: OMIM:1 {source="MONDO:equivalentTo"}
xref: Orphanet:3 {source="MONDO:equivalentTo"}
xref: OMIM:2 {source="MONDO:relatedTo"}

[Term]
id: MONDO:0000002
name: gone
is_obsolete: true
xref: OMIM:4 {source="MONDO:equivalentTo"}
""", encoding="utf-8")
    groups, labels, info = M.parse(obo)
    assert groups == {"MONDO:0000001": ["OMIM:1", "ORPHA:3"]} and info["mode"] == "obo-equivalentTo"
    rows, stats = M.disease_groups(groups, {"OMIM:1", "OMIM:2", "ORPHA:3"})
    assert rows == [("MONDO:0000001", "OMIM:1"), ("MONDO:0000001", "ORPHA:3")] and stats["groups_omim_and_orpha"] == 1


def test_npz_is_byte_reproducible(tmp_path):
    a, b = tmp_path / "a.npz", tmp_path / "b.npz"
    npzio.save(a, x=np.arange(5), s=np.array(["HP:1", "HP:2"]))
    npzio.save(b, s=np.array(["HP:1", "HP:2"]), x=np.arange(5))
    assert a.read_bytes() == b.read_bytes() and list(npzio.load(a)["s"]) == ["HP:1", "HP:2"]


def _build_script():
    spec = importlib.util.spec_from_file_location("build_hpo_kb", ROOT / "scripts/build_hpo_kb.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_build_script_end_to_end_and_reproducible(src, tmp_path, monkeypatch):
    module = _build_script()
    rows = [{"gene_symbol": "GENEA", "OMIM_disease": "D1 (000001)(AD)"}]
    monkeypatch.setattr(module, "_platform_omim_rows", lambda: (rows, {"source": "test"}))
    (src / "mondo.obo").write_text("[Term]\nid: MONDO:1\nname: m\nxref: OMIM:1 {source=\"MONDO:equivalentTo\"}\n",
                                   encoding="utf-8")
    args = ["--hp-obo", str(src / "hp.obo"), "--hpoa", str(src / "phenotype.hpoa"), "--genes-to-disease",
            str(src / "genes_to_disease.txt"), "--mondo", str(src / "mondo.obo"), "--allow-missing", "orphanet"]
    assert module.main(args + ["--kb-root", str(tmp_path / "r1")]) == 0
    assert module.main(args + ["--kb-root", str(tmp_path / "r2")]) == 0
    s1, s2 = (next((tmp_path / r).glob("2026-*")) for r in ("r1", "r2"))
    m1, m2 = (json.loads((s / "manifest.json").read_text(encoding="utf-8")) for s in (s1, s2))
    assert m1["kb"]["content_sha256"] == m2["kb"]["content_sha256"]           # same versions → identical files
    assert m1["kb"]["qc"]["passed"] and "hpo_gene.tsv" in m1["kb"]["content_sha256"]
    for name in ("terms.tsv", "ancestors.npz", "disease_closure.npz", "disease_baseline.tsv", "disease_gene.tsv",
                 "disease_gene_discordant.tsv", "disease_groups.tsv"):
        assert (s1 / name).is_file()
    assert "## 知識庫" in (s1 / "build_qc.md").read_text(encoding="utf-8")
    kb, mu = K.load_kb(s1)
    assert kb.N == 4 and len(mu) == 4
    # missing sources refuse unless allowed
    with pytest.raises(SystemExit):
        module.main(args[:-2] + ["--kb-root", str(tmp_path / "r3")])
    assert not list((tmp_path / "r3").glob("2026-*"))                        # nothing half-written
