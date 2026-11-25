"""ACMG rule engine v1.4: PM5 ≥1★, recessive PM2 threshold, ClinGen VCEP specs, collagen Gly-X-Y PM1."""
from app.services.acmg_engine import Context, evaluate_variant
from app.services.acmg_engine import cspec
from app.services.acmg_engine import ref_store

from test_acmg_engine import by_code, lof_row, met, missense_row, refs, row


# ── PM5 / PS1 star levels ─────────────────────────────────────────────────
def test_pm5_accepts_one_star_but_ps1_still_needs_two():
    one_star_other = refs(missense_by_residue={("GENE1", 34): [
        {"alt_aa": "Gln", "variant_key": "chr1-1001-G-A", "stars": 1, "category": "P"}]})
    pm5 = met(evaluate_variant(missense_row(), Context(refs=one_star_other)), "PM5")
    assert pm5 and pm5[0]["strength"] == "moderate" and "1★" in pm5[0]["reason"]
    zero_star = refs(missense_by_residue={("GENE1", 34): [
        {"alt_aa": "Gln", "variant_key": "chr1-1001-G-A", "stars": 0, "category": "P"}]})
    assert not met(evaluate_variant(missense_row(), Context(refs=zero_star)), "PM5")
    one_star_same = refs(missense_by_residue={("GENE1", 34): [
        {"alt_aa": "Trp", "variant_key": "chr1-1000-C-A", "stars": 1, "category": "P"}]})
    result = evaluate_variant(missense_row(), Context(refs=one_star_same))
    assert not met(result, "PS1") and not met(result, "PM5")


# ── PM2 for recessive genes ───────────────────────────────────────────────
def test_pm2_recessive_threshold():
    variant = missense_row(GNOMAD_E_AF="0.00015")
    rec = met(evaluate_variant(variant, Context(moi={"AR"})), "PM2")
    assert rec and rec[0]["strength"] == "moderate"                                   # ≤ 0.02%: extremely low (v1.10)
    assert not met(evaluate_variant(variant, Context(moi={"AD", "AR"})), "PM2")   # mixed → 0.0001
    assert not met(evaluate_variant(variant, Context()), "PM2")                     # unknown → 0.0001
    assert not met(evaluate_variant(missense_row(GNOMAD_E_AF="0.0003"), Context(moi={"AR"})), "PM2")
    low = met(evaluate_variant(missense_row(GNOMAD_E_AF="0.00005"), Context()), "PM2")
    assert low and low[0]["strength"] == "supporting"                                 # present but low, not recessive


def test_pm2_absent_is_moderate_but_supporting_next_to_pvs1():
    absent = met(evaluate_variant(missense_row(), Context(moi={"AD"})), "PM2")
    assert absent and absent[0]["strength"] == "moderate" and absent[0]["points"] == 2
    lof = evaluate_variant(lof_row(), Context(refs=refs()))
    pm2 = met(lof, "PM2")[0]
    assert pm2["strength"] == "supporting" and "PVS1" in pm2["reason"]
    assert lof["engine_points"] == 9 and lof["classification"] == "Likely pathogenic"
    capped = met(evaluate_variant(missense_row(), _ctx()), "PM2")                     # VCEP caps to Supporting
    assert capped and capped[0]["strength"] == "supporting"


# ── ClinGen VCEP specs ────────────────────────────────────────────────────
SPEC = {
    "svi": "GN999", "version": "1.0", "vcep": "Test Variant Curation Expert Panel", "url": "https://x",
    "criteria": {"PP2": [], "PP3": ["supporting"], "PVS1": [], "PM2": ["supporting"], "PM1": ["moderate", "strong"],
                 "BS1": ["strong"], "BA1": ["stand_alone"]},
    "thresholds": {"ba1": 0.002, "bs1": 0.0003, "bs1_ar": 0.003, "pm2": 0.00001},
    "pm1": [{"strength": "moderate", "ranges": [[20, 60]]}, {"strength": "strong", "ranges": [[34, 34]]}],
}


def _ctx(spec=SPEC, **kwargs):
    return Context(refs=refs(vcep_specs={"GENE1": spec},
                             missense_genes={"GENE1": (50, 1)}), **kwargs)


def test_vcep_not_applicable_and_strength_cap():
    result = evaluate_variant(missense_row(PKNN_EVIDENCE="PP3_Strong"), _ctx())
    pp3 = met(result, "PP3")
    assert pp3 and pp3[0]["strength"] == "supporting" and "規格最高" in pp3[0]["reason"]
    assert not met(result, "PP2") and "規格不使用 PP2" in by_code(result, "PP2")[0]["reason"]
    assert result["vcep"]["svi"] == "GN999" and any(f["code"] == "vcep" for f in result["flags"])
    lof = evaluate_variant(lof_row(), _ctx())
    assert not met(lof, "PVS1")


def test_vcep_pm1_regions_strongest_wins_and_outside_is_not_met():
    assert met(evaluate_variant(missense_row(), _ctx()), "PM1")[0]["strength"] == "strong"      # p.34
    inside = evaluate_variant(missense_row(HGVS_P="ENSP1:p.Arg50Trp", PROTEIN_POSITION="50/500"), _ctx())
    assert met(inside, "PM1")[0]["strength"] == "moderate"
    outside = evaluate_variant(missense_row(HGVS_P="ENSP1:p.Arg90Trp", PROTEIN_POSITION="90/500"), _ctx())
    assert not met(outside, "PM1") and "不在" in by_code(outside, "PM1")[0]["reason"]
    other = dict(SPEC, pm1=[], pm1_other_gene=True)
    assert not met(evaluate_variant(missense_row(), _ctx(other)), "PM1")


def test_vcep_frequency_thresholds_follow_inheritance():
    af = missense_row(GNOMAD_E_AF="0.0005")
    assert met(evaluate_variant(af, _ctx()), "BS1")                       # 0.0005 > 0.0003
    assert not met(evaluate_variant(af, _ctx(moi={"AR"})), "BS1")         # AR threshold 0.003
    ba1 = evaluate_variant(missense_row(GNOMAD_E_AF="0.003"), _ctx())
    assert met(ba1, "BA1") and "GN999" in met(ba1, "BA1")[0]["reason"]
    rare = missense_row(GNOMAD_E_AF="0.00002")
    assert not met(evaluate_variant(rare, _ctx()), "PM2")                 # VCEP PM2 ≤ 0.00001
    assert met(evaluate_variant(missense_row(GNOMAD_E_AF="0.000005"), _ctx()), "PM2")
    no_pm2 = dict(SPEC, criteria=dict(SPEC["criteria"], PM2=[]))
    assert "規格不使用 PM2" in by_code(evaluate_variant(row(), _ctx(no_pm2)), "PM2")[0]["reason"]


# ── collagen Gly-X-Y ──────────────────────────────────────────────────────
def test_collagen_glycine_gets_pm1_moderate():
    collagen = {"GENE1": {"accession": "P00000", "length": 500, "positions": frozenset({34, 37})}}
    ctx = Context(refs=refs(collagen=collagen))
    gly = missense_row(HGVS_P="ENSP1:p.Gly34Glu")
    hit = met(evaluate_variant(gly, ctx), "PM1")
    assert hit and hit[0]["strength"] == "moderate" and "Gly-X-Y" in hit[0]["reason"]
    assert not met(evaluate_variant(missense_row(HGVS_P="ENSP1:p.Arg34Glu"), ctx), "PM1")
    assert not met(evaluate_variant(missense_row(HGVS_P="ENSP1:p.Gly35Glu", PROTEIN_POSITION="35/500"), ctx), "PM1")
    other_isoform = missense_row(HGVS_P="ENSP1:p.Gly34Glu", PROTEIN_POSITION="34/480")
    assert not met(evaluate_variant(other_isoform, ctx), "PM1")


def test_package_tables_load():
    specs = ref_store.load_vcep_specs()
    assert "MYH7" in specs and specs["MYH7"]["criteria"]["PM2"] == ["supporting"]
    assert "RYR1" not in specs                                     # several rule sets → skipped
    collagen = ref_store.load_collagen()
    assert 816 in collagen["COL4A4"]["positions"] and 1752 in collagen["COL7A1"]["positions"]
    assert collagen["COL4A4"]["length"] == 1690


# ── CSpec parsing ─────────────────────────────────────────────────────────
def test_frequency_and_residue_parsing():
    assert cspec.frequencies("Use if variant is present at ≥0.000083 (0.0083%) in any sub-population") == [8.3e-05]
    assert cspec.frequencies("cutoff of  ≥8 x 10<sup>-3</sup>") == [0.008]
    assert cspec.frequencies("Allele frequency above 0.05%. at least 2,000 alleles") == [0.0005]
    assert cspec.frequencies("MAF cutoff of greater or equal to 0.0333% (or 0.000333).") == [0.000333]
    assert cspec.residues("ATP binding region: aa 19-43; TEY site: aa 169-171") == [[19, 43], [169, 171]]
    assert cspec.residues("residues: Tyr138, Arg158 and R391-S392; p.S1344") == [[138, 138], [158, 158], [391, 392], [1344, 1344]]
    assert cspec.residues("codons using transcript NM_00546.4: 175, 245, 248") == [[175, 175], [245, 245], [248, 248]]
    assert cspec.residues("within exons 3, 4 or 5 (c.676 to c.841)") == []


def _doc(genes, criteria, status="Released", sid="GN900"):
    return {"@id": f"https://cspec/api/SequenceVariantInterpretation/id/{sid}", "cspecStatus": status,
            "version": "1.0", "affiliation": {"label": "Test Variant Curation Expert Panel"},
            "ruleSets": [{"genes": [{"label": g} for g in genes], "criteriaCodes": criteria}]}


def _code(label, strengths):
    return {"label": label, "evidenceStrengths": [
        {"label": name, "applicability": app, "description": text} for name, app, text in strengths]}


def test_build_vcep_specs():
    criteria = [
        _code("PP2", [("Supporting", "Not Applicable for this VCEP", "")]),
        _code("PM2", [("Supporting", "Applicable with VCEP specification",
                       "≤0.00007 for autosomal recessive, ≤0.00002 for autosomal dominant")]),
        _code("BS1", [("Strong", "Applicable", "GnomAD FAF ≥ 0.0001 and < 0.001 (0.01-0.1%)")]),
        _code("BA1", [("Stand Alone", "Applicable as originally described", "Allele frequency ≥ 0.001")]),
        _code("PM1", [("Moderate", "Applicable", "GENEA: residues 10-20. GENEB: aa 30-40")]),
    ]
    specs = cspec.build_vcep_specs([_doc(["GENEA", "GENEB", "GENEC"], criteria),
                                    _doc(["GENEX"], criteria, status="Pilot Rules In Prep", sid="GN901"),
                                    _doc(["GENEA"], criteria, sid="GN902")])
    assert "GENEX" not in specs["genes"] and "GENEA" in specs["skipped"]
    b = specs["genes"]["GENEB"]
    assert b["criteria"]["PP2"] == [] and b["criteria"]["BA1"] == ["stand_alone"]
    assert b["thresholds"]["pm2_ar"] == 7e-05 and b["thresholds"]["pm2_ad"] == 2e-05 and b["thresholds"]["pm2"] == 2e-05
    assert b["thresholds"]["bs1"] == 0.0001 and b["thresholds"]["ba1"] == 0.001
    assert b["pm1"] == [{"strength": "moderate", "ranges": [[30, 40]]}]
    assert specs["genes"]["GENEC"]["pm1_other_gene"] is True


def test_gly_xy_positions():
    seq = "MAAA" + "GPP" * 8 + "MKLM"
    assert cspec.gly_xy_positions(seq) == [5 + 3 * i for i in range(8)]
    assert cspec.gly_xy_positions(seq, helix=[(10, 20)]) == [11, 14, 17, 20]
    assert cspec.gly_xy_positions("MAAA" + "GPP" * 5 + "MKLM") == []     # run too short


# ── PM2 for X-linked recessive genes (v1.5) ───────────────────────────────
def test_pm2_xlinked_recessive_allows_low_frequency():
    rare = missense_row(GNOMAD_E_AF="0.00005")
    xlr = evaluate_variant(rare, Context(moi={"XLR"}))
    assert met(xlr, "PM2")                                                    # female carriers: ≤ 0.0001
    assert by_code(xlr, "PM2")[0]["strength"] == "moderate"                   # ACMG 2015: recessive, extremely low
    assert met(evaluate_variant(rare, Context(moi={"XL"})), "PM2")
    assert not met(evaluate_variant(rare, Context(moi={"XLD"})), "PM2")      # dominant stays 0
    assert not met(evaluate_variant(rare, Context(moi={"XLR", "AD"})), "PM2")
    too_common = evaluate_variant(missense_row(GNOMAD_E_AF="0.0002"), Context(moi={"XLR"}))
    assert not met(too_common, "PM2")
    assert any("X 性聯隱性" in (e.get("reason") or "") for e in by_code(too_common, "PM2"))


# ── PM1 parsing fixes and hand corrections (v1.9) ─────────────────────────
def test_pm1_residues_skip_domain_cys_numbering_coordinates_and_disulfide_pairs():
    assert cspec.residues("critical Gly between Cys2-Cys3 in cbEGF-like domain") == []
    assert cspec.residues("C350 (C350-C382 disulfide in ZP-N domain)") == [[350, 350], [382, 382]]
    assert cspec.residues("substitutions between positions p.873-951; CRX binding site at position "
                          "17:8002593-8002598 (hg38)") == [[873, 951]]
    assert cspec.residues("codons R391-S392 and aa 19-43") == [[19, 43], [391, 392]]


def test_pm1_overrides_apply_only_to_their_version():
    fbn1, note, stale = cspec._apply_pm1_override("FBN1", "1.0", [{"strength": "moderate", "ranges": [[2, 4]]}])
    assert fbn1 == [{"strength": "moderate", "residue_change": "Cys"}] and note and not stale
    kept, note, stale = cspec._apply_pm1_override("FBN1", "2.0", [{"strength": "moderate", "ranges": [[2, 4]]}])
    assert kept == [{"strength": "moderate", "ranges": [[2, 4]]}] and note is None and stale
    pah, _, _ = cspec._apply_pm1_override("PAH", "2.1", [{"strength": "moderate", "ranges": [[245, 245], [268, 268]]}])
    assert pah[0]["ranges"] == [[245, 266], [268, 268], [280, 283], [322, 326], [377, 379]]
    assert cspec._apply_pm1_override("OTOF", "1.0", [])[0] == []


def test_vcep_pm1_cysteine_rule():
    spec = dict(SPEC, pm1=[{"strength": "moderate", "residue_change": "Cys"}])
    lost = evaluate_variant(missense_row(HGVS_P="ENSP1:p.Cys614Arg", PROTEIN_POSITION="614/2871"), _ctx(spec))
    assert met(lost, "PM1")[0]["strength"] == "moderate" and "Cys" in met(lost, "PM1")[0]["reason"]
    gained = evaluate_variant(missense_row(HGVS_P="ENSP1:p.Arg100Cys", PROTEIN_POSITION="100/2871"), _ctx(spec))
    assert met(gained, "PM1")
    other = evaluate_variant(missense_row(HGVS_P="ENSP1:p.Ile614Thr", PROTEIN_POSITION="614/2871"), _ctx(spec))
    assert not met(other, "PM1") and "人工判斷" in by_code(other, "PM1")[0]["reason"]
