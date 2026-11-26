"""ACMG rule engine (plan §5.2) — synthetic rows only."""
import pytest

from app.services import manual_acmg
from app.services.acmg_engine import Context, StaticReferences, evaluate_variant
from app.services.acmg_engine import parse as P
from app.services.acmg_engine.engine import load_rulebook


def row(**overrides):
    base = {
        "CHROM": "chr1", "POS": "1000", "REF": "C", "ALT": "T",
        "GENE": "GENE1", "TRANSCRIPT": "ENST00000000001.1", "TRANSCRIPT_TYPE": "MANE_SELECT",
        "HGVS_C": "ENST00000000001.1:c.100C>T", "HGVS_P": "ENSP00000000001.1:p.Arg34Trp",
        "CONSEQUENCE": "missense_variant", "IMPACT": "MODERATE",
        "EXON": "2/10", "INTRON": ".", "PROTEIN_POSITION": "34/500", "NMD": "",
        "LOFTEE": ".", "LOFTEE_FILTER": ".",
        "GNOMAD_E_AF": ".", "GNOMAD_E_EAS_AF": ".", "GNOMAD_G_AF": ".",
        "GNOMAD41_JOINT_AF": ".", "GNOMAD41_JOINT_EAS_AF": ".",
        "PKNN_EVIDENCE": ".", "PKNN_LLR": ".", "BAYESDEL_NOAF": ".",
        "PANGOLIN_SCORE": "0.01", "PHYLOP100": "1.0",
        "CLINVAR_SIG": ".", "CLINVAR_STARS": ".", "CLINVAR_SIGCONF": ".",
        "GIAB_STRATA": "",
    }
    base.update(overrides)
    return base


def lof_row(**overrides):
    values = dict(CONSEQUENCE="stop_gained", IMPACT="HIGH", HGVS_C="ENST00000000001.1:c.100C>T",
                  HGVS_P="ENSP00000000001.1:p.Arg34Ter", LOFTEE="HC")
    values.update(overrides)
    return row(**values)


def refs(**kwargs):
    kwargs.setdefault("lof_genes", {"GENE1": ("yes", "ClinVar 有 3 個致病截斷變異")})
    return StaticReferences(**kwargs)


def by_code(result, code):
    return [e for e in result["evidence"] if e["code"] == code]


def met(result, code):
    return [e for e in by_code(result, code) if e["status"] == "met"]


# ── parsers ────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("hgvs,offsets,hit", [
    ("ENST1:c.1097+6_1097+11del", [6, 7, 8, 9, 10, 11], False),
    ("ENST1:c.456+1G>C", [1], True),
    ("ENST1:c.1430-4_1430-2delinsTTC", [-4, -3, -2], True),
    ("ENST1:c.100_100+2del", [1, 2], True),
    ("ENST1:c.1096_1097+3del", [1, 2, 3], True),
    ("ENST1:c.100_102del", [0], False),
    ("ENST1:c.147+77_147+121del", list(range(77, 122)), False),
])
def test_hgvs_offsets_and_canonical_hit(hgvs, offsets, hit):
    assert P.hgvs_c_offsets(hgvs) == offsets
    assert P.canonical_splice_hit(hgvs) is hit


def test_hgvs_parsers_handle_missing_and_protein_changes():
    assert P.canonical_splice_hit(".") is None
    assert P.hgvs_p("ENSP1:p.Arg34Trp") == ("Arg", 34, "Trp")
    assert P.hgvs_p("ENSP1:p.Arg34%3D") == ("Arg", 34, "=")
    assert P.deep_intronic("ENST1:c.100+7A>G", 7, -21)
    assert not P.deep_intronic("ENST1:c.100+6A>G", 7, -21)
    assert P.deep_intronic("ENST1:c.101-25A>G", 7, -21)


@pytest.mark.parametrize("value,category", [
    ("Pathogenic", "P"), ("Likely_pathogenic", "LP"), ("Pathogenic/Likely_pathogenic", "LP"),
    ("Uncertain_significance", "VUS"), ("Benign/Likely_benign", "LB"),
    ("Conflicting_classifications_of_pathogenicity", "CONFLICT"),
    ("Pathogenic,_low_penetrance", None), ("risk_factor", None), (".", None),
])
def test_clinvar_category(value, category):
    assert P.clinvar_category(value) == category


# ── PVS1 ───────────────────────────────────────────────────────────────────
def test_recessive_gene_nonsense_gets_pvs1_through_clinvar_mechanism():
    result = evaluate_variant(lof_row(), Context(refs=refs()))
    assert met(result, "PVS1")[0]["strength"] == "very_strong"
    assert result["engine_points"] == 9  # PVS1 8 + PM2_Supporting 1
    assert result["classification"] == "Likely pathogenic"
    assert result["criteria_text"] == "PVS1,PM2_Supporting"


def test_unknown_mechanism_is_not_evaluated_not_failed():
    result = evaluate_variant(lof_row(GENE="GENE9"), Context(refs=refs()))
    pvs1 = by_code(result, "PVS1")[0]
    assert pvs1["status"] == "not_evaluated"
    assert "機轉未確立" in pvs1["reason"]
    assert any(item["code"] == "PVS1" for item in result["not_evaluated"])


def test_missing_reference_table_is_reported_as_not_loaded():
    result = evaluate_variant(lof_row(), Context())
    assert "未載入" in by_code(result, "PVS1")[0]["reason"]


def test_gene_where_lof_is_not_the_mechanism():
    ctx = Context(refs=refs(lof_genes={"GENE1": ("no", "ClinGen HI=40（劑量不敏感）")}))
    assert by_code(evaluate_variant(lof_row(), ctx), "PVS1")[0]["status"] == "not_met"


def test_splice_indel_whose_hgvs_misses_canonical_site_uses_splicing_prediction():
    ctu2 = lof_row(CONSEQUENCE="splice_donor_variant&intron_variant",
                   HGVS_C="ENST00000453996.7:c.1097+6_1097+11del", HGVS_P=".",
                   REF="AGGTGTG", ALT="A", PANGOLIN_SCORE="0.02")
    result = evaluate_variant(ctu2, Context(refs=refs()))
    assert by_code(result, "PVS1")[0]["status"] == "not_met"
    assert "不在 ±1、2" in by_code(result, "PVS1")[0]["reason"]
    assert met(result, "BP4")[0]["reason"].startswith("剪接")


def test_canonical_splice_with_frameshift_exon_skip():
    variant = lof_row(CONSEQUENCE="splice_donor_variant", HGVS_C="ENST00000000001.1:c.456+1G>C",
                      HGVS_P=".", INTRON="3/9")
    in_frame = refs(exons={("ENST00000000001", 3): {"frameshift": False, "fraction": 0.05, "last_exon": False}})
    shift = refs(exons={("ENST00000000001", 3): {"frameshift": True, "fraction": 0.2, "last_exon": False}})
    assert met(evaluate_variant(variant, Context(refs=shift)), "PVS1")[0]["strength"] == "very_strong"
    assert met(evaluate_variant(variant, Context(refs=in_frame)), "PVS1")[0]["strength"] == "moderate"
    assert met(evaluate_variant(variant, Context(refs=refs())), "PVS1")[0]["strength"] == "strong"


def test_nmd_escape_strength_uses_downstream_clinvar_then_fraction():
    escaping = lof_row(NMD="NMD_escaping_variant", PROTEIN_POSITION="480/500")
    critical = refs(truncating_positions={"GENE1": [495]})
    plain = refs(truncating_positions={"GENE1": [10]})
    assert met(evaluate_variant(escaping, Context(refs=critical)), "PVS1")[0]["strength"] == "strong"
    assert met(evaluate_variant(escaping, Context(refs=plain)), "PVS1")[0]["strength"] == "moderate"
    early = lof_row(NMD="NMD_escaping_variant", PROTEIN_POSITION="100/500")
    assert met(evaluate_variant(early, Context(refs=plain)), "PVS1")[0]["strength"] == "strong"


def test_old_tsv_without_nmd_column_uses_exon_position():
    variant = lof_row(EXON="10/10", PROTEIN_POSITION="490/500")
    del variant["NMD"]
    result = evaluate_variant(variant, Context(refs=refs()))
    assert met(result, "PVS1")[0]["strength"] == "moderate"
    assert "缺 NMD 欄位" in met(result, "PVS1")[0]["reason"]


def test_loftee_low_confidence_reasons():
    anc = lof_row(LOFTEE="LC", LOFTEE_FILTER="ANC_ALLELE")
    assert by_code(evaluate_variant(anc, Context(refs=refs())), "PVS1")[0]["status"] == "not_met"
    end = lof_row(LOFTEE="LC", LOFTEE_FILTER="END_TRUNC", PROTEIN_POSITION="495/500")
    assert met(evaluate_variant(end, Context(refs=refs())), "PVS1")[0]["strength"] == "moderate"


def test_lof_only_on_non_reporting_transcript():
    mane = row(CONSEQUENCE="intron_variant", HGVS_P=".")
    other = lof_row(TRANSCRIPT="ENST00000009999.1", TRANSCRIPT_TYPE="BEST_CONSEQUENCE")
    result = evaluate_variant([other, mane], Context(refs=refs()))
    pvs1 = by_code(result, "PVS1")[0]
    assert pvs1["status"] == "not_met" and "ENST00000009999.1" in pvs1["reason"]
    assert result["transcript"] == "ENST00000000001.1"


def test_pvs1_suppresses_pp3_and_pm4():
    variant = lof_row(PANGOLIN_SCORE="0.5")
    result = evaluate_variant(variant, Context(refs=refs()))
    assert met(result, "PVS1")
    assert not met(result, "PP3")


# ── recurrence guard (plan §5.2.5 step 0) ─────────────────────────────────
def test_wes_recurrent_variant_blocks_pvs1_and_pm2():
    variant = lof_row()
    key = P.variant_key(variant)
    ctx = Context(test_type="WES", refs=refs(recurrence={key: 220}, recurrence_patients=250))
    result = evaluate_variant(variant, ctx)
    assert by_code(result, "PVS1")[0]["status"] == "not_evaluated"
    assert by_code(result, "PM2")[0]["status"] == "not_met"
    assert result["flags"][0]["code"] == "recurrent"


def test_wes_recurrence_waits_for_enough_patients_and_respects_clinvar():
    variant = lof_row()
    key = P.variant_key(variant)
    small = Context(test_type="WES", refs=refs(recurrence={key: 40}, recurrence_patients=50))
    result = evaluate_variant(variant, small)
    assert met(result, "PVS1")
    assert any(item["code"] == "RECURRENCE" for item in result["not_evaluated"])
    protected = lof_row(CLINVAR_SIG="Pathogenic", CLINVAR_STARS="2")
    big = Context(test_type="WES", refs=refs(recurrence={key: 220}, recurrence_patients=250))
    assert not evaluate_variant(protected, big)["flags"]


def test_wgs_recurrence_uses_inhouse_af():
    variant = lof_row(INHOUSE_AF="0.03", INHOUSE_AC="84", INHOUSE_AN="2800")
    result = evaluate_variant(variant, Context(test_type="TITAN-WGS", refs=refs()))
    assert result["flags"][0]["code"] == "recurrent"
    assert met(result, "BS1")  # in-house AF also counts as a population frequency


# ── frequency ──────────────────────────────────────────────────────────────
def test_ba1_and_clinvar_exception():
    common = row(GNOMAD41_JOINT_EAS_AF="0.08")
    assert evaluate_variant(common, Context())["classification"] == "Benign"
    known = row(GNOMAD41_JOINT_EAS_AF="0.08", CLINVAR_SIG="Pathogenic", CLINVAR_STARS="2")
    # v1.12: ClinVar P/LP above 5 % protects only on the ClinGen BA1 exception list
    assert met(evaluate_variant(known, Context()), "BA1")
    listed = Context(refs=refs(ba1_exceptions={"chr1-1000-C-T"}))
    result = evaluate_variant(known, listed)
    assert not met(result, "BA1")
    assert any(flag["code"] == "frequency_conflict" for flag in result["flags"])
    below = evaluate_variant(row(GNOMAD41_JOINT_EAS_AF="0.03", CLINVAR_SIG="Pathogenic", CLINVAR_STARS="2"), Context())
    assert not met(below, "BS1") and not met(below, "BA1")


def test_bs1_and_pm2_depend_on_inheritance():
    af = row(GNOMAD_E_AF="0.002")
    assert met(evaluate_variant(af, Context(moi={"AD"})), "BS1")
    recessive = evaluate_variant(af, Context(moi={"AR"}))
    assert not met(recessive, "BS1") and not met(recessive, "PM2")
    rare = row(GNOMAD_E_AF="0.00005")
    assert not met(evaluate_variant(rare, Context(moi={"AD"})), "PM2")
    assert met(evaluate_variant(rare, Context(moi={"AR"})), "PM2")
    assert met(evaluate_variant(row(), Context(moi={"AD"})), "PM2")


def test_small_datasets_and_low_count_inhouse_are_ignored():
    small = row(GNOMAD_G_EAS_AF="0.2", TG_EAS_AF="0.3")
    assert not met(evaluate_variant(small, Context()), "BA1")
    few = row(INHOUSE_AF="0.01", INHOUSE_AC="3", INHOUSE_AN="2800")
    assert evaluate_variant(few, Context())["frequency"]["inhouse_used"] is False


# ── missense / splicing / in-frame ─────────────────────────────────────────
def test_pknn_indeterminate_does_not_fall_through_to_other_tools():
    variant = row(PKNN_LLR="0.3", BAYESDEL_NOAF="0.6")
    result = evaluate_variant(variant, Context())
    assert not met(result, "PP3")
    assert "不換工具" in by_code(result, "PP3")[0]["reason"]


def test_bayesdel_fallback_and_multi_values():
    assert met(evaluate_variant(row(BAYESDEL_NOAF="0.45"), Context()), "PP3")[0]["strength"] == "strong"
    mixed = evaluate_variant(row(BAYESDEL_NOAF="0.45&0.1"), Context())
    assert not met(mixed, "PP3") and not met(mixed, "BP4")


def test_missense_bp4_needs_benign_splicing():
    benign = row(PKNN_EVIDENCE="BP4_Moderate", PANGOLIN_SCORE="0.02")
    assert met(evaluate_variant(benign, Context()), "BP4")[0]["strength"] == "moderate"
    splice = row(PKNN_EVIDENCE="BP4_Moderate", PANGOLIN_SCORE="-0.35")
    result = evaluate_variant(splice, Context())
    assert not met(result, "BP4")
    assert met(result, "PP3")[0]["reason"].startswith("剪接")


def test_pp3_counted_once_and_pm1_cap():
    variant = row(PKNN_EVIDENCE="PP3_Strong", PANGOLIN_SCORE="0.4")
    ctx = Context(refs=refs(vcep_specs={"GENE1": {"pm1_regions": [[30, 40]]}}))
    result = evaluate_variant(variant, ctx)
    assert len(met(result, "PP3")) == 1 and met(result, "PP3")[0]["strength"] == "strong"
    assert not met(result, "PM1")


def test_bp7_takes_precedence_over_splicing_bp4():
    syn = row(CONSEQUENCE="synonymous_variant", HGVS_P="ENSP1:p.Arg34=", PANGOLIN_SCORE="0.01", PHYLOP100="0.5")
    result = evaluate_variant(syn, Context())
    assert met(result, "BP7") and not met(result, "BP4")
    conserved = evaluate_variant(row(CONSEQUENCE="synonymous_variant", PHYLOP100="6.0"), Context())
    assert not met(conserved, "BP7")


def test_inframe_indel_repeat_handling():
    base = dict(CONSEQUENCE="inframe_deletion", REF="CAGG", ALT="C", HGVS_P=".")
    assert met(evaluate_variant(row(**base), Context()), "PM4")
    repeat = evaluate_variant(row(GIAB_STRATA="tandem_repeat", PHYLOP100="0.2", **base), Context())
    assert met(repeat, "BP3") and not met(repeat, "PM4")
    no_giab = row(**base)
    del no_giab["GIAB_STRATA"]
    assert "未確認" in met(evaluate_variant(no_giab, Context()), "PM4")[0]["reason"]


def test_ps1_and_pm5_from_clinvar_residues():
    variant = row()
    key = P.variant_key(variant)
    same = refs(missense_by_residue={("GENE1", 34): [{"alt_aa": "Trp", "variant_key": "chr1-1000-C-A", "stars": 2}]})
    other = refs(missense_by_residue={("GENE1", 34): [{"alt_aa": "Gln", "variant_key": "chr1-1001-G-A", "stars": 3}]})
    self_only = refs(missense_by_residue={("GENE1", 34): [{"alt_aa": "Trp", "variant_key": key, "stars": 4}]})
    weak = refs(missense_by_residue={("GENE1", 34): [{"alt_aa": "Trp", "variant_key": "chr1-1000-C-A", "stars": 1}]})
    assert met(evaluate_variant(variant, Context(refs=same)), "PS1")
    assert met(evaluate_variant(variant, Context(refs=other)), "PM5")
    assert not met(evaluate_variant(variant, Context(refs=self_only)), "PS1")
    assert not met(evaluate_variant(variant, Context(refs=weak)), "PS1")
    assert by_code(evaluate_variant(variant, Context()), "PS1")[0]["status"] == "not_evaluated"


# ── ClinVar rule ───────────────────────────────────────────────────────────
def _mode(mode):
    import copy
    rb = copy.deepcopy(load_rulebook())
    rb["clinvar"]["rank_mode"] = mode
    return rb


def test_clinvar_rank_mode_own_keeps_the_rule_call():
    """Rule v1.13: a ClinVar ≥2★ call no longer moves points or ranking."""
    own = evaluate_variant(lof_row(CLINVAR_SIG="Pathogenic", CLINVAR_STARS="3"), Context(refs=refs()))
    assert own["points"] == own["engine_points"] == 9 and own["classification"] == "Likely pathogenic"
    assert own["clinvar_rule"]["action"] == "none" and own["clinvar_rule"]["eligible"] == "adopt"
    assert own["clinvar_rule"]["mode"] == "own" and own["rank_points"] == 9
    benign = evaluate_variant(lof_row(CLINVAR_SIG="Benign", CLINVAR_STARS="2"), Context(refs=refs()))
    assert benign["points"] == 9 and benign["clinvar_rule"]["eligible"] == "adopt"
    vcep = evaluate_variant(lof_row(CLINGEN_VCEP_CLASS="Benign"), Context(refs=refs()))
    assert vcep["clinvar_rule"]["action"] == "adopt" and vcep["points"] <= -7   # VCEP still adopted


def test_clinvar_adopt_and_floor():
    adopt = _mode("adopt")
    lp_engine = lof_row()  # engine 9 points
    adopt_p = evaluate_variant(lof_row(CLINVAR_SIG="Pathogenic", CLINVAR_STARS="3"), Context(refs=refs()), adopt)
    assert adopt_p["points"] == 10 and adopt_p["classification"] == "Pathogenic"
    assert adopt_p["clinvar_rule"]["action"] == "adopt"
    over = row(PKNN_EVIDENCE="PP3_Strong", CLINVAR_SIG="Likely_pathogenic", CLINVAR_STARS="4")
    over_result = evaluate_variant(over, Context(refs=refs(missense_by_residue={("GENE1", 34): [
        {"alt_aa": "Trp", "variant_key": "chr1-1000-C-A", "stars": 3}]})), adopt)
    assert over_result["engine_points"] == 10 and over_result["points"] == 9   # ClinVar LP caps at 9
    two_star = evaluate_variant(row(CLINVAR_SIG="Likely_pathogenic", CLINVAR_STARS="2"), Context(), adopt)
    assert two_star["points"] == 6 and two_star["clinvar_rule"]["action"] == "adopt"   # 2★ adopted (v1.1)
    assert two_star["classification"] == "Likely pathogenic"
    one_star = evaluate_variant(row(CLINVAR_SIG="Pathogenic", CLINVAR_STARS="1"), Context(), adopt)
    assert one_star["clinvar_rule"]["action"] == "none" and one_star["clinvar_rule"]["eligible"] is None
    conflict = evaluate_variant(row(CLINVAR_SIG="Likely_pathogenic", CLINVAR_STARS="2",
                                    CLINVAR_SIGCONF="Pathogenic(1)|Uncertain_significance(1)"), Context(), adopt)
    assert conflict["clinvar_rule"]["action"] == "none" and conflict["clinvar_rule"]["conflict"]
    assert evaluate_variant(lp_engine, Context(refs=refs()), adopt)["clinvar_rule"]["action"] == "none"
    # floor (validation option B): P/LP ≥2★ raises to 6 one-way, never lowers
    floor = _mode("floor")
    low = evaluate_variant(row(CLINVAR_SIG="Pathogenic", CLINVAR_STARS="2"), Context(), floor)
    assert low["points"] == 6 and low["rank_points"] == 6 and low["clinvar_rule"]["action"] == "floor"
    high = evaluate_variant(over, Context(refs=refs(missense_by_residue={("GENE1", 34): [
        {"alt_aa": "Trp", "variant_key": "chr1-1000-C-A", "stars": 3}]})), floor)
    assert high["points"] == 10 and high["clinvar_rule"]["action"] == "none"
    benign = evaluate_variant(lof_row(CLINVAR_SIG="Benign", CLINVAR_STARS="3"), Context(refs=refs()), floor)
    assert benign["points"] == 9


def test_clinvar_adjust_reranks_stored_results():
    from app.services.acmg_engine.engine import clinvar_adjust
    rb = load_rulebook()
    stored = evaluate_variant(row(CLINVAR_SIG="Likely_pathogenic", CLINVAR_STARS="2"), Context())
    rule, rank = stored["clinvar_rule"], stored["rank_points"]
    assert clinvar_adjust(rank, rule, "own", rb) == rank
    assert clinvar_adjust(rank, rule, "floor", rb) == max(rank, 6)
    assert clinvar_adjust(12.0, rule, "adopt", rb) == 9 and clinvar_adjust(12.0, rule, "floor", rb) == 12.0
    old = {"action": "adopt", "source": "ClinVar", "category": "P", "stars": 2}      # stored before v1.13
    assert clinvar_adjust(3.0, old, "adopt", rb) == 10
    vcep = {"action": "adopt", "source": "ClinGen VCEP", "category": "B", "stars": 0}
    assert clinvar_adjust(-8.0, vcep, "floor", rb) == -8.0


def test_latest_clinvar_preferred_and_vcep_adopted():
    latest = row(CLINVAR_SIG="Uncertain_significance", CLINVAR_STARS="1",
                 CLINVAR_LATEST_SIG="Pathogenic", CLINVAR_LATEST_STARS="2")
    assert evaluate_variant(latest, Context())["clinvar_rule"]["source"] == "ClinVar latest"
    vcep = evaluate_variant(lof_row(CLINGEN_VCEP_CLASS="Benign"), Context(refs=refs()))
    assert vcep["classification"] == "Benign" and vcep["points"] <= -7
    override = evaluate_variant(row(CLINVAR_SIG="Pathogenic", CLINVAR_STARS="3"),
                                Context(clinvar={"category": None, "stars": 0, "conflict": False}))
    assert override["clinvar_rule"]["action"] == "none"


# ── compatibility with the manual ACMG editor ─────────────────────────────
def test_criteria_text_round_trips_through_manual_acmg():
    ctx = Context(refs=refs(truncating_positions={"GENE1": [495]}))
    result = evaluate_variant(lof_row(NMD="NMD_escaping_variant", PROTEIN_POSITION="480/500"), ctx)
    parsed, unknown = manual_acmg.parse_criteria_text(result["criteria_text"])
    assert not unknown
    assert parsed == result["criteria"]
    assert manual_acmg.calculate(parsed)["score"] == result["engine_points"]


# ── v1.1: PM5 grade, ClinVar-derived PM1 / PP2, BP7 without phyloP ────────
def missense_row(**values):
    base = dict(CONSEQUENCE="missense_variant", IMPACT="MODERATE", HGVS_P="ENSP1:p.Arg34Trp", LOFTEE=".")
    base.update(values)
    return row(**base)


def test_pm5_strength_follows_reference_grade():
    lp_only = refs(missense_by_residue={("GENE1", 34): [
        {"alt_aa": "Gln", "variant_key": "chr1-1001-G-A", "stars": 2, "category": "LP"}]})
    with_p = refs(missense_by_residue={("GENE1", 34): [
        {"alt_aa": "Gln", "variant_key": "chr1-1001-G-A", "stars": 2, "category": "LP"},
        {"alt_aa": "Gly", "variant_key": "chr1-1001-G-C", "stars": 2, "category": "P"}]})
    assert met(evaluate_variant(missense_row(), Context(refs=lp_only)), "PM5")[0]["strength"] == "supporting"
    assert met(evaluate_variant(missense_row(), Context(refs=with_p)), "PM5")[0]["strength"] == "moderate"


def _hotspot(positions, benign=()):
    residues = {("GENE1", p): [{"alt_aa": "Gln", "variant_key": f"chr1-{p}-G-A", "stars": 2, "category": "P"}]
                for p in positions}
    return refs(missense_by_residue=residues, benign_missense={"GENE1": [
        {"protein_pos": p, "variant_key": f"chr1-{p}-A-G", "stars": 2} for p in benign]})


def test_pm1_from_clinvar_hotspot():
    variant = missense_row()                       # residue 34; own residue never counts
    assert not met(evaluate_variant(variant, Context(refs=_hotspot([30, 34]))), "PM1")
    supporting = met(evaluate_variant(variant, Context(refs=_hotspot([30, 32, 40]))), "PM1")
    assert supporting and supporting[0]["strength"] == "supporting"
    moderate = met(evaluate_variant(variant, Context(refs=_hotspot([25, 28, 30, 32, 40, 44]))), "PM1")
    assert moderate and moderate[0]["strength"] == "moderate"
    blocked = evaluate_variant(variant, Context(refs=_hotspot([30, 32, 40], benign=[36])))
    assert not met(blocked, "PM1") and by_code(blocked, "PM1")[0]["status"] == "not_met"
    far = evaluate_variant(variant, Context(refs=_hotspot([1, 2, 3])))
    assert not met(far, "PM1")
    vcep_outside = refs(missense_by_residue=_hotspot([30, 32, 40]).missense_by_residue, benign_missense={},
                        vcep_specs={"GENE1": {"pm1_regions": [[100, 120]]}})
    assert not met(evaluate_variant(variant, Context(refs=vcep_outside)), "PM1")


def test_pp2_from_gene_missense_counts():
    yes = refs(missense_genes={"GENE1": (20, 1)})
    too_benign = refs(missense_genes={"GENE1": (20, 7)})     # 26% benign > 20%
    edge = refs(missense_genes={"GENE1": (20, 5)})           # exactly 20% (v1.8)
    too_few = refs(missense_genes={"GENE1": (5, 0)})
    assert met(evaluate_variant(missense_row(), Context(refs=yes)), "PP2")
    assert not met(evaluate_variant(missense_row(), Context(refs=too_benign)), "PP2")
    assert met(evaluate_variant(missense_row(), Context(refs=edge)), "PP2")
    assert not met(evaluate_variant(missense_row(), Context(refs=too_few)), "PP2")
    assert not met(evaluate_variant(lof_row(), Context(refs=yes)), "PP2")      # not missense
    assert "PP2" in evaluate_variant(missense_row(), Context(refs=yes))["criteria_text"]


def test_bp7_without_phylop_uses_splice_prediction():
    syn = row(CONSEQUENCE="synonymous_variant", HGVS_P="ENSP1:p.Arg34=", PANGOLIN_SCORE="0.01", PHYLOP100=".")
    bp7 = met(evaluate_variant(syn, Context()), "BP7")
    assert bp7 and "缺 phyloP" in bp7[0]["reason"]
    no_splice = row(CONSEQUENCE="synonymous_variant", HGVS_P="ENSP1:p.Arg34=", PANGOLIN_SCORE=".", PHYLOP100=".")
    assert by_code(evaluate_variant(no_splice, Context()), "BP7")[0]["status"] == "not_evaluated"
