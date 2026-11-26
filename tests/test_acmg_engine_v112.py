"""ACMG rule engine v1.12: NMD from MANE exons, OMIM disease-line inheritance, ClinVar above BA1."""
from app.services import inhouse_common
from app.services.acmg_engine import Context, evaluate_variant
from app.services.acmg_engine import parse as P
from app.services.acmg_engine import runner

from test_acmg_engine import by_code, lof_row, met, refs


# STRC-like MANE transcript: 29 exons, CDS 5,325 nt (1,775 aa), last exon 86 coding nt
STRC_EXONS = [180] * 27 + [148, 86]
_JUNCTION = sum(STRC_EXONS[:-1])          # 5,008: last coding base before the last junction


def _exon_refs(lengths=STRC_EXONS, **kwargs):
    exons = {("ENST00000000001", i + 1): {"cds_length": n, "frameshift": n % 3 != 0,
                                          "fraction": n / sum(lengths), "last_exon": i == len(lengths) - 1}
             for i, n in enumerate(lengths)}
    return refs(exons=exons, **kwargs)


def _stop(codon, **values):
    base = dict(HGVS_P=f"ENSP00000000001.1:p.Arg{codon}Ter", PROTEIN_POSITION=f"{codon}/1775",
                EXON="28/29", NMD="NMD_escaping_variant", LOFTEE="HC")
    base.update(values)
    return lof_row(**base)


def test_ptc_codon():
    assert P.ptc_codon("ENSP1:p.Arg1730Ter") == 1730
    assert P.ptc_codon("ENSP1:p.(Arg1730*)") == 1730
    assert P.ptc_codon("ENSP1:p.Leu20ProfsTer5") == 24
    assert P.ptc_codon("ENSP1:p.Leu20fs") is None


def test_nmd_52_nt_from_the_last_junction_is_very_strong_despite_vep_plugin():
    codon = (_JUNCTION - 51 + 2) // 3                 # first stop base 52 nt before the junction
    assert _JUNCTION - (3 * codon - 2) + 1 == 52
    pvs1 = met(evaluate_variant(_stop(codon), Context(refs=_exon_refs())), "PVS1")
    assert pvs1 and pvs1[0]["strength"] == "very_strong" and "52 nt" in pvs1[0]["reason"]


def test_nmd_escape_in_last_50_nt_and_last_exon():
    codon = (_JUNCTION - 30 + 2) // 3
    pvs1 = by_code(evaluate_variant(_stop(codon), Context(refs=_exon_refs())), "PVS1")[0]
    assert pvs1["status"] == "met" and pvs1["strength"] != "very_strong" and "逃過" in pvs1["reason"]
    last = by_code(evaluate_variant(_stop(1760, EXON="29/29"), Context(refs=_exon_refs())), "PVS1")[0]
    assert "最後一個 exon" in last["reason"]


def test_frameshift_uses_the_new_stop():
    # the variant is 60 nt upstream but the new stop lands within 50 nt of the junction
    codon = (_JUNCTION - 60 + 2) // 3
    row = _stop(codon, CONSEQUENCE="frameshift_variant", HGVS_P=f"ENSP00000000001.1:p.Arg{codon}GlyfsTer8", NMD="")
    pvs1 = by_code(evaluate_variant(row, Context(refs=_exon_refs())), "PVS1")[0]
    assert "逃過" in pvs1["reason"]


def test_start_proximal_stop_is_not_escape():
    pvs1 = met(evaluate_variant(_stop(20, EXON="1/29"), Context(refs=_exon_refs())), "PVS1")
    assert pvs1 and pvs1[0]["strength"] == "very_strong"


def test_falls_back_to_vep_column_without_exon_table():
    pvs1 = by_code(evaluate_variant(_stop(1000), Context(refs=refs())), "PVS1")[0]
    assert "預測逃過 NMD" in pvs1["reason"]


def test_omim_disease_line_inheritance(monkeypatch):
    from app.services import omim_store
    monkeypatch.setattr(omim_store, "lookup", lambda gene=None, **_: {
        "Inheritance": "", "OMIM_disease": "Deafness, autosomal recessive 16 (603720)(AR)\n"
                                           "[Blood group trait] (111000)(AD)"})
    assert runner.omim_moi("STRC") == {"AR"}
    assert runner.disease_line_moi("Some syndrome (123456)(AD, AR)") == {"AD", "AR"}
    assert runner.disease_line_moi("Old line without MIM (AR)") == set()


def _common_plp(**values):
    base = dict(CLINVAR_SIG="Pathogenic", CLINVAR_STARS="2", GNOMAD_G_AF="0.34",
                CONSEQUENCE="frameshift_variant", HGVS_P="ENSP00000000001.1:p.Leu20fs")
    base.update(values)
    return lof_row(**base)


def test_clinvar_plp_above_ba1_is_not_adopted():
    result = evaluate_variant(_common_plp(), Context(refs=refs()))
    assert met(result, "BA1") and result["classification"] == "Benign"
    assert result["clinvar_rule"]["action"] == "overridden"
    assert result["rank_points"] == -7
    assert any(f["code"] == "clinvar_overridden" for f in result["flags"])


def test_clinvar_plp_on_ba1_exception_list_stays_protected():
    key = "chr1-1000-C-T"
    result = evaluate_variant(_common_plp(), Context(refs=refs(ba1_exceptions={key})))
    assert not met(result, "BA1") and result["clinvar_rule"]["action"] == "adopt"
    low = evaluate_variant(_common_plp(GNOMAD_G_AF="0.03"), Context(refs=refs()))
    assert result["rank_points"] >= 10 and low["clinvar_rule"]["action"] == "adopt"


def test_inhouse_common_penalty_not_exempt_when_common():
    variant = {"id": "chr4-102901325-AAC-A", "CLNSIG": "Pathogenic", "clinvar_stars": 2, "AF": 0.34,
               "inhouse_common": {"level": 2, "target": 8}}
    assert inhouse_common.apply_penalty(variant, scope="in_house", score=10, criteria_text="") == 2
    rare = dict(variant, AF=0.001, inhouse_common={"level": 2, "target": 8})
    assert inhouse_common.apply_penalty(rare, scope="in_house", score=10, criteria_text="") == 10
    v37i = dict(variant, id="chr13-20189473-C-T", AF=0.08, inhouse_common={"level": 2, "target": 8})
    assert inhouse_common.apply_penalty(v37i, scope="in_house", score=10, criteria_text="") == 10


def test_inheritance_uses_the_primary_row_gene():
    lnc = lof_row(GENE="LNC1", TRANSCRIPT="ENST00000000009.1", TRANSCRIPT_TYPE="", IMPACT="MODIFIER",
                  CONSEQUENCE="non_coding_transcript_exon_variant")
    mane = lof_row(GENE="GENE1")
    seen = []
    def moi(gene):
        seen.append(gene)
        return {"AR"} if gene == "GENE1" else None
    [(vid, result)] = list(runner.evaluate_rows({"v": [lnc, mane]}, test_type="WES", refs=refs(), moi_fn=moi))
    assert seen == ["GENE1"]
    assert by_code(result, "PS4")[0]["reason"] != "遺傳模式未知"
