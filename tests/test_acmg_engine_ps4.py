"""ACMG rule engine v1.6: PS4 by proband counting (in-house Causative cases + ClinVar affected cases)."""
import gzip
import io
import sqlite3

from app.services.acmg_engine import Context, evaluate_variant
from app.services.acmg_engine import clinvar_cases as CC, ref_store, runner
from app.services.acmg_engine.engine import load_rulebook

from test_acmg_engine import by_code, met, missense_row, refs


def _ps4(result):
    return by_code(result, "PS4")[0]


def test_ps4_inhouse_counts_and_strengths():
    for n, strength in ((1, "supporting"), (2, "moderate"), (3, "moderate"), (4, "strong"), (7, "strong")):
        r = evaluate_variant(missense_row(), Context(moi={"AD"}, refs=refs(probands={"chr1-1000-C-T": n})))
        e = _ps4(r)
        assert e["status"] == "met" and e["strength"] == strength and f"{n} 位" in e["reason"]
        assert r["ps4_inhouse"] == n
    none = evaluate_variant(missense_row(), Context(moi={"AD"}, refs=refs(probands={}, clinvar_cases={})))
    assert _ps4(none)["status"] == "not_met" and none["ps4_inhouse"] == 0


def test_ps4_needs_dominant_gene_and_pm2():
    data = refs(probands={"chr1-1000-C-T": 5})
    recessive = evaluate_variant(missense_row(), Context(moi={"AR"}, refs=data))
    assert _ps4(recessive)["status"] == "not_met" and "PM3" in _ps4(recessive)["reason"]
    assert recessive["ps4_inhouse"] is None
    mixed = evaluate_variant(missense_row(), Context(moi={"AD", "AR"}, refs=data))
    assert met(mixed, "PS4")
    common = evaluate_variant(missense_row(GNOMAD_E_AF="0.002"), Context(moi={"AD"}, refs=data))
    assert _ps4(common)["status"] == "not_met" and "PM2" in _ps4(common)["reason"]
    unknown = evaluate_variant(missense_row(), Context(refs=data))
    assert _ps4(unknown)["status"] == "not_evaluated"
    nothing = evaluate_variant(missense_row(), Context(moi={"AD"}, refs=refs()))
    assert _ps4(nothing)["status"] == "not_evaluated"


def test_ps4_clinvar_capped_at_supporting_and_stronger_source_wins():
    cv = {"123": {"submitters": 6, "individuals": 20}}
    variant = missense_row(CLINVAR_VARIATION_ID="123")
    only_cv = evaluate_variant(variant, Context(moi={"AD"}, refs=refs(clinvar_cases=cv, probands={})))
    assert met(only_cv, "PS4")[0]["strength"] == "supporting" and "ClinVar 6 個提交者" in _ps4(only_cv)["reason"]
    both = evaluate_variant(variant, Context(moi={"AD"}, refs=refs(clinvar_cases=cv,
                                                                    probands={"chr1-1000-C-T": 2})))
    assert met(both, "PS4")[0]["strength"] == "moderate"            # max, not the sum
    latest = evaluate_variant(missense_row(CLINVAR_LATEST_VARIATION_ID="123"),
                              Context(moi={"AD"}, refs=refs(clinvar_cases=cv)))
    assert met(latest, "PS4")
    assert "PS4_Supporting" in only_cv["criteria_text"]


def test_ps4_respects_vcep_spec():
    spec = {"svi": "GN1", "criteria": {"PS4": []}}
    r = evaluate_variant(missense_row(), Context(moi={"AD"}, refs=refs(probands={"chr1-1000-C-T": 4},
                                                                        vcep_specs={"GENE1": spec})))
    assert not met(r, "PS4") and "規格不使用 PS4" in _ps4(r)["reason"]
    capped = {"svi": "GN1", "criteria": {"PS4": ["supporting", "moderate"]}}
    r = evaluate_variant(missense_row(), Context(moi={"AD"}, refs=refs(probands={"chr1-1000-C-T": 4},
                                                                        vcep_specs={"GENE1": capped})))
    assert met(r, "PS4")[0]["strength"] == "moderate"


XML = """<?xml version="1.0"?>
<ClinVarVariationRelease>
 <VariationArchive VariationID="225368" Accession="VCV000225368">
  <ClassifiedRecord>
   <Classifications><GermlineClassification><Description>Pathogenic</Description></GermlineClassification></Classifications>
   <ClinicalAssertionList>
    <ClinicalAssertion ID="1">
     <ClinVarAccession Accession="SCV1" OrgID="500031" SubmitterName="Lab A"/>
     <RecordStatus>current</RecordStatus>
     <Classification><GermlineClassification>Pathogenic</GermlineClassification></Classification>
     <ObservedInList><ObservedIn>
      <Sample><Origin>germline</Origin><AffectedStatus>yes</AffectedStatus><NumberOfIndividuals>12</NumberOfIndividuals></Sample>
      <Method><MethodType>clinical testing</MethodType></Method>
     </ObservedIn></ObservedInList>
    </ClinicalAssertion>
    <ClinicalAssertion ID="2">
     <ClinVarAccession Accession="SCV2" OrgID="1006" SubmitterName="Lit"/>
     <Classification><GermlineClassification>Pathogenic</GermlineClassification></Classification>
     <ObservedInList><ObservedIn>
      <Sample><Origin>germline</Origin><AffectedStatus>yes</AffectedStatus></Sample>
      <Method><MethodType>literature only</MethodType></Method>
     </ObservedIn></ObservedInList>
    </ClinicalAssertion>
    <ClinicalAssertion ID="3">
     <ClinVarAccession Accession="SCV3" OrgID="77" SubmitterName="Lab B"/>
     <Classification><GermlineClassification>Likely pathogenic</GermlineClassification></Classification>
     <ObservedInList><ObservedIn>
      <Sample><Origin>germline</Origin><AffectedStatus>yes</AffectedStatus></Sample>
      <Method><MethodType>clinical testing</MethodType></Method>
     </ObservedIn></ObservedInList>
    </ClinicalAssertion>
   </ClinicalAssertionList>
  </ClassifiedRecord>
 </VariationArchive>
 <VariationArchive VariationID="9">
  <ClassifiedRecord><ClinicalAssertionList>
   <ClinicalAssertion ID="4">
    <ClinVarAccession Accession="SCV4" OrgID="5"/>
    <Classification><GermlineClassification>Benign</GermlineClassification></Classification>
    <ObservedInList><ObservedIn><Sample><Origin>germline</Origin><AffectedStatus>yes</AffectedStatus></Sample>
     <Method><MethodType>clinical testing</MethodType></Method></ObservedIn></ObservedInList>
   </ClinicalAssertion>
   <ClinicalAssertion ID="5">
    <ClinVarAccession Accession="SCV5" OrgID="6"/>
    <Classification><GermlineClassification>Pathogenic</GermlineClassification></Classification>
    <ObservedInList><ObservedIn><Sample><Origin>somatic</Origin><AffectedStatus>yes</AffectedStatus></Sample>
     <Method><MethodType>clinical testing</MethodType></Method></ObservedIn>
     <ObservedIn><Sample><Origin>germline</Origin><AffectedStatus>no</AffectedStatus></Sample>
     <Method><MethodType>clinical testing</MethodType></Method></ObservedIn></ObservedInList>
   </ClinicalAssertion>
  </ClinicalAssertionList></ClassifiedRecord>
 </VariationArchive>
</ClinVarVariationRelease>
"""


def test_clinvar_cases_parser_and_references(tmp_path):
    found = dict(CC.iter_cases(io.BytesIO(XML.encode())))
    assert found == {"225368": {"500031": 12, "77": 1}}       # literature-only, benign, somatic, unaffected dropped
    xml = tmp_path / "vcv.xml.gz"
    with gzip.open(xml, "wb") as handle:
        handle.write(XML.encode())
    out = tmp_path / CC.CASES_DB_NAME
    stats = CC.build(xml, out, release="2026-11-03", progress=lambda *_: None)
    assert stats["variants"] == 1
    refs_db = tmp_path / ref_store.REFS_DB_NAME          # not built: only the cases DB is present
    r = ref_store.SqliteReferences(refs_db, cases_db=out, ba1_exceptions=set())
    assert r.clinvar_affected("225368") == {"submitters": 2, "individuals": 13}
    assert r.clinvar_affected("1") == {} and "clinvar_cases" in r.available
    assert r.version["clinvar_cases_release"] == "2026-11-03"
    assert ref_store.SqliteReferences(None, ba1_exceptions=set()).clinvar_affected("225368") is None


def _observations(db, rows):
    conn = sqlite3.connect(db)
    conn.execute("CREATE TABLE observations (genome_build TEXT, variant_id TEXT, sample_id TEXT, status TEXT, "
                 "reviewer_user_id INTEGER, reviewer_username TEXT, updated_at TEXT)")
    conn.executemany("INSERT INTO observations VALUES ('hg38', ?, ?, ?, NULL, '', '')", rows)
    conn.commit()
    conn.close()


def test_inhouse_probands_dedupe_patients_and_skip_own_case(tmp_path):
    db = tmp_path / "manual.sqlite"
    _observations(db, [("chr1-1000-C-T", "S1-dragen", "1"), ("chr1-1000-C-T", "S1-nckuh", "1"),
                       ("chr1-1000-C-T", "S2", "1"), ("chr1-1000-C-T", "S3", "2"),
                       ("chr1-1000-C-T", "SELF", "1"), ("chr1-1000-C-T", "SELF_OTHER_RUN", "1")])
    mrn = {"S1-dragen": "111", "S1-nckuh": "111", "S2": "222", "SELF": "999", "SELF_OTHER_RUN": "999"}
    probands = runner.InhouseProbands.from_db(db, patient_of=lambda s: mrn.get(s, ""))
    # S1 twice → one patient; S3 is status 2 (Other); SELF and its other run are the case's own patient
    assert probands.count("chr1-1000-C-T", "SELF") == 2
    assert probands.count("1-1000-c-t", "SELF") == 2                  # key normalised like the registry
    assert probands.count("chr1-2-A-G", "SELF") == 0
    case = runner.CaseRecurrence(refs(), None, "SELF", "", probands)
    r = evaluate_variant(missense_row(), Context(moi={"AD"}, refs=case))
    assert met(r, "PS4")[0]["strength"] == "moderate" and r["ps4_inhouse"] == 2
    assert runner.summary(r)["ps4_inhouse"] == 2
    assert runner.InhouseProbands.from_db(tmp_path / "missing.sqlite") is None


def test_attach_refreshes_ps4_when_inhouse_count_changes(tmp_path, monkeypatch):
    from test_acmg_engine_runner import MISSENSE, write_review
    review = write_review(tmp_path / "S1.snv_indel.review.tsv", [MISSENSE])
    out = tmp_path / "S1.snv_acmg_engine.sqlite"
    shared = refs()
    vid = runner.ui_variant_id(MISSENSE)
    before = runner.InhouseProbands({}, patient_of=lambda s: "")
    runner.build_case(review_tsv=review, out=out, sample_id="S1", test_type="WGS", refs=shared,
                      moi_fn=lambda gene: {"AD"}, probands=before)
    after = runner.InhouseProbands({runner.ui_variant_id(MISSENSE): {"A", "B"}}, patient_of=lambda s: "")
    monkeypatch.setattr(runner, "shared_references", lambda: shared)
    monkeypatch.setattr(runner, "shared_recurrence", lambda: None)
    monkeypatch.setattr(runner, "omim_moi", lambda gene: {"AD"})
    monkeypatch.setattr(runner, "shared_probands", lambda: before)
    variants = {vid: {}}
    runner.attach(variants, results_path=out, review_tsv=review, sample_id="S1", test_type="WGS")
    assert variants[vid]["engine_acmg"]["ps4_inhouse"] == 0 and "ps4_refreshed" not in variants[vid]["engine_acmg"]
    monkeypatch.setattr(runner, "shared_probands", lambda: after)
    variants = {vid: {}}
    counts = runner.attach(variants, results_path=out, review_tsv=review, sample_id="S1", test_type="WGS")
    engine = variants[vid]["engine_acmg"]
    assert counts.get("ps4_refreshed") == 1 and engine["ps4_inhouse"] == 2
    assert "PS4_Moderate" in engine["criteria_text"] and engine["status"] == "current"


# ── v1.7: PM3 from ClinVar biallelic reports ───────────────────────────────
def _bi(n, hom=1, comphet=0):
    return {"123": {"submitters": n, "hom": hom, "comphet": comphet,
                    "detail": [{"scv": "SCV000389235", "snippet": "including in a homozygous state in one individual"}]}}


def test_pm3_from_clinvar_biallelic():
    variant = missense_row(CLINVAR_VARIATION_ID="123")
    one = evaluate_variant(variant, Context(moi={"AR"}, refs=refs(biallelic=_bi(1))))
    pm3 = met(one, "PM3")
    assert pm3 and pm3[0]["strength"] == "supporting"
    assert "SCV000389235" in pm3[0]["reason"] and "homozygous state" in pm3[0]["reason"]
    two = evaluate_variant(variant, Context(moi={"AR"}, refs=refs(biallelic=_bi(3))))
    assert met(two, "PM3")[0]["strength"] == "moderate" and "PM3" in two["criteria_text"]
    assert not by_code(evaluate_variant(variant, Context(moi={"AD"}, refs=refs(biallelic=_bi(3)))), "PM3")
    common = evaluate_variant(missense_row(CLINVAR_VARIATION_ID="123", GNOMAD_E_AF="0.02"),
                              Context(moi={"AR"}, refs=refs(biallelic=_bi(3))))
    assert by_code(common, "PM3")[0]["status"] == "not_met"
    assert by_code(evaluate_variant(variant, Context(moi={"AR"}, refs=refs())), "PM3")[0]["status"] == "not_evaluated"
    none = evaluate_variant(variant, Context(moi={"AR"}, refs=refs(biallelic={})))
    assert by_code(none, "PM3")[0]["status"] == "not_met"


def test_pm3_and_ps4_keep_only_the_stronger():
    variant = missense_row(CLINVAR_VARIATION_ID="123")
    r = evaluate_variant(variant, Context(moi={"AD", "AR"}, refs=refs(biallelic=_bi(2), probands={},
                                                                         clinvar_cases={"123": {"submitters": 3}})))
    assert met(r, "PM3") and not met(r, "PS4")                         # Moderate PM3 beats Supporting PS4
    assert "只計較強者" in by_code(r, "PS4")[0]["reason"]
    r = evaluate_variant(variant, Context(moi={"AD", "AR"}, refs=refs(biallelic=_bi(1),
                                                                         probands={"chr1-1000-C-T": 4})))
    assert met(r, "PS4") and not met(r, "PM3")                         # Strong in-house PS4 wins


XML_PM3 = """<ClinVarVariationRelease>
 <VariationArchive VariationID="225368">
  <ClassifiedRecord><ClinicalAssertionList>
   <ClinicalAssertion ID="1">
    <ClinVarAccession Accession="SCV000389235" OrgID="500031"/>
    <Classification><GermlineClassification>Pathogenic</GermlineClassification>
     <Comment>The GALC c.1901T&gt;C (p.Leu634Ser) variant has been identified in 12 individuals affected with Krabbe disease, including in a homozygous state in one individual and in a compound heterozygous state in eight individuals. It is reported in gnomAD with no homozygotes.</Comment>
    </Classification>
    <ObservedInList><ObservedIn><Sample><Origin>germline</Origin><AffectedStatus>yes</AffectedStatus></Sample>
     <Method><MethodType>literature only</MethodType></Method></ObservedIn></ObservedInList>
   </ClinicalAssertion>
   <ClinicalAssertion ID="2">
    <ClinVarAccession Accession="SCV2" OrgID="77"/>
    <Classification><GermlineClassification>Likely pathogenic</GermlineClassification></Classification>
    <ObservedInList><ObservedIn><Sample><Origin>germline</Origin><AffectedStatus>yes</AffectedStatus></Sample>
     <Method><MethodType>clinical testing</MethodType></Method>
     <ObservedData><Attribute Type="CompoundHeterozygote" integerValue="2"/></ObservedData></ObservedIn></ObservedInList>
   </ClinicalAssertion>
   <ClinicalAssertion ID="3">
    <ClinVarAccession Accession="SCV3" OrgID="88"/>
    <Classification><GermlineClassification>Uncertain significance</GermlineClassification>
     <Comment>No homozygous individuals were observed in population databases. Unaffected carrier parents were heterozygous.</Comment></Classification>
   </ClinicalAssertion>
  </ClinicalAssertionList></ClassifiedRecord>
 </VariationArchive>
</ClinVarVariationRelease>"""


def test_parser_collects_biallelic_text_and_counts(tmp_path):
    records = list(CC.iter_records(io.BytesIO(XML_PM3.encode())))
    assert len(records) == 1
    vid, ps4_found, biallelic, functional = records[0]
    assert vid == "225368" and ps4_found == {"77": 1}                 # literature-only is not a PS4 case
    assert set(biallelic) == {"500031", "77"}                          # SCV3 only negations / population text
    assert biallelic["500031"]["hom"] and biallelic["500031"]["comphet"]
    assert "homozygous state" in biallelic["500031"]["snippet"]
    assert biallelic["77"]["comphet"] and "compound heterozygote" in biallelic["77"]["snippet"]
    xml = tmp_path / "v.xml"
    xml.write_text(XML_PM3, encoding="utf-8")
    out = tmp_path / CC.CASES_DB_NAME
    assert CC.build(xml, out, progress=lambda *_: None)["pm3_variants"] == 1
    r = ref_store.SqliteReferences(tmp_path / ref_store.REFS_DB_NAME, cases_db=out, ba1_exceptions=set())
    got = r.clinvar_biallelic("225368")
    assert got["submitters"] == 2 and got["hom"] == 1 and got["comphet"] == 2
    assert got["detail"][0]["scv"] and r.clinvar_biallelic("1") == {}
    assert r.clinvar_affected("225368") == {"submitters": 1, "individuals": 1}


# ── v1.7: ranking uses the P-KNN LLR ──────────────────────────────────────
def test_rank_points_use_pknn_llr():
    pp3 = evaluate_variant(missense_row(PKNN_EVIDENCE="PP3_Moderate", PKNN_LLR="3.442"), Context(refs=refs()))
    assert met(pp3, "PP3")[0]["points"] == 2
    assert pp3["rank_points"] == pp3["engine_points"] - 2 + 3.442 and pp3["classification"] == pp3["engine_classification"]
    weak = evaluate_variant(missense_row(PKNN_EVIDENCE=".", PKNN_LLR="0.6"), Context(refs=refs()))
    assert weak["rank_points"] == weak["engine_points"] + 0.6
    bp4 = evaluate_variant(missense_row(PKNN_EVIDENCE="BP4_Supporting", PKNN_LLR="-1.5"), Context(refs=refs()))
    assert met(bp4, "BP4") and bp4["rank_points"] == bp4["engine_points"] + 1 - 1.5
    no_llr = evaluate_variant(missense_row(PKNN_EVIDENCE="PP3_Moderate"), Context(refs=refs()))
    assert no_llr["rank_points"] == no_llr["engine_points"]
    clinvar_p = missense_row(PKNN_EVIDENCE="PP3_Moderate", PKNN_LLR="3.442", CLINVAR_SIG="Pathogenic", CLINVAR_STARS="2")
    own = evaluate_variant(clinvar_p, Context(refs=refs()))                 # v1.13: ClinVar does not move ranking
    assert own["rank_points"] == pp3["rank_points"] and own["points"] == own["engine_points"]
    adopt = load_rulebook().copy()
    adopt["clinvar"] = dict(adopt["clinvar"], rank_mode="adopt")
    adopted = evaluate_variant(clinvar_p, Context(refs=refs()), adopt)
    assert adopted["points"] >= 10 and adopted["rank_points"] >= 10
    assert runner.summary(pp3)["rank_points"] == pp3["rank_points"]
