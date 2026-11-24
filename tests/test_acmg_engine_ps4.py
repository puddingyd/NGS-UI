"""ACMG rule engine v1.6: PS4 by proband counting (in-house Causative cases + ClinVar affected cases)."""
import gzip
import io
import sqlite3

from app.services.acmg_engine import Context, evaluate_variant
from app.services.acmg_engine import clinvar_cases as CC, ref_store, runner

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
