"""ACMG rule v1.13: ranking by the rule's own call, PM3 Strong on explicit
in-trans text, PS3 Supporting from ClinVar functional-study text."""
import io
import sqlite3

from app.services.acmg_engine import Context, evaluate_variant
from app.services.acmg_engine import clinvar_cases as CC, ref_store, runner

from test_acmg_engine import by_code, lof_row, met, missense_row, refs


# ── text parsing ────────────────────────────────────────────────────────────
def test_explicit_trans_sentences():
    yes = [
        "This variant was observed in trans with a pathogenic variant in a patient with deafness.",
        "In at least one individual the data is consistent with being in trans (on the opposite chromosome) "
        "from a pathogenic variant (PMID: 123456).",
        "Observed with a pathogenic variant on the opposite allele (in trans) in a patient with hearing loss.",
        "It was found in trans with a likely pathogenic variant in an affected sibling.",
    ]
    no = [
        "This variant was identified in a patient who was compound heterozygous, phase unknown.",
        "It was found in a patient in trans with a variant of uncertain significance.",
        "This variant was presumed to be in trans with a second pathogenic variant in an affected individual.",
        "It is not known whether these variants were in trans in this patient.",
        "This variant has been observed in the compound heterozygous state in individuals with Usher syndrome.",
        "The two variants were found in cis in the proband.",
    ]
    assert all(CC.explicit_trans(t) for t in yes)
    assert not any(CC.explicit_trans(t) for t in no)
    # a ";" clause inherits the hedge of its sentence
    assert not CC.explicit_trans("observed with a pathogenic variant in trans in one patient.",
                                 "Parental testing was not performed; observed with a pathogenic variant in trans "
                                 "in one patient.")


def test_functional_sentences():
    yes = {
        "Experimental studies have shown that this missense change affects GJB2 function (PMID: 12345678, 23456789).":
            ["12345678", "23456789"],
        "Published functional studies demonstrate a damaging effect: reduced enzyme activity (PMID: 9876543).":
            ["9876543"],
        "Studies have shown that this variant results in skipping of exon 5 (PMID: 31234567).": ["31234567"],
        "RNA analysis showed aberrant splicing resulting in exon 12 skipping (PMID 22222222).": ["22222222"],
        "In vitro functional studies showed that the variant protein had 5% residual activity compared to wild type.":
            [],
    }
    no = [
        "Experimental studies are conflicting or provide insufficient evidence to determine the effect of this variant.",
        "Experimental studies have shown that this missense change does not substantially affect BRCA1 function.",
        "Published functional studies demonstrate no damaging effect: normal enzyme activity (PMID: 9876543).",
        "Advanced modeling of protein sequence and biophysical properties (such as structural, functional, and "
        "spatial information) performed at Invitae indicates that this missense variant is expected to disrupt "
        "protein function.",
        "Functional studies of a different missense change at this codon showed reduced activity.",
        "To our knowledge, functional studies have not been performed for this variant.",
        "In silico analysis supports that this missense variant has a deleterious effect on protein function.",
    ]
    for text, pmids in yes.items():
        found = CC.functional_sentences(text)
        assert found and found[0][2] == pmids, text
    assert not any(CC.functional_sentences(t) for t in no)


def test_snippet_is_centred_on_the_match():
    sentence = "x" * 300 + " in trans " + "y" * 300
    snippet = CC._snippet(sentence, 301, 309)
    assert len(snippet) <= CC.SNIPPET_CHARS and "in trans" in snippet
    assert snippet.startswith("…") and snippet.endswith("…")


XML = """<ClinVarVariationRelease>
 <VariationArchive VariationID="165298">
  <ClassifiedRecord><ClinicalAssertionList>
   <ClinicalAssertion ID="1">
    <ClinVarAccession Accession="SCV000200470" OrgID="500031"/>
    <Classification><GermlineClassification>Pathogenic</GermlineClassification>
     <Comment>This variant has been observed in the compound heterozygous state in individuals with deafness. It was observed in trans with a pathogenic variant in one affected individual. Experimental studies have shown that this missense change affects STRC function (PMID: 12345678).</Comment>
    </Classification>
   </ClinicalAssertion>
   <ClinicalAssertion ID="2">
    <ClinVarAccession Accession="SCV2" OrgID="77"/>
    <Classification><GermlineClassification>Likely pathogenic</GermlineClassification>
     <Comment>This variant was identified in a compound heterozygous patient; phase unknown. Functional studies have not been performed.</Comment>
    </Classification>
   </ClinicalAssertion>
   <ClinicalAssertion ID="3">
    <ClinVarAccession Accession="SCV3" OrgID="88"/>
    <Classification><GermlineClassification>Benign</GermlineClassification>
     <Comment>Functional studies showed reduced activity in patient fibroblasts (PMID: 99999999). Observed in trans with a pathogenic variant in a patient.</Comment>
    </Classification>
   </ClinicalAssertion>
  </ClinicalAssertionList></ClassifiedRecord>
 </VariationArchive>
</ClinVarVariationRelease>"""


def test_parser_and_references_schema3(tmp_path):
    (vid, _cases, biallelic, functional), = list(CC.iter_records(io.BytesIO(XML.encode())))
    assert vid == "165298" and set(biallelic) == {"500031", "77"}            # benign submitter dropped
    assert biallelic["500031"]["trans"] == 1 and "in trans" in biallelic["500031"]["snippet"]
    assert biallelic["77"]["trans"] == 0
    assert set(functional) == {"500031"} and functional["500031"]["pmids"] == ["12345678"]
    xml = tmp_path / "v.xml"
    xml.write_text(XML, encoding="utf-8")
    out = tmp_path / CC.CASES_DB_NAME
    stats = CC.build(xml, out, progress=lambda *_: None)
    assert stats["pm3_variants"] == 1 and stats["ps3_variants"] == 1
    r = ref_store.SqliteReferences(tmp_path / ref_store.REFS_DB_NAME, cases_db=out, ba1_exceptions=set())
    got = r.clinvar_biallelic("165298")
    assert got["submitters"] == 2 and got["trans"] == 1 and got["detail"][0]["trans"]
    assert got["detail"][0]["scv"] == "SCV000200470"                          # trans submitter listed first
    ps3 = r.clinvar_functional("165298")
    assert ps3["submitters"] == 1 and ps3["detail"][0]["pmids"] == ["12345678"]
    assert r.clinvar_functional("1") == {} and r.clinvar_functional("") == {}


def test_schema2_table_still_reads_pm3_and_reports_ps3_not_evaluated(tmp_path):
    out = tmp_path / CC.CASES_DB_NAME
    con = sqlite3.connect(out)
    con.executescript("""
        CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE cases (variation_id TEXT PRIMARY KEY, submitters INTEGER NOT NULL,
            individuals INTEGER NOT NULL, submitter_ids TEXT NOT NULL, pm3_submitters INTEGER NOT NULL DEFAULT 0,
            pm3_hom INTEGER NOT NULL DEFAULT 0, pm3_comphet INTEGER NOT NULL DEFAULT 0,
            pm3_detail TEXT NOT NULL DEFAULT '[]') WITHOUT ROWID;
        INSERT INTO cases VALUES ('5', 0, 0, '', 2, 0, 2, '[]');
    """)
    con.commit()
    con.close()
    r = ref_store.SqliteReferences(tmp_path / ref_store.REFS_DB_NAME, cases_db=out, ba1_exceptions=set())
    assert r.clinvar_biallelic("5") == {"submitters": 2, "hom": 0, "comphet": 2, "trans": 0, "detail": []}
    assert r.clinvar_functional("5") is None


# ── engine ─────────────────────────────────────────────────────────────────
def _bi(n, trans=0):
    detail = [{"scv": "SCV2", "snippet": "compound heterozygous patient", "trans": False},
              {"scv": "SCV000200470", "snippet": "observed in trans with a pathogenic variant", "trans": True}]
    return {"123": {"submitters": n, "hom": 0, "comphet": 1, "trans": trans, "detail": detail}}


def test_pm3_strong_on_explicit_trans():
    variant = missense_row(CLINVAR_VARIATION_ID="123")
    strong = met(evaluate_variant(variant, Context(moi={"AR"}, refs=refs(biallelic=_bi(1, trans=1)))), "PM3")
    assert strong and strong[0]["strength"] == "strong" and strong[0]["points"] == 4
    assert "明確寫" in strong[0]["reason"] and strong[0]["reason"].index("SCV000200470") < strong[0]["reason"].index("SCV2")
    plain = met(evaluate_variant(variant, Context(moi={"AR"}, refs=refs(biallelic=_bi(2)))), "PM3")
    assert plain[0]["strength"] == "moderate"
    common = evaluate_variant(missense_row(CLINVAR_VARIATION_ID="123", GNOMAD_E_AF="0.02"),
                              Context(moi={"AR"}, refs=refs(biallelic=_bi(1, trans=1))))
    assert by_code(common, "PM3")[0]["status"] == "not_met"


def _fn(n=1):
    return {"123": {"submitters": n, "detail": [
        {"scv": "SCV000200470", "snippet": "Experimental studies have shown that this missense change affects function",
         "pmids": ["12345678", "23456789"]}]}}


def test_ps3_supporting_from_functional_text():
    variant = missense_row(CLINVAR_VARIATION_ID="123")
    r = evaluate_variant(variant, Context(refs=refs(functional=_fn())))
    ps3 = met(r, "PS3")
    assert ps3 and ps3[0]["strength"] == "supporting" and ps3[0]["points"] == 1
    assert "12345678" in ps3[0]["reason"] and "PS3_Supporting" in r["criteria_text"]
    assert by_code(evaluate_variant(variant, Context(refs=refs())), "PS3")[0]["status"] == "not_evaluated"
    assert by_code(evaluate_variant(variant, Context(refs=refs(functional={}))), "PS3")[0]["status"] == "not_met"
    common = evaluate_variant(missense_row(CLINVAR_VARIATION_ID="123", GNOMAD_E_AF="0.02"),
                              Context(refs=refs(functional=_fn())))
    assert by_code(common, "PS3")[0]["status"] == "not_met"


def test_ps3_not_with_pvs1_and_no_bs3():
    lof = evaluate_variant(lof_row(CLINVAR_VARIATION_ID="123"), Context(refs=refs(functional=_fn())))
    assert met(lof, "PVS1") and by_code(lof, "PS3")[0]["status"] == "not_met"
    assert "PVS1" in by_code(lof, "PS3")[0]["reason"]
    assert not by_code(lof, "BS3")


def test_ps3_respects_vcep_spec():
    spec = {"GENE1": {"criteria": {"PS3": []}, "vcep": "Test VCEP", "svi": "GN000"}}
    r = evaluate_variant(missense_row(CLINVAR_VARIATION_ID="123"),
                         Context(refs=refs(functional=_fn(), vcep_specs=spec)))
    assert by_code(r, "PS3")[0]["status"] == "not_met" and "不使用 PS3" in by_code(r, "PS3")[0]["reason"]


def test_case_recurrence_delegates_functional():
    base = refs(functional=_fn())
    case = runner.CaseRecurrence(base, None, "S1", None, None)
    assert case.clinvar_functional("123")["submitters"] == 1
