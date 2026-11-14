"""Independent truth coverage and presentation must not manufacture certainty."""
import copy
import gzip
import json
import sqlite3
from decimal import Decimal

import pytest

from app.services import prs_multi, prs_presentation as view
from research.prs import validation as audit
from research.prs.multidisease import provenance
from test_prs_multi import summary,payload


def snapshot(age=20, **history):
    pack=prs_multi.bundle()[0];s=summary();p=payload();p['age']=age;p['history'].update(history)
    result=prs_multi.calculator.evaluate(pack,s,prs_multi.calculator.person_for(s,p))
    return dict(result=result,inputs=dict(bundle=pack,summary=s))


@pytest.mark.parametrize('age,coverage',[(20,'all'),(20.001,'partial'),(25,'partial'),(30,'none'),(70,'none')])
def test_af_zero_annotation_uses_exact_full_horizon_and_preserves_results(age,coverage):
    saved=snapshot(age);before=copy.deepcopy(saved);p=view.assessment_presentation(saved)
    af=p['diseases']['AF'];a=af['risks']['risk_10_year']
    assert a['source_zero']==coverage
    assert saved==before
    if coverage=='all':assert a['ui_text']=='來源零值*' and a['message']==view.ZERO_NOTE
    if coverage=='partial':assert saved['result']['diseases']['AF']['risk_10_year']['risk']>0 and a['ui_text'].endswith('%')
    assert af['curve'][0]['ui_text']=='評估起點'
    assert af['risks']['risk_to_80']['source_zero']==('partial' if age<30 else 'none')


def test_unavailable_reasons_take_precedence_and_period_zero_is_not_rounded_zero():
    p=view.assessment_presentation(snapshot(AF=True))['diseases']['AF']
    assert p['risks']=={} and p['notices']==[]
    assert view.probability_text(.0000001)=='<0.01%'
    assert view.probability_text(.0000001,1)=='<0.1%'
    assert view.risk_annotation(prs_multi.bundle()[0]['models']['AF'],'AF','female',30,40,0)['source_zero']=='none'
    p=view.assessment_presentation(snapshot())['diseases']
    assert view.CANCER_NOTE in p['ALL_CANCERS']['notices']
    assert any('未依百分位遞增' in n for n in p['EPILEPSY']['notices'])


def test_new_annotations_reach_docx_and_snapshot_is_not_rewritten():
    from docx import Document
    from app.services import docx_export
    saved=snapshot();before=copy.deepcopy(saved);doc=Document()
    docx_export._render_health_prs_section(doc,saved)
    text='\n'.join(p.text for p in doc.paragraphs)
    assert '來源零值*' in text and view.ZERO_NOTE in text and view.PARTIAL_ZERO_NOTE in text
    assert view.CANCER_NOTE in text and '數值越高，代表' not in text
    assert saved==before


def test_independent_genotype_reader_does_not_turn_absence_or_failed_gt_into_reference(tmp_path):
    p=tmp_path/'truth.vcf';p.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPUBLIC\n'
        'chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n1\t2\t.\tC\tT\t.\tFAIL\t.\tGT\t0/0\n')
    rows=dict(audit.genotype_rows(p))
    assert rows['1',1,'A','G']==1 and rows['1',2,'C','T'] is None
    assert ('1',3,'A','G') not in rows


def test_source_direction_audit_rejects_model_sign_and_effect_allele_changes(tmp_path):
    model=tmp_path/'model';model.mkdir();sources=tmp_path/'sources';sources.mkdir()
    db=sqlite3.connect(model/'models.sqlite')
    db.executescript('CREATE TABLE locus(id,chrom,pos,ref,alt); CREATE TABLE model(disease,id,ea,weight,freq,row_no);')
    db.execute('INSERT INTO locus VALUES(?,?,?,?,?)',('x',1,1,'A','G'))
    db.execute('INSERT INTO model VALUES(?,?,?,?,?,?)',('CHD','x','A','0.1','0.8',1));db.commit();db.close()
    chd=sources/'chd.gz'
    with gzip.open(chd,'wt') as f:f.write('#pgs_id=PGS004879\nhm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n1\t1\tA\tG\t0.1\n')
    # First CHD row is checked before the other model source files are needed.
    key,row=next(audit.retained_sources(model,sources,chd))
    assert key=='CHD' and row['source_ea']=='A' and row['source_weight']==Decimal('.1')
    db=sqlite3.connect(model/'models.sqlite');db.execute("UPDATE model SET weight='-0.1'");db.commit();db.close()
    with pytest.raises(ValueError,match='direction or weight'):next(audit.retained_sources(model,sources,chd))


def test_missing_runtime_tools_are_unknown_not_assumed(monkeypatch,tmp_path):
    monkeypatch.setattr(provenance.shutil,'which',lambda _:None)
    value=provenance.collect({},None,'score',tmp_path)
    assert value['caller']['status']=='unknown'
    assert all(v['status']=='unavailable' for v in value['tools'].values())


def test_six_disease_independent_score_and_mask_frequency_effects(tmp_path,monkeypatch):
    from research.prs.absolute_risk.registry import DISEASES
    from research.prs.multidisease import common
    model=tmp_path/'model';model.mkdir();sources=tmp_path/'sources';sources.mkdir()
    db=sqlite3.connect(model/'models.sqlite')
    db.executescript('CREATE TABLE locus(id,chrom,pos,ref,alt);CREATE TABLE model(disease,id,ea,weight,freq,row_no);')
    db.executemany('INSERT INTO locus VALUES(?,?,?,?,?)',[('x',1,1,'A','G'),('y',1,2,'C','T')])
    expected={};manifest={'models':{}};distributions={}
    for j,(key,d) in enumerate(DISEASES.items(),1):
        a,b=Decimal(j)/10,-Decimal(j)/5
        db.executemany('INSERT INTO model VALUES(?,?,?,?,?,?)',[(key,'x','A',str(a),'.8',1),(key,'y','T',str(b),'.2',2)])
        path=sources/(d['pgs_id']+'_hmPOS_GRCh38.txt.gz')
        with gzip.open(path,'wt') as f:f.write(f'#pgs_id={d["pgs_id"]}\nhm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n1\t1\tA\tG\t{a}\n1\t2\tT\tC\t{b}\n')
        manifest['models'][key]={'n_loci':2};distributions[key]={str(i):j*(i/367-1) for i in range(735)}
        expected[key]={'total_sum':float(-Decimal(j)*Decimal('.24'))}
    db.commit();db.close();(model/'reference_scores.json').write_text(json.dumps(distributions))
    monkeypatch.setattr(common,'load_models',lambda _: (manifest,'a'*64))
    monkeypatch.setattr(audit,'source_audit',lambda *a:{'status':'PASS','scope':'synthetic test only'})
    vcf=tmp_path/'truth.vcf';vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPUBLIC\n1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\n1\t2\t.\tC\tT\t.\tPASS\t.\tGT\t1/1\n')
    mask=tmp_path/'mask.gz'
    with gzip.open(mask,'wt') as f:f.write('disease\tid\n'+''.join(k+'\tx\n' for k in DISEASES))
    freq=tmp_path/'eas.afreq';freq.write_text('ID\tREF\tALT\tALT_FREQS\tOBS_CT\nx\tA\tG\t.5\t1470\ny\tC\tT\t.2\t1470\n')
    target=tmp_path/'target.json';target.write_text(json.dumps({'diseases':expected}))
    args=(model,sources,vcf,sources/'PGS004879_hmPOS_GRCh38.txt.gz',mask,freq,target)
    rows=audit.benchmark(*args)['diseases']
    for j,(key,row) in enumerate(rows.items(),1):
        assert row['sums']['full']==pytest.approx(-.2*j)
        assert row['sums']['current']==pytest.approx(-.24*j)
        assert row['sums']['eas']==pytest.approx(-.3*j)
        assert row['numerical_match'] and row['truth_missing']==0
    bad={'diseases':copy.deepcopy(expected)};bad['diseases']['CHD']['total_sum']+=1
    target.write_text(json.dumps(bad))
    assert audit.benchmark(*args)['status']=='FAIL'
    vcf.write_text(vcf.read_text().replace('1/1','./.'))
    rows=audit.benchmark(*args)['diseases']
    assert all(r['status']=='PARTIAL_TRUTH_COVERAGE' and 'percentiles' not in r for r in rows.values())


def test_aggregate_excludes_partial_and_rejects_duplicate_comparisons(tmp_path):
    from research.prs.absolute_risk.registry import DISEASES
    full={'status':'COMPLETE_GENOTYPE_COMPARISON','percentiles':{'full':50,'current':60},
          'delta_percentile_pp':{'current':10},'changed_stratum':{'current':True},'numerical_match':False}
    p=tmp_path/'one.json';q=tmp_path/'two.json'
    data={'schema':'PRS_INDEPENDENT_BENCHMARK_V1','model_sha256':'a'*64,'inputs':{'vcf':'a'},'diseases':{k:full for k in DISEASES}}
    p.write_text(json.dumps(data));data['inputs']['vcf']='b';data['diseases']={k:{'status':'PARTIAL_TRUTH_COVERAGE'} for k in DISEASES};q.write_text(json.dumps(data))
    rows=audit.aggregate([p,q])['diseases']
    assert all(r['complete_comparisons']==1 and r['partial_comparisons']==1 and r['numerical_failures']==1 for r in rows.values())
    assert rows['CHD']['metrics']['current']['mean_signed_delta_pp']==10
    with pytest.raises(ValueError,match='Duplicate'):audit.aggregate([p,p])
