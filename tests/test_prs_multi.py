"""Synthetic patients only; real published parameters are independently pinned."""
import copy
import datetime
import hashlib
import io
import json
from pathlib import Path

import pytest

from app.services import prs, prs_multi, sample_layout
from research.prs.absolute_risk import multi, engine
from research.prs.absolute_risk.registry import DISEASES

ROOT = Path(__file__).resolve().parents[1]


def summary():
    diseases = {}
    for key, d in DISEASES.items():
        n = 610305 if key == 'CHD' else 1000
        diseases[key] = dict(pgs_id=d['pgs_id'], n_loci=n, observed_loci=n, missing_loci=0,
            total_sum=.25, observed_sum=.25, imputed_sum=0., eas_z=1., eas_percentile=97.,
            intervene_group='>=95', qc_ready=True)
    return dict(schema=multi.SCHEMA, sample='SYNTHETIC', reference_id=multi.REFERENCE,
        model_manifest_sha256='a'*64, diseases=diseases,
        ancestry_status=dict(sample='SYNTHETIC', ancestry='EAS', ancestry_low_confidence='False',
                             qc_status='PASS', pca_ld_pruned_markers=20000))


def payload(**changes):
    return dict(age=50.5, sex='female', assessment_date=datetime.date.today().isoformat(),
                history={key: False for key in DISEASES}, **changes)


def test_all_parameters_reproduce_and_risks_are_conditional():
    pack, _ = prs_multi.bundle()
    proof = json.loads((ROOT/'research/prs/releases/tw_multi_v1/reproduction.json').read_text())
    assert sum(p['rows'] for p in proof['diseases'].values()) == 1232
    for key, evidence in proof['diseases'].items():
        assert evidence['max_error_pp'] <= .0050001
        assert len(evidence['comparisons']) == evidence['rows']
        assert max(abs(r['expected']-r['actual']) for r in evidence['comparisons']) == evidence['max_error_pp']
        b = pack['models'][key]
        for sex in DISEASES[key]['sexes']:
            for percentile in (0,20,40,60,80,90,95,100):
                r = engine.integrate(b, sex, percentile, 50.5, 80)
                assert 0 <= r['risk'] <= 1
                assert r['curve'][0]['risk'] == 0
                assert r['curve'][-1]['age'] == 80
                assert [p['risk'] for p in r['curve']] == sorted(p['risk'] for p in r['curve'])
    person = multi.person_for(summary(), payload())
    result = multi.evaluate(pack, summary(), person)
    assert all(r['available'] for r in result['diseases'].values())
    assert all(r['risk_10_year']['risk'] <= r['risk_20_year']['risk'] <= r['risk_to_80']['risk'] for r in result['diseases'].values())


@pytest.mark.parametrize('age', [20,50.5,60,60.001,70,70.001,79.99])
def test_twenty_year_horizon_uses_full_competing_risk_interval(age):
    pack=prs_multi.bundle()[0];s=summary();p=payload();p['age']=age
    result=multi.evaluate(pack,s,multi.person_for(s,p))
    for key,r in result['diseases'].items():
        assert (r['risk_20_year'] is not None) == (age<=60)
        assert (r['risk_10_year'] is not None) == (age<=70)
        assert (r['reference']['risk_20_year'] is not None) == (age<=60)
        if age<=60:
            assert r['risk_20_year']==engine.integrate(pack['models'][key],'female',97,age,age+20)
            assert r['risk_10_year']['risk']<=r['risk_20_year']['risk']<=r['risk_to_80']['risk']
            assert r['risk_20_year']['curve'][0]['risk']==0
            assert r['risk_20_year']['curve'][-1]['age']==age+20
        if age==60:assert r['risk_20_year']==r['risk_to_80']


def test_per_disease_history_sex_and_qc_do_not_block_others():
    s = summary(); p = payload(); p['sex'] = 'male'
    p['history'].update(AF=True, T2D=None)
    s['diseases']['EPILEPSY'].update(observed_loci=990, missing_loci=10, qc_ready=False)
    result = multi.evaluate(prs_multi.bundle()[0], s, multi.person_for(s, p))['diseases']
    assert {k for k,r in result.items() if r['available']} == {'CHD','ALL_CANCERS'}
    assert '女性' in result['BREAST']['reason']
    assert '復發' in result['AF']['reason']
    assert '確認' in result['T2D']['reason']
    assert '品質' in result['EPILEPSY']['reason']
    p['history'].update(BREAST=True, ALL_CANCERS=False)
    with pytest.raises(ValueError, match='不一致'): multi.person_for(s,p)


@pytest.mark.parametrize('change', [dict(age=20-1e-6),dict(age=80),dict(age=True),dict(sex='unknown'),
    dict(assessment_date='9999-01-01'),dict(history={}),dict(percentile=100)])
def test_bad_clinical_inputs(change):
    p = payload(); p.update(change)
    with pytest.raises((ValueError, TypeError)): multi.person_for(summary(),p)


def test_ancestry_and_model_guards_and_older_horizon():
    s = summary(); s['ancestry_status']['sample'] = 'OTHER'
    assert not any(r['available'] for r in multi.evaluate(prs_multi.bundle()[0],s,multi.person_for(s,payload()))['diseases'].values())
    s = summary(); s['diseases']['CHD']['pgs_id'] = 'PGS004870'
    with pytest.raises(ValueError): multi.score_rows(s)
    s = summary(); p = payload(); p['age']=75
    r = multi.evaluate(prs_multi.bundle()[0],s,multi.person_for(s,p))
    assert all(v['risk_10_year'] is None and v['risk_to_80']['risk']>=0 for v in r['diseases'].values())


@pytest.fixture
def linked(tmp_path, monkeypatch):
    from app.services import prs_jobs
    monkeypatch.setattr(sample_layout,'state_dir',lambda sid,**k: tmp_path/sid)
    monkeypatch.setattr(sample_layout,'state_file',lambda sid,name,**k: tmp_path/sid/f'{sid}.{name}')
    monkeypatch.setattr(prs_jobs,'overview',lambda sid: dict(job=None,can_start=False,candidates=[]))
    prs.atomic_json(tmp_path/'SYNTHETIC/SYNTHETIC.sample_metadata.json',dict(lis_id='SYNTHETIC'))
    s = summary(); model_dir=tmp_path/'models';model_dir.mkdir()
    prs.atomic_json(model_dir/'manifest.json',dict(models={k:dict(pgs_id=r['pgs_id'],n_loci=r['n_loci']) for k,r in s['diseases'].items()}))
    s['model_manifest_sha256']=prs.read_json(model_dir/'manifest.json')[1]
    cfg=tmp_path/'cfg.json';prs.atomic_json(cfg,dict(multi_model_dir=str(model_dir)))
    monkeypatch.setenv('NGS_UI_PRS_CONFIG',str(cfg))
    prs.atomic_json(prs.score_path('SYNTHETIC'),dict(schema='PRS_SCORE_LINK_V1',sample_id='SYNTHETIC',
        source_sample_id='SYNTHETIC',pipeline_source_sha256=None,summary=s))
    return tmp_path


def test_saved_multi_assessment_report_integrity_and_staleness(linked):
    assert prs.status('SYNTHETIC')['available']
    saved = prs.assess('SYNTHETIC',payload(),'test')
    assert len(saved['result']['diseases']) == 6
    assert prs.report_assessment('SYNTHETIC') == saved
    assert prs.report_readiness('SYNTHETIC')['available']
    path = prs.assessment_dir('SYNTHETIC')/(saved['assessment_id']+'.json')
    original = path.read_bytes()
    # Older immutable results keep their source name and remain reproducible;
    # current presentation is supplied independently of that saved string.
    assert saved['result']['diseases']['ALL_CANCERS']['name']=='整體癌症（All cancers）'
    assert next(r for r in prs.status('SYNTHETIC')['diseases'] if r['id']=='ALL_CANCERS')['name']=='全部癌症（All cancers）'
    assert '整體癌症' not in json.dumps(prs_multi.explanation(),ensure_ascii=False)
    assert path.read_bytes()==original
    corrupted=copy.deepcopy(saved);corrupted['result']['diseases']['CHD']['risk_10_year']['risk']=.99
    prs.atomic_json(path,corrupted)
    assert not prs.report_readiness('SYNTHETIC')['available']
    prs.atomic_json(path,saved)
    prs.atomic_json(linked/'models/manifest.json',dict(models={}))
    assert not prs.report_readiness('SYNTHETIC')['available']
    assert not prs.status('SYNTHETIC')['available']


def test_multi_health_text_table_keeps_all_percentiles_and_unavailable_reasons(linked, monkeypatch):
    from docx import Document
    from app.services import docx_export
    p=payload();p['sex']='male';p['history']['AF']=None;p['history']['T2D']=True
    saved=prs.assess('SYNTHETIC',p,'test')
    doc=Document();docx_export._render_health_prs_section(doc,saved)
    assert not doc.tables
    text='\n'.join(p.text for p in doc.paragraphs)
    assert '乳癌' in text and '僅適用女性' in text and '心房顫動' in text
    table_rows=[p.text for p in doc.paragraphs if p.text.startswith('     ') and not p.text.strip().startswith('=')]
    assert len(table_rows)==6
    assert [r.strip().split()[0] for r in table_rows]==['冠狀動脈心臟病','第二型糖尿病','全部癌症','乳癌','心房顫動','癲癇']
    assert '已有此疾病，不估算風險' in table_rows[1]
    assert f'第{docx_export._health_prs_number(saved["result"]["diseases"]["T2D"]["percentile"])}百分位' in table_rows[1]
    header=next(p for p in doc.paragraphs if '多基因風險分數百分位' in p.text)
    assert header.text.splitlines()[1].split('\t')==[
        '     疾病','多基因風險分數百分位','未來10年風險','未來20年風險','至80歲累積風險']
    table_paragraphs=[header]+[p for p in doc.paragraphs if p.text in table_rows]
    for paragraph in table_paragraphs:
        assert [(stop.position.pt,stop.alignment) for stop in paragraph.paragraph_format.tab_stops]==[
            (126.5,0),(247.5,0),(324.5,0),(401.5,0)]
        assert paragraph._p.xpath('./w:pPr/w:autoSpaceDN/@w:val')==['0']
    assert len(table_rows[0].splitlines()[0].split('\t'))==5
    assert len(table_rows[1].splitlines()[0].split('\t'))==3
    assert '本次未提供的風險估計' not in text and '各疾病為獨立估計' not in text
    assert not any('：'+r['limitations'][0] in text for r in saved['result']['diseases'].values() if r['available'])


def test_pre_twenty_year_snapshot_is_preserved_but_needs_reassessment_for_report(linked):
    p=payload();p['history']['AF']=True
    saved=prs.assess('SYNTHETIC',p,'test')
    path=prs.assessment_dir('SYNTHETIC')/(saved['assessment_id']+'.json')
    for r in saved['result']['diseases'].values():
        if r['available']:
            r.pop('risk_20_year');r['reference'].pop('risk_20_year')
    prs.atomic_json(path,saved)
    original=path.read_bytes()
    status=prs.status('SYNTHETIC')
    assert status['latest']==saved and status['available']
    assert status['latest']['result']['person']['history']['AF'] is True
    readiness=prs.report_readiness('SYNTHETIC')
    assert not readiness['available'] and '20年' in readiness['reason']
    assert path.read_bytes()==original
    new=prs.assess('SYNTHETIC',p,'test')
    assert new['assessment_id']!=saved['assessment_id']
    assert prs.report_readiness('SYNTHETIC')['available']
    assert path.read_bytes()==original
