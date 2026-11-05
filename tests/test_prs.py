import copy
import datetime
import hashlib
import json
import io
from pathlib import Path

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient
from starlette.middleware.sessions import SessionMiddleware

from app import config
from app.auth import current_user
from app.routers import prs as router
from app.services import prs, sample_layout

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def setup(tmp_path, monkeypatch):
    # Deliberately artificial parameter values, isolated in pytest's temporary dir.
    # No sample/model fixtures are made available through the production router.
    fixtures = ROOT / 'research/prs/examples'
    bundle = json.loads((fixtures/'synthetic_bundle.json').read_text())
    bundle.update(status='unvalidated_external', model_id='intervene_estbb_tw2019_A_exact_v1',
                  population='Taiwan', age_min=20,
                  reproduction={'status':'PASS','rows':224,'max_error_pp':.0049})
    bundle['bins'] = [r for r in bundle['bins'] if r['age_start'] >= 20]
    for r in bundle['bins']:
        r['hr']['>=95'] = 2
    summary = json.loads((fixtures/'synthetic_summary.json').read_text())
    summary['sample'] = 'S1'
    summary['scores'][0].update(eas_percentile=97.0, intervene_group='>=95')
    model = tmp_path/'bundle.json'
    model.write_text(json.dumps(bundle))
    monkeypatch.setattr(config, 'PRS_BUNDLE', model)
    monkeypatch.setattr(prs, 'SCORE_DEFINITION', bundle['score_definition'])
    from app.services import prs_jobs
    monkeypatch.setattr(prs_jobs, 'overview', lambda sid: {'job':None,'can_start':False,'candidates':[]})
    monkeypatch.setattr(prs, 'RELEASED_BUNDLE_SHA256', hashlib.sha256(model.read_bytes()).hexdigest())
    monkeypatch.setattr(sample_layout, 'state_dir', lambda sid, **kw: tmp_path/sid)
    monkeypatch.setattr(sample_layout, 'state_file', lambda sid,name,**kw: tmp_path/sid/f'{sid}.{name}')
    prs.atomic_json(tmp_path/'S1/S1.sample_metadata.json', {'lis_id':'S1'})
    score = tmp_path/'summary.json'
    score.write_text(json.dumps(summary))
    prs.import_score('S1', score)
    app = FastAPI()
    app.add_middleware(SessionMiddleware, secret_key='test-only')
    app.include_router(router.router)
    client = TestClient(app)
    payload = {'age':50.5, 'sex':'male', 'assessment_date':datetime.date.today().isoformat(), 'prior_chd':False}
    return client, app, tmp_path, bundle, score, payload


def login(app):
    app.dependency_overrides[current_user] = lambda: {'username':'researcher'}


def test_auth_and_unknown_case(setup):
    client, app, *_ = setup
    assert client.get('/api/samples/S1/prs').status_code == 401
    assert client.post('/api/samples/S1/prs/assessments', json={}).status_code == 401
    assert client.post('/api/samples/S1/prs/jobs', json={}).status_code == 401
    assert client.post('/api/samples/S1/prs/jobs/unknown/cancel').status_code == 401
    assert client.get('/api/samples/S1/prs/jobs/unknown/log').status_code == 401
    login(app)
    assert client.get('/api/samples/S2/prs').status_code == 404
    with pytest.raises(ValueError): prs.status('../S1')


def test_calculation_reference_and_immutable_snapshot(setup):
    client, app, root, _, _, payload = setup
    login(app)
    assert client.get('/api/samples/S1/prs').json()['available']
    response = client.post('/api/samples/S1/prs/assessments',json=payload)
    assert response.status_code == 200, response.text
    first = response.json()
    result = first['result']
    assert result['percentile'] == 97
    assert result['risk_10_year']['risk'] > result['reference']['risk_10_year']['risk']
    assert result['risk_to_80']['curve'][0] == {'age':50.5,'risk':0,'event_free_survival':1}
    assert result['reference']['risk_to_80']['curve'][-1]['age'] == 80
    assert result['clinical_validation'] == 'NOT_VALIDATED'
    assert result['uncertainty_interval'] is None
    path = prs.assessment_dir('S1')/f"{first['assessment_id']}.json"
    before = path.read_bytes()
    second = client.post('/api/samples/S1/prs/assessments',json={**payload,'age':60}).json()
    assert first['assessment_id'] != second['assessment_id']
    assert path.read_bytes() == before
    assert len(list(prs.assessment_dir('S1').glob('*.json'))) == 2
    assert client.get('/api/samples/S1/prs').json()['latest']['assessment_id'] == second['assessment_id']
    assert first['inputs']['person']['age'] == 50.5


@pytest.mark.parametrize('change',[
    {'prior_chd':True}, {'prior_chd':None}, {'age':19.99}, {'age':80}, {'age':True},
    {'sex':'unknown'}, {'assessment_date':'9999-01-01'}, {'percentile':99}, {'allow_synthetic':True},
])
def test_invalid_inputs_do_not_publish(setup,change):
    client, app, *_, payload = setup
    login(app)
    response=client.post('/api/samples/S1/prs/assessments',json={**payload,**change})
    assert response.status_code == 400
    assert not prs.assessment_dir('S1').exists()


def test_horizon_not_truncated(setup):
    client, app, *_, payload = setup
    login(app)
    for age, available in [(20,True),(50.5,True),(60,True),(60.01,True),(70,True),(70.01,False),(79.99,False)]:
        result=client.post('/api/samples/S1/prs/assessments',json={**payload,'age':age}).json()['result']
        assert (result['risk_10_year'] is not None) == available
        assert (result['reference']['risk_10_year'] is not None) == available
        assert (result['risk_20_year'] is not None) == (age <= 60)
        assert (result['reference']['risk_20_year'] is not None) == (age <= 60)
        if age <= 60:
            assert result['risk_10_year']['risk'] <= result['risk_20_year']['risk'] <= result['risk_to_80']['risk']
            assert result['risk_20_year']['curve'][-1]['age'] == age+20
        if age == 60:
            assert result['risk_20_year'] == result['risk_to_80']
        assert result['risk_to_80']['risk'] > 0


def test_old_chd_assessment_keeps_values_and_requires_twenty_year_recalculation(setup):
    client, app, *_, payload = setup
    login(app)
    saved=client.post('/api/samples/S1/prs/assessments',json=payload).json()
    result=saved['result'];result.pop('risk_20_year');result.pop('risk_20_year_unavailable_reason')
    result['reference'].pop('risk_20_year');result['risk_difference_pp'].pop('twenty_year')
    path=prs.assessment_dir('S1')/(saved['assessment_id']+'.json');prs.atomic_json(path,saved)
    original=path.read_bytes()
    status=client.get('/api/samples/S1/prs').json()
    assert status['latest']==saved and status['available']
    readiness=client.get('/api/samples/S1/prs/report-readiness').json()
    assert not readiness['available'] and '20年' in readiness['reason']
    assert path.read_bytes()==original


def test_missing_or_synthetic_bundle_never_uses_fallback(setup):
    client, app, root, bundle, _, payload = setup
    login(app)
    for value in [None, {**bundle,'status':'synthetic'}, {**bundle,'reproduction':{}}]:
        if value is None: config.PRS_BUNDLE.unlink()
        else: config.PRS_BUNDLE.write_text(json.dumps(value))
        status=client.get('/api/samples/S1/prs').json()
        assert not status['available']
        assert status['latest'] is None
        assert client.post('/api/samples/S1/prs/assessments',json=payload).status_code == 409
    assert not prs.assessment_dir('S1').exists()


def test_source_and_model_changes_hide_stale_results(setup):
    client, app, root, bundle, score, payload = setup
    login(app)
    client.post('/api/samples/S1/prs/assessments',json=payload).raise_for_status()
    data=json.loads(score.read_text());data['scores'][0].update(eas_percentile=50,intervene_group='40-60')
    score.write_text(json.dumps(data));prs.import_score('S1',score)
    state=client.get('/api/samples/S1/prs').json()
    assert state['stale'] and state['latest'] is None
    assert len(list(prs.assessment_dir('S1').glob('*.json'))) == 1
    client.post('/api/samples/S1/prs/assessments',json=payload).raise_for_status()
    bundle['bins'][-1]['hr']['>=95']=3
    config.PRS_BUNDLE.write_text(json.dumps(bundle))
    assert client.get('/api/samples/S1/prs').json()['stale']


def test_summary_cannot_be_linked_to_wrong_patient(setup):
    _, _, root, _, score, _ = setup
    data=json.loads(score.read_text());data['sample']='OTHER';score.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='source_sample_id'):prs.import_score('S1',score)
    prs.atomic_json(root/'S1/S1.pipeline_source.json',{'source_sample_id':'OTHER'})
    prs.import_score('S1',score)
    assert prs.load_score('S1')[0]['sample']=='OTHER'
    prs.atomic_json(root/'S1/S1.pipeline_source.json',{'source_sample_id':'NEW_SOURCE'})
    assert not prs.status('S1')['available']


def test_low_confidence_ancestry_and_wrong_lock_rejected(setup):
    client, app, _, _, score, payload = setup
    login(app)
    initial=json.loads(score.read_text())
    for key,value in [('ancestry_status',{}),('model_lock_sha256','incorrect')]:
        data={**initial,key:value};score.write_text(json.dumps(data));prs.import_score('S1',score)
        assert not prs.status('S1')['available']
        assert client.post('/api/samples/S1/prs/assessments',json=payload).status_code == 400


def test_concurrent_input_change_rejected(setup,monkeypatch):
    client, app, root, _, _, payload = setup
    login(app)
    evaluate=prs.calculator.evaluate
    def changed(*args,**kw):
        result=evaluate(*args,**kw)
        prs.atomic_json(root/'S1/S1.sample_metadata.json',{'lis_id':'S1','updated':True})
        return result
    monkeypatch.setattr(prs.calculator,'evaluate',changed)
    assert client.post('/api/samples/S1/prs/assessments',json=payload).status_code == 409
    assert not prs.assessment_dir('S1').exists()


def test_score_visible_without_bundle_or_ancestry(setup):
    client, app, _, _, score, _ = setup
    login(app)
    data=json.loads(score.read_text());data['ancestry_status']={'status':'NOT_ASSESSED'}
    score.write_text(json.dumps(data));prs.import_score('S1',score)
    config.PRS_BUNDLE.unlink()
    result=client.get('/api/samples/S1/prs').json()
    assert result['score']['percentile']==97
    assert result['score']['raw_score']==data['scores'][0]['total_sum']
    assert not result['available'] and result['latest'] is None
    assert result['diseases'][0]['id']=='CHD'
    assert result['score']['qc_ready'] is False
    text=json.dumps(result['explanation'],ensure_ascii=False)
    assert 'PGS004879' in text and '610,305' in text and '735' in text
    assert '不穩定心絞痛' in text and '穩定心絞痛' in text and '心跳停止' in text
    assert 'Model A' not in text and 'outcome validation' not in text


def test_explicit_failed_ancestry_qc_rejects_risk(setup):
    client, app, _, _, score, payload=setup
    login(app)
    data=json.loads(score.read_text());data['ancestry_status']['qc_status']='FAIL'
    score.write_text(json.dumps(data));prs.import_score('S1',score)
    result=client.get('/api/samples/S1/prs').json()
    assert result['score'] and not result['available']
    assert client.post('/api/samples/S1/prs/assessments',json=payload).status_code==400


def test_health_readiness_requires_current_verified_saved_assessment(setup):
    client, app, root, _, score, payload = setup
    url='/api/samples/S1/prs/report-readiness'
    assert client.get(url).status_code == 401
    login(app)
    assert client.get('/api/samples/S2/prs/report-readiness').status_code == 404
    assert not client.get(url).json()['available']  # Score alone is insufficient.
    snapshot=client.post('/api/samples/S1/prs/assessments',json=payload).json()
    assert client.get(url).json()['assessment_id'] == snapshot['assessment_id']
    path=prs.assessment_dir('S1')/f"{snapshot['assessment_id']}.json"
    before=path.read_bytes()
    assert prs.report_assessment('S1') == snapshot and path.read_bytes() == before
    altered=copy.deepcopy(snapshot); altered['result']['risk_10_year']['risk']=.999
    prs.atomic_json(path,altered)
    assert not client.get(url).json()['available']
    path.write_bytes(before)
    prs.atomic_json(root/'S1/S1.sample_metadata.json',{'lis_id':'S1','Test':'WGS'})
    assert not client.get(url).json()['available']


@pytest.fixture
def health_report(setup, monkeypatch):
    from app.services import docx_export
    from app.routers import samples
    client, app, root, bundle, score, payload = setup
    login(app); app.include_router(samples.router)
    (root/'archive').mkdir(); monkeypatch.setattr(config,'REPORT_OUTPUT_DIR',root/'archive')
    monkeypatch.setattr(docx_export.sample_loader,'load_sample',lambda *a,**k: {'meta':{'LIS_ID':'S1','Test':'WGS'},'variants':{}})
    monkeypatch.setattr(docx_export.sample_loader,'load_sample_secondary_snv',lambda *a,**k: {'variants':{},'categories':{}})
    monkeypatch.setattr(docx_export.sample_loader,'load_sample_pgx',lambda *a,**k: {})
    monkeypatch.setattr(docx_export.sample_loader,'_secondary_panel_genes',lambda: {})
    monkeypatch.setattr(docx_export.report_store,'load',lambda sid: {})
    monkeypatch.setattr(docx_export,'_health_panel_gene_sections',lambda requested: [])
    monkeypatch.setattr(docx_export.acmg_sf_education,'render_acmg_sf_education',lambda *a,**k: None)
    return setup


def test_health_docx_prs_selection_order_methods_rounding_and_omission(health_report):
    from docx import Document
    from app.services import docx_export
    client, app, root, _, _, payload = health_report
    url='/api/samples/S1/health-report.docx?sections=acmg_sf,prs,pgx'
    assert client.get(url).status_code == 400  # Cannot force an unavailable checkbox through the API.
    assert not list((root/'archive').iterdir())
    snapshot=client.post('/api/samples/S1/prs/assessments',json=payload).json()
    response=client.get(url); assert response.status_code == 200, response.text
    doc=Document(io.BytesIO(response.content)); paragraphs=[p.text for p in doc.paragraphs]
    summary=paragraphs.index('  綜合說明:')
    assert paragraphs[summary+1] == 'ACMG疾病風險基因'
    assert paragraphs.index('ACMG疾病風險基因') < paragraphs.index('多基因風險評估') < paragraphs.index('藥物基因體學')
    assert '多基因風險評估' in next(p for p in paragraphs if p.startswith('二、檢驗套組'))
    methods=paragraphs[paragraphs.index('四、檢測方法說明'):]
    assert next(i for i,p in enumerate(methods) if 'PGS004879' in p) < next(i for i,p in enumerate(methods) if '某些藥物基因' in p)
    text='\n'.join(paragraphs+[cell.text for table in doc.tables for row in table.rows for cell in row.cells])
    assert [cell.text for cell in doc.tables[0].rows[0].cells] == ['疾病','PRS 百分位','未來10年風險','未來20年風險','至80歲累積風險']
    assert '第97.0百分位' in text
    expected=docx_export._health_prs_number(snapshot['result']['risk_10_year']['risk'],risk=True)+'%'
    assert expected in text and '不應單獨作為疾病診斷' in text
    expected_twenty=docx_export._health_prs_number(snapshot['result']['risk_20_year']['risk'],risk=True)+'%'
    assert doc.tables[0].rows[1].cells[3].text == expected_twenty
    assert 'Model A' not in text and '原始 PRS score' not in text
    omitted=Document(io.BytesIO(client.get(url.replace('acmg_sf,prs,pgx','acmg_sf,pgx')).content))
    omitted_text='\n'.join(p.text for p in omitted.paragraphs)
    assert '多基因風險評估' not in omitted_text and 'PGS004879' not in omitted_text
    assert docx_export._health_prs_number(60.55) == '60.6'
    assert docx_export._health_prs_number(.0675,risk=True) == '6.8'


def test_health_report_over_70_has_explicit_unavailable_ten_year_horizon(health_report):
    from docx import Document
    client, _, _, _, _, payload = health_report
    client.post('/api/samples/S1/prs/assessments',json={**payload,'age':75}).raise_for_status()
    assert client.get('/api/samples/S1/prs/report-readiness').json()['available']
    response=client.get('/api/samples/S1/health-report.docx?sections=prs')
    assert response.status_code == 200
    text='\n'.join(p.text for p in Document(io.BytesIO(response.content)).paragraphs)
    assert '不適用' in text and '完整十年期間超出' in text
    assert '完整二十年期間超出' in text
    assert '二、檢驗套組: 多基因風險評估' in text
    assert '藥物基因體學' not in text and 'ACMG疾病風險基因' not in text


def test_changed_prs_inputs_during_docx_export_never_archive_a_stale_report(health_report, monkeypatch):
    from app.services import docx_export
    client, _, root, _, _, payload = health_report
    client.post('/api/samples/S1/prs/assessments',json=payload).raise_for_status()
    def change(*a,**k):
        prs.atomic_json(root/'S1/S1.sample_metadata.json',{'lis_id':'S1','changed':True})
    monkeypatch.setattr(docx_export,'_section_health_annotations',change)
    response=client.get('/api/samples/S1/health-report.docx?sections=prs')
    assert response.status_code == 400
    assert not list((root/'archive').iterdir())
