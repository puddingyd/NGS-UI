import copy
import datetime
import hashlib
import json
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
    for age, available in [(70,True),(70.01,False),(79.99,False)]:
        result=client.post('/api/samples/S1/prs/assessments',json={**payload,'age':age}).json()['result']
        assert (result['risk_10_year'] is not None) == available
        assert (result['reference']['risk_10_year'] is not None) == available
        assert result['risk_to_80']['risk'] > 0


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
