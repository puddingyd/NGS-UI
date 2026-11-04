import os
from pathlib import Path
import sys
import time

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient
from starlette.middleware.sessions import SessionMiddleware

from app import config
from app.auth import current_user
from app.routers import dragen as router
from app.services import dragen_jobs, prs, prs_jobs, postprocessing_log, container_runtime
from app.workers import dragen_run


@pytest.mark.parametrize('parent_exists', [True, False])
def test_bind_cleanup_is_by_destination_and_drops_only_managed_job_sources(tmp_path, parent_exists):
    root = tmp_path/'nf_work'
    parent = root/'contexts/shared-dragen/.runtime/jobs/old/tmp'
    if parent_exists: parent.mkdir(parents=True)
    raw = f'/reference:/ref:ro,{parent}:{parent},{parent}:/tmp,/admin:/tmp,/data:/data'
    env = container_runtime.job_bind_environment(dict(APPTAINER_BIND=raw, SINGULARITY_BINDPATH=raw), root,
                                                ['/case/tmp:/case/tmp', '/case/tmp:/tmp', '/case/tmp:/var/tmp'])
    value = env['APPTAINER_BINDPATH']
    assert str(parent) not in value
    assert '/reference:/ref:ro' in value and '/data:/data' in value
    assert env['APPTAINER_BIND'] == env['SINGULARITY_BIND'] == ''
    destinations = [v.split(':')[1] for v in value.split(',')]
    assert len(destinations) == len(set(destinations))
    assert value.count(':/tmp') == 1


def test_raw_annotsv_output_is_saved_but_summary_keeps_timing_and_warning(tmp_path, capsys):
    raw = tmp_path/'postprocessing.log'
    code = ('print("[post-processing-step] cnv-rescue start"); '
            'print("AnnotSV huge configuration dump\\n"*80); '
            'print("WARNING: optional annotations unavailable"); '
            'print("[cnv-rescue] integrated_PASS=2 B=1 annotated=3"); '
            'print("[post-processing-step] cnv-rescue done")')
    dragen_run._run([sys.executable, '-c', code], label='CNV annotation', post_sample='S1',
                    post_step='cnv-rescue', raw_log=raw)
    summary = capsys.readouterr().out
    assert 'huge configuration dump' not in summary
    assert 'huge configuration dump' in raw.read_text()
    assert summary.count('step=cnv-rescue status=start') == 1
    assert 'step=cnv-rescue status=done elapsed=' in summary
    assert 'optional annotations unavailable' in summary and 'annotated=3' in summary


def test_failed_tool_retains_raw_diagnostics_and_surfaces_actual_cause(tmp_path, capsys):
    raw = tmp_path/'postprocessing.log'
    with pytest.raises(RuntimeError, match='synthetic missing input'):
        dragen_run._run([sys.executable, '-c',
            'print("[post-processing-step] spliceai start"); print("verbose dump"); '
            'print("FATAL: synthetic missing input"); raise SystemExit(3)'],
            label='postprocessing S1', post_sample='S1', raw_log=raw)
    log = capsys.readouterr().out
    assert 'step=spliceai status=failed elapsed=' in log
    assert 'verbose dump' not in log and 'verbose dump' in raw.read_text()
    assert 'FATAL: synthetic missing input' in log


def test_summary_reports_gpn_progress_without_large_signature_dump():
    events = []
    summary = postprocessing_log.Summary('S1', 'postprocessing', events.append)
    summary.feed('[perf] snv_review.gpn_start elapsed=0s provider=gpn_star signature_diff=BIG_DUMP')
    summary.feed('[perf] snv_review.gpn_done elapsed=12.5s provider=gpn_star status=complete')
    assert len(events) == 2
    assert 'step=gpn_star status=done elapsed=12.5s' in events[-1]
    assert 'BIG_DUMP' not in ''.join(events)


def test_parent_log_and_api_surface_prs_failure_without_failing_tertiary(tmp_path, monkeypatch):
    jobs = tmp_path/'tertiary_jobs'
    monkeypatch.setattr(dragen_jobs, 'TERTIARY_JOBS_DIR', jobs)
    monkeypatch.setattr(dragen_jobs, 'LEGACY_DRAGEN_JOBS_DIR', tmp_path/'legacy')
    monkeypatch.setattr(config, 'DATA_ROOT', tmp_path/'data')
    jid = 'a'*32
    parent = dict(job_id='parent-1', state='done', sample_id='S1',
                  prs_jobs=dict(S1=dict(job_id=jid, status='queued')))
    dragen_jobs.save_state('parent-1', parent)
    (jobs/'parent-1/log.txt').write_text('tertiary finished\n')
    raw = jobs/'parent-1/postprocessing.log'; raw.write_text('full AnnotSV output\n')
    child = dict(job_id=jid, sample_id='S1', parent_job_id='parent-1', status='failed',
                 created=time.time(), updated=time.time(), step='容器暫存檢查', error='exit 255')
    prs.atomic_json(prs_jobs.job_dir(jid)/'state.json', child)
    (prs_jobs.job_dir(jid)/'log.txt').write_text('FATAL: parent tmp does not exist\nTraceback:\nRuntimeError: exit 255\n')
    view = dragen_jobs.read_log('parent-1')
    assert 'status=failed' in view and 'parent tmp does not exist' in view
    assert dragen_jobs.load_state('parent-1')['state'] == 'done'
    app = FastAPI(); app.add_middleware(SessionMiddleware, secret_key='test')
    app.include_router(router.router); client = TestClient(app)
    assert client.get('/api/dragen/jobs/parent-1/postprocessing-log').status_code == 401
    app.dependency_overrides[current_user] = lambda: dict(username='tester')
    state = client.get('/api/dragen/jobs/parent-1').json()
    assert state['state'] == 'done' and state['prs_active'] is False
    assert 'parent tmp does not exist' in state['log_tail']
    assert client.get('/api/dragen/jobs/parent-1/postprocessing-log').text == raw.read_text()
    child.update(status='running', error=None, pid=os.getpid(), updated=time.time())
    prs.atomic_json(prs_jobs.job_dir(jid)/'state.json', child)
    assert client.get('/api/dragen/jobs/parent-1').json()['prs_active'] is True
    assert dragen_jobs.list_jobs()[0]['prs_active'] is True


def test_log_download_rejects_invalid_job_path():
    with pytest.raises(ValueError): dragen_jobs.postprocessing_log_path('..')
