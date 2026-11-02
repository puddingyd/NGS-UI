"""Queue-to-publish integration using synthetic gVCF and isolated fake tools.

Actual WGS/scoring resources are never used or made available by these tests.
"""
import gzip
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from app import config
from app.services import prs, prs_jobs as jobs, sample_layout, dragen_jobs
from app.workers import prs_run, dragen_run


@pytest.fixture
def env(tmp_path, monkeypatch):
    monkeypatch.setattr(config, 'DATA_ROOT', tmp_path/'data')
    monkeypatch.setattr(dragen_run, 'TERTIARY_JOBS_DIR', tmp_path/'tertiary_jobs')
    monkeypatch.setattr(sample_layout, 'state_dir', lambda sid, **kw: tmp_path/'cases'/sid)
    monkeypatch.setattr(sample_layout, 'state_file', lambda sid,name,**kw: tmp_path/'cases'/sid/f'{sid}.{name}')
    monkeypatch.setattr(dragen_jobs, 'active_sample_ids', jobs.active_ids)
    monkeypatch.setattr(jobs.subprocess, 'Popen', lambda *a,**kw: SimpleNamespace(pid=os.getpid()))
    monkeypatch.setattr(prs, 'SCORE_DEFINITION', dict(model='A',model_lock_sha256='0'*64))
    anchor=tmp_path/'run/vcf.gz/SOURCE.hard-filtered.vcf.gz'
    anchor.parent.mkdir(parents=True);anchor.write_bytes(b'anchor')
    gvcf=anchor.with_name('SOURCE.hard-filtered.gvcf.gz')
    with gzip.open(gvcf,'wt') as f:
        f.write('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n##ALT=<ID=NON_REF,Description="reference">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSOURCE\n')
    Path(str(gvcf)+'.tbi').write_bytes(b'synthetic-index')
    marker=dict(source_sample_id='SOURCE',source_vcf_path=str(anchor),pipeline_type='dragen',seq_type='WGS')
    prs.atomic_json(sample_layout.state_file('CASE','pipeline_source.json'),marker)
    cfg=dict(project=str(tmp_path/'project'))
    for key in ('raw_script','original_scorefile','fasta','dbsnp','plink_img','raw_nextflow_config','ancestry_nextflow_config','reference_panel','method2_lock'):
        p=tmp_path/key;p.write_text('synthetic');cfg[key]=str(p)
    monkeypatch.setattr(jobs,'settings',lambda:(cfg,'cfg-hash'))
    summary=json.loads((jobs.SCORING.parent/'examples/synthetic_summary.json').read_text())
    summary['sample']='SOURCE';summary['ancestry_status']={'status':'NOT_ASSESSED'}
    calls=[]
    def command(worker,args,step):
        calls.append(args)
        run=args[args.index('--run-id')+1]
        dest=Path(cfg['project'])/'integrated_runs'/run/'SOURCE'
        if '--gvcf' in args:
            prs.atomic_json(dest/'method2/summary.json',summary)
        else:
            prs.atomic_json(dest/'qc/qc.json',dict(sample='SOURCE',qc_status='PASS',ancestry='EAS',ancestry_low_confidence=False))
    monkeypatch.setattr(prs_run.Worker,'command',command)
    return SimpleNamespace(root=tmp_path,anchor=anchor,gvcf=gvcf,cfg=cfg,summary=summary,calls=calls,marker=marker)


def test_tertiary_to_score_and_ancestry_without_clinical_metadata(env):
    output=jobs.enqueue_tertiary([{'sample_id':'CASE'}],'tertiary-parent',lambda text:None)
    jid=output['CASE']['job_id']
    assert jobs.read_job(jid)['parent_job_id']=='tertiary-parent'
    assert not sample_layout.state_file('CASE','sample_metadata.json').exists()
    prs_run.Worker(jid).run()
    record,_=prs.read_json(prs.score_path('CASE'))
    assert record['summary']['ancestry_status']['ancestry']=='EAS'
    assert record['ancestry_qc_sha256']
    assert record['summary']['scores']==env.summary['scores']
    assert jobs.read_job(jid)['status']=='done'
    assert not prs.assessment_dir('CASE').exists()  # No clinical risk automatically.
    assert len(env.calls)==2 and '--resume' in env.calls[0]
    assert '--with-legacy-ancestry' not in env.calls[0]
    assert '--raw-source' in env.calls[1]
    run=record['pipeline_run_id']
    original=Path(env.cfg['project'])/'integrated_runs'/run/'SOURCE/method2/summary.json'
    assert json.loads(original.read_text())['ancestry_status']=={'status':'NOT_ASSESSED'}


def test_ancestry_failure_keeps_percentile_and_retry_reuses_raw_run(env,monkeypatch):
    original=prs_run.Worker.command
    def fail_ancestry(worker,args,step):
        if '--raw-source' in args: raise RuntimeError('fake ancestry failure')
        original(worker,args,step)
    monkeypatch.setattr(prs_run.Worker,'command',fail_ancestry)
    jid=jobs.start('CASE')['job_id'];prs_run.Worker(jid).run()
    assert jobs.read_job(jid)['status']=='done_with_warnings'
    score,_=prs.load_score('CASE')
    assert prs.calculator.score_row(score,prs.SCORE_DEFINITION,False)['eas_percentile']==50
    with pytest.raises(ValueError): prs.calculator.score_row(score,prs.SCORE_DEFINITION)
    retry=jobs.start('CASE')['job_id']
    assert jobs.read_job(retry)['pipeline_run_id']==jobs.read_job(jid)['pipeline_run_id']


def test_duplicate_cancel_and_cross_sample_access(env):
    jid=jobs.start('CASE')['job_id']
    with pytest.raises(ValueError):jobs.start('CASE')
    with pytest.raises(FileNotFoundError):jobs.cancel('OTHER',jid)
    with pytest.raises(FileNotFoundError):jobs.log('OTHER',jid)
    jobs.cancel('CASE',jid)
    worker=prs_run.Worker(jid)
    with pytest.raises(prs_run.Cancelled):worker.run()
    assert not env.calls and not prs.score_path('CASE').exists()


def test_input_change_rejects_publish_and_new_run_does_not_reuse(env):
    jid=jobs.start('CASE')['job_id']
    env.anchor.write_bytes(b'changed anchor')
    with pytest.raises(ValueError,match='來源已更新'):prs_run.Worker(jid).run()
    worker=prs_run.Worker(jid);worker.update(status='failed')
    next_id=jobs.start('CASE')['job_id']
    assert jobs.read_job(next_id)['pipeline_run_id']!=jobs.read_job(jid)['pipeline_run_id']
    assert not prs.score_path('CASE').exists()


def test_changed_gvcf_invalidates_published_score(env):
    jid=jobs.start('CASE')['job_id'];prs_run.Worker(jid).run()
    Path(str(env.gvcf)+'.tbi').write_bytes(b'new-index')
    with pytest.raises(prs.NotReady):prs.load_score('CASE')


@pytest.mark.parametrize('problem',['wrong_sample','multi_sample','not_gvcf','no_index','wrong_build'])
def test_gvcf_header_and_index_rejection(env,problem):
    if problem=='no_index':Path(str(env.gvcf)+'.tbi').unlink()
    else:
        names={'wrong_sample':'OTHER','multi_sample':'SOURCE\tOTHER','not_gvcf':'SOURCE','wrong_build':'SOURCE'}
        with gzip.open(env.gvcf,'wt') as f:
            f.write('##contig=<ID=chr1,length='+('249250621' if problem=='wrong_build' else '248956422')+'>\n')
            if problem!='not_gvcf':f.write('##ALT=<ID=NON_REF,Description="reference">\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+names[problem]+'\n')
    with pytest.raises(ValueError):jobs.start('CASE')
    assert not jobs.jobs()


def test_ambiguous_gvcf_requires_selection_and_rejects_arbitrary_path(env):
    other=env.gvcf.with_name('SOURCE.gvcf.gz');other.write_bytes(env.gvcf.read_bytes())
    Path(str(other)+'.tbi').write_bytes(b'index')
    choices=jobs.candidates(jobs.source('CASE'),env.cfg)
    assert len(choices)==2
    with pytest.raises(ValueError):jobs.start('CASE')
    with pytest.raises(ValueError):jobs.start('CASE','/etc/passwd')
    jid=jobs.start('CASE',choices[0]['id'])['job_id']
    assert jobs.read_job(jid)['plan']['gvcf']==choices[0]['path']


def test_missing_setup_isolated_from_tertiary_and_wes_rejected(env,monkeypatch):
    env.marker['seq_type']='WES'
    prs.atomic_json(sample_layout.state_file('CASE','pipeline_source.json'),env.marker)
    result=jobs.enqueue_tertiary([{'sample_id':'CASE'}],'parent',lambda text:None)
    assert result['CASE']['status']=='blocked'
    assert not jobs.jobs()
    assert 'WGS' in result['CASE']['error']


def test_server_launch_argv_is_fixed(env,monkeypatch):
    captured=[]
    def spawn(args,**kw):captured.append((args,kw));return SimpleNamespace(pid=os.getpid())
    monkeypatch.setattr(jobs.subprocess,'Popen',spawn)
    job=jobs.start('CASE')
    assert captured[0][0][-2:]==['app.workers.prs_run',job['job_id']]
    assert captured[0][1]['start_new_session'] is True
    assert 'shell' not in captured[0][1]


def test_explicit_new_run_keeps_old_raw_run(env):
    jid=jobs.start('CASE')['job_id'];prs_run.Worker(jid).run()
    next_id=jobs.start('CASE',new_run=True)['job_id']
    assert jobs.read_job(next_id)['pipeline_run_id']!=jobs.read_job(jid)['pipeline_run_id']
    assert prs.score_path('CASE').exists()


def test_tertiary_checkbox_reaches_worker(env,monkeypatch):
    captured=[]
    monkeypatch.setattr(dragen_jobs,'TERTIARY_JOBS_DIR',env.root/'tertiary_jobs')
    monkeypatch.setattr(dragen_jobs.subprocess,'Popen',lambda args,**kw:(captured.append(args) or SimpleNamespace(pid=os.getpid())))
    jid=dragen_jobs.start_job(str(env.anchor),'CASE',source_sample_id='SOURCE',mode='dragen',seq_type='WGS',with_prs=True)
    assert '--with-prs' in captured[0]
    assert dragen_jobs.load_state(jid)['with_prs'] is True


def test_cancel_stops_owned_process_group(env,monkeypatch):
    jid=jobs.start('CASE')['job_id'];worker=prs_run.Worker(jid)
    stopped=[]
    process=SimpleNamespace(pid=12345,poll=lambda:None,wait=lambda **kw:0)
    monkeypatch.setattr(prs_run.subprocess,'Popen',lambda *a,**kw:process)
    monkeypatch.setattr(prs_run.os,'killpg',lambda pid,sig:stopped.append((pid,sig)))
    # Exercise the real command monitor instead of the fixture's fake tools.
    # Original method saved separately below; no external executable is launched.
    jobs.cancel('CASE',jid)
    with pytest.raises(prs_run.Cancelled): REAL_COMMAND(worker,['fake-tool'],'test')
    assert stopped and stopped[0][0]==12345


REAL_COMMAND = prs_run.Worker.command
