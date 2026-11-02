#!/usr/bin/env python3
"""Execution-only v2.1 add-on: locks, guarded resume, fingerprints, batch summary.

Requires installation beside the unchanged pipeline.py, pipeline_v2.py and
method2_score.py from the validated integrated v2 release.
"""
import argparse
import contextlib
import csv
import datetime as dt
import fcntl
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path

import pipeline as v1
import pipeline_v2 as v2
import method2_score as m2
import v21_raw_resume as repair

ROOT=Path(__file__).resolve().parent
BASELINE=ROOT/'BASELINE_V2_SHA256SUMS.txt'
FP_SCHEMA='PGS004879_V21_ORCHESTRATION_20261031'
FIELDS=['sample','execution_status','raw_status','method2_status',
        'model','n_loci','observed_loci','missing_loci','observed_sum',
        'imputed_sum','total_sum','eas_z','eas_percentile','intervene_group',
        'method2_lock_sha256','report','message']


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')


def safe_id(s):
    if s in ('.','..') or not v1.SID.fullmatch(s):
        raise RuntimeError('Unsafe sample/run ID: '+str(s))
    return s


def baseline():
    if not BASELINE.is_file():
        raise RuntimeError('Missing BASELINE_V2_SHA256SUMS.txt; install whole add-on')
    checked=[]
    for line in BASELINE.read_text().splitlines():
        if not line.strip():continue
        digest,rel=line.split(maxsplit=1)
        rel=rel.strip()
        path=(ROOT/rel).resolve()
        if not path.is_relative_to(ROOT.resolve()) or not path.is_file():
            raise RuntimeError('Missing/unsafe validated v2 component: '+rel)
        if m2.sha256(path)!=digest:
            raise RuntimeError('Validated v2 file changed: '+rel)
        checked.append(rel)
    if len(checked)!=10:
        raise RuntimeError('Incomplete baseline component manifest')
    return checked


@contextlib.contextmanager
def guard_sample(project,sample):
    """Separate orchestrator lock; v2 itself keeps its original v1 sample lock."""
    lockdir=Path(project)/'integrated_runs'/'.v21_locks'
    lockdir.mkdir(parents=True,exist_ok=True)
    with open(lockdir/f'{sample}.lock','a+') as handle:
        try:
            fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f'Sample {sample} is already running under v2.1')
        handle.seek(0);handle.truncate()
        handle.write(f'pid={os.getpid()} start={utc()}\n');handle.flush()
        try:yield
        finally:fcntl.flock(handle,fcntl.LOCK_UN)


def fingerprint(c,args):
    base=v1.fingerprint(c,args.sample,args.mode,args.gvcf,args.raw_source)
    return {'schema':FP_SCHEMA,'v1_inputs':base,
            'method2_lock_sha256':m2.sha256(c['method2_lock']),
            'model_dir':str(Path(c['method2_model_dir']).resolve()),
            'v2_core_sha256':{p:m2.sha256(ROOT/p) for p in
                              ('pipeline.py','pipeline_v2.py','method2_score.py')},
            'v21_code_sha256':{p:m2.sha256(ROOT/p) for p in
                               ('pipeline_v21.py','v21_raw_resume.py')},
            'raw_source':str(Path(args.raw_source).resolve()) if args.raw_source else None,
            'with_legacy_ancestry':args.with_legacy_ancestry}


def completed_outputs(dest,sample,expected_lock):
    p=Path(dest)/'method2'
    status=p/'status.json'
    summary=p/'summary.json'
    scores=p/'scores.tsv'
    missing=p/'missing_loci.tsv.gz'
    report=p/'report.md'
    if not all(x.is_file() for x in (status,summary,scores,missing,report,p/'inputs.lock.json')):
        return False
    try:
        if json.loads(status.read_text()).get('status')!='PASS':return False
        d=json.loads(summary.read_text())
        if (d.get('sample')!=sample or d.get('stage')!='METHOD2_TECHNICAL_PASS' or
            d.get('model_lock_sha256')!=expected_lock):return False
        with scores.open(newline='') as stream:
            rows=list(csv.DictReader(stream,delimiter='\t'))
        if len(rows)!=2 or {r['model'] for r in rows}!={'A','B'} or \
            any(r['sample']!=sample or int(r['observed_loci'])+int(r['missing_loci'])!=int(r['n_loci']) for r in rows):
            return False
        return True
    except (ValueError,KeyError,TypeError,json.JSONDecodeError):
        return False


def output_fingerprints(dest):
    root=Path(dest)
    paths=[root/'method2'/p for p in (
        'status.json','summary.json','scores.tsv','missing_loci.tsv.gz',
        'report.md','inputs.lock.json')]
    paths += [root/'raw/audit/final_reconciliation.json',
              root/'raw/audit/all_score_locus_decisions.tsv']
    result={}
    for path in paths:
        if path.is_file():result[str(path.relative_to(root))]=m2.sha256(path)
    return result


def write_state(dest,stage,message=None):
    d={'schema':FP_SCHEMA,'stage':stage,'updated_utc':utc()}
    if message:d['message']=str(message)
    m2.json_write_atomic(Path(dest)/'orchestration/v21_status.json',d)



def verify_finished_output_without_manifest(dest,sample,args,c):
    """Recover crash after v2 PASS but before v2.1 output digest was written.

    Independently regenerate Method2 into disposable *new* staging and check
    its five deterministic output files against the existing files. Never
    overwrite or delete prior raw/Method2 data.
    """
    dest=Path(dest)
    if args.mode=='fresh':
        raw=dest/'raw'
        v2.verify_fresh_raw_origin(raw,args.gvcf,c['original_scorefile'])
    else:
        raw=Path(args.raw_source)
    v1.raw_integrity(raw,sample)
    qc=dest/'qc/qc.json'
    qc=qc if qc.is_file() else None
    compare=('summary.json','scores.tsv','missing_loci.tsv.gz','report.md')
    with tempfile.TemporaryDirectory(prefix='.v21_adoption_check_',dir=dest/'orchestration') as tmp:
        m2.score(sample,raw,c['method2_model_dir'],c['method2_lock'],
                 c['original_scorefile'],tmp,qc)
        for name in compare:
            before=dest/'method2'/name
            regenerated=Path(tmp)/name
            if not before.is_file() or m2.sha256(before)!=m2.sha256(regenerated):
                raise RuntimeError('Completed pre-manifest output differs from independent re-score: '+name)
    return True


def run_one(args):
    baseline()
    c=v1.cfg_read(args.config)
    sample=safe_id(args.sample)
    run_id=safe_id(args.run_id)
    if args.mode=='fresh':
        if not args.gvcf or args.raw_source:
            raise RuntimeError('Fresh requires --gvcf and forbids --raw-source')
        args.gvcf=str(Path(args.gvcf).resolve())
        if not Path(args.gvcf).is_file():raise RuntimeError('gVCF not found')
    else:
        if args.gvcf or not args.raw_source:
            raise RuntimeError('Reuse requires --raw-source and forbids --gvcf')
        args.raw_source=str(Path(args.raw_source).resolve())
    m2.validate_locked_model(c['method2_model_dir'],c['method2_lock'],c['original_scorefile'])
    dest=Path(c['project'])/'integrated_runs'/run_id/sample
    with guard_sample(c['project'],sample):
        stage_dir=dest/'orchestration'
        fp_path=stage_dir/'v21_inputs.lock.json'
        fp=fingerprint(c,args)
        if fp_path.exists():
            if not args.resume:
                raise RuntimeError('Existing v2.1 run requires explicit --resume; do not overwrite')
            if json.loads(fp_path.read_text())!=fp:
                raise RuntimeError('Immutable input fingerprint changed. Use NEW --run-id')
        else:
            if dest.exists():
                raise RuntimeError('Existing run without v2.1 input manifest: refusing automatic adoption')
            stage_dir.mkdir(parents=True,exist_ok=False)
            m2.json_write_atomic(fp_path,fp)
        try:
            sha=m2.sha256(c['method2_lock'])
            is_complete=completed_outputs(dest,sample,sha)
            digest_path=stage_dir/'v21_output_sha256.json'
            if is_complete:
                # Integrity must be checked before returning a prior PASS.
                if not digest_path.is_file():
                    if not args.resume:
                        raise RuntimeError('Completed v2 output missing v2.1 digest; explicit --resume required')
                    verify_finished_output_without_manifest(dest,sample,args,c)
                    m2.json_write_atomic(digest_path,output_fingerprints(dest))
                    print('V21_COMPLETED_STAGE_RECOVERED_AND_SEALED',flush=True)
                if json.loads(digest_path.read_text())!=output_fingerprints(dest):
                    raise RuntimeError('Completed output content changed; use a new run ID')
                if args.mode=='fresh':
                    v2.verify_fresh_raw_origin(dest/'raw',args.gvcf,c['original_scorefile'])
                write_state(dest,'PASS','Unchanged verified prior results; no scoring rerun')
                print(f'V21_ALREADY_COMPLETE: {sample}',flush=True)
                return
            if digest_path.exists():
                raise RuntimeError('Previously completed outputs now invalid; refusing overwrite')
            if args.mode=='fresh':
                raw=dest/'raw'
                if raw.exists():
                    try:
                        v1.raw_integrity(raw,sample)
                        v2.verify_fresh_raw_origin(raw,args.gvcf,c['original_scorefile'])
                        write_state(dest,'RAW_VALID','Existing complete raw stage validated')
                    except RuntimeError:
                        if not args.resume:
                            raise RuntimeError('Raw incomplete; explicit --resume required')
                        write_state(dest,'RAW_REPAIR','Validate roundtrip then resume same Nextflow work')
                        with v1.sample_lock(c['project'],sample):
                            repair.resume_raw_after_roundtrip(raw,sample,args.gvcf,c)
                        write_state(dest,'RAW_VALID','Validated Nextflow recovery')
            cmd=[sys.executable,str(ROOT/'pipeline_v2.py'),'run',
                 '--config',str(Path(args.config).resolve()),
                 '--sample',sample,'--run-id',run_id,'--mode',args.mode]
            if args.mode=='fresh':cmd+=['--gvcf',args.gvcf]
            else:cmd+=['--raw-source',args.raw_source]
            if args.with_legacy_ancestry:cmd+=['--with-legacy-ancestry']
            write_state(dest,'V2_RUNNING')
            v1.log_run(cmd,dest/'logs/v21_v2_stage.log')
            if not completed_outputs(dest,sample,sha):
                raise RuntimeError('v2 ended without validated Method 2 results')
            m2.json_write_atomic(digest_path,output_fingerprints(dest))
            write_state(dest,'PASS')
            print(f'V21 SAMPLE COMPLETE: {sample}',flush=True)
        except Exception as exc:
            write_state(dest,'FAILED',str(exc))
            raise


def summarize(project,run_id,write=True):
    safe_id(run_id)
    root=Path(project)/'integrated_runs'/run_id
    if not root.is_dir():raise RuntimeError('Run directory missing: '+str(root))
    result=[]
    for child in sorted(root.iterdir()):
        if not child.is_dir() or child.name.startswith('.'):continue
        if not (child/'orchestration/v21_inputs.lock.json').is_file():continue
        sid=child.name
        if not v1.SID.fullmatch(sid):continue
        state=child/'orchestration/v21_status.json'
        status='UNKNOWN';message=''
        if state.is_file():
            data=json.loads(state.read_text())
            status=data.get('stage','UNKNOWN')
            message=data.get('message','')
        mp=child/'method2'
        summary=mp/'summary.json'
        rows=[]
        summary_data={}
        if summary.is_file():
            try:
                summary_data=json.loads(summary.read_text())
                if summary_data.get('sample')==sid:
                    rows=summary_data.get('scores',[])
            except (ValueError,TypeError):pass
        digest_path=child/'orchestration/v21_output_sha256.json'
        if status=='PASS':
            if not digest_path.is_file():
                status='INVALID_OUTPUT';message='Missing output SHA256 manifest'
            elif json.loads(digest_path.read_text())!=output_fingerprints(child):
                status='INVALID_OUTPUT';message='Saved outputs changed after PASS'
        if not rows:rows=[{'model':''}]
        for row in rows:
            result.append({'sample':sid,'execution_status':status,
                           'raw_status':'PASS' if (child/'raw/audit/final_reconciliation.json').is_file() else 'NOT_VERIFIED',
                           'method2_status':'PASS' if status=='PASS' and (mp/'status.json').is_file() else 'INCOMPLETE',
                           'model':row.get('model',''),
                           **{key:row.get(key,'') for key in (
                               'n_loci','observed_loci','missing_loci','observed_sum',
                               'imputed_sum','total_sum','eas_z','eas_percentile','intervene_group')},
                           'method2_lock_sha256':summary_data.get('model_lock_sha256','') if rows!=[{'model':''}] else '',
                           'report':str(mp/'report.md') if (mp/'report.md').is_file() else '',
                           'message':message})
    if write:
        target=root/'batch_summary_v21.tsv'
        tmp=target.with_name(target.name+f'.tmp.{os.getpid()}')
        with tmp.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',extrasaction='ignore')
            w.writeheader();w.writerows(result)
        os.replace(tmp,target)
        print('BATCH SUMMARY:',target,flush=True)
    return result


def batch(a):
    config=v1.cfg_read(a.config)
    path=Path(a.manifest).resolve()
    if not path.is_file():raise RuntimeError('Batch manifest does not exist')
    with path.open(newline='') as f:
        rd=csv.DictReader(f,delimiter='\t')
        if not {'sample','mode','gvcf','raw_source'}<=set(rd.fieldnames or []):
            raise RuntimeError('Manifest requires sample,mode,gvcf,raw_source')
        rows=list(rd)
    if not rows:raise RuntimeError('Empty batch manifest')
    seen=set();jobs=[]
    for row in rows:
        sid=safe_id(row['sample']);mode=row['mode']
        if sid in seen:raise RuntimeError('Duplicate batch sample ID '+sid)
        seen.add(sid)
        if mode not in ('fresh','reuse'):raise RuntimeError('Invalid batch mode for '+sid)
        jobs.append(argparse.Namespace(config=a.config,sample=sid,run_id=a.run_id,
                                       mode=mode,gvcf=row['gvcf'] or None,
                                       raw_source=row['raw_source'] or None,
                                       with_legacy_ancestry=a.with_legacy_ancestry,
                                       resume=a.resume))
    failures=[]
    # Default serial; concurrent runs are opt-in because each can use Nextflow.
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        work={pool.submit(run_one,x):x.sample for x in jobs}
        for future in as_completed(work):
            sid=work[future]
            try:future.result()
            except Exception as e:
                failures.append((sid,str(e)))
                print(f'BATCH FAILURE {sid}: {e}',file=sys.stderr,flush=True)
    records=summarize(config['project'],a.run_id)
    expected={r.sample for r in jobs}
    passed={r['sample'] for r in records if r['execution_status']=='PASS' and r['model'] in ('A','B')}
    if failures or passed!=expected:
        raise RuntimeError(f'Batch unfinished: failures={failures}, passed={len(passed)}/{len(expected)}')
    print(f'V21 BATCH PASS: {len(expected)} samples',flush=True)



def verify_pilot(project,run_id):
    """Validate the 6 frozen Method 2 A/B outputs from the three known cases."""
    root=Path(project)/'integrated_runs'/safe_id(run_id)
    golden_counts={
        '26T00970':{'A':(609032,1273),'B':(602628,810)},
        '26T00969':{'A':(609162,1143),'B':(602768,670)},
        '26T00968':{'A':(609272,1033),'B':(602869,569)},
    }
    results=[]
    for sid,models in golden_counts.items():
        dest=root/sid
        state=dest/'orchestration/v21_status.json'
        digests=dest/'orchestration/v21_output_sha256.json'
        if not (state.is_file() and digests.is_file()):
            raise RuntimeError(f'Pilot run not complete: {sid}')
        if json.loads(state.read_text()).get('stage')!='PASS' or json.loads(digests.read_text())!=output_fingerprints(dest):
            raise RuntimeError(f'Pilot output invalid or changed: {sid}')
        scores=dest/'method2/scores.tsv'
        with scores.open(newline='') as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        if len(rows)!=2 or {r['model'] for r in rows}!={'A','B'}:
            raise RuntimeError(f'Pilot models missing: {sid}')
        for row in rows:
            model=row['model']
            e_sum,e_pct=v2.PILOT_TARGET[sid][model]
            e_obs,e_missing=models[model]
            diff_sum=abs(float(row['total_sum'])-e_sum)
            diff_pct=abs(float(row['eas_percentile'])-e_pct)
            obs=int(row['observed_loci']); missing=int(row['missing_loci'])
            ok=(row['sample']==sid and diff_sum<=2e-6 and diff_pct<=1e-6 and
                (obs,missing)==(e_obs,e_missing))
            results.append((sid,model,str(obs),str(missing),
                            f'{diff_sum:.12g}',f'{diff_pct:.12g}',
                            'PASS' if ok else 'FAIL'))
    out=root/'pilot_regression_v21.tsv'
    temp=out.with_name(out.name+f'.tmp.{os.getpid()}')
    with temp.open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t')
        w.writerow(['sample','model','observed','missing','sum_abs_diff','percentile_abs_diff','status'])
        w.writerows(results)
    os.replace(temp,out)
    if len(results)!=6 or any(r[-1]!='PASS' for r in results):
        raise RuntimeError(f'PILOT REGRESSION FAILED: {out}')
    print(f'V21 PILOT REGRESSION: PASS (6/6); {out}')
    return results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    run=sub.add_parser('run',help='Safe per-sample execution and restricted post-roundtrip resume')
    run.add_argument('--config',required=True)
    run.add_argument('--sample',required=True)
    run.add_argument('--run-id',required=True)
    run.add_argument('--mode',choices=('fresh','reuse'),default='fresh')
    run.add_argument('--gvcf')
    run.add_argument('--raw-source')
    run.add_argument('--resume',action='store_true')
    run.add_argument('--with-legacy-ancestry',action='store_true')
    b=sub.add_parser('batch',help='Run tab-separated sample manifest; continue other samples on failure')
    b.add_argument('--config',required=True)
    b.add_argument('--manifest',required=True)
    b.add_argument('--run-id',required=True)
    b.add_argument('--jobs',type=int,default=1)
    b.add_argument('--resume',action='store_true')
    b.add_argument('--with-legacy-ancestry',action='store_true')
    pilot=sub.add_parser('verify-pilot',help='Compare three fresh/reuse samples to six frozen benchmarks')
    pilot.add_argument('--config',required=True)
    pilot.add_argument('--run-id',required=True)
    sm=sub.add_parser('summarize',help='Summarize completed or failed runs without execution')
    sm.add_argument('--config',required=True)
    sm.add_argument('--run-id',required=True)
    check=sub.add_parser('verify-baseline',help='Check immutable v2 program components')
    check.add_argument('--config')
    a=parser.parse_args()
    try:
        if a.command=='verify-baseline':
            paths=baseline()
            print(f'V2 GOLDEN PROGRAM BASELINE: PASS ({len(paths)} files)')
            if a.config:
                c=v1.cfg_read(a.config)
                m2.validate_locked_model(c['method2_model_dir'],c['method2_lock'],c['original_scorefile'])
                print('METHOD 2 MODEL LOCK: PASS')
        elif a.command=='run':run_one(a)
        elif a.command=='batch':
            if a.jobs<1 or a.jobs>8:raise RuntimeError('--jobs must be between 1 and 8')
            batch(a)
        elif a.command=='verify-pilot':
            c=v1.cfg_read(a.config)
            verify_pilot(c['project'],a.run_id)
        else:
            c=v1.cfg_read(a.config)
            records=summarize(c['project'],a.run_id)
            print('SAMPLE ROWS:',len(records))
    except Exception as exc:
        parser.exit(1,'ERROR: '+str(exc)+'\n')


if __name__=='__main__':main()
