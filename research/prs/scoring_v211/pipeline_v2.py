#!/usr/bin/env python3
"""One entrypoint: raw gVCF -> validated v1 raw reconstruction -> fixed Method 2 A/B PRS.

Default raw-only execution does not rerun target-dependent pgsc_calc ancestry.
Optional --with-legacy-ancestry runs the unmodified integrated v1 ancestry/QC path too.
"""
import argparse
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

import pipeline as v1
import method2_score as m2

ROOT=Path(__file__).resolve().parent



def verify_fresh_raw_origin(raw_dir,gvcf,source):
    """Fail closed if a prior completed raw output came from another gVCF."""
    info=Path(raw_dir)/'audit/source_sha256.txt'
    if not info.is_file():m2.die(f'Raw input identity missing: {info}')
    hashes={}
    for line in info.read_text().splitlines():
        fields=line.split(maxsplit=1)
        if len(fields)!=2:m2.die('Malformed raw source SHA256 list')
        hashes[str(Path(fields[1].strip()).resolve())]=fields[0]
    for p in (gvcf,source):
        p=str(Path(p).resolve())
        if p not in hashes:m2.die(f'Raw input path not registered in source manifest: {p}')
        if hashes[p]!=m2.sha256(p):m2.die(f'Raw source checksum changed: {p}')


def run(args):
    c=v1.cfg_read(args.config)
    project=Path(c['project'])
    sample=args.sample
    run_id=args.run_id
    if not (v1.SID.fullmatch(sample) and v1.SID.fullmatch(run_id)) or sample in ('.','..') or run_id in ('.','..'):
        m2.die('Unsafe sample or run ID')
    modeldir=Path(c.get('method2_model_dir','')).resolve()
    lock=Path(c.get('method2_lock','')).resolve()
    if not c.get('method2_model_dir') or not c.get('method2_lock'):
        m2.die('Configure method2_model_dir and method2_lock first')
    # Preflight immutable model before spending hours on a raw run.
    m2.validate_locked_model(modeldir,lock,c['original_scorefile'])
    dest=project/'integrated_runs'/run_id/sample
    if args.mode=='fresh':
        if not args.gvcf: m2.die('Fresh mode requires --gvcf')
        if args.raw_source:m2.die('--raw-source only allowed with reuse')
        gvcf=str(Path(args.gvcf).resolve())
        raw=dest/'raw'
    else:
        if args.gvcf:m2.die('--gvcf only allowed with fresh')
        if not args.raw_source:m2.die('Reuse mode requires explicit --raw-source')
        gvcf=None
        raw=Path(args.raw_source).resolve()
    out=dest/'method2'
    # Unmodified v1 full runner takes its own lock; do not double-lock.
    if args.with_legacy_ancestry:
        cmd=[sys.executable,str(ROOT/'pipeline.py'),'run','--config',str(Path(args.config).resolve()),
             '--sample',sample,'--run-id',run_id,'--mode',args.mode]
        if gvcf:cmd += ['--gvcf',gvcf]
        if args.raw_source:cmd += ['--raw-source',str(raw)]
        if args.recheck_ancestry:cmd += ['--recheck-ancestry']
        v1.log_run(cmd,dest/'logs/v1_full_runner.log')
        v1.raw_integrity(raw,sample)
        if args.mode=='fresh':verify_fresh_raw_origin(raw,gvcf,c['original_scorefile'])
        with v1.sample_lock(project,sample):
            _run_scoring(sample,raw,modeldir,lock,c,dest,out,gvcf,args.mode)
    else:
        if args.recheck_ancestry:m2.die('--recheck-ancestry requires --with-legacy-ancestry')
        with v1.sample_lock(project,sample):
            if args.mode=='fresh':
                # Exact same argument convention as unmodified integrated_v1 pipeline.py.
                v1.check_paths(c,'fresh',gvcf,None,sample)
                try:
                    v1.raw_integrity(raw,sample)
                    print('RAW_ALREADY_VALID; reuse validated raw output',flush=True)
                except RuntimeError:
                    if raw.exists():
                        m2.die('Raw output directory exists but incomplete; use new --run-id, never overwrite raw audit')
                    cmd=['bash',c['raw_script'],'--work',c['raw_work'],'--gvcf',gvcf,
                         '--out',str(raw),'--scorefile',c['original_scorefile'],'--ref',c['fasta'],
                         '--dbsnp',c['dbsnp'],'--plink-img',c['plink_img'],
                         '--container-config',c['raw_nextflow_config'],
                         '--min-dp',str(c['min_dp']),'--min-gq',str(c['min_gq']),
                         '--min-overlap',str(c['min_overlap']),
                         '--palindromic-mode',c['palindromic_mode'],
                         '--threads',str(c['raw_threads']),'--memory-mb',str(c['raw_memory_mb'])]
                    v1.log_run(cmd,dest/'logs/raw_stage.log')
            v1.raw_integrity(raw,sample)
            if args.mode=='fresh':verify_fresh_raw_origin(raw,gvcf,c['original_scorefile'])
            _run_scoring(sample,raw,modeldir,lock,c,dest,out,gvcf,args.mode)


def _run_scoring(sample,raw,modeldir,lock,c,dest,out,gvcf,mode):
    inp={'schema':m2.SCHEMA,'sample':sample,'mode':mode,
         'raw_dir':str(Path(raw).resolve()),
         'raw_reconciliation_sha256':m2.sha256(Path(raw)/'audit/final_reconciliation.json'),
         'raw_locus_audit_sha256':m2.sha256(Path(raw)/'audit/all_score_locus_decisions.tsv'),
         'corrected_scorefile_sha256':m2.sha256(Path(raw)/'scorefiles/PGS004879_hmPOS_GRCh38.txt.gz'),
         'model_lock_sha256':m2.sha256(lock),
         'source_scorefile_sha256':m2.sha256(c['original_scorefile']),
         'gvcf':v1.stat_stamp(gvcf) if gvcf else None,
         'runtime_config':{k:c[k] for k in ('min_dp','min_gq','min_overlap','palindromic_mode')},
         'raw_vendor_script_sha256':m2.sha256(c['raw_script']),
         'v1_runner_sha256':m2.sha256(ROOT/'pipeline.py'),
         'v2_runner_sha256':m2.sha256(ROOT/'pipeline_v2.py'),
         'method2_scoring_code_sha256':m2.sha256(ROOT/'method2_score.py')}
    out.mkdir(parents=True,exist_ok=True)
    path=out/'inputs.lock.json'
    if path.is_file() and json.loads(path.read_text())!=inp:
        m2.die('Method 2 existing input lock differs; use another run ID or restore validated inputs')
    if not path.is_file():m2.json_write_atomic(path,inp)
    try:
        ancestryqc=dest/'qc/qc.json'
        if not ancestryqc.is_file():ancestryqc=None
        result=m2.score(sample,raw,modeldir,lock,c['original_scorefile'],out,ancestryqc)
        for r in result['scores']:
            print(f"MODEL {r['model']}: SUM {r['total_sum']:.12f}, EAS percentile {r['eas_percentile']:.6f}, group {r['intervene_group']}",flush=True)
        print('METHOD2 SCORING: TECHNICAL PASS',flush=True)
        m2.json_write_atomic(out/'status.json',{'status':'PASS','report':str(out/'report.md')})
    except Exception as exc:
        m2.json_write_atomic(out/'status.json',{'status':'FAILED','error':str(exc)})
        raise


PILOT_TARGET = {
    '26T00970': {'A':(0.366424148908984,35.10204081632653), 'B':(0.363472617655232,35.374149659863946)},
    '26T00969': {'A':(0.47483856276928,60.54421768707483), 'B':(0.47580405611060006,61.49659863945578)},
    '26T00968': {'A':(0.556804458044592,79.45578231292517), 'B':(0.55261041849784,79.45578231292517)},
}


def regression(a):
    c=v1.cfg_read(a.config)
    base=Path(c['project'])/'integrated_runs'/a.raw_run
    modeldir=Path(c['method2_model_dir']).resolve()
    lock=Path(c['method2_lock']).resolve()
    m2.validate_locked_model(modeldir,lock,c['original_scorefile'])
    summary=[]
    for sample in PILOT_TARGET:
        raw=base/sample/'raw'
        out=Path(c['project'])/'integrated_runs'/a.out_run/sample/'method2'
        qc=base/sample/'qc/qc.json'
        result=m2.score(sample,raw,modeldir,lock,c['original_scorefile'],out,qc if qc.is_file() else None)
        for row in result['scores']:
            exp_s,exp_p=PILOT_TARGET[sample][row['model']]
            ds=abs(row['total_sum']-exp_s)
            dp=abs(row['eas_percentile']-exp_p)
            if ds>2e-6 or dp>1e-6:
                m2.die(f'REGRESSION FAIL {sample} {row["model"]}: score diff {ds} / percentile diff {dp}')
            summary.append((sample,row['model'],f'{row["total_sum"]:.12f}',f'{row["eas_percentile"]:.9f}','PASS'))
    path=Path(c['project'])/'integrated_runs'/a.out_run/'method2_pilot_regression.tsv'
    import csv
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        wr=csv.writer(f,delimiter='\t')
        wr.writerow(['sample','model','total_sum','eas_percentile','status'])
        wr.writerows(summary)
    print('METHOD2 THREE-SAMPLE REGRESSION: PASS (6/6);',path)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    fr=sub.add_parser('freeze-model',help='One-time model lock, after validated v0 server inventory')
    fr.add_argument('--config',required=True)
    runp=sub.add_parser('run',help='One sample raw gVCF -> Method2 A/B')
    runp.add_argument('--config',required=True)
    runp.add_argument('--sample',required=True)
    runp.add_argument('--run-id',required=True)
    runp.add_argument('--mode',choices=['fresh','reuse'],default='fresh')
    runp.add_argument('--gvcf')
    runp.add_argument('--raw-source')
    runp.add_argument('--with-legacy-ancestry',action='store_true',help='also run unchanged v1 ancestry and QC')
    runp.add_argument('--recheck-ancestry',action='store_true')
    reg=sub.add_parser('regression',help='Validate 6 frozen pilot target scores from completed raw audits; never rerun gVCF')
    reg.add_argument('--config',required=True)
    reg.add_argument('--raw-run',default='pilot_fresh_02')
    reg.add_argument('--out-run',default='method2_v2_regression_01')
    a=p.parse_args()
    try:
        if a.command=='freeze-model':
            c=v1.cfg_read(a.config)
            m2.freeze_model(c['method2_model_dir'],c['method2_lock'],c['original_scorefile'],c['project'])
        elif a.command=='regression':regression(a)
        else:run(a)
    except Exception as exc:
        p.exit(1, f'ERROR: {exc}\n')


if __name__=='__main__': main()
