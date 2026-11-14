"""Research-only cohort driver; never registers cases or publishes UI results.

The production six-disease scorer remains the only scoring implementation.
Headers, inputs and checkpoint digests are checked before isolated batch runs.
Frequency comparisons measure sensitivity, not accuracy against genotype truth.
"""
import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
from contextlib import contextmanager
import csv
from datetime import datetime, timezone
from decimal import Decimal
import fcntl
import gzip
import hashlib
from itertools import combinations, zip_longest
import json
import math
import os
from pathlib import Path
import re
import signal
import statistics
import subprocess
import sys
import tempfile
import time

from .absolute_risk.registry import DISEASES
from .absolute_risk.multi import score_rows, ancestry_ready
from .absolute_risk.engine import GROUPS, group
from .multidisease.common import load_models, rows, sha, write_json
from .multidisease.run import checkpoint, fingerprint

SCHEMA = 'PRS_COHORT_STUDY_V1'
REPO = Path(__file__).resolve().parents[2]


def platform():
    from app import config
    from app.services import prs_jobs, prs_paths, nextflow_work, container_runtime, prs_presentation
    return config, prs_jobs, prs_paths, nextflow_work, container_runtime, prs_presentation


def now():
    return datetime.now(timezone.utc).isoformat()


def read(path):
    with Path(path).open() as f:
        return json.load(f)


def code_signature():
    return {str(p.relative_to(REPO)):sha(p) for p in
            [Path(__file__), REPO/'research/prs/validation.py',
             REPO/'backend/app/services/prs_paths.py',
             REPO/'scripts/inventory_prs_gvcfs.py',REPO/'research/prs/multidisease_cli.py',
             *sorted((REPO/'research/prs/absolute_risk').glob('*.py')),
             *sorted((REPO/'research/prs/multidisease').glob('*.py'))]}


def read_inventory(directory):
    from scripts.inventory_prs_gvcfs import classify, EXPECTED_VAL_IDS
    directory = Path(directory).resolve()
    data = read(directory/'inventory.json')
    if data.get('schema') != 'PRS_GVCF_INVENTORY_V1' or data.get('discovery_complete') is not True:
        raise ValueError('需要完整且無搜尋錯誤的 inventory.json')
    with (directory/'samples.tsv').open(newline='') as f:
        samples = list(csv.DictReader(f, delimiter='\t'))
    source_rows = {r['sample_id']:r for r in data['samples']}
    if len(source_rows) != len(data['samples']) or len(samples) != len(source_rows):
        raise ValueError('Inventory JSON/TSV row identity differs')
    seen, paths, result = set(), set(), []
    raw_root = Path(data['root']).resolve()
    for row in samples:
        sid = row['sample_id']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}',sid) or '..' in sid or sid in seen:
            raise ValueError('Invalid or duplicate inventory sample ID')
        seen.add(sid)
        kind = classify(sid)
        if kind is None or kind != (row['cohort'],row['standard_id']):
            raise ValueError('Inventory cohort/standard identity differs')
        original = source_rows.get(sid,{})
        for field in ('selection','cohort','standard_id','file_count',
                      'nonempty_matching_file_count','distinct_physical_files'):
            if str(original.get(field)) != row[field]:
                raise ValueError('Inventory JSON/TSV differs: '+sid+'/'+field)
        selected = json.loads(row['paths'])
        if (row['selection'] != 'ONE_FILE' or row['file_count'] != '1'
                or row['nonempty_matching_file_count'] != '1' or len(selected) != 1
                or selected != original.get('paths')):
            raise ValueError('每個樣本必須先唯一確認非空來源：'+sid)
        path = Path(selected[0])
        if not path.is_absolute() or not path.resolve().is_relative_to(raw_root) or str(path.resolve()) in paths:
            raise ValueError('來源不在盤點root或重複指向同一檔案：'+sid)
        paths.add(str(path.resolve()))
        result.append(dict(sample_id=sid,cohort=row['cohort'],standard_id=row['standard_id'],
                           gvcf=str(path),run_id=path.relative_to(Path(data['root'])).parts[0]))
    if any(s not in seen for s in EXPECTED_VAL_IDS):
        raise ValueError('本輪需要完整 VAL-45～54 標準品清單')
    return result, dict(directory=str(directory),root=str(raw_root),
        json_sha256=sha(directory/'inventory.json'),tsv_sha256=sha(directory/'samples.tsv'))


def header(path, sample_id, standard_id, fasta_lengths):
    """Read a bounded header only; do not infer caller version from filename."""
    declarations, formats, caller, references, symbolic = {}, set(), [], [], False
    count = 0
    opener = gzip.open if str(path).endswith(('.gz','.bgz')) else open
    with opener(path,'rt') as f:
        for _ in range(100000):
            line=f.readline(1_000_000);count+=len(line)
            if count>32_000_000 or not line.endswith('\n'):
                raise ValueError('gVCF header missing or exceeds bounds')
            if line.startswith('##contig=<'):
                c=re.search(r'ID=([^,>]+)',line);n=re.search(r'length=(\d+)',line,re.I)
                if c and n:
                    ident=c[1].removeprefix('chr')
                    if ident in declarations:raise ValueError('Duplicate/ambiguous gVCF contig')
                    declarations[ident]=int(n[1])
            m=re.match(r'##FORMAT=<ID=([^,>]+)',line)
            if m:formats.add(m[1])
            symbolic |= bool(re.match(r'##ALT=<ID=(?:NON_REF|\*),',line))
            if re.match(r'##(?:source|DRAGENVersion|GATKVersion|DeepVariantVersion)=',line,re.I):
                caller.append(line.strip()[:2000])
            if line.lower().startswith('##reference='):references.append(line.strip()[:2000])
            if line.startswith('#CHROM\t'):
                fields=line.rstrip('\r\n').split('\t')
                aliases={'HG001':{'HG001','NA12878'},'HG002':{'HG002','NA24385'}}
                if len(fields)!=10 or fields[9] not in ({sample_id}|aliases.get(standard_id,set())):
                    raise ValueError('gVCF single-sample identity differs')
                if declarations.get('1')!=248956422 or any(
                        declarations.get(str(c))!=fasta_lengths.get(str(c)) for c in range(1,23)):
                    raise ValueError('gVCF autosomal contigs do not match configured GRCh38 FASTA')
                if not symbolic or not {'GT','DP','GQ'}<=formats:
                    raise ValueError('gVCF reference-block allele or GT/DP/GQ declarations missing')
                indexes=[Path(str(path)+s) for s in ('.tbi','.csi') if Path(str(path)+s).is_file()]
                if not indexes or any(p.stat().st_size==0 or not os.access(p,os.R_OK) for p in indexes):
                    raise ValueError('gVCF readable nonempty index missing')
                return dict(status='HEADER_PASS',header_sample=fields[9],caller=caller,
                    caller_version_status='observed' if any(re.search(r'\d+\.\d+',s) for s in caller) else 'unknown',
                    reference=references,input_signature=[fingerprint(path),*[fingerprint(p) for p in indexes]],
                    coverage_status='NOT_CHECKED',genotype_identity_status='NOT_CHECKED')
    raise ValueError('gVCF header not found')


def settings(config_path):
    _, jobs, *_=platform()
    if not Path(config_path).is_file():
        raise ValueError('找不到PRS設定：'+str(config_path)+'；可用 prepare --config 指定平台設定檔')
    previous=os.environ.get('NGS_UI_PRS_CONFIG')
    try:
        os.environ['NGS_UI_PRS_CONFIG']=str(config_path)
        return jobs.settings()
    finally:
        if previous is None:os.environ.pop('NGS_UI_PRS_CONFIG',None)
        else:os.environ['NGS_UI_PRS_CONFIG']=previous


def study_read(study):
    study=Path(study).resolve();plan=read(study/'plan.json')
    if plan.get('schema')!=SCHEMA or plan.get('study')!=str(study):
        raise ValueError('Invalid/moved cohort study directory')
    return plan


def verify_current(plan):
    _,jobs,*_=platform()
    cfg,h=settings(plan['config_path'])
    if (h!=plan['config_sha256'] or jobs.resource_signature(cfg)!=plan['resources']
            or code_signature()!=plan['code_sha256']):
        raise ValueError('設定、模型、資源或程式已變更；請保留舊study，使用新的study目錄')
    return cfg


def prepare(inventory_dir, study, config_path, score_dir):
    config,jobs,paths,guard,*_=platform()
    print('[cohort] config='+str(config_path)+'；score sources='+str(score_dir),flush=True)
    cfg,cfg_hash=settings(config_path)
    if not cfg.get('multi_model_dir'):raise ValueError('需要已啟用的六疾病 multi_model_dir')
    if (cfg.get('min_dp'),cfg.get('min_gq'),cfg.get('palindromic_mode'))!=(10,20,'assume-direct'):
        raise ValueError('研究批次沿用固定DP10/GQ20/assume-direct，不變更計分政策')
    manifest,mhash=load_models(cfg['multi_model_dir'])
    samples,source=read_inventory(inventory_dir)
    study=Path(study).absolute()
    _,checked,_=paths.validation_config(study/'runtime',cfg)
    study=Path(checked['project']).parent  # Freeze the canonical existing work-root alias.
    inputs=dict(inventory=source,config_sha256=cfg_hash,resources=jobs.resource_signature(cfg),
                model_sha256=mhash,code_sha256=code_signature())
    if study.exists():
        plan=study_read(study)
        if plan['inputs']!=inputs:raise ValueError('Study inputs changed; use a new --study')
        for item in plan['samples']:
            if item['header']['status']=='HEADER_PASS':verify_sample(item)
        print('[cohort] 沿用已核對study：'+str(study),flush=True)
        return plan
    from .validation import model_audit
    with guard.acquire_lease(config.TERTIARY_NF_WORK_ROOT,'prs'):
        study.mkdir(parents=True,mode=0o700)
        (study/'frozen').mkdir()
        for name in ('manifest.json','reference_scores.json'):
            (study/'frozen'/name).write_bytes((Path(cfg['multi_model_dir'])/name).read_bytes())
        audit=model_audit(cfg['multi_model_dir'],score_dir,cfg['original_scorefile'])
        if (audit['model_sha256']!=mhash or sha(study/'frozen/manifest.json')!=mhash
                or sha(study/'frozen/reference_scores.json')!=manifest['sha256']['reference_scores.json']):
            raise ValueError('Model changed while preparing frozen inputs')
        write_json(study/'model_source_audit.json',audit)
        lengths={}
        for line in Path(cfg['fasta']+'.fai').read_text().splitlines():
            fields=line.split('\t');key=fields[0].removeprefix('chr')
            if key in lengths:raise ValueError('Ambiguous FASTA contigs')
            lengths[key]=int(fields[1])
        for i,item in enumerate(samples,1):
            try:item['header']=header(item['gvcf'],item['sample_id'],item['standard_id'],lengths)
            except (ValueError,OSError,EOFError) as exc:
                item['header']=dict(status='HEADER_FAILED',error=str(exc))
            if i%50==0 or i==len(samples):
                print('[cohort] header '+str(i)+'/'+str(len(samples)),flush=True)
        current,current_hash=settings(config_path)
        if current_hash!=cfg_hash or jobs.resource_signature(current)!=inputs['resources']:
            raise ValueError('Resources/config changed during preflight; use a new study')
        plan=dict(schema=SCHEMA,study=str(study),prepared_at=now(),inputs=inputs,samples=samples,
                  config_path=str(Path(config_path).resolve()),config_sha256=cfg_hash,
                  model_sha256=mhash,resources=inputs['resources'],code_sha256=inputs['code_sha256'],
                  score_dir=str(Path(score_dir).resolve()),
                  source_audit_sha256=sha(study/'model_source_audit.json'))
        write_json(study/'plan.json',plan)
        write_json(study/'preflight_summary.json',dict(schema='PRS_COHORT_PREFLIGHT_V1',
            model_sha256=mhash,source_model_audit=audit['status'],
            counts=dict(Counter(r['header']['status'] for r in samples)),
            cohort_counts=dict(Counter(r['cohort'] for r in samples)),
            caller_headers=dict(Counter(' | '.join(r['header'].get('caller',[])) or 'unknown' for r in samples)),
            limitations=['Header PASS does not establish WGS coverage, genotype identity, ancestry or PRS accuracy.']))
    print('[cohort] 前置檢查完成：'+str(study/'preflight_summary.json'),flush=True)
    return plan


def verify_sample(item):
    h=item['header']
    if h['status']!='HEADER_PASS':raise ValueError('Header未通過：'+item['sample_id'])
    path=Path(item['gvcf'])
    indexes=[Path(str(path)+s) for s in ('.tbi','.csi') if Path(str(path)+s).is_file()]
    actual=[fingerprint(path),*[fingerprint(p) for p in indexes]]
    if actual!=h['input_signature']:raise ValueError('gVCF/index changed after preflight: '+item['sample_id'])


def select_samples(plan,cohort,per_run=0,all_healthy=False):
    kind='standard' if cohort=='standards' else 'healthy_candidate'
    result=[r for r in plan['samples'] if r['cohort']==kind]
    if cohort=='standards' and (per_run or all_healthy):
        raise ValueError('健康人抽樣參數不可套用到標準品')
    if cohort=='healthy':
        if bool(per_run)==bool(all_healthy):raise ValueError('健康人請選 --per-run N 或 --all-healthy')
        if per_run<0:raise ValueError('--per-run must be positive')
        if per_run:
            runs=defaultdict(list)
            for r in result:runs[r['run_id']].append(r)
            result=[]
            for key,items in sorted(runs.items()):
                result.extend(sorted(items,key=lambda r:hashlib.sha256(
                    ('prs-pilot-v1/'+r['sample_id']).encode()).hexdigest())[:per_run])
    if not result:raise ValueError('No samples selected')
    bad=[r['sample_id'] for r in result if r['header']['status']!='HEADER_PASS']
    if bad:raise ValueError('所選樣本header未通過，請看plan.json：'+', '.join(bad[:10]))
    return result


@contextmanager
def study_lock(study):
    with (Path(study)/'execution.lock').open('a') as handle:
        try:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('此study已有工作執行中')
        yield


def command(args,cwd,env,log_path,lease,slot=None):
    """Own the child process group so Ctrl-C/SIGTERM cannot orphan tools."""
    *_,presentation=platform()
    with Path(log_path).open('a') as log:
        log.write('\n[cohort] '+now()+' '+json.dumps(list(map(str,args)),ensure_ascii=False)+'\n');log.flush()
        cursor=Path(log_path).stat().st_size
        p=subprocess.Popen(list(map(str,args)),cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,
                           start_new_session=True,pass_fds=tuple(sorted({lease.fileno()}|
                               ({slot.fileno()} if slot is not None else set()))))
        recent='';failed=False
        try:
            while True:
                status=p.poll()
                with Path(log_path).open('rb') as stream:
                    stream.seek(cursor);chunk=stream.read(1024*1024);cursor=stream.tell()
                recent=(recent+chunk.decode(errors='replace'))[-50000:]
                if presentation.nextflow_failed(recent):raise RuntimeError('Nextflow失敗，請看個案log')
                if status is not None:
                    # Drain all output, including errors printed just before exit.
                    if len(chunk)==1024*1024:continue
                    break
                time.sleep(.5)
            if p.returncode:raise RuntimeError('工具失敗 exit='+str(p.returncode)+'；'+str(log_path))
        except BaseException:
            failed=True;raise
        finally:
            if failed or p.poll() is None:
                try:
                    os.killpg(p.pid,signal.SIGTERM)
                    try:p.wait(timeout=5)
                    except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
                except ProcessLookupError:pass


def runtime(study,item,cfg):
    _,_,paths,*_=platform()
    project=Path(study)/'cases'/item['sample_id']
    effective,layout,wrappers=paths.validation_config(project,cfg)
    paths.prepare(layout,wrappers)
    config_path=project/'config.json'
    if config_path.exists() and read(config_path)!=effective:raise ValueError('Case runtime config changed')
    write_json(config_path,effective)
    return project,effective,layout,config_path


def run_batch(study,cohort,stage,per_run=0,all_healthy=False):
    config,jobs,paths,guard,containers,_=platform()
    plan=study_read(study);cfg=verify_current(plan)
    selected=select_samples(plan,cohort,per_run,all_healthy)
    with study_lock(study),guard.acquire_lease(config.TERTIARY_NF_WORK_ROOT,'prs') as lease:
        for index,item in enumerate(selected,1):
            sid=item['sample_id'];verify_sample(item);verify_current(plan)
            jobs.root().mkdir(parents=True,exist_ok=True)
            print('[cohort] '+str(index)+'/'+str(len(selected))+' '+sid+' 等待PRS名額',flush=True)
            with (jobs.root()/'.execution.lock').open('a') as slot:
                fcntl.flock(slot,fcntl.LOCK_EX)
                verify_sample(item);verify_current(plan)
                project,effective,layout,config_path=runtime(study,item,cfg)
                previous=read(project/'state.json') if (project/'state.json').is_file() else {}
                started=time.monotonic()
                state=dict(sample_id=sid,status='RUNNING',started_at=now(),updated_at=now(),pid=os.getpid(),stage=stage,
                           verified_stages=previous.get('verified_stages',[]))
                write_json(project/'state.json',state)
                try:
                    cwd,env=paths.environment(layout,'validation','score')
                    env.update(guard.lease_env(config.TERTIARY_NF_WORK_ROOT,'prs',lease))
                    token=__import__('uuid').uuid4().hex
                    with tempfile.NamedTemporaryFile(mode='w',dir=paths.container_tmp(layout)) as marker:
                        marker.write(token);marker.flush()
                        clean=containers.without_bind_destinations(env,
                            {str(paths.container_tmp(layout)),'/tmp','/var/tmp'})
                        command(paths.container_tmp_probe(effective,layout,Path(marker.name),token),
                                cwd,clean,project/'tools.log',lease,slot)
                    stages=['score','ancestry'] if stage=='both' else [stage]
                    for phase in stages:
                        state.update(phase=phase,updated_at=now());write_json(project/'state.json',state)
                        cwd,env=paths.environment(layout,'validation',phase)
                        env.update(guard.lease_env(config.TERTIARY_NF_WORK_ROOT,'prs',lease))
                        command([sys.executable,REPO/'research/prs/multidisease_cli.py',
                            '--config',config_path,'--sample',sid,'--gvcf',item['gvcf'],
                            '--output',project/'multi','--stage',phase],cwd,env,project/'tools.log',lease,slot)
                        verify_sample(item);verify_current(plan)
                        state['verified_stages']=sorted(set(state['verified_stages'])|{phase})
                        write_json(project/'state.json',state)
                    verify_sample(item);verify_current(plan)
                    state.update(status='COMPLETE',updated_at=now(),elapsed_seconds=time.monotonic()-started)
                except BaseException as exc:
                    state.update(status='INTERRUPTED' if isinstance(exc,KeyboardInterrupt) else 'FAILED',
                                 error=str(exc),updated_at=now(),elapsed_seconds=time.monotonic()-started)
                    write_json(project/'state.json',state);raise
                write_json(project/'state.json',state)
            print('[cohort] '+sid+' '+stage+'完成',flush=True)
            time.sleep(1.1)  # Release the slot between cases, including for UI workers.


def completed(plan):
    """Read only checksum-verified completions, not early summary.json files."""
    reference=Path(plan['study'])/'frozen/reference_scores.json'
    manifest_path=Path(plan['study'])/'frozen/manifest.json'
    if sha(manifest_path)!=plan['model_sha256']:raise ValueError('Frozen model manifest changed')
    manifest=read(manifest_path)
    if sha(reference)!=manifest['sha256']['reference_scores.json']:
        raise ValueError('Frozen reference distribution changed')
    for item in plan['samples']:
        root=Path(plan['study'])/'cases'/item['sample_id']/'multi'
        state_path=root.parent/'state.json'
        if not state_path.is_file() or 'score' not in read(state_path).get('verified_stages',[]):continue
        if not (root/'score.complete.json').is_file():continue
        expected=read(root/'inputs.lock.json')
        if (expected['sample']!=item['sample_id'] or expected['model']!=plan['model_sha256']
                or expected['gvcf']!=item['header']['input_signature'][0]):
            raise ValueError('Score completion does not belong to this study')
        checkpoint(root/'score.complete.json',expected)
        summary=read(root/'summary.json');score_rows(summary)
        if summary['sample']!=item['sample_id'] or summary['model_manifest_sha256']!=plan['model_sha256']:
            raise ValueError('Summary identity differs')
        qc=None
        if ('ancestry' in read(state_path).get('verified_stages',[])
                and checkpoint(root/'ancestry.complete.json',expected)):qc=read(root/'qc.json')
        yield item,summary,qc,root


def distance(a,b):
    a=sorted(a);b=sorted(b)
    return max(abs(bisect_right(a,v)/len(a)-bisect_right(b,v)/len(b)) for v in set(a+b))


def describe(values,reference):
    return dict(n=len(values),mean=statistics.mean(values),median=statistics.median(values),
                sd=statistics.stdev(values) if len(values)>1 else None,
                minimum=min(values),maximum=max(values),
                mean_shift_reference_sd=(statistics.mean(values)-statistics.mean(reference))/statistics.stdev(reference),
                max_ecdf_distance=distance(values,reference),p_value=None)


def concordance(a,b):
    counts=Counter({k:0 for k in ('union_loci','joint_observed','discordant','both_missing','one_missing')})
    for chrom in range(1,23):
        for left,right in zip_longest(rows(a/f'chr{chrom}.audit.tsv.gz'),rows(b/f'chr{chrom}.audit.tsv.gz')):
            if left is None or right is None or left['id']!=right['id']:
                raise ValueError('Replicate audit membership differs')
            x,y=left['dosage'],right['dosage']
            counts['union_loci']+=1
            if x!='' and y!='':
                counts['joint_observed']+=1;counts['discordant']+=x!=y
            elif x=='' and y=='':counts['both_missing']+=1
            else:counts['one_missing']+=1
    return {**dict(counts),'observed_concordance':1-counts['discordant']/counts['joint_observed']
            if counts['joint_observed'] else None}


def report(study,genotypes=False):
    config,_,_,guard,*_=platform()
    plan=study_read(study)
    with study_lock(study),guard.acquire_lease(config.TERTIARY_NF_WORK_ROOT,'prs'):
        data=list(completed(plan))
        if not data:raise ValueError('尚無已完成且checksum通過的計分結果')
        reference=read(Path(study)/'frozen/reference_scores.json')
        output=Path(tempfile.mkdtemp(prefix='report_',dir=study))
        scores=[];replicates=[];populations={}
        for disease in DISEASES:
            ref=list(reference[disease].values())
            eligible=[];all_health=[];runs=defaultdict(list);all_runs=defaultdict(list);standards=defaultdict(list)
            for item,s,qc,root in data:
                row=s['diseases'][disease]
                eas=ancestry_ready({**s,'ancestry_status':qc or {}})
                scores.append(dict(sample_id=item['sample_id'],cohort=item['cohort'],standard_id=item['standard_id'],
                    run_id=item['run_id'],disease=disease,**row,ancestry_checked=qc is not None,
                    confident_eas=eas,genotype_truth_compared=False))
                if item['cohort']=='standard':standards[item['standard_id']].append(row)
                else:
                    all_health.append(row)
                    all_runs[item['run_id']].append(row)
                    if eas and row['qc_ready']:eligible.append(row);runs[item['run_id']].append(row['total_sum'])
            for standard,values in standards.items():
                source_runs={item['run_id'] for item,_,_,_ in data if item['standard_id']==standard}
                replicates.append(dict(standard=standard,disease=disease,n_replicate_files=len(values),
                    sequencing_run_count=len(source_runs),
                    repeatability_status='DESCRIPTIVE_REPEATS' if len(values)>1 else 'INSUFFICIENT_REPLICATES',
                    score_range=max(r['total_sum'] for r in values)-min(r['total_sum'] for r in values) if len(values)>1 else None,
                    percentile_range_pp=max(r['eas_percentile'] for r in values)-min(r['eas_percentile'] for r in values) if len(values)>1 else None,
                    min_observed_fraction=min(r['observed_loci']/r['n_loci'] for r in values),
                    stratum_count=len({r['intervene_group'] for r in values})))
            populations[disease]=dict(scored_health_ids=len(all_health),confident_eas_quality_ids=len(eligible),
                stratum_expected_proportions=[.2,.2,.2,.2,.1,.05,.05],
                all_scored_descriptive_not_ancestry_validated=describe([r['total_sum'] for r in all_health],ref) if all_health else None,
                all_scored_by_run={run:dict(score=describe([r['total_sum'] for r in v],ref),
                    mean_percentile=statistics.mean(r['eas_percentile'] for r in v),
                    mean_missing_fraction=statistics.mean(r['missing_loci']/r['n_loci'] for r in v),
                    max_high_weight_missing=max(r.get('high_weight_missing',0) for r in v)) for run,v in all_runs.items()},
                observed_locus_quality_pass_ids=sum(r['qc_ready'] for r in all_health),
                all_scored_mean_percentile=statistics.mean(r['eas_percentile'] for r in all_health) if all_health else None,
                all_scored_percentile_deciles=[sum(min(9,int(r['eas_percentile']//10))==i for r in all_health) for i in range(10)],
                eligible_score_distribution=describe([r['total_sum'] for r in eligible],ref) if eligible else None,
                eligible_percentile_deciles=[sum(min(9,int(r['eas_percentile']//10))==i for r in eligible) for i in range(10)],
                eligible_stratum_counts={g:sum(r['intervene_group']==g for r in eligible) for g in GROUPS},
                eligible_by_run={run:describe(v,ref) for run,v in runs.items()})
        genotype_rows=[]
        if genotypes:
            for standard in ('HG001','HG002'):
                items=[r for r in data if r[0]['standard_id']==standard]
                for a,b in combinations(items,2):
                    print('[cohort] genotype concordance '+a[0]['sample_id']+'/'+b[0]['sample_id'],flush=True)
                    genotype_rows.append(dict(standard=standard,a=a[0]['sample_id'],b=b[0]['sample_id'],**concordance(a[3],b[3])))
        write_json(output/'summary.json',dict(schema='PRS_COHORT_REPORT_V1',created_at=now(),
            planned_sample_ids=len(plan['samples']),scored_sample_ids=len(data),
            not_scored_sample_ids=len(plan['samples'])-len(data),
            health_scored_ids=sum(r[0]['cohort']=='healthy_candidate' for r in data),
            standards_scored_ids=sum(r[0]['cohort']=='standard' for r in data),model_sha256=plan['model_sha256'],
            standard_score_repeatability=replicates,standard_genotype_repeatability=genotype_rows,
            population=populations,limitations=['No GIAB genotype-truth comparison in this report.',
                'Percentiles of standards are technical descriptors, not EAS clinical interpretations.',
                'Counts are sample IDs; relatedness and biological-person deduplication not performed.',
                'Empirical distribution distance is descriptive; no uniformity P-value or clinical validation.',
                'Case scores are never published to UI or used to replace a reference model.']))
        tsv(output/'scores.tsv',scores);tsv(output/'standard_repeatability.tsv',replicates)
        if genotype_rows:tsv(output/'standard_genotypes.tsv',genotype_rows)
        print('[cohort] 報告：'+str(output/'summary.json'),flush=True)
        return output


def tsv(path,data):
    if not data:return
    with Path(path).open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');writer.writeheader();writer.writerows(data)


def frequency_sensitivity(study):
    """Same-case missing-locus frequency substitution, NOT a truth benchmark."""
    config,jobs,paths,guard,*_=platform()
    plan=study_read(study);cfg=verify_current(plan)
    model=Path(cfg['multi_model_dir']);prefix=model/'reference'
    for suffix in ('.pgen','.pvar','.psam'):
        if not Path(str(prefix)+suffix).is_file():raise ValueError('缺少模型準備保留的 reference'+suffix)
    with study_lock(study),guard.acquire_lease(config.TERTIARY_NF_WORK_ROOT,'prs') as lease:
        data=list(completed(plan))
        if not data:raise ValueError('尚無完成計分可比較')
        output=Path(tempfile.mkdtemp(prefix='frequency_',dir=study))
        effective,layout,wrappers=paths.validation_config(output/'runtime',cfg);paths.prepare(layout,wrappers)
        cwd,env=paths.environment(layout,'frequency','reference')
        env.update(guard.lease_env(config.TERTIARY_NF_WORK_ROOT,'prs',lease))
        from .multidisease.common import plink
        genotypes={s:fingerprint(str(prefix)+s) for s in ('.pgen','.pvar','.psam')}
        jobs.root().mkdir(parents=True,exist_ok=True)
        with (jobs.root()/'.execution.lock').open('a') as slot:
            fcntl.flock(slot,fcntl.LOCK_EX)
            command(plink(effective,'--pfile',prefix,'--keep',model/'reference_calibration_eas_735.keep',
                          '--freq','--out',output/'eas'),cwd,env,output/'tools.log',lease,slot)
        # A matching frequency file alone is insufficient evidence of membership.
        psam=list(rows(Path(str(prefix)+'.psam')))
        selected={r.get('#IID') or r.get('IID') for r in psam}
        keep=members(model/'reference_calibration_eas_735.keep')
        all_keep=members(model/'reference_calibration_3330.keep')
        if len(keep)!=735 or not keep<=selected or len(psam)!=3330 or selected!=all_keep:
            raise ValueError('EAS735/reference3330 membership differs')
        frequency_rows=list(rows(output/'eas.afreq'));af={r['ID']:r for r in frequency_rows}
        baseline=list(rows(model/'reference.afreq'))
        if len(af)!=len(frequency_rows) or set(af)!={r['ID'] for r in baseline}:
            raise ValueError('EAS frequency union differs')
        for r in baseline:
            f=af[r['ID']];p=Decimal(f['ALT_FREQS']);obs=int(f['OBS_CT'])
            if (f['REF'],f['ALT'])!=(r['REF'],r['ALT']) or not 0<=p<=1 or not 0<obs<=1470:
                raise ValueError('EAS frequency allele/count differs')
        if genotypes!={s:fingerprint(str(prefix)+s) for s in genotypes}:
            raise ValueError('Reference genotypes changed during frequency calculation')
        import sqlite3
        db=sqlite3.connect((model/'models.sqlite').resolve().as_uri()+'?mode=ro',uri=True);db.row_factory=sqlite3.Row
        result=[];reference=read(Path(study)/'frozen/reference_scores.json')
        try:
            for item,summary,qc,root in data:
                terms={d:[] for d in DISEASES};seen=set()
                for r in rows(root/'missing_loci.tsv.gz'):
                    key=(r['disease'],r['id'])
                    if key in seen:raise ValueError('Duplicate mask locus')
                    seen.add(key)
                    m=db.execute('SELECT m.*,l.ref,l.alt FROM model m JOIN locus l USING(id) WHERE m.disease=? AND m.id=?',key).fetchone()
                    if not m or Decimal(r['effect_weight'])!=Decimal(m['weight']) or Decimal(r['imputed_dosage'])!=2*Decimal(m['freq']):
                        raise ValueError('Missing-mask model/dosage mismatch')
                    f=af[r['id']]
                    p=Decimal(f['ALT_FREQS']);obs=int(f['OBS_CT'])
                    if (f['REF'],f['ALT'])!=(m['ref'],m['alt']) or not 0<=p<=1 or not 0<obs<=1470:
                        raise ValueError('EAS allele/frequency/count differs')
                    effect=p if m['ea']==m['alt'] else 1-p
                    terms[r['disease']].append(float(Decimal(m['weight'])*2*(effect-Decimal(m['freq']))))
                for d,r in summary['diseases'].items():
                    if sum(k[0]==d for k in seen)!=r['missing_loci']:raise ValueError('Mask count differs from completed summary')
                    delta=math.fsum(terms[d]);alternative=math.fsum([r['total_sum'],delta])
                    ref=sorted(reference[d].values());pct=100*bisect_right(ref,alternative)/735
                    result.append(dict(sample_id=item['sample_id'],cohort=item['cohort'],standard_id=item['standard_id'],
                        disease=d,current_score=r['total_sum'],eas_fill_score=alternative,delta_score=delta,
                        delta_score_sd=delta/statistics.stdev(ref),current_percentile=r['eas_percentile'],
                        eas_fill_percentile=pct,delta_percentile_pp=pct-r['eas_percentile'],
                        changed_stratum=group(pct)!=r['intervene_group'],missing_loci=r['missing_loci']))
        finally:db.close()
        verify_current(plan)
        if genotypes!={s:fingerprint(str(prefix)+s) for s in genotypes}:
            raise ValueError('Reference genotypes changed during sensitivity calculation')
        tsv(output/'sensitivity.tsv',result)
        write_json(output/'summary.json',dict(schema='PRS_FREQUENCY_SENSITIVITY_V1',model_sha256=plan['model_sha256'],
            frequency_sha256=sha(output/'eas.afreq'),keep_sha256=sha(model/'reference_calibration_eas_735.keep'),
            reference_genotype_signatures=genotypes,
            sample_ids=len(data),diseases={d:dict(n=sum(r['disease']==d for r in result),
                max_absolute_percentile_shift_pp=max(abs(r['delta_percentile_pp']) for r in result if r['disease']==d),
                stratum_changes=sum(r['changed_stratum'] for r in result if r['disease']==d)) for d in DISEASES},
            limitations=['Frequency-choice sensitivity on the same observed sample; not genotype-truth accuracy.',
                         'Fixed original EAS735 ECDF; no reference/model activation.',
                         'Reference-masking population benchmark and GIAB accuracy remain separate tasks.']))
        print('[cohort] 補值頻率敏感度：'+str(output/'summary.json'),flush=True)
        return output


def members(path):
    ids=[]
    for r in rows(path):
        iid=r.get('#IID') or r.get('IID')
        if not iid:raise ValueError('Reference keep/PSAM is missing IID')
        ids.append(iid)
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate reference IID')
    return set(ids)


def exported_sample(path):
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('#CHROM\t'):
                fields=line.rstrip('\r\n').split('\t')
                if len(fields)!=10:raise ValueError('Reference export must contain one sample')
                return fields[9]
    raise ValueError('Reference export has no sample header')


def reference_benchmark(study,frequency_dir,output,reference_limit,mask_limit):
    """Bounded reference-genotype x actual-mask pilot; preserves partial truth."""
    config,jobs,paths,guard,*_=platform()
    plan=study_read(study);cfg=verify_current(plan);model=Path(cfg['multi_model_dir'])
    frequency_dir=Path(frequency_dir).resolve();frequency=read(frequency_dir/'summary.json')
    prefix=model/'reference';output=Path(output).absolute()
    _,checked,_=paths.validation_config(output/'runtime',cfg)
    output=Path(checked['project']).parent
    if not 1<=reference_limit<=735 or not 1<=mask_limit<=50:
        raise ValueError('Reference pilot limits: references 1..735, masks 1..50')
    if (frequency.get('schema')!='PRS_FREQUENCY_SENSITIVITY_V1'
            or frequency['model_sha256']!=plan['model_sha256']
            or sha(frequency_dir/'eas.afreq')!=frequency['frequency_sha256']
            or sha(model/'reference_calibration_eas_735.keep')!=frequency['keep_sha256']
            or {s:fingerprint(str(prefix)+s) for s in ('.pgen','.pvar','.psam')}!=frequency['reference_genotype_signatures']):
        raise ValueError('Frequency/reference inputs no longer match this study')
    reference_ids=sorted(members(model/'reference_calibration_eas_735.keep'),
                         key=lambda iid:hashlib.sha256(('reference-pilot-v1/'+iid).encode()).hexdigest())[:reference_limit]
    from . import validation as audit
    from .multidisease.common import plink
    with study_lock(study),guard.acquire_lease(config.TERTIARY_NF_WORK_ROOT,'prs') as lease:
        data=list(completed(plan))
        by_run=defaultdict(list)
        for item,s,qc,root in data:by_run[item['run_id']].append((item,root))
        ordered=[]
        for run in sorted(by_run,key=lambda x:hashlib.sha256(('mask-pilot-v1/'+x).encode()).hexdigest()):
            items=sorted(by_run[run],key=lambda x:hashlib.sha256(x[0]['sample_id'].encode()).hexdigest())
            ordered.append(items.pop(0));by_run[run]=items
        for run in sorted(by_run):ordered.extend(by_run[run])
        masks=ordered[:mask_limit]
        if not masks:raise ValueError('No completed cases supply actual missing masks')
        frozen=dict(schema='PRS_REFERENCE_MASK_PILOT_V1',model_sha256=plan['model_sha256'],
            references=reference_ids,masks=[dict(sample_id=i['sample_id'],sha256=sha(r/'missing_loci.tsv.gz')) for i,r in masks],
            frequency_sha256=frequency['frequency_sha256'],code_sha256=code_signature(),
            reference_signatures=frequency['reference_genotype_signatures'])
        if output.exists():
            if read(output/'inputs.json')!=frozen:raise ValueError('Reference pilot changed; use a new --output')
        else:
            output.mkdir(parents=True,mode=0o700);write_json(output/'inputs.json',frozen)
        effective,layout,wrappers=paths.validation_config(output/'runtime',cfg);paths.prepare(layout,wrappers)
        cwd,env=paths.environment(layout,'reference-pilot','export')
        env.update(guard.lease_env(config.TERTIARY_NF_WORK_ROOT,'prs',lease))
        reference=read(Path(study)/'frozen/reference_scores.json');reports=[]
        for index,iid in enumerate(reference_ids,1):
            token=hashlib.sha256(iid.encode()).hexdigest()[:20]
            dest=output/token;dest.mkdir(exist_ok=True)
            vcf=dest/'truth.vcf.gz';marker=dest/'export.json'
            if marker.exists():
                if read(marker)!=dict(iid=iid,sha256=sha(vcf)):raise ValueError('Reference export changed')
            else:
                keep=dest/'sample.keep';keep.write_text('#IID\n'+iid+'\n')
                with (jobs.root()/'.execution.lock').open('a') as slot:
                    fcntl.flock(slot,fcntl.LOCK_EX)
                    command(plink(effective,'--pfile',prefix,'--keep',keep,'--export','vcf','bgz',
                                  'id-paste=iid','--out',dest/'truth'),cwd,env,output/'tools.log',lease,slot)
                if exported_sample(vcf)!=iid:raise ValueError('Reference export sample differs from requested IID')
                write_json(marker,dict(iid=iid,sha256=sha(vcf)))
            for item,root in masks:
                path=dest/(item['sample_id']+'.json')
                receipt=dest/(item['sample_id']+'.complete.json')
                expected=dict(model=plan['model_sha256'],code=frozen['code_sha256'],
                    vcf=sha(vcf),mask=sha(root/'missing_loci.tsv.gz'),frequency=frequency['frequency_sha256'])
                print('[cohort] reference '+str(index)+'/'+str(len(reference_ids))+' mask '+item['sample_id'],flush=True)
                if checkpoint(receipt,expected):
                    value=read(path)
                    if (value['model_sha256']!=plan['model_sha256'] or value['inputs'].get('vcf')!=sha(vcf)
                            or value['inputs'].get('mask')!=sha(root/'missing_loci.tsv.gz')
                            or value['inputs'].get('eas_frequency')!=frequency['frequency_sha256']):
                        raise ValueError('Saved reference comparison changed')
                else:
                    value=audit.benchmark(model,plan['score_dir'],vcf,chd=cfg['original_scorefile'],
                        mask=root/'missing_loci.tsv.gz',eas_frequency=frequency_dir/'eas.afreq')
                    value['reference_calibration']={}
                    for d,r in value['diseases'].items():
                        if r['status']=='COMPLETE_GENOTYPE_COMPARISON':
                            delta=r['sums']['full']-reference[d][iid]
                            value['reference_calibration'][d]=dict(delta=delta,match=abs(delta)<=1e-5)
                    write_json(path,value)
                    write_json(receipt,dict(input=expected,files={path.name:sha(path)}))
                if any(not r['match'] for r in value.get('reference_calibration',{}).values()):
                    raise ValueError('Independent full reference score differs from frozen calibration')
                reports.append(path)
        verify_current(plan)
        if frozen['reference_signatures']!={s:fingerprint(str(prefix)+s) for s in ('.pgen','.pvar','.psam')}:
            raise ValueError('Reference genotypes changed during pilot')
        # Identical mask bytes can occur across cases; count each unique genotype/mask once.
        unique={}
        for p in reports:
            value=read(p);unique.setdefault(tuple(value['inputs'].get(k) for k in ('vcf','mask','eas_frequency')),p)
        result=audit.aggregate(list(unique.values()))
        result.update(reference_people=len(reference_ids),mask_files=len(masks),
            exported_comparisons=len(reports),unique_comparisons=len(unique),
            meaning='Technical reference-genotype masking; no clinical outcomes, no GIAB comparison or reference activation')
        write_json(output/'summary.json',result)
        print('[cohort] reference遮罩試驗：'+str(output/'summary.json'),flush=True)
        return output


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--inventory-dir',type=Path,required=True)
    a.add_argument('--study',type=Path,required=True);a.add_argument('--config',type=Path);a.add_argument('--score-dir',type=Path)
    a=sub.add_parser('run');a.add_argument('--study',type=Path,required=True)
    a.add_argument('--group',choices=('standards','healthy'),default='standards')
    a.add_argument('--stage',choices=('score','ancestry','both'))
    a.add_argument('--per-run',type=int,default=0);a.add_argument('--all-healthy',action='store_true')
    a=sub.add_parser('report');a.add_argument('--study',type=Path,required=True);a.add_argument('--genotypes',action='store_true')
    a=sub.add_parser('frequency');a.add_argument('--study',type=Path,required=True)
    a=sub.add_parser('reference');a.add_argument('--study',type=Path,required=True)
    a.add_argument('--frequency-dir',type=Path,required=True);a.add_argument('--output',type=Path,required=True)
    a.add_argument('--reference-limit',type=int,default=10);a.add_argument('--mask-limit',type=int,default=10)
    args=p.parse_args(argv)
    def cancelled(*_):raise KeyboardInterrupt('使用者終止；保留已完成checkpoint')
    old=signal.signal(signal.SIGTERM,cancelled)
    try:
        config,*_=platform()
        if args.action=='prepare':
            config_path=args.config or Path(os.environ.get('NGS_UI_PRS_CONFIG',config.DATA_ROOT/'prs_config.json'))
            score_dir=args.score_dir or config.PRS_WORK_ROOT/'_reference_prepare_v1/scorefiles'
            prepare(args.inventory_dir,args.study,config_path,score_dir)
        elif args.action=='run':
            run_batch(args.study,args.group,args.stage or 'score',args.per_run,args.all_healthy)
        elif args.action=='report':report(args.study,args.genotypes)
        elif args.action=='frequency':frequency_sensitivity(args.study)
        else:reference_benchmark(args.study,args.frequency_dir,args.output,args.reference_limit,args.mask_limit)
        return 0
    except KeyboardInterrupt:
        print('[cohort] 已終止，保留checkpoint；可用相同命令接續。',file=sys.stderr);return 130
    except (ValueError,OSError,RuntimeError,KeyError,ArithmeticError) as exc:
        print('[cohort] ERROR: '+str(exc),file=sys.stderr);return 1
    finally:signal.signal(signal.SIGTERM,old)
