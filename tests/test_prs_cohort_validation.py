"""Cohort studies must preserve scientific gates and isolate batch outputs."""
import gzip
import json
import os
from bisect import bisect_right
from decimal import Decimal
from pathlib import Path
import sqlite3
import sys
from types import SimpleNamespace

import pytest

from scripts import inventory_prs_gvcfs as inventory
from research.prs import cohort_validation as cohort
from research.prs.multidisease.common import sha,write_json
from research.prs.absolute_risk.registry import DISEASES
from research.prs.absolute_risk.engine import group


def fake_gvcf(path,sample='S1',caller='##source=DRAGENVersion4.3.6\n',chr1=248956422,
              symbolic=True,index=True):
    path.parent.mkdir(parents=True,exist_ok=True)
    text='##fileformat=VCFv4.2\n'+caller
    text+=''.join(f'##contig=<ID=chr{c},length={chr1 if c==1 else 100}>\n' for c in range(1,23))
    if symbolic:text+='##ALT=<ID=NON_REF,Description="test">\n'
    text+=''.join(f'##FORMAT=<ID={k},Number=1,Type=String,Description="test">\n' for k in ('GT','DP','GQ'))
    text+='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n'
    with gzip.open(path,'wt') as f:f.write(text)
    if index:Path(str(path)+'.tbi').write_bytes(b'not-a-real-index')
    return path


@pytest.fixture
def lengths():return {str(c):248956422 if c==1 else 100 for c in range(1,23)}


def test_header_identity_reference_policy_and_unknown_caller(tmp_path,lengths):
    p=fake_gvcf(tmp_path/'DRAGEN4.9.9'/'S1.gvcf.gz',caller='')
    value=cohort.header(p,'S1','',lengths)
    assert value['status']=='HEADER_PASS' and value['caller_version_status']=='unknown'
    assert value['coverage_status']==value['genotype_identity_status']=='NOT_CHECKED'
    q=fake_gvcf(tmp_path/'VAL-45.gvcf.gz',sample='NA12878')
    assert cohort.header(q,'VAL-45','HG001',lengths)['header_sample']=='NA12878'
    with pytest.raises(ValueError,match='identity'):cohort.header(q,'26T00001','',lengths)


@pytest.mark.parametrize('options,match',[
    ({'chr1':249250621},'GRCh38'),({'symbolic':False},'reference-block'),
    ({'sample':'S1\tS2'},'single-sample'),({'index':False},'index')])
def test_bad_header_is_not_promoted_to_ready(tmp_path,lengths,options,match):
    p=fake_gvcf(tmp_path/'S1.gvcf.gz',**options)
    with pytest.raises(ValueError,match=match):cohort.header(p,'S1','',lengths)


def test_input_signature_refuses_changed_index(tmp_path,lengths):
    p=fake_gvcf(tmp_path/'S1.gvcf.gz')
    item=dict(sample_id='S1',gvcf=str(p),header=cohort.header(p,'S1','',lengths))
    cohort.verify_sample(item)
    Path(str(p)+'.tbi').write_bytes(b'changed')
    with pytest.raises(ValueError,match='changed'):cohort.verify_sample(item)


def test_inventory_requires_complete_unique_consistent_sources(tmp_path):
    raw=tmp_path/'raw'
    for n in range(45,55):
        sid='VAL-'+str(n);fake_gvcf(raw/'run_a'/'other'/sid/'germline_seq'/(sid+'.gvcf.gz'),sample=sid)
    fake_gvcf(raw/'run_b'/'other'/'26T00001'/'germline_seq'/'26T00001.gvcf.gz',sample='26T00001')
    directory=inventory.write_reports(inventory.scan(raw),tmp_path/'inventory')
    items,_=cohort.read_inventory(directory)
    assert len(items)==11 and sum(i['standard_id']=='HG001' for i in items)==5
    data=cohort.read(directory/'inventory.json');data['discovery_complete']=False
    write_json(directory/'inventory.json',data)
    with pytest.raises(ValueError,match='完整'):cohort.read_inventory(directory)


def test_health_selection_is_bounded_deterministic_and_not_score_based():
    samples=[dict(sample_id='26T0'+str(i),cohort='healthy_candidate',run_id='run'+str(i%3),
                  header=dict(status='HEADER_PASS')) for i in range(12)]
    plan={'samples':samples}
    a=cohort.select_samples(plan,'healthy',per_run=1)
    b=cohort.select_samples({'samples':list(reversed(samples))},'healthy',per_run=1)
    assert a==b and len(a)==3 and len({r['run_id'] for r in a})==3
    with pytest.raises(ValueError,match='--per-run'):cohort.select_samples(plan,'healthy')
    assert len(cohort.select_samples(plan,'healthy',all_healthy=True))==12
    with pytest.raises(ValueError,match='標準品'):cohort.select_samples(plan,'standards',all_healthy=True)


@pytest.fixture
def runtime(tmp_path,monkeypatch):
    from app import config
    from app.services import prs_paths
    root=tmp_path.resolve();legacy=root/'resources';legacy.mkdir()
    source=legacy/'nextflow.config';source.write_text('singularity.enabled=true\n')
    cfg={k:str(legacy/k) for k in ('raw_script','original_scorefile','fasta','dbsnp','plink_img',
        'reference_panel','method2_lock','method2_model_dir','multi_model_dir')}
    cfg.update(project=str(legacy),raw_work=str(legacy),raw_nextflow_config=str(source),
               ancestry_nextflow_config=str(source),min_dp=10,min_gq=20,palindromic_mode='assume-direct')
    monkeypatch.setattr(config,'PRS_WORK_ROOT',root/'nf_work/prs')
    monkeypatch.setattr(config,'TERTIARY_NF_WORK_ROOT',root/'nf_work')
    monkeypatch.setattr(config,'DATA_ROOT',root/'data')
    return config,prs_paths,cfg


def test_validation_paths_reuse_tmp_rules_and_cannot_alias_cases_or_resources(runtime):
    config,paths,cfg=runtime
    study=config.PRS_WORK_ROOT/'_validation_20261114/cohort_v1'
    effective,layout,wrappers=paths.validation_config(study/'cases/VAL-45',cfg)
    assert effective['multi_model_dir']==cfg['multi_model_dir']
    assert effective['project']==str(study/'cases/VAL-45')
    assert all(str(study/'cases/VAL-45/.runtime/tmp') in text for text in wrappers.values())
    paths.prepare(layout,wrappers)
    _,env=paths.environment(layout,'case','score')
    assert str(study/'cases/VAL-45/.runtime/tmp')==env['TMPDIR']
    for bad in (config.PRS_WORK_ROOT/'REAL_CASE',config.PRS_WORK_ROOT/'_validation_x',Path(cfg['project'])/'outside'):
        with pytest.raises(ValueError):paths.validation_config(bad,cfg)
    alias=study/'cases/alias';alias.symlink_to(Path(cfg['project']),target_is_directory=True)
    with pytest.raises(ValueError,match='symlink'):paths.validation_config(alias,cfg)
    normal=paths.execution_config('CASE',cfg)
    assert normal[0]['project']==str(config.PRS_WORK_ROOT/'CASE')


def test_repeat_genotypes_separate_observed_discordance_and_missingness(tmp_path):
    a=tmp_path/'a';b=tmp_path/'b';a.mkdir();b.mkdir()
    for chrom in range(1,23):
        for root,values in [(a,('0','')),(b,('1' if chrom==1 else '0','1'))]:
            with gzip.open(root/f'chr{chrom}.audit.tsv.gz','wt') as f:
                f.write('id\tdosage\tstatus\n')
                for i,d in enumerate(values):f.write(f'{chrom}:{i}\t{d}\tTEST\n')
    r=cohort.concordance(a,b)
    assert r['union_loci']==44 and r['joint_observed']==22 and r['discordant']==1
    assert r['one_missing']==22 and r['both_missing']==0
    assert r['observed_concordance']==pytest.approx(21/22)


def test_validation_accepts_existing_work_root_ssd_alias_but_rejects_child_aliases(runtime,monkeypatch):
    config,paths,cfg=runtime
    physical=config.PRS_WORK_ROOT;physical.mkdir(parents=True)
    logical=physical.parent.parent/'ssd-work';logical.symlink_to(physical,target_is_directory=True)
    monkeypatch.setattr(config,'PRS_WORK_ROOT',logical)
    _,layout,wrappers=paths.validation_config(logical/'_validation_x/study/runtime',cfg)
    assert layout['project']==str(physical/'_validation_x/study/runtime')
    paths.prepare(layout,wrappers)
    assert paths.validation_config(physical/'_validation_x/study/runtime',cfg)[1]==layout
    (physical/'_validation_x/study/other').symlink_to(physical/'_validation_x/study/runtime',target_is_directory=True)
    with pytest.raises(ValueError,match='symlink'):
        paths.validation_config(logical/'_validation_x/study/other',cfg)


def study_results(study):
    ref={d:{str(i):i/735 for i in range(735)} for d in DISEASES}
    write_json(study/'frozen/reference_scores.json',ref)
    write_json(study/'frozen/manifest.json',dict(sha256={'reference_scores.json':sha(study/'frozen/reference_scores.json')}))
    model_hash=sha(study/'frozen/manifest.json')
    item=dict(sample_id='VAL-45',cohort='standard',standard_id='HG001',run_id='ONE_RUN',
              header=dict(input_signature=[['/not-read',10,20,30]]))
    root=study/'cases/VAL-45/multi';root.mkdir(parents=True)
    summary=dict(schema='PRS_MULTI_SCORE_V1',sample='VAL-45',reference_id='unrelated3330_EAS735_v1',
                 model_manifest_sha256=model_hash,diseases={})
    for d,v in DISEASES.items():
        n=610305 if d=='CHD' else 1;pct=100*74/735
        summary['diseases'][d]=dict(pgs_id=v['pgs_id'],n_loci=n,observed_loci=n,missing_loci=0,
            observed_sum=.1,imputed_sum=0.,total_sum=.1,eas_percentile=pct,eas_z=0.,
            intervene_group=group(pct),qc_ready=True,high_weight_missing=0,boundary_review_2pp=False)
    write_json(root/'summary.json',summary)
    expected=dict(sample='VAL-45',model=model_hash,gvcf=item['header']['input_signature'][0])
    write_json(root/'inputs.lock.json',expected)
    write_json(root/'score.complete.json',dict(input=expected,files={'summary.json':sha(root/'summary.json')}))
    plan=dict(schema=cohort.SCHEMA,study=str(study),model_sha256=model_hash,samples=[item])
    write_json(study/'plan.json',plan)
    return plan,root


def test_reports_ignore_early_or_not_postverified_scores_and_detect_tampering(tmp_path):
    study=tmp_path/'study';(study/'frozen').mkdir(parents=True)
    plan,root=study_results(study)
    assert list(cohort.completed(plan))==[]
    write_json(root.parent/'state.json',dict(verified_stages=['score'],status='FAILED_ANCESTRY'))
    values=list(cohort.completed(plan))
    assert len(values)==1 and values[0][2] is None
    with (root/'summary.json').open('a') as f:f.write(' ')
    with pytest.raises(ValueError,match='artifact changed'):list(cohort.completed(plan))


def test_standard_report_labels_replicate_files_not_sequencing_runs(tmp_path,runtime,monkeypatch):
    config,_,_=runtime;study=config.PRS_WORK_ROOT/'_validation_20261114/cohort_v1'
    (study/'frozen').mkdir(parents=True);_,root=study_results(study)
    write_json(root.parent/'state.json',dict(verified_stages=['score']))
    output=cohort.report(study)
    result=cohort.read(output/'summary.json')
    assert result['standards_scored_ids']==1 and result['health_scored_ids']==0
    assert all(r['n_replicate_files']==1 and r['sequencing_run_count']==1 for r in result['standard_score_repeatability'])
    assert all(v['confident_eas_quality_ids']==0 and v['eligible_score_distribution'] is None for v in result['population'].values())
    assert 'No GIAB' in result['limitations'][0]


def test_keep_files_read_header_and_reject_duplicate_ids(tmp_path):
    p=tmp_path/'keep';p.write_text('#IID\nE1\nE2\n')
    assert cohort.members(p)=={'E1','E2'}
    p.write_text('#IID\nE1\nE1\n')
    with pytest.raises(ValueError,match='Duplicate'):cohort.members(p)


def test_distribution_distance_handles_ties_without_uniformity_pvalue():
    assert cohort.distance([1,1,2],[1,2,2])==pytest.approx(1/3)
    assert cohort.distance([1,2],[1,2])==0
    assert cohort.describe([1,2],[0,1,2])['p_value'] is None


def test_prepare_freezes_sources_and_reuses_only_unchanged_headers(runtime,monkeypatch):
    config,_,cfg=runtime
    from app.services import prs_jobs
    from research.prs import validation as audit
    raw=Path(cfg['project'])/'raw'
    for n in range(45,55):
        sid='VAL-'+str(n)
        fake_gvcf(raw/'run_a'/'other'/sid/'germline_seq'/(sid+'.gvcf.gz'),sample=sid)
    directory=inventory.write_reports(inventory.scan(raw),Path(cfg['project'])/'inventory')
    model=Path(cfg['multi_model_dir']);model.mkdir()
    write_json(model/'reference_scores.json',{'synthetic':True})
    manifest=dict(sha256={'reference_scores.json':sha(model/'reference_scores.json')})
    write_json(model/'manifest.json',manifest);mhash=sha(model/'manifest.json')
    Path(cfg['fasta']+'.fai').write_text(''.join(f'chr{c}\t{248956422 if c==1 else 100}\t0\t0\t0\n' for c in range(1,23)))
    monkeypatch.setattr(cohort,'settings',lambda _:(cfg,'frozen-config'))
    monkeypatch.setattr(cohort,'load_models',lambda _:(manifest,mhash))
    monkeypatch.setattr(prs_jobs,'resource_signature',lambda _:{'synthetic':'stable'})
    monkeypatch.setattr(audit,'model_audit',lambda *a,**k:dict(status='PASS',model_sha256=mhash))
    study=config.PRS_WORK_ROOT/'_validation_20261114/cohort_v1'
    plan=cohort.prepare(directory,study,Path(cfg['project'])/'cfg.json',Path(cfg['project'])/'scores')
    assert len(plan['samples'])==10
    assert cohort.read(study/'preflight_summary.json')['counts']=={'HEADER_PASS':10}
    assert cohort.prepare(directory,study,Path(cfg['project'])/'cfg.json',Path(cfg['project'])/'scores')==plan
    Path(plan['samples'][0]['gvcf']+'.tbi').write_bytes(b'changed')
    with pytest.raises(ValueError,match='changed'):
        cohort.prepare(directory,study,Path(cfg['project'])/'cfg.json',Path(cfg['project'])/'scores')


def test_command_does_not_replay_old_failure_and_terminates_new_failed_process(tmp_path,monkeypatch):
    detector=SimpleNamespace(nextflow_failed=lambda text:'FAIL_NOW' in text)
    monkeypatch.setattr(cohort,'platform',lambda:(None,None,None,None,None,detector))
    log=tmp_path/'tool.log';log.write_text('FAIL_NOW from an earlier attempt\n')
    with (tmp_path/'lease').open('a') as lease:
        cohort.command([sys.executable,'-c','print("safe restart")'],tmp_path,os.environ.copy(),log,lease)
        with pytest.raises(RuntimeError,match='Nextflow'):
            cohort.command([sys.executable,'-c','import time;print("FAIL_NOW",flush=True);time.sleep(20)'],
                           tmp_path,os.environ.copy(),log,lease)


def test_batch_preserves_verified_score_when_ancestry_fails_then_resumes(runtime,monkeypatch):
    config,_,cfg=runtime;study=config.PRS_WORK_ROOT/'_validation_20261114/cohort_v1'
    (study/'frozen').mkdir(parents=True);plan,root=study_results(study)
    p=fake_gvcf(Path(cfg['project'])/'VAL-45.gvcf.gz',sample='VAL-45')
    item=plan['samples'][0];item['gvcf']=str(p)
    item['header']=dict(status='HEADER_PASS',input_signature=[cohort.fingerprint(p),cohort.fingerprint(str(p)+'.tbi')])
    write_json(study/'plan.json',plan)
    saved=cohort.read(root/'summary.json')
    monkeypatch.setattr(cohort,'verify_current',lambda _:cfg)
    monkeypatch.setattr(cohort.time,'sleep',lambda _:None)
    phases=[]
    def run(args,cwd,env,log,lease,slot=None):
        assert env['NGS_UI_NF_GUARD_SCOPE']=='prs'
        assert slot is not None
        if '--stage' not in args:return
        phase=args[-1];phases.append(phase)
        if phase=='ancestry':raise RuntimeError('synthetic ancestry failure')
        expected=dict(sample='VAL-45',model=plan['model_sha256'],gvcf=cohort.fingerprint(p))
        write_json(root/'inputs.lock.json',expected);write_json(root/'summary.json',saved)
        write_json(root/'score.complete.json',dict(input=expected,files={'summary.json':sha(root/'summary.json')}))
    monkeypatch.setattr(cohort,'command',run)
    with pytest.raises(RuntimeError,match='ancestry failure'):cohort.run_batch(study,'standards','both')
    state=cohort.read(root.parent/'state.json')
    assert state['status']=='FAILED' and state['verified_stages']==['score']
    assert len(list(cohort.completed(plan)))==1
    cohort.run_batch(study,'standards','score')
    assert cohort.read(root.parent/'state.json')['status']=='COMPLETE'
    assert phases==['score','ancestry','score']


@pytest.fixture
def frequency_case(runtime,monkeypatch):
    config,_,cfg=runtime;study=config.PRS_WORK_ROOT/'_validation_20261114/cohort_v1'
    (study/'frozen').mkdir(parents=True);plan,root=study_results(study)
    write_json(root.parent/'state.json',dict(verified_stages=['score']))
    model=Path(cfg['multi_model_dir']);model.mkdir()
    (model/'reference.pgen').write_bytes(b'synthetic-not-native-pgen')
    (model/'reference.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n1\t1\tx\tA\tG\n')
    (model/'reference.psam').write_text('#IID\n'+''.join('E'+str(i)+'\n' for i in range(3330)))
    for n in (735,3330):
        (model/f'reference_calibration_{"eas_735" if n==735 else "3330"}.keep').write_text(
            '#IID\n'+''.join('E'+str(i)+'\n' for i in range(n)))
    (model/'reference.afreq').write_text('#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n1\tx\tA\tG\t0.25\t6660\n')
    db=sqlite3.connect(model/'models.sqlite')
    db.executescript('CREATE TABLE locus(id,ref,alt);CREATE TABLE model(disease,id,ea,weight,freq);')
    db.execute('INSERT INTO locus VALUES(?,?,?)',('x','A','G'))
    summary=cohort.read(root/'summary.json')
    with gzip.open(root/'missing_loci.tsv.gz','wt') as f:
        f.write('disease\tid\tstatus\teffect_weight\timputed_dosage\n')
        for j,(d,r) in enumerate(summary['diseases'].items(),1):
            beta=Decimal(j)/10;ref_effect=j%2==0;p=Decimal('.75' if ref_effect else '.25')
            db.execute('INSERT INTO model VALUES(?,?,?,?,?)',(d,'x','A' if ref_effect else 'G',str(beta),str(p)))
            f.write(f'{d}\tx\tMISSING\t{beta}\t{2*p}\n')
            r.update(n_loci=610305 if d=='CHD' else 2,missing_loci=1,observed_sum=float(beta),
                     imputed_sum=float(2*p*beta),total_sum=float(beta+2*p*beta))
            r['observed_loci']=r['n_loci']-1;r['qc_ready']=r['observed_loci']/r['n_loci']>=.995
            ref=[i/735 for i in range(735)];pct=100*bisect_right(ref,r['total_sum'])/735
            r['eas_percentile']=pct;r['intervene_group']=group(pct)
    db.commit();db.close()
    write_json(root/'summary.json',summary)
    complete=cohort.read(root/'score.complete.json')
    complete['files']={p.name:sha(p) for p in (root/'summary.json',root/'missing_loci.tsv.gz')}
    write_json(root/'score.complete.json',complete)
    monkeypatch.setattr(cohort,'verify_current',lambda _:cfg)
    def fake_command(args,*_):
        assert '--keep' in args and str(model/'reference_calibration_eas_735.keep') in list(map(str,args))
        out=Path(args[args.index('--out')+1])
        Path(str(out)+'.afreq').write_text('#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n1\tx\tA\tG\t0.5\t1470\n')
    monkeypatch.setattr(cohort,'command',fake_command)
    return study,cfg,root


def test_frequency_sensitivity_preserves_effect_allele_direction_and_is_not_truth(frequency_case):
    study,_,_=frequency_case
    out=cohort.frequency_sensitivity(study)
    from research.prs.multidisease.common import rows
    data=list(rows(out/'sensitivity.tsv'))
    for j,r in enumerate(data,1):
        assert float(r['delta_score'])==pytest.approx((-.05 if j%2==0 else .05)*j)
    summary=cohort.read(out/'summary.json')
    assert summary['sample_ids']==1 and 'not genotype-truth' in summary['limitations'][0]
    assert summary['reference_genotype_signatures'].keys()=={'.pgen','.pvar','.psam'}


def test_frequency_generation_rejects_reference_change_during_command(frequency_case,monkeypatch):
    study,cfg,_=frequency_case;original=cohort.command
    def mutation(args,*rest):
        original(args,*rest)
        (Path(cfg['multi_model_dir'])/'reference.pgen').write_bytes(b'changed')
    monkeypatch.setattr(cohort,'command',mutation)
    with pytest.raises(ValueError,match='genotypes changed'):cohort.frequency_sensitivity(study)


def test_reference_mask_pilot_reconciles_unmasked_scores_and_checks_resume_receipts(frequency_case,monkeypatch):
    study,cfg,root=frequency_case
    from research.prs import validation as audit
    from research.prs.multidisease import common
    from research.prs.multidisease.common import rows
    plan=cohort.study_read(study)
    db=sqlite3.connect(Path(cfg['multi_model_dir'])/'models.sqlite');db.row_factory=sqlite3.Row
    models={r['disease']:dict(r) for r in db.execute('SELECT * FROM model')};db.close()
    ref={d:{'E'+str(i):float(Decimal(m['weight'])*(i%3 if m['ea']=='G' else 2-i%3))
            for i in range(735)} for d,m in models.items()}
    write_json(study/'frozen/reference_scores.json',ref)
    write_json(Path(cfg['multi_model_dir'])/'reference_scores.json',ref)
    write_json(study/'frozen/manifest.json',dict(sha256={'reference_scores.json':sha(study/'frozen/reference_scores.json')}))
    plan['model_sha256']=sha(study/'frozen/manifest.json');plan['score_dir']=str(study/'unused-scorefiles')
    write_json(study/'plan.json',plan)
    summary=cohort.read(root/'summary.json');summary['model_manifest_sha256']=plan['model_sha256']
    for d,r in summary['diseases'].items():
        pct=100*bisect_right(sorted(ref[d].values()),r['total_sum'])/735
        r['eas_percentile']=pct;r['intervene_group']=group(pct)
    write_json(root/'summary.json',summary)
    complete=cohort.read(root/'score.complete.json');complete['input']['model']=plan['model_sha256']
    complete['files']['summary.json']=sha(root/'summary.json')
    write_json(root/'inputs.lock.json',complete['input']);write_json(root/'score.complete.json',complete)
    frequency=cohort.frequency_sensitivity(study)
    original=cohort.command
    def native_export(args,*rest):
        if '--export' not in args:return original(args,*rest)
        keep=Path(args[args.index('--keep')+1]);iid=next(iter(cohort.members(keep)))
        dosage=int(iid[1:])%3;gt=['0/0','0/1','1/1'][dosage]
        out=Path(args[args.index('--out')+1])
        with gzip.open(str(out)+'.vcf.gz','wt') as f:
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+iid+'\n')
            f.write('1\t1\tx\tA\tG\t.\tPASS\t.\tGT\t'+gt+'\n')
    monkeypatch.setattr(cohort,'command',native_export)
    monkeypatch.setattr(common,'load_models',lambda _:(dict(models={d:dict(n_loci=1) for d in DISEASES}),plan['model_sha256']))
    monkeypatch.setattr(audit,'source_audit',lambda *a,**k:{'status':'PASS'})
    def retained(*args):
        for d,m in models.items():
            yield d,dict(m,chrom=1,pos=1,ref='A',alt='G',source_ea=m['ea'],
                         source_weight=Decimal(m['weight']),orientation='direct')
    monkeypatch.setattr(audit,'retained_sources',retained)
    output=study/'reference_pilot'
    cohort.reference_benchmark(study,frequency,output,2,1)
    result=cohort.read(output/'summary.json')
    assert result['reference_people']==2 and result['mask_files']==1
    assert all(r['complete_comparisons']==2 for r in result['diseases'].values())
    targets=list(output.glob('*/VAL-45.json'));assert len(targets)==2
    assert all(all(r['match'] for r in cohort.read(p)['reference_calibration'].values()) for p in targets)
    cohort.reference_benchmark(study,frequency,output,2,1)
    with targets[0].open('a') as f:f.write(' ')
    with pytest.raises(ValueError,match='artifact changed'):
        cohort.reference_benchmark(study,frequency,output,2,1)
