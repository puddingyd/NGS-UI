#!/usr/bin/env python3
"""Strictly limited recovery of a post-roundtrip failed v2 raw stage.

No genotype reconstruction, no modification of the vendor script, no cleaning,
no deletion of outputs. Revalidates the validated vendor-generated roundtrip
before resuming pgsc_calc with the *same* Nextflow work directory.
"""
import csv
import os
from pathlib import Path
import subprocess
import shlex

import pipeline as v1
import pipeline_v2 as v2


class EarlyStageFailure(RuntimeError):
    """Raw run did not reach a verifiable roundtrip checkpoint."""


def generated_verify_script(vendor):
    source = Path(vendor).read_text()
    prefix = 'cat > "$OUT/scripts/verify_and_reconcile.py" <<\'PY_VERIFY\'\n'
    start = source.find(prefix)
    if start < 0:
        raise EarlyStageFailure('Vendor verify-script heredoc absent; baseline changed')
    start += len(prefix)
    end = source.find('\nPY_VERIFY\n', start)
    if end < 0:
        raise EarlyStageFailure('Vendor verify-script heredoc terminator absent')
    return (source[start:end] + '\n').encode('utf-8')


def run_logged_env(cmd, logfile, cwd, env):
    """Stream a subprocess to a log while supplying the vendor-required env."""
    logfile = Path(logfile)
    logfile.parent.mkdir(parents=True,exist_ok=True)
    print('COMMAND:',shlex.join(map(str,cmd)),flush=True)
    with logfile.open('a') as log:
        log.write('\nCOMMAND: '+shlex.join(map(str,cmd))+'\n')
        with subprocess.Popen(cmd,cwd=cwd,env=env,stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,text=True,errors='replace') as proc:
            for line in proc.stdout:
                print(line,end='',flush=True)
                log.write(line)
            ret=proc.wait()
    if ret:
        raise RuntimeError(f'Command returned {ret}; see {logfile}')


def verify_existing_roundtrip(raw, sample, gvcf, c):
    raw = Path(raw).resolve()
    gvcf = str(Path(gvcf).resolve())
    required = [
        raw / 'audit/source_sha256.txt',
        raw / 'audit/upstream_reconstruction_summary.json',
        raw / 'audit/all_score_locus_decisions.tsv',
        raw / 'prep/score_ready.vcf.gz',
        raw / 'scripts/verify_and_reconcile.py',
        raw / 'scorefiles/PGS004879_hmPOS_GRCh38.txt.gz',
        raw / 'local_containers.config',
        raw / 'samplesheet.csv',
    ]
    absent = [str(p) for p in required if not p.is_file() or p.stat().st_size == 0]
    if absent:
        raise EarlyStageFailure('Cannot resume pre-roundtrip/incomplete raw stage; new --run-id required. Missing: ' + ', '.join(absent))
    expected = generated_verify_script(c['raw_script'])
    actual = (raw/'scripts/verify_and_reconcile.py').read_bytes()
    if actual != expected:
        raise EarlyStageFailure('Generated verification script differs from frozen vendor version')
    if v1.sha(raw/'local_containers.config') != v1.sha(c['raw_nextflow_config']):
        raise EarlyStageFailure('Raw-stage Nextflow configuration was changed')
    # This reads both source files, including the original large gVCF, and
    # rejects either path replacement or modified content.
    v2.verify_fresh_raw_origin(raw, gvcf, c['original_scorefile'])
    rows = list(csv.DictReader((raw/'samplesheet.csv').open(newline='')))
    if not rows or any(row.get('sampleset') != 'case'+sample or
                       row.get('format') != 'pfile' or
                       row.get('chrom') not in {str(i) for i in range(1,23)}
                       for row in rows):
        raise EarlyStageFailure('Wrong or incomplete samplesheet metadata')
    chroms = set()
    for r in rows:
        pref = Path(r['path_prefix'])
        expected_parent = raw/'roundtrip_22chr/pgen'
        if pref.parent.resolve() != expected_parent.resolve() or pref.name != 'chr'+r['chrom']:
            raise EarlyStageFailure('Samplesheet points outside the expected raw stage')
        if r['chrom'] in chroms:
            raise EarlyStageFailure('Duplicate samplesheet chromosome')
        chroms.add(r['chrom'])
        if any(not Path(str(pref)+'.'+ext).is_file() for ext in ('pgen','pvar','psam')):
            raise EarlyStageFailure('Missing roundtrip PGEN component for '+str(pref))
        if not (raw/f'roundtrip_22chr/export/chr{r["chrom"]}.vcf.gz').is_file():
            raise EarlyStageFailure('Missing exported roundtrip VCF for '+r['chrom'])
    env = os.environ.copy()
    env.update(OUT=str(raw), SAMPLE=sample, SAMPLESET='case'+sample)
    # The verifier compares GT dosage, chr/pos, REF/ALT and PGEN roundtrip.
    run_logged_env(['python3',str(raw/'scripts/verify_and_reconcile.py'),'roundtrip'],
                   raw/'logs/v21_roundtrip_recheck.log',cwd=raw,env=env)
    return env


def resume_raw_after_roundtrip(raw, sample, gvcf, c):
    raw = Path(raw).resolve()
    env = verify_existing_roundtrip(raw,sample,gvcf,c)
    cmd = ['nextflow','run','pgscatalog/pgsc_calc',
           '-r','v2.3.0','-profile','singularity','-resume',
           '-c',str(raw/'local_containers.config'),'-w',str(raw/'nextflow_work'),
           '--input',str(raw/'samplesheet.csv'),
           '--target_build','GRCh38',
           '--scorefile',str(raw/'scorefiles/PGS004879_hmPOS_GRCh38.txt.gz'),
           '--outdir',str(raw/'results'),'--skip_ancestry']
    verification = ['python3',str(raw/'scripts/verify_and_reconcile.py'),'reconcile']
    # If all downstream results already exist, verify before paying for another
    # Nextflow invocation. Otherwise resume the exact same Nextflow work tree.
    try:
        run_logged_env(verification,raw/'logs/v21_reconcile.log',cwd=raw,env=env)
        print('RAW_DOWNSTREAM_ALREADY_COMPLETE: reused after reconciliation',flush=True)
    except RuntimeError:
        print('RAW_PGSC_CALC_RESUME: invoking Nextflow -resume with existing cache',flush=True)
        v1.log_run(cmd,raw/'logs/v21_nextflow_resume.log',cwd=raw)
        run_logged_env(verification,raw/'logs/v21_reconcile.log',cwd=raw,env=env)
    # Same checksum ledger produced by the unmodified vendor script.
    with (raw/'audit/output_sha256.txt').open('w') as f:
        for p in (raw/'prep/score_ready.vcf.gz',raw/'scorefiles/PGS004879_hmPOS_GRCh38.txt.gz'):
            f.write(f'{v1.sha(p)}  {p}\n')
    v1.raw_integrity(raw,sample)
    print('RAW POST-ROUNDTRIP RESUME: PASS',flush=True)
