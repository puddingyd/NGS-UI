#!/usr/bin/env python3
"""Version-locked canonical PGS004879 scoring from validated per-locus gVCF audit.

No server-only package modules, and no changes to the legacy pipeline. Stdlib only.
This is independent of pgsc_calc's target-dependent MATCH_COMBINE step.
"""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
import os
import re
import statistics
import tempfile
from bisect import bisect_right
from collections import Counter
from decimal import Decimal
from pathlib import Path

EXPECTED = {"score_a_sha256": "358dc695afaee2387b14da0641f23a46b8d7cf33312f6ee14225eb9018d37f96",
            "score_b_sha256": "980541e5f1df2d5c5cd439f1f7020b68ce40a0af403bd6a6e44c02a62e675a42",
            "original_source_decompressed_sha256": "0ee9852635bc831bfbcbb144cb308157bed658e0e26fc4243737b3957aea4d3f"}
MODEL_FILES = [
    "reference_qc_master_v0.tsv.gz",
    "PGS004879_modelA_v0.score.tsv",
    "PGS004879_modelB_v0.score.tsv",
    "reference_qc_unrelated3330.afreq",
    "reference_calibration_3330.keep",
    "reference_calibration_eas_735.keep",
    "reference_modelA_v0.sscore",
    "reference_modelB_v0.sscore",
    "provenance/provenance_issues_v0.tsv",
]
COMP = str.maketrans('ACGT', 'TGCA')
SCHEMA = "PGS004879_METHOD2_V0_20261031"


def die(msg):
    raise RuntimeError(msg)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def gzip_content_sha256(path):
    h = hashlib.sha256()
    with gzip.open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def json_write_atomic(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, prefix='.tmp_', delete=False, encoding='utf8') as f:
        temp = Path(f.name)
        json.dump(data, f, indent=2, ensure_ascii=False, sort_keys=True)
        f.write('\n')
    os.replace(temp, path)


def gz_rows(path):
    with gzip.open(path, 'rt', newline='') as f:
        yield from csv.DictReader(f, delimiter='\t')


def uncompressed_tsv(path):
    with open(path, 'r', newline='') as f:
        yield from csv.DictReader(f, delimiter='\t')


def key(chrom, pos):
    chrom = str(chrom).removeprefix('chr')
    if not chrom.isdigit() or int(chrom) not in range(1, 23):
        die(f'Unexpected autosomal chromosome: {chrom}')
    return chrom, int(pos)


def read_ids(path):
    ids = set()
    for row in uncompressed_tsv(path):
        iid = row.get('#IID') or row.get('IID')
        if not iid or iid in ids: die(f'Missing/duplicate IID in {path}: {iid}')
        ids.add(iid)
    return ids


def sscore(path):
    values = {}
    with open(path) as f:
        header = f.readline().split()
        if '#IID' in header:
            iid_col = header.index('#IID')
        elif 'IID' in header:
            iid_col = header.index('IID')
        else:
            die(f'No IID column in .sscore {path}')
        if 'PGS_SUM' not in header:die(f'No PGS_SUM column in .sscore {path}')
        score_col = header.index('PGS_SUM')
        for line in f:
            fields = line.split()
            if len(fields)!=len(header):die(f'Malformed .sscore row in {path}')
            iid = fields[iid_col]
            if not iid or iid in values:die(f'Duplicate IID in {path}: {iid}')
            value=float(fields[score_col])
            if not math.isfinite(value):die(f'Nonfinite PRS for {iid}')
            values[iid]=value
    return values


def validate_locked_model(model_dir, lock_path, source):
    model_dir = Path(model_dir).resolve()
    lock_path = Path(lock_path).resolve()
    source = Path(source).resolve()
    if not lock_path.is_file():
        die(f'Model lock missing: {lock_path}; run freeze-model AFTER reviewing validated server files')
    lock = json.loads(lock_path.read_text())
    if lock.get('schema') != SCHEMA: die('Model lock schema mismatch')
    if Path(lock['model_dir']).resolve() != model_dir: die('Model directory differs from lock; do not silently redirect')
    if Path(lock['source_path']).resolve() != source: die('Original scorefile path differs from lock')
    if gzip_content_sha256(source) != EXPECTED['original_source_decompressed_sha256']:
        die('Official original scoring data checksum changed')
    for name in MODEL_FILES:
        p = model_dir/name
        if not p.is_file() or sha256(p) != lock['sha256'].get(name):
            die(f'LOCK VIOLATION: file is missing or modified: {p}')
    if sha256(model_dir/'PGS004879_modelA_v0.score.tsv') != EXPECTED['score_a_sha256']:
        die('Model A scorefile differs from validated version')
    if sha256(model_dir/'PGS004879_modelB_v0.score.tsv') != EXPECTED['score_b_sha256']:
        die('Model B scorefile differs from validated version')
    return lock


def complement(a):
    if not a or set(a) - set('ACGT'): die('Invalid allele: '+repr(a))
    return a.translate(COMP)[::-1]


def load_master(path, *, strict=True):
    entries = {}; rows = set(); n_b = 0; flip_count = 0
    for r in gz_rows(path):
        k = key(r['chr'], r['pos'])
        if k in entries: die(f'Duplicate coordinate in canonical master {k}')
        rownum = int(r['source_row_index'])
        if rownum in rows: die(f'Duplicate original source row {rownum}')
        rows.add(rownum)
        if r['qc_strict_1pct'] not in ('PASS', 'REVIEW'):die('Invalid B QC flag')
        n_b += r['qc_strict_1pct'] == 'PASS'
        ea, oa = r['effect_allele'], r['other_allele']
        ref, alt = r['reference_REF'], r['reference_ALT']
        if len({ea, oa, ref, alt}) != 2 or set((ea, oa)) != set((ref, alt)):
            die(f'Malformed effect/ref alleles: {k}')
        orient = 'REF' if ea == ref else 'ALT'
        if r['effect_allele_is'] != orient:die(f'Effect allele orientation mismatch: {k}')
        p = float(r['effect_freq_unrel'])
        if not math.isfinite(p) or not 0 <= p <= 1:die(f'Invalid effect frequency at {k}')
        weight = Decimal(r['effect_weight'])
        entries[k] = {'row':rownum,'id':r['reference_ID'],'ea':ea,'oa':oa,'ref':ref,'alt':alt,
                      'weight':weight,'p':Decimal(r['effect_freq_unrel']),
                      'model_b':r['qc_strict_1pct']=='PASS'}
    if strict and (len(entries)!=610305 or n_b!=603438):
        die(f'Unexpected frozen model size: A={len(entries)}, B={n_b}')
    return entries


def model_check(model_dir, source):
    model_dir = Path(model_dir)
    if sha256(model_dir/'PGS004879_modelA_v0.score.tsv') != EXPECTED['score_a_sha256']:
        die('Known Model A SHA256 mismatch')
    if sha256(model_dir/'PGS004879_modelB_v0.score.tsv') != EXPECTED['score_b_sha256']:
        die('Known Model B SHA256 mismatch')
    if gzip_content_sha256(source) != EXPECTED['original_source_decompressed_sha256']:
        die('Original scoring content does not match official 2026-10-31 audit')
    manifest = load_master(model_dir/'reference_qc_master_v0.tsv.gz')
    # Guard that the master effect-allele frequencies still correspond to the
    # immutable unrelated-3330 allele frequency table.
    by_id={r['id']:r for r in manifest.values()}
    af_seen=set()
    with (model_dir/'reference_qc_unrelated3330.afreq').open() as af:
        header=af.readline().split()
        needed={'ID','REF','ALT','ALT_FREQS','OBS_CT'}
        if not needed <= set(header):die('Reference AFREQ header lacks required fields')
        for line in af:
            fields=line.split()
            r=dict(zip(header,fields))
            vid=r['ID']
            if vid not in by_id or vid in af_seen:die(f'Unexpected or duplicate AFREQ ID {vid}')
            af_seen.add(vid)
            m=by_id[vid]
            if (r['REF'],r['ALT'])!=(m['ref'],m['alt']):die(f'AFREQ allele mismatch {vid}')
            if int(r['OBS_CT'])<=0:die(f'AFREQ missing calls {vid}')
            af_alt=float(r['ALT_FREQS'])
            af_effect=af_alt if m['ea']==m['alt'] else 1-af_alt
            if not math.isfinite(af_effect) or abs(af_effect-float(m['p']))>1e-6:
                die(f'AFREQ effect allele frequency mismatch {vid}')
    if len(af_seen)!=610305:die('Incomplete frequency source')
    for letter in 'AB':
        path = model_dir/f'PGS004879_model{letter}_v0.score.tsv'
        expected = {r['id']:(r['ea'],r['weight']) for r in manifest.values() if letter=='A' or r['model_b']}
        seen = set()
        for r in uncompressed_tsv(path):
            vid = r['ID']
            if vid not in expected or vid in seen or (r['A1'],Decimal(r['PGS'])) != expected[vid]:
                die(f'Model {letter} scoring file disagrees with master: {vid}')
            seen.add(vid)
        if len(seen)!=len(expected): die(f'Model {letter} has missing score rows')
    ref_ids = read_ids(model_dir/'reference_calibration_3330.keep')
    eas_ids = read_ids(model_dir/'reference_calibration_eas_735.keep')
    if len(ref_ids)!=3330 or len(eas_ids)!=735 or not eas_ids <= ref_ids:
        die('Unexpected reference population membership')
    for model in 'AB':
        scores = sscore(model_dir/f'reference_model{model}_v0.sscore')
        if set(scores)!=ref_ids: die(f'Unmatched reference IID set, Model {model}')
    issues = list(uncompressed_tsv(model_dir/'provenance/provenance_issues_v0.tsv'))
    if issues: die(f'Unresolved source-to-canonical provenance issues: {len(issues)}')
    return manifest


def freeze_model(model_dir, lock_path, source, project=None):
    model_dir = Path(model_dir).resolve(); source = Path(source).resolve()
    lock_path = Path(lock_path).resolve()
    if lock_path.exists(): die(f'Lock already exists, refusing overwrite: {lock_path}')
    if project is not None:
        inventory = model_dir/'provenance/method2_v0_SHA256SUMS.txt'
        if not inventory.is_file(): die('Missing historical SHA256 inventory')
        import subprocess
        cp = subprocess.run(['sha256sum','-c',str(inventory)],cwd=project,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        if cp.returncode: die('Historical server inventory failed: '+cp.stdout+' '+cp.stderr)
        print('HISTORICAL SHA256 INVENTORY: PASS')
    model_check(model_dir, source)
    lock = {'schema':SCHEMA,'created_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
            'model_dir':str(model_dir),'source_path':str(source),
            'original_source_decompressed_sha256':EXPECTED['original_source_decompressed_sha256'],
            'variants_A':610305,'variants_B':603438,'unrelated_n':3330,'eas_n':735,
            'imputation':'2*p_effect_unrelated_3330',
            'percentile':'ECDF_100_times_count_reference_LE_target_over_735',
            'sha256':{name:sha256(model_dir/name) for name in MODEL_FILES}}
    json_write_atomic(lock_path,lock)
    print('MODEL FREEZE: PASS; lock:', lock_path)
    return lock


def raw_allele_check(r,m,k):
    # The validated v1 audit refers to the original scorefile orientation, not the corrected canonical one.
    e, o = r['EFFECT_ALLELE'].upper(), r['OTHER_ALLELE'].upper()
    if r['WEIGHT'] and Decimal(r['WEIGHT']) != m['weight']:
        die(f'Source vs canonical effect weight mismatch: {k}')
    if (e,o) == (m['ea'],m['oa']):return 'DIRECT_METADATA'
    if (complement(e),complement(o)) == (m['ea'],m['oa']):return 'REVERSE_COMPLEMENT_METADATA'
    die(f'Source to canonical allele orientation mismatch: {k}: {e}/{o} -> {m["ea"]}/{m["oa"]}')


def score_audit(audit_path, master, reference, eas_ids, *, strict=True):
    seen=set(); statuses=Counter(); flip_rows=0; metadata_orientation=Counter()
    obs_terms={'A':[], 'B':[]}; imp_terms={'A':[], 'B':[]}
    obs_ct=Counter(); missing_ct=Counter(); missing_weight={'A':Decimal(0),'B':Decimal(0)}
    missing_records=[]; valid_count=0
    strand_warnings=Counter()
    excluded_accepted_count=0; excluded_accepted_sum=Decimal(0); canonical_accepted_sum=Decimal(0)
    for r in uncompressed_tsv(audit_path):
        k=key(r['CHROM'],r['POS'])
        if k not in master:
            # Original raw v1 can score reference-incompatible loci. Preserve
            # their contributions in a separate raw-to-canonical reconciliation.
            if r['STATUS'] in ('CANDIDATE_DIRECT','FLIP_CONFIRMED'):
                exdosage=Decimal(r['RECALC_EFFECT_DOSAGE'])
                exweight=Decimal(r['WEIGHT'])
                excon=Decimal(r['RECALC_CONTRIBUTION'])
                if exdosage*exweight != excon:
                    die(f'Excluded raw locus contribution mismatch at {k}')
                excluded_accepted_count+=1
                excluded_accepted_sum+=excon
            continue
        if k in seen:die(f'Duplicate audit coordinate {k}')
        seen.add(k); m=master[k]
        orientation=raw_allele_check(r,m,k)
        metadata_orientation[orientation]+=1
        if orientation=='REVERSE_COMPLEMENT_METADATA':flip_rows+=1
        status=r['STATUS'];statuses[status]+=1
        if r.get('STRAND_AMBIGUITY','').strip():
            strand_warnings[r['STRAND_AMBIGUITY']]+=1
        d=r.get('RECALC_EFFECT_DOSAGE','').strip()
        pass_qc = status in ('CANDIDATE_DIRECT','FLIP_CONFIRMED')
        if pass_qc:
            if d not in ('0','1','2','0.0','1.0','2.0'):
                die(f'Invalid observed genotype dosage at {k}: {d}')
            dosage=Decimal(d)
            # Independent reconstruction from v1 VCF GT, REF and ALT, not just
            # the precomputed contribution. This catches orientation mistakes.
            vcf_ref=r.get('VCF_REF','').upper()
            vcf_alt=r.get('VCF_ALT','').upper()
            if {vcf_ref,vcf_alt} != {m['ref'],m['alt']}:
                die(f'Raw VCF and canonical reference alleles differ: {k}: {vcf_ref}/{vcf_alt} vs {m["ref"]}/{m["alt"]}')
            gtparts=re.split(r'[/|]',r.get('VCF_GT',''))
            if len(gtparts)!=2 or any(x not in ('0','1') for x in gtparts):
                die(f'Invalid accepted biallelic GT at {k}: {r.get("VCF_GT")}')
            altcount=sum(x=='1' for x in gtparts)
            if r.get('ALT_DOSAGE','') and int(r['ALT_DOSAGE'])!=altcount:
                die(f'Raw audit ALT dosage disagrees with GT at {k}')
            recalc=altcount if m['ea']==vcf_alt else 2-altcount
            if Decimal(recalc)!=dosage:
                die(f'Effect allele dosage disagrees with canonical GT at {k}: {recalc} vs {dosage}')
            if orientation=='REVERSE_COMPLEMENT_METADATA' and status!='FLIP_CONFIRMED':
                die(f'Accepted reverse-complement allele did not pass dbSNP flip QC: {k}')
            contribution = dosage*m['weight']
            observed_text = r.get('RECALC_CONTRIBUTION','').strip()
            if not observed_text or Decimal(observed_text) != contribution:
                die(f'Audit effect-weight contribution mismatch: {k}')
            # source scorefile effect allele can be complemented to canonical allele but dosage count is unchanged
            valid_count+=1
            canonical_accepted_sum+=contribution
            obs_terms['A'].append(float(contribution));obs_ct['A']+=1
            if m['model_b']:
                obs_terms['B'].append(float(contribution));obs_ct['B']+=1
        else:
            if d not in ('','.'):die(f'Failed-QC locus has dosage {k} {status} {d}')
            contribution = (Decimal(2)*m['p'])*m['weight']
            imp_terms['A'].append(float(contribution));missing_ct['A']+=1;missing_weight['A']+=abs(m['weight'])
            missing_records.append((m['row'],m['id'],status,str(m['weight']),str(m['p']),str(Decimal(2)*m['p']),str(contribution),'A'))
            if m['model_b']:
                imp_terms['B'].append(float(contribution));missing_ct['B']+=1;missing_weight['B']+=abs(m['weight'])
                missing_records.append((m['row'],m['id'],status,str(m['weight']),str(m['p']),str(Decimal(2)*m['p']),str(contribution),'B'))
    if seen != set(master):die(f'Canonical audit coverage mismatch: seen {len(seen)} expected {len(master)}, missing {len(set(master)-seen)}')
    if strict and (len(master)!=610305 or flip_rows!=331):die(f'Canonical master orientation audit mismatch: {len(master)} {flip_rows}')
    out=[]
    for model in 'AB':
        ref=reference[model]
        vals=sorted(ref[i] for i in eas_ids)
        if len(vals)!=735 and strict:die('EAS score set size mismatch')
        observed_sum=math.fsum(obs_terms[model]); imputed_sum=math.fsum(imp_terms[model])
        total_sum=math.fsum([observed_sum,imputed_sum])
        if not math.isfinite(total_sum):die('Nonfinite resulting SUM')
        percentile=100*bisect_right(vals,total_sum)/len(vals)
        z=(total_sum-statistics.mean(vals))/statistics.stdev(vals)
        limits=[(20,'<20'),(40,'20-40'),(60,'40-60'),(80,'60-80'),(90,'80-90'),(95,'90-95')]
        group=next((name for bound,name in limits if percentile<bound),'>=95')
        boundary_distance=min(abs(percentile-v) for v in (20,40,60,80,90,95))
        out.append({'model':model,'n_loci':len(master) if model=='A' else sum(x['model_b'] for x in master.values()),
                    'observed_loci':obs_ct[model],'missing_loci':missing_ct[model],
                    'observed_sum':observed_sum,'imputed_sum':imputed_sum,'total_sum':total_sum,
                    'missing_absolute_effect_weight':float(missing_weight[model]),
                    'eas_mean':statistics.mean(vals),'eas_sd_sample':statistics.stdev(vals),
                    'eas_z':z,'eas_percentile':percentile,'intervene_group':group,
                    'boundary_distance_pp':boundary_distance,
                    'boundary_review_2pp':boundary_distance<=2.0})
        if obs_ct[model]+missing_ct[model]!=out[-1]['n_loci']:
            die('Observed/missing total inconsistency')
    if strict and obs_ct['A']/610305 < 0.995:
        die(f'Target observed fraction {obs_ct["A"]/610305:.3%} below Method2 research guard 99.5%')
    return {'scores':out,'raw_reconciliation':{'canonical_observed_exact':str(canonical_accepted_sum), 'outside_canonical_observed_exact':str(excluded_accepted_sum), 'outside_canonical_observed_n':excluded_accepted_count, 'reconstructed_raw_sum':str(canonical_accepted_sum+excluded_accepted_sum)},'qc_status_counts':dict(statuses),'metadata_orientation_counts':dict(metadata_orientation),
            'missing':missing_records,'canonical_audit_rows':len(seen),'observed':valid_count,'strand_ambiguity_warnings':dict(strand_warnings)}


def write_outputs(out_dir, sample, audit, results, model_lock, ancestry_qc=None):
    out_dir=Path(out_dir); out_dir.mkdir(parents=True,exist_ok=True)
    summary={'schema':SCHEMA,'sample':sample,'stage':'METHOD2_TECHNICAL_PASS',
             'warning':'Research-use scoring, not CHD clinical prediction/absolute risk; sample genotype QC from v1 reconstruction',
             'model_lock_sha256':sha256(model_lock),
             'ancestry_status': ancestry_qc if ancestry_qc is not None else {'status':'NOT_ASSESSED','note':'EAS percentile is a conditional comparison, not a verified target ancestry assignment'},
             'audit_input':str(Path(audit).resolve()), 'audit_sha256':sha256(audit),
             'canonical_audit_rows':results['canonical_audit_rows'],
             'observed_loci':results['observed'],
             'raw_reconciliation':results['raw_reconciliation'],
             'high_weight_missing_modelA':sorted(({
                 'reference_ID':r[1],'status':r[2],'effect_weight':r[3],
                 'imputed_contribution':r[6]
             } for r in results['missing'] if r[7]=='A' and abs(float(r[3]))>=0.001),
                 key=lambda r:abs(float(r['effect_weight'])),reverse=True),
             'qc_status_counts':results['qc_status_counts'],
             'strand_ambiguity_warnings':results['strand_ambiguity_warnings'],
             'metadata_orientation_counts':results['metadata_orientation_counts'],
             'scores':results['scores']}
    out_file=out_dir/'summary.json'; json_write_atomic(out_file,summary)
    with (out_dir/'scores.tsv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=['sample']+list(results['scores'][0]),delimiter='\t')
        wr.writeheader()
        for r in results['scores']:wr.writerow({'sample':sample,**r})
    # Deterministic gzip avoids content-hash drift across repeated runs.
    with (out_dir/'missing_loci.tsv.gz').open('wb') as fileobj:
        with gzip.GzipFile(fileobj=fileobj,mode='wb',filename='',mtime=0) as gz:
            with io.TextIOWrapper(gz,encoding='utf-8',newline='') as text:
                wr=csv.writer(text,delimiter='\t')
                wr.writerow(['source_row_index','reference_ID','qc_status','effect_weight','effect_freq_unrelated3330','imputed_dosage','imputed_contribution','model'])
                wr.writerows(results['missing'])
    report=['# PGS004879 Method 2 fixed-universe PRS','',f'- Sample: `{sample}`',
            f'- Canonical audit: {results["canonical_audit_rows"]} loci',
            '- Model version: Method 2 v0; unrelated 3330, EAS 735',
            '- Missing: 2 × unrelated3330 effect-allele frequency; no QC failed locus treated as 0/0','',
            '| Model | Observed / total | Observed SUM | Imputed SUM | Total SUM | EAS Z | EAS percentile | Group |',
            '|---|---:|---:|---:|---:|---:|---:|---|']
    for r in results['scores']:
        report.append(f'| {r["model"]} | {r["observed_loci"]:,} / {r["n_loci"]:,} | {r["observed_sum"]:.9f} | {r["imputed_sum"]:+.9f} | {r["total_sum"]:.9f} | {r["eas_z"]:+.6f} | {r["eas_percentile"]:.6f} | {r["intervene_group"]} |')
    if ancestry_qc is None:
        report.append('**Ancestry not assessed in this run:** EAS percentiles are reference comparisons only; validate target ancestry separately before interpretation.')
    elif ancestry_qc.get('ancestry')!='EAS' or str(ancestry_qc.get('ancestry_low_confidence'))!='False':
        report.append('**Ancestry review:** Target is not confidently assigned EAS; EAS percentile should NOT be interpreted as ancestry-matched.')
    if any(r['boundary_review_2pp'] for r in results['scores']):
        report.append('**Boundary review:** At least one percentile is within 2 percentage points of a reporting cutoff; display continuous percentile and QC flags.')
    if summary['high_weight_missing_modelA']:
        report.append('**High-weight missing (|beta| >= 0.001):** ' + '; '.join(
            f"{r['reference_ID']} ({r['status']}, beta={r['effect_weight']})"
            for r in summary['high_weight_missing_modelA'][:10]
        ))
    if results['strand_ambiguity_warnings']:
        report.append('**Strand ambiguity:** v1 raw reconstruction flagged palindromic allele orientation issues; review QC metadata and do not regard these as independently phased.')
    report+=['','**Limitations:** The scorer reuses validated v1 gVCF genotype reconstruction; it is not independent genotype calling. The reference percentile is empirical EAS and is not a CHD event risk. Clinical calibration remains incomplete.','']
    (out_dir/'report.md').write_text('\n'.join(report),encoding='utf-8')
    return summary


def score(sample, raw_dir, model_dir, lock_path, source, out_dir, ancestry_qc_path=None):
    from pipeline import raw_integrity
    raw_dir=Path(raw_dir).resolve()
    raw_integrity(raw_dir,sample)
    validate_locked_model(model_dir,lock_path,source)
    audit=raw_dir/'audit/all_score_locus_decisions.tsv'
    if not audit.is_file():die(f'No per-locus audit: {audit}')
    manifest=load_master(Path(model_dir)/'reference_qc_master_v0.tsv.gz')
    eas=read_ids(Path(model_dir)/'reference_calibration_eas_735.keep')
    ref={m:sscore(Path(model_dir)/f'reference_model{m}_v0.sscore') for m in 'AB'}
    if any(len(s)!=3330 for s in ref.values()) or any(not eas <= s.keys() for s in ref.values()):
        die('Reference sample membership corrupted')
    result=score_audit(audit,manifest,ref,eas)
    raw_json=json.loads((raw_dir/'audit/final_reconciliation.json').read_text())
    canonical_raw=Decimal(result['raw_reconciliation']['reconstructed_raw_sum'])
    historical_raw=Decimal(raw_json['observed_sum'])
    if abs(canonical_raw-historical_raw)>Decimal('0.000002'):
        die(f'Independent raw-to-canonical observed SUM mismatch: {canonical_raw} vs {historical_raw}')
    ancestry_qc=None
    if ancestry_qc_path is not None:
        qpath=Path(ancestry_qc_path)
        if not qpath.is_file():die(f'Ancestry QC file missing: {qpath}')
        q=json.loads(qpath.read_text())
        if q.get('sample')!=sample:die('Ancestry QC sample does not match target')
        ancestry_qc={k:q.get(k) for k in ('ancestry','ancestry_low_confidence','rf_p_eas','qc_status')}
    summary=write_outputs(out_dir,sample,audit,result,lock_path,ancestry_qc)
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    fr=sub.add_parser('freeze-model',help='Validate server model assets and create immutable SHA256 model lock')
    fr.add_argument('--model-dir',required=True);fr.add_argument('--lock',required=True);fr.add_argument('--source',required=True)
    fr.add_argument('--project',help='Verify existing historical 17-file sha256sum inventory (recommended)')
    sc=sub.add_parser('score',help='Score prevalidated raw gVCF reconstruction')
    sc.add_argument('--sample',required=True);sc.add_argument('--raw-dir',required=True);sc.add_argument('--model-dir',required=True)
    sc.add_argument('--lock',required=True);sc.add_argument('--source',required=True);sc.add_argument('--out',required=True)
    sc.add_argument('--ancestry-qc',help='Previously generated v1 qc/qc.json to confirm target ancestry')
    a=p.parse_args()
    try:
        if a.command=='freeze-model':freeze_model(a.model_dir,a.lock,a.source,a.project)
        else:
            summary=score(a.sample,a.raw_dir,a.model_dir,a.lock,a.source,a.out,a.ancestry_qc)
            for r in summary['scores']:
                print(f"METHOD2 {a.sample} Model {r['model']}: SUM={r['total_sum']:.12f} Z={r['eas_z']:+.6f} EAS={r['eas_percentile']:.6f} GROUP={r['intervene_group']}")
            print('METHOD2 SCORING: TECHNICAL PASS')
    except Exception as exc:
        p.exit(1, f'ERROR: {exc}\n')


if __name__=='__main__':main()
