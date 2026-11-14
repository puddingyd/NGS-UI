"""Read-only source and independent genotype/mean-fill audits (stdlib).

Never imports production reconstruction/scoring. Inputs are frozen source PGS,
model membership, a separate single-sample genotype VCF, and optional case mask.
Benchmarks without complete truth coverage are explicitly partial.
"""
import argparse
from bisect import bisect_right
import csv
from decimal import Decimal
import gzip
import hashlib
from itertools import zip_longest
import json
import math
from pathlib import Path
import sqlite3
import statistics

from .absolute_risk.registry import DISEASES

ROOT = Path(__file__).resolve().parent


def digest(path, decompress=False):
    h = hashlib.sha256()
    with (gzip.open(path,'rb') if decompress else open(path,'rb')) as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()


def score_rows(path, pgs):
    with gzip.open(path,'rt') as f:
        metadata = []
        for line in f:
            if not line.startswith('#'): break
            metadata.append(line.rstrip())
        if '#pgs_id='+pgs not in metadata:raise ValueError('PGS identity mismatch')
        yield from csv.DictReader(f,fieldnames=line.rstrip().split('\t'),delimiter='\t')


def source_audit(score_dir, chd=None, cancer_original=None):
    pins = {r['pgs_id']:r for r in json.loads((ROOT/'multidisease/sources.json').read_text())}
    result = dict(schema='PRS_SOURCE_AUDIT_V1', status='PASS', diseases={})
    for key,d in DISEASES.items():
        path = Path(chd) if key=='CHD' and chd else Path(score_dir)/(d['pgs_id']+'_hmPOS_GRCh38.txt.gz')
        if not path.is_file():
            result['diseases'][key] = dict(status='WAITING_FOR_DATA');continue
        h = digest(path); pin = pins.get(d['pgs_id'])
        if pin and h != pin['sha256']:raise ValueError('Source checksum mismatch: '+d['pgs_id'])
        lock = json.loads((ROOT/'scoring_v211/method2_model_lock.json').read_text())
        decompressed = digest(path,True) if key=='CHD' else None
        if key=='CHD' and decompressed != lock['original_source_decompressed_sha256']:
            raise ValueError('CHD frozen source differs')
        n=pal=0
        for r in score_rows(path,d['pgs_id']):
            n+=1;pal += {r['effect_allele'].upper(),r['other_allele'].upper()} in ({'A','T'},{'C','G'})
        result['diseases'][key] = dict(status='PASS',pgs_id=d['pgs_id'],sha256=h,
            decompressed_sha256=decompressed,source_rows=n,palindromic=pal,
            scope='Source file only; retained model and genotype truth still require benchmark')
    if cancer_original:
        hm=Path(score_dir)/'PGS004875_hmPOS_GRCh38.txt.gz'
        a=score_rows(cancer_original,'PGS004875');b=score_rows(hm,'PGS004875')
        n=bad=0
        for left,right in zip_longest(a,b):
            n+=1
            bad += left is None or right is None or any(left.get(k)!=right.get(k) for k in ('rsID','effect_allele','other_allele','effect_weight'))
        result['cancer_original_comparison']=dict(rows=n,mismatched_rows=bad,sha256=digest(cancer_original),
            status='PASS' if bad==0 else 'FAIL',scope='Row identity/alleles/weight preservation; not HR direction validation')
        if bad:result['status']='FAIL'
    if any(r['status']=='WAITING_FOR_DATA' for r in result['diseases'].values()) and result['status']=='PASS':result['status']='PARTIAL'
    return result


def genotype_rows(vcf):
    opener=gzip.open if str(vcf).endswith('.gz') else open
    with opener(vcf,'rt') as f:
        header=False
        for line in f:
            if line.startswith('##'):continue
            if line.startswith('#CHROM'):
                if len(line.rstrip().split('\t'))!=10:raise ValueError('Independent VCF must contain one sample')
                header=True;continue
            if not header:raise ValueError('Missing VCF header')
            v=line.rstrip().split('\t')
            if len(v)!=10:raise ValueError('Invalid VCF row')
            if ',' in v[4] or len(v[3])!=1 or len(v[4])!=1:continue
            fields=dict(zip(v[8].split(':'),v[9].split(':')));gt=fields.get('GT','.').replace('|','/').split('/')
            dosage=sum(int(a) for a in gt) if len(gt)==2 and set(gt)<={'0','1'} and v[6] in ('PASS','.') else None
            yield (v[0].removeprefix('chr'),int(v[1]),v[3],v[4]),dosage


def retained_sources(model_dir, score_dir, chd=None):
    db=sqlite3.connect((Path(model_dir)/'models.sqlite').resolve().as_uri()+'?mode=ro',uri=True);db.row_factory=sqlite3.Row
    try:
        for key,d in DISEASES.items():
            models={r['row_no']:dict(r) for r in db.execute('SELECT m.*,l.chrom,l.pos,l.ref,l.alt FROM model m JOIN locus l USING(id) WHERE m.disease=?',(key,))}
            path=Path(chd) if key=='CHD' and chd else Path(score_dir)/(d['pgs_id']+'_hmPOS_GRCh38.txt.gz')
            seen=set()
            for n,raw in enumerate(score_rows(path,d['pgs_id']),1):
                if n not in models:continue
                m=models[n];ea,oa=raw['effect_allele'].upper(),raw['other_allele'].upper()
                if (int(raw['hm_chr']),int(raw['hm_pos'])) != (m['chrom'],m['pos']):raise ValueError('Model/source coordinate mismatch')
                orientation='direct'
                if {ea,oa}!={m['ref'],m['alt']}:
                    ea,oa=ea.translate(str.maketrans('ACGT','TGCA')),oa.translate(str.maketrans('ACGT','TGCA'));orientation='complement'
                    if {ea,oa}!={m['ref'],m['alt']}:raise ValueError('Model/source allele mismatch')
                if ea!=m['ea'] or Decimal(raw['effect_weight'])!=Decimal(m['weight']):raise ValueError('Model/source effect direction or weight mismatch')
                if not Decimal(0)<=Decimal(m['freq'])<=Decimal(1):raise ValueError('Invalid effect frequency')
                seen.add(n)
                yield key,{**m,'source_weight':Decimal(raw['effect_weight']),'source_ea':ea,'orientation':orientation}
            if seen != set(models):raise ValueError('Missing source rows for fixed model')
    finally:db.close()


def stratum(p):
    return sum(p>=b for b in (20,40,60,80,90,95))


def model_audit(model_dir,score_dir,chd=None):
    from .multidisease.common import load_models
    manifest,mhash=load_models(model_dir)
    sources=source_audit(score_dir,chd)
    if sources['status']!='PASS':raise ValueError('All six pinned source files are required')
    counts={key:dict(n=0,ref_effect=0,complement=0,palindromic=0) for key in DISEASES}
    for key,row in retained_sources(model_dir,score_dir,chd):
        r=counts[key];r['n']+=1;r['ref_effect']+=row['source_ea']==row['ref']
        r['complement']+=row['orientation']=='complement'
        r['palindromic']+={row['ref'],row['alt']} in ({'A','T'},{'C','G'})
    for key,r in counts.items():
        if r['n']!=manifest['models'][key]['n_loci']:raise ValueError('Manifest/model count differs')
    return dict(schema='PRS_MODEL_SOURCE_AUDIT_V1',status='PASS',model_sha256=mhash,source_audit=sources,diseases=counts,
        scope='Raw source vs retained model only; no independent genotype or clinical validation')


def aggregate(paths):
    values=[json.loads(Path(p).read_text()) for p in paths]
    if not values or any(v.get('schema')!='PRS_INDEPENDENT_BENCHMARK_V1' for v in values):raise ValueError('Expected benchmark reports')
    if len({v['model_sha256'] for v in values})!=1:raise ValueError('Do not combine model versions')
    identities=[(v['inputs']['vcf'],v['inputs'].get('mask'),v['inputs'].get('eas_frequency')) for v in values]
    if len(set(identities))!=len(identities):raise ValueError('Duplicate benchmark comparisons')
    result={}
    for key in DISEASES:
        full=[v['diseases'][key] for v in values if v['diseases'][key]['status']=='COMPLETE_GENOTYPE_COMPARISON']
        metrics={}
        for scenario in ('current','eas'):
            rows=[r for r in full if scenario in r['percentiles']]
            if rows:
                delta=[r['delta_percentile_pp'][scenario] for r in rows]
                metrics[scenario]=dict(comparisons=len(rows),mean_signed_delta_pp=statistics.mean(delta),
                    median_absolute_delta_pp=statistics.median(map(abs,delta)),max_absolute_delta_pp=max(map(abs,delta)),
                    stratum_changes=sum(r['changed_stratum'][scenario] for r in rows))
        result[key]=dict(complete_comparisons=len(full),partial_comparisons=len(values)-len(full),metrics=metrics,
            numerical_failures=sum(r.get('numerical_match') is False for r in full))
    return dict(schema='PRS_BENCHMARK_AGGREGATE_V1',model_sha256=values[0]['model_sha256'],diseases=result,
        unit='Comparison files, not independent people; repeated truth samples/masks are not population validation',
        input_sha256=[digest(p) for p in paths])


def benchmark(model_dir,score_dir,vcf,chd=None,mask=None,eas_frequency=None,expected=None):
    # Model files identify the retained set; raw source determines beta/alleles.
    from .multidisease.common import load_models
    manifest,mhash=load_models(model_dir)
    audit=source_audit(score_dir,chd)
    if audit['status']!='PASS':raise ValueError('All six pinned source files are required')
    genotypes={}
    for k,d in genotype_rows(vcf):
        if k in genotypes:raise ValueError('Duplicate independent VCF allele')
        genotypes[k]=d
    masks={key:set() for key in DISEASES}
    if mask:
        with gzip.open(mask,'rt') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                if r['disease'] not in masks or r['id'] in masks[r['disease']]:raise ValueError('Invalid/duplicate mask')
                masks[r['disease']].add(r['id'])
    eas={}
    if eas_frequency:
        with open(eas_frequency) as f:
            for r in csv.DictReader(f,delimiter='\t'):
                if r['ID'] in eas:raise ValueError('Duplicate EAS frequency')
                freq=float(r['ALT_FREQS']);obs=int(r['OBS_CT'])
                if not 0<=freq<=1 or not 0<obs<=1470:raise ValueError('Invalid EAS allele frequency/count')
                eas[r['ID']]=(r['REF'],r['ALT'],freq)
    terms={key:dict(full=[],current=[],eas=[],missing_truth=0,masked=0,complements=0,ids=set()) for key in DISEASES}
    for key,m in retained_sources(model_dir,score_dir,chd):
        t=terms[key];t['ids'].add(m['id']);t['complements']+=m['orientation']=='complement'
        gt=genotypes.get((str(m['chrom']),m['pos'],m['ref'],m['alt']))
        if gt is None:t['missing_truth']+=1;continue
        dosage=Decimal(gt if m['source_ea']==m['alt'] else 2-gt);weight=m['source_weight']
        t['full'].append(float(weight*dosage))
        masked=m['id'] in masks[key];t['masked']+=masked
        t['current'].append(float(weight*(2*Decimal(m['freq']) if masked else dosage)))
        if eas_frequency:
            if m['id'] not in eas or eas[m['id']][:2]!=(m['ref'],m['alt']):raise ValueError('EAS frequencies do not cover the fixed allele set')
            p=eas[m['id']][2] if m['source_ea']==m['alt'] else 1-eas[m['id']][2]
            t['eas'].append(float(weight*(Decimal(str(2*p)) if masked else dosage)))
    distributions=json.loads((Path(model_dir)/'reference_scores.json').read_text())
    target=json.loads(Path(expected).read_text()) if expected else None
    results={}
    for key,t in terms.items():
        if masks[key]-t['ids']:raise ValueError('Mask contains non-model loci')
        n=manifest['models'][key]['n_loci'];values=sorted(distributions[key].values())
        if len(t['ids']) != n:raise ValueError('Retained model count differs from manifest')
        if len(values)!=735:raise ValueError('Expected EAS735 distribution')
        sums={name:math.fsum(t[name]) for name in ('full','current','eas') if name!='eas' or eas_frequency}
        complete=t['missing_truth']==0
        row=dict(status='COMPLETE_GENOTYPE_COMPARISON' if complete else 'PARTIAL_TRUTH_COVERAGE',n_loci=n,
            truth_observed=n-t['missing_truth'],truth_missing=t['missing_truth'],masked=t['masked'],complement_rows=t['complements'],sums=sums)
        if complete:
            row['percentiles']={name:100*bisect_right(values,v)/735 for name,v in sums.items()}
            row['delta_percentile_pp']={k:v-row['percentiles']['full'] for k,v in row['percentiles'].items()}
            row['changed_stratum']={k:stratum(v)!=stratum(row['percentiles']['full']) for k,v in row['percentiles'].items()}
            row['delta_score_sd']={k:(v-sums['full'])/statistics.stdev(values) for k,v in sums.items()}
            if target:
                delta=sums['current']-target['diseases'][key]['total_sum'];row['pipeline_delta']=delta
                row['numerical_match']=abs(delta)<=1e-5
        results[key]=row
    failed=any(r.get('numerical_match') is False for r in results.values())
    partial=any(r['status']=='PARTIAL_TRUTH_COVERAGE' for r in results.values())
    status='FAIL' if failed else 'PARTIAL' if partial else 'PASS' if expected else 'CALCULATED_WITHOUT_PIPELINE_COMPARISON'
    return dict(schema='PRS_INDEPENDENT_BENCHMARK_V1',status=status, model_sha256=mhash,source_audit=audit,diseases=results,
        inputs={k:digest(p) for k,p in dict(vcf=vcf,mask=mask,eas_frequency=eas_frequency,expected=expected).items() if p},
        limitations=['Complete truth requires explicit diploid genotypes at every retained locus; absent VCF positions are unknown, not 0/0.',
          'Source/complement membership checks are not independent phenotype validation.',
          'Fixed EAS735 ECDF is used only to isolate masking/frequency sensitivity; no new model is activated.'])


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    for name in ('sources','model','benchmark'):
        p=sub.add_parser(name);p.add_argument('--score-dir',type=Path,required=True);p.add_argument('--chd',type=Path)
        p.add_argument('--output',type=Path,required=True)
        if name=='sources':p.add_argument('--cancer-original',type=Path)
        elif name=='model':p.add_argument('--model-dir',type=Path,required=True)
        else:
            p.add_argument('--model-dir',type=Path,required=True);p.add_argument('--vcf',type=Path,required=True)
            for field in ('mask','eas-frequency','expected'):p.add_argument('--'+field,type=Path)
    p=sub.add_parser('aggregate');p.add_argument('--inputs',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True)
    a=parser.parse_args()
    if a.output.exists():parser.error('Output exists; choose a new audit filename')
    if a.command=='sources':value=source_audit(a.score_dir,a.chd,a.cancer_original)
    elif a.command=='model':value=model_audit(a.model_dir,a.score_dir,a.chd)
    elif a.command=='aggregate':value=aggregate(a.inputs)
    else:value=benchmark(a.model_dir,a.score_dir,a.vcf,a.chd,a.mask,a.eas_frequency,a.expected)
    with a.output.open('x') as f:json.dump(value,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
    if value.get('status')=='FAIL':raise SystemExit(1)


if __name__=='__main__':main()
