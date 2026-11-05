"""Exercise the shared reconstruction against real bcftools/samtools/PLINK.

Tiny synthetic genomes deliberately bypass release-size validation ONLY in the
test. Production retains fixed CHD and reference membership requirements.
"""
import gzip
import json
from pathlib import Path
import shutil
import sqlite3
import subprocess

import pytest

from research.prs.multidisease import run as pipeline
from research.prs.multidisease.common import write_json
from research.prs.multidisease.reconstruction import rules, scan_records
from research.prs.absolute_risk.registry import DISEASES


def test_conflicting_same_position_and_non_target_allele(tmp_path):
    rule=rules(tmp_path,dict(fasta='/unused',dbsnp='/unused'))
    targets=[dict(id='x',chrom=1,pos=10,ref='A',alt='G')];bases={('1',10):'A'}
    def record(gt): return rule.source_record(['chr1','10','.','A','G,T','.','PASS','.','GT:DP:GQ',gt+':20:40'])
    assert list(scan_records([record('0/2')],targets,bases,rule))[0][1] is None
    with pytest.raises(RuntimeError,match='Conflicting'):
        list(scan_records([record('0/1'),record('1/1')],targets,bases,rule))


@pytest.mark.skipif(not all(shutil.which(t) for t in ('bcftools','samtools','plink2')),reason='Native genomic tools unavailable')
def test_shared_six_column_score_real_tools_and_resume(tmp_path,monkeypatch):
    model=tmp_path/'model';model.mkdir();root=tmp_path/'run'
    fasta=tmp_path/'ref.fa';fasta.write_text(''.join('>chr'+str(c)+'\n'+'A'*100+'\n' for c in range(1,23)))
    subprocess.run(['samtools','faidx',str(fasta)],check=True)
    gvcf=tmp_path/'synthetic.gvcf'
    with gvcf.open('w') as f:
        f.write('##fileformat=VCFv4.2\n'+''.join(f'##contig=<ID=chr{c},length=100>\n' for c in range(1,23)))
        f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n')
        for key,typ in [('GT','String'),('DP','Integer'),('GQ','Integer'),('MIN_DP','Integer')]:
            f.write(f'##FORMAT=<ID={key},Number=1,Type={typ},Description="Synthetic">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n')
        for c in range(1,23):
            f.write(f'chr{c}\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=9\tGT:DP:GQ:MIN_DP\t0/0:20:40:20\n')
            f.write(f'chr{c}\t10\t.\tA\tG,<NON_REF>\t.\tPASS\t.\tGT:DP:GQ\t0/1:20:40\n')
            f.write(f'chr{c}\t11\t.\tA\t<NON_REF>\t.\tPASS\tEND=25\tGT:DP:GQ:MIN_DP\t0/0:20:40:3\n')
            f.write(f'chr{c}\t30\t.\tA\tG,T,<NON_REF>\t.\tPASS\t.\tGT:DP:GQ\t0/2:20:40\n')
    subprocess.run(['bcftools','view','-Oz','-o',str(gvcf)+'.gz',str(gvcf)],check=True)
    subprocess.run(['bcftools','index','-t',str(gvcf)+'.gz'],check=True)
    db=sqlite3.connect(model/'models.sqlite')
    db.executescript('CREATE TABLE locus(id TEXT PRIMARY KEY,chrom INTEGER,pos INTEGER,ref TEXT,alt TEXT,freq TEXT);'
        'CREATE TABLE model(disease TEXT,id TEXT,ea TEXT,weight TEXT,freq TEXT,row_no INTEGER); CREATE INDEX model_id ON model(id);')
    offsets={d:0. for d in DISEASES}
    with (model/'weights.tsv').open('w') as w,(model/'reference.afreq').open('w') as af:
        w.write('ID\tALT\t'+'\t'.join(DISEASES)+'\n');af.write('#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n')
        for c in range(1,23):
            for i,(pos,alt) in enumerate([(10,'G'),(20,'G'),(30,'T'),(30,'G')]):
                ident=f'{c}:{pos}:A:{alt}'
                db.execute('INSERT INTO locus VALUES(?,?,?,?,?,?)',(ident,c,pos,'A',alt,'0.25'))
                af.write(f'{c}\t{ident}\tA\t{alt}\t0.25\t6660\n')
                weights=[]
                for j,d in enumerate(DISEASES):
                    weight=(j+1)*.001
                    ea='A' if j%2 else alt;freq='0.75' if j%2 else '0.25'
                    db.execute('INSERT INTO model VALUES(?,?,?,?,?,?)',(d,ident,ea,str(weight),freq,4*c+i))
                    weights.append(str(-weight if j%2 else weight))
                    if j%2: offsets[d]+=2*weight
                w.write(ident+'\t'+alt+'\t'+'\t'.join(weights)+'\n')
    db.commit();db.close()
    write_json(model/'reference_scores.json',{d:{str(i):i/735 for i in range(735)} for d in DISEASES})
    manifest=dict(models={d:dict(n_loci=88) for d in DISEASES},offsets=offsets)
    monkeypatch.setattr(pipeline,'load_models',lambda _: (manifest,'a'*64))
    monkeypatch.setattr(pipeline,'score_rows',lambda s:s['diseases'])
    monkeypatch.setattr(pipeline,'plink',lambda cfg,*args:['plink2','--threads','1','--memory','640',*map(str,args)])
    cfg=dict(multi_model_dir=str(model),fasta=str(fasta),dbsnp='/unused',min_dp=10,min_gq=20,palindromic_mode='assume-direct')
    pipeline.run(cfg,'SYNTHETIC',Path(str(gvcf)+'.gz'),root,'score')
    result=json.loads((root/'summary.json').read_text())
    for j,(key,r) in enumerate(result['diseases'].items()):
        assert r['observed_loci']==44 and r['missing_loci']==44
        expected=22*(j+1)*.001*(5 if j%2 else 3)
        assert r['total_sum']==pytest.approx(expected,abs=1e-12)
    times={p.name:p.stat().st_mtime_ns for p in root.glob('chr*.json')}
    pipeline.run(cfg,'SYNTHETIC',Path(str(gvcf)+'.gz'),root,'score')
    assert times=={p.name:p.stat().st_mtime_ns for p in root.glob('chr*.json')}
    with (root/'chr1.audit.tsv.gz').open('ab') as f:f.write(b'tampered')
    # The full score checkpoint validates published files; remove it to exercise
    # failed-run chromosome recovery and reject altered genotype audit artifacts.
    (root/'score.complete.json').unlink()
    with pytest.raises(ValueError,match='artifact changed'):
        pipeline.run(cfg,'SYNTHETIC',Path(str(gvcf)+'.gz'),root,'score')
