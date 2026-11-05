"""One gVCF reconstruction, one six-column scoring job and one ancestry run."""
import argparse
from bisect import bisect_right
import csv
from decimal import Decimal
import gzip
import json
import math
from pathlib import Path
import sqlite3
import statistics

from ..absolute_risk.engine import group
from ..absolute_risk.registry import DISEASES
from ..absolute_risk.multi import score_rows
from .common import command, load_models, plink, rows, score_values, sha, write_json
from .reconstruction import rules, scan_chromosome
from .ancestry_input import stage_score, preflight


def fingerprint(path):
    p = Path(path); s = p.stat()
    return [str(p.resolve()), s.st_size, s.st_mtime_ns, s.st_ctime_ns]


def checkpoint(path, expected):
    if not path.is_file(): return False
    saved = json.loads(path.read_text())
    if saved.get('input') != expected: raise ValueError('PRS checkpoint input changed')
    for name, digest in saved['files'].items():
        if sha(path.parent/name) != digest: raise ValueError('PRS checkpoint artifact changed')
    return True


def reconstruct(cfg, sample, gvcf, root, db, expected):
    rule = rules(root, cfg)
    for chrom in range(1, 23):
        targets = [dict(r) for r in db.execute('SELECT * FROM locus WHERE chrom=? ORDER BY pos,ref,alt', (chrom,))]
        if not targets: raise ValueError('Missing chromosome in model union')
        marker = root/f'chr{chrom}.json'
        if checkpoint(marker, expected):
            print(f'[PRS] chr{chrom} 共用重建：沿用已核對結果', flush=True); continue
        audit = root/f'chr{chrom}.audit.tsv.gz'; vcf = root/f'chr{chrom}.vcf'
        with gzip.open(str(audit)+'.tmp', 'wt') as a, Path(str(vcf)+'.tmp').open('w') as f:
            w = csv.writer(a, delimiter='\t'); w.writerow(['id', 'dosage', 'status'])
            f.write('##fileformat=VCFv4.2\n'+''.join(f'##contig=<ID={c}>\n' for c in range(1,23)))
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n')
            for target, dosage, status in scan_chromosome(gvcf, chrom, targets, rule):
                w.writerow([target['id'], '' if dosage is None else dosage, status])
                gt = {None: './.', 0: '0/0', 1: '0/1', 2: '1/1'}[dosage]
                f.write(f"{chrom}\t{target['pos']}\t{target['id']}\t{target['ref']}\t{target['alt']}\t.\tPASS\t.\tGT\t{gt}\n")
        Path(str(audit)+'.tmp').replace(audit); Path(str(vcf)+'.tmp').replace(vcf)
        command(['bcftools', 'view', '-Oz', '-o', str(vcf)+'.gz', vcf])
        command(['bcftools', 'index', '-f', '-t', str(vcf)+'.gz'])
        files = [audit, Path(str(vcf)+'.gz'), Path(str(vcf)+'.gz.tbi')]
        write_json(marker, dict(input=expected, files={p.name: sha(p) for p in files}))
        vcf.unlink()
        print(f'[PRS] chr{chrom} 共用重建完成：{len(targets):,} 個位點', flush=True)


def roundtrip(root, db):
    expected = {}
    for chrom in range(1, 23):
        for r in rows(root/f'chr{chrom}.audit.tsv.gz'):
            if r['id'] in expected: raise ValueError('Duplicate shared audit locus')
            expected[r['id']] = None if r['dosage'] == '' else int(r['dosage'])
    seen = set()
    with gzip.open(root/'roundtrip.vcf.gz', 'rt') as stream:
        for line in stream:
            if line.startswith('#'): continue
            fields = line.rstrip().split('\t'); ident = fields[2]
            if len(fields) != 10 or ident not in expected or ident in seen:
                raise ValueError('Roundtrip locus identity differs')
            gt = fields[9].split(':')[fields[8].split(':').index('GT')].replace('|', '/').split('/')
            dosage = None if gt == ['.', '.'] else gt.count('1') if len(gt) == 2 and set(gt) <= {'0', '1'} else -1
            if dosage != expected[ident]: raise ValueError('Roundtrip genotype differs: '+ident)
            allele = db.execute('SELECT chrom,pos,ref,alt FROM locus WHERE id=?', (ident,)).fetchone()
            if tuple(allele) != (int(fields[0].removeprefix('chr')), int(fields[1]), fields[3], fields[4]):
                raise ValueError('Roundtrip reference allele differs')
            seen.add(ident)
    if seen != set(expected): raise ValueError('Roundtrip lost loci')


def summarize(db, root, manifest, model_dir, sample, model_hash):
    terms = {d: {'observed': [], 'imputed': [], 'high_weight_missing': 0} for d in DISEASES}
    missing_path = root/'missing_loci.tsv.gz'
    with gzip.open(missing_path, 'wt') as stream:
        writer = csv.writer(stream, delimiter='\t'); writer.writerow(['disease', 'id', 'status', 'effect_weight', 'imputed_dosage'])
        for chrom in range(1, 23):
            for r in rows(root/f'chr{chrom}.audit.tsv.gz'):
                dosage = None if r['dosage'] == '' else int(r['dosage'])
                for m in db.execute('SELECT m.*,l.alt FROM model m JOIN locus l USING(id) WHERE m.id=?', (r['id'],)):
                    observed = dosage is not None
                    effect = Decimal(dosage if m['ea'] == m['alt'] else 2-dosage) if observed else 2*Decimal(m['freq'])
                    weight = Decimal(m['weight']); key = 'observed' if observed else 'imputed'
                    terms[m['disease']][key].append(float(effect*weight))
                    if not observed:
                        writer.writerow([m['disease'], m['id'], r['status'], str(weight), str(effect)])
                        terms[m['disease']]['high_weight_missing'] += abs(weight) >= Decimal('.001')
    distributions = json.loads((model_dir/'reference_scores.json').read_text())
    plink_result = score_values(root/'scored.sscore', DISEASES)
    if set(plink_result) != {sample}: raise ValueError('Scoring sample identity mismatch')
    results = {}
    for key, t in terms.items():
        n = manifest['models'][key]['n_loci']; observed = len(t['observed']); missing = len(t['imputed'])
        if observed+missing != n: raise ValueError('Incomplete per-disease genotype coverage')
        a, b = math.fsum(t['observed']), math.fsum(t['imputed']); total = math.fsum([a, b])
        check = plink_result[sample][key]+float(manifest['offsets'][key])
        if abs(check-total) > 1e-5: raise ValueError(f'{key}: independent PLINK score reconciliation failed')
        reference = sorted(distributions[key].values())
        if len(reference) != 735 or len(set(reference)) < 2 or not all(math.isfinite(x) for x in reference):
            raise ValueError('Invalid fixed reference distribution')
        pct = 100*bisect_right(reference, total)/735
        results[key] = dict(pgs_id=DISEASES[key]['pgs_id'], n_loci=n, observed_loci=observed, missing_loci=missing,
            observed_sum=a, imputed_sum=b, total_sum=total, eas_percentile=pct,
            eas_z=(total-statistics.mean(reference))/statistics.stdev(reference), intervene_group=group(pct),
            qc_ready=observed/n >= .995, high_weight_missing=t['high_weight_missing'],
            boundary_review_2pp=min(abs(pct-x) for x in (20, 40, 60, 80, 90, 95)) <= 2)
        print(f'[PRS] {DISEASES[key]["name"]}：score={total:.6f}；percentile={pct:.2f}；可靠位點={observed}/{n}', flush=True)
    summary = dict(schema='PRS_MULTI_SCORE_V1', sample=sample, reference_id='unrelated3330_EAS735_v1',
                   model_manifest_sha256=model_hash, diseases=results)
    score_rows(summary)
    write_json(root/'summary.json', summary)
    with (root/'scores.tsv').open('w') as f:
        f.write('disease\tpgs_id\tscore\tpercentile\tqc_ready\n')
        for d, r in results.items(): f.write(f"{d}\t{r['pgs_id']}\t{r['total_sum']}\t{r['eas_percentile']}\t{r['qc_ready']}\n")


def ancestry(cfg, sample, root, model_dir):
    dest = root/'ancestry'; dest.mkdir(exist_ok=True)
    # Old installed manifests remain valid. Never rewrite their auxiliary file
    # or the calibrated scores just to complete required PGS header metadata.
    score = stage_score(model_dir/'ancestry_score.txt.gz', dest/'ancestry_score.txt.gz', 610305)
    preflight(cfg, score)
    command(plink(cfg, '--pfile', root/'union', '--extract', model_dir/'ancestry.ids', '--geno', '0',
                  '--make-pgen', '--out', dest/'target'))
    sheet = dest/'samplesheet.csv'
    with sheet.open('w') as f:
        w = csv.writer(f); w.writerow(['sampleset', 'path_prefix', 'chrom', 'format'])
        w.writerow(['case'+sample, str(dest/'target'), 'ALL', 'pfile'])
    command(['nextflow', 'run', 'pgscatalog/pgsc_calc', '-r', 'v2.3.0', '-profile', 'singularity',
        '-c', cfg['ancestry_nextflow_config'], '-resume', '-w', dest/'work', '--input', sheet,
        '--target_build', 'GRCh38', '--scorefile', score,
        '--run_ancestry', cfg['reference_panel'], '--projection_method', cfg['projection_method'],
        '--ref_label', cfg['reference_label'], '--pca_maf_target', str(cfg['pca_maf_target']),
        '--pca_geno_miss_target', str(cfg['pca_geno_miss_target']), '--max_cpus', str(cfg['max_cpus']),
        '--max_memory', cfg['max_memory'], '--outdir', dest/'results'], cwd=dest)
    pop = list(rows(dest/'results'/('case'+sample)/'score'/f'case{sample}_popsimilarity.txt.gz'))
    selected = [r for r in pop if r.get('sampleset') == 'case'+sample and r.get('IID') == sample]
    if len(selected) != 1: raise ValueError('Ancestry sample identity mismatch')
    r = selected[0]
    files = list((dest/'work').rglob('*thinned.prune.in.gz'))
    counts = set()
    for f in files:
        with gzip.open(f, 'rt') as s: counts.add(sum(bool(line.strip()) for line in s))
    if len(counts) != 1 or min(counts) <= 0: raise ValueError('Missing or inconsistent PCA marker count')
    qc = dict(sample=sample, qc_status='PASS', ancestry=r['MostSimilarPop'],
        ancestry_low_confidence=r['MostSimilarPop_LowConfidence'], pca_ld_pruned_markers=min(counts),
        marker_policy='fixed_CHD_canonical_v0', ancestry_input_sha256=sha(score),
        formatter_preflight='PASS', flags=[])
    if qc['ancestry'] != 'EAS' or qc['ancestry_low_confidence'] != 'False':
        qc.update(qc_status='PASS_WITH_REVIEW_FLAGS', flags=['祖源不符合東亞風險估算條件'])
    write_json(root/'qc.json', qc)


def run(cfg, sample, gvcf, root, stage):
    if cfg.get('min_dp') != 10 or cfg.get('min_gq') != 20 or cfg.get('palindromic_mode') != 'assume-direct':
        raise ValueError('Six-disease release requires DP10/GQ20 and assume-direct palindrome policy')
    model_dir = Path(cfg['multi_model_dir'])
    manifest, model_hash = load_models(model_dir)
    root.mkdir(parents=True, exist_ok=True)
    expected = dict(sample=sample, gvcf=fingerprint(gvcf), model=model_hash,
                    code={p.name: sha(p) for p in Path(__file__).parent.glob('*.py')},
                    fasta=fingerprint(cfg['fasta']), config=cfg)
    lock = root/'inputs.lock.json'
    if lock.exists() and json.loads(lock.read_text()) != expected: raise ValueError('PRS inputs changed; start a new run')
    write_json(lock, expected)
    db = sqlite3.connect((model_dir/'models.sqlite').resolve().as_uri()+'?mode=ro', uri=True); db.row_factory = sqlite3.Row
    try:
        if stage == 'score':
            if checkpoint(root/'score.complete.json', expected): return
            reconstruct(cfg, sample, gvcf, root, db, expected)
            command(['bcftools', 'concat', '-Oz', '-o', root/'union.vcf.gz', *[root/f'chr{c}.vcf.gz' for c in range(1, 23)]])
            command(plink(cfg, '--vcf', root/'union.vcf.gz', '--make-pgen', '--out', root/'union'))
            command(plink(cfg, '--pfile', root/'union', '--export', 'vcf', 'bgz', 'id-paste=iid', '--out', root/'roundtrip'))
            roundtrip(root, db)
            command(plink(cfg, '--pfile', root/'union', '--read-freq', model_dir/'reference.afreq',
                '--score', model_dir/'weights.tsv', '1', '2', 'header-read', 'cols=+scoresums',
                '--score-col-nums', '3-8', '--out', root/'scored'))
            summarize(db, root, manifest, model_dir, sample, model_hash)
            files = ['summary.json', 'scores.tsv', 'missing_loci.tsv.gz', 'union.pgen', 'union.pvar', 'union.psam',
                     *[f'chr{c}.audit.tsv.gz' for c in range(1,23)]]
            write_json(root/'score.complete.json', dict(input=expected, files={f: sha(root/f) for f in files}))
        else:
            if not checkpoint(root/'score.complete.json', expected): raise ValueError('Scoring has not completed')
            if checkpoint(root/'ancestry.complete.json', expected): return
            ancestry(cfg, sample, root, model_dir)
            artifacts = ('qc.json', 'ancestry/ancestry_score.txt.gz',
                         'ancestry/ancestry_score.txt.source.json', 'ancestry/format_check.json')
            write_json(root/'ancestry.complete.json', dict(input=expected, files={f: sha(root/f) for f in artifacts}))
    finally: db.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('config', 'gvcf', 'output'): p.add_argument('--'+key, type=Path, required=True)
    p.add_argument('--sample', required=True); p.add_argument('--stage', choices=('score', 'ancestry'), required=True)
    a = p.parse_args()
    run(json.loads(a.config.read_text()), a.sample, a.gvcf.resolve(), a.output.resolve(), a.stage)


if __name__ == '__main__': main()
