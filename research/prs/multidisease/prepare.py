"""Build all six fixed models from the full GRCh38 HGDP+1kGP PGEN.

One-time operator command, never run on case loading. Existing CHD resources
are read-only. A manifest is written only after CHD reference regression passes.
"""
import argparse
import csv
from decimal import Decimal
import gzip
import importlib.util
import json
from pathlib import Path
import shutil
import sqlite3
import subprocess

from ..absolute_risk.registry import DISEASES, NEW_DISEASES
from .common import SCORING, command, plink, rows, score_values, sha, write_json
from .reconstruction import rules


def legacy():
    spec = importlib.util.spec_from_file_location('fixed_chd_model', SCORING/'method2_score.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def score_rows(path, pgs_id):
    with gzip.open(path, 'rt') as f:
        metadata = []
        for line in f:
            if not line.startswith('#'): break
            metadata.append(line.strip())
        if '#pgs_id='+pgs_id not in metadata or '#HmPOS_build=GRCh38' not in metadata:
            raise ValueError('Score identity/build mismatch')
        reader = csv.DictReader(f, fieldnames=line.rstrip().split('\t'), delimiter='\t')
        for n, r in enumerate(reader, 1):
            chrom, pos = int(r['hm_chr']), int(r['hm_pos'])
            ea, oa = r['effect_allele'].upper(), r['other_allele'].upper()
            weight = Decimal(r['effect_weight'])
            if (not 1 <= chrom <= 22 or pos <= 0 or len(ea) != 1 or len(oa) != 1
                    or len({ea, oa}) != 2 or not {ea, oa} <= set('ACGT') or not weight.is_finite()
                    or r.get('is_interaction') not in ('False', '') or r.get('is_recessive') not in ('False', '')):
                raise ValueError(f'Unsupported score row {pgs_id}:{n}')
            yield n, chrom, pos, ea, oa, str(weight), r.get('rsID', ''), r.get('hm_rsID', '')


def pvar_rows(prefix):
    plain, compressed = Path(str(prefix)+'.pvar'), Path(str(prefix)+'.pvar.zst')
    if plain.is_file():
        yield from rows(plain)
    elif compressed.is_file():
        p = subprocess.Popen(['zstd', '-dc', str(compressed)], stdout=subprocess.PIPE, text=True)
        try:
            yield from csv.DictReader((s for s in p.stdout if not s.startswith('##')), delimiter='\t')
            if p.wait(): raise RuntimeError('Reference PVAR decompression failed')
        finally:
            if p.poll() is None: p.terminate(); p.wait()
    else: raise FileNotFoundError('Reference PVAR not found')


def candidates(source, reference):
    """Direct forward alleles take precedence; complements need dbSNP proof."""
    ea, oa = source['ea'], source['oa']
    ref, alt = reference['REF'], reference['ALT']
    if {ea, oa} == {ref, alt}: return ea, 0
    comp = str.maketrans('ACGT', 'TGCA')
    if {ea.translate(comp), oa.translate(comp)} == {ref, alt}: return ea.translate(comp), 1
    return None


def prepare(cfg, score_dir, reference, output):
    if output.exists(): raise ValueError('Output already exists; use a new versioned model directory')
    m2 = legacy()
    locked = Path(cfg['method2_model_dir'])
    m2.validate_locked_model(locked, cfg['method2_lock'], cfg['original_scorefile'])
    output.mkdir(parents=True)
    rule = rules(output, cfg)
    db = sqlite3.connect(output/'models.sqlite'); db.row_factory = sqlite3.Row
    db.executescript('''
        CREATE TABLE source(disease TEXT,row_no INTEGER,chrom INTEGER,pos INTEGER,ea TEXT,oa TEXT,
          weight TEXT,rsid TEXT,hm_rsid TEXT, PRIMARY KEY(disease,row_no), UNIQUE(disease,chrom,pos,ea,oa));
        CREATE INDEX source_position ON source(chrom,pos);
        CREATE TABLE candidate(disease TEXT,row_no INTEGER,id TEXT,chrom INTEGER,pos INTEGER,
          ref TEXT,alt TEXT,ea TEXT,rank INTEGER, UNIQUE(disease,row_no,id));
        CREATE TABLE locus(id TEXT PRIMARY KEY,chrom INTEGER,pos INTEGER,ref TEXT,alt TEXT,freq TEXT,
          UNIQUE(chrom,pos,ref,alt));
        CREATE TABLE model(disease TEXT,id TEXT,ea TEXT,weight TEXT,freq TEXT,row_no INTEGER,
          PRIMARY KEY(disease,id), UNIQUE(disease,row_no));
        CREATE INDEX model_locus ON model(id);
    ''')
    pins = {r['pgs_id']: r for r in json.loads((Path(__file__).parent/'sources.json').read_text())}
    for key in NEW_DISEASES:
        pid = DISEASES[key]['pgs_id']; path = score_dir/f'{pid}_hmPOS_GRCh38.txt.gz'
        if sha(path) != pins[pid]['sha256']: raise ValueError('Official score checksum mismatch: '+pid)
        db.executemany('INSERT INTO source VALUES(?,?,?,?,?,?,?,?,?)', ((key, *r) for r in score_rows(path, pid)))
    master = m2.load_master(locked/'reference_qc_master_v0.tsv.gz')
    # CHD uses exactly the old canonical allele, weight and row identities.
    db.executemany('INSERT INTO source VALUES(?,?,?,?,?,?,?,?,?)',
        (('CHD', r['row'], int(c), p, r['ea'], r['oa'], str(r['weight']), r['id'], r['id'])
         for (c, p), r in master.items()))
    db.commit()
    coordinates = {(r[0], r[1]) for r in db.execute('SELECT DISTINCT chrom,pos FROM source')}
    print(f'[PRS] 共用參考匹配：{len(coordinates):,} 個座標', flush=True)
    for r in pvar_rows(reference):
        c = r['#CHROM'].removeprefix('chr')
        if not c.isdigit() or (int(c), int(r['POS'])) not in coordinates: continue
        if len(r['REF']) != 1 or len(r['ALT']) != 1 or not {r['REF'], r['ALT']} <= set('ACGT'): continue
        if not r['ID'] or r['ID'] == '.': continue
        for s in db.execute('SELECT * FROM source WHERE chrom=? AND pos=?', (int(c), int(r['POS']))).fetchall():
            match = candidates(s, r)
            if match:
                db.execute('INSERT INTO candidate VALUES(?,?,?,?,?,?,?,?,?)',
                    (s['disease'], s['row_no'], r['ID'], int(c), int(r['POS']), r['REF'], r['ALT'], *match))
    db.commit(); del coordinates
    # Refuse equally plausible reference alleles rather than choose a row order.
    selection = db.execute('''SELECT c.*,s.oa,s.ea AS source_ea,s.weight,s.rsid,s.hm_rsid FROM candidate c
        JOIN source s USING(disease,row_no)
        WHERE c.rank=(SELECT MIN(b.rank) FROM candidate b WHERE b.disease=c.disease AND b.row_no=c.row_no)
        ORDER BY c.chrom,c.pos''').fetchall()
    counts = {}
    for r in selection: counts[r['disease'], r['row_no']] = counts.get((r['disease'], r['row_no']), 0)+1
    chosen = [r for r in selection if counts[r['disease'], r['row_no']] == 1]
    flips = {i: dict(CHROM=str(r['chrom']), POS=r['pos'], ORIENTATION='FLIP', STATUS='CANDIDATE_FLIP',
        FASTA_REF=r['ref'], COMP_EFFECT=r['ea'], COMP_OTHER=r['oa'].translate(rule.COMP),
        RSID=r['rsid'], HM_RSID=r['hm_rsid']) for i, r in enumerate(chosen) if r['rank'] == 1}
    rule.check_flips(flips)
    for i, r in enumerate(chosen):
        if i in flips and flips[i]['STATUS'] != 'FLIP_CONFIRMED': continue
        if r['disease'] == 'CHD':
            old = master[str(r['chrom']), r['pos']]
            if (r['id'], r['ref'], r['alt'], r['ea']) != (old['id'], old['ref'], old['alt'], old['ea']):
                raise ValueError('CHD reference identity changed')
        existing = db.execute('SELECT chrom,pos,ref,alt FROM locus WHERE id=?', (r['id'],)).fetchone()
        if existing and tuple(existing) != (r['chrom'], r['pos'], r['ref'], r['alt']):
            raise ValueError('Duplicate reference ID at different alleles')
        if not existing:
            db.execute('INSERT INTO locus VALUES(?,?,?,?,?,NULL)', (r['id'], r['chrom'], r['pos'], r['ref'], r['alt']))
        db.execute('INSERT INTO model VALUES(?,?,?,?,NULL,?)', (r['disease'], r['id'], r['ea'], r['weight'], r['row_no']))
    models = {}
    for key in DISEASES:
        total = db.execute('SELECT count(*) FROM source WHERE disease=?', (key,)).fetchone()[0]
        matched = db.execute('SELECT count(*) FROM model WHERE disease=?', (key,)).fetchone()[0]
        if matched/total < .995 or (key == 'CHD' and matched != 610305):
            raise ValueError(f'{key}: insufficient fixed reference match {matched}/{total}')
        models[key] = dict(pgs_id=DISEASES[key]['pgs_id'], source_loci=total, n_loci=matched)
    ids = output/'union.ids'
    ids.write_text(''.join(r[0]+'\n' for r in db.execute('SELECT id FROM locus ORDER BY chrom,pos,ref,alt')))
    for name in ('reference_calibration_3330.keep', 'reference_calibration_eas_735.keep'):
        shutil.copyfile(locked/name, output/name)
    prefix_args = ['--pfile', str(reference)]
    if not Path(str(reference)+'.pvar').is_file(): prefix_args += ['vzs']
    command(plink(cfg, *prefix_args, '--extract', ids, '--keep', output/'reference_calibration_3330.keep',
        '--make-pgen', '--out', output/'reference'))
    command(plink(cfg, '--pfile', output/'reference', '--freq', '--out', output/'reference'))
    frequency = {}
    for r in rows(output/'reference.afreq'):
        f = Decimal(r['ALT_FREQS'])
        if not f.is_finite() or not 0 <= f <= 1 or int(r['OBS_CT']) <= 0: raise ValueError('Invalid reference frequency')
        if r['ID'] in frequency: raise ValueError('Duplicate frequency row')
        frequency[r['ID']] = str(f)
        loc = db.execute('SELECT ref,alt FROM locus WHERE id=?', (r['ID'],)).fetchone()
        if not loc or tuple(loc) != (r['REF'], r['ALT']): raise ValueError('Frequency allele mismatch')
        db.execute('UPDATE locus SET freq=? WHERE id=?', (str(f), r['ID']))
    if len(frequency) != db.execute('SELECT count(*) FROM locus').fetchone()[0]: raise ValueError('Missing reference frequencies')
    for r in db.execute('SELECT m.*,l.chrom,l.pos,l.ref,l.alt,l.freq AS af FROM model m JOIN locus l USING(id)').fetchall():
        f = Decimal(r['af']) if r['ea'] == r['alt'] else 1-Decimal(r['af'])
        if r['disease'] == 'CHD':
            fixed = master[str(r['chrom']), r['pos']]['p']
            if abs(f-fixed) > Decimal('0.000001'): raise ValueError('CHD reference frequency changed')
            f = fixed
        db.execute('UPDATE model SET freq=? WHERE disease=? AND id=?', (str(f), r['disease'], r['id']))
    db.commit()
    offsets = {key: Decimal(0) for key in DISEASES}
    with (output/'weights.tsv').open('w') as f:
        f.write('ID\tALT\t'+'\t'.join(DISEASES)+'\n')
        for loc in db.execute('SELECT * FROM locus ORDER BY chrom,pos,ref,alt'):
            weights = {key: '0' for key in DISEASES}
            for r in db.execute('SELECT * FROM model WHERE id=?', (loc['id'],)):
                w = Decimal(r['weight'])
                if r['ea'] == loc['ref']: offsets[r['disease']] += 2*w; w = -w
                weights[r['disease']] = str(w)
            f.write(loc['id']+'\t'+loc['alt']+'\t'+'\t'.join(weights.values())+'\n')
    command(plink(cfg, '--pfile', output/'reference', '--read-freq', output/'reference.afreq',
        '--score', output/'weights.tsv', '1', '2', 'header-read', 'cols=+scoresums',
        '--score-col-nums', '3-8', '--out', output/'calibration'))
    reference_scores = score_values(output/'calibration.sscore', DISEASES)
    members = m2.read_ids(locked/'reference_calibration_3330.keep')
    eas = m2.read_ids(locked/'reference_calibration_eas_735.keep')
    if set(reference_scores) != members or len(members) != 3330 or len(eas) != 735 or not eas <= members:
        raise ValueError('Reference membership changed')
    old_scores = m2.sscore(locked/'reference_modelA_v0.sscore')
    delta = max(abs(reference_scores[i]['CHD']+float(offsets['CHD'])-old_scores[i]) for i in members)
    if delta > 1e-5: raise ValueError(f'CHD calibration regression failed: {delta}')
    distribution = {key: {i: (old_scores[i] if key == 'CHD' else reference_scores[i][key]+float(offsets[key])) for i in sorted(eas)} for key in DISEASES}
    write_json(output/'reference_scores.json', distribution)
    # A small, fixed CHD-only ancestry input prevents other score additions from
    # changing PCA markers. Scoring still uses all six models in one matrix.
    (output/'ancestry.ids').write_text(''.join(r[0]+'\n' for r in db.execute("SELECT id FROM model WHERE disease='CHD' ORDER BY id")))
    with gzip.open(output/'ancestry_score.txt.gz', 'wt') as f:
        f.write('#pgs_id=PGS004879\n#genome_build=GRCh38\n#variants_number=610305\n')
        f.write('chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n')
        for r in db.execute("SELECT l.*,m.ea,m.weight FROM locus l JOIN model m USING(id) WHERE disease='CHD' ORDER BY chrom,pos"):
            f.write(f"{r['chrom']}\t{r['pos']}\t{r['ea']}\t{r['ref'] if r['ea']==r['alt'] else r['alt']}\t{r['weight']}\n")
    # Keep rejected/ambiguous records for a preparation audit, outside case data.
    db.commit(); db.close()
    files = ('models.sqlite', 'weights.tsv', 'reference.afreq', 'reference_scores.json', 'ancestry.ids',
             'ancestry_score.txt.gz', 'reference_calibration_3330.keep', 'reference_calibration_eas_735.keep')
    manifest = dict(schema='PRS_MULTI_MODELS_V1', reference_id='unrelated3330_EAS735_v1',
        reference_n=3330, eas_n=735, models=models, reference_pfile=str(reference),
        chd_reference_regression='PASS', chd_max_score_difference=delta,
        offsets={k: str(v) for k, v in offsets.items()}, source_pins=pins,
        sha256={name: sha(output/name) for name in files})
    write_json(output/'manifest.json', manifest)
    print('[PRS] 六疾病參考模型完成；CHD 3330 人回歸檢查通過。', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('config', 'score-dir', 'reference-pfile', 'output'): p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args()
    prepare(json.loads(a.config.read_text()), a.score_dir.resolve(), a.reference_pfile.resolve(), a.output.resolve())


if __name__ == '__main__': main()
