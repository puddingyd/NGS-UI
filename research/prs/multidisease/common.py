import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SCORING = ROOT / 'scoring_v211'
MODEL_FILES = frozenset(('models.sqlite', 'weights.tsv', 'reference.afreq', 'reference_scores.json',
    'ancestry.ids', 'ancestry_score.txt.gz', 'reference_calibration_3330.keep', 'reference_calibration_eas_735.keep'))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(2**20), b''): h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name('.'+path.name+'.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    os.replace(temporary, path)


def command(args, **kwargs):
    print('[PRS] '+str(args[0])+': '+str(args[-1]), flush=True)
    subprocess.run(list(map(str, args)), check=True, **kwargs)


def plink(cfg, *args):
    return ['apptainer', 'exec', '--bind', '/home:/home', cfg['plink_img'], 'plink2',
            '--threads', str(cfg.get('raw_threads', 2)), '--memory', str(cfg.get('raw_memory_mb', 8192)), *map(str, args)]


def rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        yield from csv.DictReader((line for line in f if not line.startswith('##')), delimiter='\t')


def score_values(path, disease_ids):
    result = {}
    for row in rows(path):
        iid = row.get('#IID') or row.get('IID')
        if not iid or iid in result: raise ValueError('Missing/duplicate reference identity')
        values = {d: float(row[d+'_SUM']) for d in disease_ids}
        if not all(math.isfinite(x) for x in values.values()): raise ValueError('Invalid reference score')
        result[iid] = values
    return result


def load_models(directory, verify=True):
    directory = Path(directory)
    manifest = json.loads((directory/'manifest.json').read_text())
    from ..absolute_risk.registry import DISEASES
    if (manifest.get('schema') != 'PRS_MULTI_MODELS_V1' or set(manifest.get('models', {})) != set(DISEASES)
            or manifest.get('reference_id') != 'unrelated3330_EAS735_v1'
            or manifest.get('reference_n') != 3330 or manifest.get('eas_n') != 735
            or set(manifest.get('sha256', {})) != MODEL_FILES
            or set(manifest.get('offsets', {})) != set(DISEASES)
            or manifest.get('chd_reference_regression') != 'PASS'):
        raise ValueError('Incomplete or unverified six-disease calibration')
    for key, model in manifest['models'].items():
        n, total = model.get('n_loci'), model.get('source_loci')
        if (model.get('pgs_id') != DISEASES[key]['pgs_id'] or type(n) is not int or type(total) is not int
                or not 0 < n <= total or n/total < .995 or (key == 'CHD' and n != 610305)
                or not math.isfinite(float(manifest['offsets'][key]))):
            raise ValueError('Invalid fixed model identity or coverage')
    for name, digest in manifest['sha256'].items():
        path = directory/name
        if Path(name).is_absolute() or '..' in Path(name).parts or path.is_symlink():
            raise ValueError('Unsafe model artifact')
        if verify and sha(path) != digest: raise ValueError('Model artifact checksum mismatch: '+name)
    return manifest, sha(directory/'manifest.json')
