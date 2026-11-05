"""Versioned, case-local adapter for old and new canonical ancestry inputs.

The installed model stays immutable: only metadata is completed, and variant
rows are copied unchanged. pgscatalog-format checks the result before Nextflow
can start expensive reference extraction.
"""
from decimal import Decimal
import gzip
from pathlib import Path
import tempfile

from .common import command, sha, write_json

COLUMNS = 'chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n'


def header(n):
    return ('###PGS CATALOG SCORING FILE - see https://www.pgscatalog.org/downloads/#scoring_header\n'
            '#format_version=2.0\n#pgs_id=PGS004879\n'
            '#pgs_name=NGS_UI_fixed_CHD_ancestry_subset\n'
            '#trait_reported=Coronary heart disease\n'
            '#trait_mapped=coronary artery disorder\n#trait_efo=MONDO_0005010\n'
            '#genome_build=GRCh38\n'
            f'#variants_number={n}\n#weight_type=beta\n'
            '#pgp_id=PGP000618\n#citation=Jermy B et al. Nat Commun (2024); '
            'NGS-UI fixed GRCh38 subset for ancestry QC only\n')


def stage_score(source, destination, expected_n):
    """Accept the checksum-verified legacy header; validate and preserve its rows."""
    temporary = destination.with_suffix('.tmp')
    try:
        with gzip.open(source, 'rt') as src, gzip.open(temporary, 'wt') as out:
            metadata = {}
            for line in src:
                if not line.startswith('#'): break
                if line.startswith('#') and not line.startswith('##') and '=' in line:
                    key, value = line[1:].rstrip().split('=', 1)
                    if key in metadata: raise ValueError('Duplicate ancestry score metadata: '+key)
                    metadata[key] = value
            else: raise ValueError('Empty ancestry score')
            if (metadata.get('pgs_id') != 'PGS004879' or metadata.get('genome_build') != 'GRCh38'
                    or metadata.get('variants_number') != str(expected_n) or line != COLUMNS):
                raise ValueError('Canonical ancestry score identity/build/columns mismatch')
            out.write(header(expected_n)); out.write(COLUMNS)
            count = 0
            for line in src:
                fields = line.rstrip('\n').split('\t')
                if len(fields) != 5: raise ValueError('Invalid ancestry score row')
                chrom, pos, ea, oa, weight = fields
                if (not 1 <= int(chrom) <= 22 or int(pos) <= 0 or len(ea) != 1 or len(oa) != 1
                        or ea == oa or not {ea, oa} <= set('ACGT') or not Decimal(weight).is_finite()):
                    raise ValueError('Invalid ancestry score allele/weight')
                out.write(line); count += 1
            if count != expected_n: raise ValueError('Canonical ancestry score variant count mismatch')
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    write_json(destination.with_suffix('.source.json'), dict(schema='PRS_ANCESTRY_INPUT_V1',
        source=str(source), source_sha256=sha(source), staged_sha256=sha(destination),
        variants=count, variant_rows='unchanged', adapter='complete_pgs_header_v1'))
    return destination


def preflight(cfg, score):
    image = Path(cfg.get('pgscatalog_utils_img',
        str(Path(cfg['plink_img']).parent/'ghcr.io-pgscatalog-pygscatalog-pgscatalog-utils-2.0.4.img')))
    if not image.is_absolute() or not image.is_file():
        raise ValueError('找不到祖源格式檢查容器；請設定 pgscatalog_utils_img 絕對路徑：'+str(image))
    print('[PRS] 祖源輸入格式檢查（pgscatalog-format）', flush=True)
    with tempfile.TemporaryDirectory(prefix='.format-check-', dir=score.parent) as work:
        command(['apptainer', 'exec', '--bind', '/home:/home', image, 'pgscatalog-format',
            '-s', score, '-t', 'GRCh38', '--threads', str(cfg.get('raw_threads', 2)),
            '-o', work, '-l', 'log_scorefiles.json'])
    write_json(score.parent/'format_check.json', dict(score_sha256=sha(score),
        formatter='pgscatalog-format', container=str(image), status='PASS'))
    print('[PRS] 祖源輸入格式檢查通過', flush=True)
