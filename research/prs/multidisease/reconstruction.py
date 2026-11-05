"""Reuse the frozen, tested gVCF observation rules; scan each chromosome once.

Only function/import definitions are loaded from the pinned vendor heredoc.
No vendor main block, shell code or environment initialization is executed.
"""
import ast
from collections import defaultdict
from pathlib import Path
import subprocess
import types

from .common import SCORING, sha


def rules(root, cfg):
    vendor = SCORING/'vendor/PGS004879_FINAL_gVCF_resume.sh'
    pins = dict(line.split(maxsplit=1) for line in (SCORING/'BASELINE_V2_SHA256SUMS.txt').read_text().splitlines() if line.strip())
    if not any(digest == sha(vendor) and name.strip() == 'vendor/'+vendor.name for digest, name in pins.items()):
        raise ValueError('Frozen genotype reconstruction changed')
    source = vendor.read_text().split("<<'PY_GVCF'\n", 1)[1].split('\nPY_GVCF\n', 1)[0]
    tree = ast.parse(source)
    tree.body = [node for node in tree.body if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef))]
    module = types.ModuleType('frozen_gvcf_rules')
    exec(compile(tree, str(vendor), 'exec'), module.__dict__)
    module.ROOT = root
    module.SRC = root/'discovery'; module.AUD = root/'audit'
    for p in (module.SRC, module.AUD, root/'logs'): p.mkdir(parents=True, exist_ok=True)
    module.REF = Path(cfg['fasta']); module.DB = Path(cfg['dbsnp'])
    module.VALID = set('ACGT'); module.COMP = str.maketrans('ACGT', 'TGCA')
    module.MIN_DP = 10; module.MIN_GQ = 20; module.PALINDROMIC_MODE = 'assume-direct'
    return module


def scan_records(records, targets, bases, rule):
    """Yield every target, including missing sites; preserve same-position ties."""
    by_pos = defaultdict(list)
    for target in targets: by_pos[target['pos']].append(target)
    ordered = sorted(by_pos); index = 0; active = []; last = 0

    def at(pos):
        for target in by_pos[pos]:
            key = (str(target['chrom']), pos)
            score = dict(effect_allele=target['alt'], other_allele=target['ref'], effect_weight='0')
            observation = rule.choose_observation(key, active, bases, score)
            # The preparation manifest already proved forward orientation.
            # A target reference discrepancy must never trigger a new flip.
            if bases[key] != target['ref']:
                observation['STATUS'] = 'REFERENCE_FASTA_MISMATCH'
            accepted = observation['STATUS'] == 'CANDIDATE_DIRECT'
            dosage = int(observation['ALT_DOSAGE']) if accepted else None
            yield target, dosage, observation['STATUS']

    for source in records:
        if source['pos'] < last: raise ValueError('Unsorted gVCF')
        last = source['pos']
        while index < len(ordered) and ordered[index] < source['pos']:
            yield from at(ordered[index]); index += 1
        active = [r for r in active if r['end'] >= source['pos']]
        active.append(source)
    while index < len(ordered):
        yield from at(ordered[index]); index += 1


def scan_chromosome(gvcf, chromosome, targets, rule):
    rule.GVCF = Path(gvcf)
    contig = rule.get_contigs()[str(chromosome)]
    bases = rule.fasta_at_targets({(str(t['chrom']), t['pos']): {} for t in targets})
    with (rule.ROOT/f'logs/gvcf_chr{chromosome}.log').open('w') as log:
        proc = subprocess.Popen(['bcftools', 'view', '-H', '-r', contig, str(gvcf)],
                                stdout=subprocess.PIPE, stderr=log, text=True)
        try:
            def records():
                for line in proc.stdout:
                    r = rule.source_record(line.rstrip('\n').split('\t'))
                    if r['c'] != contig: raise ValueError('Unexpected chromosome from bcftools')
                    yield r
            yield from scan_records(records(), targets, bases, rule)
            if proc.wait(): raise RuntimeError('gVCF reading failed; see '+str(log.name))
        finally:
            if proc.poll() is None: proc.terminate(); proc.wait()
