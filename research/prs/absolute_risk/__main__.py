"""CLI: validate bundles or evaluate explicitly supplied research scenarios."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

from .engine import group, integrate, number, validate


def read(path):
    return json.loads(Path(path).read_text())


def evaluate(bundle, summary, person, allow_synthetic=False):
    validate(bundle)
    if not isinstance(summary, dict) or not isinstance(person, dict):
        raise ValueError('Summary and person must be JSON objects')
    if bundle['status'] == 'synthetic' and not allow_synthetic:
        raise ValueError('Synthetic fixture requires --allow-synthetic; never use for real cases')
    if bundle['status'] == 'synthetic' and person.get('synthetic') is not True:
        raise ValueError('Synthetic model requires an explicitly synthetic person')
    if person.get('prior_chd') is not False:
        raise ValueError('First-event risk requires confirmed absence of prior CHD')
    datetime.date.fromisoformat(person['assessment_date'])
    age = number(person['age'], 'age', 0, 80)
    if age == 80:
        raise ValueError('No supported future horizon at age 80')
    if summary.get('sample') != person.get('sample') or not person.get('sample'):
        raise ValueError('Sample identity mismatch')
    if summary.get('schema') != 'PGS004879_METHOD2_V0_20261031' or summary.get('stage') != 'METHOD2_TECHNICAL_PASS':
        raise ValueError('Verified Method 2 output required')
    definition = bundle['score_definition']
    if summary.get('model_lock_sha256') != definition['model_lock_sha256']:
        raise ValueError('Score model lock mismatch')
    ancestry = summary.get('ancestry_status', {})
    if not isinstance(ancestry, dict):
        raise ValueError('Invalid ancestry evidence')
    if ancestry.get('ancestry') != 'EAS' or ancestry.get('ancestry_low_confidence') not in (False, 'False'):
        raise ValueError('Confident, linked EAS ancestry evidence required')
    scores = summary.get('scores', [])
    if not isinstance(scores, list) or any(not isinstance(r, dict) for r in scores):
        raise ValueError('Invalid score records')
    rows = [r for r in scores if r.get('model') == definition['model']]
    if len(rows) != 1:
        raise ValueError('Missing or duplicate selected PRS model')
    row = rows[0]
    p = number(row['eas_percentile'], 'eas_percentile', 0, 100)
    if group(p) != row['intervene_group']:
        raise ValueError('Inconsistent percentile group')
    for field in ('total_sum', 'eas_z'):
        number(row[field], field, -float('inf'))
    expected = 610305 if definition['model'] == 'A' else 603438
    if row['n_loci'] != expected or row['observed_loci'] + row['missing_loci'] != expected:
        raise ValueError('Inconsistent locus counts')
    for field in ('observed_loci', 'missing_loci'):
        value = row[field]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError('Invalid locus count')
    ten = integrate(bundle, person['sex'], p, age, age + 10) if age <= 70 else None
    remaining = integrate(bundle, person['sex'], p, age, 80)
    return {'schema': 'CHD_RESEARCH_RESULT_V1', 'research_only': True,
            'clinical_validation': 'NOT_VALIDATED', 'parameter_status': bundle['status'],
            'model_id': bundle['model_id'], 'method': bundle['method'],
            'population': bundle['population'], 'baseline_year': bundle['baseline_year'],
            'score_definition': definition, 'person': person, 'percentile': p,
            'risk_10_year': ten, 'risk_10_year_unavailable_reason': None if ten else 'Horizon exceeds age 80',
            'risk_to_80': remaining, 'uncertainty_interval': None,
            'uncertainty_status': 'NOT_ESTIMATED',
            'limitations': bundle['limitations'],
            'score_qc': {k: summary.get(k) for k in ('high_weight_missing_modelA', 'strand_ambiguity_warnings')},
            'boundary_review_2pp': min(abs(p-b) for b in (20,40,60,80,90,95)) <= 2}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('validate', 'evaluate'))
    parser.add_argument('--bundle', required=True)
    parser.add_argument('--summary')
    parser.add_argument('--person')
    parser.add_argument('--output')
    parser.add_argument('--allow-synthetic', action='store_true')
    args = parser.parse_args()
    try:
        bundle = validate(read(args.bundle))
        if args.command == 'validate':
            result = {'structure': 'PASS', 'model_id': bundle['model_id'],
                      'parameter_status': bundle['status'], 'scientific_validation': 'NOT_VERIFIED'}
        else:
            if not args.summary or not args.person:
                raise ValueError('--summary and --person are required')
            result = evaluate(bundle, read(args.summary), read(args.person), args.allow_synthetic)
        inputs = {'bundle': args.bundle}
        if args.command == 'evaluate':
            inputs.update(summary=args.summary, person=args.person)
        result['input_sha256'] = {k: hashlib.sha256(Path(v).read_bytes()).hexdigest() for k,v in inputs.items()}
        rendered = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n'
        if args.output:
            # Exclusive creation prevents replacing previously reviewed results.
            with open(args.output, 'x') as stream:
                stream.write(rendered)
        else:
            print(rendered, end='')
    except (ValueError, KeyError, TypeError, OSError) as exc:
        parser.exit(2, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
