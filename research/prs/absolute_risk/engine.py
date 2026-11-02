"""Versioned research calculator; no implicit population data or clinical defaults.

intervene_discrete_v1 reproduces the published interval approximation.
competing_exact_v1 integrates constant cause-specific hazards within each bin.
Neither method has yet been validated against published CHD output tables.
"""
import math

GROUPS = ('<20', '20-40', '40-60', '60-80', '80-90', '90-95', '>=95')
PROPORTIONS = (.2, .2, .2, .2, .1, .05, .05)
METHODS = ('intervene_discrete_v1', 'competing_exact_v1')


def number(value, name, low=0, high=math.inf):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f'{name}: expected numeric value')
    if not math.isfinite(value) or not low <= value <= high:
        raise ValueError(f'{name}: outside [{low}, {high}]')
    return float(value)


def group(percentile):
    p = number(percentile, 'percentile', 0, 100)
    return GROUPS[next((i for i, bound in enumerate((20, 40, 60, 80, 90, 95)) if p < bound), 6)]


def validate(bundle):
    if not isinstance(bundle, dict):
        raise ValueError('Bundle must be a JSON object')
    if bundle.get('schema') != 'CHD_RESEARCH_BUNDLE_V1':
        raise ValueError('Unsupported bundle schema')
    for key in ('model_id', 'baseline_source', 'hr_source', 'baseline_year',
                'baseline_release', 'population', 'limitations', 'score_definition'):
        if not bundle.get(key):
            raise ValueError(f'Missing provenance: {key}')
    if bundle.get('endpoint') != 'first_CHD':
        raise ValueError('Only first_CHD is supported')
    if bundle.get('status') not in ('synthetic', 'unvalidated_external'):
        raise ValueError('This implementation cannot declare clinical validation')
    if bundle.get('method') not in METHODS:
        raise ValueError('Explicit supported numerical method required')
    if bundle.get('units') != 'rates_per_person_year_prevalence_fraction':
        raise ValueError('Normalize units before importing; per-100000 is not accepted')
    definition = bundle['score_definition']
    if not isinstance(definition, dict):
        raise ValueError('score_definition must be an object')
    if definition.get('model') not in ('A', 'B') or not definition.get('reference_id'):
        raise ValueError('Explicit A/B model and reference required')
    digest = definition.get('model_lock_sha256', '')
    if len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
        raise ValueError('Expected model-lock SHA256')
    if definition.get('ancestry') != 'EAS':
        raise ValueError('This adapter currently supports EAS only')
    bins = bundle.get('bins', [])
    if not isinstance(bins, list) or any(not isinstance(r, dict) for r in bins):
        raise ValueError('bins must be a list of objects')
    for row in bins:
        number(row['age_start'], 'age_start', 0, 80)
    for sex in ('male', 'female'):
        rows = sorted((r for r in bins if r.get('sex') == sex), key=lambda r: r['age_start'])
        end = 0
        for r in rows:
            start = number(r['age_start'], 'age_start', 0, 80)
            stop = number(r['age_end'], 'age_end', 0, 80)
            if start != end or stop <= start:
                raise ValueError('Age bins must be contiguous, non-overlapping and positive')
            if bundle['method'] == 'intervene_discrete_v1' and stop - start != 5:
                raise ValueError('Reproduction method requires five-year bins')
            end = stop
            number(r['incidence'], 'incidence')
            number(r['prevalence'], 'prevalence', 0, 1)
            if r['prevalence'] == 1:
                raise ValueError('No disease-free risk set at prevalence=1')
            deaths = number(r['all_cause_mortality'], 'all_cause_mortality')
            number(r['chd_mortality'], 'chd_mortality', 0, deaths)
            hr = r['hr']
            if set(hr) != set(GROUPS) or hr['40-60'] != 1:
                raise ValueError('All seven HRs required; 40-60 reference HR must equal 1')
            for value in hr.values():
                if number(value, 'HR') == 0:
                    raise ValueError('HR must be positive')
        if end != 80:
            raise ValueError(f'{sex}: complete coverage [0,80) required')
    if any(r.get('sex') not in ('male', 'female') for r in bins):
        raise ValueError('Unknown sex in model bins')
    return bundle


def integrate(bundle, sex, percentile, start_age, end_age):
    """Conditional risk given alive and CHD-free at start_age; forward integration."""
    validate(bundle)
    start = number(start_age, 'start_age', 0, 80)
    end = number(end_age, 'end_age', start, 80)
    if sex not in ('male', 'female'):
        raise ValueError('Explicit male/female model sex required')
    selected = group(percentile)
    survival, risk = 1., 0.
    rows = sorted((r for r in bundle['bins'] if r['sex'] == sex), key=lambda r: r['age_start'])
    curve = [{'age': start, 'risk': risk, 'event_free_survival': survival}]
    for row in rows:
        lower, upper = max(start, row['age_start']), min(end, row['age_end'])
        if upper <= lower:
            continue
        dt = upper - lower
        if bundle['method'] == 'intervene_discrete_v1' and dt != 5:
            raise ValueError('Partial bins require competing_exact_v1; do not silently change reproduction')
        weighted_hr = math.fsum(p * row['hr'][g] for p, g in zip(PROPORTIONS, GROUPS))
        hazard = row['incidence'] * row['hr'][selected] / weighted_hr / (1 - row['prevalence'])
        mortality = row['all_cause_mortality'] - row['chd_mortality']
        total = hazard + mortality
        if bundle['method'] == 'intervene_discrete_v1':
            increment = -math.expm1(-hazard * dt)
        else:
            increment = hazard / total * -math.expm1(-total * dt) if total else 0.
        risk += survival * increment
        survival *= math.exp(-total * dt)
        if not 0 <= risk <= 1 + 1e-12:
            raise ValueError('Calculated risk outside probability bounds')
        curve.append({'age': upper, 'risk': risk, 'event_free_survival': survival})
    return {'risk': risk, 'curve': curve, 'group': selected}
