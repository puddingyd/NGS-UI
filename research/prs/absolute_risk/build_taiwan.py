"""Build the adult Taiwan bundle only after reproducing published EstBB CHD curves.

Reads pinned source snapshots; does not download or estimate replacement HRs.
Run from the repository root: python -m research.prs.absolute_risk.build_taiwan --help
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from .engine import GROUPS, integrate, validate

ROOT = Path(__file__).resolve().parents[1]
LABELS = dict(zip(GROUPS, ('< 20%', '20-40%', '40-60%', '60-80%', '80-90%', '90-95%', '> 95%')))
PERCENTILES = dict(zip(GROUPS, (10, 30, 50, 70, 85, 92, 97)))
PINS = {
    'Sex_Stratified_Incidence.csv': '0a5da0a297137547002c90e152a71472da6843d2848eb33628197e34edf1194b',
    'Sex_Stratified_Prevalence.csv': '231e5ba00a052161a28f30e0a63c722c11c473de4427fc90476425624999caeb',
    'Sex_Stratified_Mortality.csv': '286155ea6dae1b6904c7529f5cde5c83fc38a24c9ff38aeab37ba0e8dfeba2f3',
    'HR_AgeandSexStratified_EstBB_211122.csv': 'a94d73b38c229955b2ec67c3e7af6e9bd06c829ee98ab28c85d56e46787ea23c',
    'HR_FemaleSample_EstBB.csv': 'd62ad6f71733f7d95b32e919ca386d3957fd8acea8ec456f54f153ef8d44815b',
    'chd_published_cumulative_risk.csv': 'df04cb4a66c402fb567202819e31759ad3d3d3bc1a7fce000e27b574aea295f3',
    'method2_model_lock.json': 'a3b2bb7a52762d3232cd7e72040b26470c882e5afce66d8e574edc9faeda1083',
}
SOURCE_URLS = {
    'male_hr': 'https://drive.google.com/file/d/18ylwVsjMsaaQ-FASwHat_XuhTyItmI3i/view',
    'female_hr': 'https://drive.google.com/file/d/14IQVWEdxTvp2OloAji2637VDmI8VaCRC/view',
    'paper': 'https://doi.org/10.1038/s41467-024-48938-2',
    'gbd': 'https://github.com/intervene-EU-H2020/flagship/tree/179986db2f9af3b6eb32e4e16dc9201a23d2a29c/AbsoluteRiskEstimation',
}


def pinned(path):
    path = Path(path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != PINS[path.name]:
        raise ValueError(f'Source checksum mismatch: {path.name}')
    return path


def association_rows(path, columns):
    # This exact source has V1... headers. Pin bytes BEFORE the explicit mapping.
    with pinned(path).open(newline='') as stream:
        reader = csv.reader(stream)
        if next(reader) != [''] + [f'V{i}' for i in range(1, len(columns) + 1)]:
            raise ValueError('Unexpected EstBB positional schema')
        result = []
        for row in reader:
            if len(row) != len(columns) + 1:
                raise ValueError('Malformed association row')
            item = dict(zip(columns, row[1:]))
            if item['Phenotype'] == 'I9_CHD':
                if item['PRS'] != 'CHD':
                    raise ValueError('Unexpected score')
                for key in ('HR', 'SE'):
                    if not math.isfinite(float(item[key])) or float(item[key]) <= 0:
                        raise ValueError('Invalid HR/SE')
                result.append(item)
    return result


def male_hr(rows, group, age):
    selected = [r for r in rows if r['Sex'] == 'male' and r['Group'] == LABELS[group]]
    selected.sort(key=lambda r: float(r['MinAge']))
    if len(selected) != 4:
        raise ValueError('Four age strata required for each male PRS group')
    xs = [float(r['MedianAAO']) for r in selected]
    ys = [math.log(float(r['HR'])) for r in selected]
    weights = [1 / float(r['SE']) ** 2 for r in selected]
    if xs != sorted(set(xs)):
        raise ValueError('Distinct ordered MedianAAO values required')
    mean_x = math.fsum(w*x for w, x in zip(weights, xs)) / math.fsum(weights)
    mean_y = math.fsum(w*y for w, y in zip(weights, ys)) / math.fsum(weights)
    slope = math.fsum(w*(x-mean_x)*(y-mean_y) for w, x, y in zip(weights, xs, ys)) / math.fsum(w*(x-mean_x)**2 for w, x in zip(weights, xs))
    return math.exp(mean_y + slope * (min(max(age, xs[0]), xs[-1]) - mean_x))


def source_bins(gbd_dir, hr_dir, country, minimum):
    male = association_rows(hr_dir / 'HR_AgeandSexStratified_EstBB_211122.csv',
                            'Phenotype PRS Sex MinAge MaxAge MedianAAO Group Controls Cases Beta SE Pval HR Cipos Cineg'.split())
    female = association_rows(hr_dir / 'HR_FemaleSample_EstBB.csv',
                              'Phenotype PRS Group Controls Cases Beta SE Pval HR Cipos Cineg'.split())
    values = {}
    for measure in ('Incidence', 'Prevalence', 'Mortality'):
        with pinned(gbd_dir / f'Sex_Stratified_{measure}.csv').open(newline='') as stream:
            for row in csv.DictReader(stream):
                if row['location'] != country or row['metric'] != 'Rate' or row['cause'] not in ('Ischemic heart disease', 'All causes'):
                    continue
                key = (measure, row['sex'].lower(), row['age'], row['cause'])
                low, value, high = (float(row[k]) for k in ('lower', 'val', 'upper'))
                if key in values or row['year'] != '2019' or not (0 <= low <= value <= high < math.inf):
                    raise ValueError('Invalid or duplicate GBD record')
                values[key] = value / 100000
    result = []
    for sex in ('male', 'female'):
        for start in range(minimum, 80, 5):
            # 1–4 → [0,5) is ONLY the author reproduction convention, never
            # used by the deployed adult [20,80) parameter bundle.
            label = f'{start or 1}-{start+4} years'
            hrs = {'40-60': 1.0}
            for group in GROUPS:
                if group == '40-60':
                    continue
                if sex == 'male':
                    hrs[group] = male_hr(male, group, start + 2.5)
                else:
                    matches = [r for r in female if r['Group'] == LABELS[group]]
                    if len(matches) != 1:
                        raise ValueError('Unique female group HR required')
                    hrs[group] = float(matches[0]['HR'])
            result.append(dict(sex=sex, age_start=start, age_end=start+5, hr=hrs,
                               incidence=values['Incidence', sex, label, 'Ischemic heart disease'],
                               prevalence=values['Prevalence', sex, label, 'Ischemic heart disease'],
                               all_cause_mortality=values['Mortality', sex, label, 'All causes'],
                               chd_mortality=values['Mortality', sex, label, 'Ischemic heart disease']))
    return result


def build(sources, hr_dir):
    pinned(ROOT / 'scoring_v211/method2_model_lock.json')
    common = dict(schema='CHD_RESEARCH_BUNDLE_V1', endpoint='first_CHD', status='unvalidated_external',
                  model_id='intervene_estbb_tw2019_A_exact_v1', method='competing_exact_v1',
                  age_min=20, population='Taiwan', baseline_year=2019, baseline_release='GBD2019',
                  baseline_source=SOURCE_URLS['gbd'], hr_source=SOURCE_URLS,
                  hr_cohort='INTERVENE EstBB European ancestry; male age/sex, female sex',
                  units='rates_per_person_year_prevalence_fraction', source_sha256=PINS,
                  score_definition=dict(model='A', reference_id='canonical_v0_EAS735', ancestry='EAS',
                                        model_lock_sha256=PINS['method2_model_lock.json'], pgs_id='PGS004879'),
                  limitations=[
                      '研究用途；未完成台灣 outcome validation，不作治療或篩檢決策。',
                      '歐洲 EstBB 分組 HR 外推至台灣；Model A/EAS735 分位組與原研究分位組相容性未經 outcome validation。',
                      'GBD2019 IHD 與 INTERVENE CHD 定義有差異；固定2019年率不預測未來世代變化。',
                      '以評估時仍存活且無 CHD 為條件；20至79歲，計算終點80歲。',
                      'UI 採五年區間內固定 hazard 的連續競爭風險積分；原文五年離散近似另行驗算。',
                      '尚未估計個人風險的不確定性區間；未納入血壓、血脂、吸菸或治療效果。'])
    reproduction = {**common, 'age_min': 0, 'population': 'Estonia', 'method': 'intervene_discrete_v1',
                    'bins': source_bins(sources/'gbd2019', hr_dir, 'Estonia', 0)}
    expected_path = pinned(sources/'extracted/chd_published_cumulative_risk.csv')
    comparisons, seen = [], set()
    with expected_path.open(newline='') as stream:
        for row in csv.DictReader(stream):
            if row['col_2: Country'] != 'Estonia':
                continue
            sex = {'Males': 'male', 'Females': 'female'}[row['col_3: Sex']]
            group = next(g for g in GROUPS if LABELS[g] == row['col_6: PGS Strata'])
            age = float(row['col_7: Age'])
            key = sex, group, age
            if key in seen or row['col_5: Optimal Model'] != ('Age and Sex' if sex == 'male' else 'Sex'):
                raise ValueError('Duplicate or mismatched published result')
            seen.add(key)
            actual = integrate(reproduction, sex, PERCENTILES[group], 0, age)['risk'] * 100
            expected = float(row['col_8: Cumulative Absolute Risk  (%)'])
            comparisons.append(dict(sex=sex, group=group, age=age, actual_percent=actual,
                                    expected_percent=expected, error_pp=abs(actual-expected)))
    required = {(s,g,float(a)) for s in ('male','female') for g in GROUPS for a in range(5,81,5)}
    maximum = max(r['error_pp'] for r in comparisons)
    if seen != required or maximum > .0050001:
        raise ValueError(f'Published reproduction failed: {len(seen)} rows; max error {maximum} pp')
    evidence = dict(status='PASS', cohort='EstBB', rows=224, max_error_pp=maximum,
                    tolerance_pp=.0050001, expected_source='Supplementary Data 4, Estonia CHD',
                    algorithm='intervene_discrete_v1', uncertainty_reproduced=False)
    bundle = validate({**common, 'reproduction': evidence,
                       'bins': source_bins(sources/'gbd2019', hr_dir, 'Taiwan (Province of China)', 20)})
    deltas = []
    discrete = {**bundle, 'method': 'intervene_discrete_v1'}
    for sex in ('male', 'female'):
        for percentile in PERCENTILES.values():
            for age in range(20, 80, 5):
                for stop in sorted({80, min(age+10,80)}):
                    deltas.append(abs(integrate(bundle, sex, percentile, age, stop)['risk']-
                                      integrate(discrete, sex, percentile, age, stop)['risk'])*100)
    report = dict(**evidence, comparisons=comparisons,
                  taiwan_exact_vs_discrete_max_difference_pp=max(deltas),
                  source_sha256=PINS, source_urls=SOURCE_URLS,
                  clinical_validation='NOT_VALIDATED', score_transfer_validation='NOT_VALIDATED')
    return bundle, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--hr-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    bundle, report = build(args.sources, args.hr_dir)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    for name, value in [('bundle.json', bundle), ('reproduction.json', report)]:
        (args.output_dir/name).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ('status','rows','max_error_pp','taiwan_exact_vs_discrete_max_difference_pp')}))


if __name__ == '__main__':
    main()
