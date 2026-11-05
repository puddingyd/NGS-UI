"""Strict six-disease score and clinical-input adapter; no client scores."""
import datetime
import math

from .engine import group, integrate, number, validate
from .registry import DISEASES

SCHEMA = 'PRS_MULTI_SCORE_V1'
REFERENCE = 'unrelated3330_EAS735_v1'


def score_rows(summary, require_ancestry=False):
    if summary.get('schema') != SCHEMA or not summary.get('sample'):
        raise ValueError('Invalid multi-disease score schema')
    if summary.get('reference_id') != REFERENCE:
        raise ValueError('Unrecognized reference population')
    digest = summary.get('model_manifest_sha256', '')
    if len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
        raise ValueError('Missing model release checksum')
    rows = summary.get('diseases')
    if not isinstance(rows, dict) or set(rows) != set(DISEASES):
        raise ValueError('All six disease results required')
    for key, r in rows.items():
        if r.get('pgs_id') != DISEASES[key]['pgs_id']:
            raise ValueError('Disease/PGS mismatch')
        for field in ('total_sum', 'observed_sum', 'imputed_sum', 'eas_z'):
            number(r.get(field), field, -math.inf)
        p = number(r.get('eas_percentile'), 'percentile', 0, 100)
        if group(p) != r.get('intervene_group'):
            raise ValueError('Inconsistent risk stratum')
        for field in ('n_loci', 'observed_loci', 'missing_loci'):
            v = r.get(field)
            if type(v) is not int or v < 0:
                raise ValueError('Invalid locus counts')
        if r['n_loci'] <= 0 or r['observed_loci'] + r['missing_loci'] != r['n_loci']:
            raise ValueError('Inconsistent locus counts')
        if key == 'CHD' and r['n_loci'] != 610305:
            raise ValueError('Frozen CHD locus set changed')
        if not math.isclose(r['total_sum'], math.fsum([r['observed_sum'], r['imputed_sum']]), abs_tol=1e-12):
            raise ValueError('Score reconciliation failed')
        if r.get('qc_ready') is not (r['observed_loci'] / r['n_loci'] >= .995):
            raise ValueError('Incorrect per-disease quality gate')
    if require_ancestry and not ancestry_ready(summary):
        raise ValueError('Confident EAS ancestry evidence required')
    return rows


def ancestry_ready(summary):
    a = summary.get('ancestry_status', {})
    return (a.get('sample') == summary.get('sample') and a.get('ancestry') == 'EAS'
            and a.get('ancestry_low_confidence') in (False, 'False')
            and a.get('qc_status') in ('PASS', 'PASS_WITH_REVIEW_FLAGS')
            and type(a.get('pca_ld_pruned_markers')) is int and a['pca_ld_pruned_markers'] > 0)


def person_for(summary, payload):
    if not isinstance(payload, dict) or set(payload) != {'age', 'sex', 'assessment_date', 'history'}:
        raise ValueError('請提供年齡、性別、評估日期與各疾病病史')
    number(payload['age'], '評估年齡', 20, 79.999999)
    if payload['sex'] not in ('male', 'female'):
        raise ValueError('請選擇模型使用的性別')
    if datetime.date.fromisoformat(payload['assessment_date']) > datetime.date.today():
        raise ValueError('評估日期不能晚於今天')
    history = payload['history']
    if not isinstance(history, dict) or set(history) != set(DISEASES):
        raise ValueError('各疾病病史必須明確填入是、否或尚未確認')
    if any(v is not None and type(v) is not bool for v in history.values()):
        raise ValueError('病史選項格式不正確')
    if history['BREAST'] is True and history['ALL_CANCERS'] is False:
        raise ValueError('乳癌屬於癌症；乳癌病史與整體癌症病史不一致')
    return {**payload, 'sample': summary['sample']}


def validate_pack(pack):
    if pack.get('schema') != 'PRS_RISK_PACK_V2' or set(pack.get('models', {})) != set(DISEASES):
        raise ValueError('Complete six-disease parameter pack required')
    for key, bundle in pack['models'].items():
        validate(bundle)
        p = bundle.get('reproduction', {})
        if (bundle.get('disease_id') != key or bundle.get('status') != 'unvalidated_external'
                or bundle.get('population') != 'Taiwan' or bundle.get('age_min') != 20
                or bundle.get('method') != 'competing_exact_v1'
                or p.get('status') != 'PASS' or p.get('rows') != (112 if key == 'BREAST' else 224)
                or number(p.get('max_error_pp'), 'reproduction error') > .0050001):
            raise ValueError('Unverified parameter pack')
    return pack


def evaluate(pack, summary, person):
    validate_pack(pack)
    rows = score_rows(summary)
    expected_person = person_for(summary, {k: v for k, v in person.items() if k != 'sample'})
    if person != expected_person:
        raise ValueError('Sample identity mismatch')
    results = {}
    age, sex = person['age'], person['sex']
    for key, disease in DISEASES.items():
        r = rows[key]
        reason = None
        if sex not in disease['sexes']:
            reason = '此風險模型僅適用女性'
        elif person['history'][key] is True:
            reason = '已有此疾病；不估算復發風險'
        elif person['history'][key] is None:
            reason = '請確認此疾病病史後計算'
        elif not r['qc_ready']:
            reason = '此疾病計分品質未通過，無法估算風險'
        elif not ancestry_ready(summary):
            reason = '祖源檢查尚未通過，無法估算風險'
        base = dict(id=key, name=disease['name'], pgs_id=disease['pgs_id'],
                    available=reason is None, reason=reason, percentile=r['eas_percentile'])
        if reason:
            results[key] = base
            continue
        b = pack['models'][key]
        p = r['eas_percentile']
        ten = integrate(b, sex, p, age, age+10) if age <= 70 else None
        twenty = integrate(b, sex, p, age, age+20) if age <= 60 else None
        remaining = integrate(b, sex, p, age, 80)
        results[key] = dict(**base, person=person, risk_10_year=ten, risk_20_year=twenty, risk_to_80=remaining,
            reference=dict(group='40-60',
                risk_10_year=integrate(b, sex, 50, age, age+10) if age <= 70 else None,
                risk_20_year=integrate(b, sex, 50, age, age+20) if age <= 60 else None,
                risk_to_80=integrate(b, sex, 50, age, 80)),
            model_id=b['model_id'], limitations=b['limitations'],
            conditioning='alive_and_free_of_this_disease_at_assessment', uncertainty_interval=None)
    return dict(schema='PRS_MULTI_RESULT_V1', person=person, diseases=results,
                research_only=True, clinical_validation='NOT_VALIDATED')
