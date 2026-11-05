"""Sample-bound CHD research estimates from an existing Method2 score.

No genotype execution, client-provided percentile, or clinical report mutation.
"""
import datetime as dt
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import re
import sys
import tempfile
import uuid

from .. import config
from . import sample_layout

# The deployment PYTHONPATH contains backend/, not necessarily the repo root.
# Load the single shared implementation without mutating sys.path or copying it.
_package = '_ngs_prs_risk'
if _package not in sys.modules:
    _directory = config.REPO_ROOT / 'research/prs/absolute_risk'
    _spec = importlib.util.spec_from_file_location(_package, _directory/'__init__.py',
                                                  submodule_search_locations=[str(_directory)])
    _module = importlib.util.module_from_spec(_spec)
    sys.modules[_package] = _module
    _spec.loader.exec_module(_module)
engine = importlib.import_module(_package + '.engine')
calculator = importlib.import_module(_package + '.__main__')

# Released, locally reproduced parameter artifact. A renamed/edited/synthetic
# bundle cannot become active merely by claiming that reproduction passed.
SCORE_DEFINITION = dict(model='A', model_lock_sha256='a3b2bb7a52762d3232cd7e72040b26470c882e5afce66d8e574edc9faeda1083')

RELEASED_BUNDLE_SHA256 = '322d4e68d49a2047d7326172a609b3898b8754f658f800f50806f4e2276a0a27'


class NotReady(ValueError):
    pass


def validate_sid(sid):
    if not isinstance(sid, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}', sid) or '..' in sid:
        raise ValueError('個案識別碼格式不正確')
    return sid


def read_json(path):
    if path.stat().st_size > 2_000_000:
        raise ValueError('PRS 檔案超出允許大小')
    raw = path.read_bytes()
    value = json.loads(raw)
    if not isinstance(value, dict):
        raise ValueError('PRS 檔案必須是 JSON object')
    return value, hashlib.sha256(raw).hexdigest()


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n'
    fd, temporary = tempfile.mkstemp(prefix='.prs-', dir=path.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def score_path(sid):
    validate_sid(sid)
    return sample_layout.state_dir(sid) / f'{sid}.prs.summary.json'


def assessment_dir(sid):
    validate_sid(sid)
    return sample_layout.state_dir(sid) / f'{sid}.prs' / 'assessments'


def case_metadata(sid):
    validate_sid(sid)
    path = sample_layout.state_file(sid, 'sample_metadata.json')
    if not path.is_file():
        raise FileNotFoundError('找不到已載入個案')
    return read_json(path)


def load_bundle():
    if not config.PRS_BUNDLE.is_file():
        raise NotReady('尚未安裝 CHD 研究風險參數包')
    try:
        bundle, digest = read_json(config.PRS_BUNDLE)
        if digest != RELEASED_BUNDLE_SHA256:
            raise ValueError('Parameter release checksum mismatch')
        engine.validate(bundle)
        proof = bundle.get('reproduction', {})
        if (bundle['status'] != 'unvalidated_external'
                or bundle['model_id'] != 'intervene_estbb_tw2019_A_exact_v1'
                or bundle['method'] != 'competing_exact_v1'
                or bundle.get('age_min') != 20 or bundle['population'] != 'Taiwan'
                or bundle['score_definition']['model'] != 'A'
                or proof.get('status') != 'PASS' or proof.get('rows') != 224
                or not 0 <= engine.number(proof.get('max_error_pp'), 'reproduction error') <= .0050001):
            raise ValueError('Unsupported or unverified parameter bundle')
    except (ValueError, TypeError, KeyError) as exc:
        raise NotReady('CHD 參數包未通過版本或重現檢查') from exc
    return bundle, digest


def load_score_record(sid):
    path = score_path(sid)
    if not path.is_file():
        raise NotReady('尚未匯入此個案的 PRS 計分結果')
    record, digest = read_json(path)
    if record.get('schema') != 'PRS_SCORE_LINK_V1' or record.get('sample_id') != sid:
        raise NotReady('PRS 計分結果與個案連結不符')
    summary = record.get('summary')
    if not isinstance(summary, dict) or summary.get('sample') != record.get('source_sample_id'):
        raise NotReady('PRS 來源樣本識別碼不符')
    if record.get('pipeline_source_sha256') != pipeline_source_hash(sid):
        raise NotReady('個案分析來源已更新，請重新執行 PRS 分析')
    if record.get('gvcf_signature'):
        from .prs_jobs import signature
        try:
            if any(signature(stamp[0]) != stamp for stamp in record['gvcf_signature']):
                raise NotReady('原始 gVCF 已變更，請重新分析 PRS')
        except OSError as exc:
            raise NotReady('原始 gVCF 已移動或不存在，請重新確認來源') from exc
    return record, digest


def load_score(sid):
    record, digest = load_score_record(sid)
    return record['summary'], digest


def pipeline_source_hash(sid):
    path = sample_layout.state_file(sid, 'pipeline_source.json')
    return read_json(path)[1] if path.is_file() else None


def source_identity_matches(sid, source):
    if source == sid:
        return True
    path = sample_layout.state_file(sid, 'pipeline_source.json')
    if path.is_file():
        marker, _ = read_json(path)
        return bool(source) and marker.get('source_sample_id') == source
    return False


def import_score(sid, summary_path):
    """Operator CLI import; identity must match LIS or its pipeline source marker."""
    case_metadata(sid)
    summary, digest = read_json(Path(summary_path))
    if not source_identity_matches(sid, summary.get('sample')):
        raise ValueError('summary.sample 不符合個案或 pipeline_source.source_sample_id')
    if summary.get('schema') != 'PGS004879_METHOD2_V0_20261031' or summary.get('stage') != 'METHOD2_TECHNICAL_PASS':
        raise ValueError('只接受完成技術檢查的 Method2 summary.json')
    record = dict(schema='PRS_SCORE_LINK_V1', sample_id=sid, source_sample_id=summary['sample'],
                  original_sha256=digest, pipeline_source_sha256=pipeline_source_hash(sid),
                  imported_at=dt.datetime.now(dt.timezone.utc).isoformat(), summary=summary)
    atomic_json(score_path(sid), record)
    return {'sample_id': sid, 'imported': True, 'original_sha256': digest}


def person_for(summary, payload):
    if set(payload) != {'age', 'sex', 'assessment_date', 'prior_chd'}:
        raise ValueError('請提供評估年齡、性別、日期與 CHD 病史')
    date = dt.date.fromisoformat(payload['assessment_date'])
    if date > dt.date.today():
        raise ValueError('評估日期不能晚於今天')
    if payload['prior_chd'] is not False:
        raise ValueError('僅適用確認尚未罹患 CHD 的個案')
    engine.number(payload['age'], '評估年齡', 20, 79.999999)
    if payload['sex'] not in ('male', 'female'):
        raise ValueError('請選擇模型使用的性別')
    return {**payload, 'sample': summary['sample']}


def latest_snapshot(sid):
    directory = assessment_dir(sid)
    # Filenames are generated UUIDs; no request path is interpolated here.
    candidates = sorted(directory.glob('*.json'), key=lambda p: p.stat().st_mtime_ns, reverse=True) if directory.is_dir() else []
    for path in candidates:
        if not re.fullmatch(r'[0-9a-f]{32}\.json', path.name):
            continue
        value, _ = read_json(path)
        if value.get('sample_id') == sid and value.get('schema') == 'PRS_ASSESSMENT_V1':
            return value
    return None


def status(sid):
    try:
        candidate, candidate_hash = load_score_record(sid)
        if candidate['summary'].get('schema') == 'PRS_MULTI_SCORE_V1':
            from . import prs_multi
            return prs_multi.status(sid, candidate, candidate_hash)
    except (OSError, ValueError, KeyError, TypeError):
        pass
    _, meta_hash = case_metadata(sid)
    issues, bundle, summary, score = [], None, None, None
    bundle_hash = score_hash = None
    try:
        record, score_hash = load_score_record(sid)
        summary = record['summary']
        if not source_identity_matches(sid, summary['sample']):
            raise NotReady('PRS 來源與目前個案來源不符，請重新分析')
        row = calculator.score_row(summary, SCORE_DEFINITION, require_ancestry=False)
        try:
            calculator.score_row(summary, SCORE_DEFINITION)
            qc_ready = True
        except (ValueError, KeyError, TypeError):
            qc_ready = False
        score = dict(raw_score=row['total_sum'], percentile=row['eas_percentile'],
                     group=row['intervene_group'], model='A', job_id=record.get('job_id'),
                     qc_ready=qc_ready,
                     ancestry=summary.get('ancestry_status', {}),
                     observed_loci=row['observed_loci'], missing_loci=row['missing_loci'])
    except (OSError, ValueError, KeyError, TypeError):
        summary = None
        issues.append('尚無有效且與此個案相符的 PRS 結果')
    try:
        bundle, bundle_hash = load_bundle()
    except (OSError, ValueError) as exc:
        issues.append(str(exc))
    if bundle and summary:
        try:
            calculator.score_row(summary, bundle['score_definition'])
        except (ValueError, KeyError, TypeError):
            issues.append('PRS 祖源或品質檢查尚未符合風險模型要求')
    current = dict(bundle=bundle_hash, score=score_hash, case_metadata=meta_hash)
    latest = latest_snapshot(sid)
    stale = bool(latest and latest.get('input_sha256') != current)
    from . import prs_jobs, prs_presentation
    return dict(available=not issues, issues=issues, score=score, research_only=True,
                scope_note='此個案目前僅有 CHD 計分；其他疾病需完成六疾病 PRS 分析。',
                diseases=[dict(id='CHD', name='冠心病（CHD）', score=score)],
                pipeline=prs_jobs.overview(sid), explanation=prs_presentation.explanation(),
                age_min=20, age_max_exclusive=80,
                model_id=bundle['model_id'] if bundle else None,
                latest=latest if not stale and not issues else None, stale=stale)


def assess(sid, payload, username):
    candidate, _ = load_score(sid)
    if candidate.get('schema') == 'PRS_MULTI_SCORE_V1':
        from . import prs_multi
        return prs_multi.assess(sid, payload, username)
    if isinstance(payload, dict) and 'history' in payload:
        from . import prs_multi
        person = prs_multi.calculator.person_for(candidate, payload)
        payload = {k: person[k] for k in ('age', 'sex', 'assessment_date')}
        payload['prior_chd'] = person['history']['CHD']
    _, meta_hash = case_metadata(sid)
    bundle, bundle_hash = load_bundle()
    summary, score_hash = load_score(sid)
    if not source_identity_matches(sid, summary['sample']):
        raise NotReady('PRS 來源與目前個案來源不符')
    person = person_for(summary, payload)
    try:
        result = calculator.evaluate(bundle, summary, person)
    except (ValueError, KeyError, TypeError) as exc:
        raise ValueError('PRS 品質、祖源或計分版本不符合風險模型，無法估算') from exc
    hashes = dict(bundle=bundle_hash, score=score_hash, case_metadata=meta_hash)
    # Refuse to publish when a concurrent import/metadata update changed inputs.
    if (read_json(config.PRS_BUNDLE)[1] != bundle_hash or load_score(sid)[1] != score_hash
            or case_metadata(sid)[1] != meta_hash or not source_identity_matches(sid, summary['sample'])):
        raise NotReady('輸入資料已更新，請重新估算')
    record = dict(schema='PRS_ASSESSMENT_V1', sample_id=sid, assessment_id=uuid.uuid4().hex,
                  created_at=dt.datetime.now(dt.timezone.utc).isoformat(), created_by=username,
                  input_sha256=hashes, inputs=dict(person=person, summary=summary, bundle=bundle), result=result)
    atomic_json(assessment_dir(sid)/f"{record['assessment_id']}.json", record)
    return record


def report_assessment(sid):
    """Require a current saved assessment; never create one during export."""
    candidate, _ = load_score(sid)
    if candidate.get('schema') == 'PRS_MULTI_SCORE_V1':
        from . import prs_multi
        return prs_multi.report_assessment(sid)
    _, meta_hash = case_metadata(sid)
    bundle, bundle_hash = load_bundle()
    summary, score_hash = load_score(sid)
    if not source_identity_matches(sid, summary['sample']):
        raise NotReady('PRS 來源與目前個案來源不符，請重新分析')
    saved = latest_snapshot(sid)
    if not saved:
        raise NotReady('尚未計算並保存 PRS 風險，請先至 PRS 頁面完成風險計算')
    if saved.get('input_sha256') != dict(bundle=bundle_hash, score=score_hash, case_metadata=meta_hash):
        raise NotReady('PRS 計分、模型或個案資料已更新，請重新計算風險後再匯出')
    try:
        original = saved['inputs']['person']
        person = person_for(summary, {key: original[key] for key in ('age', 'sex', 'assessment_date', 'prior_chd')})
        expected = calculator.evaluate(bundle, summary, person)
        if expected != saved['result']:
            raise ValueError('Saved assessment differs from validated inputs')
    except (ValueError, KeyError, TypeError) as exc:
        raise NotReady('PRS 風險評估未通過完整性或品質檢查，請重新計算風險') from exc
    return saved


def report_readiness(sid):
    try:
        saved = report_assessment(sid)
        return dict(available=True, reason='', assessment_id=saved['assessment_id'])
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return dict(available=False, reason=str(exc))
