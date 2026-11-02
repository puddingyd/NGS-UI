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
        with os.fdopen(fd, 'w') as stream:
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


def load_score(sid):
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
        raise NotReady('個案分析來源已更新，請重新確認並匯入 PRS 結果')
    return summary, digest


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
    _, meta_hash = case_metadata(sid)
    issues, bundle, summary = [], None, None
    bundle_hash = score_hash = None
    try:
        bundle, bundle_hash = load_bundle()
    except (OSError, ValueError) as exc:
        issues.append(str(exc))
    try:
        summary, score_hash = load_score(sid)
        if not source_identity_matches(sid, summary['sample']):
            raise NotReady('PRS 來源與目前個案來源不符，請重新匯入')
    except (OSError, ValueError, KeyError, TypeError):
        summary = None
        issues.append('尚未匯入有效且與此個案相符的 PRS 結果')
    score = None
    if bundle and summary:
        try:
            # Probe validates score/QC only. Never return this placeholder age's risk.
            probe = calculator.evaluate(bundle, summary, dict(sample=summary['sample'], age=50,
                                        sex='female', assessment_date=dt.date.today().isoformat(), prior_chd=False))
            score = {'percentile': probe['percentile'], 'group': probe['risk_to_80']['group']}
        except (ValueError, KeyError, TypeError):
            issues.append('PRS 品質、祖源或計分版本尚未符合風險模型要求')
    current = dict(bundle=bundle_hash, score=score_hash, case_metadata=meta_hash)
    latest = latest_snapshot(sid)
    stale = bool(latest and latest.get('input_sha256') != current)
    return dict(available=not issues, issues=issues, score=score, research_only=True,
                age_min=20, age_max_exclusive=80,
                model_id=bundle['model_id'] if bundle else None,
                latest=latest if not stale and not issues else None, stale=stale)


def assess(sid, payload, username):
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
