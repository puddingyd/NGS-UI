"""Platform adapter for complete, sample-bound multi-disease PRS releases."""
import importlib
import os
from pathlib import Path
import uuid

from .. import config
from . import prs, prs_presentation

calculator = importlib.import_module(prs._package+'.multi')
registry = importlib.import_module(prs._package+'.registry').DISEASES
RELEASE_SHA256 = '66949cab5117319f23db6dc7406f3d5a00e819e7d9de8e24b86600f47fa2e823'
BUNDLE = config.REPO_ROOT/'research/prs/releases/tw_multi_v1/bundle.json'


def bundle():
    value, digest = prs.read_json(BUNDLE)
    if digest != RELEASE_SHA256: raise prs.NotReady('六疾病風險參數包版本不符')
    calculator.validate_pack(value)
    return value, digest


def check_model(summary):
    path = Path(os.environ.get('NGS_UI_PRS_CONFIG', config.DATA_ROOT/'prs_config.json'))
    cfg, _ = prs.read_json(path)
    model_dir = cfg.get('multi_model_dir')
    if not model_dir: raise prs.NotReady('尚未設定六疾病參考計分模型')
    manifest, digest = prs.read_json(Path(model_dir)/'manifest.json')
    if summary.get('model_manifest_sha256') != digest:
        raise prs.NotReady('PRS 計分模型已更新，請重新分析')
    if set(manifest.get('models', {})) != set(registry): raise prs.NotReady('PRS 參考模型不完整')
    for key, row in calculator.score_rows(summary).items():
        model = manifest['models'][key]
        if row['pgs_id'] != model['pgs_id'] or row['n_loci'] != model['n_loci']:
            raise prs.NotReady('PRS 計分位點與參考模型不一致')


def status(sid, record, score_hash):
    _, meta_hash = prs.case_metadata(sid)
    summary = record['summary']; issues = []; pack_hash = None
    rows = calculator.score_rows(summary)
    if not prs.source_identity_matches(sid, summary['sample']): raise prs.NotReady('PRS 個案來源不符')
    try:
        check_model(summary); _, pack_hash = bundle()
    except (OSError, ValueError, KeyError, TypeError) as exc: issues.append(str(exc))
    ancestry = calculator.ancestry_ready(summary)
    if not ancestry: issues.append('祖源與品質檢查尚未符合風險估算要求')
    if not any(r['qc_ready'] for r in rows.values()): issues.append('尚無通過計分品質檢查的疾病')
    diseases = []
    for key, disease in registry.items():
        r = rows[key]
        score = dict(raw_score=r['total_sum'], percentile=r['eas_percentile'], group=r['intervene_group'],
            job_id=record.get('job_id'), qc_ready=ancestry and r['qc_ready'],
            observed_loci=r['observed_loci'], missing_loci=r['missing_loci'])
        diseases.append(dict(id=key, name=prs_presentation.disease_name(key, disease['name']), pgs_id=disease['pgs_id'], score=score))
    current = dict(bundle=pack_hash, score=score_hash, case_metadata=meta_hash)
    latest = prs.latest_snapshot(sid)
    stale = bool(latest and latest.get('input_sha256') != current)
    from . import prs_jobs
    return dict(available=not issues, issues=issues, diseases=diseases, score=diseases[0]['score'],
        assessment_schema='multi', latest=latest if not stale and not issues else None, stale=stale,
        research_only=True, pipeline=prs_jobs.overview(sid), explanation=explanation(),
        age_min=20, age_max_exclusive=80)


def assess(sid, payload, username):
    _, meta_hash = prs.case_metadata(sid)
    pack, pack_hash = bundle(); summary, score_hash = prs.load_score(sid)
    if not prs.source_identity_matches(sid, summary['sample']): raise prs.NotReady('PRS 個案來源不符')
    check_model(summary)
    person = calculator.person_for(summary, payload)
    result = calculator.evaluate(pack, summary, person)
    if not any(r['available'] for r in result['diseases'].values()):
        raise ValueError('尚無可計算風險的疾病；請確認病史、性別與計分品質')
    check_model(summary)
    if (prs.load_score(sid)[1] != score_hash or prs.case_metadata(sid)[1] != meta_hash
            or bundle()[1] != pack_hash or not prs.source_identity_matches(sid, summary['sample'])):
        raise prs.NotReady('輸入資料已更新，請重新計算')
    saved = dict(schema='PRS_ASSESSMENT_V1', sample_id=sid, assessment_id=uuid.uuid4().hex,
        created_at=prs.dt.datetime.now(prs.dt.timezone.utc).isoformat(), created_by=username,
        input_sha256=dict(bundle=pack_hash, score=score_hash, case_metadata=meta_hash),
        inputs=dict(person=person, summary=summary, bundle=pack), result=result)
    prs.atomic_json(prs.assessment_dir(sid)/f"{saved['assessment_id']}.json", saved)
    return saved


def report_assessment(sid):
    _, meta_hash = prs.case_metadata(sid)
    pack, pack_hash = bundle(); summary, score_hash = prs.load_score(sid)
    if not prs.source_identity_matches(sid, summary['sample']): raise prs.NotReady('PRS 個案來源不符')
    check_model(summary)
    saved = prs.latest_snapshot(sid)
    if not saved or saved.get('input_sha256') != dict(bundle=pack_hash, score=score_hash, case_metadata=meta_hash):
        raise prs.NotReady('請先完成並保存目前資料的 PRS 風險計算')
    prs.require_twenty_year_assessment(saved)
    expected = calculator.evaluate(pack, summary, saved['inputs']['person'])
    if expected != saved['result'] or not any(r['available'] for r in expected['diseases'].values()):
        raise prs.NotReady('PRS 風險評估完整性檢查未通過')
    return saved


def explanation():
    common = prs_presentation.explanation()
    return dict(method=[
        'gVCF → PRS score：共用一次基因型重建，依各疾病模型的固定變異位點與權重加總。無法可靠判讀的位點，以固定 3,330 名參考者的平均基因型補值；品質不足的疾病不估算風險。',
        '使用的模型：'+ '；'.join(prs_presentation.disease_name(k, d['name'])+' '+d['pgs_id'] for k,d in registry.items())+'。',
        'score → percentile：各疾病分別與相同 735 名東亞參考者的該疾病分數比較；不能跨疾病比較原始分數。'+prs_presentation.RANK_EXPLANATION,
        'percentile → risk：依 INTERVENE 愛沙尼亞資料庫的百分位分組相對風險，結合台灣 GBD 2019 年齡、性別別的發生率、盛行率與死亡率；估算首次發病，並考慮先因其他原因死亡的可能。',
    ], limitations=[
        common['limitations'][0], common['limitations'][1],
        prs_presentation.CANCER_NOTE,
        '只估算評估時仍存活、確認尚未罹患該疾病的 20 至未滿 80 歲個案；各疾病病史分別確認。超過 60 歲不提供完整二十年風險，超過 70 歲不提供完整十年風險，仍可查看至 80 歲風險。乳癌風險模型僅適用女性。',
        *[prs_presentation.disease_name(k, d['name'])+'：'+d['limitation'].replace('整體癌症', '全部癌症') for k,d in registry.items()],
        '全部癌症與乳癌的事件範圍重疊，兩個風險不能相加；此模型也不能用來預測特定癌別、腫瘤位置或癌症復發。',
        common['limitations'][4], common['limitations'][5],
        '尚未納入血壓、血脂、血糖、BMI、吸菸、家族史或治療效果，也沒有個人風險誤差範圍；相同百分位分組、年齡及性別會得到相同估計。',
    ], sources=[dict(label=d['pgs_id'], url='https://www.pgscatalog.org/score/'+d['pgs_id']+'/') for d in registry.values()]
       +[dict(label='INTERVENE 研究及疾病定義', url='https://pmc.ncbi.nlm.nih.gov/articles/PMC11169548/')])
