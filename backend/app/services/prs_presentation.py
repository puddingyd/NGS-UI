"""Human-readable presentation only; never modify locked scores or risk inputs."""
import re
from decimal import Decimal, ROUND_HALF_UP

RANK_EXPLANATION = ('PRS 百分位表示分數在參考族群中的排名，不是發病機率。風險依七個百分位區間估算，'
    '跨到另一區間時數字可能跳動；來源研究的相鄰組估計也不一定依序增加。')
BOUNDARY_NOTE = '此百分位接近風險分組邊界；排名小幅變動可能改變分組與風險數字。這是分組穩定性提醒，不是個人風險的誤差區間。'
CANCER_NOTE = '全部癌症：權重分數與來源相對風險分組的方向對應仍待確認；目前風險屬研究性重建，不能把高百分位直接解讀為已確認的高發病風險。'
ZERO_NOTE = '來源資料在此年齡設為零，無法據此確認零風險。'
PARTIAL_ZERO_NOTE = '本次期間有部分年齡使用來源資料設定的零發生率，估計受此限制；其餘年齡仍依模型計算，不代表零值年齡不會發病。'


def probability_text(value, digits=2):
    """Display precision is explicit; a positive risk must not look exactly zero."""
    number = Decimal(str(value)) * 100
    unit = Decimal(1).scaleb(-digits)
    rounded = number.quantize(unit, rounding=ROUND_HALF_UP)
    return ('<' + str(unit) if number > 0 and rounded == 0 else str(rounded)) + '%'


def _zero_intervals(bundle, disease, sex, start, end):
    # Only the audited GBD2019 AF source exception is a source-zero label.
    if disease != 'AF' or bundle.get('baseline_release') != 'GBD2019':
        return []
    return [[max(start, r['age_start']), min(end, r['age_end'])]
        for r in bundle['bins'] if r['sex'] == sex and 20 <= r['age_start'] < r['age_end'] <= 30
        and r['incidence'] == 0 and min(end, r['age_end']) > max(start, r['age_start'])]


def risk_annotation(bundle, disease, sex, start, endpoint, value):
    intervals = _zero_intervals(bundle, disease, sex, start, endpoint)
    coverage = 'none'
    if intervals:
        cursor = start
        for a,b in sorted(intervals):
            if a > cursor: break
            cursor = max(cursor,b)
        coverage = 'all' if cursor >= endpoint else 'partial'
    # Initial curve point is the conditioning origin, not a zero-risk forecast.
    origin = endpoint == start
    return dict(source_zero=coverage, intervals=intervals,
        ui_text='評估起點' if origin else '來源零值*' if coverage == 'all' else probability_text(value),
        report_text='來源零值*' if coverage == 'all' else probability_text(value, 1),
        message=ZERO_NOTE if coverage == 'all' else PARTIAL_ZERO_NOTE if coverage == 'partial' else '')


def assessment_presentation(snapshot):
    """Derived view of a validated snapshot; never write into its inputs/result."""
    from . import prs
    result = snapshot['result']; person = result['person']; pack = snapshot['inputs']['bundle']
    rows = result.get('diseases') or {'CHD': {**result, 'available': True}}
    views = {}
    for key, r in rows.items():
        if not r.get('available'):
            views[key] = dict(notices=[], risks={}, curve=[])
            continue
        bundle = pack['models'][key] if 'models' in pack else pack
        start, sex = person['age'], person['sex']
        notices = [CANCER_NOTE] if key == 'ALL_CANCERS' else []
        if min(abs(r.get('percentile', 50)-b) for b in (20,40,60,80,90,95)) <= 2:
            notices.append(BOUNDARY_NOTE)
        risks = {}; inversions = []
        for field, end in (('risk_10_year', start+10), ('risk_20_year', start+20), ('risk_to_80', 80)):
            value = r.get(field)
            if value is None: continue
            annotation = risk_annotation(bundle, key, sex, start, end, value['risk'])
            risks[field] = annotation
            if annotation['message'] and annotation['message'] not in notices:
                notices.append(annotation['message'])
            comparison = [prs.engine.integrate(bundle, sex, p, start, end)['risk'] for p in (10,30,50,70,85,92,97)]
            if any(a > b+1e-12 for a,b in zip(comparison,comparison[1:])):
                inversions.append({'risk_10_year':'10年','risk_20_year':'20年','risk_to_80':'至80歲'}[field])
        if inversions:
            notices.append('此疾病的'+ '／'.join(inversions)+'估計在部分相鄰分組未依百分位遞增；不能解讀為較高 PRS 具有保護作用。')
        curve = []
        for point, reference in zip(r['risk_to_80']['curve'], r['reference']['risk_to_80']['curve']):
            a = risk_annotation(bundle, key, sex, start, point['age'], point['risk'])
            a.update(age=point['age'], reference_text=risk_annotation(bundle, key, sex, start, reference['age'], reference['risk'])['ui_text'])
            curve.append(a)
        reference = {field:risk_annotation(bundle,key,sex,start,end,r['reference'][field]['risk'])
            for field,end in (('risk_10_year',start+10),('risk_20_year',start+20),('risk_to_80',80)) if r['reference'].get(field) is not None}
        views[key] = dict(notices=notices, risks=risks, curve=curve, reference=reference)
    return dict(version='PRS_PRESENTATION_V2', diseases=views)


def presented_snapshot(snapshot):
    return {**snapshot, 'presentation': assessment_presentation(snapshot)}


def disease_name(key, source_name):
    """Current display label without rewriting reproducible assessment inputs."""
    return '全部癌症（All cancers）' if key == 'ALL_CANCERS' else source_name


def nextflow_failed(raw):
    # A task can be retried by Nextflow; only its terminal workflow error is
    # grounds for stopping the owned process group.
    raw = re.sub(r'\x1b\[[0-9;?]*[A-Za-z]', '', raw)
    return bool(re.search(r'^\s*ERROR\s*~\s*Error executing process\s*>', raw, re.MULTILINE))


def failure_reason(raw, fallback):
    if 'ScoreHeader' in raw and re.search(r'\btrait_reported\s*\n\s*Field required', raw):
        return '祖源計分檔格式錯誤：缺少必要的疾病名稱欄位 trait_reported'
    from . import postprocessing_log
    detail = postprocessing_log.failure_detail(raw, fallback)
    if detail != fallback: return detail
    found = re.findall(r'^\s*(ERROR\s*~[^\n]+)', raw, re.MULTILINE)
    return found[-1][:1500] if found else fallback


def explanation():
    return dict(method=[
        'gVCF → PRS score：使用 PGS004879（INTERVENE_MegaPRS_CHD）的權重，固定計算 610,305 個已核對位點；無法可靠判讀的位點，以 3,330 名參考者的平均基因型補值後加總。',
        'PRS score → percentile：與 HGDP／1000 Genomes 中固定的 735 名東亞參考者比較排名。第 60 百分位表示分數約高於六成參考者，不是有 60% 的發病機率。',
        'percentile → 風險：通過祖源與品質檢查後，依百分位所屬的七個分組，套用 INTERVENE 愛沙尼亞生物資料庫（EstBB）的相對風險，再結合台灣 GBD 2019 各年齡、性別的發生率、盛行率與死亡率，估算首次 CHD 風險；同時考慮可能先因其他原因死亡。',
    ], limitations=[
        '僅供研究：尚未用台灣個案的實際追蹤結果，驗證預測與真實發病人數是否一致，不能單憑這個結果決定用藥或篩檢。',
        '跨族群推估：相對風險來自歐洲祖源的 EstBB；百分位來自 735 名東亞參考者，並非台灣人口常模。固定計分位點也與原始 PGS 的全部位點略有不同，這些差異對台灣預測的影響仍待驗證。',
        '疾病範圍不完全相同：INTERVENE 的 CHD 以診斷碼辨識，包含不穩定心絞痛、心肌梗塞、其他急性或慢性缺血性心臟病及心跳停止；GBD 的缺血性心臟病（IHD）則依急性心肌梗塞與慢性 IHD（如穩定心絞痛、心肌梗塞後無症狀 IHD）估計。納入事件與統計方式不同，因此這不是只預測心肌梗塞、心跳停止或心血管死亡的風險。',
        '適用對象與期間：僅適用評估時仍存活、確認尚未罹患 CHD 的 20 至未滿 80 歲個案，不估算已患病者的復發風險。超過 60 歲不提供完整二十年風險，超過 70 歲不提供完整十年風險；仍可查看至 80 歲風險，不把不足期間當成完整十年或二十年。',
        '固定年度假設：使用台灣 2019 年的人口統計率，假設未來仍維持相同年齡別情況；未預測醫療進步或人口發病率的改變。',
        '計算近似：每個五歲年齡區間內假設發病與死亡速率固定，再從目前年齡逐段累計；與原論文每五年計算一次的近似方式可能有小幅差異。',
        '尚未納入個人血壓、血脂、吸菸、糖尿病或治療效果；也尚未提供風險估計的誤差範圍。同一百分位分組、年齡與性別會得到相同風險，不是完整的個人心血管臨床評估。',
        '分數可先顯示，風險須等檢查：計分完成就能顯示 score 與 percentile；祖源不適合東亞參考比較，或品質檢查未符合要求時，不能用本流程計算風險。接近分組邊界或重要位點缺失時，仍需查看品質提醒。',
    ], sources=[
        dict(label='PGS004879 模型', url='https://www.pgscatalog.org/score/PGS004879/'),
        dict(label='INTERVENE 原研究與疾病定義附件 8、9', url='https://pmc.ncbi.nlm.nih.gov/articles/PMC11169548/'),
        dict(label='GBD 2019 IHD 定義', url='https://pubmed.ncbi.nlm.nih.gov/34922374/'),
    ])


def log_summary(raw):
    """Keep progress and diagnostics, not the legacy pipeline's other scores.

    Original stdout/stderr remains downloadable byte-for-byte. In particular,
    qc_report JSON and COMPLETED SUM/PCT must not look like the table result.
    """
    lines = ['流程摘要；表格計分與祖源／品質檢查分開執行。完整工具輸出可另行下載。']
    error_context = False
    for line in raw.splitlines():
        text = re.sub(r'\x1b\[[0-9;?]*[A-Za-z]', '', line).strip()
        if re.match(r'(?:Traceback|.*ValidationError:)', text): error_context = True
        if re.match(r'^(?:Work dir:|Container:|Tip:|executor >)', text): error_context = False
        if (re.search(r'\[postprocessing\]|\[PRS\]', text)
                or re.match(r'(?:executor >|\[(?:[0-9a-f]{2}/[0-9a-f]{6}|-\s+)\])', text)
                or re.search(r'Pipeline completed successfully|^(?:Completed at:|Duration\s*:|Succeeded\s*:)', text)
                or re.search(r'\b(?:FATAL|ERROR|WARNING)\b|(?:Error|Exception):|Traceback', text)
                or text == 'Execution cancelled -- Finishing pending tasks before exit'
                or error_context):
            lines.append(text)
    return '\n'.join(lines)


def score_messages(summary, row):
    messages = [f'固定計分位點：{row["n_loci"]:,}；可靠判讀：{row["observed_loci"]:,}；參考平均補值：{row["missing_loci"]:,}。']
    if row.get('boundary_review_2pp'):
        messages.append('提醒：表格百分位距離風險分組邊界不超過 2 個百分點，接近邊界，需留意分組穩定性。')
    missing = summary.get('high_weight_missing_modelA') or []
    if missing:
        messages.append(f'提醒：固定計分中有 {len(missing)} 個較高權重位點需補值，請查看完整計分稽核。')
    if summary.get('strand_ambiguity_warnings'):
        messages.append('提醒：計分包含正反股方向的研究假設，請查看完整計分稽核。')
    qc = summary.get('ancestry_status', {})
    messages.append(f'祖源／品質檢查：{qc.get("qc_status") or "尚未完成"}。')
    if qc.get('flags'):
        messages.append(f'祖源工具另有 {len(qc["flags"])} 項檢查提醒，詳見完整工具 Log；其輔助分數不取代表格結果。')
    return messages
