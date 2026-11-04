"""Local PRS queue: original gVCF -> frozen v2.1.1 score -> ancestry QC.

No clinical input belongs to these jobs. Only administrator configuration can
supply tools/resource paths; browser requests select a server-discovered gVCF.
"""
import contextlib
import fcntl
import gzip
import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import uuid

from .. import config
from . import prs, prs_paths, sample_layout, postprocessing_log, prs_presentation

ACTIVE = {'queued', 'running', 'cancelling'}
SCORING = config.REPO_ROOT / 'research/prs/scoring_v211'
DEFAULT_TEMPLATES = {
    'dragen': ['{anchor_dir}/{source}.hard-filtered.gvcf.gz',
               '{run_dir}/gvcf.gz/{source}.hard-filtered.gvcf.gz',
               '{run_dir}/other/{source}/germline_seq/{source}.hard-filtered.gvcf.gz',
               '{anchor_dir}/{source}.gvcf.gz', '{run_dir}/gvcf.gz/{source}.gvcf.gz'],
    # In-house callers/layouts must be explicitly configured and validated.
    'inhouse': [],
}


def root():
    return config.DATA_ROOT / 'jobs/prs'


def job_dir(jid):
    if not isinstance(jid, str) or not re.fullmatch('[0-9a-f]{32}', jid):
        raise ValueError('PRS 工作識別碼不正確')
    return root() / jid


def signature(path):
    p = Path(path).resolve()
    st = p.stat()
    return [str(p), st.st_size, st.st_mtime_ns, st.st_ctime_ns]


def settings():
    path = Path(os.environ.get('NGS_UI_PRS_CONFIG', config.DATA_ROOT / 'prs_config.json'))
    if not path.is_file():
        raise ValueError('尚未設定 PRS 執行環境（NGS_UI_PRS_CONFIG）')
    cfg, digest = prs.read_json(path)
    for key in ('project', 'raw_work', 'raw_script', 'original_scorefile', 'fasta', 'dbsnp',
                'plink_img', 'raw_nextflow_config', 'ancestry_nextflow_config', 'reference_panel',
                'method2_model_dir', 'method2_lock'):
        if not isinstance(cfg.get(key), str) or not Path(cfg[key]).is_absolute():
            raise ValueError(f'PRS 設定需要絕對路徑：{key}')
        if key not in ('project', 'raw_work') and not Path(cfg[key]).exists():
            raise ValueError(f'PRS 資源不存在：{key}')
    if hashlib.sha256(Path(cfg['method2_lock']).read_bytes()).hexdigest() != prs.SCORE_DEFINITION['model_lock_sha256']:
        raise ValueError('PRS model lock 與固定版本不符')
    if Path(cfg['raw_script']).read_bytes() != (SCORING/'vendor/PGS004879_FINAL_gVCF_resume.sh').read_bytes():
        raise ValueError('PRS raw script 與已驗證的版本不符')
    for tool in ('nextflow', 'apptainer', 'bcftools', 'samtools'):
        if not shutil.which(tool):
            raise ValueError(f'PRS 執行環境缺少 {tool}')
    return cfg, digest


def resource_signature(cfg):
    return {key: signature(cfg[key]) for key in (
        'raw_script', 'original_scorefile', 'fasta', 'dbsnp', 'plink_img',
        'raw_nextflow_config', 'ancestry_nextflow_config', 'reference_panel', 'method2_lock')}


def source(sid):
    prs.validate_sid(sid)
    path = sample_layout.state_file(sid, 'pipeline_source.json')
    if not path.is_file():
        raise ValueError('缺少三級分析來源紀錄，無法核對原始 gVCF')
    marker, digest = prs.read_json(path)
    source_id = prs.validate_sid(marker.get('source_sample_id'))
    anchor = Path(marker.get('source_vcf_path', ''))
    if not anchor.is_file():
        raise ValueError('三級分析的原始 VCF 來源不存在')
    seq_type = marker.get('seq_type')
    if not seq_type:
        try:
            metadata, _ = prs.case_metadata(sid)
            seq_type = metadata.get('seq_type') or metadata.get('test_type')
        except FileNotFoundError:
            pass
    if str(seq_type).upper() not in ('WGS', 'TITAN-WGS'):
        raise ValueError('PRS 目前限已確認為 WGS 的個案；舊來源缺少類型時需先確認')
    mode = marker.get('pipeline_type')
    if mode == 'nckuh': mode = 'inhouse'
    if mode not in DEFAULT_TEMPLATES:
        raise ValueError('PRS 不支援此分析來源')
    return dict(sample_id=sid, source_sample_id=source_id, anchor=str(anchor.resolve()),
                anchor_signature=signature(anchor), pipeline_source_sha256=digest, mode=mode)


def candidates(info, cfg):
    anchor = Path(info['anchor'])
    mapping = cfg.get('gvcf_templates', DEFAULT_TEMPLATES)
    if not isinstance(mapping, dict):
        raise ValueError('gvcf_templates 必須是來源對應表')
    templates = mapping.get(info['mode'], [])
    if not isinstance(templates, list) or any(not isinstance(t, str) for t in templates):
        raise ValueError('gvcf_templates 必須是各來源的路徑陣列')
    paths = set()
    for template in templates:
        # No glob or request-provided arbitrary path; precise source filenames only.
        rendered = template.format(anchor_dir=anchor.parent, run_dir=anchor.parent.parent,
                                   source=info['source_sample_id'])
        path = Path(rendered)
        if path.is_file(): paths.add(path.resolve())
    return [dict(id=hashlib.sha256(str(p).encode()).hexdigest(), path=str(p)) for p in sorted(paths)]


def validate_gvcf(path, source_id):
    path = Path(path)
    indexes = [Path(str(path)+s) for s in ('.tbi', '.csi') if Path(str(path)+s).is_file()]
    if not indexes:
        raise ValueError('原始 gVCF 缺少 tbi/csi index')
    symbolic = False
    grch38 = False
    with gzip.open(path, 'rt') as stream:
        for _ in range(100000):
            line = stream.readline(1_000_000)
            if not line: break
            if line.startswith('##contig=<') and re.search(r'ID=(?:chr)?1[,>]', line):
                grch38 = bool(re.search(r'length=248956422[,>]', line, re.IGNORECASE))
            if line.startswith('##ALT=<ID=NON_REF,') or line.startswith('##ALT=<ID=*,'):
                symbolic = True
            if line.startswith('#CHROM\t'):
                if line.rstrip('\n\r').split('\t')[9:] != [source_id]:
                    raise ValueError('gVCF header 樣本與三級來源不符')
                if not grch38:
                    raise ValueError('gVCF chr1 長度不符合 GRCh38，或未提供可核對的組裝資訊')
                if not symbolic:
                    raise ValueError('輸入未宣告 gVCF reference-block allele')
                return [signature(path), *[signature(p) for p in indexes]]
    raise ValueError('找不到合格的單樣本 gVCF header')


def read_job(jid):
    value, _ = prs.read_json(job_dir(jid)/'state.json')
    if value['status'] in ACTIVE and time.time()-value.get('updated', 0) > 30:
        pid = value.get('pid')
        if not pid:
            try: pid = int((job_dir(jid)/'pid').read_text())
            except (OSError, ValueError): pass
        from .dragen_jobs import _pid_alive
        if not pid or not _pid_alive(pid):
            value.update(status='failed', error='PRS worker 已停止；可重試，既有結果保留', updated=time.time())
            prs.atomic_json(job_dir(jid)/'state.json', value)
    return value


def jobs(sid=None):
    out = []
    for path in root().glob('*/state.json'):
        try:
            job = read_job(path.parent.name)
            if sid is None or job['sample_id'] == sid: out.append(job)
        except (OSError, ValueError, KeyError):
            continue
    return sorted(out, key=lambda j: j['created'], reverse=True)


def active_ids():
    return {sid for j in jobs() if j['status'] in ACTIVE for sid in (j['sample_id'], j['source_sample_id'])}


def public_job(job):
    if not job: return None
    keys = ('job_id', 'status', 'step', 'created', 'updated', 'error', 'warning', 'parent_job_id')
    return {k: job.get(k) for k in keys}


def overview(sid):
    history = jobs(sid)
    latest = history[0] if history else None
    choices, error = [], None
    try:
        cfg, _ = settings()
        info = source(sid)
        choices = candidates(info, cfg)
        if not choices: error = '找不到與三級來源相符的原始 gVCF'
    except (OSError, ValueError, KeyError) as exc:
        error = str(exc)
    return dict(job=public_job(latest), candidates=choices, setup_issue=error,
                can_start=not error and not (latest and latest['status'] in ACTIVE))


@contextlib.contextmanager
def submission_lock():
    root().mkdir(parents=True, exist_ok=True)
    with (root()/'.submit.lock').open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def start(sid, candidate_id=None, username='', parent_job_id=None, new_run=False, *, launch=True):
    # Called after tertiary promotion; metadata is intentionally NOT required.
    prs.validate_sid(sid)
    if not launch and not parent_job_id:
        raise ValueError('同步 PRS 必須由三級後處理啟動')
    cfg, cfg_hash = settings()
    info = source(sid)
    choices = candidates(info, cfg)
    if candidate_id is None and len(choices) == 1: candidate_id = choices[0]['id']
    selected = next((c for c in choices if c['id'] == candidate_id), None)
    if not selected: raise ValueError('請在 PRS 頁面確認原始 gVCF（無候選或有多筆來源）')
    stamps = validate_gvcf(selected['path'], info['source_sample_id'])
    effective, paths, wrappers = prs_paths.execution_config(sid, cfg)
    plan = dict(**info, gvcf=selected['path'], gvcf_signature=stamps,
                config_sha256=cfg_hash, resources=resource_signature(cfg),
                path_layout=paths, runtime_config_sha256=prs_paths.digest_json(effective))
    with submission_lock():
        from . import dragen_jobs
        if not parent_job_id and {sid, info['source_sample_id']} & dragen_jobs.active_sample_ids():
            raise ValueError('此個案已有分析工作執行中')
        if {sid, info['source_sample_id']} & active_ids():
            raise ValueError('此個案 PRS 已在排程或執行中')
        prior = None if new_run else next((j for j in jobs(sid) if j.get('plan') == plan and j['status'] not in ACTIVE), None)
        prs_paths.prepare(paths, wrappers)
        jid = uuid.uuid4().hex
        directory = job_dir(jid); directory.mkdir()
        record = dict(job_id=jid, sample_id=sid, source_sample_id=info['source_sample_id'],
                      status='queued', step='等待執行', created=time.time(), updated=time.time(),
                      requested_by=username, parent_job_id=parent_job_id, plan=plan,
                      pipeline_run_id=prior['pipeline_run_id'] if prior else 'ui_prs_'+jid)
        prs.atomic_json(directory/'state.json', record)
        prs.atomic_json(directory/'config.json', effective)
        if not launch:
            record.update(pid=os.getpid(), execution='tertiary_postprocessing')
            prs.atomic_json(directory/'state.json', record)
            (directory/'log.txt').touch()
            return public_job(record)
        try:
            _, worker_env = prs_paths.environment(paths, record['pipeline_run_id'], 'worker')
            worker_env['PYTHONPATH'] = str(config.REPO_ROOT/'backend')
            with (directory/'log.txt').open('a') as log:
                proc = subprocess.Popen([sys.executable, '-m', 'app.workers.prs_run', jid],
                    cwd=config.REPO_ROOT, env=worker_env,
                    stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            # Worker owns state.json after spawn. Separate PID file avoids races.
            (directory/'pid').write_text(str(proc.pid))
        except Exception as exc:
            record.update(status='failed', error=str(exc));prs.atomic_json(directory/'state.json', record)
            raise
        return public_job(record)


def cancel(sid, jid):
    job = read_job(jid)
    if job['sample_id'] != sid: raise FileNotFoundError('找不到 PRS 工作')
    if job['status'] in ACTIVE:
        # Polling cancellation avoids signalling an unrelated reused PID.
        (job_dir(jid)/'cancel').touch()
    return public_job(job)


def log(sid, jid):
    job = read_job(jid)
    if job['sample_id'] != sid: raise FileNotFoundError('找不到 PRS 工作')
    path = job_dir(jid)/'log.txt'
    with path.open('rb') as stream:
        stream.seek(max(0, path.stat().st_size-50000))
        tail = stream.read().decode('utf-8', errors='replace')
    text = prs_presentation.log_summary(tail)
    try:
        record, _ = prs.load_score_record(sid)
        if record.get('job_id') == jid:
            row = prs.calculator.score_row(record['summary'], prs.SCORE_DEFINITION, require_ancestry=False)
            text += (f'\n\n[表格採用的結果] 原始 PRS score={row["total_sum"]:.6f}；'
                     f'PRS percentile=第 {row["eas_percentile"]:.2f} 百分位\n'
                     '祖源／品質檢查的另外一組 SUM、PCT 僅供工具稽核，不用於表格或風險計算。')
            text += '\n'+'\n'.join(prs_presentation.score_messages(record['summary'], row))
        else:
            text += '\n\n此工作不是目前表格所採用的計分來源。'
    except (OSError, ValueError, KeyError, TypeError):
        text += '\n\n此工作尚無可供表格使用的有效計分結果。'
    text += f'\n[目前狀態] {job["status"]}；{job.get("step", "")} '
    if job.get('error') or job.get('warning'):
        text += '\n'+(job.get('error') or job['warning'])
    return dict(job=public_job(job), log=text, raw_log_available=True)


def raw_log_path(sid, jid):
    job = read_job(jid)
    if job['sample_id'] != sid: raise FileNotFoundError('找不到 PRS 工作')
    path = job_dir(jid)/'log.txt'
    if not path.is_file(): raise FileNotFoundError('找不到 PRS Log')
    return path


def failure_message(jid, fallback):
    return postprocessing_log.failure_detail(postprocessing_log.read_tail(job_dir(jid)/'log.txt'), fallback)


def enqueue_tertiary(samples, parent_job_id, logger=print):
    results = {}
    for sample in samples:
        sid = sample['sample_id']
        logger(postprocessing_log.event(sid, 'prs', 'start', detail='排入獨立背景工作'))
        try:
            results[sid] = start(sid, parent_job_id=parent_job_id)
            logger(postprocessing_log.event(sid, 'prs', 'queued',
                detail=f'job={results[sid]["job_id"]}；等待三級釋放個案鎖，計分結果另行發布'))
        except Exception as exc:
            results[sid] = dict(status='blocked', error=str(exc))
            logger(postprocessing_log.event(sid, 'prs', 'blocked', detail=f'{exc}；可由個案PRS區補跑'))
    return results
