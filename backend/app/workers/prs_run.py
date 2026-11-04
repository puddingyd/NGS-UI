"""PRS worker for tertiary postprocessing or a standalone card run."""
import fcntl
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import traceback
import uuid

from .. import config
from ..services import prs, prs_jobs as jobs, prs_paths, nextflow_work, postprocessing_log, container_runtime
from .dragen_run import _acquire_sample_locks, _release_sample_locks, _load_secrets


class Cancelled(Exception):
    pass


class Worker:
    def __init__(self, jid, *, log_stream=None, on_update=None):
        self.directory = jobs.job_dir(jid)
        self.job = jobs.read_job(jid)
        self.cfg, _ = prs.read_json(self.directory/'config.json')
        self.plan = self.job['plan']
        self.process = None
        self.work_lease = None
        self.log_stream, self.on_update = log_stream, on_update

    def update(self, **values):
        self.job.update(values, updated=time.time(), pid=os.getpid())
        prs.atomic_json(self.directory/'state.json', self.job)
        if self.on_update:
            self.on_update(dict(self.job))

    def message(self, text):
        print(text, file=self.log_stream, flush=True)

    def event(self, status, **fields):
        self.message(postprocessing_log.timestamped(postprocessing_log.event(self.job['sample_id'], 'prs', status, **fields)))

    def check_cancel(self):
        if (self.directory/'cancel').exists(): raise Cancelled('使用者已取消 PRS 分析')

    def verify_inputs(self):
        self.check_cancel()
        info = jobs.source(self.job['sample_id'])
        if any(info[k] != self.plan[k] for k in info):
            raise ValueError('三級分析來源已更新，請重新啟動 PRS')
        if jobs.validate_gvcf(self.plan['gvcf'], self.job['source_sample_id']) != self.plan['gvcf_signature']:
            raise ValueError('原始 gVCF 或 index 已變更，請重新啟動 PRS')
        current_config, current_hash = jobs.settings()
        effective, paths, wrappers = prs_paths.execution_config(self.job['sample_id'], current_config)
        if (current_hash != self.plan['config_sha256']
                or jobs.resource_signature(current_config) != self.plan['resources']
                or paths != self.plan.get('path_layout')
                or prs_paths.digest_json(effective) != self.plan.get('runtime_config_sha256')
                or effective != self.cfg):
            raise ValueError('PRS 執行設定已變更，請重新啟動')
        if any(not p.is_file() or p.is_symlink() or p.read_text() != content for p, content in wrappers.items()):
            raise ValueError('PRS 平台 Nextflow 設定已變更，請重新啟動')

    def stop_process(self):
        if self.process is None or self.process.poll() is not None: return
        try:
            os.killpg(self.process.pid, signal.SIGTERM)
            try: self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(self.process.pid, signal.SIGKILL);self.process.wait(timeout=5)
        except ProcessLookupError:
            pass

    def check_container_tmp(self):
        paths = self.plan['path_layout']
        project = Path(paths['project'])
        directory = prs_paths.owned_directory(project, prs_paths.container_tmp(paths).relative_to(project))
        token = uuid.uuid4().hex
        with tempfile.NamedTemporaryFile(mode='w', prefix='.ngs-ui-prs-', dir=directory) as marker:
            marker.write(token); marker.flush()
            self.command(prs_paths.container_tmp_probe(self.cfg, paths, Path(marker.name), token),
                         '容器暫存掛載與寫入檢查')

    def command(self, args, step):
        self.update(status='running', step=step)
        started = time.monotonic()
        self.event('start', detail=step)
        cwd, env = prs_paths.environment(self.plan['path_layout'], self.job['pipeline_run_id'],
                                         'ancestry' if '--raw-source' in args else 'score')
        if args[0] == 'apptainer':
            env = container_runtime.without_bind_destinations(env,
                {str(prs_paths.container_tmp(self.plan['path_layout'])), '/tmp', '/var/tmp'})
        fds = ()
        if self.work_lease is not None:
            env.update(nextflow_work.lease_env(config.TERTIARY_NF_WORK_ROOT, 'prs', self.work_lease))
            fds = (self.work_lease.fileno(),)
        output = dict(stdout=self.log_stream, stderr=subprocess.STDOUT) if self.log_stream is not None else {}
        self.process = subprocess.Popen(args, cwd=cwd, env=env, pass_fds=fds, start_new_session=True, **output)
        try:
            while self.process.poll() is None:
                self.check_cancel()
                self.update()
                time.sleep(1)
            if self.process.returncode:
                reason = jobs.failure_message(self.job['job_id'], f'exit={self.process.returncode}')
                self.event('failed', elapsed=time.monotonic()-started, detail=f'{step}；{reason}')
                raise RuntimeError(f'{step}失敗（exit {self.process.returncode}）：{reason}')
            self.event('done', elapsed=time.monotonic()-started, detail=step)
        finally:
            self.stop_process()
            self.process = None

    def publish(self, original, qc=None):
        self.verify_inputs()
        summary, digest = prs.read_json(original)
        if summary.get('sample') != self.job['source_sample_id']:
            raise ValueError('計分結果樣本識別碼不符')
        prs.calculator.score_row(summary, prs.SCORE_DEFINITION, require_ancestry=False)
        qc_hash = None
        if qc:
            evidence, qc_hash = prs.read_json(qc)
            if evidence.get('sample') != summary['sample'] or evidence.get('qc_status') not in ('PASS', 'PASS_WITH_REVIEW_FLAGS'):
                raise ValueError('祖源 QC 未通過')
            # Derived envelope only; never edit the frozen upstream summary.
            summary = {**summary, 'ancestry_status': evidence}
        artifacts = dict(score=prs_paths.publish_files(
            self.job['sample_id'], self.job['job_id'], 'score', original.parent, prs_paths.SCORE_FILES))
        if artifacts['score']['files']['summary.json'] != digest:
            raise ValueError('計分結果在發布期間已變更')
        if qc:
            artifacts['ancestry'] = prs_paths.publish_files(
                self.job['sample_id'], self.job['job_id'], 'ancestry', qc.parent, prs_paths.QC_FILES)
            if artifacts['ancestry']['files']['qc.json'] != qc_hash:
                raise ValueError('祖源結果在發布期間已變更')
        self.verify_inputs()
        record = dict(schema='PRS_SCORE_LINK_V1', sample_id=self.job['sample_id'],
                      source_sample_id=summary['sample'], original_sha256=digest,
                      pipeline_source_sha256=self.plan['pipeline_source_sha256'],
                      imported_at=prs.dt.datetime.now(prs.dt.timezone.utc).isoformat(),
                      job_id=self.job['job_id'], pipeline_run_id=self.job['pipeline_run_id'],
                      gvcf_signature=self.plan['gvcf_signature'], ancestry_qc_sha256=qc_hash,
                      path_layout=self.plan['path_layout'], artifacts=artifacts,
                      runtime_config_sha256=self.plan['runtime_config_sha256'],
                      summary=summary)
        prs.atomic_json(prs.score_path(self.job['sample_id']), record)

    def run(self, *, sample_locks_held=False):
        # Use the same cleanup admission/lease protocol as tertiary Nextflow.
        with nextflow_work.acquire_lease(config.TERTIARY_NF_WORK_ROOT, 'prs') as lease:
            self.work_lease = lease
            try:
                if sample_locks_held:
                    self._run(sample_locks_held=True)
                else:
                    self._run()
            finally:
                self.work_lease = None

    def _run(self, *, sample_locks_held=False):
        self.update(status='queued', step='等待個案分析完成' if not sample_locks_held else '等待 PRS 執行名額')
        locks = []
        slot = None
        try:
            # All callers take sample ownership before the shared PRS slot.
            # A manual job must not hold that slot while waiting for a batch
            # whose parent already owns the same sample locks.
            if not sample_locks_held:
                while True:
                    self.check_cancel()
                    try:
                        locks = _acquire_sample_locks([self.job['sample_id'], self.job['source_sample_id']]);break
                    except RuntimeError:
                        self.update();time.sleep(1)
            self.update(step='等待 PRS 執行名額')
            slot = (jobs.root()/'.execution.lock').open('a')
            while True:
                self.check_cancel()
                try:
                    fcntl.flock(slot, fcntl.LOCK_EX | fcntl.LOCK_NB);break
                except BlockingIOError:
                    self.update();time.sleep(1)
            self.verify_inputs()
            self.message('[PRS] 個案工作目錄：'+self.cfg['project'])
            self.message('[PRS] 正式結果：'+str(prs.score_path(self.job['sample_id']).parent))
            self.message('[PRS] 容器暫存：'+str(prs_paths.container_tmp(self.plan['path_layout'])))
            self.check_container_tmp()
            sid, run = self.job['source_sample_id'], self.job['pipeline_run_id']
            config_path = str(self.directory/'config.json')
            dest = Path(self.cfg['project'])/'integrated_runs'/run/sid
            self.command([sys.executable, str(jobs.SCORING/'pipeline_v21.py'), 'run',
                          '--config', config_path, '--sample', sid, '--run-id', run,
                          '--mode', 'fresh', '--gvcf', self.plan['gvcf'], '--resume'],
                         'gVCF → 原始計分 → percentile')
            original = dest/'method2/summary.json'
            self.publish(original)
            # The legacy ancestry path runs against these SAME raw outputs.
            # Its normalized percentile does NOT replace Method2 EAS ECDF.
            ancestry_run = run+'_ancestry'
            try:
                self.command([sys.executable, str(jobs.SCORING/'pipeline.py'), 'run',
                              '--config', config_path, '--sample', sid, '--run-id', ancestry_run,
                              '--mode', 'reuse', '--raw-source', str(dest/'raw')],
                             '祖源與品質檢查（計分已完成）')
                qc = Path(self.cfg['project'])/'integrated_runs'/ancestry_run/sid/'qc/qc.json'
                self.publish(original, qc)
            except RuntimeError as exc:
                self.update(status='done_with_warnings', step='計分完成；祖源檢查未完成', warning=str(exc))
                return
            self.update(status='done', step='計分與祖源檢查完成；等待臨床資訊')
        finally:
            self.stop_process()
            if slot is not None: slot.close()
            _release_sample_locks(locks)


def main(jid):
    _load_secrets()
    worker = Worker(jid)
    def terminate(*_): raise Cancelled('PRS 工作已終止')
    signal.signal(signal.SIGTERM, terminate)
    signal.signal(signal.SIGINT, terminate)
    try:
        worker.run()
        return 0
    except Cancelled as exc:
        worker.stop_process();worker.update(status='cancelled', step='已取消', error=str(exc));return 1
    except Exception as exc:
        traceback.print_exc()
        worker.stop_process();worker.update(status='failed', step='失敗', error=str(exc));return 1


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1]))
