#!/usr/bin/env python3
"""Prepare all five additional PRS models together; optionally activate them.

Run on the analysis host. No patient jobs, services or old results are changed.
Requires the existing CHD configuration and the configured reference archive.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import urllib.request

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO/'backend'))
from app import config
from app.services import prs, prs_paths, nextflow_work
from research.prs.multidisease.common import load_models, sha, write_json
from research.prs.multidisease.prepare import prepare


def download_scores(directory):
    directory.mkdir(parents=True, exist_ok=True)
    for source in json.loads((REPO/'research/prs/multidisease/sources.json').read_text()):
        path = directory/(source['pgs_id']+'_hmPOS_GRCh38.txt.gz')
        if path.is_file():
            if sha(path) != source['sha256']: raise ValueError('Existing score checksum mismatch: '+str(path))
            continue
        partial = path.with_suffix('.download')
        print('[PRS] 下載 '+source['pgs_id'], flush=True)
        with urllib.request.urlopen(source['url'], timeout=120) as response, partial.open('wb') as out:
            while block := response.read(2**20): out.write(block)
        if sha(partial) != source['sha256']: raise ValueError('Official score download checksum mismatch')
        partial.replace(path)


def extract_reference(archive, directory):
    directory.mkdir(parents=True, exist_ok=True)
    marker = directory/'complete.json'
    st = archive.stat(); stamp = [str(archive), st.st_size, st.st_mtime_ns]
    if not marker.is_file():
        print('[PRS] 解開共用參考資料；此步驟只需執行一次。', flush=True)
        proc = subprocess.Popen(['zstd', '-dc', str(archive)], stdout=subprocess.PIPE)
        try:
            with tarfile.open(fileobj=proc.stdout, mode='r|') as tar:
                for member in tar:
                    target = directory/member.name
                    if (member.issym() or member.islnk() or not (member.isdir() or member.isfile())
                            or Path(member.name).is_absolute() or not target.resolve().is_relative_to(directory.resolve())):
                        raise ValueError('Unsafe entry in reference archive')
                    if member.isdir(): target.mkdir(parents=True, exist_ok=True)
                    else:
                        target.parent.mkdir(parents=True, exist_ok=True)
                        with tar.extractfile(member) as source, target.open('wb') as out:
                            while block := source.read(2**20): out.write(block)
            if proc.wait(): raise RuntimeError('Reference decompression failed')
        finally:
            if proc.poll() is None: proc.terminate(); proc.wait()
        write_json(marker, dict(source=stamp))
    elif json.loads(marker.read_text()).get('source') != stamp:
        raise ValueError('Reference archive changed; choose a new workspace')
    candidates = [p.with_suffix('') for p in directory.rglob('*.pgen')
                  if 'GRCh38' in str(p) and Path(str(p.with_suffix(''))+'.psam').is_file()]
    if len(candidates) != 1:
        raise ValueError('無法唯一辨識 GRCh38 參考 PGEN；請用 --reference-pfile 指定完整檔名前綴')
    return candidates[0]


def activate(config_path, original_digest, output):
    # Read only state files; never terminate a job or restart a service here.
    for root in (config.DATA_ROOT/'jobs/prs', config.DATA_ROOT/'jobs/tertiary'):
        for p in root.glob('*/state.json'):
            record = json.loads(p.read_text())
            if record.get('status', record.get('state')) in ('queued','running','cancelling','starting'):
                raise ValueError('仍有三級或 PRS 工作執行中；模型已準備完成，請工作結束後再加 --activate')
    cfg, digest = prs.read_json(config_path)
    if digest != original_digest: raise ValueError('PRS config changed during preparation; retry activation')
    if cfg.get('multi_model_dir') == str(output): return
    backup = config_path.with_name(config_path.name+'.before-multi-'+datetime.datetime.now().strftime('%Y%m%d%H%M%S'))
    with backup.open('xb') as f: f.write(config_path.read_bytes())
    prs.atomic_json(config_path, {**cfg, 'multi_model_dir': str(output)})
    print('[PRS] 已啟用六疾病；新分析會使用新模型，既有個案結果保留。', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, default=Path(os.environ.get('NGS_UI_PRS_CONFIG', config.DATA_ROOT/'prs_config.json')))
    p.add_argument('--output', type=Path, default=config.NGS_UI_HOME/'biotools/prs/multi_v1')
    p.add_argument('--workspace', type=Path, default=config.PRS_WORK_ROOT/'_reference_prepare_v1')
    p.add_argument('--reference-pfile', type=Path)
    p.add_argument('--activate', action='store_true')
    a = p.parse_args(); cfg, digest = prs.read_json(a.config)
    output = a.output.resolve(); work = a.workspace.resolve()
    if not (output/'manifest.json').is_file():
        if output.exists(): raise ValueError('模型準備目錄尚未完成；請保留稽核內容，使用新的 --output 路徑重新準備')
        # Scope lease prevents the platform cleanup tool deleting preparation.
        with nextflow_work.acquire_lease(config.TERTIARY_NF_WORK_ROOT, 'prs'):
            paths = dict(project=str(work), work_root=str(config.PRS_WORK_ROOT), version=prs_paths.VERSION)
            _, env = prs_paths.environment(paths, 'prepare_v1', 'models')
            os.environ.update(env)
            download_scores(work/'scorefiles')
            reference = a.reference_pfile.resolve() if a.reference_pfile else extract_reference(Path(cfg['reference_panel']), work/'reference')
            prepare(cfg, work/'scorefiles', reference, output)
    load_models(output)
    if a.activate: activate(a.config, digest, output)
    else: print('[PRS] 模型已完成核對；加 --activate 可讓平台新工作使用六疾病。', flush=True)


if __name__ == '__main__': main()
