#!/usr/bin/env python3
"""Copy/verify the whole work root and bind it back at its original path.

Run on the Linux worker host, after stopping ngs-ui and waiting for all jobs.
Default is a read-only check. --apply installs the mount, never starts services.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
from datetime import datetime
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'backend'))
from app.services import nextflow_work


def run(argv, **kwargs):
    return subprocess.run(argv, check=True, text=True, **kwargs)


def check_paths(source: Path, target: Path) -> None:
    for path in (source, target):
        if not path.is_absolute() or not re.fullmatch(r'/[A-Za-z0-9_./-]+', str(path)):
            raise ValueError(f'路徑須為不含空白或特殊分隔符號的絕對路徑：{path}')
        if path == Path('/') or path.is_symlink() or path.resolve() != path:
            raise ValueError(f'不支援 root／符號連結／非正規路徑：{path}')
    if not source.is_dir() or source == Path.home():
        raise ValueError('來源必須是既有的獨立工作目錄')
    if source.name != 'nf_work' or target.name != 'nf_work' or target.parent == Path('/'):
        raise ValueError('此工具僅搬遷 nf_work 目錄，目的地須位於個人 scratch 目錄內')
    if source == target or source in target.parents or target in source.parents:
        raise ValueError('來源與目的地不能互相包含')
    if os.path.ismount(source):
        raise ValueError('來源已是掛載點；請先查明是否已完成搬遷，勿重複操作')
    if target.exists():
        raise ValueError(f'目的地已存在，保留待確認，不覆寫：{target}')
    ancestor = target.parent
    while not ancestor.exists():
        ancestor = ancestor.parent
    if ancestor.stat().st_dev != Path('/').stat().st_dev:
        raise ValueError('目的地必須位於已擴充的根檔案系統')
    if source.stat().st_dev == ancestor.stat().st_dev:
        raise ValueError('來源與目的地已在同一檔案系統，請先確認掛載')
    mounted = run(['findmnt', '-rn', '-o', 'TARGET'], capture_output=True).stdout.splitlines()
    if any(source in Path(path).parents for path in mounted):
        raise ValueError('工作區包含其他掛載點，需先個別盤點')


def check_idle(data_root: Path, service: str, source: Path) -> None:
    state = run(['systemctl', 'show', service, '-p', 'ActiveState', '--value'], capture_output=True).stdout.strip()
    if state not in {'inactive', 'failed'}:
        raise RuntimeError(f'{service} 尚未停止；先等所有分析完成，再停止服務')
    for kind in ('tertiary', 'dragen', 'prs', 'somatic'):
        for path in (data_root/'jobs'/kind).glob('*/state.json'):
            value = json.loads(path.read_text())
            if value.get('state') in {'queued', 'running', 'cancelling'}:
                raise RuntimeError(f'尚有未結束工作：{path}；先在 UI 處理，不修改狀態檔繞過檢查')
    scopes, unknown = nextflow_work.nextflow_process_scopes(source)
    if scopes or unknown:
        raise RuntimeError(f'仍有 Nextflow 程序，未搬移：scopes={sorted(scopes)}, PID={unknown}')


def check_service_paths(service: str, source: Path, data: Path) -> None:
    """Read only selected path settings from the live service, never secrets."""
    pid = run(['systemctl', 'show', service, '-p', 'MainPID', '--value'], capture_output=True).stdout.strip()
    if not pid or pid == '0':
        print('服務目前未執行；請確認先前已在服務執行時完成路徑檢查，包含自訂 PRS／DATA root。')
        return
    raw = (Path('/proc')/pid/'environ').read_bytes()
    env = dict(part.split(b'=', 1) for part in raw.split(b'\0') if b'=' in part)
    def setting(name, default):
        return Path(os.fsdecode(env.get(name.encode(), os.fsencode(default)))).resolve()
    home = setting('NGS_UI_HOME', source.parent)
    work = setting('NGS_UI_TERTIARY_NF_WORK_ROOT', home/'nf_work')
    prs = setting('NGS_UI_PRS_WORK_ROOT', work/'prs')
    actual_data = setting('NGS_UI_DATA_ROOT', home/'data')
    print(f'服務路徑：三級={work}；PRS={prs}；DATA={actual_data}')
    if work != source or not prs.is_relative_to(source):
        raise RuntimeError('服務使用其他三級／PRS工作區；先完成路徑盤點，不會把它誤當成已搬移')
    if actual_data != data:
        raise RuntimeError(f'請加上 --data-root {actual_data}，確保檢查的是服務實際工作狀態')


def rsync_command(source: Path, target: Path, *, verify=False) -> list[str]:
    args = ['rsync', '-aHAXS', '--numeric-ids', '--protect-args']
    if verify:
        args += ['--checksum', '--dry-run', '--itemize-changes', '--delete']
    return [*args, str(source)+'/', str(target)+'/']


def install_mount(source: Path, target: Path, service: str, stamp: str, *,
                  fstab: Path = Path('/etc/fstab'),
                  unit_root: Path = Path('/etc/systemd/system')) -> Path:
    original = fstab.read_text()
    for line in original.splitlines():
        fields = line.split()
        if fields and not fields[0].startswith('#') and len(fields) > 1 and fields[1] == str(source):
            raise RuntimeError('fstab 已有相同掛載目標，請先確認，未重複新增')
    dropin = unit_root/(service+'.service.d')/'50-ngs-scratch.conf'
    if dropin.exists():
        raise RuntimeError(f'服務掛載設定已存在，未覆寫：{dropin}')
    backup = source.with_name(source.name+'.hdd-backup-'+stamp)
    if backup.exists():
        raise RuntimeError(f'備份目錄已存在：{backup}')
    source_stat = source.stat()
    source.rename(backup)
    try:
        source.mkdir(mode=source_stat.st_mode & 0o7777)
        os.chown(source, source_stat.st_uid, source_stat.st_gid)
        run(['mount', '--bind', str(target), str(source)])
    except BaseException:
        if not os.path.ismount(source):
            source.rmdir()
            backup.rename(source)
        raise
    if not os.path.samefile(source, target):
        raise RuntimeError(f'掛載驗證未通過；HDD 備份保留在 {backup}')
    fstab_backup = fstab.with_name('fstab.ngs-scratch-'+stamp)
    shutil.copy2(fstab, fstab_backup)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=fstab.parent, prefix='.ngs-fstab-', delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(original.rstrip()+'\n\n# NGS-UI work on root NVMe; keep logical paths\n'
                         f'{target} {source} none bind,nofail 0 0\n')
            handle.flush(); os.fsync(handle.fileno())
        shutil.copystat(fstab, temporary)
        os.chown(temporary, fstab.stat().st_uid, fstab.stat().st_gid)
        temporary.replace(fstab)
        dropin.parent.mkdir(parents=True, exist_ok=True)
        dropin.write_text(f'[Unit]\nRequiresMountsFor={source}\n\n[Service]\n'
                          f'ExecStartPre=/usr/bin/mountpoint --quiet {source}\n'
                          f'ExecStartPre=/usr/bin/python3 -c "import os,sys; sys.exit(os.stat(\'/\').st_dev != os.stat(\'{source}\').st_dev)"\n')
        run(['systemctl', 'daemon-reload'])
    except BaseException:
        print(f'掛載已完成，但永久設定未完成。先勿啟動服務；fstab 備份：{fstab_backup}；HDD 備份：{backup}', file=sys.stderr)
        raise
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return backup


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('/home/n102968/NGS_UI/nf_work'))
    parser.add_argument('--target', type=Path, default=Path('/scratch/n102968/nf_work'))
    parser.add_argument('--data-root', type=Path, help='預設為 source 的同層 data；自訂 NGS_UI_DATA_ROOT 時須指定')
    parser.add_argument('--service', default='ngs-ui')
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args(argv)
    if sys.platform != 'linux':
        raise RuntimeError('僅在 Linux worker 主機執行')
    if not re.fullmatch(r'[A-Za-z0-9_-]+', args.service):
        raise ValueError('service 請使用不含 .service 的名稱')
    if args.apply and os.geteuid() != 0:
        raise RuntimeError('--apply 需要 sudo；不會自行啟動或停止服務')
    source, target = args.source, args.target
    check_paths(source, target)
    data = args.data_root or source.parent/'data'
    if not data.is_dir():
        raise ValueError(f'找不到平台資料目錄；自訂 NGS_UI_DATA_ROOT 時請使用 --data-root：{data}')
    print(f'整包搬移：{source} → {target}；保留原路徑與 HDD 備份', flush=True)
    print('複製指令：'+shlex.join(rsync_command(source, target)), flush=True)
    if not args.apply:
        check_service_paths(args.service, source, data)
        print('唯讀檢查完成；--apply 前須先等所有工作結束並停止服務。')
        return 0
    # Keep the same lock inodes in the HDD parent, throughout rsync and bind.
    os.umask(0o002)  # New guard files must remain writable by the service group.
    with nextflow_work.admission(source), ExitStack() as locks:
        for scope in ('dragen', 'nckuh', 'prs'):
            locks.enter_context(nextflow_work.open_lock(nextflow_work._lock_path(source, scope), blocking=False))
        check_idle(data, args.service, source)
        source_stat = source.stat()
        if target.parent.exists():
            if target.parent.stat().st_uid != source_stat.st_uid:
                raise RuntimeError('目的地父目錄擁有者與工作區不同；請先確認個人 scratch 目錄')
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            os.chown(target.parent, source_stat.st_uid, source_stat.st_gid)
            target.parent.chmod(0o700)
        # du accounts for sparse files and in-tree hardlinks; final rsync may
        # still fail safely if the live filesystem changes or needs more space.
        used = int(run(['du', '-s', '-B1', str(source)], capture_output=True).stdout.split()[0])
        disk = shutil.disk_usage(target.parent)
        if disk.free < used + 10*1024**3:
            raise RuntimeError('SSD 可用空間不足以複製並保留至少 10 GiB 系統空間')
        print(f'工作區 {used/1024**3:.1f} GiB；目的地可用 {disk.free/1024**3:.1f} GiB', flush=True)
        run(rsync_command(source, target))
        differences = run(rsync_command(source, target, verify=True), capture_output=True).stdout
        if differences.strip():
            raise RuntimeError('checksum／權限／連結比對不一致，來源未切換：\n'+differences)
        backup = install_mount(source, target, args.service, datetime.now().strftime('%Y%m%d-%H%M%S'))
        print(f'搬遷與永久掛載設定完成。HDD 備份：{backup}', flush=True)
        run(['findmnt', '-T', str(source)])
        run(['df', '-hT', '/', str(source)])
        print(f'服務尚未啟動；確認後執行 sudo systemctl start {args.service}。', flush=True)
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(f'停止搬遷：{error}', file=sys.stderr)
        raise SystemExit(1)
