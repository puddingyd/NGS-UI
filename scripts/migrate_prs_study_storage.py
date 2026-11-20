#!/usr/bin/env python3
"""Copy a stopped PRS study to local HDD; bind it at the unchanged study path.

Linux/root only. Default: check. --apply copies/verifies/installs the bind,
retaining the SSD backup. --release-ssd separately verifies and removes ONLY
that recorded backup. Never starts/stops services or changes scientific files.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager, ExitStack
from datetime import datetime, timezone
import errno
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import stat
import subprocess
import sys
import tempfile

# A pinned standalone copy may be run from /tmp while the checkout stays
# untouched during an active tertiary job. In that case require the repo cwd.
REPO = Path(__file__).resolve().parents[1]
if not (REPO / 'backend/app/services/nextflow_work.py').is_file(): REPO = Path.cwd()
sys.path.insert(0, str(REPO / 'backend'))
from app.services import nextflow_work as guard

SCHEMA = 'NGS_UI_PRS_STORAGE_MIGRATION_V1'
HELD_FDS = []
FSTAB = Path('/etc/fstab')


def run(argv, *, capture=False):
    child = subprocess.Popen(list(map(str, argv)), text=True, start_new_session=True,
                             stdout=subprocess.PIPE if capture else None,
                             stderr=subprocess.PIPE if capture else None,
                             pass_fds=tuple(HELD_FDS))
    try:
        out, err = child.communicate()
        if child.returncode:
            raise RuntimeError(f'{Path(str(argv[0])).name} failed ({child.returncode}): {err or ""}')
        return out or ''
    finally:
        if child.poll() is None:
            try:
                os.killpg(child.pid, signal.SIGTERM)
                try: child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL); child.wait()
            except ProcessLookupError:
                child.wait()  # Child exited between poll and the group signal.


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024**2), b''): h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def atomic(path, value, *, mode=0o600, owner=None):
    path = Path(path)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as f:
        tmp = Path(f.name)
        try:
            os.fchmod(f.fileno(), mode)
            if owner: os.fchown(f.fileno(), *owner)
            f.write(value); f.flush(); os.fsync(f.fileno())
            tmp.replace(path)
        finally:
            tmp.unlink(missing_ok=True)
    fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try: os.fsync(fd)
    finally: os.close(fd)


def save(path, value):
    value['updated_at'] = datetime.now(timezone.utc).isoformat()
    atomic(path, json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def canonical(path):
    path = Path(path)
    if (not path.is_absolute() or not re.fullmatch(r'/[A-Za-z0-9_./-]+', str(path))
            or '..' in path.parts or path.resolve() != path or path.is_symlink()):
        raise ValueError(f'路徑必須是無 symlink 的正規絕對路徑：{path}')
    return path


def shape(work, study, target):
    for p in (work, study, target): canonical(p)
    try: parts = study.relative_to(work).parts
    except ValueError: raise ValueError('study 必須位於既有 nf_work/prs 下')
    if (len(parts) != 3 or parts[0] != 'prs' or not parts[1].startswith('_validation_')
            or not re.fullmatch(r'[A-Za-z0-9_-]+', parts[2])):
        raise ValueError('只允許 nf_work/prs/_validation_*/study，不能搬整個 work 或 UI 個案')
    if (target == study or target.is_relative_to(work) or work.is_relative_to(target)
            or not target.is_relative_to(Path('/home'))
            or target.parts[1:2] != ('home',) or target.parent.name != parts[1]
            or target.name != study.name):
        raise ValueError('目的地須在 work 外的 /home，且保留 _validation_*/study 名稱')


def mount_rows(path=Path('/proc/self/mountinfo')):
    rows = []
    def unescape(value):
        return re.sub(r'\\([0-7]{3})', lambda m: chr(int(m[1], 8)), value)
    for line in path.read_text().splitlines():
        fields = line.split(); split = fields.index('-')
        rows.append(dict(target=unescape(fields[4]), fsroot=unescape(fields[3]), fstype=fields[split + 1]))
    return rows


def aliases(work, study, mounts):
    result = {study}
    parents = [r for r in mounts if study.is_relative_to(Path(r['target']))]
    if parents:
        r = max(parents, key=lambda r: len(Path(r['target']).parts))
        candidate = Path(r['fsroot']) / study.relative_to(Path(r['target']))
        if candidate.exists() and os.path.samefile(candidate, study): result.add(candidate)
    return result


def no_nested_mounts(paths, mounts):
    for r in mounts:
        p = Path(r['target'])
        if any(p != root and p.is_relative_to(root) for root in paths):
            raise RuntimeError('study／備份內仍有子掛載，停止搬移或清除：' + str(p))


def process_users(paths, proc=Path('/proc')):
    """Check paths in command/env/cwd/fds, including the physical bind alias.

    Inspect only in memory; never print command arguments or environment secrets.
    Exempt this tool and its ancestors (sudo/shell legitimately contain paths).
    """
    exempt = set(); pid = os.getpid()
    while pid and pid not in exempt:
        exempt.add(pid)
        try:
            status = (proc / str(pid) / 'status').read_text()
            pid = int(re.search(r'^PPid:\s+(\d+)', status, re.M).group(1))
        except (OSError, AttributeError): break
    needles = [str(p).encode() for p in paths]
    def mentions(raw):
        return any(re.search(re.escape(p) + rb'(?=$|[/\x00\s])', raw) for p in needles)
    def optional_read(path):
        try: return path.read_bytes()
        except OSError as error:
            # Kernel threads can have a proc entry but no userspace mm/environ
            # (ESRCH). A normal task can also exit between directory and file
            # reads (ENOENT/ESRCH). Still inspect the remaining evidence: an
            # unavailable environ must not hide a live cwd/fd in this study.
            if error.errno in {errno.ENOENT, errno.ESRCH}: return b''
            raise
    busy = []
    for entry in proc.iterdir():
        if not entry.name.isdigit() or int(entry.name) in exempt: continue
        try:
            hit = mentions(optional_read(entry / 'cmdline')) or mentions(optional_read(entry / 'environ'))
            if not hit:
                try: fds = list((entry / 'fd').iterdir())
                except OSError as error:
                    if error.errno not in {errno.ENOENT, errno.ESRCH}: raise
                    fds = []
                links = [entry / 'cwd', *fds]
                for link in links:
                    try:
                        if mentions(os.fsencode(os.readlink(link))): hit = True; break
                    except OSError as error:
                        if error.errno not in {errno.ENOENT, errno.ESRCH}: raise
            if hit: busy.append(entry.name)
        except FileNotFoundError: continue  # Process exited during the scan.
        except PermissionError:
            if entry.exists(): raise RuntimeError('無法核對程序存取；請以 sudo 執行')
    if busy: raise RuntimeError('仍有程序使用 study，未搬移／清除；PID=' + ','.join(busy))


def rsync(source, target, *, verify=False):
    args = ['rsync', '-aHAXS', '--numeric-ids', '--protect-args']
    if verify: args += ['--checksum', '--dry-run', '--itemize-changes', '--delete']
    else: args += ['--info=progress2', '--human-readable']
    return [*args, str(source) + '/', str(target) + '/']


def verify_copy(source, target):
    differences = run(rsync(source, target, verify=True), capture=True)
    if differences.strip():
        raise RuntimeError('內容／metadata 比對不同，保留來源與目的地，未刪除任何原檔：\n' + differences[:4000])


def verify_remaining(source, target):
    # After a partial rmtree, target legitimately contains more files and the
    # remaining source directories have newer mtimes. Every remaining file
    # must still match; a newly added/changed backup file must not be erased.
    args = rsync(source, target, verify=True)
    args.remove('--delete'); args.insert(1, '--omit-dir-times')
    differences = run(args, capture=True)
    if differences.strip(): raise RuntimeError('剩餘 SSD 備份與 HDD 不同，未繼續刪除：\n' + differences[:4000])


def check_plan(study, logical, expected_hash=None):
    if expected_hash and digest(study / 'plan.json') != expected_hash:
        raise RuntimeError('study plan 改變，停止操作')
    plan = read(study / 'plan.json')
    if plan.get('schema') != 'PRS_COHORT_STUDY_V1' or plan.get('study') != str(logical):
        raise ValueError('不是相符的既有 cohort study；不重寫 plan 或設定')
    manifest = study / 'frozen/manifest.json'
    if digest(manifest) != plan['model_sha256']:
        raise ValueError('凍結模型 manifest 校驗失敗')
    if digest(study / 'frozen/reference_scores.json') != read(manifest)['sha256']['reference_scores.json']:
        raise ValueError('凍結參考分布校驗失敗')
    return plan


def fstab_text(original, source, target):
    wanted = f'{target} {source} none bind,nofail 0 0'
    for line in original.splitlines():
        fields = line.split()
        if fields and not fields[0].startswith('#') and len(fields) > 1 and fields[1] == str(source):
            if fields != wanted.split(): raise RuntimeError('fstab 有不同的 study 掛載設定，未覆寫')
            return original
    return original.rstrip() + '\n\n# NGS-UI PRS study on HDD; keep frozen logical paths\n' + wanted + '\n'


def persist(source, target, fstab, stamp):
    old = fstab.read_text(); new = fstab_text(old, source, target)
    if new != old:
        backup = fstab.with_name(fstab.name + '.prs-storage-' + stamp)
        if backup.exists():
            if backup.read_text() != old: raise RuntimeError('fstab 備份與現況不同，未覆寫')
        else: shutil.copy2(fstab, backup)
        info = fstab.stat()
        atomic(fstab, new, mode=stat.S_IMODE(info.st_mode), owner=(info.st_uid, info.st_gid))
    run(['systemctl', 'daemon-reload'])  # Does not start, stop, or restart services.


@contextmanager
def lockset(work, study, target, journal, *, create=False):
    """Scope stays shared: the already-running NCKUH scope is unaffected.

    Study + global PRS execution + destination locks exclude writers. Hold BOTH
    copies of execution.lock over cutover, preventing the changed-inode gap.
    """
    global HELD_FDS
    try:
        with ExitStack() as stack:
            def hold(handle):
                stack.callback(handle.close); HELD_FDS.append(handle.fileno()); return handle
            hold(guard.acquire_lease(work, 'prs'))
            for p in (work / 'prs/.study-storage.lock', work / 'prs/.execution.lock', study / 'execution.lock'):
                if p.is_symlink(): raise ValueError('鎖檔不可為 symlink')
                hold(guard.open_lock(p, blocking=False))
            if create and not journal.exists():
                if target.exists(): raise ValueError('目的地已存在且無本工具的紀錄，未覆寫')
            yield hold
    finally:
        HELD_FDS = []


def identity(path):
    s = path.stat(); return [s.st_dev, s.st_ino]


def mkdir_owned(path, owner):
    missing = []; p = path
    while not p.exists(): missing.append(p); p = p.parent
    for p in reversed(missing):
        p.mkdir(mode=0o700); os.chown(p, *owner)


def install(source, target, record, journal, *, fstab=None):
    fstab = fstab or FSTAB
    backup = Path(record['backup'])
    if backup.exists(): raise RuntimeError('SSD 備份已存在；保留資料，先確認前次切換狀態')
    source.rename(backup)
    record['phase'] = 'switching'; save(journal, record)
    try:
        source.mkdir(mode=0o000)  # Missing bind fails closed; never silently writes to SSD.
        run(['mount', '--bind', target, source])
        if not os.path.samefile(source, target): raise RuntimeError('掛載身分核對失敗')
    except BaseException:
        if not os.path.ismount(source):
            if source.exists(): source.rmdir()
            backup.rename(source)
            record['phase'] = 'verified'; save(journal, record)
        raise
    record['phase'] = 'mounted'; save(journal, record)
    persist(source, target, fstab, record['stamp'])
    record['phase'] = 'ready'; save(journal, record)


def release(source, target, record, journal, *, fstab=None):
    fstab = fstab or FSTAB
    was_releasing = record['phase'] == 'releasing'
    backup = Path(record['backup'])
    if (not re.fullmatch(r'\d{8}T\d{12}Z', record['stamp'])
            or backup != source.with_name(source.name + '.ssd-backup-' + record['stamp'])
            or backup.resolve() != backup or backup.is_symlink() or os.path.ismount(backup)):
        raise RuntimeError('不安全的 SSD 備份路徑，拒絕刪除')
    if (record['phase'] not in {'ready', 'releasing'} or not os.path.ismount(source)
            or not os.path.samefile(source, target)):
        raise RuntimeError('尚未完成掛載與永久設定，不能釋放 SSD')
    original = fstab.read_text()
    if fstab_text(original, source, target) != original:
        raise RuntimeError('永久掛載設定不在，不能釋放 SSD')
    if backup.is_symlink() or (backup.exists() and identity(backup) != record['source_identity']):
        raise RuntimeError('SSD 備份身分不符，拒絕刪除')
    if not backup.exists():
        if record['phase'] != 'releasing': raise RuntimeError('找不到記錄中的 SSD 備份')
    elif record['phase'] == 'ready':
        print('再次逐檔核對 HDD 與 SSD 備份；核對完成前不釋放來源。', flush=True)
        verify_copy(backup, target)
        record['phase'] = 'releasing'; save(journal, record)
    # Keep a metadata snapshot of the final destination even when resuming an
    # interrupted release; a resumed analysis must not invalidate recovery.
    if tree_signature(target) != record['target_tree']:
        raise RuntimeError('HDD 內容在掛載後已變更；先不要接續分析，未繼續清除 SSD')
    if backup.exists() and was_releasing:
        verify_remaining(backup, target)
    if backup.exists(): shutil.rmtree(backup)
    record['phase'] = 'released'; save(journal, record)


def tree_signature(root):
    """Metadata only, no symlink traversal. Full contents checked by rsync."""
    h = hashlib.sha256()
    def fail(error): raise error
    for parent, dirs, files in os.walk(root, followlinks=False, onerror=fail):
        for name in sorted(dirs + files):
            path = Path(parent) / name; s = path.lstat()
            row = [str(path.relative_to(root)), s.st_mode, s.st_uid, s.st_gid,
                   s.st_size, s.st_mtime_ns, os.readlink(path) if path.is_symlink() else None]
            h.update(json.dumps(row, separators=(',', ':')).encode())
        dirs.sort()
    return h.hexdigest()


def capacity(source, target, mounts):
    ancestor = target.parent
    while not ancestor.exists(): ancestor = ancestor.parent
    disk_mounts = [r for r in mounts if ancestor.is_relative_to(Path(r['target']))]
    kind = max(disk_mounts, key=lambda r: len(Path(r['target']).parts))['fstype']
    if (kind not in {'ext4', 'xfs', 'btrfs'} or ancestor.stat().st_dev != Path('/home').stat().st_dev
            or ancestor.stat().st_dev == source.stat().st_dev):
        raise ValueError('目的地不是獨立的本機 /home 資料碟，未搬移')
    used = int(run(['du', '-sx', '-B1', source], capture=True).split()[0])
    free = shutil.disk_usage(ancestor).free
    print(f'study={used / 1024**3:.1f} GiB；HDD可用={free / 1024**3:.1f} GiB', flush=True)
    copied = int(run(['du', '-sx', '-B1', target], capture=True).split()[0]) if target.exists() else 0
    if free < max(0, used - copied) + max(100 * 1024**3, used // 10):
        raise RuntimeError('HDD 空間不足以複製並保留至少100 GiB／10%餘裕')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work-root', type=Path, default=Path('/home/n102968/NGS_UI/nf_work'))
    p.add_argument('--study', type=Path, default=Path('/home/n102968/NGS_UI/nf_work/prs/_validation_20261114/cohort_v1'))
    p.add_argument('--target', type=Path, default=Path('/home/n102968/NGS_UI/prs_validation_work/_validation_20261114/cohort_v1'))
    actions = p.add_mutually_exclusive_group()
    actions.add_argument('--apply', action='store_true')
    actions.add_argument('--release-ssd', action='store_true')
    args = p.parse_args(argv)
    if sys.platform != 'linux' or os.geteuid() != 0: raise RuntimeError('請在 Linux 伺服器以 sudo python3 執行')
    os.umask(0o002)  # Match group-writable existing work guard locks.
    work, source, target = args.work_root, args.study, args.target
    shape(work, source, target)
    journal = target.parent / ('.' + target.name + '.migration.json')
    if journal.is_symlink(): raise ValueError('搬移紀錄不可為 symlink')
    record = read(journal) if journal.exists() else None
    if record and (record.get('schema') != SCHEMA or record['source'] != str(source)
                   or record['target'] != str(target)):
        raise ValueError('搬移紀錄不符，未操作')
    if record and (not target.is_dir() or identity(target) != record['target_identity']):
        raise RuntimeError('HDD 副本目錄身分改變，停止；不覆寫其他資料')
    plan = check_plan(source, source, record['plan_sha256'] if record else None)
    mounts = mount_rows(); names = aliases(work, source, mounts)
    no_nested_mounts(names, mounts); process_users(names | {target})
    print(f'原 study 路徑保持：{source}\nHDD 實體儲存：{target}', flush=True)
    if record and record['phase'] in {'ready', 'mounted', 'released', 'releasing'}:
        with lockset(work, source, target, journal) as hold:
            if not target.is_dir() or not os.path.ismount(source) or not os.path.samefile(source, target):
                raise RuntimeError('已記錄的 study 掛載不在，停止；不回落 SSD 或重新 prepare')
            if record['phase'] == 'mounted' and args.apply:
                persist(source, target, FSTAB, record['stamp'])
                record['phase'] = 'ready'; save(journal, record)
            if args.release_ssd and record['phase'] != 'released':
                backup = Path(record['backup'])
                backup_names = aliases(work, backup, mounts) if backup.exists() else {backup}
                no_nested_mounts(backup_names, mounts); process_users(backup_names)
                if backup.exists(): hold(guard.open_lock(backup / 'execution.lock', blocking=False))
                release(source, target, record, journal)
            print('搬移狀態：' + record['phase'], flush=True)
            run(['df', '-hT', '/', source, target])
        return 0
    if args.release_ssd: raise RuntimeError('尚未完成 --apply，不能釋放 SSD')
    if record and record['phase'] not in {'copying', 'verified'}:
        raise RuntimeError('前次在切換期間中斷；保留兩份資料，先回傳搬移紀錄，不自行改 plan')
    if os.path.ismount(source): raise ValueError('study 已是掛載點，但無相符完成紀錄')
    if not record and target.exists(): raise ValueError('目的地已存在，沒有紀錄，未覆寫')
    capacity(source, target, mounts)
    fstab_text(FSTAB.read_text(), source, target)
    if not args.apply:
        print('檢查通過；未搬移、未修改 study。--apply 後仍保留 SSD 備份。')
        return 0
    with lockset(work, source, target, journal, create=True) as hold:
        process_users(names | {target})
        if record and identity(source) != record['source_identity']: raise RuntimeError('原 study 身分改變')
        if not record:
            owner = (source.stat().st_uid, source.stat().st_gid)
            mkdir_owned(target.parent, owner)
            target.mkdir(mode=0o700); os.chown(target, *owner)
            stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
            record = dict(schema=SCHEMA, source=str(source), target=str(target), phase='copying', stamp=stamp,
                          backup=str(source.with_name(source.name + '.ssd-backup-' + stamp)),
                          source_identity=identity(source), target_identity=identity(target),
                          plan_sha256=digest(source / 'plan.json'))
            save(journal, record)
        # Locks in BOTH source and destination stay held across renaming/bind.
        run(rsync(source, target))
        hold(guard.open_lock(target / 'execution.lock', blocking=False))
        print('複製完成；開始逐檔 checksum、權限與連結核對。', flush=True)
        verify_copy(source, target)
        check_plan(target, source, record['plan_sha256'])
        process_users(names | {target})
        record.update(phase='verified', target_tree=tree_signature(target)); save(journal, record)
        install(source, target, record, journal)
        print('HDD 掛載與永久設定完成；SSD 原檔仍保留。接著執行 --release-ssd，完成前不要恢復分析。', flush=True)
        run(['df', '-hT', '/', source, target])
    return 0


if __name__ == '__main__':
    def stop(*_): raise KeyboardInterrupt()
    signal.signal(signal.SIGTERM, stop)
    try: sys.exit(main())
    except KeyboardInterrupt:
        print('已停止搬移；保留來源與已複製資料。不要自行刪除或重新 prepare。', file=sys.stderr)
        sys.exit(130)
    except BlockingIOError:
        print('PRS 或此 study 仍被其他工作鎖定；請等工作正常結束後重試，不刪鎖檔。', file=sys.stderr)
        sys.exit(1)
    except (OSError, ValueError, RuntimeError, KeyError) as e:
        print('停止搬移：' + str(e), file=sys.stderr); sys.exit(1)
