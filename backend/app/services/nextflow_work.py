"""Nextflow work leases and conservative cleanup (stdlib only).

This module is also embedded verbatim in DGX2 commands: keep it independent
of app imports. Locks live beside the work root and are never removed by
cleanup. Closing a lease, rather than LOCK_UN, keeps inherited child leases
valid if a launcher exits unexpectedly.
"""
from __future__ import annotations

import fcntl
import hashlib
import os
import re
import shutil
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path

ROOT_ENV = "NGS_UI_NF_GUARD_ROOT"
SCOPE_ENV = "NGS_UI_NF_GUARD_SCOPE"
FD_ENV = "NGS_UI_NF_GUARD_FD"


def work_root(value) -> Path:
    path = Path(value)
    if not path.is_absolute() or path.is_symlink():
        raise ValueError(f"不安全的 Nextflow work root：{path}")
    path = path.resolve()
    if path == path.parent or path in {Path.home(), Path('/tmp'), Path('/private/tmp')}:
        raise ValueError(f"不安全的 Nextflow work root：{path}")
    return path


def _lock_path(root: Path, scope: str) -> Path:
    directory = root.parent / f".{root.name}.ngs-ui-locks"
    try:
        directory.mkdir(mode=0o2775, parents=True)
    except FileExistsError:
        if directory.is_symlink() or not directory.is_dir():
            raise ValueError(f"不安全的 Nextflow 鎖目錄：{directory}")
    else:
        if root.exists() and directory.stat().st_gid != root.stat().st_gid:
            os.chown(directory, -1, root.stat().st_gid)
        directory.chmod(0o2775)
    name = hashlib.sha256(scope.encode()).hexdigest() if scope else "admission"
    return directory / f"{name}.lock"


def open_lock(path: Path, *, shared=False, blocking=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    # Group-writable DGX2 batches may be managed by different users.
    fd = os.open(path, os.O_RDWR | os.O_CREAT | getattr(os, 'O_NOFOLLOW', 0), 0o666)
    handle = os.fdopen(fd, 'a+')
    try:
        fcntl.flock(fd, (fcntl.LOCK_SH if shared else fcntl.LOCK_EX)
                    | (0 if blocking else fcntl.LOCK_NB))
    except BaseException:
        handle.close()
        raise
    return handle


@contextmanager
def admission(root):
    handle = open_lock(_lock_path(work_root(root), ""))
    try:
        yield
    finally:
        handle.close()


def lease_env(root, scope, handle):
    return {ROOT_ENV: str(work_root(root)), SCOPE_ENV: scope, FD_ENV: str(handle.fileno())}


def acquire_lease(root, scope):
    root = work_root(root)
    if not scope:
        raise ValueError("Nextflow work scope is required")
    # Adoption is allowed only for the exact inherited lock inode.
    if os.environ.get(ROOT_ENV) == str(root) and os.environ.get(SCOPE_ENV) == scope:
        try:
            fd = int(os.environ[FD_ENV])
            actual, expected = os.fstat(fd), _lock_path(root, scope).stat()
            if fd > 2 and (actual.st_dev, actual.st_ino) == (expected.st_dev, expected.st_ino):
                return os.fdopen(os.dup(fd), 'a+')
        except (OSError, ValueError, KeyError):
            pass
    with admission(root):
        return open_lock(_lock_path(root, scope), shared=True)


def inherited_fds():
    try:
        fd = int(os.environ[FD_ENV])
        os.fstat(fd)
        return (fd,) if fd > 2 else ()
    except (OSError, ValueError, KeyError):
        return ()


def nextflow_process_scopes(root):
    """Unknown/manual processes fail closed; never infer safety from a PID alone."""
    try:
        result = subprocess.run(['ps', '-e', '-ww', '-o', 'pid=,comm=,args='],
                                capture_output=True, text=True, check=True)
    except (OSError, subprocess.SubprocessError) as error:
        raise RuntimeError('無法確認 Nextflow 程序，未執行清理') from error
    scopes, unknown = set(), []
    for line in result.stdout.splitlines():
        fields = line.split(None, 2)
        if len(fields) != 3:
            continue
        pid, comm, args = fields
        executable = Path(comm).name.lower()
        is_nextflow = executable == 'nextflow' or (
            executable == 'java' and re.search(r'nextflow(?:\.cli\.Launcher|[^\s]*\.jar)', args)
        ) or (executable in {'bash', 'sh'} and (
            re.match(r'\S+\s+\S*/?nextflow(?:\s|$)', args)
            or re.search(r'(?:^|[/\s])\.command\.(?:run|sh)(?:\s|$)', args)))
        if not is_nextflow:
            continue
        try:
            data = (Path('/proc') / pid / 'environ').read_bytes()
            env = dict(item.split(b'=', 1) for item in data.split(b'\0') if b'=' in item)
            scope = env.get(SCOPE_ENV.encode(), b'').decode()
            if env.get(ROOT_ENV.encode(), b'').decode() != str(root) or not scope:
                unknown.append(pid)
                continue
            # An environment marker without a live lease is not protection.
            try:
                probe = open_lock(_lock_path(root, scope), blocking=False)
            except BlockingIOError:
                scopes.add(scope)
            else:
                probe.close()
                unknown.append(pid)
        except FileNotFoundError:
            # A vanished process is safe; an unreadable /proc (e.g. macOS or
            # hidepid) is not evidence that a matching process has exited.
            try:
                os.kill(int(pid), 0)
            except ProcessLookupError:
                continue
            except PermissionError:
                pass
            unknown.append(pid)
        except (OSError, UnicodeError, ValueError):
            unknown.append(pid)
    return scopes, unknown


def _overlap(a, b):
    return a == b or a in b.parents or b in a.parents


def _references(root, directories):
    refs = set()
    def add_link(path):
        # Resolve one link at a time, including links in parent directories.
        # An alias outside work can lead through work/A to work/B; both are
        # dependencies, even though Path.resolve() would return only B.
        target = Path(os.path.abspath(path))
        for _ in range(64):
            if target == root or root in target.parents:
                refs.add(target)
            for prefix in [*reversed(target.parents), target]:
                if prefix.is_symlink():
                    suffix = target.relative_to(prefix)
                    target = Path(os.path.abspath(prefix.parent / os.readlink(prefix) / suffix))
                    break
            else:
                return
        raise RuntimeError('too many symbolic links')
    for directory in directories:
        directory = Path(directory)
        if directory.is_symlink():
            add_link(directory)
        if not directory.exists():
            continue
        if directory.resolve() == root or root in directory.resolve().parents:
            refs.add(directory.resolve())
        def fail(error):
            raise error
        for parent, folders, files in os.walk(directory, followlinks=False, onerror=fail):
            for name in folders + files:
                path = Path(parent) / name
                if path.is_symlink():
                    add_link(path)
    return refs


def _units(root, tertiary):
    units = []
    for child in sorted(root.iterdir()):
        if tertiary and child.name == 'contexts' and child.is_dir() and not child.is_symlink():
            units.extend(sorted(child.iterdir()))
        else:
            units.append(child)
    return units


def cleanup(root, *, tertiary=False, busy_scopes=(), output_roots=(), preserve_paths=()):
    root = work_root(root)
    result = dict(path=str(root), deleted_count=0, freed_bytes=0,
                  skipped_count=0, failed_count=0, deleted=[], skipped=[], failed=[],
                  existed=root.exists())
    with admission(root):
        if not root.exists():
            return result
        if not root.is_dir():
            raise ValueError(f"Nextflow work root 不是資料夾：{root}")
        units = _units(root, tertiary)
        handles = []
        try:
            scopes = {'dragen', 'nckuh'} if tertiary else {p.name for p in units}
            busy = set(busy_scopes)
            for scope in sorted(scopes):
                try:
                    handles.append(open_lock(_lock_path(root, scope), blocking=False))
                except BlockingIOError:
                    busy.add(scope)
            process_scopes, unknown = nextflow_process_scopes(root)
            busy.update(process_scopes)
            reason_all = f"無法確認 Nextflow 程序歸屬（PID {', '.join(unknown)}）" if unknown else ''
            if tertiary and any(not (root / 'contexts' / f'shared-{scope}' / 'work-dir.txt').is_file()
                                for scope in busy & {'dragen', 'nckuh'}):
                reason_all = '使用中的共用暫存尚在初始化，保留可能承接的舊暫存'
            protected = {}
            # PRS raw audits/recovery inputs are retained even after completion.
            # Reserve the default subtree across config changes; also protect an
            # administrator's custom PRS root when it lies within this work tree.
            retained = [Path(p).resolve() for p in preserve_paths]
            if tertiary:
                retained.append(root / 'prs')
            # Pointer dependencies are undirected for deletion: keep both the
            # lineage and its adopted work, including any ancestor unit.
            edges = []
            for path in units:
                scope = path.name.removeprefix('shared-') if tertiary else path.name
                if any(_overlap(path, p) for p in retained):
                    protected[path] = 'PRS 重建／稽核與接續分析資料，保留不作一般暫存清理'
                if scope in busy:
                    protected.setdefault(path, '執行中或等待中的工作正在使用此暫存')
                if path.is_symlink():
                    protected[path] = '暫存入口為符號連結，保留待確認'
                if tertiary and path.name in {'shared-dragen', 'shared-nckuh'}:
                    pointer = path / 'work-dir.txt'
                    if not pointer.exists():
                        if scope in busy:
                            reason_all = '使用中的共用暫存尚在初始化，保留可能承接的舊暫存'
                        continue
                    try:
                        value = pointer.read_text().strip()
                    except (OSError, UnicodeError):
                        reason_all = '共用暫存的 work-dir.txt 無法安全解析'
                        continue
                    target = Path(value)
                    if not value or not target.is_absolute():
                        reason_all = '共用暫存的 work-dir.txt 無法安全解析'
                        continue
                    target = target.resolve()
                    if target != root and root not in target.parents:
                        reason_all = '共用暫存指向 work root 外，保留待確認'
                        continue
                    edges.extend((path, other) for other in units if _overlap(other, target))
            def propagate():
                changed = True
                while changed:
                    changed = False
                    for a, b in edges:
                        if (a in protected) != (b in protected):
                            protected[a] = protected[b] = '與受保護的共用暫存相依'
                            changed = True
            propagate()
            try:
                # Published symlinks and inputs staged from another work tree
                # remain references even after the producing run has ended.
                pending = list(output_roots) + list(protected)
                scanned = set()
                while pending:
                    scanned.update(pending)
                    refs = _references(root, pending)
                    for path in units:
                        if any(_overlap(path, ref) for ref in refs):
                            protected.setdefault(path, '正式輸出或使用中的工作仍以連結引用此暫存')
                    propagate()
                    pending = [path for path in protected if path not in scanned]
            except (OSError, RuntimeError):
                reason_all = '無法完整檢查輸出／暫存連結，保留待確認'
            propagate()
            for path in units:
                reason = reason_all or protected.get(path)
                if path.is_symlink():
                    reason = reason or '暫存入口為符號連結，保留待確認'
                if reason:
                    result['skipped'].append(dict(path=str(path), reason=reason))
                    continue
                try:
                    size = 0
                    def fail(error):
                        raise error
                    if path.is_dir():
                        for parent, _, files in os.walk(path, followlinks=False, onerror=fail):
                            for name in files:
                                size += (Path(parent) / name).lstat().st_size
                        shutil.rmtree(path)
                    else:
                        size = path.lstat().st_size
                        path.unlink()
                    result['deleted'].append(dict(path=str(path)))
                    result['freed_bytes'] += size
                except OSError as exc:
                    result['failed'].append(dict(path=str(path), reason=str(exc)))
            for key in ('deleted', 'skipped', 'failed'):
                result[f'{key}_count'] = len(result[key])
            return result
        finally:
            for handle in handles:
                handle.close()


def main():
    os.umask(0o002)
    action, root, *args = sys.argv[1:]
    if action == 'run':
        scope, *command = args
        handle = acquire_lease(root, scope)
        try:
            env = {**os.environ, **lease_env(root, scope, handle)}
            return subprocess.call(command, env=env, pass_fds=(handle.fileno(),))
        finally:
            handle.close()
    if action != 'cleanup':
        raise ValueError('unknown action')
    print(f'只清理閒置且未被引用的批次：{root}；使用中的批次會保留。', flush=True)
    for path in sorted(work_root(root).iterdir()):
        print(f'候選批次：{path}')
    if input('確定清理以上二級分析 Nextflow 暫存？[y/N] ').lower() not in {'y', 'yes'}:
        print('已取消。')
        return 0
    result = cleanup(root, output_roots=args)
    for key, label in (('deleted', '已清理'), ('skipped', '已保留'), ('failed', '清理失敗')):
        for item in result[key]:
            print(f"{label}：{item['path']} {item.get('reason', '')}")
    print(f"清理結果：已清理 {result['deleted_count']}；保留 {result['skipped_count']}；失敗 {result['failed_count']}")
    return 1 if result['failed_count'] else 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f'停止清理／啟動：{error}', file=sys.stderr)
        sys.exit(1)
