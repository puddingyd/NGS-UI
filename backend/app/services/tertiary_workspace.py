"""Writable tertiary runtime and verified transfer to the publication disk.

Keep the logical work root stable across an SSD bind mount. This module is
stdlib-only so it can also restore runtime settings *after* the site's shell
environment has been sourced, before exec'ing Nextflow.
"""
from __future__ import annotations

import fcntl
import hashlib
import os
from pathlib import Path
import re
import shlex
import shutil
import stat
import sys
import tempfile

from . import nextflow_work


def runtime_paths(root: Path, scope: str, job_id: str) -> dict[str, Path]:
    if scope not in {'dragen', 'nckuh'} or not re.fullmatch(r'[A-Za-z0-9_-]+', job_id):
        raise ValueError('Invalid tertiary runtime scope/job')
    root = nextflow_work.work_root(root)
    if any(c in str(root) for c in (':', ',', '\n', '\r', '\0')):
        raise ValueError('三級工作路徑不能包含冒號、逗號或換行')
    runtime = root / 'contexts' / f'shared-{scope}' / '.runtime'
    job = runtime / 'jobs' / job_id
    paths = dict(root=root, runtime=runtime, job=job, tmp=job/'tmp', stage=job/'staging',
                 nextflow=runtime/'nextflow', containers=runtime/'containers',
                 apptainer=runtime/'apptainer', conda=runtime/'conda',
                 cache=runtime.parent/'launch/.nextflow/cache')
    for path in paths.values():
        if not path.resolve().is_relative_to(root):
            raise ValueError(f'三級工作子目錄不可連結至工作區外：{path}')
    return paths


def _seed_cache(source: Path, target: Path) -> None:
    """Preserve installed framework/plugins/assets without a fresh download."""
    if target.exists():
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    if not source.is_dir() or source.resolve() == target.resolve():
        return
    temporary = Path(tempfile.mkdtemp(prefix='.seed-', dir=target.parent))
    try:
        shutil.copytree(source, temporary, symlinks=True, dirs_exist_ok=True)
        temporary.replace(target)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def environment(root: Path, scope: str, job_id: str) -> dict[str, str]:
    paths = runtime_paths(root, scope, job_id)
    for name in ('runtime', 'job', 'tmp', 'nextflow', 'containers', 'apptainer', 'conda'):
        paths[name].mkdir(parents=True, exist_ok=True)
    previous_home = Path(os.environ.get('NXF_HOME', Path.home()/'.nextflow'))
    # Different jobs can initialize the same mode before taking its cache lock.
    with (paths['runtime']/'.seed.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        for source, target in (
            (previous_home/'framework', paths['nextflow']/'framework'),
            (Path(os.environ.get('NXF_ASSETS', previous_home/'assets')), paths['nextflow']/'assets'),
            (Path(os.environ.get('NXF_PLUGINS_DIR', previous_home/'plugins')), paths['nextflow']/'plugins'),
        ):
            _seed_cache(source, target)
    env = {key: str(paths[name]) for key, name in (
        ('NXF_HOME', 'nextflow'), ('NXF_TEMP', 'tmp'),
        ('NXF_SINGULARITY_CACHEDIR', 'containers'), ('NXF_APPTAINER_CACHEDIR', 'containers'),
        ('NXF_CONDA_CACHEDIR', 'conda'),
        ('APPTAINER_CACHEDIR', 'apptainer'), ('SINGULARITY_CACHEDIR', 'apptainer'),
        ('APPTAINER_TMPDIR', 'tmp'), ('SINGULARITY_TMPDIR', 'tmp'),
        ('TMPDIR', 'tmp'), ('TMP', 'tmp'), ('TEMP', 'tmp'),
    )}
    env.update(NXF_ASSETS=str(paths['nextflow']/'assets'),
               NXF_PLUGINS_DIR=str(paths['nextflow']/'plugins'),
               NXF_CACHE_DIR=str(paths['cache']))
    tmp = str(paths['tmp'])
    names = ('APPTAINER_BIND', 'APPTAINER_BINDPATH', 'SINGULARITY_BIND', 'SINGULARITY_BINDPATH')
    required = [f'{paths["root"]}:{paths["root"]}', f'{tmp}:{tmp}', f'{tmp}:/tmp', f'{tmp}:/var/tmp']
    binds = list(dict.fromkeys(b for name in names for b in os.environ.get(name, '').split(',') if b))
    for name in names:
        env[name] = ','.join([b for b in binds if b not in required] + required)
    for prefix in ('APPTAINERENV_', 'SINGULARITYENV_'):
        for name in ('TMPDIR', 'TMP', 'TEMP'):
            env[prefix+name] = tmp
    return env


def _groovy(value) -> str:
    return "'" + str(value).replace('\\', '\\\\').replace("'", "\\'").replace('\n', '\\n') + "'"


def nextflow_config(source: Path, root: Path, scope: str, job_id: str) -> Path:
    paths = runtime_paths(root, scope, job_id)
    paths['job'].mkdir(parents=True, exist_ok=True)
    tmp = str(paths['tmp'])
    options = shlex.join(['--bind', f'{paths["root"]}:{paths["root"]}',
                          '--bind', f'{tmp}:{tmp}', '--bind', f'{tmp}:/tmp',
                          '--bind', f'{tmp}:/var/tmp',
                          '--env', f'TMPDIR={tmp},TMP={tmp},TEMP={tmp}'])
    literal = _groovy(options)
    overrides = ''
    for engine in ('singularity', 'apptainer'):
        overrides += f'{engine}.cacheDir = {_groovy(paths["containers"])}\n'
        overrides += (f'{engine}.runOptions = (({engine}.runOptions ?: "").contains({literal}) '
                      f'? {engine}.runOptions : (({engine}.runOptions ?: "") + " " + {literal}).trim())\n')
    overrides += f'conda.cacheDir = {_groovy(paths["conda"])}\n'
    overrides += f'process.scratch = {_groovy(paths["tmp"])}\n'
    for key in ('TMPDIR', 'TMP', 'TEMP'):
        overrides += f'env.{key} = {_groovy(tmp)}\n'
    # The production -profile dgm must not override these writable locations.
    text = (f'// NGS-UI: execution paths only; retain source analysis settings.\n'
            f'includeConfig {_groovy(Path(source).resolve())}\n{overrides}'
            f'profiles {{\n  dgm {{\n{overrides}  }}\n}}\n')
    target = paths['job']/'nextflow.config'
    target.write_text(text)
    return target


def nextflow_command(command: list[str], root: Path, scope: str, job_id: str,
                     env_script: str, python: str, backend: Path) -> list[str]:
    # Run the stdlib-only helper after source: site scripts may reset temp/cache
    # variables or add important tool paths and container bindings.
    launcher = ['env', f'PYTHONPATH={backend}', python, '-m',
                'app.services.tertiary_workspace', str(root), scope, job_id, *command]
    shell = (f'if [ -f {shlex.quote(env_script)} ]; then source {shlex.quote(env_script)}; '
             "else echo '[nextflow] env script not found; using current environment'; fi; "
             f'exec {shlex.join(launcher)}')
    return ['bash', '-lc', shell]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def verified_copytree(source: Path, target: Path) -> dict[str, int]:
    """Materialize published files; verify bytes before any live replacement.

    Filesystem times survive transfer for raw/overlay review signatures. No
    symlink in published results may depend on disposable SSD staging/work.
    """
    if target.exists() or target.is_symlink():
        raise RuntimeError(f'publication staging already exists: {target}')
    counts = dict(files=0, bytes=0)

    def copy(directory: Path, destination: Path, ancestors: set[Path]) -> None:
        resolved = directory.resolve(strict=True)
        if resolved in ancestors:
            raise RuntimeError(f'cyclic staging directory: {directory}')
        destination.mkdir()
        for item in sorted(directory.iterdir()):
            output = destination/item.name
            if item.is_dir():
                copy(item, output, ancestors | {resolved})
                continue
            before = item.stat()
            if not stat.S_ISREG(before.st_mode):
                raise RuntimeError(f'unsupported staging file: {item}')
            shutil.copy2(item, output)
            after = item.stat()
            if ((before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns)
                    or _sha256(item) != _sha256(output)):
                raise RuntimeError(f'publication checksum mismatch: {item}')
            with output.open('rb') as handle:
                os.fsync(handle.fileno())
            counts['files'] += 1
            counts['bytes'] += before.st_size
        shutil.copystat(directory, destination)

    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        copy(source, target, set())
    except BaseException:
        shutil.rmtree(target, ignore_errors=True)
        raise
    return counts


if __name__ == '__main__':
    work_root, mode, jid, *argv = sys.argv[1:]
    if not argv:
        raise SystemExit('Nextflow command required')
    os.environ.update(environment(Path(work_root), mode, jid))
    # Preserve launcher/site PATH and inherited lease descriptors.
    os.execvpe(argv[0], argv, os.environ)
