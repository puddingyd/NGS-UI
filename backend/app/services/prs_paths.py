"""Platform-owned PRS execution paths; scientific resources remain read-only."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import shlex
import tempfile

from .. import config
from . import prs, sample_layout, container_runtime, nextflow_work

VERSION = 'NGS_UI_PRS_PATHS_V1'
SCORE_FILES = ('summary.json', 'scores.tsv', 'report.md', 'status.json',
               'inputs.lock.json', 'missing_loci.tsv.gz')
QC_FILES = ('qc.json', 'qc.tsv', 'report.md', 'missing_high_weight.tsv')


def digest_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def layout(sid, cfg):
    prs.validate_sid(sid)
    root = config.PRS_WORK_ROOT
    if not root.is_absolute():
        raise ValueError('NGS_UI_PRS_WORK_ROOT 必須是絕對路徑')
    root = root.resolve()
    project = root / sid
    # Refuse aliases into another case or a resource/research tree.
    if project.is_symlink() or project.resolve().parent != root:
        raise ValueError('PRS 個案工作路徑不可連結至其他目錄')
    protected = [Path(cfg[k]).resolve() for k in ('project', 'raw_work', 'method2_model_dir', 'multi_model_dir') if cfg.get(k)]
    resources = [Path(cfg[k]).resolve() for k in (
        'raw_script', 'original_scorefile', 'fasta', 'dbsnp', 'plink_img',
        'raw_nextflow_config', 'ancestry_nextflow_config', 'reference_panel', 'method2_lock')]
    if (any(root.is_relative_to(p) or p.is_relative_to(root) for p in protected)
            or any(p.is_relative_to(root) for p in resources)):
        raise ValueError('PRS 工作目錄不得位於原 PRS 專案／參考資源內，或包含這些資源')
    return dict(version=VERSION, work_root=str(root), project=str(project))


def groovy_string(value):
    # Single-quoted Groovy literals do not interpolate dollar signs.
    return "'" + str(value).replace('\\', '\\\\').replace("'", "\\'").replace('\n', '\\n').replace('\r', '\\r') + "'"


def container_tmp(paths):
    path = Path(paths['project']) / '.runtime/tmp'
    # Apptainer bind specifications use colon/comma separators, not shell syntax.
    if any(c in str(path) for c in (':', ',', '\n', '\r', '\0')):
        raise ValueError('PRS 容器暫存路徑不能包含冒號、逗號或換行')
    return path


def container_tmp_options(paths):
    path = str(container_tmp(paths))
    # Mount the same absolute path for Nextflow's forwarded TMPDIR, plus both
    # conventional temporary directories for software that ignores TMPDIR.
    return ['--bind', f'{path}:{path}', '--bind', f'{path}:/tmp',
            '--bind', f'{path}:/var/tmp',
            '--env', f'TMPDIR={path},TMP={path},TEMP={path}']


def execution_config(sid, cfg):
    return _execution_config(cfg, layout(sid, cfg))


def validation_config(workspace, cfg):
    """Use the same runtime rules inside an isolated validation namespace."""
    root = config.PRS_WORK_ROOT.resolve()
    workspace = Path(workspace).absolute()
    try:
        # The configured work root may itself be the existing SSD alias.
        # Accept that root alias, but never aliases below the research namespace.
        parts = workspace.relative_to(config.PRS_WORK_ROOT.absolute()).parts
    except ValueError:
        try:
            parts = workspace.relative_to(root).parts
        except ValueError:
            raise ValueError('驗證工作區必須位於既有 PRS work root')
    workspace = root.joinpath(*parts)
    if len(parts) < 3 or not parts[0].startswith('_validation_'):
        raise ValueError('驗證工作區需位於 _validation_*/study/...，不能覆寫平台個案')
    prs.validate_sid(parts[0].lstrip('_'))
    for part in parts[1:]:
        prs.validate_sid(part)
    # The leading underscore identifies research data rather than a UI sample.
    paths = layout('validation-'+parts[0].lstrip('_'), cfg)
    if workspace.resolve() != workspace:
        raise ValueError('驗證工作區不能經由 symlink 指到其他目錄')
    paths['project'] = str(workspace)
    return _execution_config(cfg, paths)


def _execution_config(cfg, paths):
    project = Path(paths['project'])
    runtime = project / '.runtime'
    effective = {**cfg, 'project': str(project), 'raw_work': str(project)}
    tmp_options = shlex.join(container_tmp_options(paths))
    wrappers = {}
    for key in ('raw_nextflow_config', 'ancestry_nextflow_config'):
        # Keep source include locations intact, then override execution-only
        # destinations, including values inside the selected singularity profile.
        overrides = '\n'.join(f'{k} = {groovy_string(v)}' for k, v in (
            ('singularity.cacheDir', runtime/'containers'),
            ('conda.cacheDir', runtime/'conda'),
            ('report.file', 'report.html'),
            ('timeline.file', 'timeline.html'),
            ('trace.file', 'trace.txt'),
            ('dag.file', 'dag.html')))
        literal = groovy_string(tmp_options)
        container_overrides = ('singularity.runOptions = '
            f"((singularity.runOptions ?: '').contains({literal}) ? singularity.runOptions : "
            f"((singularity.runOptions ?: '') + ' ' + {literal}).trim())\n")
        text = ('// NGS-UI: preserve model/process settings; own writable paths.\n'
                f'includeConfig {groovy_string(Path(cfg[key]).resolve())}\n'
                f'{overrides}\n{container_overrides}'
                f'profiles {{\n  singularity {{\n{overrides}\n{container_overrides}  }}\n}}\n')
        path = runtime/'configs'/f'{key}-{hashlib.sha256(text.encode()).hexdigest()}.config'
        wrappers[path] = text
        effective[key] = str(path)
    return effective, paths, wrappers


def owned_directory(project, relative):
    path = project / relative
    if not path.resolve().is_relative_to(project.resolve()):
        raise ValueError('PRS 工作子目錄不可連結至其他目錄')
    path.mkdir(parents=True, exist_ok=True)
    return path


def prepare(paths, wrappers):
    project = Path(paths['project'])
    project.mkdir(parents=True, exist_ok=True)
    marker = project / '.ngs-ui-prs.json'
    expected = dict(schema=VERSION, sample_id=project.name)
    if marker.exists() and prs.read_json(marker)[0] != expected:
        raise ValueError('PRS 工作目錄不屬於此個案')
    if not marker.exists():
        prs.atomic_json(marker, expected)
    for path, content in wrappers.items():
        owned_directory(project, path.parent.relative_to(project))
        if path.exists():
            if path.is_symlink() or path.read_text() != content:
                raise ValueError('PRS 平台 Nextflow 設定已變更，請檢查工作目錄')
        else:
            path.write_text(content)


def environment(paths, run_id, stage):
    project = Path(paths['project'])
    runtime = project / '.runtime'
    cwd = owned_directory(project, Path('launch') / run_id / stage)
    dirs = {
        'NXF_HOME': runtime/'nextflow', 'NXF_ASSETS': runtime/'assets',
        'NXF_PLUGINS_DIR': runtime/'plugins', 'NXF_TEMP': runtime/'tmp',
        'NXF_SINGULARITY_CACHEDIR': runtime/'containers',
        'NXF_CONDA_CACHEDIR': runtime/'conda',
        'APPTAINER_CACHEDIR': runtime/'apptainer', 'SINGULARITY_CACHEDIR': runtime/'apptainer',
        'APPTAINER_TMPDIR': runtime/'tmp', 'SINGULARITY_TMPDIR': runtime/'tmp',
        'TMPDIR': runtime/'tmp', 'TMP': runtime/'tmp', 'TEMP': runtime/'tmp',
    }
    env = dict(os.environ)
    # Keep the task database/history alongside each actual Nextflow launch.
    env.pop('NXF_CACHE_DIR', None)
    env.pop('NXF_WORK', None)  # Both frozen entry points explicitly supply -w.
    for key, path in dirs.items():
        env[key] = str(owned_directory(project, path.relative_to(project)))
    tmp = str(container_tmp(paths))
    # Direct Apptainer calls in the frozen raw script do not use Nextflow's
    # runOptions. Cover them too, preserving administrator-supplied bind paths.
    required = [f'{tmp}:{tmp}', f'{tmp}:/tmp', f'{tmp}:/var/tmp']
    env.update(container_runtime.job_bind_environment(env, config.TERTIARY_NF_WORK_ROOT, required))
    # This is a new job boundary, not an inherited tertiary lease. The PRS
    # worker acquires and explicitly passes its own lease to each command.
    for name in (nextflow_work.ROOT_ENV, nextflow_work.SCOPE_ENV, nextflow_work.FD_ENV):
        env.pop(name, None)
    for prefix in ('APPTAINERENV_', 'SINGULARITYENV_'):
        for name in ('TMPDIR', 'TMP', 'TEMP'):
            env[prefix+name] = tmp
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    return cwd, env


def container_tmp_probe(cfg, paths, marker, token):
    """Fail before scoring if the container cannot see/write the host temp mount."""
    path = container_tmp(paths)
    script = ('set -eu\n'
        'fail() { printf "[PRS] container temporary directory check failed: %s\\n" "$1" >&2; exit 72; }\n'
        f'test "$TMPDIR" = {shlex.quote(str(path))} || fail "unexpected TMPDIR: $TMPDIR"\n'
        f'for dir in {shlex.quote(str(path))} /tmp /var/tmp; do\n'
        f'  test "$(cat "$dir/{marker.name}" 2>/dev/null)" = {shlex.quote(token)} || fail "mount unavailable: $dir"\n'
        f'  probe="$dir/{marker.name}.write"\n'
        '  (umask 077; printf "PRS tmp probe\\n" > "$probe") || fail "cannot write: $dir"\n'
        '  rm "$probe" || fail "cannot remove probe: $dir"\n'
        'done\n'
        'printf "[PRS] CONTAINER_TMP_PASS: TMPDIR, /tmp, /var/tmp use case workspace\\n"\n')
    return ['apptainer', 'exec', '--no-home', *container_tmp_options(paths),
            cfg['plink_img'], '/bin/sh', '-c', script]


def publish_files(sid, jid, kind, source, names):
    """Copy complete fixed-name artifacts, then atomically publish the directory."""
    base = sample_layout.state_dir(sid, for_write=True) / f'{sid}.prs' / 'runs' / jid
    base.mkdir(parents=True, exist_ok=True)
    target = base / kind
    stage = Path(tempfile.mkdtemp(prefix=f'.{kind}-', dir=base))
    try:
        hashes = {}
        for name in names:
            src = source / name
            if not src.is_file():
                raise ValueError(f'PRS 完整結果缺少 {kind}/{name}')
            shutil.copyfile(src, stage/name)  # Dereference upstream artifacts.
            hashes[name] = file_hash(stage/name)
        manifest = dict(schema=VERSION, job_id=jid, source=str(source), files=hashes)
        prs.atomic_json(stage/'manifest.json', manifest)
        if target.exists():
            if prs.read_json(target/'manifest.json')[0] != manifest or any(
                    not (target/n).is_file() or file_hash(target/n) != h for n, h in hashes.items()):
                raise ValueError('已發布的 PRS 結果不一致')
        else:
            os.replace(stage, target)
        return dict(path=str(target), files=hashes)
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def file_hash(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def workspace_for_deletion(sid):
    """Preflight a marked platform-owned case, never legacy config.project."""
    prs.validate_sid(sid)
    path = config.PRS_WORK_ROOT.resolve() / sid
    if not path.exists():
        return None
    if path.is_symlink() or path.resolve().parent != config.PRS_WORK_ROOT.resolve():
        raise ValueError('拒絕刪除連結至其他目錄的 PRS 工作目錄')
    marker, _ = prs.read_json(path/'.ngs-ui-prs.json')
    if marker != dict(schema=VERSION, sample_id=sid):
        raise ValueError('拒絕刪除未標記的 PRS 工作目錄')
    return path
