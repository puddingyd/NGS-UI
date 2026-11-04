"""Container temp routing: portable probes plus optional real Nextflow checks.

Set NGS_UI_TEST_NEXTFLOW to a local Nextflow executable for the integration
checks. Their container engine is a stub; they do not claim a real mount/WGS test.
"""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

import pytest

from app import config
from app.services import prs_paths


@pytest.fixture
def runtime(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    legacy = root/'original_prs'
    legacy.mkdir()
    source = legacy/'nextflow.config'
    source.write_text("singularity.runOptions = '--contain'\nprofiles { singularity { singularity.runOptions = '--no-mount home'; singularity.enabled = true } }\n")
    cfg = {k: str(legacy/k) for k in ('raw_script', 'original_scorefile', 'fasta', 'dbsnp',
                                     'plink_img', 'reference_panel', 'method2_lock', 'method2_model_dir')}
    cfg.update(project=str(legacy), raw_work=str(legacy),
               raw_nextflow_config=str(source), ancestry_nextflow_config=str(source))
    monkeypatch.setattr(config, 'PRS_WORK_ROOT', root/'nf_work/prs')
    effective, paths, wrappers = prs_paths.execution_config('CASE', cfg)
    prs_paths.prepare(paths, wrappers)
    return effective, paths


def test_direct_container_calls_preserve_admin_binds_and_redirect_all_temp(runtime, monkeypatch):
    _, paths = runtime
    monkeypatch.setenv('APPTAINER_BIND', '/data:/data:ro')
    monkeypatch.setenv('SINGULARITY_BINDPATH', '/reference:/reference:ro')
    for prefix in ('APPTAINERENV_', 'SINGULARITYENV_'):
        for key in ('TMPDIR', 'TMP', 'TEMP'):
            monkeypatch.setenv(prefix+key, '/tmp/unwanted')
    _, env = prs_paths.environment(paths, 'run1', 'score')
    tmp = str(prs_paths.container_tmp(paths))
    assert env['APPTAINER_BIND'] == env['SINGULARITY_BIND'] == ''
    for key in ('APPTAINER_BINDPATH', 'SINGULARITY_BINDPATH'):
        assert '/data:/data:ro' in env[key] and '/reference:/reference:ro' in env[key]
        assert env[key].endswith(f'{tmp}:{tmp},{tmp}:/tmp,{tmp}:/var/tmp')
    for prefix in ('APPTAINERENV_', 'SINGULARITYENV_'):
        assert all(env[prefix+k] == tmp for k in ('TMPDIR', 'TMP', 'TEMP'))


@pytest.mark.parametrize('problem', [None, 'missing_mount', 'wrong_temp_env'])
def test_probe_checks_actual_marker_and_write_in_each_mount(runtime, tmp_path, problem):
    cfg, paths = runtime
    directory = prs_paths.container_tmp(paths)
    directory.mkdir(parents=True)
    marker = directory/'.test-marker'
    marker.write_text('unique-token')
    aliases = [tmp_path/'container_tmp', tmp_path/'container_var_tmp']
    for alias in aliases: alias.symlink_to(directory, target_is_directory=True)
    if problem == 'missing_mount':
        aliases[0].unlink(); aliases[0].mkdir()
    args = prs_paths.container_tmp_probe(cfg, paths, marker, 'unique-token')
    assert args[:3] == ['apptainer', 'exec', '--no-home']
    for mount in (f'{directory}:{directory}', f'{directory}:/tmp', f'{directory}:/var/tmp'):
        assert mount in args
    # Host-side model of the two mount aliases; do not write to real /tmp.
    script = args[-1].replace(' /tmp /var/tmp;', ' '+shlex.join(map(str, aliases))+';')
    result = subprocess.run(['/bin/sh', '-c', script], capture_output=True, text=True,
                            env={**os.environ, 'TMPDIR': '/tmp' if problem == 'wrong_temp_env' else str(directory)})
    assert (result.returncode == 0) == (problem is None)
    if problem is None: assert 'CONTAINER_TMP_PASS' in result.stdout
    assert not list(directory.glob('*.write'))


@pytest.mark.parametrize('suffix', ['bad,path', 'bad:path', 'bad\npath'])
def test_invalid_bind_path_rejected_before_job_can_run(runtime, monkeypatch, suffix):
    cfg, paths = runtime
    monkeypatch.setattr(config, 'PRS_WORK_ROOT', Path(paths['work_root']).parent/suffix)
    # Recover the original read-only resource config from the source include.
    original = Path(cfg['raw_script']).parent
    source_cfg = {**cfg, 'project': str(original), 'raw_work': str(original),
                  'raw_nextflow_config': str(original/'nextflow.config'),
                  'ancestry_nextflow_config': str(original/'nextflow.config')}
    with pytest.raises(ValueError, match='容器暫存路徑'):
        prs_paths.execution_config('CASE', source_cfg)


def nextflow_command():
    path = os.environ.get('NGS_UI_TEST_NEXTFLOW')
    if not path:
        pytest.skip('Set NGS_UI_TEST_NEXTFLOW for real Nextflow configuration/launcher tests')
    return str(Path(path).resolve())


@pytest.mark.parametrize('profile', [False, True])
@pytest.mark.parametrize('empty_source', [False, True])
def test_real_nextflow_config_keeps_original_options_and_adds_temp_once(runtime, tmp_path, profile, empty_source):
    nf = nextflow_command()
    cfg, paths = runtime
    if empty_source:
        Path(cfg['raw_script']).with_name('nextflow.config').write_text('')
    args = [nf, '-c', cfg['raw_nextflow_config'], 'config', '-flat']
    if profile: args += ['-profile', 'singularity']
    result = subprocess.run(args, cwd=tmp_path, capture_output=True, text=True,
        env={**os.environ, 'NXF_HOME': str(tmp_path/'nxf_home'), 'NXF_OFFLINE': 'true'}, timeout=30)
    assert result.returncode == 0, result.stderr
    line = next(x for x in result.stdout.splitlines() if x.startswith('singularity.runOptions ='))
    if not empty_source: assert ('--no-mount home' if profile else '--contain') in line
    tmp = str(prs_paths.container_tmp(paths))
    assert line.count(f'{tmp}:/tmp') == line.count(f'{tmp}:/var/tmp') == 1


def test_real_nextflow_launcher_survives_env_cleaning_and_python_uses_case_tmp(runtime, tmp_path):
    nf = nextflow_command()
    cfg, paths = runtime
    project = tmp_path/'mini_pipeline'
    project.mkdir()
    image = tmp_path/'fake.sif'; image.touch()
    # Exercise a real Nextflow-generated `env - ... singularity exec ...`
    # command, while emulating the engine's --env support without Linux mounts.
    bin_dir = tmp_path/'bin'; bin_dir.mkdir()
    engine = bin_dir/'singularity'
    engine.write_text(f'#!{sys.executable}\n'+
        'import os,sys\na=sys.argv[1:]\n'
        'if a==["--version"]:\n print("apptainer version 1.3.0"); sys.exit(0)\n'
        'assert a.pop(0)=="exec", a\n'
        'while a and a[0].startswith("-"):\n'
        ' k=a.pop(0)\n'
        ' if k=="--env":\n'
        '  for item in a.pop(0).split(","):\n'
        '   key,value=item.split("=",1); os.environ[key]=value\n'
        ' elif k in ("-B","--bind","--no-mount"): a.pop(0)\n'
        ' elif k not in ("--no-home","--pid","--contain"): raise ValueError(k)\n'
        'a.pop(0)\nos.execvpe(a[0],a,os.environ)\n')
    engine.chmod(0o755)
    project.joinpath('main.nf').write_text('''process PROBE {
    output: path 'temp.json'
    script:
    """
    python3 - <<'PY'
import json, tempfile
from pathlib import Path
with tempfile.NamedTemporaryFile() as f:
    Path('temp.json').write_text(json.dumps({'temporary_file':f.name}))
PY
    """
}
workflow { PROBE() }
''')
    project.joinpath('nextflow.config').write_text(f"process.container = '{image}'\n")
    _, env = prs_paths.environment(paths, 'integration', 'score')
    env.update(NXF_OFFLINE='true', PATH=str(bin_dir)+os.pathsep+str(Path(sys.executable).parent)+os.pathsep+env['PATH'])
    run = subprocess.run([nf, '-c', cfg['raw_nextflow_config'], 'run', str(project),
        '-profile', 'singularity', '-w', str(tmp_path/'work')], cwd=tmp_path,
        capture_output=True, text=True, env=env, timeout=30)
    assert run.returncode == 0, run.stdout+'\n'+run.stderr
    result = next((tmp_path/'work').glob('*/*/temp.json'))
    assert Path(json.loads(result.read_text())['temporary_file']).parent == prs_paths.container_tmp(paths)
    launcher = result.with_name('.command.run').read_text()
    assert 'env - PATH=' in launcher
    assert f'{prs_paths.container_tmp(paths)}:/tmp' in launcher
