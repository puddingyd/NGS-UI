"""Tertiary SSD routing, verified publication and failure containment."""
import errno
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys

import pytest

from app.services import nextflow_work, sample_layout, tertiary_workspace as workspace
from app.services import snv_overlay
from app.workers import dragen_run


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    home = root/'old-nextflow'
    for name in ('framework', 'plugins', 'assets'):
        (home/name).mkdir(parents=True)
        (home/name/'installed').write_text(name)
    monkeypatch.setenv('NXF_HOME', str(home))
    monkeypatch.delenv('NXF_ASSETS', raising=False)
    monkeypatch.delenv('NXF_PLUGINS_DIR', raising=False)
    return root/'nf_work'


def test_runtime_is_under_the_leased_mode_and_preserves_installed_nextflow(isolated, monkeypatch):
    monkeypatch.setenv('APPTAINER_BIND', '/reference:/reference:ro')
    first = workspace.environment(isolated, 'dragen', 'job1')
    second = workspace.environment(isolated, 'dragen', 'job2')
    assert first['TMPDIR'] != second['TMPDIR']
    assert first['NXF_HOME'] == second['NXF_HOME']
    for name in ('framework', 'plugins', 'assets'):
        assert (Path(first['NXF_HOME'])/name/'installed').read_text() == name
    assert first['NXF_CACHE_DIR'] == str(isolated/'contexts/shared-dragen/launch/.nextflow/cache')
    assert '/reference:/reference:ro' in first['APPTAINER_BINDPATH']
    assert first['APPTAINERENV_TMPDIR'] == first['TMPDIR']
    paths = workspace.runtime_paths(isolated, 'dragen', 'job1')
    paths['stage'].mkdir()
    (paths['stage']/'active-data').write_text('keep')
    # Existing cleanup must keep jobs, staging and caches throughout the lease.
    monkeypatch.setattr(nextflow_work, 'nextflow_process_scopes', lambda root: (set(), []))
    with nextflow_work.acquire_lease(isolated, 'dragen'):
        assert nextflow_work.cleanup(isolated, tertiary=True)['deleted_count'] == 0
    assert (paths['stage']/'active-data').read_text() == 'keep'


def test_site_script_cannot_redirect_temp_and_cache_but_keeps_tool_paths_and_binds(isolated, tmp_path):
    script = tmp_path/'site.sh'
    script.write_text('export TMPDIR=/wrong NXF_TEMP=/wrong APPTAINER_TMPDIR=/wrong\n'
                      'export NXF_CACHE_DIR=/wrong\n'
                      'export APPTAINER_BIND=/site-reference:/site-reference:ro\n'
                      'export PATH="/site-tools:$PATH"\n')
    output = tmp_path/'environment.json'
    code = 'import os,json,pathlib; pathlib.Path('+repr(str(output))+').write_text(json.dumps(dict(os.environ)))'
    command = workspace.nextflow_command([sys.executable, '-c', code], isolated, 'nckuh', 'job',
                                         str(script), sys.executable, Path(__file__).parents[1]/'backend')
    result = subprocess.run(command, capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stderr
    env = json.loads(output.read_text())
    expected = workspace.runtime_paths(isolated, 'nckuh', 'job')['tmp']
    for key in ('TMPDIR', 'NXF_TEMP', 'APPTAINER_TMPDIR', 'SINGULARITYENV_TMPDIR'):
        assert env[key] == str(expected)
    assert env['NXF_CACHE_DIR'] != '/wrong'
    assert env['PATH'].startswith('/site-tools:')
    assert '/site-reference:/site-reference:ro' in env['APPTAINER_BINDPATH']


def test_verified_copy_materializes_links_and_preserves_nanosecond_times(tmp_path):
    source, target = tmp_path/'source', tmp_path/'destination'
    source.mkdir()
    work = tmp_path/'work'; work.mkdir()
    raw = work/'raw.tsv'; raw.write_text('immutable raw\n')
    os.utime(raw, ns=(1_700_000_000_123456789, 1_700_000_000_123456789))
    (source/'linked-file').symlink_to(raw)
    (source/'linked-directory').symlink_to(work, target_is_directory=True)
    stats = workspace.verified_copytree(source, target)
    assert stats == dict(files=2, bytes=2*raw.stat().st_size)
    assert not (target/'linked-file').is_symlink()
    assert not (target/'linked-directory').is_symlink()
    assert (target/'linked-file').stat().st_mtime_ns == raw.stat().st_mtime_ns
    shutil.rmtree(source); shutil.rmtree(work)
    assert (target/'linked-directory/raw.tsv').read_text() == 'immutable raw\n'


@pytest.mark.parametrize('failure', ['corruption', 'disk_full', 'cancel'])
def test_transfer_failure_removes_partial_copy_and_never_touches_live(tmp_path, monkeypatch, failure):
    source, target = tmp_path/'source', tmp_path/'publication/job/CASE'
    source.mkdir(); (source/'raw').write_text('new')
    live = tmp_path/'live'; live.mkdir(); (live/'raw').write_text('old')
    copy = workspace.shutil.copy2
    def fail(src, dst, **options):
        copy(src, dst, **options)
        if failure == 'corruption': Path(dst).write_text('bad')
        elif failure == 'disk_full': raise OSError(errno.ENOSPC, 'disk full')
        else: raise KeyboardInterrupt()
    monkeypatch.setattr(workspace.shutil, 'copy2', fail)
    with pytest.raises((RuntimeError, OSError, KeyboardInterrupt)):
        workspace.verified_copytree(source, target)
    assert not target.exists()
    assert (live/'raw').read_text() == 'old'


def test_runtime_root_rejects_symlink_and_escaped_job_directory(tmp_path):
    real = tmp_path/'real'; real.mkdir()
    alias = tmp_path/'alias'; alias.symlink_to(real, target_is_directory=True)
    with pytest.raises(ValueError): workspace.runtime_paths(alias, 'dragen', 'job')
    with pytest.raises(ValueError): workspace.runtime_paths(real, 'dragen', '../escape')
    context = real/'contexts/shared-dragen'; context.mkdir(parents=True)
    (context/'.runtime').symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(ValueError): workspace.runtime_paths(real, 'dragen', 'job')


def test_overlay_and_rescue_signatures_remain_current_after_verified_transfer(tmp_path):
    source, copied, live = (tmp_path/name for name in ('ssd', 'publication', 'live'))
    raw = source/'03_acmg/SRC.snv_indel.acmg.tsv'
    raw.parent.mkdir(parents=True)
    raw.write_text('CHROM\tPOS\tREF\tALT\n1\t10\tA\tG\n')
    post = source/'08_postprocessing'; post.mkdir()
    overlay = post/'CASE.snv_annotations.sqlite'
    snv_overlay.build_overlay(raw, raw, overlay)
    (post/'CASE.snv_indel.review.tsv.source.json').write_text(json.dumps({'overlay': snv_overlay.overlay_signature(raw, overlay)}))
    baseline = source/'06_cnv_sv/SRC.cnv.tsv'; baseline.parent.mkdir(); baseline.write_text('cnv')
    (post/'CASE.cnv_rescue.json').write_text(json.dumps({'inputs': {'base_tsv': {'path': str(baseline)}}}))
    workspace.verified_copytree(source, copied)
    final_raw = live/raw.relative_to(source)
    final_post = live/'08_postprocessing'
    dragen_run._rebase_staged_derived_paths(sample_id='CASE', stage_post_dir=copied/'08_postprocessing',
        final_raw_tsv=final_raw, final_post_dir=final_post, source_sample_dir=source)
    copied.rename(live)
    shutil.rmtree(source)
    final_overlay = final_post/overlay.name
    assert snv_overlay.is_current(final_raw, final_overlay)
    assert json.loads((final_post/'CASE.snv_indel.review.tsv.source.json').read_text())['overlay'] == snv_overlay.overlay_signature(final_raw, final_overlay)
    assert json.loads((final_post/'CASE.cnv_rescue.json').read_text())['inputs']['base_tsv']['path'] == str(live/baseline.relative_to(source))


@pytest.mark.parametrize('failure', [None, 'second_copy', 'second_promotion', 'rollback_failure', 'prs_cancel'])
@pytest.mark.parametrize('with_prs', [False, True])
def test_worker_transfers_entire_batch_before_atomic_publication(tmp_path, monkeypatch, failure, with_prs):
    root = tmp_path.resolve()
    work, live, jobs = root/'ssd-work', root/'hdd-live', root/'jobs'
    for name, value in (('TERTIARY_NF_WORK_ROOT', work), ('PIPELINE_OUT_ROOT', live),
                        ('TERTIARY_JOBS_DIR', jobs)):
        monkeypatch.setattr(dragen_run, name, value)
    monkeypatch.setattr(sample_layout.config, 'PIPELINE_OUT_ROOT', live)
    monkeypatch.setenv('NXF_HOME', str(root/'empty-nextflow'))
    monkeypatch.setattr(sys, 'argv', ['worker', '--job-id', 'job', '--mode', 'inhouse', '--without-pgx']+(['--with-prs'] if with_prs else []))
    samples = [dict(mode='inhouse', sample_id=f'S{i}-nckuh', source_sample_id=f'S{i}',
                    seq_type='WGS', vcf_path=str(root/f'S{i}.vcf.gz'), cnv_vcf='', sv_vcf='') for i in (1, 2)]
    monkeypatch.setattr(dragen_run, '_load_samples', lambda args: samples)
    monkeypatch.setattr(dragen_run, '_load_secrets', lambda: None)
    monkeypatch.setattr(dragen_run, '_seed_shared_nextflow_context', lambda *a, **kw: None)
    monkeypatch.setattr(dragen_run, '_validate_acmg_tsv', lambda *a, **kw: None)
    monkeypatch.setattr(dragen_run, '_ensure_pipeline_annotation_versions', lambda p: p)
    monkeypatch.setattr(dragen_run.specialized, 'prepare_outputs', lambda **kw: dict(status='skipped', variants=[]))
    monkeypatch.setattr(dragen_run, '_copy_ploidy_artifacts', lambda **kw: (None, None))
    monkeypatch.setattr(dragen_run.roh, 'prepare_roh_outputs', lambda **kw: {})
    monkeypatch.setattr(dragen_run.mitomap_mito, 'annotate_mito_tsv', lambda p: False)
    monkeypatch.setattr(dragen_run, '_record_gpn_completion', lambda *a: None)
    state = {}
    monkeypatch.setattr(dragen_run, '_update', lambda jid, **kw: state.update(kw))
    monkeypatch.setattr(dragen_run, '_set_step', lambda jid, step, **kw: state.update(step=step, **kw))
    prs_calls=[]
    def prs_postprocessing(samples,jid):
        assert state['state']=='running' and state['step']=='promote-output'
        assert len(markers)==2
        with pytest.raises(RuntimeError): dragen_run._acquire_sample_locks([samples[0]['sample_id']])
        prs_calls.append(jid)
        if failure=='prs_cancel': raise SystemExit(143)
    monkeypatch.setattr(dragen_run,'_run_prs_postprocessing',prs_postprocessing)
    from app.services import sample_loader, unregistered_samples
    monkeypatch.setattr(sample_loader, 'invalidate_sample_cache', lambda p: None)
    monkeypatch.setattr(sample_loader, 'update_case_table_row', lambda sid: None)
    monkeypatch.setattr(unregistered_samples, 'refresh_samples', lambda ids: None)
    monkeypatch.setattr(dragen_run.dragen_jobs, 'nextflow_log_signature', lambda p: None)
    monkeypatch.setattr(dragen_run.dragen_jobs, 'snapshot_nextflow_log', lambda *a: None)
    for sample in samples:
        old = live/sample['sample_id']
        (old/'03_acmg').mkdir(parents=True)
        (old/'03_acmg/old').write_text('old')
        post = old/'08_postprocessing'; post.mkdir()
        (post/(sample['sample_id']+'.sample_metadata.json')).write_text('{"reviewer":"keep"}')
    def run(cmd, *, label, **kw):
        if label.startswith('2b/4'):
            stage = Path(state['staging_root'])
            assert stage.is_relative_to(work) and not stage.is_relative_to(live)
            for sample in samples:
                sid = sample['source_sample_id']
                for directory in dragen_run.REQUIRED_PIPELINE_STAGE_DIRS:
                    (stage/sid/directory).mkdir(parents=True)
                for directory, name in (('03_acmg', 'snv_indel.acmg'), ('04_mito', 'mito'),
                                        ('05_str', 'str'), ('06_cnv_sv', 'cnv.annotated'),
                                        ('06_cnv_sv', 'sv.annotated')):
                    (stage/sid/directory/f'{sid}.{name}.tsv').write_text('CHROM\tPOS\tREF\tALT\n')
        elif label.startswith('post-processing'):
            assert Path(cmd[cmd.index('--work-tsv')+1]).is_relative_to(work)
    monkeypatch.setattr(dragen_run, '_run', run)
    copy = workspace.verified_copytree
    copied = []
    def transfer(source, target):
        if failure == 'second_copy' and copied: raise OSError(errno.ENOSPC, 'full')
        result = copy(source, target); copied.append(target)
        return result
    monkeypatch.setattr(workspace, 'verified_copytree', transfer)
    replace = os.replace
    markers = []
    interrupted = False
    def replace_on_hdd(src, dst):
        nonlocal interrupted
        src, dst = Path(src), Path(dst)
        if dst.is_relative_to(live) and not dst.is_relative_to(live/'.staging') and not dst.is_relative_to(live/'.rollback'):
            if src.is_relative_to(work): raise OSError(errno.EXDEV, 'cross-device rename')
            assert len(copied) == 2
            if (failure in {'second_promotion', 'rollback_failure'}
                    and (not interrupted or failure == 'rollback_failure') and dst == live/'S2-nckuh/03_acmg'):
                interrupted = True
                raise OSError('publication interrupted')
            if dst.name.endswith('.layout.json'): markers.append(dst)
        return replace(src, dst)
    monkeypatch.setattr(os, 'replace', replace_on_hdd)
    publication_completed = failure in (None, 'prs_cancel')
    cancelled = with_prs and failure=='prs_cancel'
    assert dragen_run.main() == (143 if cancelled else (0 if publication_completed else 1))
    assert prs_calls == (['job'] if with_prs and publication_completed else [])
    if cancelled: assert state['state']=='cancelled'
    for sample in samples:
        directory = live/sample['sample_id']
        assert json.loads((directory/'08_postprocessing'/f'{sample["sample_id"]}.sample_metadata.json').read_text()) == {'reviewer': 'keep'}
        if publication_completed:
            assert (directory/'03_acmg'/f'{sample["source_sample_id"]}.snv_indel.acmg.tsv').exists()
        elif failure == 'rollback_failure' and sample['sample_id'] == 'S2-nckuh':
            assert (live/'.rollback/job/S2-nckuh/03_acmg/old').read_text() == 'old'
            assert 'rollback incomplete' in state['error']
        else:
            assert (directory/'03_acmg/old').read_text() == 'old'
    assert len(markers) == (2 if publication_completed else (1 if failure in {'second_promotion', 'rollback_failure'} else 0))
    assert (live/'.staging/job').exists() == (failure == 'rollback_failure' or cancelled)
    assert (live/'.rollback/job').exists() == (failure == 'rollback_failure' or cancelled)
    assert not workspace.runtime_paths(work, 'nckuh', 'job')['job'].exists()


@pytest.mark.parametrize('profile', [False, True])
def test_real_nextflow_config_preserves_site_options_and_routes_both_engines(isolated, tmp_path, profile):
    executable = os.environ.get('NGS_UI_TEST_NEXTFLOW')
    if not executable:
        pytest.skip('Set NGS_UI_TEST_NEXTFLOW for actual Nextflow configuration checks')
    source = tmp_path/'site.config'
    source.write_text("singularity.runOptions = '--contain'\n"
                      "profiles { dgm { singularity.runOptions = '--no-home'; params.clinicalSetting = 'keep' } }\n")
    env = {**os.environ, **workspace.environment(isolated, 'dragen', 'job'), 'NXF_OFFLINE': 'true'}
    wrapper = workspace.nextflow_config(source, isolated, 'dragen', 'job')
    cmd = [executable, '-c', str(wrapper), 'config', '-flat']
    if profile: cmd += ['-profile', 'dgm']
    result = subprocess.run(cmd, cwd=tmp_path, env=env, capture_output=True, text=True, timeout=40)
    assert result.returncode == 0, result.stderr
    options = next(line for line in result.stdout.splitlines() if line.startswith('singularity.runOptions ='))
    assert ('--no-home' if profile else '--contain') in options
    assert options.count(env['TMPDIR']+':/tmp') == 1
    assert options.count(env['TMPDIR']+':/var/tmp') == 1
    if profile: assert "params.clinicalSetting = 'keep'" in result.stdout


def test_real_nextflow_resumes_after_job_tmp_and_publication_paths_change(isolated, tmp_path):
    executable = os.environ.get('NGS_UI_TEST_NEXTFLOW')
    if not executable:
        pytest.skip('Set NGS_UI_TEST_NEXTFLOW for actual Nextflow resume checks')
    project = tmp_path/'pipeline'; project.mkdir()
    (project/'main.nf').write_text('''process COPY {
    input: path seed
    output: path 'result.txt'
    script:
    """
    cat ${seed} > result.txt
    """
}
workflow { COPY(file(params.seed)) }
''')
    source = project/'nextflow.config'; source.write_text('')
    seed = tmp_path/'seed'; seed.write_text('same source\n')
    launch = isolated/'contexts/shared-dragen/launch'; launch.mkdir(parents=True)
    task_work = isolated/'contexts/shared-dragen/work'
    for number in (1, 2):
        jid = f'job{number}'
        env = {**os.environ, **workspace.environment(isolated, 'dragen', jid), 'NXF_OFFLINE': 'true'}
        wrapper = workspace.nextflow_config(source, isolated, 'dragen', jid)
        trace = tmp_path/f'trace{number}.txt'
        cmd = [executable, '-c', str(wrapper), 'run', str(project), '-profile', 'dgm',
               '-work-dir', str(task_work), '--seed', str(seed), '-with-trace', str(trace)]
        if number == 2: cmd.append('-resume')
        result = subprocess.run(cmd, cwd=launch, env=env, capture_output=True, text=True, timeout=40)
        assert result.returncode == 0, result.stdout+'\n'+result.stderr
        assert ('CACHED' in trace.read_text()) == (number == 2)
        # This is disposable after each job; cached outputs must live in work.
        shutil.rmtree(workspace.runtime_paths(isolated, 'dragen', jid)['job'])
    result = next(task_work.glob('*/*/result.txt'))
    assert result.read_text() == 'same source\n'


def test_real_nextflow_container_launcher_forwards_tmp_after_env_cleaning(isolated, tmp_path):
    executable = os.environ.get('NGS_UI_TEST_NEXTFLOW')
    if not executable:
        pytest.skip('Set NGS_UI_TEST_NEXTFLOW for actual Nextflow launcher checks')
    tools = tmp_path/'bin'; tools.mkdir()
    engine = tools/'singularity'
    engine.write_text(f'#!{sys.executable}\n'+
        'import os,sys\na=sys.argv[1:]\n'
        'if a==["--version"]:\n print("apptainer version 1.3.0"); sys.exit(0)\n'
        'assert a.pop(0)=="exec",a\n'
        'while a and a[0].startswith("-"):\n'
        ' k=a.pop(0)\n'
        ' if k=="--env":\n'
        '  for entry in a.pop(0).split(","):\n'
        '   key,value=entry.split("=",1); os.environ[key]=value\n'
        ' elif k in ("-B","--bind","--no-mount"): a.pop(0)\n'
        ' elif k not in ("--no-home","--pid","--contain"): raise ValueError(k)\n'
        'a.pop(0)\nos.execvpe(a[0],a,os.environ)\n')
    engine.chmod(0o755)
    image = tmp_path/'fake.sif'; image.touch()
    project = tmp_path/'pipeline'; project.mkdir()
    (project/'main.nf').write_text('''process PROBE {
    output: path 'temp.json'
    script:
    """
    python3 - <<'PY'
import json,tempfile
from pathlib import Path
with tempfile.NamedTemporaryFile() as f:
    Path('temp.json').write_text(json.dumps({'temporary_file':f.name}))
PY
    """
}
workflow { PROBE() }
''')
    source = project/'nextflow.config'
    source.write_text(f"process.container = '{image}'\nprofiles {{ dgm {{ singularity.enabled = true }} }}\n")
    env_script = tmp_path/'site.sh'
    env_script.write_text(f'export PATH="{tools}:{Path(sys.executable).parent}:$PATH"\nexport TMPDIR=/wrong NXF_TEMP=/wrong\n')
    env = {**os.environ, **workspace.environment(isolated, 'dragen', 'job'), 'NXF_OFFLINE': 'true'}
    wrapper = workspace.nextflow_config(source, isolated, 'dragen', 'job')
    work = isolated/'contexts/shared-dragen/work'
    cmd = workspace.nextflow_command([executable, '-c', str(wrapper), 'run', str(project), '-profile', 'dgm', '-w', str(work)],
        isolated, 'dragen', 'job', str(env_script), sys.executable, Path(__file__).parents[1]/'backend')
    result = subprocess.run(cmd, cwd=tmp_path, env=env, capture_output=True, text=True, timeout=40)
    assert result.returncode == 0, result.stdout+'\n'+result.stderr
    output = next(work.glob('*/*/temp.json'))
    assert Path(json.loads(output.read_text())['temporary_file']).parent == Path(env['TMPDIR'])
    launcher = output.with_name('.command.run').read_text()
    assert 'env - PATH=' in launcher
    assert f'{isolated}:{isolated}' in launcher
    assert f'{env["TMPDIR"]}:/tmp' in launcher and f'{env["TMPDIR"]}:/var/tmp' in launcher
