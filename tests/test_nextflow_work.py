import subprocess
import sys
import threading
from pathlib import Path
from types import SimpleNamespace

import pytest

from app.services import dragen_jobs, nextflow_work as work


@pytest.fixture(autouse=True)
def no_external_nextflow(monkeypatch):
    monkeypatch.setattr(work, 'nextflow_process_scopes', lambda root: (set(), []))


def task(path):
    path.mkdir(parents=True)
    (path / '.command.out').write_text('output')
    return path


def context(root, scope, target=None):
    path = task(root / 'contexts' / f'shared-{scope}')
    target = target or task(path / 'work')
    (path / 'work-dir.txt').write_text(str(target.resolve()))
    return path


def test_busy_secondary_batch_keeps_entire_tree_and_cleans_idle(tmp_path):
    root = tmp_path / 'work'
    busy, idle = task(root / 'busy'), task(root / 'idle')
    with work.acquire_lease(root, 'busy'):
        result = work.cleanup(root)
    assert busy.is_dir() and not idle.exists()
    assert result['deleted_count'] == result['skipped_count'] == 1
    assert result['failed_count'] == 0
    assert work._lock_path(root, 'busy').exists()
    assert work.cleanup(root)['deleted_count'] == 1


def test_inherited_child_lease_survives_parent_close(tmp_path):
    root = tmp_path / 'work'
    batch = task(root / 'batch')
    lease = work.acquire_lease(root, 'batch')
    proc = subprocess.Popen([sys.executable, '-c', 'print("ready", flush=True); input()'],
                            pass_fds=(lease.fileno(),), stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE, text=True)
    try:
        assert proc.stdout.readline().strip() == 'ready'
        lease.close()
        assert work.cleanup(root)['skipped_count'] == 1
        assert batch.exists()
    finally:
        proc.communicate('\n', timeout=10)
        lease.close()
    assert work.cleanup(root)['deleted_count'] == 1


def test_cleanup_admission_blocks_new_launch_until_deletion_finishes(monkeypatch, tmp_path):
    root = tmp_path / 'work'
    batch = task(root / 'batch')
    entered, acquired, release = threading.Event(), threading.Event(), threading.Event()
    def launch():
        entered.set()
        with work.acquire_lease(root, 'batch'):
            acquired.set()
            release.wait(5)
    thread = threading.Thread(target=launch)
    def scan(*args):
        thread.start()
        assert entered.wait(2)
        assert not acquired.wait(.1)
        return set()
    monkeypatch.setattr(work, '_references', scan)
    try:
        result = work.cleanup(root, output_roots=[tmp_path / 'output'])
        assert result['deleted_count'] == 1 and not batch.exists()
        assert acquired.wait(2)
    finally:
        release.set()
        thread.join(5)
    assert not thread.is_alive()


def test_active_tertiary_protects_adopted_work_and_dependent_lineage(tmp_path):
    root = tmp_path / 'work'
    old = task(root / 'legacy' / 'work')
    dragen = context(root, 'dragen', old)
    nckuh = context(root, 'nckuh')
    idle = task(root / 'unreferenced')
    with work.acquire_lease(root, 'dragen'):
        result = work.cleanup(root, tertiary=True)
    assert dragen.exists() and old.exists()
    assert not nckuh.exists() and not idle.exists()
    assert result['skipped_count'] == 2
    assert result['deleted_count'] == 2


def test_two_lineages_pointing_at_same_work_are_kept_together(tmp_path):
    root = tmp_path / 'work'
    old = task(root / 'legacy')
    a, b = context(root, 'dragen', old), context(root, 'nckuh', old)
    with work.acquire_lease(root, 'nckuh'):
        result = work.cleanup(root, tertiary=True)
    assert result['deleted_count'] == 0
    assert all(p.exists() for p in (old, a, b))


@pytest.mark.parametrize('initialized', [False, True])
def test_queued_initialization_keeps_possible_seed_sources(tmp_path, initialized):
    root = tmp_path / 'work'
    old = task(root / 'legacy')
    if initialized:
        task(root / 'contexts' / 'shared-dragen')
    with work.acquire_lease(root, 'dragen'):
        result = work.cleanup(root, tertiary=True)
    assert old.exists() and result['deleted_count'] == 0


def test_output_symlink_preserves_work_and_cache_metadata(tmp_path):
    root = tmp_path / 'work'
    old = task(root / 'legacy')
    cache = context(root, 'dragen', old)
    output = tmp_path / 'output'
    output.mkdir()
    (output / 'result.bam').symlink_to(old / '.command.out')
    idle = task(root / 'idle')
    result = work.cleanup(root, tertiary=True, output_roots=[output])
    assert old.exists() and cache.exists() and not idle.exists()
    assert (output / 'result.bam').read_text() == 'output'
    assert result['skipped_count'] == 2


def test_active_work_input_symlink_protects_another_batch(tmp_path):
    root = tmp_path / 'work'
    producer, consumer = task(root / 'producer'), task(root / 'consumer')
    (consumer / 'input').symlink_to(producer / '.command.out')
    with work.acquire_lease(root, 'consumer'):
        result = work.cleanup(root)
    assert result['deleted_count'] == 0
    assert producer.exists()


def test_chained_published_symlinks_keep_all_intermediate_batches(tmp_path):
    root = tmp_path / 'work'
    a, b, c = [task(root / name) for name in ('a', 'b', 'c')]
    (a / 'linkdir').symlink_to(b, target_is_directory=True)
    (b / 'linkfile').symlink_to(c / '.command.out')
    output = tmp_path / 'output'
    output.mkdir()
    (output / 'result').symlink_to(a / 'linkdir' / 'linkfile')
    result = work.cleanup(root, output_roots=[output])
    assert result['deleted_count'] == 0
    assert result['skipped_count'] == 3
    assert (output / 'result').read_text() == 'output'


def test_mount_alias_does_not_hide_intermediate_work_reference(tmp_path):
    root = tmp_path / 'real' / 'work'
    a, b = task(root / 'a'), task(root / 'b')
    alias = tmp_path / 'alias'
    alias.symlink_to(root.parent, target_is_directory=True)
    (a / 'linkdir').symlink_to(b, target_is_directory=True)
    output = tmp_path / 'output'
    output.mkdir()
    (output / 'result').symlink_to(alias / 'work' / 'a' / 'linkdir' / '.command.out')
    result = work.cleanup(root, output_roots=[output])
    assert result['deleted_count'] == 0
    assert result['skipped_count'] == 2
    assert (output / 'result').read_text() == 'output'


def test_unknown_nextflow_keeps_everything(monkeypatch, tmp_path):
    root = tmp_path / 'work'
    batch = task(root / 'batch')
    monkeypatch.setattr(work, 'nextflow_process_scopes', lambda root: (set(), ['1234']))
    result = work.cleanup(root)
    assert batch.exists() and result['skipped_count'] == 1
    assert '1234' in result['skipped'][0]['reason']


def test_process_inspection_failure_never_deletes(monkeypatch, tmp_path):
    root = tmp_path / 'work'
    batch = task(root / 'batch')
    def denied(root):
        raise PermissionError('cannot inspect processes')
    monkeypatch.setattr(work, 'nextflow_process_scopes', denied)
    with pytest.raises(PermissionError):
        work.cleanup(root)
    assert batch.exists()


def test_root_symlink_is_rejected(tmp_path):
    target = task(tmp_path / 'real' / 'batch')
    alias = tmp_path / 'alias'
    alias.symlink_to(target.parent, target_is_directory=True)
    with pytest.raises(ValueError):
        work.cleanup(alias)
    assert target.exists()


def test_invalid_work_pointer_preserves_all_candidates(tmp_path):
    root = tmp_path / 'work'
    shared = context(root, 'dragen')
    idle = task(root / 'idle')
    (shared / 'work-dir.txt').write_text('/somewhere/outside/root')
    result = work.cleanup(root, tertiary=True)
    assert result['deleted_count'] == 0
    assert idle.exists()


def test_unreadable_references_keep_everything(monkeypatch, tmp_path):
    root = tmp_path / 'work'
    task(root / 'batch')
    def denied(*args):
        raise PermissionError('cannot inspect outputs')
    monkeypatch.setattr(work, '_references', denied)
    assert work.cleanup(root, output_roots=[tmp_path / 'output'])['skipped_count'] == 1


def test_delete_failure_is_reported_without_false_success(monkeypatch, tmp_path):
    root = tmp_path / 'work'
    task(root / 'batch')
    def denied(*args):
        raise PermissionError('cannot delete')
    monkeypatch.setattr(work.shutil, 'rmtree', denied)
    result = work.cleanup(root)
    assert result['deleted_count'] == result['freed_bytes'] == 0
    assert result['failed_count'] == 1


@pytest.mark.parametrize('value', ['', '../work', '/'])
def test_unsafe_root_rejected(value):
    with pytest.raises(ValueError):
        work.cleanup(value)


def test_tertiary_wrapper_respects_old_cache_lock(monkeypatch, tmp_path):
    root = tmp_path / 'work'
    active, idle = context(root, 'dragen'), context(root, 'nckuh')
    monkeypatch.setattr(dragen_jobs, 'TERTIARY_NF_WORK_ROOT', root)
    monkeypatch.setattr(dragen_jobs, 'TERTIARY_JOBS_DIR', tmp_path / 'jobs')
    monkeypatch.setattr(dragen_jobs, '_active_tertiary_jobs', lambda: [])
    monkeypatch.setattr(dragen_jobs, 'PIPELINE_OUT_ROOT', tmp_path / 'output')
    monkeypatch.setattr(dragen_jobs, 'LEGACY_PIPELINE_OUT_ROOT', tmp_path / 'old-output')
    with work.open_lock(tmp_path / 'jobs' / '.nextflow_locks' / 'dragen.lock'):
        result = dragen_jobs.cleanup_nf_work()
    assert active.exists() and not idle.exists()
    assert result['skipped_count'] == result['deleted_count'] == 1


def test_tertiary_wrapper_preserves_waiting_or_postprocessing_jobs(monkeypatch, tmp_path):
    root = tmp_path / 'work'
    idle, active = context(root, 'dragen'), context(root, 'nckuh')
    monkeypatch.setattr(dragen_jobs, 'TERTIARY_NF_WORK_ROOT', root)
    monkeypatch.setattr(dragen_jobs, 'TERTIARY_JOBS_DIR', tmp_path / 'jobs')
    monkeypatch.setattr(dragen_jobs, '_active_tertiary_jobs', lambda: [{'mode': 'inhouse'}])
    monkeypatch.setattr(dragen_jobs, 'PIPELINE_OUT_ROOT', tmp_path / 'output')
    monkeypatch.setattr(dragen_jobs, 'LEGACY_PIPELINE_OUT_ROOT', tmp_path / 'old-output')
    result = dragen_jobs.cleanup_nf_work()
    assert active.exists() and not idle.exists()
    assert result['skipped_count'] == result['deleted_count'] == 1


def test_start_job_passes_lifetime_lease_to_detached_worker(monkeypatch, tmp_path):
    root = tmp_path / 'work'
    ctx = context(root, 'dragen')
    vcf = tmp_path / 'S1.hard-filtered.vcf.gz'
    vcf.touch()
    monkeypatch.setattr(dragen_jobs, 'TERTIARY_NF_WORK_ROOT', root)
    monkeypatch.setattr(dragen_jobs, 'TERTIARY_JOBS_DIR', tmp_path / 'jobs')
    monkeypatch.setattr(dragen_jobs, 'LEGACY_DRAGEN_JOBS_DIR', tmp_path / 'legacy-jobs')
    monkeypatch.setattr(dragen_jobs, 'active_sample_ids', lambda: set())
    popen = subprocess.Popen
    processes = []
    def spawn(command, **kwargs):
        child = popen([sys.executable, '-c',
                       'import os; os.fstat(int(os.environ["NGS_UI_NF_GUARD_FD"])); '
                       'print("ready", flush=True); input()'],
                      **{**kwargs, 'stdin': subprocess.PIPE, 'stdout': subprocess.PIPE, 'text': True})
        processes.append(child)
        return child
    monkeypatch.setattr(dragen_jobs.subprocess, 'Popen', spawn)
    try:
        job_id = dragen_jobs.start_job(str(vcf), 'S1', mode='dragen')
        child = processes[0]
        assert child.stdout.readline().strip() == 'ready'
        assert dragen_jobs.load_state(job_id)['state'] == 'queued'
        assert work.cleanup(root, tertiary=True)['deleted_count'] == 0
        assert ctx.exists()
    finally:
        for child in processes:
            child.communicate('\n', timeout=10)
    assert work.cleanup(root, tertiary=True)['deleted_count'] == 1


def test_worker_commands_inherit_work_lease(monkeypatch, tmp_path):
    from app.workers import dragen_run
    root = tmp_path / 'work'
    with work.acquire_lease(root, 'dragen') as lease:
        for key, value in work.lease_env(root, 'dragen', lease).items():
            monkeypatch.setenv(key, value)
        dragen_run._run([sys.executable, '-c',
                         'import os; os.fstat(int(os.environ["NGS_UI_NF_GUARD_FD"]))'], label='lease')


def test_process_scanner_does_not_mistake_embedded_python_for_nextflow(monkeypatch, tmp_path):
    # Restore the real scanner under the autouse override.
    import importlib.util
    spec = importlib.util.spec_from_file_location('guard_scanner', work.__file__)
    scanner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scanner)
    root = tmp_path / 'work'
    monkeypatch.setattr(scanner.subprocess, 'run', lambda *a, **kw: SimpleNamespace(stdout=(
        '99 python3 python3 -c source_with_nextflow.cli.Launcher cleanup\n'
        '123 java java -cp nextflow.jar nextflow.cli.Launcher run main.nf\n'
        '124 bash bash -ue /work/ab/xyz/.command.run\n')))
    monkeypatch.setattr(scanner.os, 'kill', lambda *a: None)
    assert scanner.nextflow_process_scopes(root) == (set(), ['123', '124'])


def test_managed_process_maps_to_live_lease(monkeypatch, tmp_path):
    import importlib.util
    spec = importlib.util.spec_from_file_location('guard_scanner', work.__file__)
    scanner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scanner)
    root = tmp_path / 'work'
    monkeypatch.setattr(scanner.subprocess, 'run', lambda *a, **kw: SimpleNamespace(
        stdout='123 java java nextflow.cli.Launcher run main.nf\n'))
    with work.acquire_lease(root, 'batch') as lease:
        data = b'\0'.join(f'{k}={v}'.encode() for k, v in work.lease_env(root, 'batch', lease).items())
        monkeypatch.setattr(Path, 'read_bytes', lambda path: data)
        assert scanner.nextflow_process_scopes(root) == ({'batch'}, [])
