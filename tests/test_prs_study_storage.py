"""Storage cutover safeguards: synthetic trees, real rsync, no actual mounts."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

spec = importlib.util.spec_from_file_location('prs_storage', Path(__file__).parents[1] / 'scripts/migrate_prs_study_storage.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
STAMP = '20261120T010203123456Z'


@pytest.fixture
def study(tmp_path):
    root = tmp_path.resolve()
    source = root / 'work/prs/_validation_test/cohort_v1'
    source.mkdir(parents=True)
    (source / 'execution.lock').touch()
    (source / 'frozen').mkdir()
    (source / 'frozen/reference_scores.json').write_text('{"synthetic": true}')
    manifest = {'sha256': {'reference_scores.json': m.digest(source / 'frozen/reference_scores.json')}}
    (source / 'frozen/manifest.json').write_text(json.dumps(manifest))
    (source / 'plan.json').write_text(json.dumps(dict(schema='PRS_COHORT_STUDY_V1', study=str(source),
                                                    model_sha256=m.digest(source / 'frozen/manifest.json'))))
    target = root / 'hdd/_validation_test/cohort_v1'
    target.parent.mkdir(parents=True)
    journal = target.parent / '.cohort_v1.migration.json'
    record = dict(schema=m.SCHEMA, source=str(source), target=str(target), phase='verified', stamp=STAMP,
                  backup=str(source.with_name(source.name + '.ssd-backup-' + STAMP)),
                  source_identity=m.identity(source), plan_sha256=m.digest(source / 'plan.json'))
    return source, target, journal, record


def copied(study):
    source, target, journal, record = study
    shutil.copytree(source, target, symlinks=True)
    record['target_tree'] = m.tree_signature(target)
    return study


def fake_mount(monkeypatch, source, target, commands, *, fail=False):
    def run(argv, **kwargs):
        commands.append(list(map(str, argv)))
        if argv[0] == 'mount':
            if fail: raise RuntimeError('mount failed')
            source.rmdir(); source.symlink_to(target, target_is_directory=True)
        return ''
    monkeypatch.setattr(m, 'run', run)
    monkeypatch.setattr(m.os.path, 'ismount', lambda p: Path(p) == source and source.is_symlink())


def mounted(study, monkeypatch):
    source, target, journal, record = copied(study)
    commands = []
    fake_mount(monkeypatch, source, target, commands)
    fstab = target.parent / 'fstab'; fstab.write_text('# original\n')
    m.install(source, target, record, journal, fstab=fstab)
    return source, target, journal, record, fstab, commands


def test_paths_reject_ui_cases_and_work_root_replacement(monkeypatch):
    monkeypatch.setattr(m, 'canonical', lambda p: p)
    work = Path('/home/user/nf_work')
    good = work / 'prs/_validation_test/cohort_v1'
    target = Path('/home/user/prs_validation_work/_validation_test/cohort_v1')
    m.shape(work, good, target)
    for bad in (work, work / 'prs/CASE', work / 'prs/_validation_test/cohort_v1/cases'):
        with pytest.raises(ValueError): m.shape(work, bad, target)
    for bad in (good, work / 'other/_validation_test/cohort_v1', Path('/tmp/_validation_test/cohort_v1'),
                target.with_name('different')):
        with pytest.raises(ValueError): m.shape(work, good, bad)


def test_alias_and_unusual_paths_are_rejected(tmp_path):
    root = tmp_path.resolve(); real = root / 'real'; real.mkdir()
    alias = root / 'alias'; alias.symlink_to(real, target_is_directory=True)
    for path in (alias, alias / 'child', root / 'has space', root / 'real/../other'):
        with pytest.raises(ValueError): m.canonical(path)


def test_frozen_inputs_and_logical_path_must_remain_identical(study):
    source, target, journal, record = copied(study)
    m.check_plan(target, source, record['plan_sha256'])
    with pytest.raises(ValueError): m.check_plan(target, target)
    (target / 'frozen/reference_scores.json').write_text('changed')
    with pytest.raises(ValueError, match='分布'): m.check_plan(target, source)


def test_rsync_preserves_hardlinks_symlinks_sparse_and_hidden_files(study):
    if not shutil.which('rsync'): pytest.skip('rsync unavailable')
    source, target, _, _ = study
    (source / '.hidden').write_bytes(b'fixed study metadata')
    os.link(source / '.hidden', source / 'hardlink')
    (source / 'external').symlink_to('/unavailable/reference.pgen')
    with (source / 'sparse').open('wb') as f:
        f.seek(4 * 1024**2); f.write(b'end')
    m.run(m.rsync(source, target))
    m.verify_copy(source, target)
    assert os.path.samefile(target / '.hidden', target / 'hardlink')
    assert os.readlink(target / 'external') == '/unavailable/reference.pgen'
    assert (target / 'sparse').stat().st_size == 4 * 1024**2 + 3
    assert (target / '.hidden').read_bytes() == (source / '.hidden').read_bytes()
    (target / 'extra').write_text('must not be silently removed')
    with pytest.raises(RuntimeError, match='比對不同'): m.verify_copy(source, target)
    assert (target / 'extra').exists() and source.exists()
    assert '--delete' not in m.rsync(source, target)
    assert '--delete' in m.rsync(source, target, verify=True)  # dry-run only


def test_checksum_detects_same_size_and_timestamp_corruption(study):
    if not shutil.which('rsync'): pytest.skip('rsync unavailable')
    source, target, _, _ = study
    file = source / 'genotypes'; file.write_bytes(b'ACGT')
    m.run(m.rsync(source, target))
    old = (target / 'genotypes').stat()
    (target / 'genotypes').write_bytes(b'TGCA')
    os.utime(target / 'genotypes', ns=(old.st_atime_ns, old.st_mtime_ns))
    with pytest.raises(RuntimeError, match='比對不同'): m.verify_copy(source, target)


def test_partial_backup_must_be_subset_not_new_or_changed_data(study):
    if not shutil.which('rsync'): pytest.skip('rsync unavailable')
    source, target, _, _ = study
    m.run(m.rsync(source, target))
    (source / 'execution.lock').unlink()
    m.verify_remaining(source, target)  # Target-only files are expected after a partial removal.
    (source / 'new-unique-file').write_text('do not delete')
    with pytest.raises(RuntimeError, match='剩餘 SSD'): m.verify_remaining(source, target)
    assert (source / 'new-unique-file').exists()


def test_other_scope_runs_while_prs_execution_and_study_stay_locked(study):
    source, target, journal, _ = study; work = source.parents[2]
    with m.guard.acquire_lease(work, 'nckuh'):
        with m.lockset(work, source, target, journal) as hold:
            for path in (work / 'prs/.execution.lock', source / 'execution.lock'):
                with pytest.raises(BlockingIOError): m.guard.open_lock(path, blocking=False)
            target.mkdir(); (target / 'execution.lock').touch()
            hold(m.guard.open_lock(target / 'execution.lock', blocking=False))
            with pytest.raises(BlockingIOError): m.guard.open_lock(target / 'execution.lock', blocking=False)
            # No admission lock is held for the lengthy copy.
            lock = m.guard.open_lock(m.guard._lock_path(work, ''), blocking=False); lock.close()
    assert not m.HELD_FDS


def test_existing_writer_refuses_cutover_and_closes_guard_fds(study):
    source, target, journal, _ = study; work = source.parents[2]
    busy = m.guard.open_lock(source / 'execution.lock', blocking=False)
    try:
        with pytest.raises(BlockingIOError):
            with m.lockset(work, source, target, journal): pass
    finally: busy.close()
    assert source.exists() and not target.exists() and not m.HELD_FDS


def test_mount_failure_restores_original_and_never_changes_fstab(study, monkeypatch):
    source, target, journal, record = copied(study)
    fstab = target.parent / 'fstab'; fstab.write_text('# unchanged\n')
    commands = []; fake_mount(monkeypatch, source, target, commands, fail=True)
    with pytest.raises(RuntimeError, match='mount failed'):
        m.install(source, target, record, journal, fstab=fstab)
    assert m.identity(source) == record['source_identity']
    assert (source / 'plan.json').read_bytes() == (target / 'plan.json').read_bytes()
    assert not Path(record['backup']).exists() and record['phase'] == 'verified'
    assert fstab.read_text() == '# unchanged\n'


def test_mount_keeps_backup_and_installs_only_one_fstab_entry(study, monkeypatch):
    source, target, journal, record, fstab, commands = mounted(study, monkeypatch)
    assert os.path.samefile(source, target) and Path(record['backup']).is_dir()
    assert record['phase'] == m.read(journal)['phase'] == 'ready'
    assert str(target) + ' ' + str(source) + ' none bind,nofail 0 0' in fstab.read_text()
    m.persist(source, target, fstab, STAMP)
    assert fstab.read_text().count(' none bind,nofail 0 0') == 1
    assert all(c[0] != 'systemctl' or c == ['systemctl', 'daemon-reload'] for c in commands)


def test_persistence_failure_retains_mounted_copy_and_ssd_backup(study, monkeypatch):
    source, target, journal, record = copied(study)
    fake_mount(monkeypatch, source, target, [])
    monkeypatch.setattr(m, 'persist', lambda *args: (_ for _ in ()).throw(OSError('fstab write failed')))
    with pytest.raises(OSError): m.install(source, target, record, journal)
    assert os.path.samefile(source, target) and Path(record['backup']).is_dir()
    assert m.read(journal)['phase'] == 'mounted'


def test_release_only_removes_exact_verified_backup(study, monkeypatch):
    source, target, journal, record, fstab, _ = mounted(study, monkeypatch)
    checked = []; monkeypatch.setattr(m, 'verify_copy', lambda *args: checked.append(args))
    m.release(source, target, record, journal, fstab=fstab)
    assert checked == [(Path(record['backup']), target)]
    assert not Path(record['backup']).exists() and record['phase'] == 'released'
    assert (source / 'plan.json').exists() and (target / 'plan.json').exists()


@pytest.mark.parametrize('problem', ['wrong_backup', 'backup_mount', 'not_mounted', 'fstab', 'changed_target', 'checksum'])
def test_release_refuses_unsafe_or_unverified_data(study, monkeypatch, problem):
    source, target, journal, record, fstab, _ = mounted(study, monkeypatch)
    backup = Path(record['backup'])
    monkeypatch.setattr(m, 'verify_copy', lambda *args: None)
    if problem == 'wrong_backup': record['backup'] = str(source.parent)
    elif problem == 'backup_mount': monkeypatch.setattr(m.os.path, 'ismount', lambda p: Path(p) in {source, backup})
    elif problem == 'not_mounted': monkeypatch.setattr(m.os.path, 'ismount', lambda p: False)
    elif problem == 'fstab': fstab.write_text('# removed\n')
    elif problem == 'changed_target': (target / 'new-result').write_text('resumed too early')
    else: monkeypatch.setattr(m, 'verify_copy', lambda *args: (_ for _ in ()).throw(RuntimeError('checksum differs')))
    with pytest.raises(RuntimeError): m.release(source, target, record, journal, fstab=fstab)
    assert backup.exists() and (target / 'plan.json').exists()


def test_interrupted_backup_release_can_resume_without_erasing_hdd(study, monkeypatch):
    source, target, journal, record, fstab, _ = mounted(study, monkeypatch)
    monkeypatch.setattr(m, 'verify_copy', lambda *args: None)
    real_delete = m.shutil.rmtree
    def partial(path):
        (path / 'execution.lock').unlink()
        raise KeyboardInterrupt()
    monkeypatch.setattr(m.shutil, 'rmtree', partial)
    with pytest.raises(KeyboardInterrupt): m.release(source, target, record, journal, fstab=fstab)
    assert m.read(journal)['phase'] == 'releasing'
    monkeypatch.setattr(m.shutil, 'rmtree', real_delete)
    m.release(source, target, record, journal, fstab=fstab)
    assert record['phase'] == 'released' and (target / 'execution.lock').exists()


def fake_process(proc, pid, *, args=b'', env=b'', cwd=None, fd=None):
    p = proc / str(pid); (p / 'fd').mkdir(parents=True)
    (p / 'cmdline').write_bytes(args); (p / 'environ').write_bytes(env)
    if cwd: (p / 'cwd').symlink_to(cwd, target_is_directory=True)
    if fd: (p / 'fd/3').symlink_to(fd)


@pytest.mark.parametrize('where', ['args', 'env', 'cwd', 'fd'])
def test_busy_study_detected_without_disclosing_arguments_or_secrets(tmp_path, where):
    root = tmp_path.resolve(); proc = root / 'proc'; proc.mkdir()
    source = root / 'cohort_v1'; source.mkdir()
    kwargs = {'args': b'python\0--study\0' + os.fsencode(source)} if where == 'args' else {}
    if where == 'env': kwargs['env'] = b'SECRET=do-not-print\0TMPDIR=' + os.fsencode(source / '.runtime/tmp')
    if where in {'cwd', 'fd'}: kwargs[where] = source
    fake_process(proc, 777777, **kwargs)
    with pytest.raises(RuntimeError, match='777777') as e: m.process_users({source}, proc)
    assert 'SECRET' not in str(e.value) and str(source) not in str(e.value)


def test_unrelated_nckuh_and_similar_prefix_are_not_busy(tmp_path):
    proc = tmp_path.resolve() / 'proc'; proc.mkdir()
    source = tmp_path.resolve() / 'cohort_v1'
    fake_process(proc, 777777, args=b'java nextflow -work-dir /work/contexts/shared-nckuh\0')
    fake_process(proc, 777778, args=os.fsencode(str(source) + '0/file'))
    m.process_users({source}, proc)


def test_bind_alias_and_nested_mount_detection(study):
    source, _, _, _ = study
    # Samefile alias stands in for the kernel bind mapping in this unit test.
    physical = source.parent / 'physical'; physical.symlink_to(source, target_is_directory=True)
    rows = [dict(target=str(source), fsroot=str(physical), fstype='ext4')]
    assert physical in m.aliases(source.parents[2], source, rows)
    with pytest.raises(RuntimeError, match='子掛載'):
        m.no_nested_mounts({source}, [dict(target=str(source / 'cases/mount'))])


def test_mountinfo_reads_bind_roots_and_escaped_names(tmp_path):
    file = tmp_path / 'mountinfo'
    file.write_text('1 0 1:1 / / rw - ext4 /dev/root rw\n'
                    '2 1 1:1 /scratch/n102968/nf_work /home/user/nf_work rw shared:4 - ext4 /dev/root rw\n'
                    '3 1 1:2 / /has\\040space rw - nfs server:/export rw\n')
    rows = m.mount_rows(file)
    assert rows[1] == dict(target='/home/user/nf_work', fsroot='/scratch/n102968/nf_work', fstype='ext4')
    assert rows[2]['target'] == '/has space' and rows[2]['fstype'] == 'nfs'


def test_capacity_refuses_network_destination_before_any_copy(study):
    source, target, _, _ = study
    with pytest.raises(ValueError, match='本機'):
        m.capacity(source, target, [dict(target='/', fsroot='/', fstype='nfs')])


def main_environment(study, monkeypatch):
    source, target, journal, record = study
    commands = []; real_run = m.run
    fstab = target.parent / 'fstab'; fstab.write_text('# original\n')
    monkeypatch.setattr(m, 'FSTAB', fstab)
    monkeypatch.setattr(m.sys, 'platform', 'linux')
    monkeypatch.setattr(m.os, 'geteuid', lambda: 0)
    monkeypatch.setattr(m, 'shape', lambda *args: None)  # /home layout already tested separately.
    monkeypatch.setattr(m, 'capacity', lambda *args: None)  # No separate block devices in test.
    monkeypatch.setattr(m, 'mount_rows', lambda: [])
    monkeypatch.setattr(m, 'process_users', lambda *args: None)
    monkeypatch.setattr(m.os.path, 'ismount', lambda p: Path(p) == source and source.is_symlink())
    def run(argv, **kwargs):
        commands.append(list(map(str, argv)))
        if argv[0] == 'mount':
            source.rmdir(); source.symlink_to(target, target_is_directory=True); return ''
        if argv[0] in {'systemctl', 'df'}: return ''
        return real_run(argv, **kwargs)
    monkeypatch.setattr(m, 'run', run)
    args = ['--work-root', str(source.parents[2]), '--study', str(source), '--target', str(target)]
    return args, commands


def test_whole_driver_check_copy_release_and_idempotence(study, monkeypatch):
    source, target, journal, _ = study
    args, commands = main_environment(study, monkeypatch)
    original_plan = (source / 'plan.json').read_bytes()
    assert m.main(args) == 0 and not journal.exists() and not target.exists()
    assert m.main(args + ['--apply']) == 0
    assert (source / 'plan.json').read_bytes() == original_plan
    record = m.read(journal); assert record['phase'] == 'ready' and Path(record['backup']).exists()
    assert m.main(args + ['--apply']) == 0  # Doesn't copy/rename the mounted study a second time.
    assert m.main(args + ['--release-ssd']) == 0
    assert m.main(args + ['--release-ssd']) == 0
    assert m.read(journal)['phase'] == 'released' and not Path(record['backup']).exists()
    assert (source / 'plan.json').read_bytes() == original_plan
    assert len([c for c in commands if c[0] == 'mount']) == 1
    assert all(c[0] != 'systemctl' or c == ['systemctl', 'daemon-reload'] for c in commands)


def test_copy_interruption_can_resume_original_plan_and_partial_copy(study, monkeypatch):
    source, target, journal, _ = study
    args, _ = main_environment(study, monkeypatch); real_run = m.run
    def interrupted(argv, **kwargs):
        if argv[0] == 'rsync' and '--dry-run' not in argv:
            real_run(argv, **kwargs); raise KeyboardInterrupt()
        return real_run(argv, **kwargs)
    monkeypatch.setattr(m, 'run', interrupted)
    with pytest.raises(KeyboardInterrupt): m.main(args + ['--apply'])
    assert m.read(journal)['phase'] == 'copying' and source.is_dir() and not source.is_symlink()
    monkeypatch.setattr(m, 'run', real_run)
    assert m.main(args + ['--apply']) == 0 and m.read(journal)['phase'] == 'ready'


def test_capacity_failure_never_creates_target_or_journal(study, monkeypatch):
    source, target, journal, _ = study
    args, commands = main_environment(study, monkeypatch)
    monkeypatch.setattr(m, 'capacity', lambda *a: (_ for _ in ()).throw(RuntimeError('容量不足')))
    with pytest.raises(RuntimeError, match='容量不足'): m.main(args + ['--apply'])
    assert not target.exists() and not journal.exists() and not commands


def test_fstab_conflict_is_not_overwritten(study):
    source, target, _, _ = study
    old = f'/different/device {source} ext4 defaults 0 0\n'
    with pytest.raises(RuntimeError, match='fstab'): m.fstab_text(old, source, target)
