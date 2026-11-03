"""Migration safeguards, without a real mount or changes to /etc."""
import importlib.util
import json
from pathlib import Path
import subprocess

import pytest


spec = importlib.util.spec_from_file_location('scratch_migration', Path(__file__).parents[1]/'scripts/migrate_nf_work_to_scratch.py')
migration = importlib.util.module_from_spec(spec)
spec.loader.exec_module(migration)


@pytest.mark.parametrize('problem', ['service', 'queued', 'running', 'manual_nextflow'])
def test_migration_rejects_active_work(tmp_path, monkeypatch, problem):
    root = tmp_path/'data'
    job = root/'jobs/prs/job'; job.mkdir(parents=True)
    (job/'state.json').write_text(json.dumps({'state': problem if problem in {'queued', 'running'} else 'done'}))
    monkeypatch.setattr(migration, 'run', lambda *a, **kw: subprocess.CompletedProcess([], 0, stdout='active' if problem == 'service' else 'inactive'))
    monkeypatch.setattr(migration.nextflow_work, 'nextflow_process_scopes', lambda root: (set(), ['123'] if problem == 'manual_nextflow' else []))
    with pytest.raises(RuntimeError): migration.check_idle(root, 'ngs-ui', tmp_path/'nf_work')


def test_mount_failure_restores_source_and_leaves_fstab_intact(tmp_path, monkeypatch):
    source = tmp_path/'home/nf_work'; source.mkdir(parents=True)
    (source/'.nextflow').write_text('session metadata')
    target = tmp_path/'scratch/nf_work'; target.mkdir(parents=True)
    fstab = tmp_path/'fstab'; fstab.write_text('# original\n')
    monkeypatch.setattr(migration, 'run', lambda *a, **kw: (_ for _ in ()).throw(OSError('mount failed')))
    with pytest.raises(OSError, match='mount failed'):
        migration.install_mount(source, target, 'ngs-ui', 'stamp', fstab=fstab, unit_root=tmp_path/'systemd')
    assert (source/'.nextflow').read_text() == 'session metadata'
    assert fstab.read_text() == '# original\n'
    assert not source.with_name('nf_work.hdd-backup-stamp').exists()


def test_success_installs_bind_dependency_but_never_starts_service(tmp_path, monkeypatch):
    source = tmp_path/'home/nf_work'; source.mkdir(parents=True)
    (source/'.nextflow').write_text('session metadata')
    target = tmp_path/'scratch/nf_work'; target.mkdir(parents=True)
    fstab = tmp_path/'fstab'; fstab.write_text('# original\n')
    commands = []
    monkeypatch.setattr(migration, 'run', lambda argv, **kw: commands.append(argv))
    samefile = migration.os.path.samefile
    monkeypatch.setattr(migration.os.path, 'samefile', lambda a, b: True if (a, b) == (source, target) else samefile(a, b))
    backup = migration.install_mount(source, target, 'ngs-ui', 'stamp', fstab=fstab, unit_root=tmp_path/'systemd')
    assert (backup/'.nextflow').read_text() == 'session metadata'
    assert f'{target} {source} none bind,nofail 0 0' in fstab.read_text()
    assert (tmp_path/'fstab.ngs-scratch-stamp').read_text() == '# original\n'
    settings = (tmp_path/'systemd/ngs-ui.service.d/50-ngs-scratch.conf').read_text()
    assert f'RequiresMountsFor={source}' in settings
    assert f'ExecStartPre=/usr/bin/mountpoint --quiet {source}' in settings
    assert "os.stat('/').st_dev != os.stat(" in settings
    assert commands == [['mount', '--bind', str(target), str(source)], ['systemctl', 'daemon-reload']]


def test_checksum_verification_is_read_only_and_keeps_links(tmp_path):
    args = migration.rsync_command(tmp_path/'source', tmp_path/'target', verify=True)
    assert '--checksum' in args and '--dry-run' in args and '--itemize-changes' in args
    assert '-aHAXS' in args  # Includes hardlinks, sparse files, ACLs and xattrs.
