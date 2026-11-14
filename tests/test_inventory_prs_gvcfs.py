"""Inventory counts must distinguish candidates, repeats, controls and failures."""
import csv
import json
from pathlib import Path
import subprocess
import sys

import pytest

from scripts import inventory_prs_gvcfs as inventory


def gvcf(root, sid, run='run_a', filename=None, data=b'not-read-by-inventory'):
    p = root/run/'other'/sid/'germline_seq'/(filename or sid+'.hard-filtered.gvcf.gz')
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(data)
    return p


def test_counts_ids_separately_from_repeat_files_and_filters_other_samples(tmp_path):
    root = tmp_path/'raw';a = gvcf(root,'26T00001');gvcf(root,'26T00001',run='run_b')
    gvcf(root,'26T00002',filename='26T00002.gvcf.gz')
    gvcf(root,'26T10001');gvcf(root,'25T00001');gvcf(root,'VAL-57')
    gvcf(root,'26T00003',filename='26T00003.vcf.gz')
    before = (a.stat().st_mtime_ns, a.read_bytes())
    r = inventory.scan(root)
    assert r['discovery_complete']
    assert r['counts']['healthy_sample_ids'] == 2
    assert r['counts']['healthy_gvcf_paths'] == 3
    assert r['counts']['healthy_ids_with_multiple_files'] == 1
    assert r['samples'][0]['selection'] == 'MULTIPLE_FILES_NEED_REVIEW'
    assert (a.stat().st_mtime_ns,a.read_bytes()) == before
    assert all(not x['header_checked'] for x in r['files'])


def test_ten_val_runs_are_two_biological_standards(tmp_path):
    for n in range(45,55):gvcf(tmp_path,'VAL-'+str(n))
    r = inventory.scan(tmp_path)
    assert r['counts']['healthy_sample_ids'] == 0
    assert r['counts']['standard_sample_ids'] == 10
    assert r['counts']['standard_biological_ids'] == 2
    assert r['counts']['standard_gvcf_paths'] == 10
    assert r['expected_val_ids_not_found'] == []


def test_filename_folder_conflict_is_visible_and_not_usable(tmp_path):
    gvcf(tmp_path,'26T00001',filename='26T00002.gvcf.gz')
    r = inventory.scan(tmp_path)
    assert r['counts']['healthy_sample_ids'] == 1
    assert r['counts']['identity_conflicts'] == 1
    assert r['counts']['healthy_sample_ids_with_nonempty_matching_files'] == 0
    assert r['samples'][0]['selection'] == 'CHECK_IDENTITY'
    assert not r['files'][0]['ready_for_header_check']


def test_empty_files_and_missing_index_are_not_header_ready(tmp_path):
    gvcf(tmp_path,'26T00001',data=b'')
    good = gvcf(tmp_path,'26T00002')
    r = inventory.scan(tmp_path)
    assert r['counts']['healthy_sample_ids'] == 2
    assert r['counts']['healthy_sample_ids_with_nonempty_matching_files'] == 1
    assert r['counts']['paths_ready_for_header_check'] == 0
    Path(str(good)+'.tbi').write_bytes(b'index')
    assert inventory.scan(tmp_path)['counts']['paths_ready_for_header_check'] == 1


def test_directory_symlinks_are_not_followed_and_discovery_is_partial(tmp_path):
    raw = tmp_path/'raw';gvcf(raw,'26T00001')
    outside = tmp_path/'elsewhere';gvcf(outside,'26T00002')
    (raw/'linked_run').symlink_to(outside,target_is_directory=True)
    r = inventory.scan(raw)
    assert not r['discovery_complete']
    assert r['counts']['healthy_sample_ids'] == 1
    assert r['counts']['skipped_directory_symlinks'] == 1


def test_file_symlink_aliases_and_cross_id_aliases_are_reported(tmp_path):
    p = gvcf(tmp_path,'26T00001')
    alias = tmp_path/'run_b'/'other'/'26T00001'/'germline_seq'/p.name
    alias.parent.mkdir(parents=True);alias.symlink_to(p)
    r = inventory.scan(tmp_path)
    assert r['counts']['healthy_gvcf_paths'] == 2
    assert r['samples'][0]['distinct_physical_files'] == 1
    assert r['samples'][0]['selection'] == 'SAME_FILE_MULTIPLE_PATHS'
    cross = tmp_path/'run_c'/'other'/'26T00002'/'germline_seq'/'26T00002.gvcf.gz'
    cross.parent.mkdir(parents=True);cross.symlink_to(p)
    r = inventory.scan(tmp_path)
    assert r['counts']['physical_files_shared_across_ids'] == 1
    assert all(s['selection']=='CHECK_IDENTITY' for s in r['samples'])


def test_io_errors_produce_partial_counts_instead_of_false_zero_success(tmp_path,monkeypatch):
    root = tmp_path/'raw';root.mkdir()
    def denied(path,**kwargs):
        kwargs['onerror'](PermissionError(13,'permission denied',str(path/'run_a')))
        return iter(())
    monkeypatch.setattr(inventory.os,'walk',denied)
    r = inventory.scan(root)
    assert r['counts']['healthy_sample_ids'] == 0
    assert not r['discovery_complete'] and r['counts']['errors'] == 1
    assert inventory.main(['--root',str(root),'--output-dir',str(tmp_path/'reports')]) == 2


def test_reports_do_not_overwrite_and_cannot_be_written_under_raw_root(tmp_path):
    root=tmp_path/'raw';gvcf(root,'26T00001')
    r=inventory.scan(root)
    a=inventory.write_reports(r,tmp_path/'reports');b=inventory.write_reports(r,tmp_path/'reports')
    assert a != b
    assert json.loads((a/'inventory.json').read_text())['counts'] == r['counts']
    with (a/'files.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==1 and rows[0]['sample_id']=='26T00001'
    with pytest.raises(ValueError,match='原始資料'):
        inventory.write_reports(r,root/'reports')
    assert not (root/'reports').exists()


def test_generic_and_uncompressed_filenames_preserve_eligible_folder_identity(tmp_path):
    gvcf(tmp_path,'26T00001',filename='unidentified.gvcf.gz')
    p=tmp_path/'run_b'/'26T00002.g.vcf';p.parent.mkdir();p.write_bytes(b'not read')
    r=inventory.scan(tmp_path)
    assert r['counts']['healthy_sample_ids']==2
    assert r['counts']['identity_conflicts']==1
    assert any(f['identity_basis']=='filename_only' for f in r['files'])


def test_broken_file_link_makes_discovery_partial(tmp_path):
    gvcf(tmp_path,'26T00001')
    broken=tmp_path/'run_b'/'26T00002.gvcf.gz'
    broken.parent.mkdir();broken.symlink_to(tmp_path/'not_present')
    r=inventory.scan(tmp_path)
    assert not r['discovery_complete'] and r['counts']['errors']==1
    assert r['counts']['healthy_sample_ids']==2
    assert r['counts']['healthy_sample_ids_with_nonempty_matching_files']==1


def test_standalone_cli_writes_reports_and_reports_missing_root_as_failure(tmp_path):
    root=tmp_path/'raw';gvcf(root,'26T00001');gvcf(root,'VAL-45')
    script=Path(inventory.__file__).resolve()
    result=subprocess.run([sys.executable,str(script),'--root',str(root),
        '--output-dir',str(tmp_path/'reports')],capture_output=True,text=True)
    assert result.returncode==0 and '26T0 健康人候選檢體 ID：1' in result.stdout
    reports=list((tmp_path/'reports').glob('inventory_*/inventory.json'))
    assert len(reports)==1
    assert json.loads(reports[0].read_text())['counts']['standard_biological_ids']==1
    result=subprocess.run([sys.executable,str(script),'--root',str(tmp_path/'absent'),
        '--output-dir',str(tmp_path/'unused')],capture_output=True,text=True)
    assert result.returncode==1 and 'ERROR' in result.stderr
    assert not (tmp_path/'unused').exists()
