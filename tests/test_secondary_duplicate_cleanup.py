import importlib.util
import json
import os
from contextlib import nullcontext
from datetime import datetime
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location("secondary_cleanup", Path(__file__).resolve().parents[1]/"scripts/cleanup_secondary_duplicates.py")
cleanup = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cleanup)
CUTOFF = datetime(2026,11,5,tzinfo=cleanup.TAIPEI).timestamp()


@pytest.fixture
def tree(tmp_path,monkeypatch):
    root=tmp_path/'nextflow_output';root.mkdir()
    tertiary=tmp_path/'tertiary_output';tertiary.mkdir()
    dates={}
    def add(batch, sid, day, complete=True):
        folder=root/batch/sid
        for rel in ('02_alignment/'+sid+'.aligned.sorted.bam','02_alignment/'+sid+'.aligned.sorted.bam.bai',
                    '04_snv_indel/'+sid+'.ensemble.fixed.vcf.gz','04_snv_indel/'+sid+'.ensemble.fixed.vcf.gz.tbi'):
            path=folder/rel;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(b'SYNTHETIC')
        if not complete:(folder/'04_snv_indel'/f'{sid}.ensemble.fixed.vcf.gz.tbi').unlink()
        dates[str(folder)]=datetime(2026,11,day,tzinfo=cleanup.TAIPEI).timestamp()
        os.utime(folder,(dates[str(folder)],dates[str(folder)]))
        return folder
    return root,tertiary,add,dates


def actions(plan,sid='S1'):
    return {Path(r['path']).parent.name:r for r in plan['groups'][sid]}


def test_keep_latest_and_all_november_5_onward_versions(tree):
    root,tertiary,add,_=tree
    for batch,day in [('old_WES',1),('new_WES',4),('protected_WES',5),('latest_WES',6)]:add(batch,'S1',day)
    add('only_WES','S10',1)
    plan=cleanup.audit(root,tertiary,CUTOFF)
    rows=actions(plan)
    assert rows['old_WES']['action']==rows['new_WES']['action']=='DELETE'
    assert rows['protected_WES']['action']==rows['latest_WES']['action']=='KEEP'
    assert 'S10' not in plan['groups']


def test_equal_latest_dates_are_both_retained(tree):
    root,tertiary,add,_=tree
    add('a_WES','S1',1);add('b_WES','S1',4);add('c_WES','S1',4)
    rows=actions(cleanup.audit(root,tertiary,CUTOFF))
    assert rows['a_WES']['action']=='DELETE'
    assert rows['b_WES']['action']==rows['c_WES']['action']=='KEEP'


@pytest.mark.parametrize('condition',['incomplete_latest','wes_wgs'])
def test_uncertain_groups_are_not_deleted(tree,condition):
    root,tertiary,add,_=tree
    add('old_WES','S1',1)
    add('new_WGS' if condition=='wes_wgs' else 'new_WES','S1',4,complete=condition!='incomplete_latest')
    assert all(r['action']=='KEEP' for r in actions(cleanup.audit(root,tertiary,CUTOFF)).values())


def test_sample_folder_mtime_at_cutoff_is_protected(tree):
    root,tertiary,add,_=tree
    old=add('old_WES','S1',1);add('new_WES','S1',6)
    os.utime(old,(CUTOFF,CUTOFF))
    row=actions(cleanup.audit(root,tertiary,CUTOFF))['old_WES']
    assert row['action']=='KEEP' and '2026/11/5' in row['reason']


def test_no_creation_time_lookup_and_ctime_does_not_block_mtime_policy(tree,monkeypatch):
    root,tertiary,add,_=tree
    add('old_WES','S1',1);add('new_WES','S1',4)
    def unexpected(*args,**kwargs):
        raise AssertionError('audit must not call stat for creation time')
    monkeypatch.setattr(cleanup.subprocess,'run',unexpected)
    original=cleanup.signature
    def recent_ctime(path):
        result=original(path)
        result[5]=int((CUTOFF+86400)*1e9)
        return result
    monkeypatch.setattr(cleanup,'signature',recent_ctime)
    plan=cleanup.audit(root,tertiary,CUTOFF)
    assert plan['date_basis']=='sample_directory_mtime'
    assert actions(plan)['old_WES']['action']=='DELETE'


def test_child_file_dates_do_not_replace_sample_folder_mtime(tree):
    root,tertiary,add,_=tree
    old=add('old_WES','S1',1);add('new_WES','S1',4)
    child=old/'02_alignment/S1.aligned.sorted.bam'
    os.utime(child,(CUTOFF+86400,CUTOFF+86400))
    assert old.stat().st_mtime<CUTOFF
    assert actions(cleanup.audit(root,tertiary,CUTOFF))['old_WES']['action']=='DELETE'


@pytest.mark.parametrize('reference',['symlink','source_marker'])
def test_referenced_old_outputs_are_retained(tree,reference):
    root,tertiary,add,_=tree
    old=add('old_WES','S1',1);new=add('new_WES','S1',4)
    if reference=='symlink':
        (new/'old_source').symlink_to(old/'04_snv_indel',target_is_directory=True)
    else:
        marker=tertiary/'S1-nckuh/08_postprocessing/S1-nckuh.pipeline_source.json'
        marker.parent.mkdir(parents=True);marker.write_text(json.dumps({'source_vcf_path':str(old/'04_snv_indel/S1.ensemble.fixed.vcf.gz')}))
    row=actions(cleanup.audit(root,tertiary,CUTOFF))['old_WES']
    assert row['action']=='KEEP' and row['references']


@pytest.mark.parametrize('mount_index',[0,1],ids=['dgm','dgx'])
@pytest.mark.parametrize('foreign',[False,True],ids=['same-host-reference','other-host-reference'])
@pytest.mark.parametrize('reference',['symlink','source_path','source_vcf_path'])
def test_references_protect_old_output_on_either_shared_mount(tree,tmp_path,monkeypatch,mount_index,foreign,reference):
    root,tertiary,add,_=tree
    old=add('old_WES','S1',1);add('new_WES','S1',4)
    mounts=(tmp_path/'home/datalake_Intermediate',tmp_path/'datalake_Intermediate')
    monkeypatch.setattr(cleanup,'SHARED_MOUNTS',mounts)
    local=mounts[mount_index]
    new_root=local/'pipeline/nextflow_output'
    new_root.parent.mkdir(parents=True)
    root.rename(new_root)
    new_tertiary=local/'pipeline/tertiary_output'
    tertiary.rename(new_tertiary)
    target_mount=mounts[1-mount_index] if foreign else local
    target=target_mount/'pipeline/nextflow_output'/old.relative_to(root)/'04_snv_indel/S1.ensemble.fixed.vcf.gz'
    state=new_tertiary/'S1-nckuh/08_postprocessing'
    state.mkdir(parents=True)
    if reference=='symlink':
        # The other host's absolute path is deliberately unmounted here.
        (state/'original.vcf.gz').symlink_to(target)
    else:
        (state/'S1-nckuh.pipeline_source.json').write_text(json.dumps({reference:str(target)}))
    plan=cleanup.audit(new_root,new_tertiary,CUTOFF)
    row=actions(plan)['old_WES']
    assert not plan['dependency_error']
    assert row['action']=='KEEP' and row['reason']=='仍被其他資料引用'


def test_similar_mount_name_does_not_falsely_protect_unreferenced_old_output(tree,tmp_path,monkeypatch):
    root,tertiary,add,_=tree
    add('old_WES','S1',1);add('new_WES','S1',4)
    monkeypatch.setattr(cleanup,'SHARED_MOUNTS',(tmp_path,tmp_path.parent/'other-mount'))
    state=tertiary/'S1-nckuh/08_postprocessing'
    state.mkdir(parents=True)
    target=tmp_path.parent/'other-mount-extra/nextflow_output/old_WES/S1'
    (state/'S1-nckuh.pipeline_source.json').write_text(json.dumps({'source_path':str(target)}))
    assert actions(cleanup.audit(root,tertiary,CUTOFF))['old_WES']['action']=='DELETE'


def test_permission_denied_still_prevents_deletion(tree,monkeypatch):
    root,tertiary,add,_=tree
    add('old_WES','S1',1);add('new_WES','S1',4)
    walk=cleanup.os.walk
    def unreadable(top,**kwargs):
        if Path(top)==tertiary:
            kwargs['onerror'](PermissionError(13,'Permission denied',str(tertiary/'case/08_postprocessing/case.prs/runs/job/ancestry')))
        yield from walk(top,**kwargs)
    monkeypatch.setattr(cleanup.os,'walk',unreadable)
    plan=cleanup.audit(root,tertiary,CUTOFF)
    assert 'Permission denied' in plan['dependency_error']
    assert all(row['action']=='KEEP' for row in plan['groups']['S1'])


def test_report_only_dgm_invocation_never_prompts_or_deletes(tree,tmp_path,monkeypatch,capsys):
    root,tertiary,add,_=tree
    old=add('old_WES','S1',1);new=add('new_WES','S1',4)
    report=tmp_path/'report-dgm'
    monkeypatch.setattr(cleanup.sys,'argv',['cleanup','--root',str(root),'--reference-root',str(tertiary),'--report-dir',str(report)])
    def unexpected(*args,**kwargs):
        raise AssertionError('report-only entry must never ask for deletion or delete')
    monkeypatch.setattr(cleanup,'delete_reviewed',unexpected)
    monkeypatch.setattr('builtins.open',unexpected)
    assert cleanup.main()==0
    plan=json.loads((report/'plan.json').read_text())
    assert actions(plan)['old_WES']['action']=='DELETE'
    assert old.exists() and new.exists() and not (report/'actions.jsonl').exists()
    output=capsys.readouterr().out
    assert str(tertiary) in output and '本次只有盤點，沒有刪除' in output


def test_changed_retained_latest_aborts_before_deletion(tree,tmp_path,monkeypatch):
    root,tertiary,add,_=tree
    old=add('old_WES','S1',1);new=add('new_WES','S1',4)
    plan=cleanup.audit(root,tertiary,CUTOFF)
    (new/'02_alignment/S1.aligned.sorted.bam').write_bytes(b'CHANGED')
    monkeypatch.setattr(cleanup,'idle_guard',lambda p:nullcontext())
    with pytest.raises(RuntimeError,match='狀態已變動'):
        cleanup.delete_reviewed(plan,tmp_path,tmp_path)
    assert old.exists() and new.exists()


def test_delete_only_reviewed_sample_and_preserve_external_link_targets(tree,tmp_path,monkeypatch):
    root,tertiary,add,dates=tree
    old=add('old_WES','S1',1);new=add('new_WES','S1',4)
    other=add('old_WES','OTHER',1)
    external=tmp_path/'original.fastq.gz';external.write_bytes(b'RAW INPUT')
    (old/'external-link').symlink_to(external)
    os.utime(old,(dates[str(old)],dates[str(old)]))
    sheet=old.parent/'samplesheet.csv';sheet.write_text('sample\nS1\nOTHER\n')
    plan=cleanup.audit(root,tertiary,CUTOFF)
    monkeypatch.setattr(cleanup,'idle_guard',lambda p:nullcontext())
    monkeypatch.setattr(cleanup,'no_nextflow',lambda:None)
    monkeypatch.setattr(cleanup,'nested_mounts',lambda p:None)
    cleanup.delete_reviewed(plan,tmp_path,tmp_path)
    assert not old.exists() and new.exists() and other.exists() and sheet.exists()
    assert external.read_bytes()==b'RAW INPUT'
    events=[json.loads(line) for line in (tmp_path/'actions.jsonl').read_text().splitlines()]
    assert [e['status'] for e in events]==['staging','deleting','deleted']


def test_enter_does_not_delete(tree,tmp_path,monkeypatch,capsys):
    root,_,add,_=tree
    old=add('old_WES','S1',1);add('new_WES','S1',4)
    class Terminal:
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def write(self,*args):pass
        def flush(self):pass
        def readline(self):return '\n'
    real_open=open
    monkeypatch.setattr('builtins.open',lambda p,*args,**kwargs:Terminal() if str(p)=='/dev/tty' else real_open(p,*args,**kwargs))
    monkeypatch.setattr(cleanup.sys,'argv',['cleanup','--root',str(root),'--review','--report-dir',str(tmp_path/'report')])
    assert cleanup.main()==0
    assert old.exists() and '未確認，沒有刪除' in capsys.readouterr().out


def test_busy_platform_lease_blocks_cleanup(tmp_path,monkeypatch):
    work=tmp_path/'work';work.mkdir()
    locks=tmp_path/'.work.ngs-ui-locks';locks.mkdir()
    (locks/'admission.lock').touch()
    lease=locks/'scope.lock';lease.touch()
    monkeypatch.setattr(cleanup,'no_nextflow',lambda:None)
    with lease.open('a') as handle:
        cleanup.fcntl.flock(handle,cleanup.fcntl.LOCK_SH)
        with pytest.raises(BlockingIOError):
            with cleanup.idle_guard(work):raise AssertionError('must not enter')


def test_missing_reference_mount_still_reports_duplicates_but_deletes_nothing(tree):
    root,tertiary,add,_=tree
    add('old_WES','S1',1);add('new_WES','S1',4)
    tertiary.rmdir()
    plan=cleanup.audit(root,tertiary,CUTOFF)
    assert plan['dependency_error']
    assert all(r['action']=='KEEP' for r in plan['groups']['S1'])


@pytest.mark.parametrize('filename,expected_root,review',[
    ('paste_cleanup_secondary_duplicates.sh','/datalake_Intermediate/pipeline/nextflow_output',True),
    ('paste_cleanup_secondary_duplicates_dgm.sh','/home/datalake_Intermediate/pipeline/nextflow_output',False),
])
def test_paste_bundle_contains_the_tested_current_tool(filename,expected_root,review):
    root=Path(__file__).resolve().parents[1]
    source=(root/'scripts/cleanup_secondary_duplicates.py').read_text()
    bundle=(root/'scripts'/filename).read_text()
    assert "<<'NGS_SECONDARY_REVIEW_PY'\n"+source+'\nNGS_SECONDARY_REVIEW_PY\n' in bundle
    invocation=bundle.split('\nNGS_SECONDARY_REVIEW_PY\n',1)[1]
    assert '--root '+expected_root in invocation
    assert ('--review' in invocation)==review
    if not review:
        assert '--reference-root /home/datalake_Intermediate/pipeline/tertiary_output' in invocation
