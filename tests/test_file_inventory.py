from concurrent.futures import ThreadPoolExecutor
import json
import sqlite3
from pathlib import Path
import threading
import time
from types import SimpleNamespace

from fastapi import FastAPI
from fastapi.testclient import TestClient
import pytest
from starlette.middleware.sessions import SessionMiddleware

from app import config
from app.auth import current_user
from app.routers import files, secondary as secondary_api, dragen as tertiary_api
from app.services import file_index, secondary_analysis as secondary, dragen_jobs as tertiary
from app.workers import file_index_update, secondary_index_update, tertiary_index_update


def touch(path, content=b"1234567890"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return path


@pytest.fixture
def tree(tmp_path, monkeypatch):
    raw, wes, output, published = [tmp_path/n for n in ("raw","wes","output","tertiary")]
    for p in (raw, wes, output, published): p.mkdir()
    monkeypatch.setattr(config, "FILE_INDEX_PATH", tmp_path/"state/index.sqlite")
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", published)
    monkeypatch.setattr(config, "LEGACY_PIPELINE_OUT_ROOT", tmp_path/"old")
    monkeypatch.setattr(config, "LEGACY_TERTIARY_OUTPUT_ROOT", tmp_path/"legacy")
    monkeypatch.setattr(secondary, "SECONDARY_WES_FASTQ_ROOTS", [wes])
    monkeypatch.setattr(secondary, "SECONDARY_WGS_FASTQ_ROOTS", [raw])
    monkeypatch.setattr(tertiary, "DRAGEN_VCF_ROOTS", [raw])
    monkeypatch.setattr(tertiary, "INHOUSE_VCF_ROOTS", [output])
    roots = dict(wes=[wes], wgs=[raw], dragen=[raw], inhouse=[output], dragen_bam=[raw], inhouse_bam=[output])
    monkeypatch.setattr(file_index, "data_roots", lambda:roots)
    return raw, wes, output, published


def wgs(raw, sid="S1", lanes=8):
    for lane in range(1,lanes+1):
        for read in (1,2):
            touch(raw/f"run/fastq.gz/{sid}_S50_L{lane:03d}_R{read}_001.fastq.gz")


def client(auth=True):
    app = FastAPI()
    app.add_middleware(SessionMiddleware, secret_key="synthetic-test-only")
    for router in (files.router, secondary_api.router, tertiary_api.router): app.include_router(router)
    if auth: app.dependency_overrides[current_user] = lambda:{"username":"test"}
    return TestClient(app)


def test_six_columns_and_analysis_views_share_one_generation(tree):
    raw, wes, output, published = tree
    wgs(raw)
    touch(raw/"run/bam/S1.bam"); touch(raw/"run/bam/S1.bam.bai")
    touch(raw/"run/bam/S1.repeats.bam"); touch(raw/"run/bam/S1.repeats.bam.bai")
    for suffix in ("hard-filtered.vcf.gz","cnv.vcf.gz","sv.vcf.gz","cnv_sv.vcf.gz","ploidy.vcf.gz"):
        touch(raw/f"run/vcf.gz/S1.{suffix}")
        touch(raw/f"run/vcf.gz/S1.{suffix}.tbi")
    touch(raw/"run/other/S1/germline_seq/S1.hard-filtered.gvcf.gz")
    touch(raw/"run/other/S1/germline_seq/S1.targeted.json")
    touch(raw/"run/vcf.gz/VCF_ONLY.hard-filtered.vcf.gz")
    touch(wes/"run/W1_S1_R1_001.fastq.gz"); touch(wes/"run/W1_S1_R2_001.fastq.gz")
    touch(wes/"run/ONLY_S2_R1_001.fastq.gz")
    touch(output/"batch/W1/04_snv_indel/W1.ensemble.fixed.vcf.gz")
    touch(output/"batch/W1/02_alignment/W1.aligned.sorted.bam")
    touch(output/"batch/W1/02_alignment/W1.aligned.sorted.bai")
    for name in ("W1.cnv_sv.vcf.gz","W1.gvcf.gz","W1.targeted.json"):
        touch(output/f"batch/W1/05_cnv_sv/{name}")
    for name in ("W1.gcnv.vcf.gz","W1.delly.vcf.gz"):
        touch(output/f"batch/W1/05_cnv_sv/{name}")
    touch(published/"S1-dragen/08_postprocessing/S1-dragen.pipeline_source.json",
          json.dumps({"source_sample_id":"S1"}).encode())
    views = file_index.refresh()
    data = file_index.list_samples()
    samples = {s["sample_id"]:s for s in data["samples"]}
    assert set(samples) == {"S1","W1","ONLY"}
    assert len(samples["S1"]["fastq"]) == 16
    assert len(samples["S1"]["bam"]) == 2
    assert len(samples["S1"]["other"]) == 5
    assert not any("ploidy" in f["name"] or "repeats" in f["name"] or f["name"].endswith(".tbi") for s in samples.values() for col in ("bam","vcf","other") for f in s[col])
    assert len(samples["W1"]["other"]) == 2
    assert len(samples["W1"]["bam"]) == 2
    assert samples["ONLY"]["bam"] == samples["ONLY"]["vcf"] == samples["ONLY"]["other"] == []
    assert samples["S1"]["tertiary"][0]["name"] == "S1-dragen"
    assert views["secondary"]["wgs"][0]["lane_count"] == 8
    assert len(views["secondary"]["wes"]) == 1  # incomplete pair stays inventory-only
    assert {r["sample_id"] for r in views["tertiary"]["dragen"]} == {"S1","VCF_ONLY"}
    assert secondary.load_index()["generation"] == tertiary.load_index()["generation"] == data["meta"]["generation"]
    with client() as c:
        assert c.get('/api/files?q=S1').json()["total"] == 1
        assert c.get('/api/files?q=%').json()["total"] == 0
        assert c.get('/api/secondary/fastqs').json()["wgs"][0]["fastq_file_count"] == 16
        assert len(c.get('/api/dragen/vcfs').json()["dragen"]) == 2


def test_reanalysis_unknown_names_and_test_fastqs_do_not_pollute_samples(tree):
    raw, wes, output, _ = tree
    touch(wes/"Reanalysis/RENAMED/ORIGINAL.R1.clean.fastq.gz")
    touch(wes/"Reanalysis/RENAMED/ORIGINAL.R2.clean.fastq.gz")
    touch(wes/"run/undecodable.fastq.gz")
    touch(output/"pgx_getrm/_links/REFERENCE_S1_R1_001.fastq.gz")
    (output/"pgx_getrm/_links/BROKEN_S2_R1_001.fastq.gz").symlink_to(output/"absent")
    touch(output/"batch/RENAMED/04_snv_indel/RENAMED.ensemble.fixed.vcf.gz")
    touch(output/"batch/RENAMED/02_alignment/RENAMED.aligned.sorted.bam")
    file_index.refresh()
    data = file_index.list_samples()
    assert data["total"] == 1
    row = data["samples"][0]
    assert row["sample_id"] == "RENAMED"
    assert len(row["bam"]) == len(row["vcf"]) == 1
    assert data["meta"]["unrecognized_fastq_count"] == 1


def test_overlap_refresh_coalesces_and_failed_scan_preserves_generation(tree, monkeypatch):
    raw, *_ = tree
    wgs(raw,lanes=1)
    original = file_index.build
    entered, release = threading.Event(), threading.Event()
    calls = []
    def slow(previous):
        calls.append(1); entered.set(); assert release.wait(5)
        return original(previous)
    monkeypatch.setattr(file_index,"build",slow)
    with ThreadPoolExecutor(2) as pool:
        a = pool.submit(file_index.refresh,force=True)
        assert entered.wait(5)
        requested = time.time_ns()
        b = pool.submit(file_index.refresh,force=True,requested_at_ns=requested)
        time.sleep(.05); release.set()
        assert a.result()["meta"]["generation"] == b.result()["meta"]["generation"]
    assert len(calls) == 1
    before = file_index.metadata()["generation"]
    monkeypatch.setattr(file_index,"build",original)
    raw.rename(raw.with_name("offline"))
    with pytest.raises(OSError): file_index.refresh(force=True)
    assert file_index.metadata()["generation"] == before
    assert file_index.list_samples()["total"] == 1


def test_download_range_auth_and_bounded_browsing(tree,tmp_path):
    raw, _, _, published = tree
    wgs(raw,lanes=1)
    secret = touch(tmp_path/"secret.txt",b"SECRET")
    touch(published/"S1-dragen/03_acmg/S1.tsv",b"abcdefghij")
    touch(published/"S1-dragen/03_acmg/<script>.txt",b"safe")
    (published/"S1-dragen/escape.txt").symlink_to(secret)
    file_index.refresh()
    item = file_index.list_samples()["samples"][0]["fastq"][0]
    link = '/api/files/download/'+item["id"]
    root = file_index.identity(published.resolve())
    with client(False) as c:
        assert c.get(link).status_code == 401
        assert c.get('/api/files/browse',params={"root":root}).status_code == 401
    with client() as c:
        assert c.get(link).content == b"1234567890"
        part=c.get(link,headers={"Range":"bytes=2-5"})
        assert part.status_code == 206 and part.content == b"3456"
        assert part.headers["content-range"] == "bytes 2-5/10"
        assert c.get(link,headers={"Range":"bytes=-3"}).content == b"890"
        assert c.head(link).headers["content-length"] == '10'
        assert c.head(link).content == b''
        assert c.get(link,headers={"Range":"bytes=20-25"}).status_code == 416
        page=c.get('/api/files/browse',params={"root":root,"path":"S1-dragen/03_acmg"})
        assert page.status_code == 200 and "上一層" in page.text
        assert '&lt;script&gt;.txt' in page.text and '<script>' not in page.text
        assert c.get('/api/files/browse',params={"root":root,"path":"../"}).status_code == 403
        assert c.get('/api/files/tertiary-download',params={"root":root,"path":"S1-dragen/escape.txt"}).status_code == 403
        response=c.get('/api/files/tertiary-download',params={"root":root,"path":"S1-dragen/03_acmg/S1.tsv"})
        assert response.content == b"abcdefghij"
        p=Path(item["path"]);p.unlink();p.symlink_to(secret)
        assert c.get(link).status_code == 403
        p.unlink()
        assert c.get(link).status_code == 404


def test_root_configuration_change_rebuilds_without_waiting_for_ttl(tree):
    raw, wes, *_ = tree
    wgs(raw,lanes=1)
    first = file_index.refresh()["meta"]["generation"]
    roots = file_index.data_roots()
    extra=wes/"extra";extra.mkdir()
    roots["wes"].append(extra)
    assert file_index.stale(file_index.metadata())
    assert file_index.refresh()["meta"]["generation"] != first


def test_legacy_daily_timers_and_view_requests_do_not_scan_again(tree, monkeypatch):
    raw, *_ = tree
    wgs(raw,lanes=1)
    file_index_update.main()
    generation=file_index.metadata()["generation"]
    def unexpected(*args,**kwargs):
        raise AssertionError("a second scan in the same daily window")
    monkeypatch.setattr(file_index,"build",unexpected)
    assert secondary_index_update.main() == tertiary_index_update.main() == 0
    assert secondary.refresh_index()["generation"] == tertiary.refresh_index()["generation"] == generation


def test_all_runs_are_preserved_and_first_two_fastqs_remain_a_pair(tree):
    raw, *_ = tree
    wgs(raw,lanes=1)
    touch(raw/"another_run/fastq.gz/S1_S50_L001_R1_001.fastq.gz")
    touch(raw/"another_run/fastq.gz/S1_S50_L001_R2_001.fastq.gz")
    touch(raw/"run/bam/S1.bam")
    touch(raw/"another_run/bam/S1.bam")
    file_index.refresh()
    sample=file_index.list_samples()["samples"][0]
    assert len(sample["fastq"]) == 4 and len(sample["bam"]) == 2
    assert Path(sample["fastq"][0]["path"]).parent == Path(sample["fastq"][1]["path"]).parent
    assert '_R1_' in sample["fastq"][0]["name"] and '_R2_' in sample["fastq"][1]["name"]


def test_merged_and_simple_fastq_names_are_inventory_only(tree):
    raw, *_ = tree
    for name in ('MERGED_R1_merged.fastq.gz','MERGED_R2_merged.fastq.gz','PLAIN.R1.fq.gz','PLAIN.R2.fq.gz'):
        touch(raw/"run/fastq.gz"/name)
    views=file_index.refresh()
    rows=file_index.list_samples()["samples"]
    assert {s["sample_id"] for s in rows} == {"MERGED","PLAIN"}
    assert all(len(s["fastq"])==2 for s in rows)
    assert views["secondary"]["wgs"] == []


def test_blank_query_never_builds_index_or_returns_all_samples(tree, monkeypatch):
    def unexpected(*args, **kwargs):
        raise AssertionError("empty query must not enumerate or build")
    monkeypatch.setattr(file_index, "refresh", unexpected)
    monkeypatch.setattr(file_index, "list_samples", unexpected)
    with client() as c:
        for query in ('', '   '):
            response=c.get('/api/files', params={'q':query})
            assert response.status_code == 200
            assert response.json()['samples'] == [] and response.json()['total'] == 0
    assert not config.FILE_INDEX_PATH.exists()


def test_existing_snapshot_gets_both_folder_dates_and_copy_paths_without_rebuild(tree, monkeypatch):
    raw, _, _, published = tree
    wgs(raw, lanes=1)
    folders=[published/'S1-dragen', published/'S1-nckuh']
    for folder in folders:
        touch(folder/'08_postprocessing'/f'{folder.name}.pipeline_source.json',json.dumps({'source_sample_id':'S1'}).encode())
    file_index.refresh()
    generation=file_index.metadata()['generation']
    monkeypatch.setattr(file_index,'_native_birth_time',lambda p:1762362000 if p==folders[0] else 1762189200)
    with client() as c:
        entries=c.get('/api/files',params={'q':'S1'}).json()['samples'][0]['tertiary']
    assert {e['absolute_path'] for e in entries} == {str(p) for p in folders}
    assert {e['created_date'] for e in entries} == {'2025-11-06','2025-11-04'}
    assert all(e['created_at'].endswith('+08:00') for e in entries)
    assert file_index.metadata()['generation'] == generation


def test_linux_birth_time_is_batched_and_unsupported_time_is_not_guessed(tmp_path,monkeypatch):
    a,b=tmp_path/'a',tmp_path/'b'
    a.mkdir();b.mkdir()
    monkeypatch.setattr(file_index,'_native_birth_time',lambda p:None)
    monkeypatch.setattr(file_index.sys,'platform','linux')
    calls=[]
    def run(args,**kwargs):
        calls.append(args)
        return SimpleNamespace(stdout=b'\0'.join((bytes(a),'1762362000'.encode(),bytes(b),b'0'))+b'\0',returncode=0)
    monkeypatch.setattr(file_index.subprocess,'run',run)
    assert file_index._birth_times([a,b]) == {a:1762362000}
    assert len(calls)==1 and calls[0][-2:] == [str(a),str(b)]
    assert r'--printf=%n\0%W\0' in calls[0]


def test_annotation_record_date_and_new_annotation_file_work_with_existing_snapshot(tree,monkeypatch):
    raw, _, _, published = tree
    wgs(raw,'26WE0142',lanes=1)
    folder=published/'26WE0142-nckuh'
    touch(folder/'08_postprocessing/26WE0142-nckuh.pipeline_source.json',json.dumps({
        'source_sample_id':'26WE0142','annotated_at':'2026-11-06 13:24:00'
    }).encode())
    file_index.refresh()
    generation=file_index.metadata()['generation']
    annotation=touch(folder/'03_acmg/26WE0142.snv_indel.acmg.tsv',b'CHROM\tPOS\nchr1\t123\n')
    monkeypatch.setattr(file_index,'_birth_times',lambda paths:{})
    with client() as c:
        data=c.get('/api/files',params={'q':'26WE0142'}).json()['samples'][0]
        entry=data['tertiary'][0]
        assert entry['created_date'] is None
        assert entry['display_date']=='2026-11-06' and entry['date_source']=='tertiary_annotation'
        assert entry['date_at']=='2026-11-06T13:24:00+08:00'
        result=data['vcf'][0]
        assert result['name']==annotation.name and result['path']==str(annotation)
        response=c.get('/api/files/tertiary-download',params={'root':result['tertiary_root'],'path':result['relative_path']})
        assert response.status_code==200 and response.content==annotation.read_bytes()
    assert file_index.metadata()['generation']==generation


def test_primary_root_only_hides_legacy_links_even_in_old_snapshot(tree):
    raw, _, _, published = tree
    wgs(raw,lanes=1)
    folders=[published/'S1-dragen',config.LEGACY_PIPELINE_OUT_ROOT/'S1-nckuh',config.LEGACY_TERTIARY_OUTPUT_ROOT/'S1-dragen']
    for folder in folders:
        touch(folder/'08_postprocessing'/f'{folder.name}.pipeline_source.json',json.dumps({'source_sample_id':'S1'}).encode())
        touch(folder/'03_acmg/S1.snv_indel.acmg.tsv')
    file_index.refresh()
    data=file_index.list_samples()['samples'][0]
    assert len(data['tertiary'])==1 and len(data['vcf'])==1
    assert data['vcf'][0]['path'].startswith(str(published))
    # Simulate a v9.99 generation with cached links to both retired roots.
    with sqlite3.connect(config.FILE_INDEX_PATH) as db:
        stored=json.loads(db.execute('SELECT payload FROM samples').fetchone()[0])
        for folder in folders[1:]:
            stored['tertiary'].append(dict(name=folder.name,path=folder.name,root=file_index.identity(folder.parent)))
        db.execute('UPDATE samples SET payload=?',(json.dumps(stored),));db.commit()
    with client() as c:
        data=c.get('/api/files',params={'q':'S1'}).json()['samples'][0]
        assert len(data['tertiary'])==1
        assert data['tertiary'][0]['absolute_path']==str(folders[0])
        for folder in folders[1:]:
            assert c.get('/api/files/browse',params={'root':file_index.identity(folder.parent),'path':folder.name}).status_code==403
            assert c.get('/api/files/tertiary-download',params={'root':file_index.identity(folder.parent),'path':folder.name+'/03_acmg/S1.snv_indel.acmg.tsv'}).status_code==403


def test_annotation_filename_is_exact_and_outside_symlinks_are_rejected(tree,tmp_path):
    _,_,_,published=tree
    folder=published/'26WE0142-nckuh';folder.mkdir()
    touch(folder/'03_acmg/26WE01420.snv_indel.acmg.tsv')
    assert file_index._annotation_entry(published,folder,'26WE0142') is None
    secret=touch(tmp_path/'secret.tsv')
    (folder/'03_acmg/26WE0142.snv_indel.acmg.tsv').symlink_to(secret)
    assert file_index._annotation_entry(published,folder,'26WE0142') is None
