"""Regression for the real ScoreHeader failure and immutable installed models."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from research.prs.multidisease import ancestry_input as adapter
from research.prs.multidisease.common import sha
from research.prs.multidisease import run as pipeline
from app.services import prs_presentation


@pytest.fixture
def legacy_score(tmp_path):
    score = tmp_path/'model/ancestry_score.txt.gz'; score.parent.mkdir()
    body = '1\t10\tG\tA\t0.000123\n2\t20\tC\tT\t-0.0015\n'
    with gzip.open(score, 'wt') as f:
        f.write('#pgs_id=PGS004879\n#genome_build=GRCh38\n#variants_number=2\n'+adapter.COLUMNS+body)
    return score, body


def test_old_model_rows_and_checksum_unchanged(legacy_score, tmp_path):
    source, body = legacy_score; before = source.read_bytes()
    target = adapter.stage_score(source, tmp_path/'ancestry_score.txt.gz', 2)
    assert source.read_bytes() == before
    with gzip.open(target, 'rt') as f:
        text = f.read()
    assert text.split(adapter.COLUMNS)[1] == body
    assert '#trait_reported=Coronary heart disease\n' in text
    provenance = json.loads(target.with_suffix('.source.json').read_text())
    assert provenance['source_sha256'] == sha(source)
    assert provenance['staged_sha256'] == sha(target)


@pytest.mark.parametrize('changed', ['#genome_build=GRCh37', '#pgs_id=PGS000001', '#variants_number=3'])
def test_identity_or_count_mismatch_is_not_repaired(legacy_score, tmp_path, changed):
    source, body = legacy_score
    with gzip.open(source, 'rt') as f: text = f.read()
    prefix = changed.split('=')[0]
    with gzip.open(source, 'wt') as f:
        f.write('\n'.join(changed if line.startswith(prefix+'=') else line for line in text.split('\n')))
    with pytest.raises(ValueError, match='mismatch'):
        adapter.stage_score(source, tmp_path/'staged.gz', 2)
    assert not (tmp_path/'staged.gz').exists()


def test_preflight_failure_never_starts_target_or_nextflow(legacy_score, tmp_path, monkeypatch):
    source, _ = legacy_score
    (tmp_path/'ancestry').mkdir()
    # Isolate only the release-size check for this two-row synthetic file.
    real = adapter.stage_score
    monkeypatch.setattr(pipeline, 'stage_score', lambda src, dst, n: real(src, dst, 2))
    def failed(*_): raise RuntimeError('synthetic formatter failure')
    monkeypatch.setattr(pipeline, 'preflight', failed)
    calls=[]; monkeypatch.setattr(pipeline, 'command', lambda *a, **k: calls.append(a))
    with pytest.raises(RuntimeError, match='formatter'):
        pipeline.ancestry({}, 'SYNTHETIC', tmp_path, source.parent)
    assert calls == []


def test_preflight_uses_existing_image_and_cleans_temporary_output(legacy_score, tmp_path, monkeypatch):
    source, _ = legacy_score
    target = adapter.stage_score(source, tmp_path/'staged.txt.gz', 2)
    image=tmp_path/'formatter.img'; image.touch(); calls=[]
    monkeypatch.setattr(adapter, 'command', lambda args: calls.append(list(map(str,args))))
    adapter.preflight(dict(pgscatalog_utils_img=str(image),plink_img='/unused/plink.img'), target)
    assert calls[0][:6]==['apptainer','exec','--bind','/home:/home',str(image),'pgscatalog-format']
    assert str(target) in calls[0] and not list(tmp_path.glob('.format-check-*'))
    assert json.loads((tmp_path/'format_check.json').read_text())['score_sha256']==sha(target)


def test_ancestry_nextflow_uses_case_input_and_keeps_model_unchanged(legacy_score, tmp_path, monkeypatch):
    source, _ = legacy_score; before=source.read_bytes(); calls=[]
    real=adapter.stage_score
    monkeypatch.setattr(pipeline,'stage_score',lambda src,dst,n:real(src,dst,2))
    monkeypatch.setattr(pipeline,'preflight',lambda cfg,score:calls.append('preflight'))
    def command(args, **kwargs):
        calls.append(args[0])
        if args[0]=='nextflow':
            staged=Path(args[args.index('--scorefile')+1])
            assert staged==tmp_path/'ancestry/ancestry_score.txt.gz'
            dest=tmp_path/'ancestry/results/caseSYNTHETIC/score'; dest.mkdir(parents=True)
            with gzip.open(dest/'caseSYNTHETIC_popsimilarity.txt.gz','wt') as f:
                f.write('sampleset\tIID\tMostSimilarPop\tMostSimilarPop_LowConfidence\ncaseSYNTHETIC\tSYNTHETIC\tEAS\tFalse\n')
            work=tmp_path/'ancestry/work/task'; work.mkdir(parents=True)
            with gzip.open(work/'ref_thinned.prune.in.gz','wt') as f:f.write('marker1\nmarker2\n')
    monkeypatch.setattr(pipeline,'command',command)
    cfg=dict(plink_img='/unused',ancestry_nextflow_config='/cfg',reference_panel='/ref',
        projection_method='oadp',reference_label='SuperPop',pca_maf_target=0,
        pca_geno_miss_target=.1,max_cpus=2,max_memory='8.GB')
    pipeline.ancestry(cfg,'SYNTHETIC',tmp_path,source.parent)
    assert calls==['preflight','apptainer','nextflow']
    assert source.read_bytes()==before
    assert json.loads((tmp_path/'qc.json').read_text())['pca_ld_pruned_markers']==2


def test_actual_pgscatalog_formatter(legacy_score, tmp_path):
    # Optional external primary-tool verification, also exercised manually in CI.
    packages = os.environ.get('NGS_UI_TEST_PGS_FORMAT_PACKAGES')
    if not packages: pytest.skip('Set NGS_UI_TEST_PGS_FORMAT_PACKAGES to pgscatalog.core installation')
    source, _ = legacy_score
    target = adapter.stage_score(source, tmp_path/'staged.txt.gz', 2)
    out = tmp_path/'formatted'; out.mkdir()
    env = {**os.environ, 'PYTHONPATH': packages}
    command = [sys.executable, '-c', 'from pgscatalog.core.cli.format_cli import main; main()',
               '-t', 'GRCh38', '--threads', '1', '-o', str(out), '-l', 'check.json', '-s']
    failed = subprocess.run([*command, str(source)], capture_output=True, text=True, env=env, timeout=30)
    assert failed.returncode and 'trait_reported' in failed.stderr
    passed = subprocess.run([*command, str(target)], capture_output=True, text=True, env=env, timeout=30)
    assert passed.returncode == 0, passed.stderr
    log = json.loads((out/'check.json').read_text())
    assert log
    assert 'PGS Catalog header not detected' not in passed.stderr


def test_summary_keeps_validation_field_and_ignores_retry_warnings():
    raw = ('ERROR ~ Error executing process > FORMAT_SCOREFILES\nTraceback (most recent call last):\n'
           '  pydantic_core._pydantic_core.ValidationError: 1 validation error for ScoreHeader\n'
           '  trait_reported\n    Field required [type=missing]\nWork dir:\n/irrelevant\n')
    assert prs_presentation.nextflow_failed(raw)
    assert not prs_presentation.nextflow_failed('WARN: Task failed -- Execution is retried (1)\n')
    assert 'trait_reported' in prs_presentation.log_summary(raw)
    assert 'Field required' in prs_presentation.log_summary(raw)
    assert '缺少必要' in prs_presentation.failure_reason(raw, 'exit=1')
