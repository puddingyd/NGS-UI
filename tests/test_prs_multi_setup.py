"""Release integrity and activation guards, using synthetic resources only."""
import json

import pytest

from app.services import prs
from research.prs.absolute_risk.registry import DISEASES
from research.prs.multidisease.common import MODEL_FILES, load_models, sha, write_json
from scripts import setup_prs_multi as setup


def test_release_rejects_missing_artifacts_wrong_identity_and_tampering(tmp_path):
    for name in MODEL_FILES:
        (tmp_path/name).write_text('synthetic artifact')
    models = {k:dict(pgs_id=d['pgs_id'], n_loci=610305 if k=='CHD' else 1000,
                     source_loci=610305 if k=='CHD' else 1000) for k,d in DISEASES.items()}
    manifest = dict(schema='PRS_MULTI_MODELS_V1', reference_id='unrelated3330_EAS735_v1',
        reference_n=3330, eas_n=735, models=models, offsets={k:'0' for k in DISEASES},
        chd_reference_regression='PASS', sha256={name:sha(tmp_path/name) for name in MODEL_FILES})
    write_json(tmp_path/'manifest.json',manifest)
    assert load_models(tmp_path)[0] == manifest
    changed=json.loads(json.dumps(manifest));changed['sha256'].pop('reference.afreq')
    write_json(tmp_path/'manifest.json',changed)
    with pytest.raises(ValueError):load_models(tmp_path)
    changed=json.loads(json.dumps(manifest));changed['models']['T2D']['pgs_id']='PGS004875'
    write_json(tmp_path/'manifest.json',changed)
    with pytest.raises(ValueError):load_models(tmp_path)
    write_json(tmp_path/'manifest.json',manifest)
    (tmp_path/'weights.tsv').write_text('changed weight')
    with pytest.raises(ValueError,match='checksum'):load_models(tmp_path)


def test_activation_refuses_active_jobs_and_concurrent_config_edit(tmp_path, monkeypatch):
    monkeypatch.setattr(setup.config,'DATA_ROOT',tmp_path/'data')
    path=tmp_path/'cfg.json';prs.atomic_json(path,dict(project='existing'))
    original=sha(path); output=tmp_path/'models'
    state=tmp_path/'data/jobs/tertiary/job/state.json'
    prs.atomic_json(state,dict(state='running'))
    with pytest.raises(ValueError,match='執行中'):setup.activate(path,original,output)
    assert sha(path)==original
    prs.atomic_json(state,dict(state='done'))
    prs.atomic_json(path,dict(project='changed'))
    with pytest.raises(ValueError,match='changed'):setup.activate(path,original,output)
    current=path.read_bytes()
    setup.activate(path,sha(path),output)
    assert json.loads(path.read_text())==dict(project='changed',multi_model_dir=str(output))
    assert next(tmp_path.glob('cfg.json.before-multi-*')).read_bytes()==current
