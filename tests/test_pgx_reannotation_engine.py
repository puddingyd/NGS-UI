"""Pinned real-tool integration uses generated calls, never patient fixtures."""
import copy
import json
import os
from pathlib import Path
import subprocess

import pytest

from app.services import pgx_reannotation as r


def test_version_mismatch_is_rejected_before_execution(tmp_path):
    with pytest.raises(ValueError, match="版本不相容"):
        r.reannotate({"pharmcatVersion": "other"}, {}, tmp_path,
                     lambda args: pytest.fail("must not execute"))


def test_changed_baseline_advice_is_never_accepted(tmp_path):
    raw = {"pharmcatVersion": r.VERSION, "dataVersion": r.DATA_VERSION,
           "genes": {"DPYD": {"sourceDiplotypes": [{"label": "c.1A>G"}, {"label": "c.2C>T"}]}}}
    def altered(args):
        result = copy.deepcopy(raw)
        result["genes"]["DPYD"]["sourceDiplotypes"].pop()
        r.atomic(tmp_path / "baseline.report.json", result)
    with pytest.raises(ValueError, match="無法完整重現"):
        r.reannotate(raw, {}, tmp_path, altered)


def test_runtime_rejects_substituted_jar(tmp_path, monkeypatch):
    fake = tmp_path / "pharmcat.jar"; fake.write_bytes(b"untrusted-version")
    monkeypatch.setenv("NGS_UI_PGX_PHARMCAT_JAR", str(fake))
    with pytest.raises(ValueError, match="校驗不符"):
        r.runtime()


@pytest.mark.skipif(not os.environ.get("NGS_UI_TEST_PHARMCAT_JAR"), reason="requires official pinned PharmCAT JAR")
def test_real_tool_recomputes_shared_guidance_preserves_calls_and_rejects_unknown(tmp_path, monkeypatch):
    monkeypatch.setenv("NGS_UI_PGX_PHARMCAT_JAR", os.environ["NGS_UI_TEST_PHARMCAT_JAR"])
    monkeypatch.delenv("NGS_UI_PGX_PHARMCAT_IMAGE", raising=False)
    engine = r.runtime()
    (tmp_path / "tmp").mkdir()
    def run(args):
        result = subprocess.run(r.command(engine, tmp_path, args), cwd=tmp_path,
                                capture_output=True, text=True, timeout=45)
        assert result.returncode == 0, result.stderr
    (tmp_path / "original.tsv").write_text("CYP2C9\t*1/*3\t\t\nHLA-B\t*40:01/*51:01\t\t\nSLCO1B1\t*1/*5\t\t\n")
    run(["-phenotyper", "-po", str(tmp_path / "original.tsv"), "-o", str(tmp_path), "-bf", "original"])
    run(["-reporter", "-ri", str(tmp_path / "original.phenotype.json"), "-reporterJson",
         "-rt", "SYNTHETIC", "-o", str(tmp_path), "-bf", "original"])
    raw = r.read(tmp_path / "original.report.json")
    new = r.reannotate(raw, {"CYP2C9": {"genotype": "*1/*1"}}, tmp_path, run)
    assert new["genes"]["CYP2C9"]["sourceDiplotypes"][0]["activityScore"] == "2.0"
    assert new["genes"]["CYP2C9"]["sourceDiplotypes"][0]["phenotypes"] == ["Normal Metabolizer"]
    for g in ("HLA-B", "SLCO1B1"):
        assert raw["genes"][g]["sourceDiplotypes"] == new["genes"][g]["sourceDiplotypes"]
    from app.adapters.pgx_tsv import _compact_drug_annotations
    old_annotations = _compact_drug_annotations(raw)
    new_annotations = _compact_drug_annotations(new)
    for genes in ({"CYP2C9", "HLA-B"}, {"CYP2C9", "SLCO1B1"}):
        before = [a for a in old_annotations if set(a["genes"]) == genes]
        after = [a for a in new_annotations if set(a["genes"]) == genes]
        assert before and after and before != after
    with pytest.raises(ValueError, match="無法辨識"):
        r.reannotate(raw, {"CYP2C9": {"genotype": "*1/*999"}}, tmp_path, run)
