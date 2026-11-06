"""Pinned real-tool integration uses generated calls, never patient fixtures."""
import copy
import importlib.util
import os
from pathlib import Path
import subprocess
import sys

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


def test_hospital_container_uses_actual_registered_hash_and_checks_version(tmp_path, monkeypatch):
    image = "/home/datalake_Intermediate/pipeline/nextflow_containers/pharmcat_3.2.0.sif"
    monkeypatch.setenv("NGS_UI_PGX_PHARMCAT_JAR", str(tmp_path / "missing.jar"))
    monkeypatch.delenv("NGS_UI_PGX_PHARMCAT_IMAGE", raising=False)
    monkeypatch.setenv("APPTAINER_BINDPATH", "/unrelated:/pharmcat")
    original_is_file = Path.is_file
    monkeypatch.setattr(Path, "is_file", lambda path: str(path) == image or original_is_file(path))
    calls = []
    def execute(args, **kwargs):
        calls.append(args)
        assert "APPTAINER_BINDPATH" not in kwargs["env"]
        if "sha256sum" in args:
            return subprocess.CompletedProcess(args, 0, r.HOSPITAL_JAR_SHA256 + "  /pharmcat/pharmcat.jar\n", "")
        return subprocess.CompletedProcess(args, 0, "PharmCAT 3.2.0\n", "")
    monkeypatch.setattr(subprocess, "run", execute)
    engine = r.runtime()
    assert engine == {"jar": "/pharmcat/pharmcat.jar", "image": image, "sha256": r.HOSPITAL_JAR_SHA256, "jar_on_host": False}
    assert image in calls[0] and image in calls[1]
    assert "-version" in calls[1]


def test_container_jar_is_not_shadowed_by_an_unrelated_host_file(tmp_path):
    jar = tmp_path / "pharmcat.jar"; jar.write_text("unrelated host copy")
    work = tmp_path / "work"; work.mkdir()
    engine = {"jar": str(jar), "image": "synthetic.sif", "jar_on_host": False}
    args = r.command(engine, work, ["-version"])
    assert args[args.index("--bind") + 1] == str(work)
    engine["jar_on_host"] = True
    args = r.command(engine, work, ["-version"])
    assert str(jar) + ":" + str(jar) + ":ro" in args[args.index("--bind") + 1]


@pytest.mark.parametrize("sha,version", [("0" * 64, "PharmCAT 3.2.0"), (r.HOSPITAL_JAR_SHA256, "PharmCAT 3.1.0")])
def test_matching_version_alone_or_matching_hash_alone_cannot_pass(tmp_path, monkeypatch, sha, version):
    monkeypatch.setenv("NGS_UI_PGX_PHARMCAT_JAR", str(tmp_path / "missing.jar"))
    monkeypatch.setenv("NGS_UI_PGX_PHARMCAT_IMAGE", "/synthetic/pharmcat.sif")
    def execute(args, **kwargs):
        return subprocess.CompletedProcess(args, 0, sha + "  /pharmcat/pharmcat.jar" if "sha256sum" in args else version, "")
    monkeypatch.setattr(subprocess, "run", execute)
    with pytest.raises(ValueError, match="校驗不符|版本不符"):
        r.runtime()


def test_baseline_with_different_data_version_cannot_publish(tmp_path):
    raw = {"pharmcatVersion": r.VERSION, "dataVersion": r.DATA_VERSION,
           "genes": {"CYP2C9": {"sourceDiplotypes": [{"label": "*1/*1"}]}}}
    def different_database(args):
        r.atomic(tmp_path / "baseline.report.json", {**raw, "dataVersion": "different-database"})
    with pytest.raises(ValueError, match="版本不相容"):
        r.reannotate(raw, {}, tmp_path, different_database)


@pytest.mark.parametrize("data_version", [r.DATA_VERSION, "old-database"])
def test_preflight_prints_actual_hospital_identity_and_rejects_other_data_version(monkeypatch, capsys, data_version):
    path = Path(__file__).resolve().parents[1] / "scripts/check_pgx_pharmcat.py"
    spec = importlib.util.spec_from_file_location("pgx_preflight", path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    engine = {"image": "/synthetic/hospital.sif", "jar": "/pharmcat/pharmcat.jar", "sha256": r.HOSPITAL_JAR_SHA256}
    monkeypatch.setattr(r, "runtime", lambda: engine)
    monkeypatch.setattr(sys, "argv", [str(path)])
    def execute(args, cwd, **kwargs):
        if "-reporter" in args:
            r.atomic(cwd / "probe.report.json", {"pharmcatVersion": r.VERSION, "dataVersion": data_version,
                "genes": {"CYP2C9": {"sourceDiplotypes": [{"label": "*1/*1", "phenotypes": ["Normal Metabolizer"], "activityScore": "2.0"}]}}})
        return subprocess.CompletedProcess(args, 0, "", "")
    monkeypatch.setattr(subprocess, "run", execute)
    if data_version != r.DATA_VERSION:
        with pytest.raises(ValueError, match="data version mismatch"):
            module.main()
        assert "PASS" not in capsys.readouterr().out
    else:
        module.main()
        output = capsys.readouterr().out
        assert "PASS" in output and r.HOSPITAL_JAR_SHA256 in output and engine["image"] in output


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
