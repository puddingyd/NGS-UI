import gzip
import importlib.util
import json
import sys
from pathlib import Path

import pytest

from app.services import sample_layout


def _script(name):
    path = Path(__file__).parents[1] / "scripts" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_read_only_scan_finds_json_and_vcf_calls_without_hard_filtered(tmp_path, monkeypatch):
    scanner = _script("scan_targeted_callers.py")
    post = tmp_path / "08_postprocessing"; post.mkdir()
    source = tmp_path / "raw" / "vcf.gz"; source.mkdir(parents=True)
    (post / "26T00054-dragen.pipeline_source.json").write_text(json.dumps({
        "source_sample_id": "SOURCE", "source_vcf_path": str(source / "SOURCE.hard-filtered.vcf.gz")
    }))
    (source / "SOURCE.targeted.json").write_text(json.dumps({
        "sampleId": "SOURCE", "genomeBuild": "hg38",
        "cyp21a2": {"totalCopyNumber": 4, "recombinantHaplotypes": ["NM_000500.9:c.710T>A", ""], "variants": []},
        "smn": {"smn1CopyNumber": 2, "smn2CopyNumber": 1,
                "variants": [{"hgvs": "NM_000344.4:c.*3+80T>G", "altCopyNumber": 1}]},
        "hba": {"genotype": "--/aa", "genotypeFilter": "PASS", "variants": []},
    }))
    with gzip.open(source / "SOURCE.targeted.vcf.gz", "wt") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSOURCE\n")
        f.write("chr6\t32038514\t.\tC\tT\t50\tPASS\tALLELE_ID=.,NM_000500.9:c.92C>T\tGT\t0/1\n")
        f.write("chr16\t173000\t.\tC\tT\t50\tPASS\tALLELE_ID=.,NM_000517.6:c.427T>C\tGT\t0/1\n")
        f.write("chr5\t70952074\t.\tT\tG\t0\tTargetedRepeatConflict\tEVENT=NM_000344.4:c.*3+80T>G\tGT\t0/1\n")
        f.write("chr6\t32038517\t.\tC\tT\t0\tTargetedWT\tALLELE_ID=.,NM_000500.9:c.713T>A\tGT\t0/0\n")
    monkeypatch.setattr(sample_layout, "state_file", lambda sample, name: post / f"{sample}.{name}")
    rows = scanner.scan_case("26T00054-dragen")
    assert {r["caller"] for r in rows} == {"cyp21a2", "smn", "hba"}
    assert all(r["status"] == "candidate" for r in rows)
    assert "SMN1_CN" not in next(r for r in rows if r["caller"] == "smn")["json_signals"]
    assert "not diagnostic alone" in next(r for r in rows if r["caller"] == "smn")["json_signals"]
    assert "713T>A" not in next(r for r in rows if r["caller"] == "cyp21a2")["vcf_signals"]
    assert "SMN2_CN" in scanner.json_signals("smn", {"smn1CopyNumber": 2, "smn2CopyNumber": 3})[0]
    (source / "SOURCE.targeted.json").unlink()
    vcf_only = scanner.scan_case("26T00054-dragen")
    assert {r["caller"] for r in vcf_only} == {"cyp21a2", "smn", "hba"}
    assert all(r["json_signals"] == "" for r in vcf_only)
    assert not (source / "SOURCE.hard-filtered.vcf.gz").exists()


def test_batch_backfill_uses_dragen_suffix_skips_complete_and_continues_errors(tmp_path, monkeypatch, capsys):
    rebuilder = _script("rebuild_specialized.py")
    existing = tmp_path / "a-dragen.specialized.manifest.json"; existing.write_text("{}")
    monkeypatch.setattr(sample_layout, "iter_sample_ids", lambda: ["c-nckuh", "b-dragen", "a-dragen"])
    monkeypatch.setattr(sample_layout, "state_file", lambda sample, name: existing if sample == "a-dragen" else tmp_path / f"{sample}.{name}")
    monkeypatch.setattr(rebuilder.specialized, "load", lambda sample: {"status": "complete"})
    assert list(rebuilder.batch_samples()) == [("a-dragen", "skipped_existing"), ("b-dragen", "pending")]
    calls = []
    def fake_rebuild(sample, *, from_source, apply, source_vcf=None):
        calls.append((sample, from_source, apply))
        raise ValueError("source missing")
    monkeypatch.setattr(rebuilder, "rebuild_one", fake_rebuild)
    monkeypatch.setattr(sys, "argv", ["rebuild_specialized.py", "--all-dragen", "--from-source"])
    with pytest.raises(SystemExit) as exc:
        rebuilder.main()
    assert exc.value.code == 1
    assert calls == [("b-dragen", True, False)]
    assert '"skipped_existing": 1' in capsys.readouterr().out
