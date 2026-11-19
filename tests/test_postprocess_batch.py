"""Batch post-processing must produce exactly the per-sample results.

Every step of scripts/postprocess_batch.py runs once for all samples (one DB
pass per resource). This builds small real resources for each step, runs the
same three overlapping samples as one batch and one at a time (through the
run_stopgaps.sh wrapper), and requires every output file to match.
"""
import csv
import gzip
import importlib.util
import json
import os
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

from app.services import clinvar_latest_store, litvar2_store

REPO = Path(__file__).resolve().parents[1]
FIELDS = ["CHROM", "POS", "REF", "ALT", "RS_ID", "GENE", "TRANSCRIPT", "HGVS_C", "HGVS_P",
          "CONSEQUENCE", "IMPACT", "GNOMAD_G_AF", "DP", "DP_DV", "VAF_DV", "GT_DV",
          "CLINVAR_SIG", "CLINVAR_STARS", "CLINVAR_DN", "CLINVAR_SIGCONF", "CLINVAR_VARIATION_ID",
          "ACMG_CLASS"]
VARIANTS = {
    "a": dict(CHROM="chr1", POS="100", REF="A", ALT="G", RS_ID="rs123", GENE="GENE1",
              TRANSCRIPT="ENST00000000001.5", HGVS_C="c.1A>G", GNOMAD_G_AF="0.001",
              CLINVAR_SIG="Uncertain significance", CLINVAR_VARIATION_ID="1001"),
    "b": dict(CHROM="chr1", POS="200", REF="C", ALT="T", GENE="GENE2", GNOMAD_G_AF="",
              CLINVAR_SIG="Pathogenic", CLINVAR_VARIATION_ID="1002"),
    "c": dict(CHROM="chr2", POS="10", REF="A", ALT="G", GENE="GENE4", GNOMAD_G_AF="0.0001"),
    "d": dict(CHROM="chr2", POS="11", REF="A", ALT="G", GENE="GENE4", GNOMAD_G_AF="0.5"),
    "e": dict(CHROM="chr3", POS="500", REF="AT", ALT="A", GENE="GENE6", GNOMAD_G_AF="0.004", DP="8", DP_DV="8"),
}
SAMPLES = {"S1": (["a", "b", "c"], "WES"), "S2": (["a", "c", "d", "e"], "WES"), "S3": (["b", "d", "e"], "WGS")}


def _resources(root: Path) -> list[str]:
    with gzip.open(root / "genebe.tsv.gz", "wt") as f:
        f.write("#chr\tpos\tref\talt\tacmg_classification\tacmg_score\tacmg_criteria\n"
                "chr1\t100\tA\tG\tLikely_benign\t-4\tBP4\nchr2\t10\tA\tG\tVUS\t2\tPM2\n")
    with gzip.open(root / "clinvar.vcf.gz", "wt") as f:
        f.write("##fileformat=VCFv4.2\n##fileDate=2026-08-03\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "1\t100\t1001\tA\tG\t.\tPASS\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;CLNDN=X\n"
                "2\t11\t2011\tA\tG\t.\tPASS\tCLNSIG=Likely_pathogenic;CLNREVSTAT=reviewed_by_expert_panel;CLNDN=Y\n")
    clinvar_latest_store.build_database(root / "clinvar.vcf.gz", root / "clinvar.sqlite")
    giab = root / "giab"
    giab.mkdir()
    (giab / "h.bed").write_text("chr1\t90\t150\n")
    (giab / "strata_manifest.json").write_text(json.dumps([{"file": "h.bed", "label": "homopolymer"}]))
    with gzip.open(root / "inhouse.vcf.gz", "wt") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t100\t.\tA\tG\t.\t.\tINHOUSE_AC=3;INHOUSE_AN=100;INHOUSE_AF=0.03\n"
                "chr3\t500\t.\tAT\tA\t.\t.\tINHOUSE_AC=2;INHOUSE_AN=100;INHOUSE_AF=0.02\n")
    with gzip.open(root / "mane.txt.gz", "wt") as f:
        f.write("Ensembl_nuc\tRefSeq_nuc\tRefSeq_prot\tMANE_status\n"
                "ENST00000000001.5\tNM_000001.1\tNP_000001.1\tMANE Select\n")
    spec = importlib.util.spec_from_file_location("litvar_fixture", REPO / "tests" / "test_litvar2.py")
    fixture = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixture)
    litvar2_store.build_database(fixture._write_bulk(root / "lit.json.gz"), root / "litvar2.sqlite",
                                 dataset_date="2026-07-01", source_url="https://example.test/l.json.gz")
    msa = root / "msa.tsv.bgz"
    msa.write_text("1\t100\tA\tG\t-3.5\n2\t10\tA\tG\t-7\n")
    Path(f"{msa}.tbi").touch()
    (root / "cds.bed").write_text("chr1\t50\t700\nchr2\t1\t50\nchr3\t400\t600\n")
    return ["--genebe-db", str(root / "genebe.tsv.gz"), "--clinvar-latest-db", str(root / "clinvar.sqlite"),
            "--giab-strat-dir", str(giab), "--inhouse-af-db", str(root / "inhouse.vcf.gz"),
            "--mane-summary", str(root / "mane.txt.gz"), "--litvar2-db", str(root / "litvar2.sqlite"),
            "--gpn-msa-db", str(msa), "--gpn-star-dir", str(root / "no-star"),
            "--candidate-bed", str(root / "cds.bed"), "--skip-spliceai"]


def _inputs(cases: Path) -> list[dict]:
    shutil.rmtree(cases, ignore_errors=True)
    entries = []
    for sid, (keys, seq) in SAMPLES.items():
        raw = cases / sid / "03_acmg" / f"{sid}.snv_indel.acmg.tsv"
        raw.parent.mkdir(parents=True)
        with raw.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for key in keys:
                writer.writerow({f: "" for f in FIELDS} | dict(
                    CONSEQUENCE="missense_variant", IMPACT="MODERATE", DP="40", DP_DV="40",
                    VAF_DV="0.5", GT_DV="0/1", ACMG_CLASS="VUS") | VARIANTS[key])
        os.utime(raw, (1_700_000_000, 1_700_000_000))
        post = cases / sid / "08_postprocessing"
        post.mkdir()
        work = post / ".snv_indel.job.working.tsv"
        shutil.copyfile(raw, work)
        entries.append(dict(sample=sid, raw_tsv=str(raw), work_tsv=str(work), post_dir=str(post), seq_type=seq))
    return entries


def _scrub(value):
    """Drop run-time values (timestamps, rebuilt-file mtimes)."""
    if isinstance(value, dict):
        return {k: ("<volatile>" if k.endswith(("_at", "mtime_ns", "ctime_ns")) else _scrub(v))
                for k, v in value.items()}
    if isinstance(value, list):
        return [_scrub(v) for v in value]
    return value


def _files(root: Path) -> dict:
    out = {}
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        rel = str(path.relative_to(root))
        if path.suffix == ".sqlite":
            with sqlite3.connect(path) as conn:
                out[rel] = [(name, sorted(map(repr, conn.execute(f"SELECT * FROM '{name}'"))))
                            for (name,) in conn.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")]
        elif path.suffix == ".json":
            out[rel] = _scrub(json.loads(path.read_text()))
        else:
            out[rel] = path.read_bytes()
    return out


def test_batch_matches_one_sample_at_a_time(tmp_path):
    flags = _resources(tmp_path)
    tabix = tmp_path / "tabix"
    tabix.write_text(f"#!{sys.executable}\nimport sys\nrows=[l.rstrip('\\n').split('\\t') for l in open(sys.argv[1])]\n"
                     "for reg in sys.argv[2:]:\n c,r=reg.split(':');p=int(r.split('-')[0])\n"
                     " [print('\\t'.join(x)) for x in rows if x[0]==c and int(x[1])==p]\n")
    tabix.chmod(0o755)
    env = {**os.environ, "NGS_UI_PYTHON_BIN": sys.executable, "TABIX_BIN": str(tabix),
           "NGS_UI_DATA_ROOT": str(tmp_path / "data"), "GENEBE_USER": "", "GENEBE_API_KEY": ""}
    env.pop("NGS_UI_CDS_CANDIDATE_BED", None)
    cases = tmp_path / "cases"

    for entry in _inputs(cases):
        proc = subprocess.run(["bash", str(REPO / "scripts/run_stopgaps.sh"), "--work-tsv", entry["work_tsv"],
                               "--raw-tsv", entry["raw_tsv"], "--post-dir", entry["post_dir"],
                               "--sample", entry["sample"], "--seq-type", entry["seq_type"], *flags],
                              env=env, capture_output=True, text=True)
        assert proc.returncode == 0, proc.stdout + proc.stderr
    single = _files(cases)

    manifest = tmp_path / "batch.json"
    manifest.write_text(json.dumps({"samples": _inputs(cases)}))
    proc = subprocess.run([sys.executable, str(REPO / "scripts/postprocess_batch.py"), "--manifest", str(manifest),
                           "--raw-log", str(tmp_path / "raw.log"), *flags], env=env, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    batch = _files(cases)

    assert batch.keys() == single.keys()
    assert [k for k in batch if batch[k] != single[k]] == []
    header = (cases / "S2/08_postprocessing/.snv_indel.job.working.tsv").read_text().splitlines()[0]
    for column in ("GENEBE_ACMG_SCORE", "CLINVAR_LATEST_SIG", "GIAB_STRATA", "INHOUSE_AF", "REFSEQ_NUC",
                   "LITVAR2_STATUS"):
        assert column in header
    review = (cases / "S1/08_postprocessing/S1.snv_indel.review.tsv").read_text()
    assert "-3.5" in review  # GPN-MSA served through the batch prefetch
    # each batch-wide step ran once; each sample step once per sample
    log = proc.stdout
    assert log.count("[post-processing-step] genebe start") == 1
    assert log.count("[sample-step] review-tsv done sample=") == 3
    assert "[gpn-msa] status=complete" in log
