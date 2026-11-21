"""docs/ops/pangolin_cache/pangolin_cache.py — persistent Pangolin score cache."""
import gzip
import importlib.util
from pathlib import Path

import pytest

TOOL = Path(__file__).resolve().parents[1] / "docs" / "ops" / "pangolin_cache" / "pangolin_cache.py"
spec = importlib.util.spec_from_file_location("pangolin_cache", TOOL)
pc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pc)

HEADER = [
    "##fileformat=VCFv4.2",
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS",
]
INFO_LINE = ('##INFO=<ID=Pangolin,Number=.,Type=String,Description="Pangolin splice scores. '
             'Format: gene|pos:score_change|pos:score_change|warnings,...">')


def _row(chrom, pos, alt, info="DP=10"):
    return f"{chrom}\t{pos}\t.\tA\t{alt}\t50\tPASS\t{info}\tGT\t0/1"


def _score(chrom, pos, alt):
    if pos % 3 == 0:
        return None  # Pangolin returned no score
    return f"ENSG{pos}|{pos % 7}:-0.{pos % 10}|1:0.0{alt == 'G' and 1 or 2}|Warnings:"


def _fake_pangolin(miss_vcf: Path, out_vcf: Path):
    """Pangolin-like output: PyVCF-style QUAL and INFO/Pangolin appended."""
    header, body = pc.read_vcf(str(miss_vcf))
    out = [line for line in header if not line.startswith("#CHROM")] + [INFO_LINE, header[-1]]
    for row in body:
        fields = row.split("\t")
        fields[5] = "50.0"
        value = _score(fields[0], int(fields[1]), fields[4])
        out.append(pc.with_pangolin("\t".join(fields), value))
    out_vcf.write_text("\n".join(out) + "\n")


def _ns_args(cache: Path):
    return ["--cache", str(cache), "--sif-name", "pangolin_cu121_1.0.0.sif",
            "--pangolin-db", "/nonexistent/gencode.v47.annotation.db",
            "--ref-fasta", "/nonexistent/hg38.fa", "--distance", "50", "--mask", "True"]


def _run_sample(tmp: Path, cache: Path, name: str, rows):
    work = tmp / name
    work.mkdir()
    vcf = work / "splice_for_pangolin.vcf"
    vcf.write_text("\n".join(HEADER + rows) + "\n")
    args = _ns_args(cache) + ["--sample", name]
    assert pc.main(["split", *args, "--in", str(vcf), "--miss-out", str(work / "miss.vcf"),
                    "--hits-out", str(work / "hits.tsv")]) == 0
    _, miss = pc.read_vcf(str(work / "miss.vcf"))
    if miss:
        _fake_pangolin(work / "miss.vcf", work / "scored.vcf")
    assert pc.main(["merge", *args, "--in", str(vcf), "--hits", str(work / "hits.tsv"),
                    "--pangolin-out", str(work / "scored.vcf"),
                    "--out", str(work / f"{name}.pangolin.vcf")]) == 0
    return work / f"{name}.pangolin.vcf", len(miss)


def _uncached_reference(tmp: Path, name: str, rows) -> Path:
    vcf = tmp / f"{name}.input.vcf"
    vcf.write_text("\n".join(HEADER + rows) + "\n")
    out = tmp / f"{name}.reference.vcf"
    _fake_pangolin(vcf, out)
    return out


def test_cached_output_matches_direct_pangolin_and_skips_known_variants(tmp_path):
    cache = tmp_path / "cache" / "pangolin.sqlite"
    s1 = [_row("chr1", p, "G") for p in range(100, 130)]
    s2 = [_row("chr1", p, "G") for p in range(115, 145)] + [_row("chr2", 7, "T", ".")]

    out1, scored1 = _run_sample(tmp_path, cache, "s1", s1)
    out2, scored2 = _run_sample(tmp_path, cache, "s2", s2)
    assert scored1 == 30
    assert scored2 == 16  # 15 shared variants came from the cache

    for name, rows, out in (("s1", s1, out1), ("s2", s2, out2)):
        ref = _uncached_reference(tmp_path, name, rows)
        assert pc.main(["compare", str(ref), str(out)]) == 0

    # Rows keep their order and the Pangolin header line is present once.
    header, body = pc.read_vcf(str(out2))
    assert [pc.row_key(r) for r in body] == [pc.row_key(r) for r in s2]
    assert sum(line.startswith("##INFO=<ID=Pangolin,") for line in header) == 1
    assert header[-1].startswith("#CHROM")
    # "no score" is cached as well: a repeat run sends nothing to Pangolin.
    _, scored_again = _run_sample(tmp_path, cache, "s1_again", s1)
    assert scored_again == 0


def test_merge_fails_when_pangolin_output_lacks_variants(tmp_path):
    cache = tmp_path / "c.sqlite"
    vcf = tmp_path / "in.vcf"
    vcf.write_text("\n".join(HEADER + [_row("chr1", 1, "G"), _row("chr1", 2, "G")]) + "\n")
    args = _ns_args(cache)
    pc.main(["split", *args, "--in", str(vcf), "--miss-out", str(tmp_path / "m.vcf"),
             "--hits-out", str(tmp_path / "h.tsv")])
    partial = tmp_path / "scored.vcf"
    partial.write_text("\n".join(HEADER + [_row("chr1", 1, "G", "Pangolin=X|1:0.1|2:0.0|")]) + "\n")
    assert pc.main(["merge", *args, "--in", str(vcf), "--hits", str(tmp_path / "h.tsv"),
                    "--pangolin-out", str(partial), "--out", str(tmp_path / "o.vcf")]) == 2
    assert not (tmp_path / "o.vcf").exists()
    # Nothing from the failed task was cached.
    assert pc.main(["split", *args, "--in", str(vcf), "--miss-out", str(tmp_path / "m2.vcf"),
                    "--hits-out", str(tmp_path / "h2.tsv")]) == 0
    assert len(pc.read_vcf(str(tmp_path / "m2.vcf"))[1]) == 2


def test_namespace_changes_with_pangolin_settings(tmp_path):
    fasta = tmp_path / "hg38.fa"
    fasta.write_text(">chr1\nACGT\n")

    def namespace(distance):
        args = pc.argparse.Namespace(sif_name="a.sif", pangolin_db=str(tmp_path / "db"),
                                     ref_fasta=str(fasta), distance=distance, mask="True")
        return pc.namespace_from_args(args)

    assert namespace("50") != namespace("100")
    assert "ref=hg38.fa:11" in namespace("50")
    assert "db=db:missing" in namespace("50")


def test_seed_skips_conflicting_keys(tmp_path):
    def gz(name, rows):
        path = tmp_path / name
        with gzip.open(path, "wt") as handle:
            handle.write("\n".join(HEADER + rows) + "\n")
        return str(path)

    a = gz("a.pangolin.vcf.gz", [_row("chr1", 1, "G", "DP=1;Pangolin=G1|1:0.1|2:0.0|"),
                                 _row("chr1", 2, "G", "DP=1;Pangolin=G1|1:0.5|2:0.0|"),
                                 _row("chr1", 3, "G")])
    b = gz("b.pangolin.vcf.gz", [_row("chr1", 1, "G", "DP=9;Pangolin=G1|1:0.1|2:0.0|"),
                                 _row("chr1", 2, "G", "DP=9;Pangolin=G1|1:0.6|2:0.0|"),
                                 _row("chr1", 3, "G")])
    cache = tmp_path / "seed.sqlite"
    conflicts = tmp_path / "conflicts.tsv"
    assert pc.main(["seed", *_ns_args(cache), "--conflicts", str(conflicts), a, b]) == 0
    lines = conflicts.read_text().splitlines()
    assert len(lines) == 2 and lines[1].startswith("chr1\t2\tA\tG\t")

    conn = pc.connect(str(cache))
    ns = pc.namespace_id(conn, pc.namespace_from_args(pc.argparse.Namespace(
        sif_name="pangolin_cu121_1.0.0.sif", pangolin_db="/nonexistent/gencode.v47.annotation.db",
        ref_fasta="/nonexistent/hg38.fa", distance="50", mask="True")))
    found = pc.lookup(conn, ns, {("chr1", 1, "A", "G"), ("chr1", 2, "A", "G"), ("chr1", 3, "A", "G")})
    assert found == {("chr1", 1, "A", "G"): "G1|1:0.1|2:0.0|", ("chr1", 3, "A", "G"): None}


@pytest.mark.parametrize("info,value,expected", [
    (".", "X", "Pangolin=X"),
    ("DP=3", None, "DP=3"),
    ("DP=3;Pangolin=old", "new", "DP=3;Pangolin=new"),
    (".", None, "."),
])
def test_with_pangolin(info, value, expected):
    assert pc.with_pangolin(_row("chr1", 1, "G", info), value).split("\t")[7] == expected


def test_apply_to_pipeline_previews_then_installs_with_backups(tmp_path):
    import subprocess
    import sys

    tool = TOOL.parent / "apply_to_pipeline.py"
    process = (TOOL.parent / "PANGOLIN_SCORE.nf").read_text(encoding="utf-8")
    process = process[process.index("process PANGOLIN_SCORE {"):]
    # A production-like module: the old (v3.6) process between two other blocks.
    old = (process.replace("${run_block}", "        pangolin \\\\\n            splice_for_pangolin.vcf \\\\\n")
                  .replace("pangolin_cache.py", "OLD"))
    code = tmp_path / "tertiary_code"
    (code / "modules").mkdir(parents=True)
    (code / "scripts").mkdir()
    module = code / "modules" / "snv_annotation.nf"
    module.write_text("process VEP_ANNOTATE {\n}\n\n" + old + "\nworkflow SNV_ANNOTATE {\n}\n",
                      encoding="utf-8")
    config = code / "nextflow_tertiary.config"
    config.write_text('params {\n    pangolin_db  = "x.db"\n}\n', encoding="utf-8")
    original = module.read_text(encoding="utf-8")
    stage = tmp_path / "stage"

    def run(*extra):
        return subprocess.run([sys.executable, str(tool), "--tertiary-code", str(code),
                               "--stage-dir", str(stage), *extra],
                              capture_output=True, text=True)

    preview = run()
    assert preview.returncode == 0, preview.stderr
    assert "PREVIEW ONLY" in preview.stdout
    assert module.read_text(encoding="utf-8") == original  # untouched

    assert run("--install").returncode == 1  # scripts/pangolin_cache.py missing
    (code / "scripts" / "pangolin_cache.py").write_bytes(TOOL.read_bytes())
    install = run("--install")
    assert install.returncode == 0, install.stderr
    new_module = module.read_text(encoding="utf-8")
    assert new_module.count("process PANGOLIN_SCORE {") == 1
    assert "pangolin_cache.py split" in new_module
    assert new_module.startswith("process VEP_ANNOTATE {") and "workflow SNV_ANNOTATE" in new_module
    assert 'pangolin_cache_db = "${params.ref_dir}/tertiary/pangolin/cache/pangolin_scores.sqlite"' \
        in config.read_text(encoding="utf-8")
    backups = sorted(p.name for p in code.rglob("*.bak-*"))
    assert len(backups) == 2
    assert (code / "modules" / backups[1]).read_text(encoding="utf-8") == original
    # Already installed: refuses instead of applying twice.
    assert run().returncode == 1
