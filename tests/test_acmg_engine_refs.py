"""ACMG rule engine reference tables and WES recurrence store — synthetic files only."""
import gzip
import json
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from app.services.acmg_engine import Context, evaluate_variant
from app.services.acmg_engine import parse as P
from app.services.acmg_engine import recurrence as R
from app.services.acmg_engine import ref_store as S

REPO = Path(__file__).resolve().parents[1]

CLINGEN = """#ClinGen dosage sensitivity curation
#Gene Symbol\tGene ID\tcytoBand\tGenomic Location\tHaploinsufficiency Score\tHaploinsufficiency Description
GENE_HI3\t1\t1p36\tchr1:1-2\t3\tSufficient evidence
GENE_HI40\t2\t1p36\tchr1:3-4\t40\tDosage sensitivity unlikely
GENE_HI1\t3\t1p36\tchr1:5-6\t1\tLittle evidence
"""

MANE_COLUMNS = ["#NCBI_GeneID", "Ensembl_Gene", "HGNC_ID", "symbol", "name", "RefSeq_nuc", "RefSeq_prot",
                "Ensembl_nuc", "Ensembl_prot", "MANE_status", "GRCh38_chr", "chr_start", "chr_end", "chr_strand"]
MANE_ROWS = [
    ["1", "ENSG1", "HGNC:1", "GENE1", "gene one", "NM_000001.3", "NP_000001.3", "ENST00000000001.1",
     "ENSP00000000001.1", "MANE Select", "1", "1", "5000", "+"],
    ["4", "ENSG4", "HGNC:4", "GENE_HI40", "gene", "NM_000040.1", "NP_000040.1", "ENST00000000040.1",
     "ENSP00000000040.1", "MANE Select", "1", "1", "5000", "+"],
]

VS_COLUMNS = ["#AlleleID", "Type", "Name", "GeneSymbol", "ClinicalSignificance", "Assembly", "Chromosome",
              "ReviewStatus", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"]
TWO = "criteria provided, multiple submitters, no conflicts"
ONE = "criteria provided, single submitter"
NONE = "no assertion criteria provided"


def vs(name, sig, review=ONE, assembly="GRCh38", chrom="1", pos="1000", ref="C", alt="T"):
    return ["1", "single nucleotide variant", name, "", sig, assembly, chrom, review, pos, ref, alt]


VS_ROWS = [
    vs("NM_000001.3(GENE1):c.100C>T (p.Arg34Ter)", "Pathogenic", TWO),
    vs("NM_000001.3(GENE1):c.400del (p.Leu134fs)", "Likely pathogenic", pos="1400", ref="CT", alt="C"),
    vs("NM_000001.3(GENE1):c.456+1G>C", "Pathogenic", pos="1500", ref="G", alt="C"),
    vs("NM_000001.3(GENE1):c.101G>A (p.Arg34Gln)", "Pathogenic", "reviewed by expert panel", pos="1001",
       ref="G", alt="A"),
    vs("NM_000001.3(GENE1):c.102G>A (p.Arg34=)", "Benign", TWO, pos="1002", ref="G", alt="A"),
    vs("NM_000001.3(GENE1):c.600C>T (p.Arg200Ter)", "Pathogenic", TWO, assembly="GRCh37", pos="9"),
    vs("NM_000001.3(GENE1):c.700C>T (p.Arg234Ter)", "Pathogenic", NONE, pos="1700"),
    vs("NM_000001.3(GENE1):c.800C>T (p.Arg267Ter)", "Conflicting classifications of pathogenicity", ONE,
       pos="1800"),
    vs("NM_999999.1(GENE2):c.10C>T (p.Arg4Ter)", "Pathogenic", ONE, chrom="2", pos="10"),
    vs("NM_999999.1(GENE2):c.20C>T (p.Arg7Ter)", "Pathogenic", ONE, chrom="2", pos="20"),
    vs("NM_888888.1(GENE3):c.30C>T (p.Arg10Ter)", "Pathogenic", TWO, chrom="3", pos="30"),
    vs("NM_000040.1(GENE_HI40):c.10C>T (p.Arg4Ter)", "Pathogenic", TWO, chrom="4", pos="10"),
    vs("NM_000040.1(GENE_HI40):c.20C>T (p.Arg7Ter)", "Pathogenic", TWO, chrom="4", pos="20"),
]


def _gtf(tx, number, kind, start, end):
    attrs = f'gene_id "ENSG1"; transcript_id "{tx}.1"; exon_number "{number}";'
    return "\t".join(["chr1", "MANE", kind, str(start), str(end), ".", "+", ".", attrs])


def write_gz(path: Path, text: str) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)
    return path


def tsv(rows):
    return "".join("\t".join(r) + "\n" for r in rows)


@pytest.fixture()
def sources(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    gtf_lines = [
        "#!genome-build GRCh38",
        _gtf("ENST00000000001", 1, "exon", 1, 200), _gtf("ENST00000000001", 1, "CDS", 141, 200),
        _gtf("ENST00000000001", 2, "exon", 301, 400), _gtf("ENST00000000001", 2, "CDS", 301, 400),
        _gtf("ENST00000000001", 3, "exon", 501, 800), _gtf("ENST00000000001", 3, "CDS", 501, 640),
    ]
    return {
        "clingen_hi": raw / "clingen.tsv",
        "mane_summary": write_gz(raw / "mane.summary.txt.gz", tsv([MANE_COLUMNS] + MANE_ROWS)),
        "clinvar_summary": write_gz(raw / "variant_summary.txt.gz", tsv([VS_COLUMNS] + VS_ROWS)),
        "mane_gtf": write_gz(raw / "MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz", "\n".join(gtf_lines) + "\n"),
        "_write": (raw / "clingen.tsv").write_text(CLINGEN, encoding="utf-8"),
    }


def build(tmp_path, sources, **kwargs):
    files = {k: v for k, v in sources.items() if not k.startswith("_")}
    files.update(kwargs)
    counts = S.build_reference_db(tmp_path / "refs", clinvar_release="2026-11-15", **files)
    return counts, tmp_path / "refs" / S.REFS_DB_NAME


# ── parsers / builder ──────────────────────────────────────────────────────
def test_parsers(sources):
    assert S.parse_clingen_hi(sources["clingen_hi"]) == {"GENE_HI3": "3", "GENE_HI40": "40", "GENE_HI1": "1"}
    mane = S.parse_mane_summary(sources["mane_summary"])
    assert mane["NM_000001"] == {"gene": "GENE1", "ensembl": "ENST00000000001", "status": "MANE Select"}
    exons = S.parse_mane_gtf(sources["mane_gtf"])
    assert exons[("ENST00000000001", 1)] == {"cds_length": 60, "frameshift": False, "fraction": 0.2,
                                              "last_exon": False}
    assert exons[("ENST00000000001", 2)]["frameshift"] is True
    assert exons[("ENST00000000001", 3)]["last_exon"] is True
    assert S.review_stars("Reviewed by expert panel") == 3
    assert S.review_stars(NONE) == 0


def test_clinvar_records_and_classification(sources):
    records = list(S.iter_clinvar_plp(sources["clinvar_summary"]))
    names = {r["hgvs_c"] for r in records}
    assert "c.102G>A" not in names          # benign
    assert "c.600C>T" not in names          # GRCh37
    assert "c.800C>T" not in names          # conflicting
    first = next(r for r in records if r["hgvs_c"] == "c.100C>T")
    assert first["variant_key"] == "chr1-1000-C-T" and first["stars"] == 2 and first["gene"] == "GENE1"
    kinds = {r["hgvs_c"]: S.classify_clinvar(r)[0] for r in records}
    assert kinds["c.100C>T"] == "truncating"
    assert kinds["c.400del"] == "truncating"
    assert kinds["c.456+1G>C"] == "splice"
    assert kinds["c.101G>A"] == "missense"


def test_build_reference_db(tmp_path, sources):
    counts, db = build(tmp_path, sources)
    assert counts["truncating"] == 4          # GENE1 ×2 + GENE_HI40 ×2 (MANE, ≥1★); 0★ and non-MANE excluded
    assert counts["missense"] == 1
    assert counts["exons"] == 3
    conn = sqlite3.connect(db)
    lof = {g: (d, n) for g, d, n in conn.execute("SELECT gene, decision, truncating_count FROM lof_genes")}
    assert lof["GENE1"] == ("yes", 3)         # two truncating + one canonical splice, ≥1★
    assert lof["GENE2"] == ("yes", 2)         # any transcript counts toward the mechanism
    assert "GENE3" not in lof                 # one variant is not enough
    assert lof["GENE_HI3"][0] == "yes"
    assert lof["GENE_HI40"][0] == "no"        # HI=40 overrides ClinVar truncating variants
    assert "GENE_HI1" not in lof
    meta = dict(conn.execute("SELECT key, value FROM meta"))
    assert meta["clinvar_release"] == "2026-11-15"
    assert json.loads(meta["sources"])["mane_gtf"]["size"] > 0
    conn.close()
    assert not list((tmp_path / "refs").glob(".*.tmp"))


def test_sqlite_references(tmp_path, sources):
    _, db = build(tmp_path, sources)
    refs = S.SqliteReferences(db, ba1_exceptions={"chr9-9-A-G"})
    assert refs.available == {"lof_genes", "truncating_positions", "missense_by_residue", "exons"}
    assert refs.lof_mechanism("GENE1")[0] == "yes"
    assert refs.lof_mechanism("GENE_HI40")[0] == "no"
    assert refs.lof_mechanism("GENE9") is None
    assert refs.truncating_downstream("GENE1", 100) is True
    assert refs.truncating_downstream("GENE1", 200) is False
    assert refs.clinvar_missense("GENE1", 34) == [{"alt_aa": "Gln", "variant_key": "chr1-1001-G-A", "stars": 3}]
    assert refs.exon_skip("ENST00000000001.1", 2) == {"frameshift": True, "fraction": pytest.approx(1 / 3),
                                                       "last_exon": False}
    assert refs.exon_skip("ENST00000000001", 9) is None
    assert refs.ba1_exception("chr9-9-A-G") and not refs.ba1_exception("chr1-1-A-G")
    assert refs.wes_recurrence("chr1-1000-C-T") is None
    assert refs.version["clinvar_release"] == "2026-11-15"


def test_without_gtf_exons_are_not_loaded(tmp_path, sources):
    _, db = build(tmp_path, sources, mane_gtf=None)
    refs = S.SqliteReferences(db)
    assert "exons" not in refs.available
    assert refs.exon_skip("ENST00000000001", 2) is None


def test_missing_db_behaves_as_not_loaded(tmp_path):
    refs = S.SqliteReferences(tmp_path / "missing.sqlite")
    assert refs.available == frozenset()
    assert refs.lof_mechanism("GENE1") is None


def test_ba1_exception_file_is_empty_by_default():
    assert S.load_ba1_exceptions() == set()


# ── engine on real (synthetic) reference database ─────────────────────────
def variant(**overrides):
    base = {
        "CHROM": "chr1", "POS": "1000", "REF": "C", "ALT": "T", "GENE": "GENE1",
        "TRANSCRIPT": "ENST00000000001.1", "TRANSCRIPT_TYPE": "MANE_SELECT",
        "HGVS_C": "ENST00000000001.1:c.100C>T", "HGVS_P": "ENSP00000000001.1:p.Arg34Ter",
        "CONSEQUENCE": "stop_gained", "IMPACT": "HIGH", "EXON": "2/3", "INTRON": ".",
        "PROTEIN_POSITION": "34/300", "NMD": "", "LOFTEE": "HC", "LOFTEE_FILTER": ".",
        "GNOMAD_E_AF": ".", "GNOMAD_G_AF": ".", "PANGOLIN_SCORE": "0.01", "PHYLOP100": "1.0",
        "CLINVAR_SIG": ".", "CLINVAR_STARS": ".", "GIAB_STRATA": "",
    }
    base.update(overrides)
    return base


def test_engine_uses_sqlite_references(tmp_path, sources):
    _, db = build(tmp_path, sources)
    ctx = Context(refs=S.SqliteReferences(db))
    nmd = evaluate_variant(variant(), ctx)
    assert [e["strength"] for e in nmd["evidence"] if e["code"] == "PVS1" and e["status"] == "met"] == ["very_strong"]
    escape = evaluate_variant(variant(NMD="NMD_escaping_variant", PROTEIN_POSITION="20/300"), ctx)
    pvs1 = [e for e in escape["evidence"] if e["code"] == "PVS1"][0]
    assert pvs1["strength"] == "strong" and "下游" in pvs1["reason"]
    splice = evaluate_variant(variant(CONSEQUENCE="splice_donor_variant", HGVS_C="ENST00000000001.1:c.456+1G>C",
                                      HGVS_P=".", INTRON="2/2", POS="1500", REF="G", ALT="C"), ctx)
    assert [e for e in splice["evidence"] if e["code"] == "PVS1"][0]["strength"] == "very_strong"
    missense = evaluate_variant(variant(CONSEQUENCE="missense_variant", IMPACT="MODERATE", LOFTEE=".",
                                        HGVS_P="ENSP00000000001.1:p.Arg34Trp"), ctx)
    assert [e for e in missense["evidence"] if e["code"] == "PM5" and e["status"] == "met"]
    not_lof = evaluate_variant(variant(GENE="GENE_HI40"), ctx)
    assert [e for e in not_lof["evidence"] if e["code"] == "PVS1"][0]["status"] == "not_met"


# ── WES recurrence ─────────────────────────────────────────────────────────
REVIEW_HEADER = ["CHROM", "POS", "REF", "ALT", "GENE"]


def make_case(root: Path, sid: str, variants, *, test_type="WES", mrn="", source=None):
    post = root / sid / "08_postprocessing"
    post.mkdir(parents=True)
    rows = [REVIEW_HEADER] + [[c, str(p), r, a, "G"] for c, p, r, a in variants]
    (post / f"{sid}.snv_indel.review.tsv").write_text(tsv(rows), encoding="utf-8")
    (post / f"{sid}.snv_indel.review.tsv.source.json").write_text(json.dumps({"test_type": test_type}))
    (post / f"{sid}.pipeline_source.json").write_text(json.dumps(source or {"seq_type": test_type}))
    if mrn:
        (post / f"{sid}.sample_metadata.json").write_text(json.dumps({"mrn": mrn}))
    return post


def test_patient_keys():
    assert R.base_sample_id("24WE0001_2-nckuh") == "24WE0001_2"
    assert R.patient_key("24WE0001_2-nckuh") == "id:24WE0001"
    assert R.patient_key("24WE0001-dragen") == "id:24WE0001"
    assert R.patient_key("24WE0001_3", "123") == "mrn:123"
    assert R.is_test_sample("23WE0044T-nckuh")
    assert not R.is_test_sample("26T00970-dragen")      # TITAN IDs contain T but do not end with it
    assert not R.is_test_sample("24WE0001_2-nckuh")


SHARED = ("chr1", 1000, "C", "T")


def test_sync_merges_patients_and_skips_test_and_wgs(tmp_path):
    root = tmp_path / "tertiary_output"
    make_case(root, "24WE0001-nckuh", [SHARED, ("chr2", 5, "A", "G"), ("chrM", 73, "A", "G")])
    make_case(root, "24WE0001_2-nckuh", [SHARED])                       # same patient by ID
    make_case(root, "24WE0002-nckuh", [SHARED], mrn="M1")
    make_case(root, "24WE0099-nckuh", [SHARED], mrn="M1")               # same patient by MRN
    make_case(root, "24WE0003-nckuh", [SHARED], mrn="M3")
    make_case(root, "23WE0044T-nckuh", [SHARED])                        # test sample
    make_case(root, "26T00970-dragen", [SHARED], test_type="WGS")       # TITAN WGS
    make_case(root, "24WG0001-dragen", [SHARED], test_type="WGS")
    (root / "24WE0500-nckuh" / "08_postprocessing").mkdir(parents=True)  # no review TSV
    (root / ".staging").mkdir()
    store = R.WesRecurrence(tmp_path / "rec.sqlite")
    summary = R.sync_cases(store, root, roster={})
    assert sorted(summary["added"]) == ["24WE0001-nckuh", "24WE0001_2-nckuh", "24WE0002-nckuh",
                                        "24WE0003-nckuh", "24WE0099-nckuh"]
    assert set(summary["skipped"]) == {"23WE0044T-nckuh", "26T00970-dragen", "24WG0001-dragen", "24WE0500-nckuh"}
    assert summary["patients"] == 3
    assert store.lookup("chr1-1000-C-T") == (3, 3)
    assert store.lookup("chr2-5-A-G") == (1, 3)
    assert store.lookup("chrM-73-A-G") == (0, 3)                        # mitochondrial rows skipped
    again = R.sync_cases(store, root, roster={})
    assert again["added"] == [] and len(again["unchanged"]) == 5
    store.close()


def test_rerun_replaces_and_removed_case_is_pruned(tmp_path):
    root = tmp_path / "tertiary_output"
    post = make_case(root, "24WE0001-nckuh", [SHARED])
    make_case(root, "24WE0002-nckuh", [SHARED])
    store = R.WesRecurrence(tmp_path / "rec.sqlite")
    R.sync_cases(store, root, roster={})
    review = post / "24WE0001-nckuh.snv_indel.review.tsv"
    review.write_text(tsv([REVIEW_HEADER, ["chr3", "7", "G", "T", "G"]]), encoding="utf-8")
    os.utime(review, ns=(1, 2_000_000_000_000_000_000))
    summary = R.sync_cases(store, root, ["24WE0001-nckuh"], roster={})
    assert summary["updated"] == ["24WE0001-nckuh"]
    assert store.lookup("chr1-1000-C-T") == (1, 2)
    assert store.lookup("chr3-7-G-T") == (1, 2)
    # case turns into a test-type mismatch → its old contribution is dropped
    (post / "24WE0001-nckuh.snv_indel.review.tsv.source.json").write_text(json.dumps({"test_type": "WGS"}))
    (post / "24WE0001-nckuh.pipeline_source.json").write_text(json.dumps({"seq_type": "WGS"}))
    summary = R.sync_cases(store, root, ["24WE0001-nckuh"], roster={})
    assert summary["removed"] == ["24WE0001-nckuh"]
    # deleted directory → pruned only by a full scan
    import shutil
    shutil.rmtree(root / "24WE0002-nckuh")
    dry = R.sync_cases(store, root, prune=True, dry_run=True, roster={})
    assert dry["removed"] == ["24WE0002-nckuh"] and store.patient_count() == 1
    R.sync_cases(store, root, prune=True, roster={})
    assert store.patient_count() == 0
    store.close()


def test_roster_mrn_merges_unregistered_cases(tmp_path):
    root = tmp_path / "tertiary_output"
    make_case(root, "24WE0001-nckuh", [SHARED])
    make_case(root, "24WE0777-nckuh", [SHARED], source={"seq_type": "WES", "source_sample_id": "24WE0777"})
    store = R.WesRecurrence(tmp_path / "rec.sqlite")
    roster = {"24WE0001": {"mrn": "P9"}, "24WE0777": {"mrn": "P9"}}
    summary = R.sync_cases(store, root, roster=roster)
    assert summary["patients"] == 1
    store.close()


def test_engine_reads_recurrence_store(tmp_path, sources):
    _, db = build(tmp_path, sources)
    root = tmp_path / "tertiary_output"
    for i in range(120):
        make_case(root, f"24WE{i:04d}-nckuh", [SHARED] if i < 10 else [("chr5", i + 1, "A", "C")])
    store = R.WesRecurrence(tmp_path / "rec.sqlite")
    R.sync_cases(store, root, roster={})
    ctx = Context(test_type="WES", refs=S.SqliteReferences(db, recurrence=store))
    result = evaluate_variant(variant(), ctx)
    assert result["flags"] and result["flags"][0]["code"] == "recurrent"     # 10/120 ≥ 5%
    rare = evaluate_variant(variant(POS="1001", REF="G", ALT="A"), ctx)
    assert not rare["flags"]
    assert ctx.refs.version["wes_recurrence"]["patients"] == 120
    store.close()


# ── CLIs ───────────────────────────────────────────────────────────────────
def _run(script, *args, env=None):
    full = dict(os.environ, PYTHONPATH=str(REPO / "backend"), **(env or {}))
    return subprocess.run([sys.executable, str(REPO / "scripts" / "acmg_engine" / script), *map(str, args)],
                          capture_output=True, text=True, env=full, check=True)


def test_build_references_cli(tmp_path, sources):
    out = tmp_path / "out"
    proc = _run("build_references.py", "--clinvar", sources["clinvar_summary"], "--mane-gtf", sources["mane_gtf"],
                "--mane-summary", sources["mane_summary"], "--clingen-hi", sources["clingen_hi"], "--out", out)
    report = json.loads(proc.stdout)
    assert report["counts"]["truncating"] == 4
    assert (out / S.REFS_DB_NAME).is_file()


def test_update_wes_recurrence_cli(tmp_path):
    root = tmp_path / "tertiary_output"
    make_case(root, "24WE0001-nckuh", [SHARED])
    make_case(root, "23WE0044T-nckuh", [SHARED])
    db = tmp_path / "rec.sqlite"
    env = {"NGS_UI_PATIENT_LIST_DIR": str(tmp_path / "no_roster")}
    dry = json.loads(_run("update_wes_recurrence.py", "--all", "--dry-run", "--root", root, "--db", db,
                          env=env).stdout)
    assert dry["added"] == 1 and dry["skipped"] == {"測試檢體（T）": 1}
    report = json.loads(_run("update_wes_recurrence.py", "--all", "--root", root, "--db", db, env=env).stdout)
    assert report["added"] == 1 and report["patients"] == 1 and report["roster_entries"] == 0
    removed = json.loads(_run("update_wes_recurrence.py", "--remove", "24WE0001-nckuh", "--db", db,
                              env=env).stdout)
    assert removed["removed"] == ["24WE0001-nckuh"] and removed["patients"] == 0
