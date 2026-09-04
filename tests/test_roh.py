from pathlib import Path

from app import config
from app.services import roh


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _inhouse_source(root: Path, sample: str = "SRC") -> Path:
    return _write(root / "04_snv_indel" / f"{sample}.ensemble.fixed.vcf.gz", "vcf")


def test_inhouse_prefers_existing_automap_without_research_flag(tmp_path):
    source_vcf = _inhouse_source(tmp_path)
    _write(
        tmp_path / "08_roh" / "SRC.HomRegions.tsv",
        "#Chr\tBegin\tEnd\tSize(Mb)\tNb_variants\tPercentage_homozygosity\n"
        "chr1\t101\t2000100\t2.0\t30\t95.0\n",
    )
    _write(
        tmp_path / "08_roh" / "SRC.roh.txt",
        "RG\tSRC\tchr2\t1\t2000000\t2000000\t40\t99\n",
    )
    post = tmp_path / "stage" / "08_postprocessing"

    summary = roh.prepare_roh_outputs(
        mode="inhouse", source_vcf=source_vcf, source_sample_id="SRC",
        sample_id="LIS", post_dir=post, research_only=False,
    )

    assert summary["source"] == "automap"
    assert summary["generated_automap"] is False
    assert summary["region_count_default"] == 1
    assert (post / "LIS.roh.source.automap.tsv").is_file()
    assert (post / "LIS.roh.source.bcftools.txt").is_file()


def test_inhouse_research_only_runs_automap_then_uses_it(tmp_path):
    source_vcf = _inhouse_source(tmp_path)
    _write(tmp_path / "04_snv_indel" / "SRC.haplotypecaller.vcf.gz", "hc")
    _write(
        tmp_path / "08_roh" / "SRC.roh.txt",
        "RG\tSRC\tchr2\t1\t2000000\t2000000\t40\t99\n",
    )
    called = []

    def runner(vcf, sample):
        called.append((vcf, sample))
        output = _write(
            tmp_path / "generated" / "SRC.HomRegions.tsv",
            "#Chr\tBegin\tEnd\tSize(Mb)\tNb_variants\tPercentage_homozygosity\n"
            "chr3\t1\t3000000\t3.0\t50\t97.0\n",
        )
        return output, None

    summary = roh.prepare_roh_outputs(
        mode="inhouse", source_vcf=source_vcf, source_sample_id="SRC",
        sample_id="LIS", post_dir=tmp_path / "post", research_only=True,
        automap_runner=runner,
    )

    assert called == [(tmp_path / "04_snv_indel" / "SRC.haplotypecaller.vcf.gz", "SRC")]
    assert summary["source"] == "automap"
    assert summary["generated_automap"] is True


def test_inhouse_without_automap_uses_bcftools_default_confidence_filter(tmp_path):
    source_vcf = _inhouse_source(tmp_path)
    _write(
        tmp_path / "08_roh" / "SRC.roh.txt",
        "# This file was produced by: bcftools roh(1.23.1)\n"
        "RG\tSRC\tchr1\t1\t2000000\t2000000\t30\t25\n"
        "RG\tSRC\tchr2\t1\t2000000\t2000000\t20\t25\n"
        "RG\tSRC\tchrM\t1\t2000000\t2000000\t30\t25\n",
    )
    post = tmp_path / "post"

    summary = roh.prepare_roh_outputs(
        mode="inhouse", source_vcf=source_vcf, source_sample_id="SRC",
        sample_id="LIS", post_dir=post, research_only=False,
    )

    assert summary["source"] == "bcftools"
    assert summary["region_count_all"] == 2  # chrM is intentionally excluded.
    assert summary["region_count_default"] == 1
    assert summary["autosomal_total_mb_default"] == 2.0


def test_failed_research_automap_is_nonfatal_and_falls_back_to_bcftools(tmp_path):
    source_vcf = _inhouse_source(tmp_path)
    _write(tmp_path / "04_snv_indel" / "SRC.haplotypecaller.vcf.gz", "hc")
    _write(
        tmp_path / "08_roh" / "SRC.roh.txt",
        "RG\tSRC\tchr1\t1\t2000000\t2000000\t30\t25\n",
    )

    def failed_runner(*_args):
        raise RuntimeError("container unavailable")

    summary = roh.prepare_roh_outputs(
        mode="inhouse", source_vcf=source_vcf, source_sample_id="SRC",
        sample_id="LIS", post_dir=tmp_path / "post", research_only=True,
        automap_runner=failed_runner,
    )

    assert summary["source"] == "bcftools"
    assert "AutoMap 執行失敗" in summary["warnings"][0]


def test_dragen_always_uses_native_roh_and_large_region_threshold(tmp_path):
    source_vcf = _write(tmp_path / "run" / "vcf.gz" / "SRC.hard-filtered.vcf.gz", "vcf")
    germline = tmp_path / "run" / "other" / "SRC" / "germline_seq"
    _write(
        germline / "SRC.roh.bed",
        "chr1\t0\t2999999\t1.0\t100\t2\n"
        "chr5\t100\t3100100\t19.88\t1419\t16\n",
    )
    _write(
        germline / "SRC.roh_metrics.csv",
        "VARIANT CALLER,,Percent SNVs in large ROH ( >= 3000000),0.070\n"
        "VARIANT CALLER,,Number of large ROH ( >= 3000000),1\n",
    )
    calls = []

    summary = roh.prepare_roh_outputs(
        mode="dragen", source_vcf=source_vcf, source_sample_id="SRC",
        sample_id="LIS", post_dir=tmp_path / "post", research_only=True,
        automap_runner=lambda *_: calls.append(True),
    )

    assert calls == []
    assert summary["source"] == "dragen"
    assert summary["region_count_all"] == 2
    assert summary["region_count_default"] == 1
    assert summary["source_details"]["metrics"]["Percent SNVs in large ROH ( >= 3000000)"] == 0.07


def test_variant_in_roh_join_uses_half_open_boundaries(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "PIPELINE_OUT_ROOT", tmp_path)
    post = tmp_path / "LIS" / "08_postprocessing"
    post.mkdir(parents=True)
    _write(
        post / "LIS.roh_regions.tsv",
        "region_id\tchrom\tstart0\tend0\tdisplay_start\tdisplay_end\tlength_bp\tlength_mb\tsource\tpasses_default_filter\tn_markers\tn_hom\tn_het\thomozygosity_pct\tquality\tscore\n"
        "roh-1\tchr1\t100\t200\t101\t200\t100\t0.0001\tautomap\t1\t30\t\t\t95\t\t\n",
    )
    variants = {
        "start": {"CHROM": "1", "POS": 101},
        "end": {"CHROM": "chr1", "POS": 200},
        "outside": {"CHROM": "chr1", "POS": 201},
    }

    roh.annotate_variants(variants, "LIS")

    assert variants["start"]["in_roh"] is True
    assert variants["end"]["in_roh"] is True
    assert variants["outside"]["in_roh"] is False
