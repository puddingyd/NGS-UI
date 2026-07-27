import gzip

from app.services.ploidy import load_sample_ploidy, read_karyotype
from app.workers.dragen_run import (
    _copy_dragen_ploidy_vcf,
    _copy_ploidy_artifacts,
    _find_dragen_ploidy_vcf,
    _find_nckuh_ploidy_vcf,
)


def _write_ploidy(path, call):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(f"##estimatedSexKaryotype={call}\n#CHROM\tPOS\n")


def test_find_and_copy_dragen_ploidy_by_source_vcf_basename(tmp_path):
    source_vcf = tmp_path / "VAL-3.hard-filtered.vcf.gz"
    source_vcf.touch()
    expected = tmp_path / "VAL-3.ploidy.vcf.gz"
    _write_ploidy(expected, "XY")
    destination_dir = tmp_path / "tertiary_output" / "LIS-001"
    destination_dir.mkdir(parents=True)

    assert _find_dragen_ploidy_vcf(source_vcf, "VAL-3") == expected
    copied = _copy_dragen_ploidy_vcf(source_vcf, "VAL-3", destination_dir)

    destination = destination_dir / "LIS-001.ploidy.vcf.gz"
    assert copied == (expected, destination)
    assert read_karyotype(destination) == "XY"


def test_find_dragen_ploidy_never_uses_another_sample(tmp_path):
    source_vcf = tmp_path / "VAL-3.hard-filtered.vcf.gz"
    source_vcf.touch()
    _write_ploidy(tmp_path / "VAL-4.ploidy.vcf.gz", "XX")

    assert _find_dragen_ploidy_vcf(source_vcf, "VAL-3") is None


def test_missing_dragen_ploidy_removes_stale_ui_sidecar(tmp_path):
    source_vcf = tmp_path / "VAL-3.hard-filtered.vcf.gz"
    source_vcf.touch()
    destination_dir = tmp_path / "tertiary_output" / "LIS-001"
    destination_dir.mkdir(parents=True)
    stale = destination_dir / "LIS-001.ploidy.vcf.gz"
    _write_ploidy(stale, "XX")

    assert _copy_dragen_ploidy_vcf(source_vcf, "VAL-3", destination_dir) is None
    assert not stale.exists()


def test_load_sample_ploidy_preserves_x_and_xxy(tmp_path):
    _write_ploidy(tmp_path / "ploidy.vcf.gz", "XXY")
    result = load_sample_ploidy(tmp_path)
    assert result["karyotype"] == "XXY"
    assert result["source"] == "ploidy.vcf.gz"
    assert result["aneuploidy_suspected"] is True
    assert result["karyotype_interpretation"] == "possible 47,XXY"

    _write_ploidy(tmp_path / "ploidy.vcf.gz", "X")
    assert read_karyotype(tmp_path / "ploidy.vcf.gz") == "X"


def test_nckuh_ploidy_is_copied_from_alignment_qc_with_qc_text(tmp_path):
    sample_root = tmp_path / "26WE0049"
    source_vcf = sample_root / "04_snv_indel" / "26WE0049.ensemble.fixed.vcf.gz"
    source_vcf.parent.mkdir(parents=True)
    source_vcf.touch()
    ploidy_vcf = sample_root / "03_alignment_qc" / "26WE0049.ploidy.vcf.gz"
    ploidy_vcf.parent.mkdir(parents=True)
    _write_ploidy(ploidy_vcf, "XY")
    qc_txt = ploidy_vcf.parent / "26WE0049.ploidy_qc.txt"
    qc_txt.write_text("future use only\n", encoding="utf-8")
    post = tmp_path / "tertiary_output" / "LIS-001" / "08_postprocessing"
    post.mkdir(parents=True)

    assert _find_nckuh_ploidy_vcf(source_vcf, "26WE0049") == ploidy_vcf
    copied_vcf, copied_qc = _copy_ploidy_artifacts(
        mode="inhouse",
        source_vcf=source_vcf,
        source_sample_id="26WE0049",
        sample_id="LIS-001",
        post_dir=post,
    )

    assert copied_vcf == (ploidy_vcf, post / "LIS-001.ploidy.vcf.gz")
    assert copied_qc == (qc_txt, post / "LIS-001.ploidy_qc.txt")


def test_ploidy_parser_flags_non_pass_nuclear_contig(tmp_path):
    sample_dir = tmp_path / "S1"
    sample_dir.mkdir()
    path = sample_dir / "S1.ploidy.vcf.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(
            "##source=NCKUH_PLOIDY_MOSDEPTH\n"
            "##seqType=WES\n"
            "##estimatedSexKaryotype=XY\n"
            "##referenceSexKaryotype=unknown\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t1\t.\tN\t.\t.\tPASS\tEND=10\tDC:NDC:RATIO\t30:1.0:1.0\n"
            "chrY\t1\t.\tN\t.\t.\tSUSPECT\tEND=10\tDC:NDC:RATIO\t15:1.6:0.8\n"
            "chrM\t1\t.\tN\t.\t.\tSUSPECT\tEND=10\tDC:NDC:RATIO\t4:.:0.2\n"
        )
    result = load_sample_ploidy(sample_dir)

    assert result["karyotype"] == "XY"
    assert result["pipeline_source"] == "NCKUH_PLOIDY_MOSDEPTH"
    assert result["seq_type"] == "WES"
    assert len(result["chromosomes"]) == 3
    assert result["chromosomes"][0]["end"] == 10
    assert [row["chrom"] for row in result["warnings"]] == ["chrY"]
    assert result["warnings"][0]["dosage_call"] == "gain"
    assert result["warnings"][0]["call_label"] == "Gain signal"
    assert result["warnings"][0]["confidence"] == "suspect"
    assert result["warnings"][0]["observed_ratio"] == 0.8
    assert result["aneuploidy_suspected"] is True


def test_dragen_pass_dup_is_an_aneuploidy_signal_and_ratio_is_derived(tmp_path):
    path = tmp_path / "ploidy.vcf.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(
            "##source=DRAGEN_PLOIDY\n"
            "##autosomeDepthOfCoverage=60\n"
            "##estimatedSexKaryotype=XY\n"
            "##referenceSexKaryotype=XY\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t1\t.\tN\t.\t30\tPASS\tEND=10\tDC:NDC\t60:1.0\n"
            "chr21\t1\t.\tN\t<DUP>\t35\tPASS\tEND=10;SVTYPE=DUP\tDC:NDC\t90:1.5\n"
            "chrX\t1\t.\tN\t.\t30\tPASS\tEND=10\tDC:NDC\t30:1.0\n"
        )

    result = load_sample_ploidy(tmp_path)
    chr21 = result["abnormal_chromosomes"][0]
    chrx = result["chromosomes"][2]

    assert result["pipeline_kind"] == "dragen"
    assert result["autosome_depth"] == 60
    assert result["aneuploidy_suspected"] is True
    assert chr21["chrom"] == "chr21"
    assert chr21["dosage_call"] == "gain"
    assert chr21["call_label"] == "Gain"
    assert chr21["confidence"] == "pass"
    assert chr21["interpretation"] == "possible trisomy 21"
    assert chr21["observed_ratio"] == 1.5
    assert chr21["ratio_source"] == "derived"
    assert chrx["expected_ratio"] == 0.5
    assert chrx["observed_ratio"] == 0.5


def test_dragen_lowqual_call_is_retained_but_lowqual_normal_is_qc_only(tmp_path):
    path = tmp_path / "ploidy.vcf.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(
            "##source=DRAGEN_PLOIDY\n"
            "##estimatedSexKaryotype=XX\n"
            "##referenceSexKaryotype=XX\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr18\t1\t.\tN\t<DEL>\t8\tLowQual\tEND=10;SVTYPE=DEL\tDC:NDC\t15:0.5\n"
            "chr20\t1\t.\tN\t.\t8\tLowQual\tEND=10\tDC:NDC\t30:1.0\n"
        )

    result = load_sample_ploidy(tmp_path)

    assert [row["chrom"] for row in result["abnormal_chromosomes"]] == ["chr18"]
    assert result["abnormal_chromosomes"][0]["call_label"] == "Loss signal"
    assert result["abnormal_chromosomes"][0]["confidence"] == "low"
    assert [row["chrom"] for row in result["qc_warnings"]] == ["chr18", "chr20"]


def test_missing_karyotype_is_not_mislabeled_as_aneuploidy(tmp_path):
    path = tmp_path / "ploidy.vcf.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1\t.\tN\t.\t.\tPASS\tEND=10\n"
        )

    result = load_sample_ploidy(tmp_path)

    assert result["exists"] is True
    assert result["karyotype"] == ""
    assert result["aneuploidy_suspected"] is False
