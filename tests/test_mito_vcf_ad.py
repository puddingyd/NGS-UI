"""Mito card AD comes from the paired VEP VCF when the TSV omits it."""
import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from app.adapters import mito_tsv
from app.services import sample_layout


class MitoVcfAdTest(unittest.TestCase):
    def test_direct_tsv_load_finds_single_sample_vcf_sibling(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tsv = root / "sample.mito.tsv"
            tsv.write_text(
                "CHROM\tPOS\tREF\tALT\tGENE\tDP\n"
                "chrM\t73\tA\tG\tMT-TF\t10\n",
                encoding="utf-8",
            )
            with gzip.open(root / "sample.mito.vep.vcf.gz", "wt", encoding="utf-8") as handle:
                handle.write(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
                    "chrM\t73\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:4,6\n"
                )
            with patch.object(mito_tsv.clinvar_mito, "lookup", return_value={}), \
                 patch.object(mito_tsv.inhouse_af_mito, "lookup", return_value={}), \
                 patch.object(mito_tsv.panel_deadzone, "canonical_panel_gene_symbol", side_effect=lambda gene: gene):
                variants, _ = mito_tsv.load_mito_tsv(tsv)
            self.assertEqual(variants["chrM-73-A-G"]["AD"], "4,6")

    def test_exact_allele_and_sample_join_with_multiallelic_depths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tsv = root / "sample.mito.tsv"
            tsv.write_text(
                "CHROM\tPOS\tREF\tALT\tGENE\tDP\tAD\n"
                "chrM\t73\tA\tG\tMT-TF\t10\t\n"
                "chrM\t100\tA\tG\tMT-TF\t21\t\n"
                "chrM\t100\tA\tT\tMT-TF\t21\t\n"
                "chrM\t200\tC\tT\tMT-TF\t9\t1,8\n"
                "chrM\t300\tG\tA\tMT-TF\t8\t\n",
                encoding="utf-8",
            )
            vcf = root / "sample.mito.vep.vcf.gz"
            with gzip.open(vcf, "wt", encoding="utf-8") as handle:
                handle.write(
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tother\tsample\n"
                    "chrM\t73\t.\tA\tG\t.\tPASS\t.\tGT:DP:AD\t0/1:10:9,1\t0/1:10:4,6\n"
                    "chrM\t100\t.\tA\tG,T\t.\tPASS\t.\tGT:AD:DP\t0/1:20,1,0:21\t1/2:5,7,9:21\n"
                    "chrM\t200\t.\tC\tT\t.\tPASS\t.\tGT:AD\t0/1:8,1\t0/1:2,7\n"
                )
            with patch.object(mito_tsv.clinvar_mito, "lookup", return_value={}), \
                 patch.object(mito_tsv.inhouse_af_mito, "lookup", return_value={}), \
                 patch.object(mito_tsv.panel_deadzone, "canonical_panel_gene_symbol", side_effect=lambda gene: gene):
                variants, categories = mito_tsv.load_mito_tsv(
                    tsv, vcf_path=vcf, vcf_sample_id="sample",
                )
                wrong, _ = mito_tsv.load_mito_tsv(
                    tsv, vcf_path=vcf, vcf_sample_id="absent",
                )
            self.assertEqual(variants["chrM-73-A-G"]["AD"], "4,6")
            self.assertEqual(variants["chrM-100-A-G"]["AD"], "5,7")
            self.assertEqual(variants["chrM-100-A-T"]["AD"], "5,9")
            self.assertEqual(variants["chrM-200-C-T"]["AD"], "1,8")
            self.assertEqual(variants["chrM-300-G-A"]["AD"], "")
            self.assertEqual(wrong["chrM-73-A-G"]["AD"], "")
            self.assertEqual(sum(map(len, categories.values())), 5)

    def test_vcf_cache_refreshes_when_file_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            vcf = Path(directory) / "sample.mito.vep.vcf.gz"
            def write(ad):
                with gzip.open(vcf, "wt", encoding="utf-8") as handle:
                    handle.write(
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
                        f"chrM\t73\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:{ad}\n"
                    )
            write("4,6")
            self.assertEqual(mito_tsv._vcf_ad_by_allele(vcf, "sample")[(73, "A", "G")], "4,6")
            write("40,60")
            self.assertEqual(mito_tsv._vcf_ad_by_allele(vcf, "sample")[(73, "A", "G")], "40,60")

    def test_derived_mitomap_tsv_still_uses_paired_pipeline_vcf(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sample_dir = root / "case"
            pipeline_dir = sample_dir / "04_mito"
            state_dir = sample_dir / "08_postprocessing"
            pipeline_dir.mkdir(parents=True)
            state_dir.mkdir()
            (state_dir / "case.layout.json").write_text(
                json.dumps({"layout_version": 3}), encoding="utf-8",
            )
            (state_dir / "case.pipeline_source.json").write_text(
                json.dumps({"source_sample_id": "source"}), encoding="utf-8",
            )
            (state_dir / "case.mito.annotated.tsv").write_text("", encoding="utf-8")
            (pipeline_dir / "source.mito.tsv").write_text("", encoding="utf-8")
            paired = pipeline_dir / "source.mito.vep.vcf.gz"
            paired.touch()
            with patch.object(sample_layout.config, "PIPELINE_OUT_ROOT", root), \
                 patch.object(sample_layout.config, "LEGACY_TERTIARY_OUTPUT_ROOT", root / "old-ui"), \
                 patch.object(sample_layout.config, "LEGACY_PIPELINE_OUT_ROOT", root / "old-pipeline"):
                self.assertEqual(sample_layout.mito_tsv("case"), state_dir / "case.mito.annotated.tsv")
                self.assertEqual(sample_layout.mito_vep_vcf("case"), paired)


if __name__ == "__main__":
    unittest.main()
