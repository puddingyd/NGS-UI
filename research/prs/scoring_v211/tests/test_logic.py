import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pipeline
import qc_report
import compare_benchmarks
import discover_gvcfs


class LogicTests(unittest.TestCase):
    def test_group_boundaries(self):
        self.assertEqual(qc_report.group(35.6462585), "20-40")
        self.assertEqual(qc_report.group(60.6803), "60-80")
        self.assertEqual(qc_report.group(79.5918), "60-80")
        self.assertEqual(qc_report.group(80), "80-90")
        self.assertEqual(qc_report.group(95), ">=95")

    def test_benchmark_compare(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            cfg=json.loads((Path(__file__).parents[1]/"benchmarks.json").read_text())
            for sid, expected in cfg["samples"].items():
                q=p/sid/"qc/qc.json"
                q.parent.mkdir(parents=True)
                q.write_text(json.dumps({**expected,"qc_status":"PASS_WITH_REVIEW_FLAGS"}))
            rows=compare_benchmarks.check(p,Path(__file__).parents[1]/"benchmarks.json")
            self.assertTrue(all(row[-1]=="PASS" for row in rows))
            sample=p/"26T00969/qc/qc.json"
            data=json.loads(sample.read_text());data["normalized_sum"]+=0.01
            sample.write_text(json.dumps(data))
            rows=compare_benchmarks.check(p,Path(__file__).parents[1]/"benchmarks.json")
            self.assertTrue(any(row[-1]=="FAIL" and row[1]=="normalized_sum" for row in rows))

    def test_derive_real_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"raw"/"audit"
            p.mkdir(parents=True)
            (p/"source_sha256.txt").write_text("a"*64+"  /real/file/26T00970.hard-filtered.gvcf.gz\n")
            self.assertEqual(discover_gvcfs.source_for(p.parent),"/real/file/26T00970.hard-filtered.gvcf.gz")

    def test_sample_id_guard(self):
        self.assertTrue(pipeline.SID.fullmatch("26T00970"))
        self.assertFalse(pipeline.SID.fullmatch("../oops"))
        self.assertFalse(pipeline.SID.fullmatch("sample;rm -rf"))

    def test_raw_integrity_fixture(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/"audit").mkdir()
            (p/"audit/final_reconciliation.json").write_text(json.dumps({"status":"PASS","iid":"26T00970"}))
            (p/"scorefiles").mkdir()
            (p/"scorefiles/PGS004879_hmPOS_GRCh38.txt.gz").write_text("fixture")
            (p/"roundtrip_22chr/pgen").mkdir(parents=True)
            prefix=p/"roundtrip_22chr/pgen/chr1"
            for ext in ("pgen","psam","pvar"):
                Path(str(prefix)+"."+ext).write_text("fixture")
            with (p/"samplesheet.csv").open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=["sampleset","path_prefix","chrom","format"])
                w.writeheader();w.writerow({"sampleset":"case26T00970","path_prefix":str(prefix),"chrom":"1","format":"pfile"})
            (p/"results/case26T00970/score").mkdir(parents=True)
            (p/"results/case26T00970/score/aggregated_scores.txt.gz").write_text("fixture")
            result=pipeline.raw_integrity(p,"26T00970")
            self.assertEqual(result["iid"],"26T00970")
            with self.assertRaises(RuntimeError): pipeline.raw_integrity(p,"26T00969")

    def test_qc_fixture(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);sample="26T00970"
            raw=p/"raw";anc=p/"ancestry"
            (raw/"audit").mkdir(parents=True)
            (raw/"audit/final_reconciliation.json").write_text(json.dumps({"status":"PASS","iid":sample,"observed_sum":"0.351907","observed_denom":1218692}))
            (raw/"audit/upstream_reconstruction_summary.json").write_text(json.dumps({"matched":609346}))
            score=anc/"results"/f"case{sample}"/"score";score.mkdir(parents=True)
            with gzip.open(score/f"case{sample}_popsimilarity.txt.gz","wt") as f:
                w=csv.DictWriter(f,fieldnames=["sampleset","IID","REFERENCE","SuperPop","MostSimilarPop","MostSimilarPop_LowConfidence","RF_P_EAS"],delimiter="\t")
                w.writeheader()
                w.writerow({"sampleset":f"case{sample}","IID":sample,"REFERENCE":"False","MostSimilarPop":"EAS","MostSimilarPop_LowConfidence":"False","RF_P_EAS":"1"})
                for i in range(735):w.writerow({"sampleset":"reference","IID":f"REF{i:04d}","REFERENCE":"True","SuperPop":"EAS"})
            with gzip.open(score/f"case{sample}_pgs.txt.gz","wt") as f:
                w=csv.DictWriter(f,fieldnames=["sampleset","IID","PGS","SUM","Z_MostSimilarPop","Z_norm1","Z_norm2","percentile_MostSimilarPop"],delimiter="\t")
                w.writeheader();w.writerow({"sampleset":f"case{sample}","IID":sample,"PGS":"PGS004879_hmPOS_GRCh38","SUM":"0.354661","Z_MostSimilarPop":"-0.3445","Z_norm1":"-0.3608","Z_norm2":"-0.3924","percentile_MostSimilarPop":"35.6462585"})
                for i in range(735):
                    w.writerow({"sampleset":"reference","IID":f"REF{i:04d}","SUM":f"{i/735:.7f}"})
            match=score.parent/"match";match.mkdir()
            with (match/f"case{sample}_summary.csv").open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=["match_status","count"]);w.writeheader()
                for k,v in (("matched",609032),("excluded",314),("unmatched",1331)):w.writerow({"match_status":k,"count":v})
            with gzip.open(match/f"case{sample}_log.csv.gz","wt") as f:
                w=csv.DictWriter(f,fieldnames=["chr_name","chr_position","effect_allele","other_allele","effect_weight","match_status","exclude","match_type","row_nr"])
                w.writeheader();w.writerow({"chr_name":"9","chr_position":"22081398","effect_allele":"G","other_allele":"T","effect_weight":"0.01395","match_status":"unmatched","row_nr":"339539"})
            markers=anc/"work/ancestry/filter";markers.mkdir(parents=True)
            with gzip.open(markers/"GRCh38_HGDP+1kGP_ALL_thinned.prune.in.gz","wt") as f:f.write("1\n2\n")
            result,missing=qc_report.calc(sample,raw,anc,"checksum","ref.tar.zst",.995,.001,2)
            self.assertEqual(result["matched"],609032)
            self.assertEqual(result["raw_sum"],.351907)
            self.assertEqual(result["intervene_group"],"20-40")
            self.assertEqual(result["high_weight_unmatched_n"],1)
            self.assertEqual(missing[0]["chr_position"],"22081398")
            self.assertEqual(result["qc_status"],"PASS_WITH_REVIEW_FLAGS")


if __name__=="__main__":unittest.main()
