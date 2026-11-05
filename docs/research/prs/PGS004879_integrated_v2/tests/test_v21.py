"""Synthetic regression for v2.1 orchestration. No medical/genomic claims."""
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pipeline_v21 as run
import v21_raw_resume as rr


class V21Tests(unittest.TestCase):
    def test_baseline_manifest(self):
        self.assertEqual(len(run.baseline()),10)

    def test_refuse_mutated_baseline(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'pipeline.py').write_text('different')
            (root/'BASELINE_V2_SHA256SUMS.txt').write_text('0'*64+'  pipeline.py\n')
            with patch.object(run,'ROOT',root),patch.object(run,'BASELINE',root/'BASELINE_V2_SHA256SUMS.txt'):
                with self.assertRaisesRegex(RuntimeError,'changed'):
                    run.baseline()

    def test_sample_lock_prevents_second_writer(self):
        with tempfile.TemporaryDirectory() as tmp:
            with run.guard_sample(tmp,'26T00970'):
                with self.assertRaisesRegex(RuntimeError,'already running'):
                    with run.guard_sample(tmp,'26T00970'):
                        pass

    def test_resumable_generated_script_matches_vendor(self):
        vendor=Path(run.ROOT)/'vendor/PGS004879_FINAL_gVCF_resume.sh'
        code=rr.generated_verify_script(vendor)
        self.assertIn(b'def roundtrip():',code)
        self.assertIn(b'def reconcile():',code)
        self.assertNotIn(b'PY_VERIFY',code)

    def test_early_failure_is_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            raw=Path(tmp)
            with self.assertRaises(rr.EarlyStageFailure):
                rr.verify_existing_roundtrip(raw,'26T00970','/abs/nonexistent',{})

    def test_complete_result_requires_correct_counts_and_lock(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)
            m=dest/'method2';m.mkdir()
            for n in ('status.json','summary.json','scores.tsv','missing_loci.tsv.gz','report.md','inputs.lock.json'):
                (m/n).write_text('fixture')
            (m/'status.json').write_text(json.dumps({'status':'PASS'}))
            (m/'summary.json').write_text(json.dumps({'sample':'26T00970','stage':'METHOD2_TECHNICAL_PASS','model_lock_sha256':'locked'}))
            with (m/'scores.tsv').open('w') as f:
                w=csv.DictWriter(f,fieldnames=['sample','model','n_loci','observed_loci','missing_loci'],delimiter='\t')
                w.writeheader()
                w.writerow({'sample':'26T00970','model':'A','n_loci':10,'observed_loci':9,'missing_loci':1})
                w.writerow({'sample':'26T00970','model':'B','n_loci':7,'observed_loci':6,'missing_loci':1})
            self.assertTrue(run.completed_outputs(dest,'26T00970','locked'))
            self.assertFalse(run.completed_outputs(dest,'26T00970','wrong'))
            self.assertFalse(run.completed_outputs(dest,'26T00969','locked'))
            (m/'scores.tsv').write_text('sample\tmodel\tn_loci\tobserved_loci\tmissing_loci\n26T00970\tA\t10\t9\t0\n')
            self.assertFalse(run.completed_outputs(dest,'26T00970','locked'))

    def test_summarize_partial_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)/'integrated_runs'/'r001'/'26T00970'
            (d/'orchestration').mkdir(parents=True)
            (d/'orchestration/v21_status.json').write_text(json.dumps({'stage':'FAILED','message':'simulated failure'}))
            (d/'orchestration/v21_inputs.lock.json').write_text('{}')
            result=run.summarize(tmp,'r001')
            self.assertEqual(len(result),1)
            self.assertEqual(result[0]['execution_status'],'FAILED')
            self.assertTrue((d.parent/'batch_summary_v21.tsv').is_file())

    def test_manifest_duplicate_sample_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'list.tsv'
            p.write_text('sample\tmode\tgvcf\traw_source\n26T00970\tfresh\t/a.gz\t\n26T00970\tfresh\t/b.gz\t\n')
            from argparse import Namespace
            args=Namespace(config='/unused',manifest=str(p),run_id='r001',jobs=1,resume=False,with_legacy_ancestry=False)
            with patch.object(run.v1,'cfg_read',return_value={'project':tmp}):
                with self.assertRaisesRegex(RuntimeError,'Duplicate batch'):
                    run.batch(args)

    def test_output_hashes_change_detection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            path=root/'method2/scores.tsv';path.parent.mkdir()
            path.write_text('pass')
            earlier=run.output_fingerprints(root)
            path.write_text('modified')
            self.assertNotEqual(earlier,run.output_fingerprints(root))


if __name__=='__main__':unittest.main()

# Note: module-level pytest/unittest discovery picks up this class as well.
class V21RunnerMocks(unittest.TestCase):
    def test_new_run_delegates_to_original_v2_without_touching_science(self):
        from argparse import Namespace
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            gvcf=root/'original.gvcf.gz';gvcf.write_text('fake')
            lock=root/'frozen.lock';lock.write_text('frozen')
            dest=root/'integrated_runs'/'run001'/'sample1'
            args=Namespace(config=str(root/'config.json'),sample='sample1',run_id='run001',
                           mode='fresh',gvcf=str(gvcf),raw_source=None,
                           with_legacy_ancestry=False,resume=False)
            cfg={'project':str(root),'method2_lock':str(lock),
                 'method2_model_dir':str(root),'original_scorefile':str(root/'source.gz')}
            def simulated_v2(cmd,logfile):
                self.assertIn('pipeline_v2.py',cmd[1])
                self.assertIn('fresh',cmd)
                md=dest/'method2';md.mkdir(parents=True)
                (md/'status.json').write_text(json.dumps({'status':'PASS'}))
                (md/'summary.json').write_text(json.dumps({'sample':'sample1',
                          'stage':'METHOD2_TECHNICAL_PASS',
                          'model_lock_sha256':run.m2.sha256(lock)}))
                with (md/'scores.tsv').open('w') as f:
                    w=csv.DictWriter(f,fieldnames=['sample','model','n_loci','observed_loci','missing_loci'],delimiter='\t')
                    w.writeheader()
                    for m in 'AB': w.writerow({'sample':'sample1','model':m,'n_loci':10,'observed_loci':9,'missing_loci':1})
                for n in ('missing_loci.tsv.gz','report.md','inputs.lock.json'):
                    (md/n).write_text('fixture')
            with patch.object(run.v1,'cfg_read',return_value=cfg), \
                 patch.object(run.m2,'validate_locked_model'), \
                 patch.object(run,'fingerprint',return_value={'immutable':'yes'}), \
                 patch.object(run.v1,'log_run',side_effect=simulated_v2):
                run.run_one(args)
            self.assertEqual(json.loads((dest/'orchestration/v21_status.json').read_text())['stage'],'PASS')
            self.assertTrue((dest/'orchestration/v21_output_sha256.json').is_file())
            # The successful result must not be overwritten by a new run.
            with patch.object(run.v1,'cfg_read',return_value=cfg), \
                 patch.object(run.m2,'validate_locked_model'), \
                 patch.object(run,'fingerprint',return_value={'immutable':'yes'}), \
                 patch.object(run.v1,'log_run',side_effect=AssertionError('should not rerun')), \
                 patch.object(run.v2,'verify_fresh_raw_origin'):
                args.resume=True
                run.run_one(args)
            # Modifying a previously completed score must fail closed.
            (dest/'method2/scores.tsv').write_text('tampered')
            with patch.object(run.v1,'cfg_read',return_value=cfg), \
                 patch.object(run.m2,'validate_locked_model'), \
                 patch.object(run,'fingerprint',return_value={'immutable':'yes'}):
                with self.assertRaisesRegex(RuntimeError,'content changed|invalid'):
                    run.run_one(args)

class V21ResumeMocks(unittest.TestCase):
    def test_resume_nextflow_uses_exact_existing_work_tree(self):
        with tempfile.TemporaryDirectory() as tmp:
            raw=Path(tmp)
            (raw/'prep').mkdir();(raw/'scorefiles').mkdir();(raw/'audit').mkdir()
            (raw/'prep/score_ready.vcf.gz').write_bytes(b'a')
            (raw/'scorefiles/PGS004879_hmPOS_GRCh38.txt.gz').write_bytes(b'b')
            first=[True]
            def reconciliation(cmd,logfile,cwd,env):
                if first[0]:
                    first[0]=False
                    raise RuntimeError('simulated incomplete Nextflow score')
            commands=[]
            def nextflow(cmd,logfile,cwd):
                commands.append(cmd)
            with patch.object(rr,'verify_existing_roundtrip',return_value={'OUT':str(raw)}), \
                 patch.object(rr,'run_logged_env',side_effect=reconciliation), \
                 patch.object(rr.v1,'log_run',side_effect=nextflow), \
                 patch.object(rr.v1,'raw_integrity',return_value={'status':'PASS'}):
                rr.resume_raw_after_roundtrip(raw,'s1','/unused',{})
            self.assertEqual(len(commands),1)
            cmd=commands[0]
            self.assertEqual(cmd[0],'nextflow')
            self.assertIn('-resume',cmd)
            self.assertEqual(cmd[cmd.index('-w')+1],str(raw/'nextflow_work'))
            self.assertEqual(cmd[cmd.index('--outdir')+1],str(raw/'results'))
            self.assertIn('--skip_ancestry',cmd)
            self.assertEqual(len((raw/'audit/output_sha256.txt').read_text().splitlines()),2)

    def test_no_nextflow_when_reconcile_already_valid(self):
        with tempfile.TemporaryDirectory() as tmp:
            raw=Path(tmp);(raw/'prep').mkdir();(raw/'scorefiles').mkdir();(raw/'audit').mkdir()
            (raw/'prep/score_ready.vcf.gz').write_bytes(b'a')
            (raw/'scorefiles/PGS004879_hmPOS_GRCh38.txt.gz').write_bytes(b'b')
            with patch.object(rr,'verify_existing_roundtrip',return_value={'OUT':str(raw)}), \
                 patch.object(rr,'run_logged_env'), \
                 patch.object(rr.v1,'log_run',side_effect=AssertionError('Nextflow should not run')), \
                 patch.object(rr.v1,'raw_integrity',return_value={'status':'PASS'}):
                rr.resume_raw_after_roundtrip(raw,'s1','/unused',{})

class V21PilotTest(unittest.TestCase):
    def test_frozen_pilot_six_combination_match(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'integrated_runs'/'pilot'
            for sample,models in run.v2.PILOT_TARGET.items():
                dest=root/sample
                (dest/'orchestration').mkdir(parents=True)
                (dest/'method2').mkdir()
                (dest/'orchestration/v21_status.json').write_text(json.dumps({'stage':'PASS'}))
                (dest/'orchestration/v21_inputs.lock.json').write_text('{}')
                fields=['sample','model','observed_loci','missing_loci','total_sum','eas_percentile']
                counts={'26T00970':{'A':(609032,1273),'B':(602628,810)},
                        '26T00969':{'A':(609162,1143),'B':(602768,670)},
                        '26T00968':{'A':(609272,1033),'B':(602869,569)}}
                with (dest/'method2/scores.tsv').open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t')
                    w.writeheader()
                    for model,(score,pct) in models.items():
                        obs,miss=counts[sample][model]
                        w.writerow({'sample':sample,'model':model,'observed_loci':obs,'missing_loci':miss,
                                    'total_sum':score,'eas_percentile':pct})
                (dest/'orchestration/v21_output_sha256.json').write_text(json.dumps(run.output_fingerprints(dest)))
            result=run.verify_pilot(tmp,'pilot')
            self.assertEqual(len(result),6)
            self.assertTrue(all(row[-1]=='PASS' for row in result))
            first=root/'26T00970/method2/scores.tsv'
            first.write_text('malformed')
            with self.assertRaisesRegex(RuntimeError,'invalid or changed'):
                run.verify_pilot(tmp,'pilot')
