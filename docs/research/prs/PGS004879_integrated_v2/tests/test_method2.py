import csv
import math
import sys
import tempfile
import unittest
from decimal import Decimal
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import method2_score as m2


class CanonicalMethod2Tests(unittest.TestCase):
    def setUp(self):
        self.d=tempfile.TemporaryDirectory()
        self.root=Path(self.d.name)
        self.audit=self.root/'all_score_locus_decisions.tsv'
        rows=[
            {'CHROM':'chr1','POS':'100','EFFECT_ALLELE':'A','OTHER_ALLELE':'G','WEIGHT':'0.5',
             'STATUS':'CANDIDATE_DIRECT','RECALC_EFFECT_DOSAGE':'1','RECALC_CONTRIBUTION':'0.5','VCF_REF':'A','VCF_ALT':'G','VCF_GT':'0/1','ALT_DOSAGE':'1'},
            {'CHROM':'chr2','POS':'200','EFFECT_ALLELE':'A','OTHER_ALLELE':'G','WEIGHT':'-0.2',
             'STATUS':'LOW_GQ','RECALC_EFFECT_DOSAGE':'','RECALC_CONTRIBUTION':'','VCF_REF':'C','VCF_ALT':'T','VCF_GT':'0/1','ALT_DOSAGE':''},
        ]
        with self.audit.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
        self.master={
            ('1',100):{'row':1,'id':'1:100:A:G','ea':'A','oa':'G','ref':'A','alt':'G',
                       'weight':Decimal('0.5'),'p':Decimal('0.25'),'model_b':True},
            ('2',200):{'row':2,'id':'2:200:C:T','ea':'T','oa':'C','ref':'C','alt':'T',
                       'weight':Decimal('-0.2'),'p':Decimal('0.25'),'model_b':False},
        }
        self.eas={f'R{x:03d}' for x in range(735)}
        self.reference={'A':{iid:i*.001 for i,iid in enumerate(sorted(self.eas))},
                        'B':{iid:i*.001 for i,iid in enumerate(sorted(self.eas))}}

    def tearDown(self):
        self.d.cleanup()

    def test_scoring_observed_and_imputed(self):
        result=m2.score_audit(self.audit,self.master,self.reference,self.eas,strict=False)
        a,b=result['scores']
        self.assertAlmostEqual(a['observed_sum'],.5)
        self.assertAlmostEqual(a['imputed_sum'],-.1)
        self.assertAlmostEqual(a['total_sum'],.4)
        self.assertEqual((a['observed_loci'],a['missing_loci']), (1,1))
        self.assertEqual((b['observed_loci'],b['missing_loci']), (1,0))
        self.assertEqual(result['metadata_orientation_counts']['REVERSE_COMPLEMENT_METADATA'],1)
        self.assertEqual(len(result['missing']),1)
        # Output artifacts must be materialized and replayable.
        model_lock=self.root/'model_lock.json';model_lock.write_text('{}')
        summary=m2.write_outputs(self.root/'out','SAMPLE',self.audit,result,model_lock)
        self.assertEqual(summary['scores'][0]['total_sum'],.4)
        self.assertTrue((self.root/'out/missing_loci.tsv.gz').is_file())
        self.assertTrue((self.root/'out/report.md').is_file())

    def test_reject_unexpected_nonmissing_failed_qc(self):
        s=self.audit.read_text(); self.audit.write_text(s.replace('LOW_GQ\t\t','LOW_GQ\t1\t'))
        with self.assertRaisesRegex(RuntimeError,'Failed-QC locus has dosage'):
            m2.score_audit(self.audit,self.master,self.reference,self.eas,strict=False)

    def test_reject_allele_orientation_mismatch(self):
        self.master[('2',200)]['ea']='A'
        with self.assertRaisesRegex(RuntimeError,'allele orientation mismatch'):
            m2.score_audit(self.audit,self.master,self.reference,self.eas,strict=False)

    def test_reject_missing_canonical_row(self):
        self.master[('3',300)]={'row':3,'id':'3:300:A:G','ea':'A','oa':'G','ref':'A','alt':'G','weight':Decimal('0.1'),'p':Decimal('.2'),'model_b':True}
        with self.assertRaisesRegex(RuntimeError,'coverage mismatch'):
            m2.score_audit(self.audit,self.master,self.reference,self.eas,strict=False)

    def test_original_raw_vs_canonical_exclusion_ledger(self):
        # A raw v1 genotype can have passed QC while a reference-incompatible
        # position was omitted from fixed Model A. Keep it in the ledger.
        with self.audit.open('a') as f:
            f.write('chr3\t300\tG\tA\t-0.1\tCANDIDATE_DIRECT\t2\t-0.2\tG\tA\t1/1\t2\n')
        result=m2.score_audit(self.audit,self.master,self.reference,self.eas,strict=False)
        rec=result['raw_reconciliation']
        self.assertEqual(rec['outside_canonical_observed_n'],1)
        self.assertEqual(rec['outside_canonical_observed_exact'],'-0.2')
        self.assertEqual(rec['reconstructed_raw_sum'],'0.3')

    def test_unexpected_content_after_model_lock(self):
        model=self.root/'model';model.mkdir()
        for name in m2.MODEL_FILES:
            path=model/name; path.parent.mkdir(parents=True,exist_ok=True);path.write_text('fixture\n')
        source=self.root/'source.gz'
        import gzip,json
        with gzip.open(source,'wb') as f:f.write(b'fixture')
        lock={'schema':m2.SCHEMA,'model_dir':str(model.resolve()),'source_path':str(source.resolve()),
              'sha256':{name:m2.sha256(model/name) for name in m2.MODEL_FILES}}
        lockfile=self.root/'lock.json';lockfile.write_text(json.dumps(lock))
        (model/m2.MODEL_FILES[0]).write_text('tampered')
        # Genuine source identity is checked first and fails closed.
        with self.assertRaisesRegex(RuntimeError,'Official original scoring data checksum changed'):
            m2.validate_locked_model(model,lockfile,source)

if __name__=='__main__': unittest.main()
