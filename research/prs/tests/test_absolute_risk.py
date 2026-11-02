import copy
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from research.prs.absolute_risk.engine import GROUPS, PROPORTIONS, group, integrate, validate
from research.prs.absolute_risk.__main__ import evaluate

ROOT = Path(__file__).resolve().parents[1]


class RiskTests(unittest.TestCase):
    def setUp(self):
        def fixture(name):
            return json.loads((ROOT / 'examples' / (name + '.json')).read_text())
        self.b = fixture('synthetic_bundle')
        self.s = fixture('synthetic_summary')
        self.p = fixture('synthetic_person')

    def test_constant_hazard_analytical_solution(self):
        result = integrate(self.b, 'female', 50, 52.5, 62.5)
        self.assertAlmostEqual(result['risk'], .01/.03*(1-math.exp(-.03*10)), places=14)
        self.assertAlmostEqual(result['curve'][-1]['event_free_survival'], math.exp(-.3))

    def test_conditional_identity(self):
        a = integrate(self.b, 'male', 50, 0, 50)
        b = integrate(self.b, 'male', 50, 0, 60)
        conditional = (b['risk']-a['risk'])/a['curve'][-1]['event_free_survival']
        self.assertAlmostEqual(conditional, integrate(self.b, 'male', 50, 50, 60)['risk'])

    def test_original_interval_formula(self):
        self.b['method'] = 'intervene_discrete_v1'
        expected = (1-math.exp(-.05))*(1+math.exp(-.15))
        self.assertAlmostEqual(integrate(self.b, 'male', 50, 50, 60)['risk'], expected)
        with self.assertRaises(ValueError):
            integrate(self.b, 'male', 50, 51, 61)

    def test_zero_incidence_and_no_competing_mortality(self):
        for r in self.b['bins']:
            r['incidence'] = 0
        self.assertEqual(integrate(self.b, 'male', 0, 0, 80)['risk'], 0)
        for r in self.b['bins']:
            r['incidence'], r['all_cause_mortality'] = .01, 0
        self.assertAlmostEqual(integrate(self.b, 'male', 50, 0, 80)['risk'], 1-math.exp(-.8))

    def test_stratum_incidence_normalization(self):
        hrs = [.5,.7,1,1.2,1.4,1.6,2]
        for r in self.b['bins']:
            r['hr'] = dict(zip(GROUPS,hrs))
        weighted = sum(p*h for p,h in zip(PROPORTIONS,hrs))
        self.assertAlmostEqual(sum(p*.01*h/weighted for p,h in zip(PROPORTIONS,hrs)), .01)
        self.assertLess(integrate(self.b,'female',10,50,60)['risk'],integrate(self.b,'female',99,50,60)['risk'])

    def test_boundary_convention(self):
        self.assertEqual([group(x) for x in (0,20,40,60,80,90,95)],list(GROUPS))
        self.assertEqual(group(100),'>=95')
        for p in (-1,101,float('nan'),float('inf'),True):
            with self.assertRaises(ValueError): group(p)

    def test_bundle_rejections(self):
        edits = [('prevalence',1),('incidence',-1),('incidence',float('nan')),
                 ('chd_mortality',.5),('age_end',0),('sex','unknown')]
        for key,value in edits:
            with self.subTest(key=key,value=value):
                b=copy.deepcopy(self.b);b['bins'][0][key]=value
                with self.assertRaises(ValueError):validate(b)
        for mutation in ('gap','duplicate','missing_hr'):
            b=copy.deepcopy(self.b)
            if mutation=='gap':b['bins'].pop(2)
            elif mutation=='duplicate':b['bins'].append(b['bins'][0])
            else:del b['bins'][0]['hr']['>=95']
            with self.assertRaises(ValueError):validate(b)

    def test_person_and_score_gates(self):
        with self.assertRaises(ValueError):evaluate(self.b,self.s,self.p)
        for key,val in [('prior_chd',True),('prior_chd',None),('sample','wrong'),('sex','unknown'),('age',80),('assessment_date','invalid'),('synthetic',False)]:
            p={**self.p,key:val}
            with self.subTest(key=key),self.assertRaises(ValueError):evaluate(self.b,self.s,p,True)
        for key,val in [('model_lock_sha256','bad'),('ancestry_status',{}),('stage','FAILED')]:
            with self.subTest(key=key),self.assertRaises(ValueError):evaluate(self.b,{**self.s,key:val},self.p,True)

    def test_horizon_and_research_status(self):
        result=evaluate(self.b,self.s,{**self.p,'age':75},True)
        self.assertIsNone(result['risk_10_year'])
        self.assertGreater(result['risk_to_80']['risk'],0)
        self.assertEqual(result['clinical_validation'],'NOT_VALIDATED')
        self.assertIsNone(result['uncertainty_interval'])

    def test_cli_and_refusal_to_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            output=Path(d)/'result.json'
            cmd=[sys.executable,'-m','research.prs.absolute_risk','evaluate','--bundle',str(ROOT/'examples/synthetic_bundle.json'),'--summary',str(ROOT/'examples/synthetic_summary.json'),'--person',str(ROOT/'examples/synthetic_person.json'),'--allow-synthetic','--output',str(output)]
            subprocess.run(cmd,check=True,capture_output=True)
            before=output.read_bytes()
            self.assertEqual(len(json.loads(before)['input_sha256']),3)
            self.assertNotEqual(subprocess.run(cmd,capture_output=True).returncode,0)
            self.assertEqual(output.read_bytes(),before)

    def test_import_checksums(self):
        manifest=json.loads((ROOT/'IMPORT_PROVENANCE.json').read_text())
        for name,digest in manifest['files'].items():
            self.assertEqual(hashlib.sha256((ROOT/'scoring_v211'/name).read_bytes()).hexdigest(),digest,name)


if __name__ == '__main__':
    unittest.main()
