"""Reproduce all six EstBB endpoints before releasing Taiwan GBD2019 parameters.

Source CSVs are immutable snapshots. Published XLSX supplies expected values,
never coefficients. Run as python -m research.prs.absolute_risk.build_multi.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from . import build_taiwan as chd
from .engine import GROUPS, integrate, validate
from .registry import DISEASES

PINS = {**chd.PINS,
    'HR_FullSample_EstBB211122.csv':'5677324b18e97b3fae26f6b68d029f7222be5e2e338a6e997a0373cdb5227389',
    'HR_AgeStratified_EstBB_211122.csv':'af990427888c354b1d216334671fc9573bd907a27ea7dafd50ff22440dd754d2',
    'BreastCancerProstateCancer_Incidence.csv':'ed20c4c45881f9c8d0faa0524a37e3b9446f1fb06e4889a00eda13f9aa3c5555',
    'BreastCancerProstateCancer_Prevalence.csv':'440a7b799d79c0a0f82d250285673e21696d1fbc1888bb853adda8838f070056',
    'BreastCancerProstateCancer_Mortality.csv':'7f9d5c402c751d3a55dc81619eab5fe93494dcbb0bd2c1dda2aaf6109c6faaaa',
    '41467_2024_48938_MOESM7_ESM.xlsx':'659e782a8f8a9beb5168f0a74cc89100247dc74b33c3ac46c07474748340d5c8',
    'taiwan_raw_gbd_audit.json':'6204fd3d6092e36708ebec8aa31ce3c7a4ff207c6e0203d525dff32898f90752',
}
ASSOCIATIONS = {
    'age_sex': ('HR_AgeandSexStratified_EstBB_211122.csv', 'Phenotype PRS Sex MinAge MaxAge MedianAAO Group Controls Cases Beta SE Pval HR Cipos Cineg'),
    'age': ('HR_AgeStratified_EstBB_211122.csv', 'Phenotype PRS MinAge MaxAge MedianAAO Group Controls Cases Beta SE Pval HR Cipos Cineg'),
    'sex': ('HR_FemaleSample_EstBB.csv', 'Phenotype PRS Group Controls Cases Beta SE Pval HR Cipos Cineg'),
    'full': ('HR_FullSample_EstBB211122.csv', 'Phenotype PRS Group Controls Cases Beta SE Pval HR Cipos Cineg'),
}
# The public EstBB cancer CSV has the opposite score-group orientation from
# Supplementary Data 4 and the PGS Catalog EB performance (HR/SD 1.15).
# This label mapping is an inferred source reconciliation, independently
# checked against ALL 224 published values. It never reciprocates/fits HRs.
CANCER_SOURCE_GROUPS = dict(zip(GROUPS,
    ('> 80%', '60-80%', '40-60%', '20-40%', '10-20%', '5-10%', '< 5%')))


def checked(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != PINS[path.name]:
        raise ValueError(f'Source changed: {path}')
    return path


def associations(path, columns):
    with checked(path).open(newline='') as stream:
        rows = csv.reader(stream)
        fields = columns.split()
        if next(rows) != ['']+[f'V{i}' for i in range(1,len(fields)+1)]:
            raise ValueError('Unexpected author association schema')
        result = []
        for row in rows:
            if len(row) != len(fields)+1:
                raise ValueError('Malformed association row')
            result.append(dict(zip(fields,row[1:])))
        return result


def bins(disease_id, sources, hr, extra, country, minimum):
    disease = DISEASES[disease_id]
    coefficients = {}
    for mode in set(disease['hr'].values()):
        filename, fields = ASSOCIATIONS[mode]
        path = hr/filename if (hr/filename).is_file() else extra/filename
        coefficients[mode] = [r for r in associations(path,fields)
                              if r['Phenotype']==disease['phenotype'] and r['PRS']==disease['prs']]
    values = {}
    for measure in ('Incidence','Prevalence','Mortality'):
        path = (extra/f'BreastCancerProstateCancer_{measure}.csv' if disease_id=='BREAST'
                else sources/'gbd2019'/f'Sex_Stratified_{measure}.csv')
        for r in csv.DictReader(checked(path).open()):
            if r['location']!=country or r['metric']!='Rate' or r['cause'] not in (disease['cause'],'All causes'):
                continue
            key = (measure,r['sex'].lower(),r['age'],r['cause'])
            lo,val,hi = (float(r[k]) for k in ('lower','val','upper'))
            if key in values or r['year']!='2019' or not 0<=lo<=val<=hi<math.inf:
                raise ValueError('Invalid GBD row')
            values[key]=val/100000
    result=[]
    for sex in disease['sexes']:
        mode=disease['hr'][sex]
        rows=coefficients[mode]
        for start in range(minimum,80,5):
            label=f'{start or 1}-{start+4} years'
            hrs={'40-60':1.0}
            for group in GROUPS:
                if group=='40-60':continue
                source_group = CANCER_SOURCE_GROUPS[group] if disease_id == 'ALL_CANCERS' else chd.LABELS[group]
                selected=[{**r, 'Group': chd.LABELS[group]} for r in rows if r['Group']==source_group
                          and (mode!='age_sex' or r['Sex']==sex)]
                for r in selected:
                    if any(not math.isfinite(float(r[k])) or float(r[k])<=0 for k in ('HR','SE')):
                        raise ValueError('Invalid HR or SE')
                if mode in ('age','age_sex'):
                    # Same weighted log-HR regression and boundary clamp as CHD.
                    hrs[group]=chd.male_hr([{**r,'Sex':'male'} for r in selected],group,start+2.5)
                else:
                    if len(selected)!=1:raise ValueError(f'{disease_id}: expected unique HR')
                    hrs[group]=float(selected[0]['HR'])
            result.append(dict(sex=sex,age_start=start,age_end=start+5,hr=hrs,
                incidence=values['Incidence',sex,label,disease['cause']],
                prevalence=values['Prevalence',sex,label,disease['cause']],
                all_cause_mortality=values['Mortality',sex,label,'All causes'],
                disease_mortality=values['Mortality',sex,label,disease['cause']]))
    return result


def build(sources,hr,extra):
    import openpyxl  # Build-time only; not required by the clinical UI.
    publication=checked(sources/'publication/41467_2024_48938_MOESM7_ESM.xlsx')
    workbook=openpyxl.load_workbook(publication,read_only=True,data_only=True)
    expected=[r for sheet in workbook for r in sheet.values if len(r)>7 and r[1]=='Estonia']
    workbook.close()
    models={}; proof={}
    raw_sources=json.loads(checked(extra/'taiwan_raw_gbd_audit.json').read_text())
    raw={}
    for source in raw_sources:
        for row in source['rows']:
            key=(row['cause'],row['measure'],row['sex'].lower(),row['age'])
            if key in raw:raise ValueError('Duplicate original Taiwan baseline')
            raw[key]=float(row['val'])/100000
    for disease_id,disease in DISEASES.items():
        common=dict(schema='PRS_DISEASE_BUNDLE_V2',disease_id=disease_id,endpoint=disease['endpoint'],
            status='unvalidated_external',model_id=f'intervene_estbb_tw2019_{disease_id.lower()}_v1',
            method='competing_exact_v1',population='Taiwan',baseline_year=2019,baseline_release='GBD2019',
            baseline_source=chd.SOURCE_URLS['gbd'],hr_source='INTERVENE EstBB original cohort',
            units='rates_per_person_year_prevalence_fraction',age_min=20,
            score_definition=dict(pgs_id=disease['pgs_id'],reference_id='unrelated3330_EAS735_v1',ancestry='EAS'),
            limitations=[disease['limitation'],'研究用途；歐洲相對風險外推至台灣，尚未完成台灣追蹤結果驗證。'],
            source_sha256=PINS,
            source_group_mapping=CANCER_SOURCE_GROUPS if disease_id=='ALL_CANCERS' else chd.LABELS,
            source_group_mapping_evidence=('INFERRED: reflection of public EstBB source labels; verified against all 224 published Estonia values; PGS Catalog PPM021187 reports EB HR/SD 1.15'
                if disease_id=='ALL_CANCERS' else 'Original published percentile labels'))
        reproduction={**common,'population':'Estonia','age_min':0,'method':'intervene_discrete_v1',
                      'bins':bins(disease_id,sources,hr,extra,'Estonia',0)}
        comparisons=[];seen=set()
        for row in expected:
            if row[0]!=disease['published']:continue
            sex={'Males':'male','Females':'female'}[row[2]]
            mode={'age_sex':'Age and Sex','age':'Age','sex':'Sex','full':'None'}[disease['hr'][sex]]
            if row[4]!=mode:raise ValueError('Published model selection differs')
            group=next(k for k,v in chd.LABELS.items() if v==row[5]); age=float(row[6])
            key=(sex,group,age)
            if key in seen:raise ValueError('Duplicate published row')
            seen.add(key)
            actual=100*integrate(reproduction,sex,chd.PERCENTILES[group],0,age)['risk']
            comparisons.append(dict(sex=sex,group=group,age=age,expected=float(row[7]),
                                    actual=actual,error_pp=abs(actual-float(row[7]))))
        required={(s,g,float(a)) for s in disease['sexes'] for g in GROUPS for a in range(5,81,5)}
        maximum=max((r['error_pp'] for r in comparisons),default=math.inf)
        if seen!=required or maximum>.0050001:
            raise ValueError(f'{disease_id}: reproduction failed ({len(seen)} rows; max={maximum} pp)')
        evidence=dict(status='PASS',rows=len(seen),max_error_pp=maximum,tolerance_pp=.0050001,
                      source='Supplementary Data 4 / Estonia',uncertainty_reproduced=False)
        models[disease_id]=validate({**common,'reproduction':evidence,
            'bins':bins(disease_id,sources,hr,extra,'Taiwan (Province of China)',20)})
        checked_rows=0; author_zero_fills=[]
        for row in models[disease_id]['bins']:
            age=f"{row['age_start']}-{row['age_end']-1} years"
            for field,measure,cause in [('incidence','Incidence',disease['cause']),
                ('prevalence','Prevalence',disease['cause']),('disease_mortality','Deaths',disease['cause']),
                ('all_cause_mortality','Deaths','All causes')]:
                key=(cause,measure,row['sex'],age)
                if key not in raw:
                    if not (disease_id=='AF' and field=='disease_mortality' and row['age_start']<30 and row[field]==0):
                        raise ValueError(f'Missing original baseline: {key}')
                    author_zero_fills.append(key)
                elif not math.isclose(row[field],raw[key],rel_tol=1e-12,abs_tol=1e-14):
                    raise ValueError(f'Original baseline mismatch: {key}')
                else:checked_rows+=1
        models[disease_id]['baseline_audit']=dict(original_rate_cells=checked_rows,
            author_missing_death_rate_zero_fills=author_zero_fills)
        proof[disease_id]={**evidence,'comparisons':comparisons}
        print(disease_id,len(seen),maximum,flush=True)
    return dict(schema='PRS_RISK_PACK_V2',models=models),dict(status='PASS',diseases=proof)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('sources','hr','extra','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();pack,proof=build(a.sources,a.hr,a.extra)
    a.output.mkdir(parents=True,exist_ok=False)
    for name,data in [('bundle.json',pack),('reproduction.json',proof)]:
        (a.output/name).write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
