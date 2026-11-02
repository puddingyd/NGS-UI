#!/usr/bin/env python3
"""Strict three-case regression against recorded pilot results; no server access required."""
import argparse
import csv
import json
from pathlib import Path


def check(run_dir,expected):
    cfg=json.loads(Path(expected).read_text())
    tolerances=cfg["tolerances"]
    rows=[]
    for sample, gold in cfg["samples"].items():
        output=Path(run_dir)/sample/"qc/qc.json"
        if not output.is_file():
            rows.append((sample,"file","MISSING",str(output),"FAIL"))
            continue
        got=json.loads(output.read_text())
        for key,value in gold.items():
            observed=got.get(key)
            if observed is None:
                rows.append((sample,key,"MISSING",value,"FAIL"));continue
            if key=="raw_sum":tol=tolerances["raw_sum_abs"]
            elif key=="normalized_sum":tol=tolerances["normalized_sum_abs"]
            elif key=="percentile":tol=tolerances["percentile_pp_abs"]
            elif key.startswith("z_"):tol=tolerances["z_abs"]
            else:tol=0
            passed=(abs(float(observed)-float(value))<=tol)
            rows.append((sample,key,observed,value,"PASS" if passed else "FAIL"))
        rows.append((sample,"qc_status",got.get("qc_status"),"PASS or PASS_WITH_REVIEW_FLAGS",
                     "FAIL" if got.get("qc_status")=="FAIL" else "PASS"))
    return rows


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-dir",required=True);ap.add_argument("--expected",required=True)
    a=ap.parse_args()
    rows=check(a.run_dir,a.expected)
    dest=Path(a.run_dir)/"regression_comparison.tsv";dest.parent.mkdir(parents=True,exist_ok=True)
    with dest.open("w",newline="") as f:
        wr=csv.writer(f,delimiter="\t");wr.writerow(["sample","metric","observed","expected","result"]);wr.writerows(rows)
    failed=[r for r in rows if r[-1]=="FAIL"]
    for r in rows:print("\t".join(map(str,r)))
    print(f"REGRESSION: {'FAIL' if failed else 'PASS'}; {len(failed)} mismatches; {dest}")
    if failed:raise SystemExit(1)

if __name__=="__main__":main()
