#!/usr/bin/env python3
"""Recover exact source gVCF paths from validated raw run source_sha256.txt files.
Does not assume NovaSeq batch paths or scan the filesystem for patient IDs.
"""
import argparse
import csv
from pathlib import Path

SAMPLES=("26T00970","26T00969","26T00968")

def source_for(raw):
    p=raw/"audit/source_sha256.txt"
    if not p.is_file():return None
    for line in p.read_text().splitlines():
        parts=line.split(maxsplit=1)
        if len(parts)!=2:continue
        name=parts[1].strip().removeprefix("*")
        if name.endswith(".gvcf.gz") and len(parts[0])==64:
            return name
    return None

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project",required=True)
    ap.add_argument("--out",required=True)
    a=ap.parse_args()
    project=Path(a.project)
    rows=[]
    for sid in SAMPLES:
        raw=project/"results"/f"{sid}_PGS004879_FINAL_gVCF_to_raw_score"
        g=source_for(raw)
        if not g:print(f"WARNING: gVCF provenance unavailable for {sid}: {raw}")
        elif not Path(g).is_file():print(f"WARNING: provenance path not found on current server for {sid}: {g}")
        rows.append({"sample":sid,"gvcf":g or "","raw_source":str(raw)})
    o=Path(a.out);o.parent.mkdir(parents=True,exist_ok=True)
    with o.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["sample","gvcf","raw_source"],delimiter="\t")
        w.writeheader();w.writerows(rows)
    print("MANIFEST:",o)
    for row in rows:print(row)

if __name__=="__main__":main()
