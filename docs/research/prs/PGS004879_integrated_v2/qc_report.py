#!/usr/bin/env python3
"""Read official pgsc_calc outputs + validated raw reconciliation; produce QC TSV/JSON/MD.
This file never alters genotype/reference/Nextflow outputs.
"""
import argparse
import csv
import gzip
import json
import math
import os
from bisect import bisect_right
from pathlib import Path

BOUNDARIES = [20.0, 40.0, 60.0, 80.0, 90.0, 95.0]
GROUPS = ["<20", "20-40", "40-60", "60-80", "80-90", "90-95", ">=95"]


def group(p):
    for i, x in enumerate(BOUNDARIES):
        if p < x:return GROUPS[i]
    return GROUPS[-1]


def open_rows(path, delimiter):
    fn=gzip.open if str(path).endswith(".gz") else open
    with fn(path,"rt",newline="") as f:
        yield from csv.DictReader(f,delimiter=delimiter)


def atom_write(path, text):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    temp=p.with_name(p.name+".tmp")
    temp.write_text(text)
    os.replace(temp,p)


def calc(sample,rawdir,ancdir,score_sha256,ref_panel,min_match,high_weight,boundary_pp):
    rawdir=Path(rawdir);ancdir=Path(ancdir)
    base=ancdir/"results"/f"case{sample}"
    scoredir=base/"score"
    reconcile=json.loads((rawdir/"audit/final_reconciliation.json").read_text())
    if reconcile.get("iid") != sample or reconcile.get("status") != "PASS":
        raise RuntimeError("raw reconciliation did not PASS for requested sample")
    raw_summary=rawdir/"audit/upstream_reconstruction_summary.json"
    up=json.loads(raw_summary.read_text()) if raw_summary.is_file() else {}
    poprows=list(open_rows(scoredir/f"case{sample}_popsimilarity.txt.gz","\t"))
    pgsrows=list(open_rows(scoredir/f"case{sample}_pgs.txt.gz","\t"))
    pop=next(r for r in poprows if r.get("sampleset")==f"case{sample}" and r.get("IID")==sample)
    target=next(r for r in pgsrows if r.get("sampleset")==f"case{sample}" and r.get("IID")==sample)
    references={r["IID"]:float(r["SUM"]) for r in pgsrows if r.get("sampleset")=="reference"}
    EAS={r["IID"] for r in poprows if r.get("REFERENCE")=="True" and r.get("SuperPop")=="EAS"}
    eas_scores=sorted(references[i] for i in EAS if i in references)
    matchrows=list(open_rows(base/"match"/f"case{sample}_summary.csv",","))
    counts={r["match_status"].strip().lower():int(r["count"]) for r in matchrows}
    matched,excluded,unmatched=[counts.get(k,0) for k in ("matched","excluded","unmatched")]
    total=matched+excluded+unmatched
    if not total:raise RuntimeError("Empty match QC summary")
    normsum=float(target["SUM"])
    pct=float(target["percentile_MostSimilarPop"])
    pct_empirical=100*bisect_right(eas_scores,normsum)/len(eas_scores) if eas_scores else float("nan")
    # Sorted proportion is a cross-check, not a replacement for pgsc_calc's own percentile.
    distance=min(abs(pct-x) for x in BOUNDARIES)
    marker=ancdir/"work/ancestry/filter/GRCh38_HGDP+1kGP_ALL_thinned.prune.in.gz"
    if marker.is_file():
        with gzip.open(marker,"rt") as f:markers=sum(1 for _ in f)
    else:
        candidates=list((ancdir/"work").rglob("*thinned.prune.in.gz"))
        markers=0
        if candidates:
            with gzip.open(candidates[0],"rt") as f:markers=sum(1 for _ in f)
    missing_high=[]
    path=base/"match"/f"case{sample}_log.csv.gz"
    for r in open_rows(path,","):
        if r.get("match_status") == "matched":continue
        try:w=float(r["effect_weight"])
        except (TypeError,ValueError):continue
        if abs(w) >= high_weight:
            missing_high.append({k:r.get(k,"") for k in ("chr_name","chr_position","effect_allele","other_allele","effect_weight","match_status","exclude","match_type","row_nr")})
    missing_high.sort(key=lambda r: abs(float(r["effect_weight"])),reverse=True)
    flags=[];errors=[]
    if total!=610677:errors.append(f"Expected exactly 610677 PGS rows, observed {total}")
    if matched/total<min_match:errors.append(f"Variant match fraction {matched/total:.4%} below {min_match:.4%}")
    if len(references)!=3330:flags.append(f"Reference count differs from pilot 3330: {len(references)}")
    if len(eas_scores)!=735:flags.append(f"EAS reference count differs from pilot 735: {len(eas_scores)}")
    if pop.get("MostSimilarPop")!="EAS":flags.append(f"Ancestry is {pop.get('MostSimilarPop')}, not EAS")
    if pop.get("MostSimilarPop_LowConfidence")!="False":flags.append("Ancestry assignment marked low confidence")
    if not math.isfinite(pct) or not (0<=pct<=100):errors.append("Invalid normalized percentile")
    if distance < boundary_pp:flags.append(f"Within {boundary_pp:g} percentile points of INTERVENE stratum boundary ({distance:.4f})")
    if markers==0:flags.append("Missing/empty LD-pruned PCA marker count")
    if missing_high:flags.append(f"{len(missing_high)} nonmatched variants with abs(weight)>={high_weight}; inspect missing_high_weight.tsv")
    if len(eas_scores) and abs(pct-pct_empirical)>0.5:flags.append(f"EAS empirical ECDF {pct_empirical:.4f} differs from reported percentile {pct:.4f}; check normalization definition")
    if reconcile.get("observed_sum") is None:errors.append("Raw SUM reconciliation absent")
    if len(pgsrows)==0:errors.append("No scoring rows")
    status="FAIL" if errors else "PASS_WITH_REVIEW_FLAGS" if flags else "PASS"
    result={
      "sample":sample,"qc_status":status,"score_id":target.get("PGS"),
      "raw_sum":float(reconcile["observed_sum"]),"raw_denom":int(reconcile["observed_denom"]),
      "raw_matched":int(up.get("matched",reconcile.get("matched",0))),
      "normalized_sum":normsum,"percentile":pct,"percentile_empirical_check":pct_empirical,
      "z_most_similar_pop":float(target["Z_MostSimilarPop"]),
      "z_norm1":float(target["Z_norm1"]),"z_norm2":float(target["Z_norm2"]),
      "ancestry":pop.get("MostSimilarPop"),"rf_p_eas":float(pop.get("RF_P_EAS","nan")),
      "ancestry_low_confidence":pop.get("MostSimilarPop_LowConfidence"),
      "reference_n":len(references),"eas_reference_n":len(eas_scores),
      "matched":matched,"excluded":excluded,"unmatched":unmatched,
      "match_fraction":matched/total,"pca_ld_pruned_markers":markers,
      "intervene_group":group(pct),"nearest_group_boundary_pp":distance,
      "boundary_review":distance<boundary_pp,"high_weight_unmatched_n":len(missing_high),
      "high_weight_abs_cutoff":high_weight,"min_match_fraction":min_match,
      "scoring_file_sha256":score_sha256,"reference_panel":ref_panel,
      "raw_dir":str(rawdir.resolve()),"ancestry_dir":str(ancdir.resolve()),
      "flags":flags,"errors":errors,
      "research_only":True,"absolute_risk_calibrated_to_taiwan":False
    }
    return result,missing_high


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("sample","raw","ancestry","out","score-sha256","ref-panel"):
        p.add_argument("--"+key,required=True)
    p.add_argument("--min-match",type=float,default=0.995)
    p.add_argument("--high-weight",type=float,default=0.001)
    p.add_argument("--boundary-pp",type=float,default=2.0)
    a=p.parse_args()
    result,missing=calc(a.sample,a.raw,a.ancestry,a.score_sha256,a.ref_panel,a.min_match,a.high_weight,a.boundary_pp)
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    atom_write(out/"qc.json",json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    keys=list(result.keys())
    with (out/"qc.tsv").open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=keys,delimiter="\t")
        wr.writeheader(); wr.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in result.items()})
    with (out/"missing_high_weight.tsv").open("w",newline="") as f:
        cols=["chr_name","chr_position","effect_allele","other_allele","effect_weight","match_status","exclude","match_type","row_nr"]
        wr=csv.DictWriter(f,fieldnames=cols,delimiter="\t");wr.writeheader();wr.writerows(missing)
    lines=[f"# PGS004879 — {a.sample} QC Report", "",f"**Result: {result['qc_status']}** | Research use only; not a clinically calibrated CHD risk estimate.","",
           "| Metric | Value |", "|---|---:|",
           f"| Raw SUM (own matched loci) | {result['raw_sum']:.6f} |",
           f"| Raw matched variants | {result['raw_matched']:,} |",
           f"| Normalized-run PGS SUM (target-specific matched loci) | {result['normalized_sum']:.6f} |",
           f"| EAS percentile (pgsc_calc) | {result['percentile']:.4f} |",
           f"| Z_MostSimilarPop | {result['z_most_similar_pop']:.6f} |",
           f"| INTERVENE stratum (exploratory) | {result['intervene_group']} |",
           f"| Nearest boundary (percentage points) | {result['nearest_group_boundary_pp']:.4f} |",
           f"| Ancestry / RF_P_EAS | {result['ancestry']} / {result['rf_p_eas']} |",
           f"| Reference / EAS reference | {result['reference_n']} / {result['eas_reference_n']} |",
           f"| Matched / excluded / unmatched | {result['matched']:,} / {result['excluded']:,} / {result['unmatched']:,} |",
           f"| Match fraction | {result['match_fraction']:.5%} |",
           f"| PCA LD-pruned marker count | {result['pca_ld_pruned_markers']:,} |",
           f"| Nonmatched high-weight SNPs (abs weight >= {a.high_weight}) | {len(missing)} |",
           "", "## Review / Fail", ""]
    lines.extend([f"- {x}" for x in result["errors"]+result["flags"]] or ["- None"])
    lines.extend(["", "## High-weight nonmatched SNPs (top 20)","","| GRCh38 chr:pos | Effect allele | Weight | Status |","|---|---|---:|---|"])
    for r in missing[:20]:
        lines.append(f"| {r['chr_name']}:{r['chr_position']} | {r['effect_allele']} | {r['effect_weight']} | {r['match_status']} |")
    lines += ["", "## Important interpretation limits", "",
              "- Raw SUM and pgsc_calc SUM use different matched-variant subsets; they are **not** expected to be identical.",
              "- In pgsc_calc v2.3.0, reference scoring can differ by target because matching excludes different loci.",
              "- RF_P_EAS=1.0 is a classifier probability, not a biological ancestry proportion.",
              "- An INTERVENE stratum derived from this pilot reference is not a validated Taiwan CHD absolute-risk category.",
              "- Variant matching, genotype quality, strand/palindromic loci, reference mean imputation and scores near group boundaries require review.",
              "", f"**Reference panel:** `{result['reference_panel']}`", f"**Corrected scoring-file SHA256:** `{result['scoring_file_sha256']}`", ""]
    atom_write(out/"report.md","\n".join(lines))
    print(json.dumps({k:result[k] for k in ("sample","qc_status","raw_sum","normalized_sum","percentile","intervene_group","matched","pca_ld_pruned_markers","high_weight_unmatched_n","flags","errors")},ensure_ascii=False,indent=2))
    if result["qc_status"]=="FAIL":raise SystemExit(1)

if __name__=="__main__":main()
