#!/usr/bin/env python3
"""Benchmark the phenotype engine on phenopacket cases (plan §11, §12.1; seed of phase 3).

Each case's non-excluded HPO terms are scored; the known disease and gene are
located in the disease and gene rankings with mid ranks for ties, and
top-1 / top-10 / MRR are printed.

    .venv/bin/python scripts/benchmark_pheno.py --snapshot SNAPSHOT --cases phenopacket_dir \
        [--case-list ids.txt] [--gene-map genes_to_disease.txt]

``--gene-map`` replaces the snapshot's disease → gene links with HPO
genes_to_disease.txt (MENDELIAN / UNKNOWN rows), which is how the pilot in
§12.1 was run.  A JSON file of ``{"case": ID}`` records or a text file with one
ID per line restricts the cases.  ``--rows OUT.json`` saves per-case ranks.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import statistics
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "backend"))

from app.services.pheno_engine import engine as E  # noqa: E402


def load_cases(folder: Path) -> list[dict]:
    cases = []
    for path in sorted(glob.glob(os.path.join(folder, "**", "*.json"), recursive=True)):
        try:
            data = json.load(open(path, encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if not isinstance(data, dict) or "phenotypicFeatures" not in data:
            continue
        obs = [f.get("type", {}).get("id", "") for f in data.get("phenotypicFeatures") or [] if not f.get("excluded")]
        dis = [d.get("term", {}).get("id", "") for d in data.get("diseases") or []]
        genes = set()
        for it in data.get("interpretations") or []:
            diag = it.get("diagnosis") or {}
            for gi in diag.get("genomicInterpretations") or []:
                vd = (gi.get("variantInterpretation") or {}).get("variationDescriptor") or {}
                for sym in ((vd.get("geneContext") or {}).get("symbol"), (gi.get("gene") or {}).get("symbol")):
                    if sym:
                        genes.add(sym)
            if (diag.get("disease") or {}).get("id"):
                dis.append(diag["disease"]["id"])
        cases.append({"id": data.get("id"), "obs": obs, "dis": sorted(set(dis)), "genes": sorted(genes)})
    return cases


def read_ids(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    try:
        data = json.loads(text)
    except ValueError:
        return [line.strip() for line in text.splitlines() if line.strip()]
    return [r["case"] if isinstance(r, dict) else str(r) for r in data]


def g2d_links(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return [{"disease_id": r["disease_id"], "gene_symbol": r["gene_symbol"], "use": 1, "source": "g2d"}
                for r in rows if r.get("association_type") in ("MENDELIAN", "UNKNOWN")]


def metrics(ranks: list[float]) -> str:
    n = len(ranks)
    if not n:
        return "n=0"
    return (f"top1={sum(r <= 1 for r in ranks) / n:.3f}  top10={sum(r <= 10 for r in ranks) / n:.3f}  "
            f"MRR={sum(1 / r for r in ranks) / n:.3f}  median={statistics.median(ranks):.1f}  n={n}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--snapshot", type=Path, required=True)
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--case-list", type=Path, default=None)
    ap.add_argument("--gene-map", type=Path, default=None)
    ap.add_argument("--rows", type=Path, default=None)
    args = ap.parse_args(argv)
    started = time.time()
    links = g2d_links(args.gene_map) if args.gene_map else None
    eng = E.Engine(args.snapshot, links=links)
    cases = load_cases(args.cases)
    if args.case_list:
        by_id = {c["id"]: c for c in cases}
        cases = [by_id[i] for i in read_ids(args.case_list) if i in by_id]
    dranks, granks, rows, skipped = [], [], [], 0
    for c in cases:
        target = next((d for d in c["dis"] if d in eng.did), None)
        res = eng.score(c["obs"])
        if res.empty or target is None:
            skipped += 1
            continue
        d = float(res.d_rank_mid[eng.did[target]])
        gene = next((g for g in c["genes"] if g in eng.gid), None)
        g = float(res.g_rank_mid[eng.gid[gene]]) if gene else None
        dranks.append(d)
        if g is not None:
            granks.append(g)
        rows.append({"case": c["id"], "disease": target, "disease_rank": d, "gene": gene, "gene_rank": g,
                     "terms": len(res.q)})
    print(f"知識庫 {eng.snapshot}：{eng.kb.N:,} 個疾病、{len(eng.genes):,} 個基因；病例 {len(rows)}（略過 {skipped}）")
    print(f"疾病  {metrics(dranks)}")
    print(f"基因  {metrics(granks)}")
    print(f"耗時 {time.time() - started:.0f} 秒")
    if args.rows:
        args.rows.write_text(json.dumps(rows, ensure_ascii=False, indent=1), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
