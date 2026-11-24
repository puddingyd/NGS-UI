#!/usr/bin/env python3
"""Literature validation of the phenotype engine (plan §11, §15 stage 6).

Full report (time-split on a frozen old knowledge base, LOPO and full on the
current one, comparators, statistics, rank-point calibration)::

    .venv/bin/python scripts/benchmark_pheno.py run \\
        --cases phenopacket-store/0.1.27 \\
        --old-snapshot hpo_kb_2023/2023-10-09_2023-10-09 --old-gene-map genes_to_disease_v2023-10-09.txt \\
        --current-snapshot hpo_kb/2026-09-01_2026-09-02 --current-gene-map genes_to_disease.txt \\
        --out benchmark_out [--lirical lirical_dir] [--exomiser exomiser_dir]

Case sets are the locked lists in ``research/pheno_benchmark/`` (``--sets``).
``select-test`` drew the locked test list once; it is kept for the record and
refuses to overwrite an existing file.

Legacy quick mode (one snapshot, one list)::

    .venv/bin/python scripts/benchmark_pheno.py --snapshot S --cases DIR [--case-list ids] [--gene-map g2d]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "backend"))

from app.services.pheno_engine import benchmark as B  # noqa: E402
from app.services.pheno_engine import engine as E  # noqa: E402

SETS = REPO / "research" / "pheno_benchmark"
SET_LABELS = {"dev": "開發資料 400 例（參數在此選定）", "fresh": "第二組 481 例（點數合併校準用）",
              "test": "鎖定測試組（未參與任何參數）"}


def g2d_links(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return [{"disease_id": r["disease_id"], "gene_symbol": r["gene_symbol"], "use": 1, "source": "g2d"}
                for r in rows if r.get("association_type") in ("MENDELIAN", "UNKNOWN")]


def load_engine(snapshot: Path, gene_map: Path | None):
    return E.Engine(snapshot, links=g2d_links(gene_map) if gene_map else None)


def progress(label):
    def fn(k, n):
        print(f"  {label} {k}/{n}", flush=True)
    return fn


def load_sets(folder: Path, by_id: dict) -> dict:
    out = {}
    for key in ("dev", "fresh", "test"):
        path = folder / f"{key}.txt"
        if path.is_file():
            ids = B.read_ids(path)
            out[key] = [by_id[i] for i in ids if i in by_id]
            missing = len(ids) - len(out[key])
            if missing:
                print(f"注意：{key} 有 {missing} 例在病例資料夾找不到（phenopacket-store 版本不同？）")
    return out


def cmd_run(args) -> int:
    t0 = time.time()
    cases = B.load_cases(args.cases)
    by_id = {c.id: c for c in cases}
    sets = load_sets(args.sets, by_id)
    print(f"病例 {len(cases):,}；" + "、".join(f"{k} {len(v)}" for k, v in sets.items()), flush=True)
    external = {}
    for tool in ("lirical", "exomiser"):
        folder = getattr(args, tool)
        if folder:
            external[tool] = B.external_ranks(folder)
            print(f"{tool}：讀到 {len(external[tool])} 例的結果")
    methods = ["engine", "nobaseline", "resnik", "platform"]
    sections, meta = [], {"cases_dir": str(args.cases), "calibration": {}}
    old = load_engine(args.old_snapshot, args.old_gene_map)
    meta["old_snapshot"] = old.snapshot
    groups_old = len({old.groups.get(d, d) for d in old.kb.diseases})
    comp = B.Comparators(old)
    ts_rows = {}
    for key, sel in sets.items():
        rows, cov = B.evaluate(old, sel, methods=methods, comparators=comp, external=external,
                               progress=progress(f"time-split {key}"))
        ts_rows[key] = rows
        sections.append({"key": f"timesplit_{key}", "title": f"Time-split（主要條件）：{SET_LABELS.get(key, key)}",
                         "note": f"知識庫凍結在 {old.snapshot}；只評估知識庫建立時尚未引用其文獻的病例（第 11.2 節）。",
                         "rows": rows, "coverage": cov, "methods": methods + list(external)})
    calib_rows = [r for k in ("dev", "fresh") for r in ts_rows.get(k, [])]
    if calib_rows:
        meta["calibration"]["基因名次（開發＋第二組合併）"] = {
            "note": f"{sum(r.get('engine|g') is not None for r in calib_rows)} 例，$f_0$＝區間寬度／{len(old.genes):,} 個基因；95% 區間為以基因群組 bootstrap 2,000 次。",
            "rows": B.calibrate([r.get("engine|g") for r in calib_rows], len(old.genes),
                                cohorts=[r["cohort"] for r in calib_rows])}
        meta["calibration"]["疾病名次（開發＋第二組合併，MONDO 同組取較佳名次）"] = {
            "note": f"{len(calib_rows)} 例，$f_0$＝區間寬度／{groups_old:,} 個疾病組（OMIM／Orphanet 同一疾病算一組）。"
                    "只用於同一基因內不同疾病的分配（第 4.8.7 節）。",
            "rows": B.calibrate([r.get("engine|dg") for r in calib_rows], groups_old,
                                cohorts=[r["cohort"] for r in calib_rows])}
    if args.current_snapshot:
        cur = load_engine(args.current_snapshot, args.current_gene_map)
        meta["current_snapshot"] = cur.snapshot
        everyone = [c for sel in sets.values() for c in sel]
        comp_cur = B.Comparators(cur)
        rows, cov = B.evaluate(cur, everyone, methods=methods, comparators=comp_cur, external=external,
                               progress=progress("full"))
        full_rows = rows
        sections.append({"key": "full", "title": "完整知識庫（樂觀上限，只當正確性檢查）",
                         "note": f"知識庫 {cur.snapshot}；病例的文獻多半已被註記引用（考古題考到原題）。",
                         "rows": rows, "coverage": cov, "methods": methods + list(external)})
        if not args.no_lopo:
            print("LOPO：計算移除後的基準分 …", flush=True)
            specs = B.lopo_specs(cur, everyone, B.disease_annotation_refs(args.current_snapshot))
            B.lopo_baselines(cur, specs, progress=progress("LOPO baseline"))
            cited_cases = [c for c in everyone if specs.get(c.id) and specs[c.id].removed]
            rows, cov = B.evaluate(cur, cited_cases, methods=["engine", "nobaseline", "resnik"], lopo=specs,
                                   comparators=comp_cur, progress=progress("LOPO"))
            removed = sorted(r["lopo_removed"] / r["lopo_total"] for r in rows if r.get("lopo_total"))
            ids = {r["case"] for r in rows}
            before = B.metrics([r.get("engine|g") for r in full_rows if r["case"] in ids])
            after = B.metrics([r.get("engine|g") for r in rows])
            note = (f"知識庫 {cur.snapshot}。{len(everyone)} 例中有 {len(cited_cases)} 例的文獻已被目前知識庫的目標疾病註記引用；"
                    f"把這些註記從目標疾病移除（中位數移除 {removed[len(removed) // 2]:.0%}；基準分隨之重算，IC 不重算），"
                    f"模擬「這個病只能從其他文獻得知」。同一批病例在完整知識庫的基因 MRR {before.get('mrr', 0):.3f}、"
                    f"top-10 {before.get('top10', 0):.3f}，移除後 {after.get('mrr', 0):.3f}、{after.get('top10', 0):.3f}。") if removed else ""
            sections.append({"key": "lopo", "title": "LOPO（悲觀壓力測試）", "note": note, "rows": rows,
                             "coverage": cov, "methods": ["engine", "nobaseline", "resnik"]})
    meta["intro"] = (f"病例資料 {args.cases}；time-split 知識庫 {meta['old_snapshot']}"
                     + (f"；目前知識庫 {meta.get('current_snapshot')}" if meta.get("current_snapshot") else "")
                     + f"。同分用平均名次。產生時間 {time.strftime('%Y-%m-%d %H:%M')}，耗時 {time.time() - t0:.0f} 秒。")
    B.write_report(args.out, sections, meta)
    print(f"完成：{args.out / 'report.md'}（{time.time() - t0:.0f} 秒）")
    for sec in sections:
        rows = sec["rows"]
        mt = B.metrics([r.get("engine|g") for r in rows])
        if mt.get("n"):
            print(f"  {sec['title']}：基因 top-1 {mt['top1']:.3f}、top-10 {mt['top10']:.3f}、MRR {mt['mrr']:.3f}（{mt['n']} 例）")
    return 0


def cmd_select_test(args) -> int:
    out = args.sets / "test.txt"
    if out.exists() and not args.force:
        print(f"{out} 已存在；測試組只抽一次並鎖定（第 11.3 節），不覆寫。")
        return 1
    cases = B.load_cases(args.cases)
    by_id = {c.id: c for c in cases}
    old = load_engine(args.old_snapshot, args.old_gene_map)
    cited = B.hpoa_citations(args.old_hpoa)
    elig = B.timesplit_eligible(cases, old, cited)
    used = set()
    for key in ("dev", "fresh"):
        if (args.sets / f"{key}.txt").is_file():
            used |= set(B.read_ids(args.sets / f"{key}.txt"))
    used_pmids = set().union(*(by_id[i].pmids for i in used if i in by_id)) if used else set()
    pool = [c for c in elig if c.id not in used and not (c.pmids & used_pmids)]
    picked = B.select(pool, per_cohort=args.per_cohort, seed=args.seed)
    header = (f"鎖定測試組（第 11.3 節）：time-split 合格（{old.snapshot}）{len(elig)} 例中，排除開發與第二組病例及其文獻，"
              f"剩 {len(pool)} 例；每個 cohort 最多 {args.per_cohort} 例，seed {args.seed}，共 {len(picked)} 例。"
              f"\n產生：{time.strftime('%Y-%m-%d')}，之後不得重抽或依結果調整參數。")
    B.write_ids(out, [c.id for c in picked], header)
    print(header)
    return 0


def legacy(argv) -> int:
    ap = argparse.ArgumentParser(description="quick mode")
    ap.add_argument("--snapshot", type=Path, required=True)
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--case-list", type=Path, default=None)
    ap.add_argument("--gene-map", type=Path, default=None)
    ap.add_argument("--rows", type=Path, default=None)
    args = ap.parse_args(argv)
    started = time.time()
    eng = load_engine(args.snapshot, args.gene_map)
    cases = B.load_cases(args.cases)
    if args.case_list:
        by_id = {c.id: c for c in cases}
        cases = [by_id[i] for i in B.read_ids(args.case_list) if i in by_id]
    rows, cov = B.evaluate(eng, cases, methods=("engine",))
    print(f"知識庫 {eng.snapshot}：{eng.kb.N:,} 個疾病、{len(eng.genes):,} 個基因；病例 {len(rows)}（略過 {cov['cases'] - len(rows)}）")
    for level, name in (("d", "疾病"), ("g", "基因")):
        mt = B.metrics([r.get(f"engine|{level}") for r in rows])
        if mt.get("n"):
            print(f"{name}  top1={mt['top1']:.3f}  top10={mt['top10']:.3f}  MRR={mt['mrr']:.3f}  "
                  f"median={mt['median']:.1f}  n={mt['n']}")
    print(f"耗時 {time.time() - started:.0f} 秒")
    if args.rows:
        args.rows.write_text(json.dumps(rows, ensure_ascii=False, indent=1), encoding="utf-8")
    return 0


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] not in ("run", "select-test"):
        return legacy(argv)
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "select-test"):
        p = sub.add_parser(name)
        p.add_argument("--cases", type=Path, required=True, help="phenopacket-store 資料夾")
        p.add_argument("--old-snapshot", type=Path, required=True, help="time-split 用的舊知識庫快照（2023-10-09）")
        p.add_argument("--old-gene-map", type=Path, default=None, help="舊版 genes_to_disease.txt")
        p.add_argument("--sets", type=Path, default=SETS)
        if name == "run":
            p.add_argument("--current-snapshot", type=Path, default=None)
            p.add_argument("--current-gene-map", type=Path, default=None)
            p.add_argument("--no-lopo", action="store_true")
            p.add_argument("--lirical", type=Path, default=None, help="每例一個 {case}.tsv 的 LIRICAL 結果")
            p.add_argument("--exomiser", type=Path, default=None)
            p.add_argument("--out", type=Path, required=True)
        else:
            p.add_argument("--old-hpoa", type=Path, required=True)
            p.add_argument("--per-cohort", type=int, default=3)
            p.add_argument("--seed", type=int, default=2026)
            p.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    return cmd_run(args) if args.cmd == "run" else cmd_select_test(args)


if __name__ == "__main__":
    sys.exit(main())
