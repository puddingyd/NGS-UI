#!/usr/bin/env python3
"""Run the phenotype engine (plan §4.3–4.8, phase 2) and print a plain summary.

For a registered sample, the four explanation files are written next to the
analysis version's analysis.json (``analyses/{ver}/{SID}.pheno_engine_*.tsv``
and ``…_run.json``); the current phenotype score and ranking are not changed.

    .venv/bin/python scripts/pheno_engine.py --sample 檢體編號
    .venv/bin/python scripts/pheno_engine.py --hpo HP:0001250,HP:0001263=2   (=2 表示權重 2)

``--snapshot DIR`` uses a specific knowledge-base snapshot instead of the
active one.  Output that names a sample stays on this server.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "backend"))

from app.services.pheno_engine import engine as E  # noqa: E402


def summary(result: E.Result, top: int) -> str:
    eng = result.engine
    lines = [f"知識庫：{eng.snapshot}（{eng.kb.N:,} 個疾病、{len(eng.genes):,} 個基因）"]
    terms = result.pre["terms"]
    lines.append(f"計分的病人 HPO：{len(terms)} 個")
    for t in terms:
        note = "" if t["how"] == "current" else f"（由 {t['input_id']} 換算）"
        weight = f"，權重 {t['weight']:g}" if t["weight"] != 1 else ""
        lines.append(f"  {t['hpo_id']} {t['name']}{note}{weight}")
    if result.pre["removed"]:
        lines.append("未計分：")
        for r in result.pre["removed"]:
            kept = f"（保留 {r['kept']}）" if r.get("kept") and r["reason"] == "redundant_ancestor" else ""
            lines.append(f"  {r['hpo_id']}：{E.REMOVED_ZH.get(r['reason'], r['reason'])}{kept}")
    if result.empty:
        lines.append("沒有可計分的表型 term，無法排序。")
        return "\n".join(lines)
    lead = result.lead_nats()
    lines.append("")
    lines.append(f"基因排名前 {top}（P＝表型點數，暫定校準）：")
    if not eng.genes:
        lines.append("  （這個快照沒有疾病–基因對應，只能排疾病）")
    for g in result.top_genes(top):
        tags = [x for x, on in (("明顯領先", g["lead"]), ("接近同分", g["near_tie"] and g["rank"] > 1)) if on]
        lines.append(f"  {g['rank']:>3}. {g['gene']:<10} G={g['G']:6.2f}  P={g['P']:+5.1f}  "
                     f"{g['best_disease']} {g['best_disease_name'][:50]}（解釋 {g['explained']:.0%}）"
                     + (f"  [{'、'.join(tags)}]" if tags else ""))
    if lead is not None:
        lines.append(f"  第 1 名領先第 2 名 {lead:.2f} nats")
    lines.append("")
    lines.append(f"疾病排名前 {top}：")
    for d in result.top_diseases(top):
        genes = "、".join(d["genes"][:4]) or "（無對應基因）"
        lines.append(f"  {d['rank']:>3}. {d['disease_id']:<13} S*={d['S_star']:6.2f}  解釋 {d['explained']:.0%}  "
                     f"{d['name'][:50]}  {genes}")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sample", help="registered sample ID")
    ap.add_argument("--version", default=None, help="analysis version (default: active)")
    ap.add_argument("--hpo", help="comma-separated HPO IDs, optional =weight")
    ap.add_argument("--snapshot", type=Path, default=None, help="knowledge-base snapshot folder")
    ap.add_argument("--top", type=int, default=10)
    args = ap.parse_args(argv)
    if not args.sample and not args.hpo:
        ap.error("請給 --sample 或 --hpo")
    engine = E.load_engine(args.snapshot)
    if engine is None:
        print("找不到表型知識庫：請先在首頁按「更新 HPO」，或用 --snapshot 指定快照資料夾", file=sys.stderr)
        return 2
    if args.sample:
        try:
            out = E.run_for_analysis(args.sample, args.version, engine=engine)
        except (FileNotFoundError, RuntimeError) as exc:
            print(exc, file=sys.stderr)
            return 1
        print(summary(out["result"], args.top))
        print(f"\n已寫入 {len(out['paths'])} 個檔案到 {out['paths'][0].parent}")
        return 0
    items = []
    for token in args.hpo.split(","):
        token = token.strip()
        if not token:
            continue
        hid, _, w = token.partition("=") if "=" in token else (token, "", "1")
        items.append((hid.strip(), w or 1))
    print(summary(engine.score(items), args.top))
    return 0


if __name__ == "__main__":
    sys.exit(main())
