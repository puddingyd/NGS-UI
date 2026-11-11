#!/usr/bin/env python3
"""Explain how the stored OMIM Update List matches the platform's content.

    PYTHONPATH=backend python scripts/omim_updates_check.py

Prints what was stored from the Update List pages, how many platform
gene / disease MIMs appear in it, and — per kind of platform content —
why an update does or does not count as "OMIM 有更新" (content date vs.
OMIM's listed date). Read-only.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))

from app.services import disease_picks, omim_capture, omim_db, omim_store, omim_updates  # noqa: E402


def main() -> int:
    if not omim_db.exists():
        print("平台 OMIM 資料庫不存在")
        return 1
    with omim_db.connect() as con:
        total = con.execute("SELECT COUNT(*) FROM omim_updates").fetchone()[0]
        kinds = con.execute("SELECT kind, COUNT(*) FROM omim_updates GROUP BY kind").fetchall()
        prefixes = con.execute("SELECT prefix, COUNT(*) FROM omim_updates GROUP BY prefix").fetchall()
        span = con.execute("SELECT MIN(listed_date), MAX(listed_date), COUNT(DISTINCT mim) FROM omim_updates").fetchone()
        months = con.execute("SELECT COUNT(*), SUM(complete) FROM omim_update_months").fetchone()
        sample = con.execute("SELECT * FROM omim_updates ORDER BY listed_date DESC LIMIT 5").fetchall()
    print(f"== 已存的 Update List：{total} 筆，{span[2]} 個不同 MIM，日期 {span[0]} ~ {span[1]}；月份 {months[0]}（完整 {months[1] or 0}）")
    print("   類型：", {k: n for k, n in kinds})
    print("   符號：", {(p or '(無)'): n for p, n in prefixes})
    for r in sample:
        print("   例：", dict(r))

    entry_latest, any_latest = omim_updates.latest_dates()
    cov = omim_updates.coverage()
    print(f"== 涵蓋：需要從 {cov['needed_from']} 起；連續涵蓋起點 {cov['coverage_start'] or '—'}；缺 {len(cov['missing'])} 個月")

    omim_store.ensure_loaded()
    rows = {id(r): r for r in list(omim_store._state.get("by_omim_id", {}).values())
            + list(omim_store._state.get("by_gene", {}).values())}.values()
    texts = omim_db.phenotype_text_sources()
    gene_caps = omim_db.gene_sync_captures()
    times = omim_capture.capture_times()
    stats = Counter()
    examples: dict[str, list[str]] = {}

    def note(key, text):
        stats[key] += 1
        if len(examples.setdefault(key, [])) < 5:
            examples[key].append(text)

    seen = set()
    for row in rows:
        gene = row.get("gene_symbol") or ""
        gene_mim = str(row.get("OMIM_id") or "").split(".")[0]
        done = (row.get("Legacy_done", row.get("Done")) or "")[:10]
        if gene_mim in entry_latest:
            stats["基因 MIM 出現在更新清單"] += 1
            if gene_mim in gene_caps:
                cid, at = gene_caps[gene_mim]
                day = omim_updates.us_date(times.get(cid) or at)
                key = "已同步基因頁：OMIM 更新晚於擷取" if entry_latest[gene_mim] > day else "已同步基因頁：擷取已含最新更新"
                note(key, f"{gene} {gene_mim} 擷取 {day} / OMIM {entry_latest[gene_mim]}")
            else:
                stats["基因頁未同步（不另加步驟）"] += 1
        meta = {m.get("slot"): m for m in row.get("Disease_meta") or []}
        for slot, text in disease_picks.disease_cells(row):
            pm = disease_picks.phenotype_mim(disease_picks.first_line(text))
            if not pm or (gene_mim, pm) in seen:
                continue
            seen.add((gene_mim, pm))
            stats["平台疾病（基因×疾病 MIM）"] += 1
            listed = any_latest.get(pm, "")
            if not listed:
                continue
            stats["疾病 MIM 出現在更新清單"] += 1
            lines = len([ln for ln in text.splitlines() if ln.strip()])
            source = (meta.get(slot) or {}).get("source", "")
            label = f"{gene} {pm} OMIM {listed}"
            if lines <= 1:
                note("只有名稱（本來就是缺說明）", label)
            elif source == "capture":
                src = texts.get(pm)
                day = max(omim_updates.us_date(times.get(src[0]) or src[1]), src[2] or "") if src else ""
                note("擷取採用：OMIM 更新晚於內容" if listed > day else "擷取採用：內容已含最新更新",
                     f"{label} / 內容 {day}")
            elif not done:
                note("OMIM.xlsx 原文：沒有 Done 日期（無法判斷）", label)
            elif listed > done:
                note("OMIM.xlsx 原文：OMIM 更新晚於 Done", f"{label} / Done {done}")
            else:
                note("OMIM.xlsx 原文：Done 之後沒有更新", f"{label} / Done {done}")
    print("== 比對結果")
    for key, n in stats.items():
        print(f"   {key}：{n}")
        for ex in examples.get(key, []):
            print(f"      例：{ex}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
