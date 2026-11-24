"""Before/after report for the cleaned HPO → gene table (plan §10.7, phase 0).

Compares the official ``phenotype_to_genes.txt`` with a cleaned snapshot from
``hpo_gene_table`` at two levels:

* table: term–gene pairs removed / added, terms that are not phenotypes but
  carried genes (inheritance, onset, …);
* cases: for every registered case's active analysis (HPO + panels), the
  phenotype-related gene set and the ranking of reviewer-marked SNVs
  (Causative / Other / Candidate) under ``total_score = geno + pheno`` with the
  live ACMG cascade (reviewer → In-house → GeneBe).  Only the phenotype part
  differs between the two runs.

Read-only.  ``report.md`` / ``summary.json`` / ``term_changes.tsv`` carry no
sample IDs; ``case_changes.tsv`` and ``marked_variants.tsv`` do (server only).
"""
from __future__ import annotations

import csv
import json
import statistics
import time
from collections import Counter, defaultdict
from pathlib import Path

from . import hpo_gene_table

STATUS_LABEL = {"1": "第一類（Causative）", "2": "第二類（Other）", "C": "Candidate"}


# ── scoring with an explicit table ─────────────────────────────────────────
class Table:
    def __init__(self, gene_table: Path, term_status: Path | None, canonical, label: str):
        self.label = label
        self.genes = hpo_gene_table.read_gene_table(gene_table, canonical)
        self.aliases, self.non_phenotype = hpo_gene_table.read_term_status(term_status)

    def norm(self, hpo_id: str) -> str:
        key = (hpo_id or "").strip()
        return self.aliases.get(key, key)

    def score(self, hpo_list, panel_pairs, panels: dict[str, set[str]]):
        """Same formula as phenotype_scorer.compute_pheno_score, with this table."""
        pairs, skipped = [], []
        for entry in hpo_list or []:
            if isinstance(entry, dict):
                hid = (entry.get("phenotype") or entry.get("hpo_id") or "").strip()
                try:
                    weight = float(entry.get("weight", 1) or 1)
                except (TypeError, ValueError):
                    weight = 1.0
            else:
                hid, weight = entry[0], float(entry[1])
            hid = self.norm(hid)
            if not hid:
                continue
            if hid in self.non_phenotype:
                skipped.append(hid)
                continue
            pairs.append((hid, weight, self.genes.get(hid, ())))
        for name, weight in panel_pairs:
            pairs.append((name, weight, panels.get(name, ())))
        total = sum(w for _, w, _ in pairs)
        accum: dict[str, float] = defaultdict(float)
        for _, weight, genes in pairs:
            for gene in genes:
                accum[gene] += weight
        scores = {g: 100.0 * w / total for g, w in accum.items() if w > 0} if total > 0 else {}
        return scores, skipped


def _panel_pairs(panels_list, panels: dict[str, set[str]], resolve) -> list[tuple[str, float]]:
    out = []
    for entry in panels_list or []:
        if isinstance(entry, dict):
            name = (entry.get("name") or "").strip()
            try:
                weight = float(entry.get("weight", 1) or 1)
            except (TypeError, ValueError):
                weight = 1.0
        elif isinstance(entry, str):
            name, weight = entry.strip(), 1.0
        else:
            name, weight = entry[0], float(entry[1])
        name = resolve(name)
        if name in panels:
            out.append((name, weight))
    return out


# ── table level ────────────────────────────────────────────────────────────
def compare_tables(old: Table, new: Table, names: dict[str, str]) -> tuple[dict, list[dict]]:
    rows = []
    removed = added = 0
    for term in sorted(set(old.genes) | set(new.genes)):
        a, b = old.genes.get(term, set()), new.genes.get(term, set())
        lost, gain = len(a - b), len(b - a)
        removed += lost
        added += gain
        if lost or gain:
            status = "non_phenotype" if term in new.non_phenotype else (
                "phenotype" if term in names else "unknown")
            rows.append({"hpo_id": term, "name": names.get(term, ""), "status": status,
                         "old_genes": len(a), "new_genes": len(b), "removed": lost, "added": gain})
    rows.sort(key=lambda r: (-r["removed"], r["hpo_id"]))
    summary = {
        "old_terms": len(old.genes), "new_terms": len(new.genes),
        "old_pairs": sum(len(g) for g in old.genes.values()),
        "new_pairs": sum(len(g) for g in new.genes.values()),
        "pairs_removed": removed, "pairs_added": added,
        "terms_changed": len(rows),
        "non_phenotype_terms_with_genes": sum(1 for r in rows if r["status"] == "non_phenotype" and r["old_genes"]),
    }
    return summary, rows


# ── case level ─────────────────────────────────────────────────────────────
def _rank(scores: dict[str, float], vid: str) -> float:
    mine = scores[vid]
    higher = sum(1 for v, s in scores.items() if v != vid and s > mine)
    ties = sum(1 for v, s in scores.items() if v != vid and s == mine)
    return 1 + higher + ties / 2.0


def _rank_summary(ranks: list[float]) -> dict:
    if not ranks:
        return {"n": 0}
    n = len(ranks)
    return {"n": n, "median": statistics.median(ranks),
            "top1": round(sum(r <= 1 for r in ranks) / n, 4), "top5": round(sum(r <= 5 for r in ranks) / n, 4),
            "top10": round(sum(r <= 10 for r in ranks) / n, 4), "top20": round(sum(r <= 20 for r in ranks) / n, 4),
            "mrr": round(sum(1.0 / r for r in ranks) / n, 4)}


def _case_phenotype(sample_id: str):
    """(test_type, hpo_list, panels_list) of the active analysis, read-only."""
    from . import analyses_store, sample_layout, sample_loader
    from .acmg_engine import runner
    try:
        meta = json.loads(sample_layout.state_file(sample_id, "sample_metadata.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(meta, dict):
        return None
    version = sample_loader._resolve_version(sample_id, None, meta.get("active_analysis"))
    if version is not None:
        analysis = analyses_store.read_version(sample_id, version) or {}
        hpo, panels = analysis.get("hpo") or [], analysis.get("selected_panels") or []
    else:
        hpo = meta.get("hpo") or meta.get("patient_phenotype") or []
        panels = meta.get("selected_panels") or []
    return runner.case_inputs(sample_id)["test_type"], meta, hpo, panels


def _case_variants(sample_id: str, meta: dict):
    """{vid: (gene, points, hgvs)} for review variants + {vid: status} of marked ones; None if no review."""
    from .acmg_engine import runner
    from .acmg_engine import parse as P
    from .acmg_engine.engine import _primary_row
    from .acmg_engine.replacement_diff import (
        _bulk_manual, _num, _sample_snapshot, effective, manual_choice, marked_statuses)
    from .acmg_engine.validation import class_code, first_value, normalize_variant_id

    case = runner.case_inputs(sample_id)
    if not case["review"].is_file():
        return None
    statuses = marked_statuses(meta)
    if not statuses:
        return {}, {}
    edits = meta.get("edits") if isinstance(meta.get("edits"), dict) else {}
    results = dict(runner.iter_results(case["existing"], detail=False)) if Path(case["existing"]).is_file() else {}
    with open(case["review"], encoding="utf-8", errors="replace", newline="") as handle:
        groups = runner.group_variant_rows(csv.DictReader(handle, delimiter="\t"))
    by_norm = {normalize_variant_id(v): v for v in groups}
    statuses = {by_norm[n]: s for vid, s in statuses.items() if (n := normalize_variant_id(vid)) in by_norm}
    edits = {by_norm.get(normalize_variant_id(k), k): v for k, v in edits.items()}
    global_manual = {normalize_variant_id(k): v for k, v in (_bulk_manual(groups.keys()) or {}).items()}
    out = {}
    for vid, rows in groups.items():
        row = _primary_row(rows)
        manual = manual_choice(_sample_snapshot(edits.get(vid)), global_manual.get(normalize_variant_id(vid)))
        result = results.get(vid)
        inhouse = ({"class": class_code(result.get("classification")), "points": _num(result.get("points"))}
                   if result else {"class": class_code(first_value(rows, "ACMG_CLASS")),
                                   "points": _num(first_value(rows, "ACMG_SCORE") or first_value(rows, "ACMG_POINTS"))})
        if inhouse.get("class") is None and inhouse.get("points") is None:
            inhouse = None
        genebe = ({"source": "genebe", "class": class_code(first_value(rows, "GENEBE_ACMG_CLASS")),
                   "points": _num(first_value(rows, "GENEBE_ACMG_SCORE"))}
                  if first_value(rows, "GENEBE_ACMG_CLASS") else None)
        chosen = effective(manual, dict(inhouse, source="inhouse") if inhouse else None, genebe, "genebe")
        out[vid] = (P.text(row.get("GENE")), chosen["points"],
                    " ".join(x for x in (P.text(row.get("HGVS_C")).split(":")[-1],
                                         P.text(row.get("HGVS_P")).split(":")[-1]) if x))
    return out, statuses


def run(*, out_dir: Path, old_table: Path, new_snapshot: Path, sample_ids=None, progress=None,
        case_phenotype=_case_phenotype, case_variants=_case_variants) -> dict:
    from . import panel_deadzone, phenotype_scorer, sample_layout
    from .acmg_engine.validation import geno_score

    started = time.time()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def canonical(gene: str) -> str:
        return panel_deadzone.canonical_panel_gene_symbol(gene) or (gene or "").strip()

    old = Table(old_table, None, canonical, "phenotype_to_genes.txt")
    new = Table(Path(new_snapshot) / hpo_gene_table.GENE_TABLE, Path(new_snapshot) / hpo_gene_table.TERM_STATUS,
                canonical, Path(new_snapshot).name)
    names = {}
    with open(Path(new_snapshot) / hpo_gene_table.TERM_STATUS, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            names[row["hpo_id"]] = row.get("name") or ""
    table_summary, term_rows = compare_tables(old, new, names)
    panels = phenotype_scorer._load_panels()

    sample_ids = sorted(set(sample_layout.iter_sample_ids())) if sample_ids is None else list(sample_ids)
    cases = Counter()
    term_use = Counter()
    skipped_use = Counter()
    gene_counts = {"old": [], "new": []}
    shrink = []
    case_rows, marked_rows = [], []
    marked = defaultdict(lambda: {"ranks_old": [], "ranks_new": [], "lost_panel": [], "gained_panel": 0,
                                  "in_panel_old": 0, "in_panel_new": 0, "move": Counter()})
    lost_detail = Counter()

    for index, sid in enumerate(sample_ids, 1):
        if progress and (index == 1 or index % 25 == 0 or index == len(sample_ids)):
            progress(f"[{index}/{len(sample_ids)}] {sid}")
        info = case_phenotype(sid)
        if info is None:
            cases["skipped_unregistered"] += 1
            continue
        tt, meta, hpo, panels_list = info
        pairs = _panel_pairs(panels_list, panels, phenotype_scorer._resolve_panel_name)
        if not hpo and not pairs:
            cases["skipped_no_phenotype"] += 1
            continue
        cases[tt] += 1
        s_old, _ = old.score(hpo, pairs, panels)
        s_new, skipped = new.score(hpo, pairs, panels)
        for entry in hpo:
            hid = new.norm((entry.get("phenotype") or entry.get("hpo_id") or "") if isinstance(entry, dict) else entry[0])
            term_use[hid] += 1
        for hid in set(skipped):
            skipped_use[hid] += 1
        if skipped:
            cases["with_excluded_terms"] += 1
        a, b = set(s_old), set(s_new)
        gene_counts["old"].append(len(a))
        gene_counts["new"].append(len(b))
        if a:
            shrink.append((len(a) - len(b)) / len(a))
        case_rows.append({"sample": sid, "test_type": tt, "hpo_terms": len(hpo), "panels": len(pairs),
                          "excluded_terms": ";".join(sorted(set(skipped))),
                          "genes_old": len(a), "genes_new": len(b), "lost": len(a - b), "gained": len(b - a)})

        variants = case_variants(sid, meta)
        if variants is None:
            cases["no_review_tsv"] += 1
            continue
        var, statuses = variants
        if not statuses:
            continue
        tot_old, tot_new = {}, {}
        for vid, (gene, points, _) in var.items():
            g = canonical(gene)
            geno = geno_score(points) if points is not None else 0
            tot_old[vid] = geno + s_old.get(g, 0.0)
            tot_new[vid] = geno + s_new.get(g, 0.0)
        for vid, status in statuses.items():
            gene, points, hgvs = var[vid]
            g = canonical(gene)
            p_old, p_new = s_old.get(g, 0.0), s_new.get(g, 0.0)
            r_old, r_new = _rank(tot_old, vid), _rank(tot_new, vid)
            bucket = marked[(tt, status)]
            bucket["ranks_old"].append(r_old)
            bucket["ranks_new"].append(r_new)
            bucket["in_panel_old"] += p_old > 0
            bucket["in_panel_new"] += p_new > 0
            bucket["move"]["better" if r_new < r_old else ("worse" if r_new > r_old else "same")] += 1
            if p_old > 0 and p_new == 0:
                terms = [new.norm(e.get("phenotype") or e.get("hpo_id") or "") for e in hpo if isinstance(e, dict)]
                causes = sorted(t for t in terms if g in old.genes.get(t, ()) and g not in new.genes.get(t, ()))
                bucket["lost_panel"].append({"gene": gene, "hgvs": hgvs, "causes": causes})
                lost_detail[(tt, status, gene, hgvs, ";".join(causes))] += 1
            elif p_old == 0 and p_new > 0:
                bucket["gained_panel"] += 1
            marked_rows.append({"sample": sid, "test_type": tt, "status": status, "gene": gene, "hgvs": hgvs,
                                "variant": vid, "pheno_old": round(p_old, 2), "pheno_new": round(p_new, 2),
                                "rank_old": r_old, "rank_new": r_new, "n_variants": len(var)})

    report = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "old_table": str(old_table), "new_snapshot": str(new_snapshot),
        "elapsed_seconds": round(time.time() - started, 1),
        "table": table_summary,
        "cases": dict(cases),
        "excluded_term_use": [{"hpo_id": t, "name": names.get(t, ""), "cases": n,
                               "old_genes": len(old.genes.get(t, ()))} for t, n in skipped_use.most_common()],
        "gene_set": {
            "old_median": statistics.median(gene_counts["old"]) if gene_counts["old"] else None,
            "new_median": statistics.median(gene_counts["new"]) if gene_counts["new"] else None,
            "shrink_median": round(statistics.median(shrink), 4) if shrink else None,
            "shrink_over_50pct": sum(1 for x in shrink if x > 0.5),
        },
        "marked": {f"{tt}|{status}": {
            "rank_old": _rank_summary(b["ranks_old"]), "rank_new": _rank_summary(b["ranks_new"]),
            "in_panel_old": b["in_panel_old"], "in_panel_new": b["in_panel_new"],
            "lost_panel": len(b["lost_panel"]), "gained_panel": b["gained_panel"], "move": dict(b["move"])}
            for (tt, status), b in sorted(marked.items())},
        "lost_panel_variants": [{"test_type": tt, "status": st, "gene": gene, "hgvs": hgvs,
                                 "causes": [f"{c} {names.get(c, '')}".strip() for c in causes.split(";") if c],
                                 "cases": n}
                                for (tt, st, gene, hgvs, causes), n in sorted(lost_detail.items())],
        "top_terms_changed": [r for r in term_rows if term_use.get(r["hpo_id"])][:30],
    }
    for r in report["top_terms_changed"]:
        r["cases"] = term_use.get(r["hpo_id"], 0)
    (out_dir / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    (out_dir / "report.md").write_text(render_markdown(report), encoding="utf-8")
    _write_tsv(out_dir / "term_changes.tsv", [dict(r, cases=term_use.get(r["hpo_id"], 0)) for r in term_rows])
    _write_tsv(out_dir / "case_changes.tsv", case_rows)
    _write_tsv(out_dir / "marked_variants.tsv", marked_rows)
    return report


def _write_tsv(path: Path, rows: list[dict]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        if not rows:
            handle.write("")
            return
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _pct(x) -> str:
    return "—" if x is None else f"{100 * x:.1f}%"


def render_markdown(report: dict) -> str:
    t = report["table"]
    out = ["# 表型分數查表修正前後差異（階段 0）", "",
           f"- 產生時間：{report['generated']}（{report['elapsed_seconds']} 秒）",
           f"- 舊表：{report['old_table']}", f"- 新表：{report['new_snapshot']}",
           "- 本檔與 summary.json、term_changes.tsv 不含檢體編號；case_changes.tsv、marked_variants.tsv 含檢體編號，只留在伺服器上。",
           "", "## 1. 查表本身", "",
           "| 項目 | 舊表 | 新表 |", "|---|---|---|",
           f"| 有基因的 HPO term | {t['old_terms']:,} | {t['new_terms']:,} |",
           f"| HPO term–基因配對 | {t['old_pairs']:,} | {t['new_pairs']:,} |", "",
           f"移除 {t['pairs_removed']:,} 個配對、新增 {t['pairs_added']:,} 個；"
           f"{t['non_phenotype_terms_with_genes']:,} 個非表型 term（遺傳方式、發病年齡等）在舊表有基因，新表不再參與計分。"]
    c = report["cases"]
    g = report["gene_set"]
    n_cases = sum(v for k, v in c.items() if k in ("WES", "WGS"))
    out += ["", "## 2. 個案的表型相關基因", "",
            f"- 比較 {n_cases} 個有 HPO／panel 的個案（" + "、".join(
                f"{k} {v}" for k, v in c.items() if k in ("WES", "WGS")) + "）；"
            f"略過未登錄 {c.get('skipped_unregistered', 0)}、沒有 HPO／panel {c.get('skipped_no_phenotype', 0)}。",
            f"- 表型相關基因數（分數 > 0）中位數：舊 {g['old_median']} → 新 {g['new_median']}；"
            f"每個個案減少比例中位數 {_pct(g['shrink_median'])}，減少超過一半的個案 {g['shrink_over_50pct']} 個。",
            f"- {c.get('with_excluded_terms', 0)} 個個案選了不參與計分的 term（新表把它們移出分母，其餘 term 的權重比例因此提高）："]
    for item in report["excluded_term_use"][:15]:
        out.append(f"  - {item['hpo_id']} {item['name']}：{item['cases']} 個個案（舊表 {item['old_genes']:,} 個基因）")
    out += ["", "## 3. 判讀者標記的 SNV", "",
            "名次＝個案 review TSV 內所有變異依 total_score（有效 ACMG 點數換算 + 表型分數）排序；有效 ACMG 依判讀者 → In-house → GeneBe，"
            "不含院內常見扣分；兩次只有表型分數不同。", "",
            "| 類別 | 個數 | 臨床表現相關（舊→新） | 失去 | 新增 | 名次中位數（舊→新） | top10（舊→新） | top20（舊→新） | MRR（舊→新） | 名次較好／相同／較差 |",
            "|---|---|---|---|---|---|---|---|---|---|"]
    for key, m in report["marked"].items():
        tt, st = key.split("|")
        ro, rn = m["rank_old"], m["rank_new"]
        if not ro.get("n"):
            continue
        mv = m["move"]
        out.append(f"| {tt} {STATUS_LABEL.get(st, st)} | {ro['n']} | {m['in_panel_old']}→{m['in_panel_new']} | "
                   f"{m['lost_panel']} | {m['gained_panel']} | {ro['median']}→{rn['median']} | "
                   f"{_pct(ro['top10'])}→{_pct(rn['top10'])} | {_pct(ro['top20'])}→{_pct(rn['top20'])} | "
                   f"{ro['mrr']}→{rn['mrr']} | {mv.get('better', 0)}／{mv.get('same', 0)}／{mv.get('worse', 0)} |")
    lost = report["lost_panel_variants"]
    out += ["", "### 失去「臨床表現相關基因」的標記變異", ""]
    if not lost:
        out.append("沒有。")
    for item in lost:
        causes = "、".join(item["causes"]) or "（其他 term 的基因清單變動）"
        out.append(f"- {item['test_type']} {STATUS_LABEL.get(item['status'], item['status'])} {item['gene']} "
                   f"{item['hgvs']}（{item['cases']} 個案）：舊表因 {causes} 帶入")
    out += ["", "## 4. 個案常用且基因變動最多的 term（前 30）", "",
            "| HPO | 名稱 | 狀態 | 使用個案 | 舊基因數 | 新基因數 | 移除 | 新增 |", "|---|---|---|---|---|---|---|---|"]
    for r in report["top_terms_changed"]:
        out.append(f"| {r['hpo_id']} | {r['name']} | {r['status']} | {r['cases']} | {r['old_genes']:,} | "
                   f"{r['new_genes']:,} | {r['removed']:,} | {r['added']:,} |")
    return "\n".join(out) + "\n"
