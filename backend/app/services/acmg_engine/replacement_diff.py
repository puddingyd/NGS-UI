"""Before/after report for replacing the pipeline ACMG with the rule engine (plan §5.2.2, step 6).

Today the platform's effective ACMG for an SNV is chosen in this order::

    reviewer (per-sample snapshot) → reviewer (global reusable) → GeneBe → pipeline ACMG_*

and the same points drive tier 1C (points ≥ 4) and ``total_score``
(``geno_score(points) + pheno_score``).  Step 7 replaces the last slot with the
rule engine.  This module simulates, for every published case, what would
change; nothing is written to any case.

Scenario A is what step 7 adopted (2026-11-24): the engine becomes the
"In-house" ACMG and moves ahead of GeneBe — reviewer → *In-house (engine)* →
GeneBe.  Class, tier 1C and total_score (the ranking V) all follow these
points, so every variant with an engine result can change, not only those
that used to fall through to the pipeline.  Scenario B (ranking from
reviewer → engine only, GeneBe never used) is kept in ``summary.json`` for
reference; with engine results present it equals A.

``report.md`` / ``summary.json`` hold aggregate numbers only.
``marked_variants.tsv`` lists the reviewer-marked variants (status 1/2/C) *with
sample IDs* so they can be checked case by case on the server; do not copy it
out of the server or paste it anywhere.
"""
from __future__ import annotations

import csv
import json
import statistics
import time
from collections import Counter, defaultdict
from pathlib import Path

from . import parse as P
from . import runner
from .engine import _primary_row
from .validation import class_code, first_value, geno_score, normalize_variant_id, _default_pheno, _canonical

CLASS_ORDER = {"B": 0, "LB": 1, "VUS": 2, "LP": 3, "P": 4}
SOURCES = ("manual_sample", "manual_global", "genebe", "pipeline", "engine", "none")
_POINT_COLUMNS = ("GENEBE_ACMG_SCORE", "ACMG_SCORE", "ACMG_POINTS")


# ── per-variant helpers ─────────────────────────────────────────────────────
def _num(value):
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return P.first_number(str(value))


def manual_choice(sample_snapshot: dict | None, global_assertion: dict | None) -> dict | None:
    """Same precedence as ``sample_loader._apply_effective_acmg`` for the two reviewer slots."""
    if sample_snapshot:
        return {"source": "manual_sample", "class": class_code(sample_snapshot.get("classification")),
                "points": _num(sample_snapshot.get("score"))}
    if global_assertion and (global_assertion.get("reusable_criteria") or not global_assertion.get("criteria")):
        return {"source": "manual_global", "class": class_code(global_assertion.get("reusable_classification")),
                "points": _num(global_assertion.get("reusable_score"))}
    return None


def effective(manual: dict | None, genebe: dict | None, last: dict | None, last_source: str) -> dict:
    if manual:
        return manual
    if genebe:
        return genebe
    if last and (last.get("class") or last.get("points") is not None):
        return dict(last, source=last_source)
    return {"source": "none", "class": None, "points": None}


def base_tier(row: dict) -> str:
    """Tier without any ACMG points: 1A (ClinVar), 1B (LOFTEE HC), 1C (predictors only) or 2."""
    from ...adapters import snv_tsv
    stripped = dict(row)
    for column in _POINT_COLUMNS:
        stripped[column] = ""
    return snv_tsv.classify_tier(stripped)


def tier_with(base: str, points) -> str:
    """Mirror of the manual-overlay tier rule: 1A/1B stay; points ≥ 4 or a predictor trigger → 1C."""
    if base in ("1A", "1B", "1C"):
        return base
    return "1C" if points is not None and points >= 4 else "2"


def total_score(points, pheno: float) -> float:
    return (geno_score(points) if points is not None else 0) + pheno


def _rank(scores: dict[str, float], vid: str) -> float:
    mine = scores[vid]
    higher = sum(1 for v, s in scores.items() if v != vid and s > mine)
    ties = sum(1 for v, s in scores.items() if v != vid and s == mine)
    return 1 + higher + ties / 2.0


def _top(scores: dict[str, float], n: int) -> set[str]:
    return {v for v, _ in sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))[:n]}


def _rank_summary(ranks: list[float]) -> dict:
    if not ranks:
        return {"n": 0}
    n = len(ranks)
    return {"n": n, "median": statistics.median(ranks),
            "top1": round(sum(r <= 1 for r in ranks) / n, 4), "top5": round(sum(r <= 5 for r in ranks) / n, 4),
            "top10": round(sum(r <= 10 for r in ranks) / n, 4), "top20": round(sum(r <= 20 for r in ranks) / n, 4),
            "mrr": round(sum(1.0 / r for r in ranks) / n, 4)}


def _direction(before: str | None, after: str | None) -> str:
    if before == after:
        return "same"
    if before in CLASS_ORDER and after in CLASS_ORDER:
        return "up" if CLASS_ORDER[after] > CLASS_ORDER[before] else "down"
    return "gained" if before is None else ("lost" if after is None else "other")


def marked_statuses(meta: dict) -> dict[str, str]:
    """{variant_id: '1'|'2'|'C'} for reviewer-marked SNVs (Causative wins over Other over Candidate)."""
    out = {}
    for vid, raw in (meta.get("status") or {}).items():
        parts = {p.strip().upper() for p in str(raw or "").split(",")}
        for code in ("1", "2", "C"):
            if code in parts:
                out[vid] = code
                break
    return out


def _load_meta(sample_id: str) -> dict:
    from .. import sample_layout
    try:
        data = json.loads(sample_layout.state_file(sample_id, "sample_metadata.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def _sample_snapshot(edit) -> dict | None:
    from .. import sample_loader
    return sample_loader._legacy_manual_snapshot(edit or {})


def _bulk_manual(variant_ids) -> dict:
    from .. import manual_acmg
    try:
        return manual_acmg.bulk_current("hg38", variant_ids)
    except Exception:
        return {}


# ── main pass ──────────────────────────────────────────────────────────────
def run(*, out_dir: Path, sample_ids=None, refs=None, meta_fn=_load_meta, manual_fn=_bulk_manual,
        snapshot_fn=_sample_snapshot, pheno_fn=_default_pheno, progress=None) -> dict:
    from .. import sample_layout

    started = time.time()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_ids = sorted(set(sample_layout.iter_sample_ids())) if sample_ids is None else list(sample_ids)

    cases = Counter()
    versions = Counter()
    by_tt = defaultdict(lambda: {
        "occurrences": 0, "source_before": Counter(), "source_after": Counter(),
        "class_transition": Counter(), "direction": Counter(), "tier_transition": Counter(),
        "engine_missing": 0})
    marked = defaultdict(lambda: {
        "n": 0, "by_status": Counter(), "class_changed": Counter(), "class_transition": Counter(),
        "tier_transition": Counter(), "report_class_changed": 0,
        "ranks": {"before": [], "A": [], "B": []}, "rank_move": {"A": Counter(), "B": Counter()}})
    top_overlap = defaultdict(lambda: {"A": [], "B": []})
    tier_variants: dict[tuple, dict] = {}
    marked_rows: list[dict] = []
    pheno_missing = 0

    for index, sid in enumerate(sample_ids, 1):
        if progress and (index == 1 or index % 25 == 0 or index == len(sample_ids)):
            progress(f"[{index}/{len(sample_ids)}] {sid}")
        case = runner.case_inputs(sid)
        if not case["review"].is_file() or not Path(case["existing"]).is_file():
            cases["skipped_no_results"] += 1
            continue
        if refs is not None and not runner.is_current(case["existing"], case["review"], refs):
            cases["skipped_stale"] += 1
            continue
        tt = case["test_type"]
        cases[tt] += 1
        meta = meta_fn(sid) or {}
        statuses = marked_statuses(meta)
        edits = meta.get("edits") if isinstance(meta.get("edits"), dict) else {}
        results = dict(runner.iter_results(case["existing"], detail=False))
        with open(case["review"], encoding="utf-8", errors="replace", newline="") as handle:
            groups = runner.group_variant_rows(csv.DictReader(handle, delimiter="\t"))
        global_manual = {normalize_variant_id(k): v for k, v in (manual_fn(groups.keys()) or {}).items()}
        by_norm = {normalize_variant_id(v): v for v in groups}
        statuses = {by_norm[n]: s for vid, s in statuses.items() if (n := normalize_variant_id(vid)) in by_norm}
        edits = {by_norm.get(normalize_variant_id(k), k): v for k, v in edits.items()}
        want_rank = bool(set(statuses) & set(groups))
        pheno = pheno_fn(sid) if want_rank else {}
        if want_rank and pheno is None:
            pheno_missing += 1
            pheno = {}
        scores = {"before": {}, "A": {}, "B": {}}
        info: dict[str, dict] = {}
        stats = by_tt[tt]

        for vid, rows in groups.items():
            row = _primary_row(rows)
            gene = P.text(row.get("GENE"))
            result = results.get(vid)
            manual = manual_choice(snapshot_fn(edits.get(vid)), global_manual.get(normalize_variant_id(vid)))
            gb_class = class_code(first_value(rows, "GENEBE_ACMG_CLASS"))
            genebe = ({"source": "genebe", "class": gb_class, "points": _num(first_value(rows, "GENEBE_ACMG_SCORE"))}
                      if first_value(rows, "GENEBE_ACMG_CLASS") else None)
            pipeline = {"class": class_code(first_value(rows, "ACMG_CLASS")),
                        "points": _num(first_value(rows, "ACMG_SCORE") or first_value(rows, "ACMG_POINTS"))}
            engine = ({"class": class_code(result.get("classification")), "points": _num(result.get("points"))}
                      if result else None)
            if result:
                versions[result.get("version") or "?"] += 1
            else:
                stats["engine_missing"] += 1
            before = effective(manual, genebe, pipeline, "pipeline")
            after = effective(manual, dict(engine, source="engine") if engine else None, genebe, "genebe")
            ranking_b = effective(manual, None, engine, "engine")
            base = base_tier(row)
            tier_before, tier_after = tier_with(base, before["points"]), tier_with(base, after["points"])

            stats["occurrences"] += 1
            stats["source_before"][before["source"]] += 1
            stats["source_after"][after["source"]] += 1
            if True:  # every occurrence can change under the adopted cascade
                stats["class_transition"][f"{before['class'] or '—'}→{after['class'] or '—'}"] += 1
                stats["direction"][_direction(before["class"], after["class"])] += 1
            stats["tier_transition"][f"{tier_before}→{tier_after}"] += 1
            if tier_before != tier_after:
                key = (tt, P.variant_key(row))
                item = tier_variants.setdefault(key, {
                    "test_type": tt, "variant": P.variant_key(row), "gene": gene,
                    "hgvs_c": P.text(row.get("HGVS_C")).split(":")[-1],
                    "tier_before": tier_before, "tier_after": tier_after,
                    "class_before": before["class"], "class_after": after["class"],
                    "points_before": before["points"], "points_after": after["points"], "cases": 0})
                item["cases"] += 1

            if want_rank:
                pheno_score = float(pheno.get(_canonical(gene), 0) or 0)
                scores["before"][vid] = total_score(before["points"], pheno_score)
                scores["A"][vid] = total_score(after["points"], pheno_score)
                scores["B"][vid] = total_score(ranking_b["points"], pheno_score)
            if vid in statuses:
                info[vid] = {"row": row, "gene": gene, "before": before, "after": after, "B": ranking_b,
                             "tier_before": tier_before, "tier_after": tier_after}

        if want_rank:
            top_before = _top(scores["before"], 20)
            for scenario in ("A", "B"):
                if top_before:
                    top_overlap[tt][scenario].append(len(top_before & _top(scores[scenario], 20)) / len(top_before))
        for vid, status in statuses.items():
            item = info.get(vid)
            if not item:
                continue
            m = marked[tt]
            m["n"] += 1
            m["by_status"][status] += 1
            before, after = item["before"], item["after"]
            if before["class"] != after["class"]:
                m["class_changed"][status] += 1
                m["class_transition"][f"{before['class'] or '—'}→{after['class'] or '—'}"] += 1
                if status in ("1", "2"):
                    m["report_class_changed"] += 1
            m["tier_transition"][f"{item['tier_before']}→{item['tier_after']}"] += 1
            ranks = {s: _rank(scores[s], vid) for s in ("before", "A", "B")}
            for s, r in ranks.items():
                m["ranks"][s].append(r)
            for s in ("A", "B"):
                move = "better" if ranks[s] < ranks["before"] else ("worse" if ranks[s] > ranks["before"] else "same")
                m["rank_move"][s][move] += 1
            row = item["row"]
            marked_rows.append({
                "sample": sid, "test_type": tt, "status": status, "variant": P.variant_key(row),
                "gene": item["gene"], "hgvs_c": P.text(row.get("HGVS_C")).split(":")[-1],
                "hgvs_p": P.text(row.get("HGVS_P")).split(":")[-1],
                "source_before": before["source"], "class_before": before["class"], "points_before": before["points"],
                "source_after": after["source"], "class_after": after["class"], "points_after": after["points"],
                "ranking_points_B": item["B"]["points"],
                "tier_before": item["tier_before"], "tier_after": item["tier_after"],
                "rank_before": ranks["before"], "rank_A": ranks["A"], "rank_B": ranks["B"],
                "class_changed": before["class"] != after["class"]})

    report = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "seconds": round(time.time() - started),
        "cases": dict(cases),
        "engine_versions": dict(versions),
        "pheno_unavailable_cases": pheno_missing,
        "by_test_type": {tt: {k: (dict(v) if isinstance(v, Counter) else v) for k, v in s.items()}
                         for tt, s in by_tt.items()},
        "marked": {tt: {
            "n": m["n"], "by_status": dict(m["by_status"]), "class_changed": dict(m["class_changed"]),
            "class_transition": dict(m["class_transition"]), "tier_transition": dict(m["tier_transition"]),
            "report_class_changed": m["report_class_changed"],
            "ranks": {s: _rank_summary(v) for s, v in m["ranks"].items()},
            "rank_move": {s: dict(c) for s, c in m["rank_move"].items()}} for tt, m in marked.items()},
        "top20_overlap": {tt: {s: (round(statistics.mean(v), 4) if v else None) for s, v in d.items()}
                          | {"cases": len(d["A"])} for tt, d in top_overlap.items()},
    }
    _write_tsv(out_dir / "tier_changes.tsv", sorted(tier_variants.values(), key=lambda r: (-r["cases"], r["variant"])))
    _write_tsv(out_dir / "marked_variants.tsv", marked_rows)
    (out_dir / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2, default=str),
                                          encoding="utf-8")
    (out_dir / "report.md").write_text(render_markdown(report), encoding="utf-8")
    return report


def _write_tsv(path: Path, rows: list[dict]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        if not rows:
            handle.write("")
            return
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if v is None else v) for k, v in row.items()})


# ── markdown ───────────────────────────────────────────────────────────────
_SOURCE_NAMES = {"manual_sample": "判讀者（本個案）", "manual_global": "判讀者（跨個案共用）", "genebe": "GeneBe",
                 "pipeline": "pipeline ACMG", "engine": "新 In-house", "none": "沒有分類"}
_SCENARIOS = {"before": "取代前（GeneBe → pipeline）", "A": "取代後（新 In-house → GeneBe）", "B": "B（只用新 In-house）"}
_STATUS = {"1": "Causative（1）", "2": "Other（2）", "C": "Candidate（C）"}
_DIRECTION = {"same": "不變", "up": "往致病方向", "down": "往良性方向", "gained": "原本沒有分類",
              "lost": "變成沒有分類", "other": "其他"}


def _pct(num, den) -> str:
    return f"{100.0 * num / den:.1f}%" if den else "—"


def _counter_line(counter: dict, total: int | None = None, limit: int = 12) -> str:
    items = sorted(counter.items(), key=lambda kv: -kv[1])[:limit]
    if not items:
        return "—"
    if total:
        return "、".join(f"{k} {v}（{_pct(v, total)}）" for k, v in items)
    return "、".join(f"{k} {v}" for k, v in items)


def render_markdown(report: dict) -> str:
    out = ["# ACMG 規則引擎取代 pipeline ACMG：前後差異（模擬）", ""]
    out.append(f"- 產生時間：{report['generated_at']}；耗時 {report['seconds']} 秒")
    out.append(f"- 個案：{report['cases']}；規則版本：{report['engine_versions']}")
    out.append("- 取代前：判讀者 → GeneBe → pipeline ACMG。取代後（步驟 7 採用）：判讀者 → 新 In-house（平台規則） → GeneBe，"
               "分類、tier 1C 與排序分數都用這個點數。")
    out.append("")

    out.append("## 1. 有效分類的來源（變異次）")
    out.append("")
    out.append("| 檢測 | 變異次 | 取代前來源 | 取代後來源 | 沒有新 In-house 結果 |")
    out.append("|---|---|---|---|---|")
    for tt, s in sorted(report["by_test_type"].items()):
        total = s["occurrences"]
        names = lambda c: {_SOURCE_NAMES.get(k, k): v for k, v in c.items()}  # noqa: E731
        out.append(f"| {tt} | {total} | {_counter_line(names(s['source_before']), total)} | "
                   f"{_counter_line(names(s['source_after']), total)} | {s['engine_missing']} |")
    out.append("")
    out.append("取代後新 In-house 排在 GeneBe 前面，所以有結果的變異都改用新 In-house；判讀者手動判讀仍優先。")
    out.append("")

    out.append("## 2. 分類改變（全部變異次）")
    out.append("")
    for tt, s in sorted(report["by_test_type"].items()):
        total = sum(s["direction"].values())
        out.append(f"### {tt}（{total} 次）")
        out.append("")
        directions = {_DIRECTION.get(k, k): v for k, v in s["direction"].items()}
        out.append(f"- 方向：{_counter_line(directions, total)}")
        out.append(f"- 分類轉換（取代前→A，前 15）：{_counter_line(s['class_transition'], None, 15)}")
        out.append("")

    out.append("## 3. Tier 改變（1A／1B 不受影響）")
    out.append("")
    for tt, s in sorted(report["by_test_type"].items()):
        changed = {k: v for k, v in s["tier_transition"].items() if k.split("→")[0] != k.split("→")[1]}
        out.append(f"- {tt}：{_counter_line(changed) if changed else '沒有改變'}（全部 {s['occurrences']} 次）")
    out.append("")
    out.append("改變 tier 的變異列在 `tier_changes.tsv`（只有變異座標與基因，沒有檢體編號）。")
    out.append("")

    out.append("## 4. 判讀者已標記的變異（1／2／C）")
    out.append("")
    for tt, m in sorted(report["marked"].items()):
        out.append(f"### {tt}（{m['n']} 個）")
        out.append("")
        by_status = {_STATUS.get(k, k): v for k, v in m["by_status"].items()}
        changed = {_STATUS.get(k, k): v for k, v in m["class_changed"].items()}
        out.append(f"- 狀態：{_counter_line(by_status)}")
        out.append(f"- 顯示分類改變：{_counter_line(changed) if changed else '沒有'}；"
                   f"其中會出現在診斷報告（1／2）者 **{m['report_class_changed']}** 個")
        if m["class_transition"]:
            out.append(f"- 分類轉換：{_counter_line(m['class_transition'])}")
        tier_changed = {k: v for k, v in m["tier_transition"].items() if k.split("→")[0] != k.split("→")[1]}
        out.append(f"- Tier 改變：{_counter_line(tier_changed) if tier_changed else '沒有'}")
        out.append("")
        out.append("| 情境 | n | 中位名次 | 第 1 名 | 前 5 | 前 10 | 前 20 | MRR | 名次較好／相同／較差 |")
        out.append("|---|---|---|---|---|---|---|---|---|")
        for s in ("before", "A"):
            r = m["ranks"][s]
            if not r.get("n"):
                continue
            move = m["rank_move"].get(s)
            move_text = (f"{move.get('better', 0)}／{move.get('same', 0)}／{move.get('worse', 0)}"
                         if move else "—")
            out.append(f"| {_SCENARIOS[s]} | {r['n']} | {r['median']} | {_pct(r['top1'], 1)} | {_pct(r['top5'], 1)} | "
                       f"{_pct(r['top10'], 1)} | {_pct(r['top20'], 1)} | {r['mrr']} | {move_text} |")
        out.append("")
    out.append("名次是該個案 review TSV 全部變異依 ACMG 分數＋表型分數排序（同平台公式，不分 tier）。"
               "逐筆明細在 `marked_variants.tsv`，**含檢體編號，只供在伺服器上核對，請勿複製或貼出**。")
    out.append("")

    out.append("## 5. 每個個案前 20 名和取代前相同的比例（有標記變異的個案）")
    out.append("")
    for tt, d in sorted(report["top20_overlap"].items()):
        out.append(f"- {tt}（{d['cases']} 例）：{_pct(d['A'] or 0, 1)}")
    out.append("")

    out.append("## 限制")
    out.append("")
    out.append("- 用三級後處理存好的規則引擎結果；規則或參考表版本已變更（過期）的個案略過，見個案統計的 skipped_stale。")
    out.append("- tier 依平台規則模擬（1A ClinVar、1B LOFTEE HC、1C 點數 ≥4 或預測工具），未含 Somatic 與 specialized caller 點位。")
    out.append("- 排序點數未扣除 PP4／BP5／隱性成對的 PM3；新 In-house 不給這些條目，只有判讀者手動判讀可能含有。")
    out.append("- 只涵蓋 review TSV 內的變異。")
    if report.get("pheno_unavailable_cases"):
        out.append(f"- {report['pheno_unavailable_cases']} 個個案讀不到表型分數，名次只依 ACMG 分數。")
    out.append("")
    return "\n".join(out)
