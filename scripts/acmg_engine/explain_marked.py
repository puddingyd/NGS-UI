#!/usr/bin/env python3
"""Why did the In-house rule give these reviewer-marked variants a different class?

Reads ``marked_variants.tsv`` from diff_replacement.py and, for the selected
rows, prints the rule's met criteria (with points and reasons), the criteria it
could not evaluate, and GeneBe's criteria and ClinVar, so the two can be
compared side by side::

    .venv/bin/python scripts/acmg_engine/explain_marked.py ~/NGS_UI/acmg_validation/diff2/marked_variants.tsv
    .venv/bin/python scripts/acmg_engine/explain_marked.py TSV --status 1 --to VUS     # default
    .venv/bin/python scripts/acmg_engine/explain_marked.py TSV --status 1,2 --to any
    .venv/bin/python scripts/acmg_engine/explain_marked.py TSV --tsv ~/NGS_UI/acmg_validation/causative_vus.tsv

Read-only.  Output has no sample IDs (the same variant in several cases is
shown once with a case count) unless ``--with-sample`` is given; output with
sample IDs is for checking on the server only.  ``--tsv`` writes a record with
one row per case *including sample IDs* (class at the time of that diff run,
the case's current In-house class and criteria, blank columns for the
reviewer); keep it on the server.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))

from app.services.acmg_engine import parse as P  # noqa: E402
from app.services.acmg_engine import runner  # noqa: E402
from app.services.acmg_engine.validation import normalize_variant_id  # noqa: E402


def _selected(path: Path, statuses: set[str], to: str):
    with open(path, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("status") not in statuses:
                continue
            before, after = row.get("class_before") or "", row.get("class_after") or ""
            if before == after:
                continue
            if to != "any" and after != to:
                continue
            yield row


def _review_rows(review: Path, wanted: set[str]) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    with open(review, encoding="utf-8", errors="replace", newline="") as handle:
        for vid, rows in runner.group_variant_rows(csv.DictReader(handle, delimiter="\t")).items():
            if normalize_variant_id(vid) in wanted:
                out[normalize_variant_id(vid)] = rows
    return out


def _first(rows, column):
    for row in rows or []:
        value = P.text(row.get(column))
        if value:
            return value
    return ""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("marked_tsv", type=Path)
    ap.add_argument("--status", default="1", help="comma-separated reviewer statuses (default 1)")
    ap.add_argument("--to", default="VUS", help="class after replacement (P/LP/VUS/LB/B) or 'any' (default VUS)")
    ap.add_argument("--with-sample", action="store_true", help="also print sample IDs (server only)")
    ap.add_argument("--tsv", type=Path, default=None,
                    help="write a per-case record WITH sample IDs to this file (server only)")
    args = ap.parse_args(argv)
    records: list[dict] = []

    by_sample: dict[str, list[dict]] = {}
    for row in _selected(args.marked_tsv, {s.strip() for s in args.status.split(",") if s.strip()}, args.to):
        by_sample.setdefault(row["sample"], []).append(row)

    merged: dict[str, dict] = {}
    for sid, rows in by_sample.items():
        case = runner.case_inputs(sid)
        wanted = {normalize_variant_id(r["variant"]) for r in rows}
        review = _review_rows(case["review"], wanted) if case["review"].is_file() else {}
        ids = [runner.ui_variant_id(group[0]) for group in review.values()]
        stored = runner.lookup_many(case["existing"], ids, detail=True)
        stored = {normalize_variant_id(k): v for k, v in stored.items()}
        for row in rows:
            key = normalize_variant_id(row["variant"])
            item = merged.setdefault(f"{row['test_type']}|{key}", {
                "row": row, "result": stored.get(key), "tsv": review.get(key), "samples": []})
            item["samples"].append(sid)
            result = stored.get(key) or {}
            tsv = review.get(key)
            records.append({
                "sample": sid, "test_type": row["test_type"], "status": row["status"],
                "gene": row["gene"], "hgvs_c": row["hgvs_c"], "hgvs_p": row["hgvs_p"], "variant": row["variant"],
                "consequence": _first(tsv, "CONSEQUENCE"),
                "class_before": f"{row['class_before']}（{row['points_before']}，{row['source_before']}）",
                "class_in_diff_run": f"{row['class_after']}（{row['points_after']}）",
                "inhouse_now": (f"{result.get('classification', '')}（{result.get('points', '')}）" if result else "找不到結果"),
                "inhouse_version": result.get("version", ""),
                "inhouse_criteria": result.get("criteria_text", ""),
                "not_evaluated": "；".join(f"{e['code']}：{e['reason']}" for e in result.get("not_evaluated") or []),
                "clinvar": (_first(tsv, "CLINVAR_SIG") or "—") + (f"（{_first(tsv, 'CLINVAR_STARS')}★）"
                                                                  if _first(tsv, "CLINVAR_STARS") else ""),
                "genebe_criteria": _first(tsv, "GENEBE_ACMG_CRITERIA"),
                "manual_acmg_done": "", "case_evidence_added": "", "reviewer_note": "",
            })

    if args.tsv:
        args.tsv.parent.mkdir(parents=True, exist_ok=True)
        with open(args.tsv, "w", encoding="utf-8-sig", newline="") as handle:
            fields = list(records[0]) if records else ["sample"]
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for rec in sorted(records, key=lambda r: (r["test_type"], r["gene"], r["sample"])):
                writer.writerow(rec)
        print(f"已寫出 {len(records)} 列（含檢體編號，只留在伺服器上）：{args.tsv}")
    if not merged:
        print("沒有符合條件的點位。")
        return 0
    for item in merged.values():
        row, result, tsv = item["row"], item["result"] or {}, item["tsv"]
        head = " ".join(x for x in (row["test_type"], row["gene"], row["hgvs_c"], row["hgvs_p"]) if x)
        print("=" * 78)
        print(f"{head}  （{len(item['samples'])} 個個案）")
        if args.with_sample:
            print(f"  檢體：{', '.join(item['samples'])}")
        print(f"  {row['variant']}  {_first(tsv, 'CONSEQUENCE')}")
        print(f"  取代前 {row['class_before']}（{row['points_before']}，{row['source_before']}）"
              f" → In-house {row['class_after']}（{row['points_after']}）")
        stars = _first(tsv, "CLINVAR_STARS")
        print(f"  ClinVar：{_first(tsv, 'CLINVAR_SIG') or '—'}{f'（{stars}★）' if stars else ''}"
              f"；最新版 {_first(tsv, 'CLINVAR_LATEST_SIG') or '—'}")
        print(f"  GeneBe criteria：{_first(tsv, 'GENEBE_ACMG_CRITERIA') or '—'}")
        if not result:
            print("  （找不到 In-house 結果）")
            continue
        freq = result.get("frequency") or {}
        print(f"  頻率：max AF {freq.get('max_af')}，gnomAD {freq.get('gnomad_max')}，院內 AF {freq.get('inhouse_af')}")
        print(f"  In-house criteria：{result.get('criteria_text') or '—'}"
              f"（規則點數 {result.get('engine_points')}，ClinVar 規則 {(result.get('clinvar_rule') or {}).get('action', 'none')}）")
        for ev in result.get("evidence") or []:
            if ev.get("status") == "met":
                print(f"    + {ev['code']} {ev.get('strength')}（{ev.get('points')}）：{ev.get('reason')}")
        for ev in result.get("evidence") or []:
            if ev.get("status") == "not_met" and ev.get("reason"):
                print(f"    - {ev['code']} 未符合：{ev['reason']}")
        for ev in result.get("not_evaluated") or []:
            print(f"    ? {ev['code']} 未評估：{ev['reason']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
