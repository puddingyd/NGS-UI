#!/usr/bin/env python3
"""Show why the In-house rule gave one variant its criteria (read-only).

Prints the stored result of one case (whether it is current or stale), every
criterion with its status and reason, the ClinVar / ClinGen VCEP rule, the
frequency columns the rule read, and — with ``--recompute`` — the result the
current code and reference tables would give now::

    .venv/bin/python scripts/acmg_engine/show_variant.py 26WE0142-nckuh chrX-50090678-C-T
    .venv/bin/python scripts/acmg_engine/show_variant.py 26WE0142-nckuh chrX-50090678-C-T --recompute

The variant may be given as chrN-pos-ref-alt (with or without "chr").
Prints the sample ID you typed; keep the output on the server.
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

FREQ_COLUMNS = ("GNOMAD_E_AF", "GNOMAD_E_EAS_AF", "GNOMAD_E_AF_DBNSFP", "GNOMAD_E_EAS_AF_DBNSFP", "GNOMAD_G_AF",
                "GNOMAD41_JOINT_AF", "GNOMAD41_JOINT_EAS_AF", "GNOMAD_G_EAS_AF", "TG_EAS_AF",
                "INHOUSE_AF", "INHOUSE_AC", "INHOUSE_AN")


def _rows(review: Path, wanted: str):
    with open(review, encoding="utf-8", errors="replace", newline="") as handle:
        for vid, rows in runner.group_variant_rows(csv.DictReader(handle, delimiter="\t")).items():
            if normalize_variant_id(vid) == wanted:
                return vid, rows
    return None, None


def _print_result(title: str, result: dict) -> None:
    rule = result.get("clinvar_rule") or {}
    own = runner.own_call(runner.summary(result))
    print(f"== {title}（規則版本 {result.get('version')}）")
    print(f"  排序用分類／點數：{result.get('classification')}（{result.get('points')}）")
    print(f"  規則自算分類／點數：{own[0]}（{own[1]}）")
    print(f"  ClinVar／專家規則：{rule.get('action')}（來源 {rule.get('source') or '—'}，"
          f"分類 {rule.get('category') or '—'}，{rule.get('stars') or 0}★）")
    if result.get("vcep"):
        print(f"  套用 ClinGen 規格：{result['vcep'].get('vcep')} {result['vcep'].get('svi')}")
    freq = result.get("frequency") or {}
    print(f"  頻率：規則採用的最大 AF {freq.get('max_af')}（gnomAD {freq.get('gnomad_max')}，"
          f"院內 AF {freq.get('inhouse_af')}，院內 AF 是否採用 {freq.get('inhouse_used')}）")
    for status, mark in (("met", "+"), ("not_met", "-"), ("not_evaluated", "?")):
        for ev in result.get("evidence") or []:
            if ev.get("status") == status:
                strength = f" {ev.get('strength')}（{ev.get('points')}）" if status == "met" else ""
                print(f"    {mark} {ev.get('code')}{strength}：{ev.get('reason') or ''}")
    for flag in result.get("flags") or []:
        print(f"    ! {flag.get('code')}：{flag.get('message')}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("sample")
    ap.add_argument("variant", help="chrN-pos-ref-alt")
    ap.add_argument("--recompute", action="store_true", help="also evaluate with the current code and tables")
    args = ap.parse_args(argv)

    case = runner.case_inputs(args.sample)
    if not case["review"].is_file():
        print(f"找不到 review TSV：{case['review']}", file=sys.stderr)
        return 2
    wanted = normalize_variant_id(args.variant if args.variant.startswith("chr") else f"chr{args.variant}")
    vid, rows = _rows(case["review"], wanted)
    if not rows:
        print("這個變異不在 review TSV（主畫面篩選後的清單）裡。", file=sys.stderr)
        return 1
    row = rows[0]
    print(f"{args.sample} {vid} {P.text(row.get('GENE'))} {P.text(row.get('HGVS_C')).split(':')[-1]} "
          f"{P.text(row.get('CONSEQUENCE'))}（{case['test_type']}）")
    print("  review TSV 頻率欄：" + "，".join(f"{c}={P.text(row.get(c)) or '.'}" for c in FREQ_COLUMNS if c in row))
    refs = runner.shared_references()
    existing = Path(case["existing"])
    stored = runner.lookup_many(existing, [vid], detail=True).get(vid) if existing.is_file() else None
    if stored:
        status = "目前有效" if runner.is_current(existing, case["review"], refs) else "已過期（需重算）"
        _print_result(f"已存結果，{status}", stored)
    else:
        print("== 沒有已存結果（個案尚未計算或結果檔不存在）")
    if args.recompute:
        recurrence = runner.shared_recurrence() if case["test_type"] == "WES" else None
        case_refs = runner.CaseRecurrence(refs, recurrence, args.sample,
                                          runner.case_patient(recurrence, args.sample, case["mrn"]))
        for _, result in runner.evaluate_rows({vid: rows}, test_type=case["test_type"], refs=case_refs):
            _print_result("目前程式重算", result)
    return 0


if __name__ == "__main__":
    sys.exit(main())
