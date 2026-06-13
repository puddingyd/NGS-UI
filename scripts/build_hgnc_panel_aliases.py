#!/usr/bin/env python3
"""Build panel-safe HGNC alias mappings from official HGNC downloads.

Inputs:
  reference/hgnc/hgnc_complete_set.txt
  reference/hgnc/withdrawn.txt
  reference/hgnc/manual_panel_aliases.tsv

Output:
  ngs_panel_deadzone/panel/panel_gene_aliases.tsv

Rules are conservative for panel gene lists:
  - previous symbols map to their current approved symbol
  - alias symbols map only when the alias has one unique target
  - withdrawn symbols map only when merged/split into one approved symbol
  - aliases that are themselves current approved symbols are not auto-mapped
  - manual aliases are kept as curated overrides when their target is current
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def _split_pipe(raw: str) -> list[str]:
    return [part.strip() for part in (raw or "").split("|") if part.strip()]


def _read_current_symbols(path: Path) -> tuple[set[str], dict[str, str], dict[str, set[str]], dict[str, set[str]]]:
    current: set[str] = set()
    by_prev: dict[str, set[str]] = defaultdict(set)
    by_alias: dict[str, set[str]] = defaultdict(set)
    symbol_to_hgnc: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if (row.get("status") or "").strip().lower() != "approved":
                continue
            symbol = (row.get("symbol") or "").strip()
            hgnc_id = (row.get("hgnc_id") or "").strip()
            if not symbol:
                continue
            current.add(symbol)
            if hgnc_id:
                symbol_to_hgnc[symbol] = hgnc_id
            for alias in _split_pipe(row.get("prev_symbol") or ""):
                if alias != symbol:
                    by_prev[alias].add(symbol)
            for alias in _split_pipe(row.get("alias_symbol") or ""):
                if alias != symbol:
                    by_alias[alias].add(symbol)
    return current, symbol_to_hgnc, by_prev, by_alias


def _read_withdrawn(path: Path, current: set[str]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    replacement_col = "MERGED_INTO_REPORT(S) (i.e HGNC_ID|SYMBOL|STATUS)"
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            symbol = (row.get("WITHDRAWN_SYMBOL") or "").strip()
            replacements = (row.get(replacement_col) or "").strip()
            if not symbol or not replacements:
                continue
            parts = [part.strip() for part in replacements.split(",") if part.strip()]
            parsed: list[str] = []
            for part in parts:
                fields = part.split("|")
                if len(fields) < 3:
                    continue
                repl_symbol = fields[1].strip()
                status = fields[2].strip().lower()
                if repl_symbol in current and status == "approved":
                    parsed.append(repl_symbol)
            if len(set(parsed)) == 1:
                target = parsed[0]
                if symbol != target:
                    out[symbol].add(target)
    return out


def _read_manual(path: Path) -> list[tuple[str, str, str]]:
    if not path.is_file():
        return []
    rows: list[tuple[str, str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            alias = (row.get("alias") or "").strip()
            target = (row.get("hgnc_symbol") or "").strip()
            source = (row.get("source") or "manual").strip() or "manual"
            if alias and target and alias != target:
                rows.append((alias, target, source))
    return rows


def _add_unique(
    aliases: dict[str, tuple[str, str]],
    conflicts: list[tuple[str, str, str, str]],
    alias: str,
    target: str,
    source: str,
) -> None:
    existing = aliases.get(alias)
    if existing and existing[0] != target:
        conflicts.append((alias, existing[0], target, source))
        return
    aliases[alias] = (target, source if not existing else existing[1])


def build_aliases(
    complete_set: Path,
    withdrawn: Path,
    manual: Path,
) -> tuple[dict[str, tuple[str, str]], list[tuple[str, str, str, str]]]:
    current, _symbol_to_hgnc, by_prev, by_alias = _read_current_symbols(complete_set)
    withdrawn_map = _read_withdrawn(withdrawn, current)
    aliases: dict[str, tuple[str, str]] = {}
    conflicts: list[tuple[str, str, str, str]] = []

    for alias, targets in sorted(by_prev.items()):
        if alias in current:
            conflicts.append((alias, "current_symbol", "|".join(sorted(targets)), "hgnc_prev_symbol"))
            continue
        if len(targets) == 1:
            _add_unique(aliases, conflicts, alias, next(iter(targets)), "hgnc_prev_symbol")
        else:
            conflicts.append((alias, "multiple_targets", "|".join(sorted(targets)), "hgnc_prev_symbol"))

    for alias, targets in sorted(withdrawn_map.items()):
        if alias in current:
            conflicts.append((alias, "current_symbol", "|".join(sorted(targets)), "hgnc_withdrawn"))
            continue
        if len(targets) == 1:
            _add_unique(aliases, conflicts, alias, next(iter(targets)), "hgnc_withdrawn")
        else:
            conflicts.append((alias, "multiple_targets", "|".join(sorted(targets)), "hgnc_withdrawn"))

    for alias, targets in sorted(by_alias.items()):
        if alias in current:
            conflicts.append((alias, "current_symbol", "|".join(sorted(targets)), "hgnc_alias_symbol"))
            continue
        if len(targets) == 1:
            _add_unique(aliases, conflicts, alias, next(iter(targets)), "hgnc_alias_symbol")
        else:
            conflicts.append((alias, "multiple_targets", "|".join(sorted(targets)), "hgnc_alias_symbol"))

    for alias, target, source in _read_manual(manual):
        if target not in current:
            conflicts.append((alias, "manual_target_not_current", target, source))
            continue
        aliases[alias] = (target, source)

    return aliases, conflicts


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--complete-set", type=Path, default=Path("reference/hgnc/hgnc_complete_set.txt"))
    parser.add_argument("--withdrawn", type=Path, default=Path("reference/hgnc/withdrawn.txt"))
    parser.add_argument("--manual", type=Path, default=Path("reference/hgnc/manual_panel_aliases.tsv"))
    parser.add_argument("--out", type=Path, default=Path("ngs_panel_deadzone/panel/panel_gene_aliases.tsv"))
    parser.add_argument("--conflicts", type=Path, default=Path("docs/ops/hgnc_alias_conflicts.tsv"))
    args = parser.parse_args()

    aliases, conflicts = build_aliases(args.complete_set, args.withdrawn, args.manual)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["alias", "hgnc_symbol", "source"])
        for alias, (target, source) in sorted(aliases.items(), key=lambda item: item[0].upper()):
            writer.writerow([alias, target, source])

    args.conflicts.parent.mkdir(parents=True, exist_ok=True)
    with args.conflicts.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["alias", "existing_or_reason", "candidate_target", "source"])
        for row in sorted(conflicts, key=lambda item: (item[0].upper(), item[3])):
            writer.writerow(row)

    print(f"aliases\t{len(aliases)}")
    print(f"conflicts\t{len(conflicts)}")
    print(f"out\t{args.out}")
    print(f"conflicts_out\t{args.conflicts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
