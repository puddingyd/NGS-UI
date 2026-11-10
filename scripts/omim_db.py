#!/usr/bin/env python3
"""Manage the platform OMIM database (NGS_UI_HOME/OMIM/omim.sqlite).

    PYTHONPATH=backend python scripts/omim_db.py init            # build from OMIM.xlsx
    PYTHONPATH=backend python scripts/omim_db.py status
    PYTHONPATH=backend python scripts/omim_db.py reviews [--kind KIND]
    PYTHONPATH=backend python scripts/omim_db.py export --out OMIM_export.xlsx

``init`` first freezes legacy SNV report disease ticks
(scripts/migrate_report_disease_picks.py --apply) while the platform still
reads OMIM.xlsx, then builds the database and installs it atomically. From
then on omim_store reads the database; OMIM.xlsx is no longer read (use
``export`` to produce a workbook). ``init --force`` rebuilds from the
workbook and DISCARDS adopted captures; the previous database is kept as a
timestamped .bak file next to it.

Remove (or rename) the database file to fall back to OMIM.xlsx.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app.config import OMIM_XLSX  # noqa: E402
from app.services import omim_db  # noqa: E402


def _migrate_picks() -> None:
    spec = importlib.util.spec_from_file_location(
        "migrate_report_disease_picks", REPO_ROOT / "scripts" / "migrate_report_disease_picks.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.main(["--apply"])


def cmd_init(args) -> int:
    target = omim_db.db_path()
    xlsx = Path(args.xlsx) if args.xlsx else OMIM_XLSX
    if omim_db.exists() and not args.force:
        print(f"{target} already exists. Use --force to rebuild from {xlsx} "
              "(adopted captures in the current database will be discarded; a .bak is kept).")
        return 1
    if not xlsx.is_file():
        print(f"workbook not found: {xlsx}")
        return 1
    if not args.skip_pick_migration:
        print("== freezing legacy SNV report disease ticks (migrate_report_disease_picks --apply)")
        _migrate_picks()
    print(f"== building OMIM database from {xlsx}")
    built = omim_db.build_from_xlsx(xlsx, target)
    backup = omim_db.install(Path(built["path"]), target)
    print(json.dumps(built["stats"], ensure_ascii=False, indent=2))
    if backup:
        print(f"previous database backed up to {backup}")
    print(f"installed {target}; restart is not required (omim_store reloads on change).")
    return 0


def cmd_status(_args) -> int:
    print(json.dumps(omim_db.status(), ensure_ascii=False, indent=2))
    return 0


def cmd_reviews(args) -> int:
    if not omim_db.exists():
        print("database not built yet")
        return 1
    data = omim_db.list_reviews("open", args.kind or "", args.limit)
    print(json.dumps(data["open_counts"], ensure_ascii=False))
    for item in data["items"]:
        detail = item["detail"]
        print(f"[{item['id']}] {item['kind']}\t{item.get('gene_symbol') or ''}\t"
              f"{detail.get('omim_disease') or detail.get('label') or detail.get('disease_cell') or ''}")
    return 0


def cmd_export(args) -> int:
    if not omim_db.exists():
        print("database not built yet")
        return 1
    Path(args.out).write_bytes(omim_db.export_xlsx())
    print(f"wrote {args.out}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Platform OMIM database")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("init", help="build the database from OMIM.xlsx")
    p.add_argument("--xlsx", help="workbook path (default: NGS_UI_OMIM_XLSX)")
    p.add_argument("--force", action="store_true", help="rebuild even if the database exists")
    p.add_argument("--skip-pick-migration", action="store_true",
                   help="do not run migrate_report_disease_picks --apply first")
    p.set_defaults(func=cmd_init)
    p = sub.add_parser("status")
    p.set_defaults(func=cmd_status)
    p = sub.add_parser("reviews")
    p.add_argument("--kind", default="")
    p.add_argument("--limit", type=int, default=500)
    p.set_defaults(func=cmd_reviews)
    p = sub.add_parser("export")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_export)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
