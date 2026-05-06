#!/usr/bin/env python3
"""Migrate from flat per-sample layout to per-analysis-version layout.

Before:
    tertiary_output/{sid}/
      sample_metadata.json     (meta + hpo + selected_panels)
      snv_indel.annotated.tsv
      pheno_score.tsv
      exomiser_results.tsv
      lirical_results.tsv
      analysis_files/...
    data/reports/{sid}.json    (status + edits + ... + clinical/comment/tags)

After:
    tertiary_output/{sid}/
      sample_metadata.json     (meta + clinical/comment/tags + status/edits/
                                panels/manual_variants + active_analysis)
      snv_indel.annotated.tsv
      analyses/default/
        analysis.json          (hpo + selected_panels + note)
        pheno_score.tsv
        exomiser_results.tsv
        lirical_results.tsv
        analysis_files/...
    data/reports/{sid}.json.bak  (old report kept for safety)

Idempotent: re-runs detect already-migrated samples (analyses/default/
analysis.json present) and skip them. Safe to run before or after the
new backend code is deployed because the loader keeps a fallback read
path against the old layout.

Usage:
    python3 scripts/migrate_to_versioned_layout.py
    NGS_UI_HOME=/home/n102968/NGS_UI python3 scripts/migrate_to_versioned_layout.py
"""
from __future__ import annotations

import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path


def _home() -> Path:
    env = os.environ.get("NGS_UI_HOME")
    if env:
        return Path(env)
    here = Path(__file__).resolve()
    # repo root is parent of scripts/, NGS_UI_HOME is parent of repo
    repo = here.parent.parent
    parent = repo.parent
    if (parent / "tertiary_output").exists():
        return parent
    if (repo / "tertiary_output").exists():
        return repo
    raise SystemExit(
        "Cannot locate NGS_UI_HOME; set the env var or run from a tree "
        "where tertiary_output/ exists next to or inside the repo."
    )


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _read_json(p: Path) -> dict:
    if not p.exists():
        return {}
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except json.JSONDecodeError:
        return {}


def _write_json(p: Path, data: dict) -> None:
    p.write_text(
        json.dumps(data, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def migrate_sample(tertiary_dir: Path, reports_dir: Path) -> str:
    """Returns 'migrated' / 'skipped' / 'no-meta'."""
    sid = tertiary_dir.name
    meta_path = tertiary_dir / "sample_metadata.json"
    if not meta_path.exists():
        return "no-meta"

    default_dir = tertiary_dir / "analyses" / "default"
    if (default_dir / "analysis.json").exists():
        return "skipped"

    meta = _read_json(meta_path)

    # Pull version-specific fields out of meta.
    hpo = meta.pop("hpo", None) or meta.pop("patient_phenotype", None) or []
    panels = meta.pop("selected_panels", None) or []

    # Pull reviewer-decision fields from the old data/reports/{sid}.json.
    report_path = reports_dir / f"{sid}.json"
    report = _read_json(report_path)

    # Merge reviewer-side state into the patient record.
    meta["clinical_description"] = report.get("clinical_description", "")
    meta["comment"] = report.get("comment", "")
    meta["tags"] = report.get("tags") or []
    meta["status"] = report.get("status") or {}
    meta["edits"] = report.get("edits") or {}
    meta["panels"] = report.get("panels") or {}
    meta["manual_variants"] = report.get("manual_variants") or []
    if not meta.get("category") and report.get("category"):
        meta["category"] = report.get("category")
    meta["active_analysis"] = "default"
    now = _now()
    meta.setdefault("created_at", now)
    meta["updated_at"] = now

    _write_json(meta_path, meta)

    # Write analysis.json for the default version.
    default_dir.mkdir(parents=True, exist_ok=True)
    _write_json(
        default_dir / "analysis.json",
        {
            "hpo": hpo,
            "selected_panels": panels,
            "note": "",
            "created_at": now,
            "updated_at": now,
        },
    )

    # Move sidecars + the analysis_files run directory into the default dir.
    for fname in ("pheno_score.tsv", "exomiser_results.tsv", "lirical_results.tsv"):
        src = tertiary_dir / fname
        if src.exists():
            shutil.move(str(src), str(default_dir / fname))
    af = tertiary_dir / "analysis_files"
    if af.exists():
        shutil.move(str(af), str(default_dir / "analysis_files"))

    # Backup the old report file (don't delete; rerunning the migration
    # safely is more important than reclaiming a few KB).
    if report_path.exists():
        bak = report_path.with_suffix(".json.bak")
        report_path.rename(bak)

    return "migrated"


def main() -> int:
    home = _home()
    tertiary = home / "tertiary_output"
    reports = home / "data" / "reports"
    if not tertiary.is_dir():
        print(f"!! {tertiary} not found", file=sys.stderr)
        return 1

    counts = {"migrated": 0, "skipped": 0, "no-meta": 0}
    for sub in sorted(tertiary.iterdir()):
        if not sub.is_dir() or sub.name.startswith("_"):
            continue
        result = migrate_sample(sub, reports)
        counts[result] += 1
        print(f"  {result:8} {sub.name}")
    print(f"\nDone: {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
