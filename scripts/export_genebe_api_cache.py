#!/usr/bin/env python3
"""Export all successful GeneBe live-API cache rows as one DB-compatible TSV."""
from __future__ import annotations

import argparse
import csv
import json
import os
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from annotate_acmg_genebe import (
    DB_EXPORT_FIELDS,
    _score_for_export,
    _variant_sort_key,
    db_classification,
)


def export_cache(cache: Path, output: Path) -> int:
    if not cache.is_file():
        raise FileNotFoundError(cache)
    rows: dict[str, tuple[str, str, str, str, str, str, str]] = {}
    with sqlite3.connect(cache) as conn:
        for chrom, pos, ref, alt, cls, score, criteria in conn.execute(
            "SELECT chrom, pos, ref, alt, acmg_classification, "
            "acmg_score, acmg_criteria FROM results WHERE status='success'"
        ):
            key = f"{chrom}:{pos}:{ref}:{alt}"
            rows[key] = (
                str(chrom), str(pos), str(ref).upper(), str(alt).upper(),
                db_classification(str(cls)) or ".",
                _score_for_export(str(score)),
                str(criteria or ".").strip() or ".",
            )

    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_suffix(output.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(DB_EXPORT_FIELDS),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for key in sorted(rows, key=_variant_sort_key):
            chrom, pos, ref, alt, cls, score, criteria = rows[key]
            writer.writerow({
                "#chr": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "acmg_classification": cls,
                "acmg_score": score,
                "acmg_criteria": criteria,
            })
    os.replace(tmp, output)

    sidecar = output.with_suffix(output.suffix + ".json")
    sidecar_tmp = sidecar.with_suffix(sidecar.suffix + ".tmp")
    sidecar_tmp.write_text(
        json.dumps({
            "schema": "genebe-api-cache-export-v1",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "cache": str(cache.resolve()),
            "rows": len(rows),
            "columns": list(DB_EXPORT_FIELDS),
        }, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    os.replace(sidecar_tmp, sidecar)
    return len(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cache",
        type=Path,
        default=Path(os.environ.get(
            "NGS_UI_GENEBE_API_CACHE",
            Path.home() / "NGS_UI" / "biotools" / "genebe"
            / "genebe_api_cache.sqlite",
        )),
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    try:
        count = export_cache(args.cache.resolve(), args.out.resolve())
    except (FileNotFoundError, sqlite3.Error, ValueError) as exc:
        print(f"ERROR: {exc}")
        return 2
    print(f"[genebe-api-export] rows={count} → {args.out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
