#!/usr/bin/env python3
"""Build the ClinVar case table for PS4 (affected cases) and PM3 (biallelic patients), rule engine v1.6–1.7.

Downloads the ClinVar VCV XML release (about 5 GB; unchanged files are
skipped) to ``<NGS_UI_ACMG_REFS_DIR>/raw/`` and writes
``<NGS_UI_ACMG_REFS_DIR>/clinvar_cases.sqlite``: per Variation ID, how many
submitters reported the variant in affected individuals (clinical testing,
research, provider interpretation, case-control; literature-only, curation,
benign submissions and somatic samples are left out), and which submissions
describe homozygous / compound-heterozygous patients (structured counts or
the comment text, with quoted snippets for PM3).  Parsing takes roughly
30–60 minutes; the old table stays in use until the new one is complete.

    .venv/bin/python scripts/acmg_engine/build_clinvar_cases.py
    .venv/bin/python scripts/acmg_engine/build_clinvar_cases.py --xml /path/ClinVarVCVRelease_00-latest.xml.gz

Afterwards recompute the stored engine results:
    .venv/bin/python scripts/acmg_engine/build_engine_results.py --all --force
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from app import config  # noqa: E402
from app.services.acmg_engine import clinvar_cases as CC  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    refs_dir = Path(config.ACMG_REFS_DIR)
    ap.add_argument("--xml", type=Path, default=None, help="use this local XML (.xml or .xml.gz), no download")
    ap.add_argument("--url", default=CC.XML_URL)
    ap.add_argument("--out", type=Path, default=refs_dir / CC.CASES_DB_NAME)
    ap.add_argument("--force", action="store_true", help="download even if unchanged")
    args = ap.parse_args(argv)
    release = ""
    if args.xml is None:
        from download_references import fetch
        raw = refs_dir / "raw"
        raw.mkdir(parents=True, exist_ok=True)
        args.xml = raw / CC.XML_NAME
        info = fetch(args.url, args.xml, args.force)
        release = info.get("last_modified", "")
    else:
        sidecar = args.xml.with_name(args.xml.name + ".json")
        if sidecar.is_file():
            release = json.loads(sidecar.read_text(encoding="utf-8")).get("last_modified", "")
    if not args.xml.is_file():
        print(f"找不到 {args.xml}", file=sys.stderr)
        return 2
    print(f"解析 {args.xml}（約 30–60 分鐘）…", flush=True)
    last = [0]

    def progress(msg):
        last[0] += 1
        if last[0] % 20 == 0:
            print(msg, flush=True)

    stats = CC.build(args.xml, args.out, release=release, progress=progress)
    print("\n" + "=" * 60)
    print(f"完成：{stats['variants']:,} 個 ClinVar 變異有病例資料（其中 {stats['pm3_variants']:,} 個有同型合子／"
          f"複合異型合子病人的描述，供 PM3；{stats.get('ps3_variants', 0):,} 個有功能研究異常的描述，供 PS3），"
          f"耗時 {stats['seconds'] / 60:.0f} 分鐘")
    print(f"寫入 {stats['path']}")
    print("接著重算 In-house 結果：.venv/bin/python scripts/acmg_engine/build_engine_results.py --all --force")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
