#!/usr/bin/env python3
"""Download public gene-disease sources and rebuild the SQLite index."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skip-download", action="store_true")
    ap.add_argument("--raw-dir", type=Path)
    ap.add_argument("--db", type=Path)
    args = ap.parse_args()

    py = sys.executable
    if not args.skip_download:
        cmd = [py, str(REPO_ROOT / "scripts" / "download_gene_disease_sources.py")]
        if args.raw_dir:
            cmd.extend(["--raw-dir", str(args.raw_dir)])
        subprocess.check_call(cmd, cwd=str(REPO_ROOT))

    cmd = [py, str(REPO_ROOT / "scripts" / "build_gene_disease_index.py")]
    if args.raw_dir:
        cmd.extend(["--raw-dir", str(args.raw_dir)])
    if args.db:
        cmd.extend(["--db", str(args.db)])
    subprocess.check_call(cmd, cwd=str(REPO_ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
