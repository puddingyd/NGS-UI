#!/usr/bin/env python3
"""Download the ACMG rule engine's public sources (plan §5.2.11).

* ClinVar ``variant_summary.txt.gz`` (GRCh38 P/LP truncating / missense);
* MANE Ensembl genomic GTF (per-exon CDS lengths).

Files go to ``<NGS_UI_ACMG_REFS_DIR>/raw/`` with a ``.json`` sidecar holding
the source URL and HTTP Last-Modified (used as the ClinVar release date).
Unchanged files (same size and Last-Modified) are skipped unless ``--force``.
Downloads are written to a temporary name and renamed when complete.

    python3 scripts/acmg_engine/download_references.py
    # time-split validation: an archived ClinVar month
    python3 scripts/acmg_engine/download_references.py --skip-mane \
      --clinvar-url https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2025-01.txt.gz
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))

from app import config  # noqa: E402

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
MANE_DIR_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
USER_AGENT = "ngs-ui-acmg-refs"


def _request(url: str, method: str = "GET"):
    return urllib.request.urlopen(urllib.request.Request(url, method=method, headers={"User-Agent": USER_AGENT}),
                                  timeout=60)


def mane_gtf_url() -> str:
    with _request(MANE_DIR_URL) as resp:
        listing = resp.read().decode("utf-8", "replace")
    names = sorted(set(re.findall(r"MANE\.GRCh38\.v([\d.]+)\.ensembl_genomic\.gtf\.gz", listing)),
                   key=lambda v: [int(x) for x in v.split(".") if x.isdigit()])
    if not names:
        raise SystemExit(f"no MANE GTF found in {MANE_DIR_URL}")
    return f"{MANE_DIR_URL}MANE.GRCh38.v{names[-1]}.ensembl_genomic.gtf.gz"


def fetch(url: str, dest: Path, force: bool = False) -> dict:
    with _request(url, "HEAD") as resp:
        size = int(resp.headers.get("Content-Length") or 0)
        modified = resp.headers.get("Last-Modified") or ""
    sidecar = dest.with_name(dest.name + ".json")
    if dest.is_file() and sidecar.is_file() and not force:
        try:
            old = json.loads(sidecar.read_text(encoding="utf-8"))
        except ValueError:
            old = {}
        if old.get("last_modified") == modified and dest.stat().st_size == size:
            print(f"[skip] {dest.name} 已是最新（{modified}）")
            return old
    tmp = dest.with_name(f".{dest.name}.{os.getpid()}.part")
    done = 0
    started = time.time()
    with _request(url) as resp, open(tmp, "wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
            done += len(chunk)
            if size and done % (64 << 20) < (1 << 20):
                print(f"  {dest.name}: {done / size:.0%}", flush=True)
    if size and done != size:
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"{url}: downloaded {done} bytes, expected {size}")
    os.replace(tmp, dest)
    info = {"url": url, "last_modified": modified, "size": done,
            "downloaded_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "seconds": round(time.time() - started)}
    sidecar.write_text(json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[done] {dest}（{done / 1e6:.1f} MB，{modified}）")
    return info


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dest", type=Path, default=Path(config.ACMG_REFS_DIR) / "raw")
    ap.add_argument("--clinvar-url", default=CLINVAR_URL)
    ap.add_argument("--clinvar-name", default="", help="output file name (default: name in the URL)")
    ap.add_argument("--mane-gtf-url", default="", help="default: newest GTF in MANE current/")
    ap.add_argument("--skip-clinvar", action="store_true")
    ap.add_argument("--skip-mane", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    args.dest.mkdir(parents=True, exist_ok=True)
    if not args.skip_clinvar:
        name = args.clinvar_name or args.clinvar_url.rsplit("/", 1)[-1]
        fetch(args.clinvar_url, args.dest / name, args.force)
    if not args.skip_mane:
        url = args.mane_gtf_url or mane_gtf_url()
        fetch(url, args.dest / url.rsplit("/", 1)[-1], args.force)
    return 0


if __name__ == "__main__":
    sys.exit(main())
