#!/usr/bin/env python3
"""Download public gene-disease resources used by the NGS-UI disease index."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tempfile
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app.config import GENE_DISEASE_RAW_DIR  # noqa: E402

SOURCES = {
    "gencc": {
        "filename": "gencc_submissions.csv",
        "urls": [
            "https://search.thegencc.org/download/action/submissions-export-csv",
            "https://thegencc.org/download/action/submissions-export-csv",
        ],
    },
    "clingen": {
        "filename": "clingen_gene_validity.csv",
        "urls": [
            "https://search.clinicalgenome.org/kb/gene-validity/download",
        ],
    },
    "mondo": {
        "filename": "mondo.json",
        "urls": [
            "https://purl.obolibrary.org/obo/mondo.json",
            "https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.json",
        ],
    },
}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _download_one(name: str, raw_dir: Path, *, timeout: int) -> dict:
    spec = SOURCES[name]
    target = raw_dir / spec["filename"]
    errors = []
    for url in spec["urls"]:
        tmp = None
        started = time.time()
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "NGS-UI gene-disease updater"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                suffix = target.suffix or ".tmp"
                with tempfile.NamedTemporaryFile(
                    "wb", delete=False, suffix=suffix, dir=str(raw_dir)
                ) as out:
                    tmp = Path(out.name)
                    while True:
                        chunk = resp.read(1024 * 1024)
                        if not chunk:
                            break
                        out.write(chunk)
                if tmp.stat().st_size == 0:
                    raise RuntimeError("downloaded empty file")
                tmp.replace(target)
                return {
                    "source": name,
                    "url": resp.geturl(),
                    "requested_url": url,
                    "path": str(target),
                    "size": target.stat().st_size,
                    "sha256": _sha256(target),
                    "downloaded_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    "elapsed_seconds": round(time.time() - started, 2),
                }
        except Exception as exc:
            errors.append(f"{url}: {exc}")
            if tmp and tmp.exists():
                tmp.unlink(missing_ok=True)
    raise RuntimeError(f"failed to download {name}: " + " | ".join(errors))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-dir", type=Path, default=GENE_DISEASE_RAW_DIR)
    ap.add_argument("--only", choices=sorted(SOURCES), action="append")
    ap.add_argument("--timeout", type=int, default=240)
    args = ap.parse_args()

    args.raw_dir.mkdir(parents=True, exist_ok=True)
    selected = args.only or list(SOURCES)
    manifest = []
    for name in selected:
        print(f"[download] {name}", flush=True)
        manifest.append(_download_one(name, args.raw_dir, timeout=args.timeout))
    out_path = args.raw_dir / "download_manifest.json"
    out_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"[download] wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
