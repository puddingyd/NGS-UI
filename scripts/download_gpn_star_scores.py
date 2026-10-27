#!/usr/bin/env python3
"""Install the pinned canonical GPN-Star-M LLR bundle (~55 GiB).

Only --download transfers score shards. Default is a manifest-only plan.
Requires curl and Polars; verifies SHA256, schema and a sentinel per shard.
An installation is activated by manifest.json only after all 24 shards pass.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))
from app.services import gpn_star

REPOSITORY = "songlab/gpn-star-scores"
REVISION = "a7b13bbf0d2338d74a7e5f0f8466e41ac0722f50"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def select_files(release: dict) -> list[dict]:
    if release["repository"]["id"] != REPOSITORY:
        raise ValueError("unexpected release repository")
    selected = [f for f in release["parquet"]["files"]
                if f["score_set"] == gpn_star.MODEL and f["score_type"] == "llr"]
    if len(selected) != 24 or {f["chrom"] for f in selected} != gpn_star.PRIMARY_CONTIGS:
        raise ValueError("release must provide exactly 24 human LLR shards")
    for item in selected:
        expected = f"data/{gpn_star.MODEL}/llr/llr_chr{item['chrom']}.parquet"
        if (item["path"] != expected or item["assembly"] != "hg38"
                or not re.fullmatch(r"[0-9a-f]{64}", item["sha256"])
                or item["size"] <= 0 or item["rows"] <= 0):
            raise ValueError(f"invalid release entry: {item['path']}")
    return selected


def verify_shard(path: Path, entry: dict) -> None:
    import polars as pl

    if path.stat().st_size != entry["size"] or sha256(path) != entry["sha256"]:
        raise ValueError(f"size/SHA256 mismatch: {path}; remove the corrupt file before retrying")
    expected = {"chrom": pl.String, "pos": pl.Int64, "ref": pl.String,
                "alt": pl.String, "llr_calibrated": pl.Float32}
    schema = pl.read_parquet_schema(path)
    if any(schema.get(k) != v for k, v in expected.items()):
        raise ValueError(f"unexpected Parquet schema: {path}")
    rows = pl.scan_parquet(path).head(3).collect().to_dicts()
    if not rows or any(row["chrom"] != entry["chrom"] or row["pos"] < 1
                       or row["ref"] == row["alt"]
                       or row["ref"] not in {"A", "C", "G", "T"}
                       or row["alt"] not in {"A", "C", "G", "T"}
                       or row["llr_calibrated"] is None
                       or not math.isfinite(row["llr_calibrated"]) for row in rows):
        raise ValueError(f"invalid Parquet sentinel: {path}")


def install(root: Path, release: dict, revision: str, base_url: str) -> None:
    import polars  # noqa: F401 — fail before starting any large transfer

    files = select_files(release)
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".download.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        marker = root / "manifest.json"
        if marker.exists():
            current = gpn_star.validate_database(root)
            if current["revision"] != revision:
                raise ValueError("active bundle has another revision; install to a new --dest and switch NGS_UI_GPN_STAR_DIR")
        remaining = 0
        for entry in files:
            dest = gpn_star.shard_path(root, entry["chrom"])
            part = root / ".downloads" / revision / f"{dest.name}.part"
            existing = dest if dest.exists() else part
            size = existing.stat().st_size if existing.exists() else 0
            remaining += max(0, entry["size"] - size)
        if shutil.disk_usage(root).free < remaining + 1024 ** 3:
            raise ValueError(f"insufficient disk space: need {remaining / 2**30:.1f} GiB plus 1 GiB reserve")
        local_files = []
        for entry in files:
            dest = gpn_star.shard_path(root, entry["chrom"])
            if dest.exists():
                verify_shard(dest, entry)
            else:
                part = root / ".downloads" / revision / f"{dest.name}.part"
                part.parent.mkdir(parents=True, exist_ok=True)
                if not part.exists() or part.stat().st_size < entry["size"]:
                    print(f"Downloading {entry['path']}", flush=True)
                    subprocess.run(["curl", "--fail", "--location", "--retry", "4",
                                    "--continue-at", "-", "--output", str(part),
                                    f"{base_url}/{entry['path']}"], check=True)
                verify_shard(part, entry)
                dest.parent.mkdir(parents=True, exist_ok=True)
                os.replace(part, dest)
            local_files.append({**entry, "path": dest.relative_to(root).as_posix()})
        payload = {"schema_version": 1, "repository": REPOSITORY,
                   "revision": revision, "model": gpn_star.MODEL,
                   "assembly": "hg38", "score_field": gpn_star.SCORE_FIELD,
                   "license": release["repository"].get("license"),
                   "files": local_files, "source_release": release}
        temporary = root / ".manifest.json.tmp"
        temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary, marker)
        print(f"Verified and activated: {root}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dest", type=Path, default=gpn_star.GPN_STAR_DIR)
    parser.add_argument("--revision", default=REVISION, help="immutable 40-character Hugging Face commit")
    parser.add_argument("--download", action="store_true", help="download, verify and activate all 24 shards")
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.revision):
        parser.error("--revision must be an immutable commit SHA")
    base_url = f"https://huggingface.co/datasets/{REPOSITORY}/resolve/{args.revision}"
    try:
        with tempfile.TemporaryDirectory(prefix="gpn-star-manifest-") as folder:
            path = Path(folder) / "release.json"
            subprocess.run(["curl", "--fail", "--location", "--retry", "4",
                            "--max-time", "60", "--output", str(path),
                            f"{base_url}/manifest/release.json"], check=True)
            release = json.loads(path.read_text(encoding="utf-8"))
        files = select_files(release)
        print(f"{gpn_star.MODEL}: {len(files)} shards, {sum(f['size'] for f in files)/2**30:.2f} GiB")
        print(f"Revision: {args.revision}\nDestination: {args.dest}")
        if args.download:
            install(args.dest, release, args.revision, base_url)
        else:
            print("Plan only. Use --download to transfer and verify the score bundle.")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
