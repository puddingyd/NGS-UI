#!/usr/bin/env python3
"""Atomically download and rebuild the weekly ClinVar GRCh38 UI database."""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import shutil
import sqlite3
import sys
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app import config  # noqa: E402
from app.services import clinvar_latest_store  # noqa: E402


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "NGS-UI ClinVar updater/1"})
    with urllib.request.urlopen(request, timeout=180) as response, destination.open("wb") as target:
        shutil.copyfileobj(response, target, length=1024 * 1024)


def update(
    *,
    vcf_file: Path | None,
    url: str,
    db_path: Path,
    vcf_path: Path,
    manifest_path: Path,
    release_date: str = "",
) -> dict[str, object]:
    output_dir = db_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="clinvar-update-", dir=output_dir) as temp_name:
        temp_dir = Path(temp_name)
        candidate_vcf = temp_dir / "clinvar.vcf.gz"
        if vcf_file is None:
            _download(url, candidate_vcf)
            source = url
        else:
            shutil.copyfile(vcf_file, candidate_vcf)
            source = str(vcf_file)
        if candidate_vcf.stat().st_size == 0:
            raise RuntimeError("downloaded ClinVar VCF is empty")
        candidate_db = temp_dir / "clinvar_latest.sqlite"
        built = clinvar_latest_store.build_database(
            candidate_vcf,
            candidate_db,
            release_date=release_date,
        )
        if not built.get("release_date"):
            raise RuntimeError("ClinVar VCF has no parseable ##fileDate")
        with sqlite3.connect(candidate_db) as connection:
            check = connection.execute("PRAGMA quick_check").fetchone()[0]
            count = int(connection.execute("SELECT count(*) FROM records").fetchone()[0])
        if check != "ok" or count <= 0:
            raise RuntimeError(f"ClinVar database validation failed: {check}; rows={count}")

        vcf_path.parent.mkdir(parents=True, exist_ok=True)
        staged_vcf = vcf_path.with_suffix(vcf_path.suffix + ".tmp")
        staged_db = db_path.with_suffix(db_path.suffix + ".tmp")
        shutil.copyfile(candidate_vcf, staged_vcf)
        shutil.copyfile(candidate_db, staged_db)
        os.replace(staged_vcf, vcf_path)
        os.replace(staged_db, db_path)

        manifest = {
            "schema_version": 1,
            "status": "complete",
            "release_date": built["release_date"],
            "record_count": count,
            "source": source,
            "source_url": url,
            "sha256": _sha256(candidate_vcf),
            "updated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "database": str(db_path),
            "vcf": str(vcf_path),
        }
        clinvar_latest_store.write_manifest(manifest_path, manifest)
        return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf-file", type=Path)
    parser.add_argument("--url", default=config.CLINVAR_LATEST_URL)
    parser.add_argument("--db", type=Path, default=config.CLINVAR_LATEST_DB)
    parser.add_argument("--vcf-output", type=Path, default=config.CLINVAR_LATEST_VCF)
    parser.add_argument("--manifest", type=Path, default=config.CLINVAR_LATEST_MANIFEST_PATH)
    parser.add_argument("--lock", type=Path, default=config.CLINVAR_LATEST_LOCK_PATH)
    parser.add_argument("--release-date", default="")
    args = parser.parse_args()

    args.lock.parent.mkdir(parents=True, exist_ok=True)
    with args.lock.open("a+", encoding="utf-8") as lock_handle:
        try:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print("ClinVar update already running", file=sys.stderr)
            return 0
        try:
            manifest = update(
                vcf_file=args.vcf_file,
                url=args.url,
                db_path=args.db,
                vcf_path=args.vcf_output,
                manifest_path=args.manifest,
                release_date=args.release_date,
            )
        except Exception as exc:
            print(f"ERROR: ClinVar update failed: {exc}", file=sys.stderr)
            return 1
    print(json.dumps(manifest, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
