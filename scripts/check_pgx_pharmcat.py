#!/usr/bin/env python3
"""Read-only preflight; --install explicitly fetches the pinned public JAR."""
import argparse
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import urllib.request

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))
from app import config
from app.services import pgx_reannotation as pgx


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--install", action="store_true", help="download the official JAR only after SHA256 verification")
    args = parser.parse_args()
    if args.install:
        target = Path(os.environ.get("NGS_UI_PGX_PHARMCAT_JAR", config.BIOTOOLS_DIR / "pharmcat" / f"pharmcat-{pgx.VERSION}-all.jar"))
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            if hashlib.sha256(target.read_bytes()).hexdigest() != pgx.JAR_SHA256:
                raise ValueError("Existing JAR has a different checksum; no file was overwritten")
        else:
            url = f"https://github.com/PharmGKB/PharmCAT/releases/download/v{pgx.VERSION}/pharmcat-{pgx.VERSION}-all.jar"
            with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as f:
                tmp = Path(f.name)
                try:
                    with urllib.request.urlopen(url, timeout=90) as response:
                        while chunk := response.read(1024 * 1024):
                            f.write(chunk)
                    f.flush()
                    if hashlib.sha256(tmp.read_bytes()).hexdigest() != pgx.JAR_SHA256:
                        raise ValueError("Downloaded JAR checksum mismatch")
                    tmp.replace(target)
                finally:
                    tmp.unlink(missing_ok=True)
    engine = pgx.runtime()
    with tempfile.TemporaryDirectory(prefix="pgx-preflight-") as name:
        directory = Path(name); (directory / "tmp").mkdir()
        result = subprocess.run(pgx.command(engine, directory, ["-version"]), cwd=directory,
                                capture_output=True, text=True, timeout=45)
        if result.returncode or pgx.VERSION not in result.stdout + result.stderr:
            raise ValueError("PharmCAT could not start with the expected version: " + result.stderr)
    print(f"PASS: PharmCAT {pgx.VERSION}; official JAR SHA256 verified. No sample data changed.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1)
