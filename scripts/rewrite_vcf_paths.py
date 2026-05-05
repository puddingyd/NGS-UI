#!/usr/bin/env python3
"""Rewrite the `vcf_path` field in every tertiary_output sample_metadata.json
from the old VCF directory to the new one.

Usage:
    python3 scripts/rewrite_vcf_paths.py \\
        --root /home/n102968/NGS_UI/tertiary_output \\
        --old  /home/n102968/vcf \\
        --new  /home/n102968/NGS_UI/vcf

Idempotent: a path that already matches `--new` is left as-is.
"""
import argparse
import json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True,
                    help="tertiary_output directory to scan")
    ap.add_argument("--old", required=True, help="old VCF prefix to replace")
    ap.add_argument("--new", required=True, help="new VCF prefix")
    args = ap.parse_args()

    root = Path(args.root)
    if not root.is_dir():
        raise SystemExit(f"not a directory: {root}")

    changed = 0
    for meta in sorted(root.glob("*/sample_metadata.json")):
        try:
            data = json.loads(meta.read_text(encoding="utf-8"))
        except Exception as exc:
            print(f"skip (parse error): {meta}: {exc}")
            continue
        vcf = data.get("vcf_path", "")
        if not vcf or not vcf.startswith(args.old):
            continue
        data["vcf_path"] = vcf.replace(args.old, args.new, 1)
        meta.write_text(
            json.dumps(data, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        changed += 1
        print(f"updated: {meta.parent.name}: {vcf} -> {data['vcf_path']}")

    print(f"\n{changed} file(s) updated under {root}")


if __name__ == "__main__":
    main()
