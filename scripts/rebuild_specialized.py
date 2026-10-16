#!/usr/bin/env python3
"""Explicit targeted-caller backfill. No Nextflow/VEP rerun, no service restart.

Default: rebuild from local snapshots. --from-source explicitly reads original
JSON/VCF using the recorded exact sample. Writes only the four bundle files and
source-size metadata, preserving reviewer state and all 00-07 outputs.
"""
from __future__ import annotations
import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))
from app.services import sample_layout, specialized
from app.workers.dragen_run import _acquire_sample_locks, _release_sample_locks


def publish(staging: Path, dest: Path, sample: str, *, names=None):
    """Same-filesystem promotion, manifest last; rollback on any failure."""
    names = list(names or specialized.MANAGED_NAMES)
    names.sort(key=lambda name: name == specialized.NAMES["manifest"])
    backup = staging / ".backup"
    backup.mkdir()
    completed = []
    try:
        for name in names:
            target = dest / f"{sample}.{name}"
            old = backup / target.name
            incoming = staging / target.name
            if target.exists():
                os.replace(target, old)
            completed.append((target, old))
            os.replace(incoming, target)
    except BaseException:
        for target, old in reversed(completed):
            if target.exists():
                target.unlink()
            if old.exists():
                os.replace(old, target)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--from-source", action="store_true", help="explicitly read original targeted JSON/VCF")
    parser.add_argument("--source-vcf", type=Path, help="override relocated VCF path, requires --from-source")
    parser.add_argument("--apply", action="store_true", help="publish; default validates in temporary staging only")
    args = parser.parse_args()
    if not args.sample or Path(args.sample).name != args.sample or args.sample in {".", ".."}:
        parser.error("invalid sample ID")
    if args.source_vcf and not args.from_source:
        parser.error("--source-vcf requires --from-source")
    locks = _acquire_sample_locks([args.sample])
    try:
        manifest = sample_layout.state_file(args.sample, specialized.NAMES["manifest"])
        source_path = sample_layout.state_file(args.sample, "pipeline_source.json")
        dest = manifest.parent
        if not dest.is_dir():
            raise ValueError("請先完成個案三級分析")
        raw = sample_layout.snv_raw_tsv(args.sample)
        if not raw.is_file():
            raise ValueError("三級 SNV 註解不存在")
        with tempfile.TemporaryDirectory(prefix=".specialized-", dir=dest) as folder:
            stage = Path(folder)
            names = set(specialized.MANAGED_NAMES)
            if args.from_source:
                source = json.loads(source_path.read_text())
                vcf = args.source_vcf or Path(source.get("source_vcf_path") or "")
                source_id = source.get("source_sample_id")
                if not source_id or not vcf.is_file():
                    raise ValueError("缺少明確 source sample 或來源 VCF")
                result = specialized.prepare_outputs(mode=source.get("pipeline_type", "dragen"), source_vcf=vcf,
                    source_sample_id=source_id, sample_id=args.sample, post_dir=stage, raw_tsv=raw)
                if result["status"] != "complete":
                    raise ValueError("未找到 targeted JSON；保留既有結果")
                source["source_vcf_size"] = vcf.stat().st_size
                (stage / f"{args.sample}.pipeline_source.json").write_text(json.dumps(source, ensure_ascii=False, indent=2))
                names.add("pipeline_source.json")
            else:
                result = specialized.rebuild_local(sample_id=args.sample, post_dir=stage, raw_tsv=raw)
            print(json.dumps({"sample": args.sample, "apply": args.apply, "findings": len(result["variants"]),
                              "loci": {k: v["status"] for k, v in result["loci"].items()}}, ensure_ascii=False))
            if args.apply:
                publish(stage, dest, args.sample, names=names)
    finally:
        _release_sample_locks(locks)


if __name__ == "__main__":
    main()
