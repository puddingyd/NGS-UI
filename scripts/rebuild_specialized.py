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


def rebuild_one(sample: str, *, from_source: bool, apply: bool, source_vcf: Path | None = None) -> dict:
    """Process one case under the same lock as tertiary analysis."""
    locks = _acquire_sample_locks([sample])
    try:
        manifest = sample_layout.state_file(sample, specialized.NAMES["manifest"])
        source_path = sample_layout.state_file(sample, "pipeline_source.json")
        dest = manifest.parent
        if not dest.is_dir():
            raise ValueError("請先完成個案三級分析")
        raw = sample_layout.snv_raw_tsv(sample)
        if not raw.is_file():
            raise ValueError("三級 SNV 註解不存在")
        with tempfile.TemporaryDirectory(prefix=".specialized-", dir=dest) as folder:
            stage = Path(folder)
            names = set(specialized.MANAGED_NAMES)
            if from_source:
                source = json.loads(source_path.read_text())
                vcf = source_vcf or Path(source.get("source_vcf_path") or "")
                source_id = source.get("source_sample_id")
                if not source_id or not vcf.is_file():
                    raise ValueError("缺少明確 source sample 或來源 VCF")
                result = specialized.prepare_outputs(mode=source.get("pipeline_type", "dragen"), source_vcf=vcf,
                    source_sample_id=source_id, sample_id=sample, post_dir=stage, raw_tsv=raw)
                if result["status"] != "complete":
                    raise ValueError("未找到 targeted JSON；保留既有結果")
                source["source_vcf_size"] = vcf.stat().st_size
                (stage / f"{sample}.pipeline_source.json").write_text(json.dumps(source, ensure_ascii=False, indent=2))
                names.add("pipeline_source.json")
            else:
                result = specialized.rebuild_local(sample_id=sample, post_dir=stage, raw_tsv=raw)
            summary = {"sample": sample, "status": "applied" if apply else "validated", "findings": len(result["variants"]),
                       "loci": {k: v["status"] for k, v in result["loci"].items()}}
            if apply:
                publish(stage, dest, sample, names=names)
            return summary
    finally:
        _release_sample_locks(locks)


def batch_samples(*, include_existing: bool = False):
    """Enumerate tertiary cases only; never glob the raw sequencing tree."""
    for sample in sorted(sample_layout.iter_sample_ids()):
        if not sample.endswith("-dragen"):
            continue
        if not include_existing and sample_layout.state_file(sample, specialized.NAMES["manifest"]).is_file():
            try:
                if specialized.load(sample).get("status") == "complete":
                    yield sample, "skipped_existing"
                    continue
            except (OSError, ValueError, KeyError):
                pass  # A damaged bundle must be rebuilt, not silently skipped.
        yield sample, "pending"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--sample")
    target.add_argument("--all-dragen", action="store_true", help="all tertiary cases with the -dragen suffix")
    parser.add_argument("--from-source", action="store_true", help="explicitly read original targeted JSON/VCF")
    parser.add_argument("--source-vcf", type=Path, help="override one relocated VCF path; requires --sample --from-source")
    parser.add_argument("--include-existing", action="store_true", help="also rebuild valid existing bundles in batch mode")
    parser.add_argument("--apply", action="store_true", help="publish; default validates in temporary staging only")
    args = parser.parse_args()
    if args.sample and (Path(args.sample).name != args.sample or args.sample in {".", ".."}):
        parser.error("invalid sample ID")
    if args.source_vcf and (not args.from_source or args.all_dragen):
        parser.error("--source-vcf requires --sample and --from-source")
    if args.include_existing and not args.all_dragen:
        parser.error("--include-existing requires --all-dragen")
    if args.all_dragen and not args.from_source:
        parser.error("--all-dragen requires --from-source for old cases")
    if args.sample:
        print(json.dumps(rebuild_one(args.sample, from_source=args.from_source, apply=args.apply,
                                     source_vcf=args.source_vcf), ensure_ascii=False), flush=True)
        return
    counts = {"validated": 0, "applied": 0, "skipped_existing": 0, "error": 0}
    for sample, status in batch_samples(include_existing=args.include_existing):
        if status == "skipped_existing":
            row = {"sample": sample, "status": status}
        else:
            try:
                row = rebuild_one(sample, from_source=True, apply=args.apply)
            except (OSError, ValueError, RuntimeError, KeyError, json.JSONDecodeError) as exc:
                row = {"sample": sample, "status": "error", "error": str(exc)}
        counts[row["status"]] += 1
        print(json.dumps(row, ensure_ascii=False), flush=True)
    print(json.dumps({"summary": counts}, ensure_ascii=False), flush=True)
    if counts["error"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
