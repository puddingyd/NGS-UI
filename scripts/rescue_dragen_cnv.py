#!/usr/bin/env python3
"""Annotate DRAGEN-rescued CNVs inside a job's private 08 staging directory."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))
from app.config import NGS_UI_HOME
from app.services.dragen_cnv_rescue import annotsv_command, build_rescue


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dragen-vcf", type=Path, required=True)
    parser.add_argument("--base-tsv", type=Path, required=True)
    parser.add_argument("--post-dir", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--source-sample", required=True)
    args = parser.parse_args()
    suffix = ".hard-filtered.vcf.gz"
    if not args.dragen_vcf.name.endswith(suffix):
        parser.error("--dragen-vcf must be the DRAGEN .hard-filtered.vcf.gz anchor")
    prefix = args.dragen_vcf.name[:-len(suffix)]

    def annotate(vcf: Path, output: Path) -> None:
        subprocess.run(annotsv_command(vcf, output, NGS_UI_HOME), check=True, cwd=vcf.parent)

    print("[post-processing-step] cnv-rescue start", flush=True)
    result = build_rescue(
        raw_cnv=args.dragen_vcf.with_name(prefix + ".cnv.vcf.gz"),
        joint_cnv=args.dragen_vcf.with_name(prefix + ".cnv_sv.vcf.gz"),
        base_tsv=args.base_tsv, post_dir=args.post_dir,
        sample_id=args.sample, source_sample=args.source_sample, annotate=annotate,
    )
    if result["status"] == "skipped":
        print("[cnv-rescue] WARNING: skipped; missing input: " + ", ".join(result["missing"]), flush=True)
    else:
        counts = result["counts"]
        print(f"[cnv-rescue] A={counts.get('rescued_A', 0)} B={counts.get('rescued_B', 0)} "
              f"added={result['added_events']}", flush=True)
        if counts.get("unmatched_or_ambiguous_original"):
            print(f"[cnv-rescue] WARNING: {counts['unmatched_or_ambiguous_original']} "
                  "integrated CNVs could not be uniquely mapped to the original CNV", flush=True)
    print("[post-processing-step] cnv-rescue done", flush=True)


if __name__ == "__main__":
    main()
