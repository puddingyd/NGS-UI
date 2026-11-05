#!/usr/bin/env bash
# Copy this entire file into the server shell. No repo update or sudo required.
bash <<'NGS_FILE_PROBE_SHELL'
probe_dir=$(mktemp -d /tmp/ngs-file-probe.XXXXXX) || exit 1
probe_args=()
probe_cfg=${NGS_UI_PRS_CONFIG:-${NGS_UI_DATA_ROOT:-${NGS_UI_HOME:-$HOME/NGS_UI}/data}/prs_config.json}
if [ -f "$probe_cfg" ]; then
  probe_args+=(--prs-config "$probe_cfg")
fi
# Add e.g. probe_args+=(--root /extra/root) for additional known data roots.
python3 - "${probe_args[@]}" > "$probe_dir/inventory.json" 2> "$probe_dir/scan.log" <<'NGS_FILE_PROBE_PY'
#!/usr/bin/env python3
"""Standalone, read-only path probe. No app imports, cache writes or file reads.

Only directory entries/stat, source markers and the optional PRS config are read.
JSON goes to stdout; progress/errors go to stderr. Run from any directory.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import sys


RAW = "/home/datalake_Raw"
OUTPUT = "/home/datalake_Intermediate/pipeline/nextflow_output"
DEFAULTS = {
    "NGS_UI_SECONDARY_WES_FASTQ_ROOTS": f"{RAW}/NextSeq2000:{RAW}/Other/Reanalysis:/datalake_Raw/Other/Reanalysis",
    "NGS_UI_SECONDARY_WGS_FASTQ_ROOTS": f"{RAW}/Novaseq",
    "NGS_UI_INHOUSE_VCF_ROOTS": OUTPUT,
    "NGS_UI_SECONDARY_OUTPUT_ROOT": OUTPUT,
    "NGS_UI_INHOUSE_BAM_ROOT": OUTPUT,
    "NGS_UI_DRAGEN_VCF_ROOTS": f"{RAW}/Novaseq:/home/datalake_Intermediate/n102968",
    "NGS_UI_DRAGEN_BAM_ROOT": f"{RAW}/Novaseq:/home/datalake_Intermediate/n102968",
    "NGS_UI_TERTIARY_ROOT": "/home/datalake_Intermediate/pipeline/tertiary_output",
}
LANE = re.compile(r"^(.+)_S\d+_L(\d{3})_R([12])_001\.fastq\.gz$", re.I)
WES = re.compile(r"^(.+)_S\d+_R([12])_001\.fastq\.gz$", re.I)
REANALYSIS = re.compile(r"^(.+)\.R([12])\.clean\.fastq\.gz$", re.I)
SUFFIXES = {
    ".hard-filtered.gvcf.gz": "gVCF", ".gvcf.gz": "gVCF", ".g.vcf.gz": "gVCF",
    ".hard-filtered.vcf.gz": "SNV_VCF", ".ensemble.fixed.vcf.gz": "SNV_VCF",
    ".cnv_sv.vcf.gz": "CNV_SV_VCF", ".cnv.vcf.gz": "CNV_VCF",
    ".sv.vcf.gz": "SV_VCF", ".gcnv.vcf.gz": "CNV_VCF",
    ".delly.vcf.gz": "SV_VCF", ".mito.vcf.gz": "MITO_VCF",
    ".call.cns": "CNVKIT_CNS", ".cnr": "CNVKIT_CNR",
    ".cnv.annotated.tsv": "CNV_TSV", ".sv.annotated.tsv": "SV_TSV",
    ".cnv.review.tsv": "CNV_REVIEW_TSV", ".cnv.rescued.annotated.tsv": "CNV_RESCUE_TSV",
    ".targeted.json": "TARGETED_JSON",
    ".ploidy.vcf.gz": "DERIVED_VCF", ".vcf_from_tsv.vcf.gz": "DERIVED_VCF",
    ".repeats.bam": "REPEATS_BAM", ".aligned.sorted.bam": "BAM", ".bam": "BAM",
}
PATTERNS = (".fastq.gz", ".fq.gz", ".bam", ".bai", ".vcf", ".vcf.gz",
            ".gvcf.gz", ".tbi", ".csi", ".call.cns", ".cnr", ".cnv.annotated.tsv",
            ".sv.annotated.tsv", ".cnv.review.tsv", ".cnv.rescued.annotated.tsv", ".targeted.json")


def configured_roots():
    roots = {}
    for key, default in DEFAULTS.items():
        fallback = default
        if key == "NGS_UI_INHOUSE_BAM_ROOT":
            fallback = os.environ.get("NGS_UI_BAM_ROOT", default)
        elif key == "NGS_UI_TERTIARY_ROOT":
            fallback = os.environ.get("NGS_UI_PIPELINE_OUT_ROOT", default)
        roots[key] = [str(Path(p).expanduser()) for p in os.environ.get(key, fallback).split(":") if p]
    return roots


def scan(roots, sample_filter="", prs_config=None):
    errors, unknown_fastqs, files = [], [], []
    samples, aliases, ui_sources = {}, defaultdict(set), {}
    seen_dirs, seen_files = set(), set()

    def error(path, exc):
        errors.append({"path": str(path), "error": str(exc)})
        print(f"WARNING {path}: {exc}", file=sys.stderr)

    def read_object(path):
        if path.stat().st_size > 2_000_000:
            raise ValueError("JSON exceeds 2 MB")
        value = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(value, dict):
            raise ValueError("JSON must be an object")
        return value

    # One walk for each physical directory, even when configured roots overlap.
    # Directory symlinks are reported instead of silently skipped/followed.
    def walk(folder):
        try:
            st = folder.stat()
            key = (st.st_dev, st.st_ino)
            if key in seen_dirs:
                return
            seen_dirs.add(key)
            with os.scandir(folder) as entries:
                for entry in entries:
                    path = Path(entry.path)
                    try:
                        if entry.is_dir(follow_symlinks=False):
                            yield from walk(path)
                        elif entry.is_symlink() and path.is_dir():
                            error(path, "directory symlink not followed; add its target using --root")
                        elif folder.name == "08_postprocessing" and entry.name in (
                            f"{folder.parent.name}.pipeline_source.json", "pipeline_source.json"
                        ):
                            preferred = folder / f"{folder.parent.name}.pipeline_source.json"
                            if entry.name == "pipeline_source.json" and preferred.exists():
                                continue
                            try:
                                source = read_object(path).get("source_sample_id")
                                if isinstance(source, str) and source:
                                    ui_sources[str(folder.parent)] = source
                            except (OSError, ValueError) as exc:
                                error(path, exc)
                        elif entry.name.lower().endswith(PATTERNS):
                            yield path
                    except OSError as exc:
                        error(path, exc)
        except OSError as exc:
            error(folder, exc)

    for root in sorted({p for paths in roots.values() for p in paths}):
        print(f"SCAN {root}", file=sys.stderr)
        for path in walk(Path(root)):
            canonical = str(path.absolute())
            if canonical in seen_files:
                continue
            seen_files.add(canonical)
            name = path.name
            try:
                st = path.stat()
            except OSError as exc:
                error(path, exc)
                continue
            record = {"path": str(path), "bytes": st.st_size, "mtime": st.st_mtime}
            if name.lower().endswith((".fastq.gz", ".fq.gz")):
                lane, wes, reanalysis = LANE.match(name), WES.match(name), REANALYSIS.match(name)
                match = lane or wes or reanalysis
                if not match:
                    unknown_fastqs.append(record)
                    continue
                source = match.group(1)
                sid = path.parent.name if reanalysis else source
                item = samples.setdefault(sid, {"sample_id": sid, "source_ids": [], "fastqs": []})
                aliases[source].add(sid)
                if source not in item["source_ids"]:
                    item["source_ids"].append(source)
                read = match.group(3) if lane else match.group(2)
                record.update(source_id=source, lane=f"L{lane.group(2)}" if lane else "", read=read)
                item["fastqs"].append(record)
                continue
            index_suffix = next((s for s in (".tbi", ".csi") if name.endswith(s)), "")
            bare = name[:-len(index_suffix)] if index_suffix else name
            # Handle both sample.bam.bai and sample.bai.
            bai = bare.endswith(".bai")
            if bai:
                bare = bare[:-4]
                if not bare.endswith(".bam"):
                    bare += ".bam"
            suffix = next((s for s in SUFFIXES if bare.endswith(s)), "")
            source = bare[:-len(suffix)] if suffix else ""
            kind = "BAI" if bai else ("VCF_INDEX" if index_suffix else SUFFIXES.get(suffix, "OTHER"))
            record.update(source_id=source, kind=kind)
            files.append(record)

    # Inspect templates without importing app.config (which creates runtime dirs).
    templates = None
    if prs_config:
        try:
            cfg = read_object(Path(prs_config))
            templates = cfg.get("gvcf_templates")
            if templates is not None and not isinstance(templates, dict):
                raise ValueError("gvcf_templates must be a mapping")
        except (OSError, ValueError, AttributeError) as exc:
            error(prs_config, exc)
    # Custom exact gVCF paths may live outside all scanned roots. Same rendering
    # as prs_jobs.candidates; custom mapping replaces defaults in the real PRS.
    for anchor in list(files):
        if templates is None or anchor["kind"] != "SNV_VCF":
            continue
        path = Path(anchor["path"])
        mode = "inhouse" if path.name.endswith(".ensemble.fixed.vcf.gz") else "dragen"
        values = templates.get(mode, [])
        if not isinstance(values, list) or any(not isinstance(t, str) for t in values):
            error(prs_config, f"invalid {mode} gvcf_templates")
            continue
        for template in values:
            try:
                candidate = Path(template.format(anchor_dir=path.parent, run_dir=path.parent.parent,
                                                 source=anchor["source_id"]))
                if not candidate.is_file() or str(candidate.absolute()) in seen_files:
                    continue
                seen_files.add(str(candidate.absolute()))
                st = candidate.stat()
                files.append(dict(path=str(candidate), source_id=anchor["source_id"], kind="gVCF",
                                  bytes=st.st_size, mtime=st.st_mtime, via="prs_config"))
                for ext in (".tbi", ".csi"):
                    idx = Path(str(candidate) + ext)
                    if idx.is_file() and str(idx.absolute()) not in seen_files:
                        seen_files.add(str(idx.absolute()))
                        files.append(dict(path=str(idx), source_id=anchor["source_id"], kind="VCF_INDEX",
                                          bytes=idx.stat().st_size, mtime=idx.stat().st_mtime, via="prs_config"))
            except (OSError, KeyError, ValueError, IndexError) as exc:
                error(template, exc)

    files_by_sample = defaultdict(list)
    for record in files:
        path = Path(record["path"])
        for parent in path.parents:
            source = ui_sources.get(str(parent))
            if source and record["source_id"] == parent.name:
                record.update(source_id=source, ui_sample_id=parent.name, via="pipeline_source.json")
                break
        for sid in aliases.get(record["source_id"], set()):
            files_by_sample[sid].append(record)
    for sid, item in sorted(samples.items()):
        pairs = defaultdict(set)
        for fastq in item["fastqs"]:
            pairs[(str(Path(fastq["path"]).parent), fastq["source_id"], fastq["lane"])].add(fastq["read"])
        item["fastq_count"] = len(item["fastqs"])
        item["fastqs"].sort(key=lambda f: f["path"])
        item["incomplete_pairs"] = [dict(directory=k[0], source_id=k[1], lane=k[2], reads=sorted(v))
                                    for k, v in sorted(pairs.items()) if v != {"1", "2"}]
        item["files"] = sorted(files_by_sample[sid], key=lambda f: f["path"])
    # An exact source ID match groups candidates only: it does NOT prove that
    # output belongs to a particular FASTQ run. Every path is preserved.
    selected = [v for sid, v in sorted(samples.items()) if not sample_filter or sid == sample_filter
                or sample_filter in v["source_ids"]]
    return {"generated_at": datetime.now(timezone.utc).isoformat(), "roots": roots,
            "association": "exact source ID candidates; run identity not verified",
            "samples": selected, "total_fastq_samples": len(samples),
            "unmatched_files": [f for f in files if not aliases.get(f["source_id"])],
            "unrecognized_fastqs": unknown_fastqs, "scan_errors": errors,
            "prs_gvcf_templates": templates,
            "limitations": ["File existence/stat only; no content, integrity or sample-header validation.",
                            "Directory symlinks not followed; use --root TARGET when reported.",
                            "Missing roots/permission errors mean an incomplete scan.",
                            "Unrecognized FASTQ names require mapping before the final inventory.",
                            "FASTQ-only samples and incomplete R1/R2 pairs are retained."]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", default="", help="Filter displayed samples; still scans all configured roots")
    parser.add_argument("--root", action="append", default=[], help="Additional root (repeatable)")
    parser.add_argument("--prs-config", default=os.environ.get("NGS_UI_PRS_CONFIG"), help="Read only gvcf_templates")
    args = parser.parse_args()
    roots = configured_roots()
    roots["extra"] = args.root
    report = scan(roots, args.sample, args.prs_config)
    json.dump(report, sys.stdout, ensure_ascii=False, indent=2)
    print()
    print(f"FASTQ samples={report['total_fastq_samples']} shown={len(report['samples'])} "
          f"unrecognized_fastqs={len(report['unrecognized_fastqs'])} "
          f"unmatched_files={len(report['unmatched_files'])} scan_errors={len(report['scan_errors'])}", file=sys.stderr)
    return 2 if report["scan_errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())

NGS_FILE_PROBE_PY
probe_status=$?
cat "$probe_dir/scan.log"
python3 - "$probe_dir/inventory.json" <<'NGS_FILE_PROBE_SUMMARY'
import json, sys
from pathlib import Path
p = Path(sys.argv[1])
r = json.loads(p.read_text())
print("\nReport:", p)
print("FASTQ samples:", r['total_fastq_samples'])
print("Unrecognized FASTQs:", len(r['unrecognized_fastqs']))
print("Files without FASTQ match:", len(r['unmatched_files']))
print("Scan errors:", len(r['scan_errors']))
print("PRS custom templates read:", r['prs_gvcf_templates'] is not None)
for s in r['samples']:
    if s['sample_id'] == '26T01345' or '26T01345' in s['source_ids']:
        print("\nSample", s['sample_id'], "FASTQs:", s['fastq_count'])
        for f in s['fastqs']:
            print(f["lane"], "R" + f["read"], f["path"])
        for f in s['files']:
            print(f["kind"], f["path"])
NGS_FILE_PROBE_SUMMARY
printf '\nProbe exit code: %s (2 means incomplete scan; inspect scan_errors).\n' "$probe_status"
printf 'Service setting locations (read only):\n'
systemctl show ngs-ui.service -p FragmentPath -p DropInPaths --no-pager 2>/dev/null || true
printf 'Current index timer status (read only):\n'
systemctl list-timers --all --no-pager 'ngs-ui-*index*' 2>/dev/null || true
NGS_FILE_PROBE_SHELL
