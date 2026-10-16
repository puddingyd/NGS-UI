#!/usr/bin/env python3
"""Read-only inventory of CYP21A2, SMN, and HBA calls in old DRAGEN cases.

Writes tab-separated results to stdout. Only the small targeted JSON and
targeted VCF are read; this tool does not read a whole hard-filtered WGS VCF.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))
from app.services import gba, sample_layout, specialized

CALLERS = ("cyp21a2", "smn", "hba")
COLUMNS = ("sample_id", "source_sample_id", "caller", "status", "json_signals",
           "vcf_signals", "json_path", "targeted_vcf_path", "error")


def _positive(value):
    return type(value) in (int, float) and value > 0


def json_signals(caller: str, node: dict | None) -> list[str]:
    """Find actual JSON calls; a normal SMN2 copy count is not a disease allele."""
    if not isinstance(node, dict):
        return []
    signals = []
    variants = node.get("variants") if isinstance(node.get("variants"), list) else []
    if caller == "cyp21a2":
        haplotypes = node.get("recombinantHaplotypes") if isinstance(node.get("recombinantHaplotypes"), list) else []
        signals.extend("haplotype:" + h for h in haplotypes if isinstance(h, str) and h)
        signals.extend(
            f"variant:{v.get('alleleId', '?')} (copies={v.get('alleleCopyNumber')}, FILTER={v.get('filter', '?')})"
            for v in variants if isinstance(v, dict) and _positive(v.get("alleleCopyNumber"))
        )
        if node.get("deletionBreakpointInGene") is True:
            signals.append("deletion_breakpoint_in_gene:true")
        elif type(node.get("totalCopyNumber")) is int and node["totalCopyNumber"] != 4:
            signals.append(f"combined_CN:{node['totalCopyNumber']} (gene/pseudogene combined; review)")
    elif caller == "smn":
        cn = node.get("smn1CopyNumber")
        if type(cn) is int and cn != 2:
            signals.append(f"SMN1_CN:{cn}")
        smn2_cn = node.get("smn2CopyNumber")
        if type(smn2_cn) is int and smn2_cn != 1:
            signals.append(f"SMN2_CN:{smn2_cn} (copy-count test case; not an SMN1 disease allele)")
        for v in variants:
            if isinstance(v, dict) and _positive(v.get("altCopyNumber")):
                signals.append(f"marker:{v.get('hgvs', '?')} (ALT_CN={v['altCopyNumber']}; not diagnostic alone)")
    else:
        gt = node.get("genotype")
        if isinstance(gt, str) and gt not in ("", "None", "aa/aa"):
            signals.append(f"genotype:{gt} (FILTER={node.get('genotypeFilter', '?')})")
        if type(node.get("totalCopyNumber")) is int and node["totalCopyNumber"] != 4 and not signals:
            signals.append(f"combined_CN:{node['totalCopyNumber']} (check full genotype)")
        for v in variants:
            if isinstance(v, dict) and _positive(v.get("alleleCopyNumber")):
                signals.append(f"variant:{v.get('alleleId', '?')} (copies={v['alleleCopyNumber']}, FILTER={v.get('filter', '?')})")
    return signals


def vcf_signals(rows: list[dict]) -> dict[str, list[str]]:
    """Inventory GT-positive targeted VCF alleles, not a clinical classification."""
    out = {caller: [] for caller in CALLERS}
    seen = {caller: set() for caller in CALLERS}
    for allele in specialized.allele_evidence(rows):
        caller = allele["locus"]
        if caller not in out:
            continue
        aid = allele.get("allele_id") or ""
        if caller == "cyp21a2" and not aid.startswith("NM_000500."):
            continue
        if caller == "hba" and not aid.startswith(("NM_000517.", "NM_000558.")):
            continue
        # SMN is dominated by copy-number records; any active targeted VCF
        # allele is listed as an observation, never as proven SMN1 disease.
        key = (allele["id"], aid, allele["filter"])
        if key in seen[caller]:
            continue
        seen[caller].add(key)
        out[caller].append(f"{aid or allele['id']} (GT={allele['gt']}, FILTER={allele['filter']})")
    return out


def scan_case(sample_id: str, *, fallback_build: str = "hg38", include_normal: bool = False) -> list[dict]:
    source_path = sample_layout.state_file(sample_id, "pipeline_source.json")
    try:
        source = json.loads(source_path.read_text(encoding="utf-8"))
        source_id = source["source_sample_id"]
        source_vcf = Path(source["source_vcf_path"])
        if not source_id or not source_vcf.name:
            raise ValueError("missing source sample or VCF path")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return [{"sample_id": sample_id, "source_sample_id": "", "caller": "*", "status": "error",
                 "json_signals": "", "vcf_signals": "", "json_path": "", "targeted_vcf_path": "", "error": f"pipeline_source.json: {exc}"}]

    base = {"sample_id": sample_id, "source_sample_id": source_id}
    targeted_vcf = source_vcf.parent / f"{source_id}.targeted.vcf.gz"
    json_path = None
    payload = {}
    errors = []
    try:
        json_path = gba.discover_json(source_vcf, source_id)
        if json_path:
            payload = json.loads(specialized._read_small(json_path))
            if not isinstance(payload, dict):
                raise ValueError("targeted JSON must be an object")
            if payload.get("sampleId") != source_id:
                raise ValueError("targeted JSON sampleId does not match source sample")
    except (OSError, ValueError, TypeError) as exc:
        errors.append(f"JSON: {exc}")
    build = {"GRCh38": "hg38", "GRCh37": "hg19"}.get(payload.get("genomeBuild"), payload.get("genomeBuild")) or fallback_build
    if build not in specialized.REGIONS:
        errors.append(f"unsupported build: {build}")
    vcf = {caller: [] for caller in CALLERS}
    if targeted_vcf.is_file() and build in specialized.REGIONS:
        try:
            _, rows = specialized._extract(targeted_vcf, "targeted", source_id, build)
            vcf = vcf_signals(rows)
        except (OSError, ValueError) as exc:
            errors.append(f"targeted VCF: {exc}")

    output = []
    for caller in CALLERS:
        js = json_signals(caller, payload.get(caller))
        vs = vcf[caller]
        if not (js or vs or errors or include_normal):
            continue
        status = "candidate" if (js or vs) else "error" if errors else "not_available" if not json_path and not targeted_vcf.is_file() else "no_signal"
        if status == "candidate" and errors:
            status = "candidate_partial"
        output.append({**base, "caller": caller, "status": status, "json_signals": "; ".join(js),
                       "vcf_signals": "; ".join(vs), "json_path": str(json_path or ""),
                       "targeted_vcf_path": str(targeted_vcf) if targeted_vcf.is_file() else "",
                       "error": "; ".join(errors)})
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", help="one exact tertiary sample ID, for a quick check")
    parser.add_argument("--include-normal", action="store_true", help="also list no-signal caller rows")
    parser.add_argument("--fallback-build", choices=("hg38", "hg19"), default="hg38",
                        help="used only when targeted JSON does not provide a genome build")
    args = parser.parse_args()
    if args.sample and (not args.sample.endswith("-dragen") or Path(args.sample).name != args.sample):
        parser.error("--sample must be an exact -dragen tertiary ID")
    samples = [args.sample] if args.sample else sorted(s for s in sample_layout.iter_sample_ids() if s.endswith("-dragen"))
    writer = csv.DictWriter(sys.stdout, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    counts = {"samples": 0, "candidate_rows": 0, "error_rows": 0}
    for sample in samples:
        counts["samples"] += 1
        for row in scan_case(sample, fallback_build=args.fallback_build, include_normal=args.include_normal):
            writer.writerow(row)
            if row["status"].startswith("candidate"):
                counts["candidate_rows"] += 1
            if row["error"]:
                counts["error_rows"] += 1
    print(json.dumps(counts, ensure_ascii=False), file=sys.stderr)


if __name__ == "__main__":
    main()
