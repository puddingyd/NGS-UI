#!/usr/bin/env python3
"""Backfill only SpliceAI onto a tertiary SNV working TSV.

dbNSFP annotation belongs to the v3.6 Nextflow pipeline.  This post-processing
step deliberately invokes VEP with only the SpliceAI plugin and writes one
derived column, ``SPLICEAI_MAX``.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

DEFAULT_REF_DIR = "/home/pipeline/reference/hg38"
DEFAULT_VEP_SIF = "/home/pipeline/nextflow_containers/vep_115.sif"
DEFAULT_NGS_UI_HOME = os.environ.get("NGS_UI_HOME") or str(Path.home() / "NGS_UI")
DEFAULT_SPLICEAI_SNV = (
    f"{DEFAULT_NGS_UI_HOME}/biotools/spliceai/"
    "spliceai_scores.raw.snv.hg38.vcf.gz"
)
DEFAULT_SPLICEAI_INDEL = (
    f"{DEFAULT_NGS_UI_HOME}/biotools/spliceai/"
    "spliceai_scores.raw.indel.hg38.vcf.gz"
)
TSV_COL_SPLICEAI = "SPLICEAI_MAX"
BedIndex = dict[str, tuple[list[int], list[tuple[int, int]]]]


def _open_vcf(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open("r")


def _normalize_chrom(chrom: str) -> str:
    value = str(chrom or "").strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    if value.upper() in {"M", "MT"}:
        return "MT"
    return value.upper() if value.upper() in {"X", "Y"} else value


def _row_max_af(row: dict[str, str], columns: list[str]) -> float:
    maximum = 0.0
    for column in columns:
        raw = str(row.get(column) or "").strip()
        if not raw or raw.upper() in {".", "NA", "N/A"}:
            continue
        try:
            maximum = max(maximum, float(raw))
        except ValueError:
            pass
    return maximum


def load_bed(path: Path) -> BedIndex:
    raw: dict[str, list[tuple[int, int]]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 3:
                continue
            try:
                start, end = int(columns[1]), int(columns[2])
            except ValueError:
                continue
            if end > start:
                raw.setdefault(_normalize_chrom(columns[0]), []).append((start, end))
    result: BedIndex = {}
    for chrom, intervals in raw.items():
        intervals.sort()
        merged: list[tuple[int, int]] = []
        for start, end in intervals:
            if not merged or start > merged[-1][1]:
                merged.append((start, end))
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        result[chrom] = ([start for start, _end in merged], merged)
    return result


def _overlaps_bed(bed: BedIndex | None, chrom: str, pos: str, ref: str) -> bool:
    if bed is None:
        return True
    try:
        start = int(pos) - 1
    except ValueError:
        return False
    end = start + max(1, len(ref or ""))
    indexed = bed.get(_normalize_chrom(chrom))
    if not indexed:
        return False
    starts, intervals = indexed
    index = bisect.bisect_right(starts, start) - 1
    if index >= 0 and intervals[index][1] > start:
        return True
    index += 1
    return index < len(intervals) and intervals[index][0] < end


def tsv_to_sites(
    tsv_in: Path,
    vcf_out: Path,
    *,
    max_af: float | None,
    af_cols: list[str],
    candidate_bed: BedIndex | None = None,
) -> tuple[int, int, int]:
    seen: set[tuple[str, str, str, str]] = set()
    rows: list[tuple[str, int, str, str]] = []
    dropped_af = dropped_bed = 0
    with tsv_in.open("r", encoding="utf-8", newline="") as source:
        for row in csv.DictReader(source, delimiter="\t"):
            chrom, pos, ref, alt = (
                str(row.get(field) or "").strip()
                for field in ("CHROM", "POS", "REF", "ALT")
            )
            if not all((chrom, pos, ref, alt)) or "*" in {ref, alt}:
                continue
            key = (chrom, pos, ref, alt)
            if key in seen:
                continue
            if max_af is not None and _row_max_af(row, af_cols) > max_af:
                dropped_af += 1
                continue
            if not _overlaps_bed(candidate_bed, chrom, pos, ref):
                dropped_bed += 1
                continue
            try:
                numeric_pos = int(pos)
            except ValueError:
                continue
            seen.add(key)
            rows.append((chrom, numeric_pos, ref, alt))
    rows.sort(key=lambda row: (_normalize_chrom(row[0]), row[1], row[2], row[3]))
    with vcf_out.open("w", encoding="utf-8") as destination:
        destination.write("##fileformat=VCFv4.2\n")
        destination.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in rows:
            destination.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
    return len(rows), dropped_af, dropped_bed


def _parse_csq_format(vcf_path: Path) -> list[str]:
    pattern = re.compile(
        r'##INFO=<ID=CSQ,.*Description="[^"]*Format:\s*([^"]+)"',
        re.IGNORECASE,
    )
    with _open_vcf(vcf_path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            match = pattern.match(line.rstrip("\n"))
            if match:
                return match.group(1).split("|")
    return []


def run_vep(args: argparse.Namespace, sites: Path, vep_out: Path) -> None:
    binds = [
        args.ref_dir,
        str(sites.parent),
        str(Path(args.spliceai_snv).parent),
        str(Path(args.spliceai_indel).parent),
    ]
    binds = list(dict.fromkeys(binds))
    cmd = [
        "apptainer", "exec",
        "--bind", ",".join(binds),
        args.vep_sif,
        "vep",
        "--input_file", str(sites),
        "--output_file", str(vep_out),
        "--vcf", "--compress_output", "bgzip",
        "--offline", "--cache",
        "--dir_cache", args.vep_cache,
        "--dir_plugins", "/opt/vep/Plugins",
        "--assembly", "GRCh38",
        "--fasta", args.ref_fasta,
        "--fork", str(args.fork),
        "--symbol", "--canonical", "--biotype", "--mane", "--flag_pick",
        "--pick_order",
        "mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length",
        "--plugin", f"SpliceAI,snv={args.spliceai_snv},indel={args.spliceai_indel}",
        "--force_overwrite", "--no_stats", "--safe",
    ]
    print(f"[spliceai] running VEP on {sites.name}", file=sys.stderr)
    proc = subprocess.run(cmd, check=False, stderr=subprocess.PIPE, text=True)
    if proc.stderr:
        noisy = re.compile(
            r"(Use of uninitialized value \$alt_allele in string eq at /plugins/SpliceAI\.pm|"
            r"^WARNING: \d+ : Use of uninitialized value \$alt_allele)"
        )
        for line in proc.stderr.splitlines():
            if not noisy.search(line):
                print(line, file=sys.stderr)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)


def parse_vep_vcf(vep_vcf: Path) -> dict[tuple[str, str, str, str], str]:
    """Return the maximum absolute SpliceAI delta score for each site."""
    fields = _parse_csq_format(vep_vcf)
    if not fields:
        return {}
    indexes = [
        fields.index(name)
        for name in (
            "SpliceAI_pred_DS_AG",
            "SpliceAI_pred_DS_AL",
            "SpliceAI_pred_DS_DG",
            "SpliceAI_pred_DS_DL",
        )
        if name in fields
    ]
    pick_index = fields.index("PICK") if "PICK" in fields else -1
    csq_re = re.compile(r"(?:^|;)CSQ=([^;]+)")
    out: dict[tuple[str, str, str, str], str] = {}
    with _open_vcf(vep_vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            match = csq_re.search(cols[7])
            if not match:
                continue
            entries = [entry.split("|") for entry in match.group(1).split(",")]
            picked = next(
                (
                    values for values in entries
                    if pick_index >= 0
                    and pick_index < len(values)
                    and values[pick_index] == "1"
                ),
                entries[0],
            )
            scores: list[float] = []
            for index in indexes:
                if index >= len(picked):
                    continue
                for token in picked[index].split("&"):
                    try:
                        scores.append(abs(float(token)))
                    except ValueError:
                        pass
            if not scores:
                continue
            for alt in cols[4].split(","):
                out[(cols[0], cols[1], cols[3], alt)] = f"{max(scores):.4f}"
    return out


def merge_into_tsv(
    in_tsv: Path,
    out_tsv: Path,
    annotations: dict[tuple[str, str, str, str], str],
) -> tuple[int, int]:
    overwriting = in_tsv.resolve() == out_tsv.resolve()
    target = Path(str(out_tsv) + ".tmp") if overwriting else out_tsv
    row_count = filled_count = 0
    with in_tsv.open("r", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if TSV_COL_SPLICEAI not in fields:
            fields.append(TSV_COL_SPLICEAI)
        with target.open("w", encoding="utf-8", newline="") as destination:
            writer = csv.DictWriter(
                destination,
                fieldnames=fields,
                delimiter="\t",
                extrasaction="ignore",
                lineterminator="\n",
            )
            writer.writeheader()
            for row in reader:
                row_count += 1
                key = tuple(
                    (row.get(name) or "").strip()
                    for name in ("CHROM", "POS", "REF", "ALT")
                )
                score = annotations.get(key)
                if score and not (row.get(TSV_COL_SPLICEAI) or "").strip():
                    row[TSV_COL_SPLICEAI] = score
                    filled_count += 1
                writer.writerow(row)
    if overwriting:
        os.replace(target, out_tsv)
    return row_count, filled_count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--out-tsv")
    parser.add_argument("--workdir")
    parser.add_argument("--ref-dir", default=DEFAULT_REF_DIR)
    parser.add_argument("--vep-sif", default=DEFAULT_VEP_SIF)
    parser.add_argument("--vep-cache")
    parser.add_argument("--ref-fasta")
    parser.add_argument("--fork", type=int, default=4)
    parser.add_argument("--max-af", type=float, default=0.01)
    parser.add_argument("--af-cols", default="GNOMAD_G_AF")
    parser.add_argument("--candidate-bed")
    parser.add_argument("--spliceai-snv", default=DEFAULT_SPLICEAI_SNV)
    parser.add_argument("--spliceai-indel", default=DEFAULT_SPLICEAI_INDEL)
    args = parser.parse_args()

    args.vep_cache = args.vep_cache or f"{args.ref_dir}/tertiary/vep_cache"
    args.ref_fasta = args.ref_fasta or f"{args.ref_dir}/Homo_sapiens_assembly38.fasta"
    in_tsv = Path(args.tsv).resolve()
    out_tsv = Path(args.out_tsv).resolve() if args.out_tsv else in_tsv
    spliceai_snv = Path(args.spliceai_snv)
    spliceai_indel = Path(args.spliceai_indel)
    required = [
        in_tsv,
        spliceai_snv,
        Path(str(spliceai_snv) + ".tbi"),
        spliceai_indel,
        Path(str(spliceai_indel) + ".tbi"),
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        print("ERROR: required SpliceAI input not found: " + ", ".join(missing), file=sys.stderr)
        return 2

    candidate_bed = None
    if args.candidate_bed:
        bed_path = Path(args.candidate_bed).resolve()
        if not bed_path.is_file():
            print(f"ERROR: --candidate-bed not found: {bed_path}", file=sys.stderr)
            return 2
        candidate_bed = load_bed(bed_path)

    if args.workdir:
        workdir = Path(args.workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        temp_context = None
    else:
        temp_context = tempfile.TemporaryDirectory(prefix="spliceai-")
        workdir = Path(temp_context.name)

    try:
        sites = workdir / "sites.vcf"
        vep_vcf = workdir / "sites.vep.vcf.gz"
        af_cols = [value.strip() for value in args.af_cols.split(",") if value.strip()]
        count, dropped_af, dropped_bed = tsv_to_sites(
            in_tsv,
            sites,
            max_af=args.max_af,
            af_cols=af_cols,
            candidate_bed=candidate_bed,
        )
        print(
            f"[spliceai] {count} sites (AF filtered {dropped_af}; BED filtered {dropped_bed})",
            file=sys.stderr,
        )
        if count:
            run_vep(args, sites, vep_vcf)
            annotations = parse_vep_vcf(vep_vcf)
        else:
            annotations = {}
        rows, filled = merge_into_tsv(in_tsv, out_tsv, annotations)
        print(f"[spliceai] backfilled {filled}/{rows} rows -> {out_tsv}", file=sys.stderr)
    finally:
        if temp_context is not None:
            temp_context.cleanup()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
