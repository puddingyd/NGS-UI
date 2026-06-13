#!/usr/bin/env python3
"""Compare GeneBe DB SpliceAI coverage against the current extra-VEP path.

This is a read-only diagnostic for deciding whether the local
genebe_hg38.tsv.gz can replace the post-processing extra-VEP SpliceAI step.
It samples variants from the GeneBe bgzip TSV, runs VEP with the same
SpliceAI plugin/data files used by annotate_extra_vep.py, and reports
where GeneBe has/misses a SpliceAI max score compared with VEP.

Typical server run:

    scripts/compare_genebe_spliceai_coverage.py \
      --genebe-db /home/n102968/NGS_UI/biotools/genebe/genebe_hg38.tsv.gz \
      --max-sites 100000 \
      --sample-mode chrom-balanced \
      --out-dir /home/n102968/NGS_UI/tmp/genebe_spliceai_check

Use --max-sites 0 only when you intentionally want to run every row in
the selected region/chromosome through VEP.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import random
import re
import secrets
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_GENEBE_DB = "/home/n102968/NGS_UI/biotools/genebe/genebe_hg38.tsv.gz"
DEFAULT_REF_DIR = "/home/pipeline/reference/hg38"
DEFAULT_VEP_SIF = "/home/pipeline/nextflow_containers/vep_115.sif"
CANONICAL_CHROMS = [*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY"]


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    consequence: str
    genebe_splice_score: str
    genebe_splice_pred: str
    genebe_splice_source: str

    @property
    def key(self) -> tuple[str, str, str, str]:
        return (self.chrom, str(self.pos), self.ref, self.alt)

    @property
    def kind(self) -> str:
        if len(self.ref) == 1 and len(self.alt) == 1:
            return "snv"
        if self.ref and self.alt and self.alt != ".":
            return "indel"
        return "other"


def _expand_path(raw: str | None) -> Path | None:
    if raw is None:
        return None
    return Path(os.path.expandvars(os.path.expanduser(raw))).resolve()


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="")


def _is_missing(value: object) -> bool:
    s = str(value or "").strip()
    return not s or s.upper() in {".", "NA", "N/A", "NONE", "NULL"}


def _to_float(value: object) -> float | None:
    if _is_missing(value):
        return None
    try:
        f = float(str(value).strip())
    except ValueError:
        return None
    return f if math.isfinite(f) else None


def _fmt_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6g}"


def _parse_region(region: str | None) -> tuple[str, int | None, int | None] | None:
    if not region:
        return None
    chrom, sep, rest = region.partition(":")
    if not sep:
        return chrom, None, None
    start_s, sep2, end_s = rest.replace(",", "").partition("-")
    if not sep2:
        pos = int(start_s)
        return chrom, pos, pos
    return chrom, int(start_s), int(end_s)


def _variant_passes(v: Variant, args) -> bool:
    if args.chrom and v.chrom not in args.chrom:
        return False
    if args.region:
        chrom, start, end = args.region
        if v.chrom != chrom:
            return False
        if start is not None and v.pos < start:
            return False
        if end is not None and v.pos > end:
            return False
    if args.variant_kind != "all" and v.kind != args.variant_kind:
        return False
    has_gb = _to_float(v.genebe_splice_score) is not None
    if args.only == "genebe-has-spliceai" and not has_gb:
        return False
    if args.only == "genebe-missing-spliceai" and has_gb:
        return False
    if not args.include_mito and v.chrom in {"chrM", "chrMT", "MT", "M"}:
        return False
    return True


def read_genebe_variants(path: Path, args) -> list[Variant]:
    if not path.is_file():
        raise FileNotFoundError(f"GeneBe DB not found: {path}")

    max_sites = None if args.max_sites == 0 else args.max_sites
    variants: list[Variant] = []
    chrom_buckets: dict[str, list[Variant]] = {}
    chrom_seen: dict[str, int] = {}
    seen: set[tuple[str, str, str, str]] = set()
    selected_seen = 0
    reservoir_seen = 0
    rng = random.Random(args.seed)
    chrom_cap = None
    if max_sites is not None and args.sample_mode == "chrom-balanced":
        chrom_count = len(args.chrom or CANONICAL_CHROMS)
        chrom_cap = max(1, math.ceil(max_sites / chrom_count))

    with _open_text(path) as fh:
        header: list[str] | None = None
        for line in fh:
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                break
        if not header:
            raise ValueError(f"No # header found in {path}")
        idx = {name: i for i, name in enumerate(header)}
        required = {"chr", "pos", "ref", "alt", "spliceai_max_score"}
        missing = sorted(required - set(idx))
        if missing:
            raise ValueError(f"GeneBe DB missing required columns: {', '.join(missing)}")

        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            try:
                chrom = fields[idx["chr"]].strip()
                pos = int(fields[idx["pos"]])
                ref = fields[idx["ref"]].strip().upper()
                alt = fields[idx["alt"]].strip().upper()
            except (IndexError, ValueError):
                continue
            if not chrom or not ref or not alt or "*" in (ref, alt):
                continue
            v = Variant(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                gene=_field(fields, idx, "gene_symbol"),
                consequence=_field(fields, idx, "consequences"),
                genebe_splice_score=_field(fields, idx, "spliceai_max_score"),
                genebe_splice_pred=_field(fields, idx, "spliceai_max_prediction"),
                genebe_splice_source=_field(fields, idx, "splice_source_selected"),
            )
            if not _variant_passes(v, args):
                continue
            selected_seen += 1
            if args.step > 1 and (selected_seen - 1) % args.step:
                continue
            if v.key in seen:
                continue
            seen.add(v.key)
            if max_sites is None:
                variants.append(v)
                continue
            if args.sample_mode == "chrom-balanced":
                bucket = chrom_buckets.setdefault(v.chrom, [])
                chrom_seen[v.chrom] = chrom_seen.get(v.chrom, 0) + 1
                assert chrom_cap is not None
                if len(bucket) < chrom_cap:
                    bucket.append(v)
                else:
                    j = rng.randrange(chrom_seen[v.chrom])
                    if j < chrom_cap:
                        bucket[j] = v
            elif args.sample_mode == "random":
                reservoir_seen += 1
                if len(variants) < max_sites:
                    variants.append(v)
                else:
                    j = rng.randrange(reservoir_seen)
                    if j < max_sites:
                        variants[j] = v
            else:
                variants.append(v)
                if len(variants) >= max_sites:
                    break
    if max_sites is not None and args.sample_mode == "chrom-balanced":
        for chrom in sorted(chrom_buckets, key=_chrom_sort_key):
            variants.extend(chrom_buckets[chrom])
        variants = variants[:max_sites]
    return variants


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if c.isdigit():
        return int(c), ""
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return order.get(c.upper(), 99), c


def _field(fields: list[str], idx: dict[str, int], name: str) -> str:
    i = idx.get(name)
    if i is None or i >= len(fields):
        return ""
    value = fields[i].strip()
    return "" if value == "." else value


def write_sites_vcf(variants: list[Variant], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for contig in [*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY", "chrM", "chrMT"]:
            fh.write(f"##contig=<ID={contig}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in sorted(variants, key=lambda x: (x.chrom, x.pos, x.ref, x.alt)):
            fh.write(f"{v.chrom}\t{v.pos}\t.\t{v.ref}\t{v.alt}\t.\t.\t.\n")


def run_vep_spliceai(args, sites_vcf: Path, out_vcf: Path) -> None:
    for label, path in (
        ("--vep-sif", args.vep_sif),
        ("--ref-fasta", args.ref_fasta),
        ("--spliceai-snv", args.spliceai_snv),
        ("--spliceai-indel", args.spliceai_indel),
    ):
        if not path or not Path(path).exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    binds = [
        str(Path(args.ref_dir).resolve()),
        str(sites_vcf.parent.resolve()),
        str(Path(args.spliceai_snv).resolve().parent),
        str(Path(args.spliceai_indel).resolve().parent),
    ]
    seen: set[str] = set()
    binds = [b for b in binds if not (b in seen or seen.add(b))]
    cmd = [
        "apptainer", "exec",
        "--bind", ",".join(binds),
        str(args.vep_sif),
        "vep",
        "--input_file", str(sites_vcf),
        "--output_file", str(out_vcf),
        "--vcf",
        "--compress_output", "bgzip",
        "--offline", "--cache",
        "--dir_cache", str(args.vep_cache),
        "--dir_plugins", "/opt/vep/Plugins",
        "--assembly", "GRCh38",
        "--fasta", str(args.ref_fasta),
        "--fork", str(args.fork),
        "--symbol", "--canonical", "--biotype",
        "--mane", "--flag_pick",
        "--pick_order", "mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length",
        "--force_overwrite", "--no_stats", "--safe",
        "--plugin", f"SpliceAI,snv={args.spliceai_snv},indel={args.spliceai_indel}",
    ]
    print(f"[compare] running VEP SpliceAI on {sites_vcf}", file=sys.stderr)
    proc = subprocess.run(cmd, check=False, stderr=subprocess.PIPE, text=True)
    if proc.stderr:
        skip_re = re.compile(
            r"(Use of uninitialized value \$alt_allele in string eq at /plugins/SpliceAI\.pm|"
            r"^WARNING: \d+ : Use of uninitialized value \$alt_allele)"
        )
        for line in proc.stderr.splitlines():
            if not skip_re.search(line):
                print(line, file=sys.stderr)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)


def _open_vcf(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def _parse_csq_format(vcf_path: Path) -> list[str]:
    pat = re.compile(r'##INFO=<ID=CSQ,.*Description="[^"]*Format:\s*([^"]+)"', re.IGNORECASE)
    with _open_vcf(vcf_path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            m = pat.match(line.rstrip("\n"))
            if m:
                return m.group(1).split("|")
    return []


def parse_vep_spliceai(vep_vcf: Path) -> dict[tuple[str, str, str, str], str]:
    fmt = _parse_csq_format(vep_vcf)
    if not fmt:
        return {}
    ds_fields = (
        "SpliceAI_pred_DS_AG",
        "SpliceAI_pred_DS_AL",
        "SpliceAI_pred_DS_DG",
        "SpliceAI_pred_DS_DL",
    )
    ds_idx = [fmt.index(name) for name in ds_fields if name in fmt]
    pick_idx = fmt.index("PICK") if "PICK" in fmt else -1
    if not ds_idx:
        return {}

    out: dict[tuple[str, str, str, str], str] = {}
    info_csq_pat = re.compile(r"(?:^|;)CSQ=([^;]+)")
    with _open_vcf(vep_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _, ref, alt = cols[:5]
            match = info_csq_pat.search(cols[7])
            if not match:
                continue
            entries = match.group(1).split(",")
            picked = None
            for entry in entries:
                vals = entry.split("|")
                if pick_idx >= 0 and pick_idx < len(vals) and vals[pick_idx] == "1":
                    picked = vals
                    break
            if picked is None:
                picked = entries[0].split("|")
            vals = []
            for i in ds_idx:
                if i < len(picked):
                    f = _to_float(picked[i])
                    if f is not None:
                        vals.append(abs(f))
            if vals:
                for alt_one in alt.split(","):
                    out[(chrom, pos, ref, alt_one)] = f"{max(vals):.4f}"
    return out


def compare_and_write(
    variants: list[Variant],
    vep_scores: dict[tuple[str, str, str, str], str],
    out_dir: Path,
    tolerance: float,
    write_all: bool,
) -> None:
    summary: dict[str, int] = {
        "total_variants": 0,
        "genebe_has_spliceai": 0,
        "extra_vep_has_spliceai": 0,
        "both_have_spliceai": 0,
        "neither_has_spliceai": 0,
        "genebe_only": 0,
        "extra_vep_only": 0,
        "both_equal_within_tolerance": 0,
        "both_different": 0,
    }
    by_kind: dict[str, dict[str, int]] = {}
    mismatch_rows: list[dict[str, str]] = []
    all_rows: list[dict[str, str]] = []

    for v in variants:
        gb = _to_float(v.genebe_splice_score)
        ev = _to_float(vep_scores.get(v.key))
        has_gb = gb is not None
        has_ev = ev is not None
        diff = None if gb is None or ev is None else ev - gb

        summary["total_variants"] += 1
        kind_counts = by_kind.setdefault(v.kind, {k: 0 for k in summary})
        kind_counts["total_variants"] += 1
        for metric, condition in (
            ("genebe_has_spliceai", has_gb),
            ("extra_vep_has_spliceai", has_ev),
            ("both_have_spliceai", has_gb and has_ev),
            ("neither_has_spliceai", not has_gb and not has_ev),
            ("genebe_only", has_gb and not has_ev),
            ("extra_vep_only", has_ev and not has_gb),
            ("both_equal_within_tolerance", diff is not None and abs(diff) <= tolerance),
            ("both_different", diff is not None and abs(diff) > tolerance),
        ):
            if condition:
                summary[metric] += 1
                kind_counts[metric] += 1

        row = {
            "chr": v.chrom,
            "pos": str(v.pos),
            "ref": v.ref,
            "alt": v.alt,
            "kind": v.kind,
            "gene": v.gene,
            "consequences": v.consequence,
            "genebe_spliceai_max_score": _fmt_float(gb),
            "extra_vep_spliceai_max_score": _fmt_float(ev),
            "delta_extra_minus_genebe": _fmt_float(diff),
            "genebe_splice_prediction": v.genebe_splice_pred,
            "genebe_splice_source": v.genebe_splice_source,
            "category": _category(has_gb, has_ev, diff, tolerance),
        }
        if row["category"] != "both_equal_or_both_missing":
            mismatch_rows.append(row)
        if write_all:
            all_rows.append(row)

    _write_metric_tsv(out_dir / "summary.tsv", summary)
    _write_by_kind_tsv(out_dir / "summary_by_kind.tsv", by_kind)
    _write_rows(out_dir / "mismatches.tsv", mismatch_rows)
    if write_all:
        _write_rows(out_dir / "all_variants.tsv", all_rows)


def _category(has_gb: bool, has_ev: bool, diff: float | None, tolerance: float) -> str:
    if has_gb and has_ev:
        return "both_different" if diff is not None and abs(diff) > tolerance else "both_equal_or_both_missing"
    if has_gb:
        return "genebe_only"
    if has_ev:
        return "extra_vep_only"
    return "both_equal_or_both_missing"


def _write_metric_tsv(path: Path, metrics: dict[str, int]) -> None:
    total = metrics.get("total_variants", 0) or 0
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "count", "fraction_of_total"])
        for key, value in metrics.items():
            frac = (value / total) if total else 0
            writer.writerow([key, value, f"{frac:.6f}"])


def _write_by_kind_tsv(path: Path, by_kind: dict[str, dict[str, int]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["kind", "metric", "count", "fraction_of_kind"])
        for kind, metrics in sorted(by_kind.items()):
            total = metrics.get("total_variants", 0) or 0
            for key, value in metrics.items():
                frac = (value / total) if total else 0
                writer.writerow([kind, key, value, f"{frac:.6f}"])


def _write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "category", "chr", "pos", "ref", "alt", "kind", "gene",
        "consequences", "genebe_spliceai_max_score",
        "extra_vep_spliceai_max_score", "delta_extra_minus_genebe",
        "genebe_splice_prediction", "genebe_splice_source",
    ]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_run_metadata(args, out_dir: Path) -> None:
    items = {
        "genebe_db": str(args.genebe_db),
        "max_sites": str(args.max_sites),
        "sample_mode": args.sample_mode,
        "seed": "" if args.seed is None else str(args.seed),
        "step": str(args.step),
        "chrom": ",".join(args.chrom or []),
        "region": "" if args.region is None else ":".join(
            [args.region[0], "-".join("" if x is None else str(x) for x in args.region[1:])]
        ),
        "variant_kind": args.variant_kind,
        "only": args.only,
        "include_mito": str(bool(args.include_mito)),
        "tolerance": str(args.tolerance),
        "spliceai_snv": str(args.spliceai_snv),
        "spliceai_indel": str(args.spliceai_indel),
        "vep_sif": str(args.vep_sif),
        "ref_dir": str(args.ref_dir),
    }
    with open(out_dir / "run_metadata.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        for key, value in items.items():
            writer.writerow([key, value])


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--genebe-db", default=DEFAULT_GENEBE_DB,
                    help=f"default {DEFAULT_GENEBE_DB}")
    ap.add_argument("--out-dir",
                    help="output directory (default: ./genebe_spliceai_compare_<UTC timestamp>)")
    ap.add_argument("--max-sites", type=int, default=100000,
                    help="maximum selected variants to send to VEP; 0 means no limit (default 100000)")
    ap.add_argument("--step", type=int, default=1,
                    help="deterministically take every Nth selected GeneBe row (default 1)")
    ap.add_argument("--sample-mode", choices=("first", "random", "chrom-balanced"),
                    default="chrom-balanced",
                    help="how to select up to --max-sites variants (default chrom-balanced)")
    ap.add_argument("--random-sample", action="store_true",
                    help="backward-compatible alias for --sample-mode random")
    ap.add_argument("--seed", type=int,
                    help="random seed for random/chrom-balanced sampling; default is generated from system entropy")
    ap.add_argument("--chrom", action="append",
                    help="restrict to a chromosome as stored in DB, e.g. chr1; may repeat")
    ap.add_argument("--region",
                    help="restrict to one region, e.g. chr7:140700000-140800000")
    ap.add_argument("--variant-kind", choices=("all", "snv", "indel"), default="all")
    ap.add_argument("--only",
                    choices=("all", "genebe-has-spliceai", "genebe-missing-spliceai"),
                    default="all")
    ap.add_argument("--include-mito", action="store_true",
                    help="include chrM/MT rows (SpliceAI usually has no mitochondrial scores)")
    ap.add_argument("--tolerance", type=float, default=0.0001,
                    help="numeric tolerance when comparing max scores (default 0.0001)")
    ap.add_argument("--write-all", action="store_true",
                    help="also write all_variants.tsv, not only mismatches.tsv")
    ap.add_argument("--sites-vcf",
                    help="reuse an existing sites VCF instead of rebuilding from GeneBe")
    ap.add_argument("--vep-vcf",
                    help="reuse an existing VEP VCF instead of running VEP")
    ap.add_argument("--ref-dir", default=DEFAULT_REF_DIR)
    ap.add_argument("--vep-sif", default=DEFAULT_VEP_SIF)
    ap.add_argument("--vep-cache", default=None,
                    help="default $REF_DIR/tertiary/vep_cache")
    ap.add_argument("--ref-fasta", default=None,
                    help="default $REF_DIR/Homo_sapiens_assembly38.fasta")
    ap.add_argument("--spliceai-snv",
                    default="$HOME/NGS_UI/biotools/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz")
    ap.add_argument("--spliceai-indel",
                    default="$HOME/NGS_UI/biotools/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz")
    ap.add_argument("--fork", type=int, default=4)
    args = ap.parse_args()

    if args.step < 1:
        print("ERROR: --step must be >= 1", file=sys.stderr)
        return 2
    if args.random_sample:
        args.sample_mode = "random"
    if args.seed is None and args.sample_mode != "first":
        args.seed = secrets.randbits(63)
    args.region = _parse_region(args.region)
    args.genebe_db = _expand_path(args.genebe_db)
    args.ref_dir = _expand_path(args.ref_dir)
    args.vep_sif = _expand_path(args.vep_sif)
    args.vep_cache = _expand_path(args.vep_cache) if args.vep_cache else args.ref_dir / "tertiary/vep_cache"
    args.ref_fasta = _expand_path(args.ref_fasta) if args.ref_fasta else args.ref_dir / "Homo_sapiens_assembly38.fasta"
    args.spliceai_snv = _expand_path(args.spliceai_snv)
    args.spliceai_indel = _expand_path(args.spliceai_indel)

    if args.out_dir:
        out_dir = _expand_path(args.out_dir)
    else:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        out_dir = Path.cwd() / f"genebe_spliceai_compare_{stamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    write_run_metadata(args, out_dir)

    print(f"[compare] reading GeneBe variants from {args.genebe_db}", file=sys.stderr)
    if args.sample_mode != "first":
        print(f"[compare] sample_mode={args.sample_mode} seed={args.seed}", file=sys.stderr)
    variants = read_genebe_variants(args.genebe_db, args)
    if not variants:
        print("ERROR: no variants selected from GeneBe DB", file=sys.stderr)
        return 1
    print(f"[compare] selected {len(variants)} variants", file=sys.stderr)

    sites_vcf = _expand_path(args.sites_vcf) if args.sites_vcf else out_dir / "genebe_sites.vcf"
    if args.sites_vcf:
        print(f"[compare] reusing sites VCF: {sites_vcf}", file=sys.stderr)
    else:
        write_sites_vcf(variants, sites_vcf)
        print(f"[compare] wrote {sites_vcf}", file=sys.stderr)

    vep_vcf = _expand_path(args.vep_vcf) if args.vep_vcf else out_dir / "genebe_sites.vep.vcf.gz"
    if args.vep_vcf:
        print(f"[compare] reusing VEP VCF: {vep_vcf}", file=sys.stderr)
    else:
        run_vep_spliceai(args, sites_vcf, vep_vcf)

    vep_scores = parse_vep_spliceai(vep_vcf)
    print(f"[compare] parsed SpliceAI scores for {len(vep_scores)} VEP sites", file=sys.stderr)
    compare_and_write(variants, vep_scores, out_dir, args.tolerance, args.write_all)
    print(f"[compare] done: {out_dir}", file=sys.stderr)
    print(out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
