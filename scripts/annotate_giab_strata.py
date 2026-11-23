#!/usr/bin/env python3
"""Annotate snv_indel.annotated.tsv with GIAB genome-stratification labels.

The GIAB / GA4GH genome stratifications are BED files that flag genomic
contexts where short-read calling is unreliable: homopolymers, tandem
repeats, segmental duplications, low-mappability regions, GC extremes,
and other difficult regions (MHC/KIR/VDJ, false duplications, ...).

This post-processing intersects every variant in the TSV against a configured
set of stratification BEDs and writes the matching short labels into a
single ``GIAB_STRATA`` column (comma-separated, e.g. ``homopolymer,segdup``).
The SNV adapter splits that column into a list and the variant card shows
one badge per label so reviewers can see at a glance when a call sits in a
difficult region.

It is fill-or-augment and idempotent: re-running overwrites the column.
Unlike the GeneBe / Extra-VEP post-processings it runs on the *whole* TSV (not
just the AF/CDS candidate set) because the labels are display-only QC
hints and interval lookup is cheap.

Strata are described by a manifest JSON in the stratification directory::

    [
      {"file": "GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
       "label": "homopolymer", "display": "Homopolymer",
       "tooltip": "GIAB: homopolymer region (low-confidence indel calls)"},
      ...
    ]

``label`` goes into the TSV; ``display`` / ``tooltip`` are carried by the
manifest only (the frontend has its own label→display map). When no
manifest is present every ``*.bed`` / ``*.bed.gz`` in the directory is
loaded and the label is derived from the filename.

Usage:
    scripts/annotate_giab_strata.py \\
        --tsv tertiary_output/<SID>/snv_indel.annotated.tsv \\
        [--strat-dir /path/to/giab_stratification] \\
        [--manifest /path/to/strata_manifest.json] \\
        [--column GIAB_STRATA]
"""
from __future__ import annotations

import argparse
import bisect
import csv
from array import array
import gzip
import json
import os
import re
import sys
import tempfile
from pathlib import Path

# Under the biotools tree like cds_combined.bed / spliceai; env / flag override.
DEFAULT_STRAT_DIR = os.environ.get(
    "NGS_UI_GIAB_STRAT_DIR",
    str(Path.home() / "NGS_UI" / "biotools" / "giab_stratification"),
)
DEFAULT_COLUMN = "GIAB_STRATA"
MANIFEST_NAME = "strata_manifest.json"

# (sorted-starts, merged-intervals) per normalised chromosome.
# chrom -> (starts, ends) of merged, sorted intervals. Plain arrays (not lists
# of tuples) keep a whole-genome stratum to a few bytes per interval and give
# forked batch workers nothing to copy-on-write through reference counts.
BedIndex = dict[str, "tuple[array, array]"]


def _norm_chrom(chrom: str) -> str:
    s = (chrom or "").strip()
    if s.lower().startswith("chr"):
        s = s[3:]
    if s in ("M", "MT", "m", "mt"):
        return "MT"
    return s.upper() if s.upper() in ("X", "Y") else s


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def load_bed(path: Path) -> BedIndex:
    """Load + merge a (optionally gzipped) BED into an interval index."""
    raw: dict[str, list[tuple[int, int]]] = {}
    with _open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if end <= start:
                continue
            raw.setdefault(_norm_chrom(parts[0]), []).append((start, end))

    out: BedIndex = {}
    for chrom, intervals in raw.items():
        intervals.sort()
        merged: list[tuple[int, int]] = []
        for start, end in intervals:
            if not merged or start > merged[-1][1]:
                merged.append((start, end))
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        out[chrom] = (array("q", (s for s, _ in merged)), array("q", (e for _, e in merged)))
    return out


def _variant_span_0based(pos: str, ref: str) -> tuple[int, int] | None:
    try:
        start = int(pos) - 1
    except (ValueError, TypeError):
        return None
    if start < 0:
        return None
    return start, start + max(1, len(ref or ""))


def _overlaps(bed: BedIndex, chrom: str, pos: str, ref: str) -> bool:
    span = _variant_span_0based(pos, ref)
    if span is None:
        return False
    start, end = span
    item = bed.get(_norm_chrom(chrom))
    if not item:
        return False
    starts, ends = item
    i = bisect.bisect_right(starts, start) - 1
    if i >= 0 and ends[i] > start:
        return True
    j = i + 1
    return j < len(starts) and starts[j] < end


def _label_from_filename(name: str) -> str:
    """Derive a short label from a GIAB BED filename when no manifest entry
    exists, e.g. ``GRCh38_AllHomopolymers_gt6bp...bed.gz`` → ``allhomopolymers``."""
    base = re.sub(r"\.bed(\.gz)?$", "", name, flags=re.IGNORECASE)
    base = re.sub(r"^GRCh3[78]_", "", base)
    base = re.sub(r"_slop\d+$", "", base)
    base = re.sub(r"[^A-Za-z0-9]+", "_", base).strip("_").lower()
    return base or name.lower()


def load_manifest(strat_dir: Path, manifest_path: Path | None) -> list[dict]:
    """Return [{file,label,...}]. Falls back to directory scan when no
    manifest file is present so a dev machine can just drop BEDs in."""
    path = manifest_path or (strat_dir / MANIFEST_NAME)
    if path.is_file():
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        if not isinstance(data, list):
            raise SystemExit(f"manifest must be a JSON list: {path}")
        entries = []
        for item in data:
            if not isinstance(item, dict) or not item.get("file"):
                continue
            label = (item.get("label") or _label_from_filename(item["file"])).strip()
            entries.append({"file": item["file"], "label": label})
        return entries

    # No manifest: scan the directory for BED files.
    entries = []
    for p in sorted(strat_dir.glob("*")):
        if p.name.lower().endswith((".bed", ".bed.gz")):
            entries.append({"file": p.name, "label": _label_from_filename(p.name)})
    return entries


def annotate(tsv: Path, strat_dir: Path, manifest_path: Path | None, column: str) -> None:
    beds = load_strata(strat_dir, manifest_path)
    if beds:
        annotate_with(tsv, beds, column)


def load_strata(strat_dir: Path, manifest_path: Path | None) -> list[tuple[str, BedIndex]]:
    """Load every stratum BED once (reused for each TSV of a batch)."""
    if not strat_dir.is_dir():
        print(f"[giab-strata] dir not found: {strat_dir} — skipping", file=sys.stderr)
        return []

    entries = load_manifest(strat_dir, manifest_path)
    if not entries:
        print(f"[giab-strata] no BEDs found under {strat_dir} — skipping", file=sys.stderr)
        return []

    # Load every BED once. Preserve manifest order for deterministic output.
    beds: list[tuple[str, BedIndex]] = []
    seen_labels: set[str] = set()
    for entry in entries:
        bed_path = strat_dir / entry["file"]
        if not bed_path.is_file():
            print(f"[giab-strata] missing BED, skipped: {bed_path}", file=sys.stderr)
            continue
        label = entry["label"]
        if label in seen_labels:
            print(f"[giab-strata] duplicate label '{label}' — keeping first", file=sys.stderr)
            continue
        seen_labels.add(label)
        beds.append((label, load_bed(bed_path)))

    if not beds:
        print(f"[giab-strata] no loadable BEDs under {strat_dir} — skipping", file=sys.stderr)
    return beds


def annotate_with(tsv: Path, beds: list[tuple[str, BedIndex]], column: str) -> None:
    """Rewrite ``tsv`` with ``column`` filled, streaming one row at a time.

    Memory does not grow with the TSV (a WGS TSV is ~6M rows; holding them as
    dicts took ~50 GB per sample). Output is identical to the former
    DictReader/DictWriter rewrite: blank lines are dropped, short rows are
    padded with empty cells and cells beyond the header are discarded.
    """
    fd, tmp_name = tempfile.mkstemp(
        dir=str(tsv.parent), prefix=tsv.name + ".", suffix=".tmp"
    )
    n_rows = n_hit = 0
    try:
        with open(tsv, "r", encoding="utf-8", newline="") as fi, \
                os.fdopen(fd, "w", encoding="utf-8", newline="") as fo:
            reader = csv.reader(fi, delimiter="\t")
            header = next(reader, None)
            if not header:
                print(f"[giab-strata] empty TSV: {tsv}", file=sys.stderr)
                fo.close()
                os.unlink(tmp_name)
                return
            width = len(header)
            out_header = list(header)
            if column in header:
                col_idx = header.index(column)
            else:
                col_idx = width
                out_header.append(column)
            ic = header.index("CHROM") if "CHROM" in header else None
            ip = header.index("POS") if "POS" in header else None
            ir = header.index("REF") if "REF" in header else None
            writer = csv.writer(fo, delimiter="\t", lineterminator="\n")
            writer.writerow(out_header)
            for row in reader:
                if not row:
                    continue
                if len(row) < width:
                    row = row + [""] * (width - len(row))
                elif len(row) > width:
                    row = row[:width]
                chrom = row[ic].strip() if ic is not None else ""
                pos = row[ip].strip() if ip is not None else ""
                ref = row[ir].strip() if ir is not None else ""
                labels = [label for label, bed in beds if _overlaps(bed, chrom, pos, ref)]
                value = ",".join(labels)
                if col_idx == width:
                    row.append(value)
                else:
                    row[col_idx] = value
                writer.writerow(row)
                n_rows += 1
                if labels:
                    n_hit += 1
        # Atomic replace so a crash can't leave a half-written TSV.
        os.replace(tmp_name, tsv)
    except BaseException:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        raise

    print(
        f"[giab-strata] {n_rows} variants, {n_hit} in >=1 stratum, "
        f"{len(beds)} BEDs: {', '.join(label for label, _ in beds)}"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True, type=Path,
                    help="snv_indel.annotated.tsv to annotate in place")
    ap.add_argument("--strat-dir", type=Path, default=Path(DEFAULT_STRAT_DIR),
                    help=f"directory of GIAB stratification BEDs (default {DEFAULT_STRAT_DIR})")
    ap.add_argument("--manifest", type=Path, default=None,
                    help=f"strata manifest JSON (default <strat-dir>/{MANIFEST_NAME})")
    ap.add_argument("--column", default=DEFAULT_COLUMN,
                    help=f"TSV column to write (default {DEFAULT_COLUMN})")
    args = ap.parse_args()

    if not args.tsv.is_file():
        raise SystemExit(f"--tsv not found: {args.tsv}")
    annotate(args.tsv, args.strat_dir, args.manifest, args.column)


if __name__ == "__main__":
    main()
