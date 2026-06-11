#!/usr/bin/env python3
"""GeneBe ACMG annotation — write a second opinion to GENEBE_* columns.

Reads `--tsv snv_indel.annotated.tsv`, looks every variant up in the
**local GeneBe database** (`genebe_hg38.tsv.gz`, a bgzip TSV) and writes
the GeneBe classification into NEW columns:
    GENEBE_ACMG_SCORE
    GENEBE_ACMG_CRITERIA
    GENEBE_ACMG_CLASS

This replaces the old live-GeneBe-API call (pygenebe via apptainer):
lookups are now fully offline, need no credentials/network, and have no
rate limit. By default the WHOLE TSV is annotated (no AF / candidate-BED
gate) — every variant present in the DB gets a second opinion, including
those only reachable via gene search. The DB is a pre-computed cache, so
a variant absent from it (e.g. a novel coding indel) simply gets no
GeneBe second opinion; the pipeline's own ACMG_CLASS still shows. There
is no API fallback by design (DB-only).

Implementation: the DB is read in a single streaming pass (the wanted
variant keys are held in a set, the DB is decompressed once and matched
against them). NGS-UI therefore does NOT use the tabix `.tbi` index at
all, so a stale/missing index can't break this step; malformed DB rows
(non-integer pos) are skipped defensively.

The pipeline's own ACMG_SCORE / ACMG_CRITERIA / ACMG_CLASS columns are
NEVER touched — the UI shows both side by side. Only the DB's
`acmg_score` / `acmg_criteria` are read (by COLUMN NAME from the '#'
header, so a slim 7-column DB and a full 55-column DB work identically);
the displayed class is derived locally from the score via classify(), so
switching from API to DB is a pure data-source swap. On a DB miss the
existing GENEBE_* cell is left as-is (re-runs keep prior values).

Score → class mapping (GeneBe acmg_score → 5-tier label):
    >= 10        Pathogenic
    6..9         Likely pathogenic
    0..5         Uncertain significance
    -1..-6       Likely benign
    <= -7        Benign

DB path (flag / env):
    --genebe-db   NGS_UI_GENEBE_DB   (default
                  $HOME/NGS_UI/biotools/genebe/genebe_hg38.tsv.gz)
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

BedIndex = dict[str, tuple[list[int], list[tuple[int, int]]]]

_MISSING = ("", ".", "NA", "N/A")


def classify(score: float | None) -> str:
    if score is None:
        return ""
    if score >= 10:  return "Pathogenic"
    if score >=  6:  return "Likely pathogenic"
    if score >=  0:  return "Uncertain significance"
    if score >= -6:  return "Likely benign"
    return "Benign"


def _to_float(s: str) -> float | None:
    s = (s or "").strip()
    if s in _MISSING:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _row_max_af(row: dict, cols: list[str]) -> float:
    """Largest numeric AF across the listed columns; missing → 0."""
    m = 0.0
    for c in cols:
        v = row.get(c)
        if v is None:
            continue
        s = str(v).strip()
        if not s or s.upper() in ("NA", "N/A", "."):
            continue
        try:
            f = float(s)
        except ValueError:
            continue
        if f > m:
            m = f
    return m


def _norm_chrom(chrom: str) -> str:
    s = chrom.strip()
    if s.lower().startswith("chr"):
        s = s[3:]
    if s in ("M", "MT", "m", "mt"):
        return "MT"
    return s.upper() if s.upper() in ("X", "Y") else s


def _chr_prefixed(chrom: str) -> str:
    """DB rows are chr-prefixed (chr1…chrX, chrM). Normalise a TSV CHROM
    to that form so keys match regardless of the TSV's prefix style."""
    s = (chrom or "").strip()
    if not s:
        return s
    body = s[3:] if s.lower().startswith("chr") else s
    if body in ("MT", "mt", "M", "m"):
        body = "M"
    return "chr" + body


def _vkey(chrom: str, pos: str, ref: str, alt: str) -> str:
    """Canonical (chr-prefixed) variant key used to match TSV ↔ DB."""
    return f"{_chr_prefixed(chrom)}:{pos}:{ref}:{alt}"


def load_bed(path: Path) -> BedIndex:
    """Load and merge a BED file into a small interval index."""
    raw: dict[str, list[tuple[int, int]]] = {}
    with open(path, "r", encoding="utf-8") as fh:
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
        out[chrom] = ([start for start, _ in merged], merged)
    return out


def _variant_span_0based(pos: str, ref: str) -> tuple[int, int] | None:
    try:
        start = int(pos) - 1
    except ValueError:
        return None
    if start < 0:
        return None
    return start, start + max(1, len(ref or ""))


def _overlaps_bed(bed: BedIndex | None, chrom: str, pos: str, ref: str) -> bool:
    if bed is None:
        return True
    span = _variant_span_0based(pos, ref)
    if span is None:
        return False
    start, end = span
    item = bed.get(_norm_chrom(chrom))
    if not item:
        return False
    starts, intervals = item
    i = bisect.bisect_right(starts, start) - 1
    if i >= 0 and intervals[i][1] > start:
        return True
    j = i + 1
    return j < len(intervals) and intervals[j][0] < end


def collect_wanted(
    tsv_in: Path,
    *,
    max_af: float | None,
    af_cols: list[str],
    candidate_bed: BedIndex | None,
) -> tuple[set[str], int, int, int]:
    """Set of canonical variant keys to look up.

    Whole-TSV by default (max_af=None, candidate_bed=None). The optional
    AF / BED gate is kept for emergencies / WES speed. Returns
    (wanted, n_dropped_by_af, n_skipped_star, n_dropped_by_bed).
    """
    wanted: set[str] = set()
    n_af = n_star = n_bed = 0
    with open(tsv_in, "r", encoding="utf-8", newline="") as fi:
        for row in csv.DictReader(fi, delimiter="\t"):
            chrom = (row.get("CHROM") or "").strip()
            pos   = (row.get("POS")   or "").strip()
            ref   = (row.get("REF")   or "").strip()
            alt   = (row.get("ALT")   or "").strip()
            if not (chrom and pos and ref and alt):
                continue
            if "*" in (ref, alt):
                n_star += 1
                continue
            if max_af is not None and _row_max_af(row, af_cols) > max_af:
                n_af += 1
                continue
            if not _overlaps_bed(candidate_bed, chrom, pos, ref):
                n_bed += 1
                continue
            wanted.add(_vkey(chrom, pos, ref, alt))
    return wanted, n_af, n_star, n_bed


# ---- local GeneBe DB (bgzip TSV, streamed) --------------------------

class GeneBeDBError(RuntimeError):
    pass


def read_db_header(db: Path) -> dict[str, int]:
    """Read the '#'-prefixed header (first line) → {col_name: index}.

    Uses Python gzip so it only touches the first block — cheap, and
    works on bgzip files (bgzip is gzip-compatible).
    """
    try:
        with gzip.open(db, "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#"):
                    cols = line.rstrip("\n").lstrip("#").split("\t")
                    return {name: i for i, name in enumerate(cols)}
                break
    except OSError as e:
        raise GeneBeDBError(f"cannot read DB header: {e}") from e
    raise GeneBeDBError("DB has no '#'-prefixed header line")


def _db_line_stream(db: Path):
    """Yield decompressed DB lines + return the subprocess (or None).

    Prefer a C decompressor (bgzip/zcat/gzip) for speed; fall back to
    Python gzip. Returns (proc_or_None, iterable_of_lines).
    """
    for exe in (["bgzip", "-dc"], ["zcat"], ["gzip", "-dc"]):
        if shutil.which(exe[0]):
            proc = subprocess.Popen(
                exe + [str(db)], stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, text=True, bufsize=1 << 20,
            )
            return proc, proc.stdout
    return None, gzip.open(db, "rt", encoding="utf-8")


def index_db_rows(lines, idx: dict[str, int], wanted: set[str]) -> dict[str, tuple[str, str]]:
    """Single pass over DB lines → {vkey: (acmg_score, acmg_criteria)} for
    the wanted keys. Pure (no subprocess) so it is unit-testable.

    Defensive: skips the header, short rows, and rows whose pos isn't an
    integer (the DB can carry malformed '.' placeholder rows).
    """
    ci_chr, ci_pos = idx["chr"], idx["pos"]
    ci_ref, ci_alt = idx["ref"], idx["alt"]
    ci_score, ci_crit = idx["acmg_score"], idx["acmg_criteria"]
    need = max(ci_chr, ci_pos, ci_ref, ci_alt, ci_score, ci_crit)
    hits: dict[str, tuple[str, str]] = {}
    for line in lines:
        if not line or line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) <= need:
            continue
        pos = f[ci_pos]
        if not pos.isdigit():
            continue
        key = f"{f[ci_chr]}:{pos}:{f[ci_ref]}:{f[ci_alt]}"
        if key in wanted:
            hits[key] = (f[ci_score], f[ci_crit])
    return hits


def scan_db(db: Path, idx: dict[str, int], wanted: set[str]) -> dict[str, tuple[str, str]]:
    """Stream the whole DB once and collect hits for the wanted keys."""
    proc, lines = _db_line_stream(db)
    try:
        hits = index_db_rows(lines, idx, wanted)
    finally:
        if proc is not None:
            err = proc.stderr.read() if proc.stderr else ""
            proc.wait()
            if proc.returncode not in (0, None):
                raise GeneBeDBError(
                    f"DB decompression failed (rc={proc.returncode}): "
                    f"{(err or '').strip()}")
        else:
            lines.close()
    return hits


def assemble_gb(db_hits: dict[str, tuple[str, str]]) -> dict[str, tuple[str, str, str]]:
    """{vkey: (score, criteria)} → {vkey: (score, criteria, class)}."""
    gb: dict[str, tuple[str, str, str]] = {}
    for key, (score_raw, crit_raw) in db_hits.items():
        score_raw = (score_raw or "").strip()
        crit_raw = (crit_raw or "").strip()
        score_out = "" if score_raw in _MISSING else score_raw
        crit_out = "" if crit_raw in _MISSING else crit_raw
        cls = classify(_to_float(score_raw))
        if score_out or crit_out or cls:
            gb[key] = (score_out, crit_out, cls)
    return gb


def preflight_db(db: Path) -> dict[str, int]:
    """Validate the DB is readable and has the required columns."""
    if not db.is_file():
        raise GeneBeDBError(f"GeneBe DB not found: {db}")
    idx = read_db_header(db)
    for col in ("chr", "pos", "ref", "alt", "acmg_score", "acmg_criteria"):
        if col not in idx:
            raise GeneBeDBError(
                f"DB header missing required column '{col}'; got {list(idx)}")
    return idx


# ---- merge back into the TSV ----------------------------------------

def merge_into_tsv(in_tsv: Path, out_tsv: Path, gb: dict) -> tuple[int, int]:
    """Write a new TSV with GENEBE_* backfilled. Returns (n_filled, n_total).

    On a DB hit the GENEBE_* cells are (over)written; on a miss the row's
    existing GENEBE_* values are left untouched. Pipeline ACMG_* is never
    touched.
    """
    in_tsv = Path(in_tsv)
    out_tsv = Path(out_tsv)
    overwriting = in_tsv.resolve() == out_tsv.resolve()
    target = Path(str(out_tsv) + ".tmp") if overwriting else out_tsv

    n_filled = 0
    n_total = 0
    with open(in_tsv, "r", encoding="utf-8", newline="") as fi:
        reader = csv.DictReader(fi, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for col in ("GENEBE_ACMG_SCORE", "GENEBE_ACMG_CRITERIA",
                    "GENEBE_ACMG_CLASS"):
            if col not in fieldnames:
                fieldnames.append(col)
        with open(target, "w", encoding="utf-8", newline="") as fo:
            writer = csv.DictWriter(fo, fieldnames=fieldnames, delimiter="\t",
                                    extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                n_total += 1
                key = _vkey((row.get("CHROM") or "").strip(),
                            (row.get("POS")   or "").strip(),
                            (row.get("REF")   or "").strip(),
                            (row.get("ALT")   or "").strip())
                gbt = gb.get(key)
                if gbt:
                    score, crit, cls = gbt
                    if score: row["GENEBE_ACMG_SCORE"]    = score
                    if crit:  row["GENEBE_ACMG_CRITERIA"] = crit
                    if cls:   row["GENEBE_ACMG_CLASS"]    = cls
                    if score or crit or cls:
                        n_filled += 1
                writer.writerow(row)
    if overwriting:
        os.replace(target, out_tsv)
    return n_filled, n_total


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--tsv", required=True,
                    help="snv_indel.annotated.tsv (updated in place unless "
                         "--out-tsv given)")
    ap.add_argument("--out-tsv",
                    help="write merged TSV here instead of overwriting --tsv")
    ap.add_argument("--genebe-db",
                    default=os.environ.get(
                        "NGS_UI_GENEBE_DB",
                        str(Path.home() / "NGS_UI" / "biotools" / "genebe"
                            / "genebe_hg38.tsv.gz")),
                    help="local GeneBe DB (bgzip TSV); env NGS_UI_GENEBE_DB")
    ap.add_argument("--max-af", type=float, default=-1.0,
                    help="optional: drop sites whose AF > this before lookup "
                         "(default -1 = whole TSV, no AF gate)")
    ap.add_argument("--af-cols", default="GNOMAD_G_AF",
                    help="comma-separated AF columns for --max-af")
    ap.add_argument("--candidate-bed",
                    help="optional: restrict lookup to these BED regions "
                         "(default: whole TSV)")
    args = ap.parse_args()

    in_tsv = Path(args.tsv).resolve()
    if not in_tsv.is_file():
        print(f"ERROR: --tsv 找不到：{in_tsv}", file=sys.stderr)
        return 2
    out_tsv = Path(args.out_tsv).resolve() if args.out_tsv else in_tsv

    db = Path(args.genebe_db).resolve()
    try:
        idx = preflight_db(db)
    except GeneBeDBError as e:
        print(f"ERROR: GeneBe DB 不可用：{e}", file=sys.stderr)
        return 2
    print(f"[genebe] DB: {db}", file=sys.stderr)

    candidate_bed = None
    if args.candidate_bed:
        bed_path = Path(args.candidate_bed).resolve()
        if not bed_path.is_file():
            print(f"ERROR: --candidate-bed 找不到：{bed_path}", file=sys.stderr)
            return 2
        candidate_bed = load_bed(bed_path)

    max_af = None if args.max_af < 0 else args.max_af
    af_cols = [c.strip() for c in args.af_cols.split(",") if c.strip()]
    wanted, n_af, n_star, n_bed = collect_wanted(
        in_tsv, max_af=max_af, af_cols=af_cols, candidate_bed=candidate_bed,
    )
    scope = "whole TSV" if (max_af is None and candidate_bed is None) else "gated"
    print(f"[genebe] {len(wanted)} unique variants to look up ({scope}; "
          f"AF-dropped {n_af}, BED-dropped {n_bed}, '*'-skipped {n_star})",
          file=sys.stderr)
    if not wanted:
        merge_into_tsv(in_tsv, out_tsv, {})
        print("[genebe] nothing to look up", file=sys.stderr)
        return 0

    t0 = time.time()
    try:
        db_hits = scan_db(db, idx, wanted)
    except GeneBeDBError as e:
        print(f"ERROR: GeneBe DB 掃描失敗：{e}", file=sys.stderr)
        return 2
    gb = assemble_gb(db_hits)
    print(f"[genebe] scanned DB in {time.time() - t0:.0f}s, "
          f"{len(gb)} variants matched", file=sys.stderr)

    n_filled, n_total = merge_into_tsv(in_tsv, out_tsv, gb)
    print(f"[genebe] backfilled ACMG for {n_filled}/{n_total} TSV rows",
          file=sys.stderr)
    print(f"[genebe] done → {out_tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
