#!/usr/bin/env python3
"""GeneBe ACMG annotation — write a second opinion to GENEBE_* columns.

Reads `--tsv snv_indel.annotated.tsv`, selects candidate sites
(skipping spanning-deletion '*' alleles, high-AF variants per --max-af,
and anything outside --candidate-bed), looks each up in the **local
GeneBe database** (`genebe_hg38.tsv.gz`, a bgzip + tabix TSV) and writes
the GeneBe classification into NEW columns named:
    GENEBE_ACMG_SCORE
    GENEBE_ACMG_CRITERIA
    GENEBE_ACMG_CLASS

This replaces the old live-GeneBe-API call (pygenebe via apptainer):
lookups are now fully offline, need no credentials/network, and have no
rate limit. The local DB is a pre-computed cache, so a variant absent
from it (e.g. a novel coding indel) simply gets no GeneBe second
opinion — the pipeline's own ACMG_CLASS still shows. There is no API
fallback by design (DB-only).

The pipeline's own ACMG_SCORE / ACMG_CRITERIA / ACMG_CLASS columns are
NEVER touched — the UI shows both side by side. Only the GeneBe DB's
`acmg_score` / `acmg_criteria` are read; the displayed class is derived
locally from the score via classify() exactly as before, so switching
from API to DB is a pure data-source swap with no display change for
variants present in both.

The DB is queried by COLUMN NAME (read from the '#'-prefixed header),
not by fixed column number, so a slim 7-column DB
(`#chr pos ref alt acmg_classification acmg_score acmg_criteria`) and a
full 55-column DB work identically.

By default updates --tsv in place; use --out-tsv to write elsewhere.

Score → class mapping (GeneBe acmg_score → 5-tier label):
    >= 10        Pathogenic
    6..9         Likely pathogenic
    0..5         Uncertain significance
    -1..-6       Likely benign
    <= -7        Benign

DB path (flag / env):
    --genebe-db   NGS_UI_GENEBE_DB   (default
                  $HOME/NGS_UI/biotools/genebe/genebe_hg38.tsv.gz)
Requires the `tabix` binary on PATH and the DB's `.tbi` index sitting
next to the `.gz` and newer than it.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import os
import shutil
import subprocess
import sys
import tempfile
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
    """Largest numeric AF across the listed columns; missing → 0.

    Matches ACMG BA1/BS1 practice: a missing AF in gnomAD isn't
    evidence of commonness, so it's kept.
    """
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
    to that form for tabix regions + key matching."""
    s = (chrom or "").strip()
    if not s:
        return s
    if s.lower().startswith("chr"):
        body = s[3:]
    else:
        body = s
    if body in ("MT", "mt", "M", "m"):
        body = "M"
    return "chr" + body


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


def collect_candidates(
    tsv_in: Path,
    *,
    max_af: float | None,
    af_cols: list[str],
    candidate_bed: BedIndex | None,
) -> tuple[list[tuple[str, str, str, str]], int, int, int]:
    """Pick the variants to look up. Same gate as the old VCF builder.

    Returns (candidates, n_dropped_by_af, n_skipped_star, n_dropped_by_bed)
    where each candidate is the TSV's own (CHROM, POS, REF, ALT).
    """
    seen: set = set()
    out: list[tuple[str, str, str, str]] = []
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
            key = (chrom, pos, ref, alt)
            if key in seen:
                continue
            if max_af is not None and _row_max_af(row, af_cols) > max_af:
                n_af += 1
                continue
            if not _overlaps_bed(candidate_bed, chrom, pos, ref):
                n_bed += 1
                continue
            seen.add(key)
            out.append(key)
    return out, n_af, n_star, n_bed


# ---- local GeneBe DB (bgzip + tabix) --------------------------------

class GeneBeDBError(RuntimeError):
    pass


def read_db_columns(db: Path) -> dict[str, int]:
    """Read the '#'-prefixed header via `tabix -H` → {col_name: index}."""
    try:
        out = subprocess.run(
            ["tabix", "-H", str(db)],
            check=True, text=True, capture_output=True,
        ).stdout
    except subprocess.CalledProcessError as e:
        raise GeneBeDBError(f"tabix -H failed: {e.stderr.strip()}") from e
    header = ""
    for line in out.splitlines():
        if line.startswith("#"):
            header = line
            break
    if not header:
        raise GeneBeDBError("DB has no '#'-prefixed header line")
    cols = header.lstrip("#").split("\t")
    return {name: i for i, name in enumerate(cols)}


def parse_db_rows(lines, idx: dict[str, int]) -> dict[tuple[str, str, str, str],
                                                      tuple[str, str]]:
    """{(chr,pos,ref,alt): (acmg_score, acmg_criteria)} from tabix output.

    Pure (no subprocess) so it can be unit-tested with synthetic lines.
    Defensive: skips short rows and rows whose pos isn't an integer (the
    DB occasionally carries malformed '.' placeholder rows).
    """
    ci_chr, ci_pos = idx["chr"], idx["pos"]
    ci_ref, ci_alt = idx["ref"], idx["alt"]
    ci_score, ci_crit = idx["acmg_score"], idx["acmg_criteria"]
    need = max(ci_chr, ci_pos, ci_ref, ci_alt, ci_score, ci_crit)
    hits: dict[tuple[str, str, str, str], tuple[str, str]] = {}
    for line in lines:
        if not line or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) <= need:
            continue
        pos = f[ci_pos]
        if not pos.isdigit():
            continue
        hits[(f[ci_chr], pos, f[ci_ref], f[ci_alt])] = (f[ci_score], f[ci_crit])
    return hits


def query_db(db: Path, regions_bed: Path, idx: dict[str, int]) -> dict:
    """Run `tabix -R <regions> <db>` and parse the rows."""
    proc = subprocess.Popen(
        ["tabix", "-R", str(regions_bed), str(db)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    assert proc.stdout is not None
    hits = parse_db_rows(proc.stdout, idx)
    _, err = proc.communicate()
    # tabix exits non-zero only on real index errors; a region whose
    # contig is absent from the index is just a stderr warning. We
    # already probed the DB in preflight, so surface stderr but don't
    # abort on the benign case (we still got rows).
    if proc.returncode not in (0, None) and not hits:
        raise GeneBeDBError(f"tabix -R failed: {(err or '').strip()}")
    return hits


def assemble_gb(
    candidates: list[tuple[str, str, str, str]],
    db_hits: dict[tuple[str, str, str, str], tuple[str, str]],
) -> dict[tuple[str, str, str, str], tuple[str, str, str]]:
    """Map DB hits back onto the TSV's own keys + derive the class.

    Returns {(tsv_chrom,pos,ref,alt): (score, criteria, class)} for the
    candidates found in the DB (with a usable score and/or criteria).
    """
    gb: dict = {}
    for (chrom, pos, ref, alt) in candidates:
        db_key = (_chr_prefixed(chrom), pos, ref, alt)
        hit = db_hits.get(db_key)
        if not hit:
            continue
        score_raw, crit_raw = hit
        score_raw = (score_raw or "").strip()
        crit_raw = (crit_raw or "").strip()
        score_out = "" if score_raw in _MISSING else score_raw
        crit_out = "" if crit_raw in _MISSING else crit_raw
        cls = classify(_to_float(score_raw))
        if score_out or crit_out or cls:
            gb[(chrom, pos, ref, alt)] = (score_out, crit_out, cls)
    return gb


def write_regions_bed(candidates: list[tuple[str, str, str, str]], path: Path) -> int:
    """One sorted, unique 1-bp BED region per candidate POS (covers the
    DB row whose tabix begin=end=pos)."""
    seen: set[tuple[str, int]] = set()
    for chrom, pos, _ref, _alt in candidates:
        if pos.isdigit():
            seen.add((_chr_prefixed(chrom), int(pos)))
    rows = sorted(seen, key=lambda x: (x[0], x[1]))
    with open(path, "w", encoding="utf-8") as fo:
        for chrom, pos in rows:
            fo.write(f"{chrom}\t{pos - 1}\t{pos}\n")
    return len(rows)


def preflight_db(db: Path) -> dict[str, int]:
    """Validate tabix + DB + index, return the column index map."""
    if shutil.which("tabix") is None:
        raise GeneBeDBError("`tabix` not found on PATH (install htslib)")
    if not db.is_file():
        raise GeneBeDBError(f"GeneBe DB not found: {db}")
    tbi = Path(str(db) + ".tbi")
    csi = Path(str(db) + ".csi")
    index = tbi if tbi.is_file() else (csi if csi.is_file() else None)
    if index is None:
        raise GeneBeDBError(
            f"missing tabix index for {db.name} — run: "
            f"tabix -s 1 -b 2 -e 2 -c '#' {db}")
    if index.stat().st_mtime < db.stat().st_mtime:
        raise GeneBeDBError(
            f"index {index.name} is OLDER than the data file — rebuild it: "
            f"tabix -s 1 -b 2 -e 2 -c '#' {db}")
    idx = read_db_columns(db)
    for col in ("chr", "pos", "ref", "alt", "acmg_score", "acmg_criteria"):
        if col not in idx:
            raise GeneBeDBError(
                f"DB header missing required column '{col}'; got {list(idx)}")
    return idx


# ---- merge back into the TSV ----------------------------------------

def merge_into_tsv(in_tsv: Path, out_tsv: Path, gb: dict) -> tuple[int, int]:
    """Write a new TSV with GENEBE_* backfilled. Returns (n_filled, n_total).

    Only blank GENEBE_* cells are populated; pipeline ACMG_* is untouched.
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
                key = (
                    (row.get("CHROM") or "").strip(),
                    (row.get("POS")   or "").strip(),
                    (row.get("REF")   or "").strip(),
                    (row.get("ALT")   or "").strip(),
                )
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
                    help="local GeneBe DB (bgzip + tabix TSV); "
                         "env NGS_UI_GENEBE_DB")
    ap.add_argument("--workdir",
                    help="keep the candidate regions BED here "
                         "(default: temp dir, removed on success)")
    ap.add_argument("--max-af", type=float, default=0.01,
                    help="drop sites whose GNOMAD_G_AF > this "
                         "(default 0.01; use -1 to disable)")
    ap.add_argument("--af-cols", default="GNOMAD_G_AF",
                    help="comma-separated AF columns to check "
                         "(default GNOMAD_G_AF)")
    ap.add_argument("--candidate-bed",
                    help="BED regions eligible for lookup "
                         "(0-based half-open; variants outside are skipped)")
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
        n_intervals = sum(len(intervals) for _, intervals in candidate_bed.values())
        print(f"[genebe] candidate BED enabled: {bed_path} "
              f"({n_intervals} merged intervals)", file=sys.stderr)

    max_af = None if args.max_af < 0 else args.max_af
    af_cols = [c.strip() for c in args.af_cols.split(",") if c.strip()]
    candidates, n_af, n_star, n_bed = collect_candidates(
        in_tsv, max_af=max_af, af_cols=af_cols, candidate_bed=candidate_bed,
    )
    af_note = (f"dropped {n_af} above AF {max_af}"
               if max_af is not None else "AF filter off")
    print(f"[genebe] {len(candidates)} unique candidate sites "
          f"({af_note}; dropped {n_bed} outside candidate BED; "
          f"skipped {n_star} with '*')", file=sys.stderr)
    if not candidates:
        print("[genebe] 0 candidates — nothing to look up", file=sys.stderr)
        # Still ensure the GENEBE_* columns exist for a stable schema.
        merge_into_tsv(in_tsv, out_tsv, {})
        return 0

    if args.workdir:
        wd = Path(args.workdir)
        wd.mkdir(parents=True, exist_ok=True)
        wd_ctx = None
    else:
        wd_ctx = tempfile.TemporaryDirectory(prefix="genebe-")
        wd = Path(wd_ctx.name)
    try:
        regions = wd / "candidates.bed"
        n_regions = write_regions_bed(candidates, regions)
        try:
            db_hits = query_db(db, regions, idx)
        except GeneBeDBError as e:
            print(f"ERROR: GeneBe DB 查詢失敗：{e}", file=sys.stderr)
            return 2
        gb = assemble_gb(candidates, db_hits)
        n_filled, n_total = merge_into_tsv(in_tsv, out_tsv, gb)
        print(f"[genebe] {len(db_hits)} DB rows over {n_regions} positions; "
              f"backfilled ACMG for {n_filled}/{n_total} TSV rows "
              f"({len(gb)} candidates matched)", file=sys.stderr)
        print(f"[genebe] done → {out_tsv}", file=sys.stderr)
    finally:
        if wd_ctx is not None:
            wd_ctx.cleanup()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
