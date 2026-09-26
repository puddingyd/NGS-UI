#!/usr/bin/env python3
"""Annotate snv_indel.annotated.tsv with in-house allele frequency (per-allele).

Adds three columns from the in-house AF sites VCF (built by
scripts/inhouse_af/publish_af.py):

    INHOUSE_AC   in-house alt allele count      (Number=A: per-ALT, comma-sep)
    INHOUSE_AN   in-house total called alleles  (Number=A)
    INHOUSE_AF   INHOUSE_AC / INHOUSE_AN        (Number=A)

MATCHING. The DB is `bcftools norm -m-` split + left-aligned, so a TSV row that
is multiallelic (`ALT="C,CA"`) or whose indel is written differently won't
exact-match. Each ALT is therefore split out and reduced to its minimal
representation (trim common suffix then prefix) before lookup. This is
reference-independent on purpose: normalising with `bcftools norm --check-ref x`
against the *local* FASTA drops every variant whose REF disagrees with it, which
on a machine whose hg38 differs from the DB's build reference lost ~16% of
matches for no speed gain. The minimal-representation join matches ~99.96%.
It misses only indels the DB left-shifted into a repeat (rare).

SPEED (a cohort backfill runs this once per sample, so it matters):
  * the TSV is streamed twice (collect keys, then rewrite) with
    `csv.reader`/`csv.writer` instead of being parsed into a list of dicts, so a
    6M-row WGS TSV never sits in memory. csv rather than raw line splitting
    because fields may be quoted and may contain embedded newlines.
  * the sites VCF is read through `bgzip -dc`/`pigz`/`zcat` when available;
    Python's gzip module is several times slower over ~57M lines.
  A SQLite index of the sites VCF was tried and REJECTED: a python-sqlite3
  indexed lookup costs ~7.5 us against ~0.45 us for a line parse, so replacing
  a 57M-line scan with 6M lookups measured SLOWER, not faster.

Per-allele values are comma-joined in ALT order (`.` for a non-matching allele).
The SNV adapter shows the first ALT's value on the card (consistent with the
existing first-ALT VAF/AD behaviour); the full per-allele data stays in the TSV.

Fill-or-augment, idempotent, atomic replace. No-op (exit 0) when the DB is
missing. Usage:

    scripts/annotate_inhouse_af.py --tsv <snv_indel.annotated.tsv> \\
        [--db <inhouse_af.hg38.vcf.gz>]
    scripts/annotate_inhouse_af.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from array import array
from pathlib import Path

DEFAULT_DB = os.environ.get(
    "NGS_UI_INHOUSE_AF_DB",
    str(Path.home() / "NGS_UI" / "biotools" / "inhouse_af" / "inhouse_af.hg38.vcf.gz"),
)
COL_AC, COL_AN, COL_AF = "INHOUSE_AC", "INHOUSE_AN", "INHOUSE_AF"
_SYMBOLIC = {"*", ".", "", "<NON_REF>", "<*>"}
_MAX_ALT = 255          # alleles per row we can pack into the key index
_CSV_MAX = 4 * 1024 * 1024


def norm_chrom(chrom: str) -> str:
    c = (chrom or "").strip()
    if not c:
        return ""
    return c if c.lower().startswith("chr") else "chr" + c


def info_get(info: str, key: str):
    kv = key + "="
    for field in info.split(";"):
        if field.startswith(kv):
            return field[len(kv):]
    return None


def info3(info: str):
    """(AC, AN, AF) from an INFO string in one pass; '.' when absent."""
    ac = an = af = "."
    for fld in info.split(";"):
        if fld.startswith("INHOUSE_AC="):
            ac = fld[11:]
        elif fld.startswith("INHOUSE_AN="):
            an = fld[11:]
        elif fld.startswith("INHOUSE_AF="):
            af = fld[11:]
    return ac, an, af


def minimal_repr(pos: int, ref: str, alt: str):
    """Parsimonious (minimal) representation: trim common suffix then prefix.

    Reference-independent, so it does not left-shift into a repeat; that is the
    known (rare) miss. Everything else canonicalises to the DB's form."""
    ref = (ref or "").upper()
    alt = (alt or "").upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def alt_alleles(alt_field: str):
    """Split a possibly-multiallelic ALT into (index, allele), skipping symbolic."""
    out = []
    for j, a in enumerate((alt_field or "").split(",")):
        a = a.strip()
        if a in _SYMBOLIC:
            continue
        out.append((j, a))
    return out


# --------------------------------------------------------------------------
# pass 1: TSV -> key index ; lookup ; pass 2: rewrite
# --------------------------------------------------------------------------

def scan_keys(tsv: Path, ic: int, ip: int, ir: int, ia: int):
    """Stream the TSV once. Returns (index, n_alts_per_row).

    index maps a normalised key to packed (row<<8 | allele) ints — a bare int
    for the common single-owner case, a list only on collision."""
    index: dict[tuple, object] = {}
    nalts = array("H")
    with open(tsv, "r", encoding="utf-8", newline="") as f:
        rdr = csv.reader(f, delimiter="\t")
        next(rdr, None)
        for i, row in enumerate(rdr):
            if len(row) <= ia:
                nalts.append(0)
                continue
            toks = (row[ia] or "").split(",")
            nalts.append(min(len(toks), 65535))
            chrom = norm_chrom(row[ic] if ic < len(row) else "")
            pos = (row[ip] or "").strip() if ip < len(row) else ""
            ref_a = (row[ir] or "").strip() if ir < len(row) else ""
            if not (chrom and pos.isdigit() and ref_a):
                continue
            base = i << 8
            for j, a in enumerate(toks):
                if j > _MAX_ALT:
                    break
                a = a.strip()
                if a in _SYMBOLIC:
                    continue
                k = minimal_repr(int(pos), ref_a, a)
                k = (chrom, k[0], k[1], k[2])
                packed = base | j
                v = index.get(k)
                if v is None:
                    index[k] = packed
                elif isinstance(v, int):
                    index[k] = [v, packed]
                else:
                    v.append(packed)
    return index, nalts


def _spread(index, k, val, hits):
    v = index[k]
    if isinstance(v, int):
        hits[v] = val
    else:
        for p in v:
            hits[p] = val


def open_db_text(db: str):
    """Open the sites VCF as text, preferring a C decompressor.

    Returns (stream, finish). `finish()` closes and reaps a child process when
    one was used. Python's gzip module is several times slower over ~57M lines,
    and this is re-read once per sample during a backfill."""
    if db.endswith(".gz"):
        for tool, args in (("bgzip", ["-dc"]), ("pigz", ["-dc"]), ("zcat", [])):
            if shutil.which(tool):
                p = subprocess.Popen([tool, *args, db], stdout=subprocess.PIPE,
                                     stderr=subprocess.DEVNULL, text=True,
                                     bufsize=1 << 20)

                def finish(p=p):
                    try:
                        p.stdout.close()
                    finally:
                        p.wait()
                return p.stdout, finish
        fh = gzip.open(db, "rt", encoding="utf-8")
        return fh, fh.close
    fh = open(db, "r", encoding="utf-8")
    return fh, fh.close


def join_scan(index, db: str):
    """One streaming pass over the sites VCF, keeping only keys we need."""
    hits = {}
    fh, finish = open_db_text(db)
    try:
        for line in fh:
            if not line or line[0] == "#":
                continue
            F = line.rstrip("\n").split("\t", 8)
            if len(F) < 8:
                continue
            try:
                k = (F[0], int(F[1]), F[3].upper(), F[4].upper())
            except ValueError:
                continue
            if k in index:
                _spread(index, k, info3(F[7]), hits)
    finally:
        finish()
    return hits


def write_out(tsv: Path, header, nalts, hits, tmp_path: str) -> int:
    """Stream the TSV again, filling the three columns. Returns rows with a hit."""
    out_header = list(header)
    for c in (COL_AC, COL_AN, COL_AF):
        if c not in out_header:
            out_header.append(c)
    i_ac, i_an, i_af = (out_header.index(c) for c in (COL_AC, COL_AN, COL_AF))
    width = len(out_header)
    n_hit = 0
    with open(tsv, "r", encoding="utf-8", newline="") as fi, \
            open(tmp_path, "w", encoding="utf-8", newline="") as fo:
        rdr = csv.reader(fi, delimiter="\t")
        wtr = csv.writer(fo, delimiter="\t", lineterminator="\n")
        next(rdr, None)
        wtr.writerow(out_header)
        for i, row in enumerate(rdr):
            if len(row) < width:
                row.extend([""] * (width - len(row)))
            n = nalts[i] if i < len(nalts) else 0
            base = i << 8
            ac_p, an_p, af_p, any_hit = [], [], [], False
            for j in range(min(n, _MAX_ALT + 1)):
                v = hits.get(base | j)
                if v:
                    any_hit = True
                    ac_p.append(v[0]); an_p.append(v[1]); af_p.append(v[2])
                else:
                    ac_p.append("."); an_p.append("."); af_p.append(".")
            if any_hit:
                row[i_ac] = ",".join(ac_p)
                row[i_an] = ",".join(an_p)
                row[i_af] = ",".join(af_p)
                n_hit += 1
            else:
                row[i_ac] = row[i_an] = row[i_af] = ""
            wtr.writerow(row)
    return n_hit


def annotate(tsv: Path, db: str) -> int:
    if not os.path.exists(db):
        print(f"[inhouse-af] DB not found: {db} — skipping (no-op)", file=sys.stderr)
        return 0

    with open(tsv, "r", encoding="utf-8", newline="") as f:
        header = next(csv.reader(f, delimiter="\t"), None)
    if not header:
        print(f"[inhouse-af] empty TSV: {tsv}", file=sys.stderr)
        return 0
    cols = {name: i for i, name in enumerate(header)}
    missing = [c for c in ("CHROM", "POS", "REF", "ALT") if c not in cols]
    if missing:
        print(f"[inhouse-af] {tsv} has no {','.join(missing)} column — skipping",
              file=sys.stderr)
        return 0

    index, nalts = scan_keys(tsv, cols["CHROM"], cols["POS"], cols["REF"], cols["ALT"])

    hits = join_scan(index, db)

    fd, tmp_name = tempfile.mkstemp(dir=str(tsv.parent), prefix=tsv.name + ".", suffix=".tmp")
    os.close(fd)
    try:
        n_hit = write_out(tsv, header, nalts, hits, tmp_name)
        os.replace(tmp_name, tsv)
    except BaseException:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        raise

    print(f"[inhouse-af] {len(nalts)} variants, {n_hit} matched in-house AF DB")
    return 0


# --------------------------------------------------------------------------

def _mk_db(path, rows):
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c, p, r, a, ac, an, af in rows:
            f.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\t.\t"
                    f"INHOUSE_AC={ac};INHOUSE_AN={an};INHOUSE_AF={af}\n")


def _read_tsv(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def selftest() -> int:
    assert norm_chrom("1") == "chr1"
    assert minimal_repr(45330228, "CAA", "C") == (45330228, "CAA", "C")
    assert minimal_repr(45330228, "CAA", "CA") == (45330228, "CA", "C")
    assert minimal_repr(100, "AT", "AG") == (101, "T", "G")
    assert alt_alleles("C,CA") == [(0, "C"), (1, "CA")]
    assert alt_alleles("A,*,G") == [(0, "A"), (2, "G")]
    assert info_get("INHOUSE_AC=5;INHOUSE_AF=0.05", "INHOUSE_AC") == "5"
    assert info3("INHOUSE_AC=5;INHOUSE_AN=10;INHOUSE_AF=0.5") == ("5", "10", "0.5")
    assert info3("X=1") == (".", ".", ".")

    d = tempfile.mkdtemp()
    db = os.path.join(d, "inhouse.vcf.gz")
    _mk_db(db, [
        ("chr1", 45330228, "CAA", "C", 1235, 2794, "0.442019"),
        ("chr1", 45330228, "CA", "C", 1136, 2794, "0.406586"),
        ("chr7", 500, "AT", "A", 10, 1000, "0.01"),
    ])
    # a quoted field with an embedded comma+quote, like the MANE_ALL column
    body = ('CHROM\tPOS\tREF\tALT\tGENE\tNOTE\n'
            'chr1\t45330228\tCAA\tC,CA\tMUTYH\t"[{""tx"": ""NM_1"", ""x"": 2}]"\n'
            'chr1\t200\tA\tG\tTP53\tplain\n'
            'chr7\t500\tAT\tA\tBRCA1\t"a,b"\n'
            'chr2\t9\tA\t*\tX\t.\n'
            'chr7\t500\tAT\t*,A\tSPAN\t.\n')

    def fresh():
        p = Path(d) / "snv.tsv"
        p.write_text(body, encoding="utf-8")
        return p

    def check(rows):
        assert rows[0][COL_AF] == "0.442019,0.406586", rows[0]
        assert rows[0][COL_AC] == "1235,1136", rows[0]
        assert rows[0]["NOTE"] == '[{"tx": "NM_1", "x": 2}]', rows[0]   # quoting survived
        assert rows[1][COL_AF] == "", rows[1]                           # miss -> blank
        assert rows[2][COL_AF] == "0.01", rows[2]
        assert rows[2]["NOTE"] == "a,b", rows[2]
        assert rows[3][COL_AF] == "", rows[3]                           # symbolic ALT only
        # Number=A: one value PER ALT TOKEN, so a spanning-deletion '*' holds
        # its slot with '.' and the real allele stays in position 2. Emitting
        # only the non-symbolic values would shift 0.01 onto '*'.
        assert rows[4][COL_AF] == ".,0.01", rows[4]
        assert rows[4][COL_AC] == ".,10", rows[4]

    t = fresh()
    annotate(t, db)
    rows = _read_tsv(t)
    check(rows)

    # the C-decompressor path and the gzip-module path must agree exactly
    real_which = shutil.which
    try:
        shutil.which = lambda _n: None      # force gzip.open
        t2 = Path(d) / "snv2.tsv"
        t2.write_text(body, encoding="utf-8")
        annotate(t2, db)
        assert _read_tsv(t2) == rows, "gzip-module path disagrees with bgzip path"
    finally:
        shutil.which = real_which

    # picks up a refreshed DB (no cached state to go stale)
    _mk_db(db, [("chr1", 45330228, "CAA", "C", 1, 2, "0.5")])
    t = fresh()
    annotate(t, db)
    assert _read_tsv(t)[0][COL_AF] == "0.5,.", _read_tsv(t)[0]

    # idempotent: re-running overwrites in place, no duplicate columns
    before = _read_tsv(t)
    annotate(t, db)
    after = _read_tsv(t)
    assert before == after, "re-run not idempotent"
    with open(t, encoding="utf-8") as f:
        hdr = next(csv.reader(f, delimiter="\t"))
    assert hdr.count(COL_AF) == 1, hdr

    # missing DB is a no-op
    assert annotate(t, os.path.join(d, "nope.vcf.gz")) == 0
    print("selftest OK — per-allele join, quoting preserved, idempotent, "
          "decompressor paths agree")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", type=Path, help="snv_indel.annotated.tsv to annotate in place")
    ap.add_argument("--db", default=DEFAULT_DB, help=f"in-house AF sites VCF (default {DEFAULT_DB})")
    ap.add_argument("--no-cache", action="store_true",
                    help="skip the SQLite lookup index and stream the DB instead "
                         "(slower per sample; use when the DB dir is read-only or tight on space)")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    csv.field_size_limit(_CSV_MAX)
    if args.selftest:
        return selftest()
    if not args.tsv:
        ap.error("--tsv required (or --selftest)")
    if not args.tsv.is_file():
        raise SystemExit(f"--tsv not found: {args.tsv}")
    return annotate(args.tsv, args.db)


if __name__ == "__main__":
    raise SystemExit(main())
