#!/usr/bin/env python3
"""Phase B (1/2) — accumulate per-sample ingest output into the in-house AF DB.

Folds new samples' per_sample/{id}/ contributions (from ingest_sample.py) into:

  counts.sqlite        per normalized site: n_hom/n_het/n_hemi + n_mt_hom/n_mt_het
                       (+ a `samples` table = manifest / dedup)
  an_track.bg.gz       cumulative genome-wide AN = Σ ploidy-weight of callable
                       samples, as a bedGraph (chrom start end AN)

AN-track update uses the event/delta method (approved design): decode the old
track + the new samples' weighted BEDs into (+w at start, -w at end) events,
sort, prefix-sum → the new step function. Only `sort`/`awk`/`bgzip` — no bedtools,
no genotype re-processing (incremental in the meaningful sense).

Idempotent: a sample already in the `samples` table is skipped (counts and the
AN track both gated by it), so re-running a batch never double-counts.

Usage:
  scripts/inhouse_af/accumulate.py --db-dir $NGS_UI_HOME/biotools/inhouse_af
    [--per-sample-dir <dir>]   # default <db-dir>/per_sample
    [--samples id1,id2,...]    # default: every per_sample/ dir not yet ingested
"""
from __future__ import annotations

import argparse
import glob
import os
import shlex
import shutil
import sqlite3
import subprocess
import sys
from datetime import datetime, timezone, timedelta

TAIPEI = timezone(timedelta(hours=8))

KLASS_COL = {
    "hom": "n_hom", "het": "n_het", "hemi": "n_hemi",
    "mt_hom": "n_mt_hom", "mt_het": "n_mt_het",
}
COLS = ["n_hom", "n_het", "n_hemi", "n_mt_hom", "n_mt_het"]

DECODE_AWK = r"""awk -v OFS='\t' '{print $1,$2,$4; print $1,$3,(-$4)}'"""
# prefix-sum sorted (chrom,pos,delta) -> bedGraph (chrom start end AN>0)
PREFIXSUM_AWK = r"""awk -v OFS='\t' '
$1!=c{c=$1; prev=-1; s=0}
{p=$2+0; d=$3+0;
 if(prev>=0 && p>prev && s>0) print c,prev,p,s;
 s+=d; prev=p}'"""


def init_db(conn):
    conn.executescript("""
    PRAGMA journal_mode=WAL;
    PRAGMA synchronous=OFF;
    CREATE TABLE IF NOT EXISTS variant_counts(
      chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
      n_hom INTEGER DEFAULT 0, n_het INTEGER DEFAULT 0, n_hemi INTEGER DEFAULT 0,
      n_mt_hom INTEGER DEFAULT 0, n_mt_het INTEGER DEFAULT 0,
      PRIMARY KEY(chrom,pos,ref,alt)
    );
    CREATE TABLE IF NOT EXISTS samples(sample_id TEXT PRIMARY KEY, added_at TEXT);
    """)


def already_ingested(conn) -> set:
    return {r[0] for r in conn.execute("SELECT sample_id FROM samples")}


def upsert_counts(conn, counts_tsv: str):
    """Add one sample's counts.tsv into variant_counts (delta +1 per row)."""
    rows = []
    with open(counts_tsv, encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            chrom, pos, ref, alt, klass = p[0], int(p[1]), p[2], p[3], p[4]
            col = KLASS_COL.get(klass)
            if not col:
                continue
            d = [0, 0, 0, 0, 0]
            d[COLS.index(col)] = 1
            rows.append((chrom, pos, ref, alt, *d))
    sql = f"""INSERT INTO variant_counts(chrom,pos,ref,alt,{','.join(COLS)})
              VALUES(?,?,?,?,?,?,?,?,?)
              ON CONFLICT(chrom,pos,ref,alt) DO UPDATE SET
              {', '.join(f'{c}={c}+excluded.{c}' for c in COLS)}"""
    conn.executemany(sql, rows)
    return len(rows)


def find_bed(sample_dir: str):
    for name in ("callable.weighted.bed.gz", "callable.weighted.bed"):
        p = os.path.join(sample_dir, name)
        if os.path.exists(p):
            return p
    return None


def rebuild_an_track(old_track, new_beds, out_path, bgzip_bin, sort_tmp):
    """events(old + new) | sort | prefix-sum | bgzip -> out_path (atomic).

    sort spills to sort_tmp (a big disk, NOT /tmp) and compresses temp runs —
    the event stream for a whole cohort is tens of GB."""
    parts = []
    if old_track and os.path.exists(old_track):
        rdr = "zcat" if old_track.endswith(".gz") else "cat"
        parts.append(f"{rdr} {shlex.quote(old_track)} | {DECODE_AWK}")
    for b in new_beds:
        rdr = "zcat" if b.endswith(".gz") else "cat"
        parts.append(f"{rdr} {shlex.quote(b)} | {DECODE_AWK}")
    if not parts:
        return
    os.makedirs(sort_tmp, exist_ok=True)
    tmp = out_path + ".tmp"
    comp = ""
    if shutil.which("gzip"):
        comp = "--compress-program=gzip"
    sort_cmd = f"LC_ALL=C sort -T {shlex.quote(sort_tmp)} -S 50% {comp} -k1,1 -k2,2n"
    pipeline = (
        "( " + " ; ".join(parts) + " ) | " + sort_cmd + " | "
        + PREFIXSUM_AWK + f" | {bgzip_bin} > {shlex.quote(tmp)}"
    )
    subprocess.run(["bash", "-c", "set -o pipefail; " + pipeline], check=True)
    os.replace(tmp, out_path)
    if shutil.which("tabix") and out_path.endswith(".gz"):
        subprocess.run(["tabix", "-f", "-p", "bed", out_path], check=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db-dir", required=True)
    ap.add_argument("--per-sample-dir")
    ap.add_argument("--samples", help="comma-separated sample ids (default: all new)")
    ap.add_argument("--sort-tmp", help="big scratch dir for sort (default <db-dir>/.sorttmp; NOT /tmp)")
    args = ap.parse_args()

    db_dir = args.db_dir
    per_sample = args.per_sample_dir or os.path.join(db_dir, "per_sample")
    os.makedirs(db_dir, exist_ok=True)
    db_path = os.path.join(db_dir, "counts.sqlite")
    an_track = os.path.join(db_dir, "an_track.bg.gz")
    bgzip_bin = shutil.which("bgzip") or "gzip"
    sort_tmp = args.sort_tmp or os.path.join(db_dir, ".sorttmp")

    conn = sqlite3.connect(db_path)
    init_db(conn)
    done = already_ingested(conn)

    if args.samples:
        want = [s.strip() for s in args.samples.split(",") if s.strip()]
    else:
        want = sorted(os.path.basename(d) for d in glob.glob(os.path.join(per_sample, "*"))
                      if os.path.isdir(d))
    new = [s for s in want if s not in done]

    if new:
        print(f"[accumulate] adding {len(new)} sample(s) (have {len(done)})", file=sys.stderr)
        now = datetime.now(TAIPEI).strftime("%Y-%m-%d %H:%M:%S")
        for sid in new:
            d = os.path.join(per_sample, sid)
            counts_tsv = os.path.join(d, "counts.tsv")
            bed = find_bed(d)
            if not (os.path.exists(counts_tsv) and bed):
                print(f"[accumulate]   skip {sid}: missing counts.tsv or callable BED", file=sys.stderr)
                continue
            conn.execute("BEGIN")
            n = upsert_counts(conn, counts_tsv)
            conn.execute("INSERT INTO samples(sample_id,added_at) VALUES(?,?)", (sid, now))
            conn.commit()
            print(f"[accumulate]   + {sid}  ({n} variant rows)", file=sys.stderr)

    # Decide what feeds the AN-track rebuild:
    #  - track present  -> incremental: old track + the BEDs of the new samples
    #  - track missing  -> (re)build from ALL ingested samples (self-healing after
    #    a crash, or a fresh DB)
    ingested = sorted(already_ingested(conn))
    conn.close()
    have_track = os.path.exists(an_track)
    if have_track:
        bed_ids = new
        old = an_track
    else:
        bed_ids = ingested
        old = None
    beds = [b for b in (find_bed(os.path.join(per_sample, s)) for s in bed_ids) if b]

    if beds:
        kind = "old + %d" % len(beds) if have_track else "fresh, %d BEDs" % len(beds)
        print(f"[accumulate] rebuilding AN track ({kind}); sort tmp={sort_tmp}…", file=sys.stderr)
        rebuild_an_track(old, beds, an_track, bgzip_bin, sort_tmp)
        print(f"[accumulate] AN track -> {an_track}", file=sys.stderr)
    else:
        print("[accumulate] AN track unchanged.", file=sys.stderr)

    # report cohort size
    conn = sqlite3.connect(db_path)
    total = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    nvar = conn.execute("SELECT COUNT(*) FROM variant_counts").fetchone()[0]
    conn.close()
    print(f"[accumulate] cohort={total} samples, {nvar} distinct variants", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
