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

The rebuild runs PER CHROMOSOME IN PARALLEL. The prefix sum resets at every
chromosome change, so chromosomes are independent and splitting is exactly
equivalent to one global sort — but a single global sort put everything
through one single-threaded awk fed by one single-threaded gzip, which took
days at 1397 WGS. Two further wins: consecutive gVCF ref blocks are exactly
adjacent with the same ploidy weight, so they are coalesced into one interval
before emitting events (orders of magnitude fewer events), and the prefix sum
merges touching segments of equal AN, so the track comes out canonical and
smaller. `--selftest` proves the AN at every base is unchanged.

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
import tempfile
import time
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


def _pick(*names):
    """First tool on PATH, or None."""
    for n in names:
        if shutil.which(n):
            return n
    return None


def _decomp_cmd():
    """Fastest available .gz reader. pigz decompresses with several threads."""
    t = _pick("pigz", "gzip", "zcat") or "gzip"
    return "zcat" if t == "zcat" else f"{t} -dc"


def _awk_bin():
    """mawk is several times faster than gawk on this kind of line loop."""
    return _pick("mawk", "gawk", "awk") or "awk"


def _run_parallel(cmds, jobs):
    """Run shell commands with at most `jobs` in flight; raise on any failure."""
    running, failed = [], []
    cmds = list(cmds)
    while cmds or running:
        while cmds and len(running) < jobs:
            c = cmds.pop(0)
            running.append((c, subprocess.Popen(["bash", "-c", "set -o pipefail; " + c])))
        time.sleep(0.2)
        for item in running[:]:
            c, p = item
            rc = p.poll()
            if rc is None:
                continue
            running.remove(item)
            if rc != 0:
                failed.append((rc, c))
    if failed:
        rc, c = failed[0]
        raise subprocess.CalledProcessError(rc, c)


def rebuild_an_track_single(old_track, new_beds, out_path, bgzip_bin, sort_tmp):
    """Reference implementation: one global sort. Correct but single-threaded
    in the prefix-sum, so it does not scale past a few hundred samples. Kept
    because --selftest checks the parallel path produces identical output."""
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
    comp = "--compress-program=gzip" if shutil.which("gzip") else ""
    sort_cmd = f"LC_ALL=C sort -T {shlex.quote(sort_tmp)} -S 50% {comp} -k1,1 -k2,2n"
    pipeline = ("( " + " ; ".join(parts) + " ) | " + sort_cmd + " | "
                + PREFIXSUM_AWK + f" | {bgzip_bin} > {shlex.quote(tmp)}")
    subprocess.run(["bash", "-c", "set -o pipefail; " + pipeline], check=True)
    os.replace(tmp, out_path)


def rebuild_an_track(old_track, new_beds, out_path, bgzip_bin, sort_tmp,
                     jobs=None, keep_work=False):
    """Per-chromosome parallel rebuild -> out_path (atomic).

    WHY: the prefix sum resets at every chromosome change
    (`$1!=c{c=$1; prev=-1; s=0}`), so chromosomes are completely independent.
    Splitting the event stream per chromosome is therefore EXACTLY equivalent
    to one global sort, but lets every chromosome sort and scan in parallel.
    That matters: at 1397 WGS the single global sort spent days in one
    single-threaded prefix-sum awk fed by one single-threaded gzip.

    Three phases:
      1. demux  — read every input once, split events into <chrom>.<worker>
                  files holding just "pos<TAB>delta" (chrom is the filename,
                  so the rows are ~45% smaller and the sort key is one
                  numeric field).
      2. reduce — per chromosome: sort -k1,1n | prefix-sum -> <chrom>.bg.
                  Chromosomes run in parallel.
      3. concat — chromosomes in LEXICOGRAPHIC order -> bgzip -> out_path.

    Output chromosome order MUST stay lexicographic: publish_af.py merge-joins
    this track against `SELECT ... ORDER BY chrom,pos` from SQLite, which is
    byte order — the same order a single `sort -k1,1` produced.
    """
    inputs = []
    if old_track and os.path.exists(old_track):
        inputs.append(old_track)
    inputs.extend(new_beds)
    if not inputs:
        return

    jobs = jobs or min(os.cpu_count() or 4, 16)
    os.makedirs(sort_tmp, exist_ok=True)
    work = tempfile.mkdtemp(prefix="an_track.", dir=sort_tmp)
    dec, awk = _decomp_cmd(), _awk_bin()
    try:
        # ---- phase 1: demux events per chromosome -------------------------
        # One awk PER INPUT FILE (not per worker) so the coalescing state below
        # resets at every file boundary; output uses >> because each awk exits
        # after its file. The work dir is fresh, so the first >> starts empty.
        #
        # COALESCING: a gVCF's consecutive ref blocks are exactly adjacent
        # (start[i+1] == end[i]) and carry the same ploidy weight, so a run of
        # callable blocks collapses into ONE interval and emits one +w/-w pair
        # instead of thousands. The AN step function is unchanged (the selftest
        # proves it) but the event volume drops by orders of magnitude — which
        # is what made the old single global sort take days.
        chunks = [inputs[i::jobs] for i in range(jobs)]
        cmds = []
        for w, chunk in enumerate(chunks):
            if not chunk:
                continue
            lst = os.path.join(work, f"in.{w}.txt")
            with open(lst, "w") as f:
                f.write("\n".join(chunk) + "\n")
            coalesce = (
                f'{awk} -v d={shlex.quote(work)} -v id={w} '
                "'{c=$1; s=$2+0; e=$3+0; wt=$4+0;"
                " if(c==pc && s==pe && wt==pw){pe=e; next}"
                " if(pc!=\"\"){f=d\"/\"pc\".\"id;"
                " print ps\"\\t\"pw >> f; print pe\"\\t\"(-pw) >> f}"
                " pc=c; ps=s; pe=e; pw=wt}"
                " END{if(pc!=\"\"){f=d\"/\"pc\".\"id;"
                " print ps\"\\t\"pw >> f; print pe\"\\t\"(-pw) >> f}}'")
            cmds.append(
                f'while IFS= read -r f; do {{ case "$f" in *.gz) {dec} "$f";;'
                f' *) cat "$f";; esac; }} | {coalesce};'
                f' done < {shlex.quote(lst)}')
        print(f"[accumulate]   phase 1/3 demux: {len(inputs)} inputs, {len(cmds)} workers "
              f"({dec}, {awk})", file=sys.stderr)
        _run_parallel(cmds, jobs)

        chroms = sorted({f.rsplit(".", 1)[0] for f in os.listdir(work)
                         if not f.startswith("in.") and not f.endswith(".bg")})
        if not chroms:
            raise RuntimeError("AN-track demux produced no events")

        # ---- phase 2: per-chromosome sort + prefix-sum --------------------
        npar = max(1, min(len(chroms), jobs))
        mem = max(2, 60 // npar)
        cmds = []
        for c in chroms:
            bg = os.path.join(work, c + ".bg")
            # prefix-sum + canonicalise: hold the pending interval and extend it
            # while the AN is unchanged, so touching segments with equal AN
            # (events that cancel across samples) collapse into one row.
            prefix = (f"{awk} -v c={shlex.quote(c)} 'BEGIN{{prev=-1; s=0; hv=-1}}"
                      "{p=$1+0; d=$2+0;"
                      " if(prev>=0 && p>prev && s>0){"
                      "  if(hv==s && he==prev){he=p}"
                      "  else {if(hv>=0) printf \"%s\\t%d\\t%d\\t%d\\n\", c, hs, he, hv;"
                      "        hs=prev; he=p; hv=s}}"
                      " s+=d; prev=p}"
                      " END{if(hv>=0) printf \"%s\\t%d\\t%d\\t%d\\n\", c, hs, he, hv}'")
            cmds.append(
                f"cat {shlex.quote(work)}/{shlex.quote(c)}.[0-9]* | "
                f"LC_ALL=C sort -k1,1n -S {mem}% -T {shlex.quote(sort_tmp)} | "
                f"{prefix} > {shlex.quote(bg)}")
        print(f"[accumulate]   phase 2/3 reduce: {len(chroms)} chromosomes, "
              f"{npar} in parallel (-S {mem}% each)", file=sys.stderr)
        _run_parallel(cmds, npar)

        # ---- phase 3: concat in lexicographic order -----------------------
        tmp = out_path + ".tmp"
        bgs = " ".join(shlex.quote(os.path.join(work, c + ".bg")) for c in chroms)
        print("[accumulate]   phase 3/3 concat + bgzip", file=sys.stderr)
        subprocess.run(["bash", "-c", "set -o pipefail; "
                        f"cat {bgs} | {bgzip_bin} > {shlex.quote(tmp)}"], check=True)
        os.replace(tmp, out_path)
    finally:
        if not keep_work:
            shutil.rmtree(work, ignore_errors=True)
    # Only index when the output is really BGZF. bgzip_bin falls back to plain
    # gzip when bgzip is absent (and --selftest passes "gzip"), and tabix would
    # then just warn "not BGZF" — noise that looks like a failure.
    if (out_path.endswith(".gz") and os.path.basename(bgzip_bin) == "bgzip"
            and shutil.which("tabix")):
        subprocess.run(["tabix", "-f", "-p", "bed", out_path], check=False)


def _write_bed(path, rows):
    import gzip as _gz
    with _gz.open(path, "wt") as f:
        for c, s, e, w in rows:
            f.write(f"{c}\t{s}\t{e}\t{w}\n")


def _read_gz(path):
    import gzip as _gz
    with _gz.open(path, "rt") as f:
        return f.read()


def _per_base(text):
    """Expand a small bedGraph to {(chrom,pos): AN} — the ground truth the AN
    track actually encodes, independent of how it is chopped into rows."""
    m = {}
    for line in text.strip().split("\n"):
        if not line:
            continue
        c, s, e, v = line.split("\t")
        for p in range(int(s), int(e)):
            m[(c, p)] = int(v)
    return m


def _has_noop_boundary(text):
    """True if any two adjacent rows touch and carry the same AN (a boundary
    that encodes nothing). The coalescing demux should leave none."""
    prev = None
    for line in text.strip().split("\n"):
        if not line:
            continue
        c, s, e, v = line.split("\t")
        cur = (c, int(s), int(e), int(v))
        if prev and prev[0] == cur[0] and prev[2] == cur[1] and prev[3] == cur[3]:
            return True
        prev = cur
    return False


def selftest():
    """The parallel per-chromosome rebuild must encode the SAME AN at every
    base as one global sort — that equivalence is the whole justification for
    splitting. It is not byte-identical: coalescing also drops no-op row
    boundaries, so the track is canonical and smaller (same AN everywhere)."""
    work = tempfile.mkdtemp(prefix="an_selftest.")
    try:
        # overlapping intervals, mixed weights, and chr1/chr2/chr10 so that
        # lexicographic vs natural chromosome order actually differ.
        samples = [
            [("chr1", 10, 20, 2), ("chr1", 30, 40, 2), ("chr2", 5, 15, 1), ("chr10", 1, 5, 2)],
            [("chr1", 15, 35, 2), ("chr2", 10, 20, 2), ("chr10", 3, 9, 1), ("chrX", 7, 11, 1)],
            [("chr1", 10, 20, 1), ("chrM", 1, 100, 1), ("chr2", 5, 15, 2)],
        ]
        beds = []
        for i, rows in enumerate(samples):
            p = os.path.join(work, f"s{i}.bed.gz")
            _write_bed(p, rows)
            beds.append(p)
        st = os.path.join(work, "sorttmp")

        # --- full rebuild: parallel == single -----------------------------
        a, b = os.path.join(work, "a.gz"), os.path.join(work, "b.gz")
        rebuild_an_track_single(None, beds, a, "gzip", st)
        rebuild_an_track(None, beds, b, "gzip", st, jobs=4)
        ta, tb = _read_gz(a), _read_gz(b)
        assert ta.strip() and tb.strip(), "empty output"
        pa, pb = _per_base(ta), _per_base(tb)
        assert pa == pb, ("FULL MISMATCH at "
                          f"{sorted(k for k in set(pa) | set(pb) if pa.get(k) != pb.get(k))[:10]}")
        # the parallel track must additionally be canonical (no no-op boundaries)
        assert not _has_noop_boundary(tb), "parallel track still has no-op boundaries"

        rows = [l.split("\t") for l in tb.strip().split("\n")]
        got = {(r[0], int(r[1]), int(r[2])): int(r[3]) for r in rows}
        # chr1 15-20 is covered by 2 (s0) + 2 (s1) + 1 (s2) = 5
        assert got.get(("chr1", 15, 20)) == 5, got
        # chr1 20-30 only s1 remains = 2
        assert got.get(("chr1", 20, 30)) == 2, got
        # chromosome order must be lexicographic (chr1, chr10, chr2, chrM, chrX)
        order = []
        for r in rows:
            if not order or order[-1] != r[0]:
                order.append(r[0])
        assert order == sorted(set(order)), f"chrom order not lexicographic: {order}"

        # --- incremental: old track + new BED, parallel == single ---------
        base_a, base_b = os.path.join(work, "ba.gz"), os.path.join(work, "bb.gz")
        rebuild_an_track_single(None, beds[:2], base_a, "gzip", st)
        shutil.copyfile(base_a, base_b)
        inc_a, inc_b = os.path.join(work, "ia.gz"), os.path.join(work, "ib.gz")
        rebuild_an_track_single(base_a, beds[2:], inc_a, "gzip", st)
        rebuild_an_track(base_b, beds[2:], inc_b, "gzip", st, jobs=4)
        ia, ib = _per_base(_read_gz(inc_a)), _per_base(_read_gz(inc_b))
        assert ia == ib, "INCREMENTAL MISMATCH (parallel vs single)"
        # incremental (2 then +1) must equal the full rebuild of all 3
        assert ia == pa, "INCREMENTAL != FULL"

        print(f"selftest OK — same AN at every base (parallel == single, "
              f"full == incremental); parallel track canonical: "
              f"{len(rows)} rows vs {len(ta.strip().splitlines())} unmerged, chroms={order}")
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db-dir")
    ap.add_argument("--per-sample-dir")
    ap.add_argument("--samples", help="comma-separated sample ids (default: all new)")
    ap.add_argument("--sort-tmp", help="big scratch dir for sort (default <db-dir>/.sorttmp; NOT /tmp)")
    ap.add_argument("--jobs", type=int, default=None,
                    help="parallel workers for the AN-track rebuild (default: min(nproc,16))")
    ap.add_argument("--rebuild-an-track", action="store_true",
                    help="force a full AN-track rebuild from every ingested sample, even when "
                         "no new samples were added (use after an interrupted rebuild left "
                         "counts.sqlite ahead of the track)")
    ap.add_argument("--selftest", action="store_true",
                    help="verify the parallel AN-track rebuild matches the single-sort "
                         "reference implementation, then exit")
    args = ap.parse_args()

    if args.selftest:
        return selftest()
    if not args.db_dir:
        ap.error("--db-dir required (or --selftest)")

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
    have_track = os.path.exists(an_track) and not args.rebuild_an_track
    if have_track:
        bed_ids = new
        old = an_track
    else:
        # full rebuild from every ingested sample: fresh DB, missing track, or
        # --rebuild-an-track after an interrupted rebuild left counts.sqlite
        # ahead of the track (re-running without this would silently keep the
        # stale track and inflate every AF).
        bed_ids = ingested
        old = None
    beds = [b for b in (find_bed(os.path.join(per_sample, s)) for s in bed_ids) if b]

    if beds:
        kind = "old + %d" % len(beds) if have_track else "full, %d BEDs" % len(beds)
        print(f"[accumulate] rebuilding AN track ({kind}); sort tmp={sort_tmp}…", file=sys.stderr)
        t0 = time.time()
        rebuild_an_track(old, beds, an_track, bgzip_bin, sort_tmp, jobs=args.jobs)
        print(f"[accumulate] AN track -> {an_track}  ({time.time()-t0:.0f}s)", file=sys.stderr)
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
