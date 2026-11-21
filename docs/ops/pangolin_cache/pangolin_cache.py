#!/usr/bin/env python3
"""Persistent Pangolin score cache for the tertiary PANGOLIN_SCORE process.

A Pangolin score depends only on the variant (CHROM, POS, REF, ALT), the
reference FASTA, the gene-model DB, the Pangolin build (container) and the
command-line options (--distance, -m). It does not depend on the sample. In
the NCKUH cohort ~98% of every sample's splice candidates were already scored
for an earlier sample, so PANGOLIN_SCORE can look them up and run the GPU
model only for variants it has never seen.

Subcommands (Python 3.6+, standard library only; runs inside the Pangolin
container or on the host):

  split    splice_for_pangolin.vcf -> pangolin_miss.vcf (rows to score) +
           pangolin_hits.tsv (cached values, frozen for this task)
  merge    original rows + cached hits + Pangolin output for the misses ->
           <sample>.pangolin.vcf (same row order as the input), and store the
           newly computed values in the cache
  seed     import existing <sample>.pangolin.vcf.gz outputs; keys whose value
           differs between files are skipped and written to a conflict report
  compare  compare two Pangolin VCFs key by key (validation)
  namespace print the cache namespace for the given settings
  stats    row counts per namespace

The cache is partitioned by a *namespace* built from the container name, the
gene-model DB and reference FASTA (basename + size) and the Pangolin options,
so changing any of them starts an empty partition automatically instead of
reusing scores computed under different conditions.

Every cached row records where it came from (sample, host, time) for audit.
A value of NULL means "Pangolin was run and returned no score" (e.g. the
variant is outside every gene model); it is cached too, so such variants are
not re-run every time.

Downstream (parse_vep_csq.load_pangolin_scores) reads only CHROM/POS/REF/ALT
and the INFO/Pangolin value, which this tool reproduces exactly.
"""
import argparse
import datetime
import gzip
import os
import socket
import sqlite3
import sys

SCHEMA_VERSION = "1"
PANGOLIN_INFO_DEFAULT = (
    '##INFO=<ID=Pangolin,Number=.,Type=String,Description="Pangolin splice scores. '
    'Format: gene|pos:score_change|pos:score_change|warnings,...">'
)
NULL_TOKEN = "\\N"


# ── namespace ─────────────────────────────────────────────────────────────

def _file_tag(path):
    try:
        size = os.path.getsize(path)
    except OSError:
        size = "missing"
    return "{}:{}".format(os.path.basename(path.rstrip("/")), size)


def namespace_from_args(args):
    return "pangolin|sif={}|db={}|ref={}|distance={}|mask={}".format(
        args.sif_name,
        _file_tag(args.pangolin_db),
        _file_tag(args.ref_fasta),
        args.distance,
        args.mask,
    )


def add_namespace_args(parser, need_cache=True):
    if need_cache:
        parser.add_argument("--cache", required=True, help="SQLite cache file")
    parser.add_argument("--sif-name", required=True,
                        help="Pangolin container file name (params.pangolin_sif)")
    parser.add_argument("--pangolin-db", required=True, help="params.pangolin_db")
    parser.add_argument("--ref-fasta", required=True, help="params.ref_fasta")
    parser.add_argument("--distance", required=True, help="pangolin --distance")
    parser.add_argument("--mask", required=True, help="pangolin -m")


# ── database ──────────────────────────────────────────────────────────────

def connect(path):
    directory = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(directory):
        os.makedirs(directory, exist_ok=True)
    conn = sqlite3.connect(path, timeout=600)
    conn.execute("PRAGMA busy_timeout = 600000")
    try:
        conn.execute("PRAGMA journal_mode = WAL")
    except sqlite3.OperationalError:
        pass  # e.g. a filesystem without shared-memory support; rollback journal still works
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS meta (
            key TEXT PRIMARY KEY, value TEXT
        );
        CREATE TABLE IF NOT EXISTS namespaces (
            id INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE, created_at TEXT
        );
        CREATE TABLE IF NOT EXISTS scores (
            ns INTEGER NOT NULL, chrom TEXT NOT NULL, pos INTEGER NOT NULL,
            ref TEXT NOT NULL, alt TEXT NOT NULL,
            pangolin TEXT,
            source TEXT, host TEXT, created_at TEXT,
            PRIMARY KEY (ns, chrom, pos, ref, alt)
        ) WITHOUT ROWID;
        """
    )
    conn.execute("INSERT OR IGNORE INTO meta VALUES ('schema_version', ?)", (SCHEMA_VERSION,))
    conn.commit()
    return conn


def namespace_id(conn, name, create=True):
    row = conn.execute("SELECT id FROM namespaces WHERE name = ?", (name,)).fetchone()
    if row:
        return row[0]
    if not create:
        return None
    conn.execute(
        "INSERT OR IGNORE INTO namespaces (name, created_at) VALUES (?, ?)",
        (name, _now()),
    )
    conn.commit()
    return conn.execute("SELECT id FROM namespaces WHERE name = ?", (name,)).fetchone()[0]


def lookup(conn, ns, keys):
    """Return {key: value-or-None} for the keys present in the cache."""
    found = {}
    if not keys:
        return found
    conn.execute(
        "CREATE TEMP TABLE IF NOT EXISTS wanted "
        "(chrom TEXT, pos INTEGER, ref TEXT, alt TEXT)"
    )
    conn.execute("DELETE FROM temp.wanted")
    conn.executemany("INSERT INTO temp.wanted VALUES (?, ?, ?, ?)", list(keys))
    for chrom, pos, ref, alt, value in conn.execute(
        "SELECT s.chrom, s.pos, s.ref, s.alt, s.pangolin FROM temp.wanted w "
        "JOIN scores s ON s.ns = ? AND s.chrom = w.chrom AND s.pos = w.pos "
        "AND s.ref = w.ref AND s.alt = w.alt",
        (ns,),
    ):
        found[(chrom, int(pos), ref, alt)] = value
    conn.execute("DELETE FROM temp.wanted")
    return found


def store(conn, ns, values, source):
    """Insert computed values; an existing row is never overwritten."""
    now, host = _now(), socket.gethostname()
    before = conn.total_changes
    conn.executemany(
        "INSERT OR IGNORE INTO scores VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [(ns, k[0], k[1], k[2], k[3], v, source, host, now) for k, v in values.items()],
    )
    conn.commit()
    return conn.total_changes - before


def _now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


# ── VCF helpers ───────────────────────────────────────────────────────────

def _open(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def read_vcf(path):
    """Return (header_lines, body_rows); rows keep their exact text (no newline)."""
    header, body = [], []
    with _open(path) as handle:
        for line in handle:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            (header if line.startswith("#") else body).append(line)
    return header, body


def row_key(row):
    fields = row.split("\t", 5)
    if len(fields) < 5:
        raise ValueError("malformed VCF row: {!r}".format(row[:200]))
    return (fields[0], int(fields[1]), fields[3], fields[4])


def pangolin_value(row):
    """INFO/Pangolin of a row, or None when the row has no Pangolin field."""
    fields = row.split("\t")
    if len(fields) < 8:
        return None
    for item in fields[7].split(";"):
        if item.startswith("Pangolin="):
            return item[len("Pangolin="):]
    return None


def with_pangolin(row, value):
    """Row with INFO/Pangolin set to value (None = leave the row unscored)."""
    fields = row.split("\t")
    info = [item for item in fields[7].split(";")
            if item and item != "." and not item.startswith("Pangolin=")]
    if value is not None:
        info.append("Pangolin=" + value)
    fields[7] = ";".join(info) if info else "."
    return "\t".join(fields)


def pangolin_info_line(header):
    for line in header:
        if line.startswith("##INFO=<ID=Pangolin,"):
            return line
    return None


def header_with_pangolin(header, info_line):
    out = [line for line in header if not line.startswith("##INFO=<ID=Pangolin,")]
    index = next((i for i, line in enumerate(out) if line.startswith("#CHROM")), len(out))
    out.insert(index, info_line)
    return out


def _write_lines(path, lines):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        for line in lines:
            handle.write(line + "\n")


# ── subcommands ───────────────────────────────────────────────────────────

def cmd_split(args):
    header, body = read_vcf(args.input)
    keys = []
    for row in body:
        keys.append(row_key(row))
    conn = connect(args.cache)
    ns = namespace_id(conn, namespace_from_args(args))
    cached = lookup(conn, ns, set(keys))
    conn.close()

    miss_rows = [row for row, key in zip(body, keys) if key not in cached]
    _write_lines(args.miss_out, header + miss_rows)
    with open(args.hits_out, "w", encoding="utf-8", newline="") as handle:
        handle.write("#chrom\tpos\tref\talt\tpangolin\n")
        for key in sorted(set(keys) & set(cached), key=lambda k: (k[0], k[1], k[2], k[3])):
            value = cached[key]
            handle.write("{}\t{}\t{}\t{}\t{}\n".format(
                key[0], key[1], key[2], key[3], NULL_TOKEN if value is None else value))
    hits = len(body) - len(miss_rows)
    print("[PANGOLIN-CACHE] {}: candidates={} cached={} to_score={} namespace={}".format(
        args.sample or "-", len(body), hits, len(miss_rows), namespace_from_args(args)),
        file=sys.stderr)
    return 0


def read_hits(path):
    hits = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            chrom, pos, ref, alt, value = line.split("\t", 4)
            hits[(chrom, int(pos), ref, alt)] = None if value == NULL_TOKEN else value
    return hits


def cmd_merge(args):
    header, body = read_vcf(args.input)
    hits = read_hits(args.hits)

    computed = {}
    info_line = None
    miss_keys = {row_key(row) for row in body} - set(hits)
    if miss_keys:
        if not args.pangolin_out or not os.path.isfile(args.pangolin_out):
            print("[PANGOLIN-CACHE] ERROR: {} variants needed scoring but the Pangolin "
                  "output is missing: {}".format(len(miss_keys), args.pangolin_out),
                  file=sys.stderr)
            return 2
        p_header, p_body = read_vcf(args.pangolin_out)
        info_line = pangolin_info_line(p_header)
        for row in p_body:
            computed[row_key(row)] = pangolin_value(row)
        absent = miss_keys - set(computed)
        if absent:
            example = sorted(absent)[0]
            print("[PANGOLIN-CACHE] ERROR: Pangolin output lacks {} of {} scored variants "
                  "(e.g. {}:{} {}>{})".format(len(absent), len(miss_keys), *example),
                  file=sys.stderr)
            return 2

    conn = connect(args.cache)
    ns = namespace_id(conn, namespace_from_args(args))
    if info_line:
        conn.execute("INSERT OR REPLACE INTO meta VALUES ('pangolin_info_line', ?)", (info_line,))
        conn.commit()
    else:
        row = conn.execute(
            "SELECT value FROM meta WHERE key = 'pangolin_info_line'").fetchone()
        info_line = row[0] if row else PANGOLIN_INFO_DEFAULT
    added = store(conn, ns, {k: computed[k] for k in miss_keys}, "sample:" + (args.sample or "-"))
    conn.close()

    out_rows = []
    from_cache = 0
    for row in body:
        key = row_key(row)
        if key in hits:
            value = hits[key]
            from_cache += 1
        else:
            value = computed[key]
        out_rows.append(with_pangolin(row, value))
    _write_lines(args.output, header_with_pangolin(header, info_line) + out_rows)
    print("[PANGOLIN-CACHE] {}: wrote {} rows ({} from cache, {} newly scored, "
          "{} new cache rows)".format(args.sample or "-", len(out_rows), from_cache,
                                      len(out_rows) - from_cache, added),
          file=sys.stderr)
    return 0


def _vcf_values(path):
    _, body = read_vcf(path)
    values = {}
    for row in body:
        values[row_key(row)] = pangolin_value(row)
    return values


def cmd_seed(args):
    min_mtime = None
    if args.since:
        min_mtime = datetime.datetime.strptime(args.since, "%Y-%m-%d").timestamp()
    seen = {}        # key -> value (first seen)
    conflicts = {}   # key -> set of values
    files = 0
    skipped_old = 0
    for path in args.files:
        if min_mtime is not None and os.path.getmtime(path) < min_mtime:
            skipped_old += 1
            continue
        files += 1
        for key, value in _vcf_values(path).items():
            if key in conflicts:
                conflicts[key].add(value)
            elif key in seen:
                if seen[key] != value:
                    conflicts[key] = {seen.pop(key), value}
            else:
                seen[key] = value
    if args.conflicts:
        with open(args.conflicts, "w", encoding="utf-8", newline="") as handle:
            handle.write("#chrom\tpos\tref\talt\tvalues\n")
            for key in sorted(conflicts):
                values = sorted(NULL_TOKEN if v is None else v for v in conflicts[key])
                handle.write("{}\t{}\t{}\t{}\t{}\n".format(*key, " || ".join(values)))
    print("[PANGOLIN-CACHE] seed: files={} skipped_older={} keys={} consistent={} "
          "conflicting={} (conflicts are not cached)".format(
              files, skipped_old, len(seen) + len(conflicts), len(seen), len(conflicts)),
          file=sys.stderr)
    if args.dry_run:
        return 0
    conn = connect(args.cache)
    ns = namespace_id(conn, namespace_from_args(args))
    added = store(conn, ns, seen, "seed")
    conn.close()
    print("[PANGOLIN-CACHE] seed: {} new cache rows".format(added), file=sys.stderr)
    return 0


def cmd_compare(args):
    left, right = _vcf_values(args.left), _vcf_values(args.right)
    only_left = sorted(set(left) - set(right))
    only_right = sorted(set(right) - set(left))
    differ = sorted(k for k in set(left) & set(right) if left[k] != right[k])
    print("compare: left={} right={} shared={} only_left={} only_right={} differ={}".format(
        len(left), len(right), len(set(left) & set(right)),
        len(only_left), len(only_right), len(differ)))
    for key in differ[: args.show]:
        print("DIFF\t{}:{} {}>{}\t{}\t{}".format(*key, left[key], right[key]))
    for key in only_left[: args.show]:
        print("ONLY_LEFT\t{}:{} {}>{}".format(*key))
    for key in only_right[: args.show]:
        print("ONLY_RIGHT\t{}:{} {}>{}".format(*key))
    return 0 if not (only_left or only_right or differ) else 1


def cmd_namespace(args):
    print(namespace_from_args(args))
    return 0


def cmd_stats(args):
    if not os.path.isfile(args.cache):
        print("no cache at {}".format(args.cache))
        return 1
    conn = sqlite3.connect("file:{}?mode=ro".format(os.path.abspath(args.cache)), uri=True)
    for ns_id, name, created in conn.execute(
            "SELECT id, name, created_at FROM namespaces ORDER BY id"):
        total, scored = conn.execute(
            "SELECT COUNT(*), COUNT(pangolin) FROM scores WHERE ns = ?", (ns_id,)).fetchone()
        print("[{}] {}\n    created {}  rows {}  with score {}  no score {}".format(
            ns_id, name, created, total, scored, total - scored))
    conn.close()
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("split", help="split candidates into cached hits and rows to score")
    add_namespace_args(p)
    p.add_argument("--in", dest="input", required=True)
    p.add_argument("--miss-out", required=True)
    p.add_argument("--hits-out", required=True)
    p.add_argument("--sample", default="")
    p.set_defaults(func=cmd_split)

    p = sub.add_parser("merge", help="assemble the sample's Pangolin VCF and store new scores")
    add_namespace_args(p)
    p.add_argument("--in", dest="input", required=True)
    p.add_argument("--hits", required=True)
    p.add_argument("--pangolin-out", default="")
    p.add_argument("--out", dest="output", required=True)
    p.add_argument("--sample", default="")
    p.set_defaults(func=cmd_merge)

    p = sub.add_parser("seed", help="import existing *.pangolin.vcf.gz outputs")
    add_namespace_args(p)
    p.add_argument("--since", default="", help="only files modified on/after YYYY-MM-DD")
    p.add_argument("--conflicts", default="", help="write conflicting keys to this TSV")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("files", nargs="+")
    p.set_defaults(func=cmd_seed)

    p = sub.add_parser("compare", help="compare two Pangolin VCFs key by key")
    p.add_argument("left")
    p.add_argument("right")
    p.add_argument("--show", type=int, default=20)
    p.set_defaults(func=cmd_compare)

    p = sub.add_parser("namespace", help="print the namespace for these settings")
    add_namespace_args(p, need_cache=False)
    p.set_defaults(func=cmd_namespace)

    p = sub.add_parser("stats", help="row counts per namespace")
    p.add_argument("--cache", required=True)
    p.set_defaults(func=cmd_stats)

    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 1
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
