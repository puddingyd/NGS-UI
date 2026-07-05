#!/usr/bin/env python3
"""Select the in-house allele-frequency cohort from a list of DRAGEN gVCFs.

The raw datalake holds the same sample under two locations within one run
(``.../<run>/other/<sample>/<sample>.hard-filtered.gvcf.gz`` and
``.../<run>/vcf.gz/<sample>.hard-filtered.gvcf.gz``), plus the occasional
broken/empty gVCF and standard-reference (control) samples that must NOT
contribute to a population allele frequency. This script turns a flat list
of gVCF paths into a clean, deduplicated cohort:

  * one path per sample id (prefer the ``/vcf.gz/`` copy, larger size wins
    the tiebreak),
  * drop gVCFs smaller than ``--min-size-mb`` (failed/empty runs),
  * drop explicit ids (``--exclude-id-file``) and numeric ranges
    (``--exclude-range PREFIX:START-END``, e.g. the ``VAL-37``..``VAL-54``
    standard references).

Stdlib only — runs on the cluster head node without any container.

INPUT (one of):
  --sizes FILE   `du -h`-style TSV: "<human-size>\\t<path>" per line
                 (this is what `du -h .../*.hard-filtered.gvcf.gz` emits).
  --paths FILE   plain list of gVCF paths, one per line (size checks then
                 rely on `os.stat`, so the paths must be reachable).

OUTPUT:
  --out-manifest FILE  audit TSV: sample_id, run, size, status, path
                       (every input sample, included AND excluded).
  --out-list FILE      included gVCF paths only, one per line — this is the
                       file you feed to build_inhouse_af.sh.

Example (reproduce the 613-sample Phase 0 cohort):
  scripts/inhouse_af/select_cohort.py \\
    --sizes hard_filtered_gvcf_sizes.txt \\
    --exclude-range 'VAL-:37-54' \\
    --out-manifest cohort_manifest.tsv \\
    --out-list cohort_gvcfs.txt
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

DEFAULT_SUFFIX = ".hard-filtered.gvcf.gz"
_RUN_RE = re.compile(r"^\d{8}_[A-Za-z0-9_]+$")  # e.g. 20251118_LH00873_0004_A237NFYLT4
_HSIZE_RE = re.compile(r"^([0-9.]+)\s*([KMGTP]?)i?B?$", re.IGNORECASE)
_UNIT = {"": 1, "K": 1024, "M": 1024 ** 2, "G": 1024 ** 3, "T": 1024 ** 4, "P": 1024 ** 5}


def human_to_bytes(s: str) -> int:
    """'1.2M' / '3.3G' / '512' → bytes. Returns -1 if unparseable."""
    m = _HSIZE_RE.match(s.strip())
    if not m:
        return -1
    return int(float(m.group(1)) * _UNIT[m.group(2).upper()])


def bytes_to_human(n: int) -> str:
    if n < 0:
        return "?"
    v = float(n)
    for u in ("", "K", "M", "G", "T"):
        if v < 1024 or u == "T":
            return (f"{v:.0f}{u}" if u == "" else f"{v:.1f}{u}")
        v /= 1024
    return f"{v:.1f}T"


def sample_id_of(path: str, suffix: str) -> str:
    name = os.path.basename(path)
    return name[: -len(suffix)] if name.endswith(suffix) else name


def run_of(path: str) -> str:
    """Best-effort run/batch label = first path segment that looks like a
    DRAGEN run folder (YYYYMMDD_...); fall back to the grandparent dir."""
    parts = Path(path).parts
    for p in parts:
        if _RUN_RE.match(p):
            return p
    return parts[-3] if len(parts) >= 3 else "unknown"


def parse_range(spec: str) -> tuple[str, int, int]:
    """'VAL-:37-54' → ('VAL-', 37, 54). Matches ids '{prefix}{N}'."""
    try:
        prefix, rng = spec.split(":", 1)
        lo, hi = rng.split("-", 1)
        return prefix, int(lo), int(hi)
    except ValueError:
        raise SystemExit(f"--exclude-range must be PREFIX:START-END, got: {spec!r}")


def in_range(sid: str, ranges: list[tuple[str, int, int]]) -> bool:
    for prefix, lo, hi in ranges:
        if sid.startswith(prefix):
            tail = sid[len(prefix):]
            if tail.isdigit() and lo <= int(tail) <= hi:
                return True
    return False


def read_inputs(args) -> list[tuple[str, int]]:
    """→ list of (path, size_bytes). size_bytes = -1 when unknown."""
    rows: list[tuple[str, int]] = []
    if args.sizes:
        with open(args.sizes, encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                # "<size>\t<path>"  (du -h). Be tolerant of multiple spaces.
                parts = line.split("\t") if "\t" in line else line.split(None, 1)
                if len(parts) != 2:
                    print(f"[warn] skip unparseable line: {line!r}", file=sys.stderr)
                    continue
                size_s, path = parts[0].strip(), parts[1].strip()
                rows.append((path, human_to_bytes(size_s)))
    else:
        with open(args.paths, encoding="utf-8") as f:
            for line in f:
                path = line.strip()
                if not path or path.startswith("#"):
                    continue
                try:
                    size = os.stat(path).st_size
                except OSError:
                    size = -1
                rows.append((path, size))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--sizes", help="du -h-style TSV: <size>\\t<path>")
    src.add_argument("--paths", help="plain list of gVCF paths")
    ap.add_argument("--suffix", default=DEFAULT_SUFFIX,
                    help=f"gVCF filename suffix to strip (default {DEFAULT_SUFFIX})")
    ap.add_argument("--min-size-mb", type=float, default=50.0,
                    help="drop gVCFs smaller than this (failed/empty); default 50")
    ap.add_argument("--exclude-id-file", action="append", default=[],
                    help="file of sample ids to exclude (one per line, # comments ok)")
    ap.add_argument("--exclude-range", action="append", default=[],
                    help="numeric id range PREFIX:START-END (repeatable), "
                         "e.g. 'VAL-:37-54' for standard references")
    ap.add_argument("--out-manifest", help="write audit TSV here")
    ap.add_argument("--out-list", help="write included gVCF paths here")
    args = ap.parse_args()

    ranges = [parse_range(s) for s in args.exclude_range]
    excl_ids: set[str] = set()
    for fn in args.exclude_id_file:
        with open(fn, encoding="utf-8") as f:
            for line in f:
                line = line.split("#", 1)[0].strip()
                if line:
                    excl_ids.add(line)
    min_bytes = int(args.min_size_mb * 1024 * 1024)

    rows = read_inputs(args)

    # ---- dedup by sample id: prefer /vcf.gz/ copy, larger size as tiebreak ----
    best: dict[str, dict] = {}
    max_seen: dict[str, int] = {}
    for path, size in rows:
        sid = sample_id_of(path, args.suffix)
        max_seen[sid] = max(max_seen.get(sid, -1), size)
        pref = 2 if "/vcf.gz/" in path else 1
        cur = best.get(sid)
        if (cur is None or pref > cur["pref"]
                or (pref == cur["pref"] and size > cur["size"])):
            best[sid] = {"path": path, "size": size, "pref": pref, "run": run_of(path)}

    # ---- classify ----
    records = []
    for sid, info in best.items():
        peak = max_seen[sid]  # judge broken-ness on the largest copy seen
        if peak >= 0 and peak < min_bytes:
            status = "EXCLUDE_broken_empty"
        elif sid in excl_ids:
            status = "EXCLUDE_id"
        elif in_range(sid, ranges):
            status = "EXCLUDE_standard_ref"
        else:
            status = "include"
        records.append((sid, info["run"], info["size"], status, info["path"]))

    records.sort(key=lambda r: (r[3] != "include", r[0]))  # include first, then by id

    included = [r for r in records if r[3] == "include"]

    if args.out_manifest:
        with open(args.out_manifest, "w", encoding="utf-8") as f:
            f.write("sample_id\trun\tgvcf_size\tstatus\tgvcf_path\n")
            for sid, run, size, status, path in records:
                f.write(f"{sid}\t{run}\t{bytes_to_human(size)}\t{status}\t{path}\n")
        print(f"[ok] manifest → {args.out_manifest}", file=sys.stderr)

    if args.out_list:
        with open(args.out_list, "w", encoding="utf-8") as f:
            for sid, run, size, status, path in included:
                f.write(path + "\n")
        print(f"[ok] gVCF list → {args.out_list}", file=sys.stderr)

    # ---- summary ----
    from collections import Counter
    by_status = Counter(r[3] for r in records)
    by_run = Counter(r[1] for r in included)
    print("", file=sys.stderr)
    print(f"[summary] input paths     : {len(rows)}", file=sys.stderr)
    print(f"[summary] unique samples  : {len(records)}", file=sys.stderr)
    for st, n in sorted(by_status.items()):
        print(f"[summary]   {st:<22} {n}", file=sys.stderr)
    print(f"[summary] INCLUDED        : {len(included)}", file=sys.stderr)
    print(f"[summary] runs/batches    : {len(by_run)}", file=sys.stderr)
    if not included:
        print("[error] no samples included — check filters", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
