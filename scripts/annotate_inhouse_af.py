#!/usr/bin/env python3
"""Annotate snv_indel.annotated.tsv with in-house allele frequency (per-allele).

Adds three columns from the in-house AF sites VCF (built by
scripts/inhouse_af/publish_af.py):

    INHOUSE_AC   in-house alt allele count      (Number=A: per-ALT, comma-sep)
    INHOUSE_AN   in-house total called alleles  (Number=A)
    INHOUSE_AF   INHOUSE_AC / INHOUSE_AN        (Number=A)

The DB is `bcftools norm -m-` split + left-aligned, so a TSV row that is
multiallelic (`ALT="C,CA"`) or whose indel is represented differently won't
exact-match. To match the DB representation we normalize the TSV variants the
SAME way before joining:

  * PRIMARY (`bcftools`): write the TSV variants to a mini VCF (ID=row_allele),
    `bcftools norm -m- -f ref` (split multiallelics + left-align), sort, then
    `bcftools annotate -a <DB>` — a fast sorted merge-join. Values map back
    per (row, allele). One-shot per sample; C-speed.
  * FALLBACK (pure Python, no bcftools/ref): split ALT on comma + trim to
    minimal representation, then a single streaming pass over the DB. Handles
    multiallelics and simple indels; misses only indels the DB left-shifted
    into a repeat (rare).

Per-allele values are comma-joined in ALT order (`.` for a non-matching
allele). The SNV adapter shows the first ALT's value on the card (consistent
with the existing first-ALT VAF/AD behavior); the full per-allele data stays
in the TSV.

Fill-or-augment, idempotent, atomic replace. No-op (exit 0) when the DB is
missing. Usage:

    scripts/annotate_inhouse_af.py --tsv <snv_indel.annotated.tsv> \\
        [--db <inhouse_af.hg38.vcf.gz>] [--ref <hg38.fasta>]
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
from pathlib import Path

DEFAULT_DB = os.environ.get(
    "NGS_UI_INHOUSE_AF_DB",
    str(Path.home() / "NGS_UI" / "biotools" / "inhouse_af" / "inhouse_af.hg38.vcf.gz"),
)
# Candidate reference FASTAs for the bcftools path (first existing one with a
# .fai wins). Override with --ref or NGS_UI_INHOUSE_AF_REF.
REF_CANDIDATES = [
    os.environ.get("NGS_UI_INHOUSE_AF_REF", ""),
    os.path.join(os.environ.get("NGS_UI_IGV_REF_DIR", "/home/pipeline/reference/hg38"),
                 "Homo_sapiens_assembly38.fasta"),
    "/home/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta",
    "/home/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta",
]
COL_AC, COL_AN, COL_AF = "INHOUSE_AC", "INHOUSE_AN", "INHOUSE_AF"
_SYMBOLIC = {"*", ".", "", "<NON_REF>", "<*>"}


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


def minimal_repr(pos: int, ref: str, alt: str):
    """Parsimonious (minimal) representation: trim common suffix then prefix.

    Does NOT reference-left-align (that needs the FASTA); the bcftools path
    handles left-shifting. Enough to match the DB for non-shifted indels and
    to canonicalize the anchor base."""
    ref = (ref or "").upper()
    alt = (alt or "").upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def alt_alleles(alt_field: str):
    """Split a possibly-multiallelic ALT into (index, allele) skipping symbolic."""
    out = []
    for j, a in enumerate((alt_field or "").split(",")):
        a = a.strip()
        if a in _SYMBOLIC:
            continue
        out.append((j, a))
    return out


# --------------------------------------------------------------------------
# tool / reference resolution
# --------------------------------------------------------------------------

def resolve_bcftools():
    sif = os.environ.get("BCFTOOLS_SIF")
    if sif and os.path.exists(sif):
        binds = os.environ.get("APPTAINER_BIND", "/home")
        return ["apptainer", "exec", "--bind", binds, sif, "bcftools"]
    b = os.environ.get("BCFTOOLS_BIN", "bcftools")
    return [b] if shutil.which(b) else None


def resolve_ref(explicit):
    for r in ([explicit] if explicit else []) + REF_CANDIDATES:
        if r and os.path.exists(r) and os.path.exists(r + ".fai"):
            return r
    return None


# --------------------------------------------------------------------------
# join implementations -> hits[(row_i, alt_j)] = (ac, an, af)
# --------------------------------------------------------------------------

def join_bcftools(rows, db, ref, bcftools):
    """norm -m- -f ref (split+left-align) | sort | annotate -a db. Robust: any
    failure raises so the caller can fall back to Python."""
    tmpd = tempfile.mkdtemp(prefix="inhouse_af.")
    mini = os.path.join(tmpd, "mini.vcf")
    try:
        # contig header from the .fai so bcftools has a contig order to sort by
        contigs = []
        with open(ref + ".fai") as f:
            for line in f:
                p = line.split("\t")
                if len(p) >= 2:
                    contigs.append(f"##contig=<ID={p[0]},length={p[1]}>")
        with open(mini, "w") as o:
            o.write("##fileformat=VCFv4.2\n")
            o.write("\n".join(contigs) + "\n")
            o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for i, r in enumerate(rows):
                chrom = norm_chrom(r.get("CHROM", ""))
                pos = (r.get("POS") or "").strip()
                ref_a = (r.get("REF") or "").strip().upper()
                if not (chrom and pos.isdigit() and ref_a):
                    continue
                for j, a in alt_alleles(r.get("ALT", "")):
                    o.write(f"{chrom}\t{pos}\t{i}_{j}\t{ref_a}\t{a.upper()}\t.\t.\t.\n")

        norm_gz = os.path.join(tmpd, "norm.vcf.gz")
        sort_gz = os.path.join(tmpd, "sorted.vcf.gz")
        subprocess.run(bcftools + ["norm", "-m-", "-f", ref, "--check-ref", "x",
                                   mini, "-Oz", "-o", norm_gz],
                       check=True, stderr=subprocess.DEVNULL)
        subprocess.run(bcftools + ["sort", "-T", tmpd, norm_gz, "-Oz", "-o", sort_gz],
                       check=True, stderr=subprocess.DEVNULL)
        # `bcftools annotate -a` uses the synced reader, which requires the MAIN
        # input to be an indexed file (even a stdin stream fails with "could not
        # load index"). Index the sorted file before annotating.
        subprocess.run(bcftools + ["index", "-t", sort_gz],
                       check=True, stderr=subprocess.DEVNULL)
        p = subprocess.Popen(
            bcftools + ["annotate", "-a", db,
                        "-c", "INFO/INHOUSE_AC,INFO/INHOUSE_AN,INFO/INHOUSE_AF", sort_gz],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        hits = {}
        for line in p.stdout:
            if not line or line[0] == "#":
                continue
            F = line.rstrip("\n").split("\t", 8)
            if len(F) < 8:
                continue
            rid = F[2]
            info = F[7]
            ac = info_get(info, "INHOUSE_AC")
            if ac is None:
                continue
            an = info_get(info, "INHOUSE_AN")
            af = info_get(info, "INHOUSE_AF")
            try:
                i_s, j_s = rid.split("_", 1)
                key = (int(i_s), int(j_s))
            except ValueError:
                continue
            hits[key] = (ac or ".", an or ".", af or ".")
        err = p.stderr.read()
        rc = p.wait()
        if rc != 0:
            raise RuntimeError(f"bcftools annotate rc={rc}: {err.strip()[:300]}")
        return hits
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)


def join_python(rows, db):
    """Pure-Python: minimal-repr keys + one streaming pass over the DB."""
    index: dict[tuple, list] = {}
    for i, r in enumerate(rows):
        chrom = norm_chrom(r.get("CHROM", ""))
        pos = (r.get("POS") or "").strip()
        ref_a = (r.get("REF") or "").strip()
        if not (chrom and pos.isdigit() and ref_a):
            continue
        for j, a in alt_alleles(r.get("ALT", "")):
            p, rr, aa = minimal_repr(int(pos), ref_a, a)
            index.setdefault((chrom, p, rr, aa), []).append((i, j))
    hits = {}
    opener = gzip.open if db.endswith(".gz") else open
    with opener(db, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            F = line.rstrip("\n").split("\t", 8)
            if len(F) < 8:
                continue
            try:
                key = (F[0], int(F[1]), F[3].upper(), F[4].upper())
            except ValueError:
                continue
            targets = index.get(key)
            if not targets:
                continue
            info = F[7]
            ac = info_get(info, "INHOUSE_AC") or "."
            an = info_get(info, "INHOUSE_AN") or "."
            af = info_get(info, "INHOUSE_AF") or "."
            for t in targets:
                hits[t] = (ac, an, af)
    return hits


# --------------------------------------------------------------------------

def assemble(rows, hits) -> int:
    """Write per-allele comma-joined columns from hits. Returns rows with >=1
    matched allele."""
    n_hit = 0
    for i, r in enumerate(rows):
        alleles = alt_alleles(r.get("ALT", ""))
        if not alleles:
            r[COL_AC] = r.get(COL_AC, "") or ""
            r[COL_AN] = r.get(COL_AN, "") or ""
            r[COL_AF] = r.get(COL_AF, "") or ""
            continue
        ac_p, an_p, af_p, any_hit = [], [], [], False
        for j, _a in alleles:
            v = hits.get((i, j))
            if v:
                any_hit = True
                ac_p.append(v[0]); an_p.append(v[1]); af_p.append(v[2])
            else:
                ac_p.append("."); an_p.append("."); af_p.append(".")
        if any_hit:
            r[COL_AC] = ",".join(ac_p)
            r[COL_AN] = ",".join(an_p)
            r[COL_AF] = ",".join(af_p)
            n_hit += 1
        else:
            r[COL_AC] = r[COL_AN] = r[COL_AF] = ""
    return n_hit


def annotate(tsv: Path, db: str, ref_arg: str | None) -> int:
    if not os.path.exists(db):
        print(f"[inhouse-af] DB not found: {db} — skipping (no-op)", file=sys.stderr)
        return 0

    with open(tsv, "r", encoding="utf-8", newline="") as fi:
        reader = csv.DictReader(fi, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            print(f"[inhouse-af] empty TSV: {tsv}", file=sys.stderr)
            return 0
        rows = list(reader)
    for col in (COL_AC, COL_AN, COL_AF):
        if col not in fieldnames:
            fieldnames.append(col)

    bcftools = resolve_bcftools()
    ref = resolve_ref(ref_arg)
    used = "python"
    hits = None
    if bcftools and ref:
        try:
            hits = join_bcftools(rows, db, ref, bcftools)
            used = "bcftools"
        except Exception as e:  # noqa: BLE001 — degrade, never block stop-gaps
            print(f"[inhouse-af] bcftools path failed ({e}); falling back to python",
                  file=sys.stderr)
            hits = None
    if hits is None:
        hits = join_python(rows, db)

    n_hit = assemble(rows, hits)

    fd, tmp_name = tempfile.mkstemp(dir=str(tsv.parent), prefix=tsv.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fo:
            writer = csv.DictWriter(fo, fieldnames=fieldnames, delimiter="\t",
                                    extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp_name, tsv)
    except BaseException:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        raise

    print(f"[inhouse-af] {len(rows)} variants, {n_hit} matched in-house AF DB "
          f"(join={used}{', ref='+os.path.basename(ref) if used=='bcftools' else ''})")
    return 0


# --------------------------------------------------------------------------

def selftest() -> int:
    assert norm_chrom("1") == "chr1"
    assert minimal_repr(45330228, "CAA", "C") == (45330228, "CAA", "C")
    assert minimal_repr(45330228, "CAA", "CA") == (45330228, "CA", "C")
    assert minimal_repr(100, "AT", "AG") == (101, "T", "G")   # SNV inside
    assert alt_alleles("C,CA") == [(0, "C"), (1, "CA")]
    assert alt_alleles("A,*,G") == [(0, "A"), (2, "G")]       # skip spanning-del *
    assert info_get("INHOUSE_AC=5;INHOUSE_AF=0.05", "INHOUSE_AC") == "5"

    # end-to-end via the python path on a multiallelic row
    import io
    d = tempfile.mkdtemp()
    tsv = Path(d) / "snv.tsv"
    tsv.write_text(
        "CHROM\tPOS\tREF\tALT\tGENE\n"
        "chr1\t45330228\tCAA\tC,CA\tMUTYH\n"   # both alleles in DB
        "chr1\t200\tA\tG\tTP53\n"              # not in DB
        "chr7\t500\tAT\tA\tBRCA1\n",           # deletion, in DB after trim
        encoding="utf-8")
    db = Path(d) / "inhouse.vcf.gz"
    with gzip.open(db, "wt", encoding="utf-8") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("chr1\t45330228\t.\tCAA\tC\t.\t.\tINHOUSE_AC=595;INHOUSE_AN=1354;INHOUSE_AF=0.439439\n")
        f.write("chr1\t45330228\t.\tCA\tC\t.\t.\tINHOUSE_AC=571;INHOUSE_AN=1354;INHOUSE_AF=0.421713\n")
        f.write("chr7\t500\t.\tAT\tA\t.\t.\tINHOUSE_AC=10;INHOUSE_AN=1000;INHOUSE_AF=0.01\n")
    annotate(tsv, str(db), ref_arg=None)   # no ref -> python path
    out = list(csv.DictReader(io.StringIO(tsv.read_text()), delimiter="\t"))
    assert out[0]["INHOUSE_AF"] == "0.439439,0.421713", out[0]
    assert out[0]["INHOUSE_AC"] == "595,571", out[0]
    assert out[1]["INHOUSE_AF"] == "", out[1]          # miss -> clean blank
    assert out[2]["INHOUSE_AF"] == "0.01", out[2]
    assert annotate(tsv, str(Path(d) / "nope.vcf.gz"), None) == 0   # missing DB = no-op
    print("selftest OK")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", type=Path, help="snv_indel.annotated.tsv to annotate in place")
    ap.add_argument("--db", default=DEFAULT_DB, help=f"in-house AF sites VCF (default {DEFAULT_DB})")
    ap.add_argument("--ref", default=None, help="hg38 FASTA (+.fai) for bcftools normalize; "
                    "auto-detected from NGS_UI_INHOUSE_AF_REF / common paths if omitted")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    if not args.tsv:
        ap.error("--tsv required (or --selftest)")
    if not args.tsv.is_file():
        raise SystemExit(f"--tsv not found: {args.tsv}")
    return annotate(args.tsv, args.db, args.ref)


if __name__ == "__main__":
    raise SystemExit(main())
