#!/usr/bin/env python3
"""Annotate snv_indel.annotated.tsv with in-house allele frequency.

Adds three columns from the in-house AF sites VCF (built by
scripts/inhouse_af/publish_af.py):

    INHOUSE_AC   in-house alt allele count
    INHOUSE_AN   in-house total called alleles (DP>=10, ploidy-weighted)
    INHOUSE_AF   INHOUSE_AC / INHOUSE_AN

The SNV adapter surfaces these next to the gnomAD AF; the card shows a
``AF_nckuh`` row as ``<AF> (AC/AN)`` so the reviewer can see local frequency
(incl. rare variants) with the cohort size visible.

Matching is a single streaming pass over the sites VCF (NOT per-variant tabix
seeks): read the TSV keys into memory, stream the sorted sites VCF once, and
fill the columns for keys that match on the normalized `(chrom,pos,ref,alt)`.
Same representation on both sides (`bcftools norm -m-` left-aligned); complex
indels whose representation disagrees simply don't match (same limitation as
the GeneBe join) — they render blank, never wrong.

Fill-or-augment, idempotent, atomic replace. **No-op (exit 0) when the DB is
missing**, so it is safe to wire into run_stopgaps.sh unconditionally.

Usage:
    scripts/annotate_inhouse_af.py \\
        --tsv tertiary_output/<SID>/snv_indel.annotated.tsv \\
        [--db  $NGS_UI_HOME/biotools/inhouse_af/inhouse_af.hg38.vcf.gz]
    scripts/annotate_inhouse_af.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
import tempfile
from pathlib import Path

DEFAULT_DB = os.environ.get(
    "NGS_UI_INHOUSE_AF_DB",
    str(Path.home() / "NGS_UI" / "biotools" / "inhouse_af" / "inhouse_af.hg38.vcf.gz"),
)
COL_AC, COL_AN, COL_AF = "INHOUSE_AC", "INHOUSE_AN", "INHOUSE_AF"


def norm_chrom(chrom: str) -> str:
    """Both TSV and sites VCF use `chrN`; normalize defensively either way."""
    c = (chrom or "").strip()
    if not c:
        return ""
    return c if c.lower().startswith("chr") else "chr" + c


def variant_key(chrom: str, pos: str, ref: str, alt: str):
    return (norm_chrom(chrom), (pos or "").strip(),
            (ref or "").strip().upper(), (alt or "").strip().upper())


def info_get(info: str, key: str):
    """Fetch a single INFO value (key=val;...). None if absent."""
    kv = key + "="
    for field in info.split(";"):
        if field.startswith(kv):
            return field[len(kv):]
    return None


def annotate(tsv: Path, db: str) -> int:
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

    # index the TSV variants: key -> [row indices] (usually one row per variant)
    index: dict[tuple, list[int]] = {}
    for i, r in enumerate(rows):
        k = variant_key(r.get("CHROM", ""), r.get("POS", ""),
                        r.get("REF", ""), r.get("ALT", ""))
        index.setdefault(k, []).append(i)
        # start blank so the columns exist on every row even without a hit
        r[COL_AC] = r.get(COL_AC, "") or ""
        r[COL_AN] = r.get(COL_AN, "") or ""
        r[COL_AF] = r.get(COL_AF, "") or ""

    for col in (COL_AC, COL_AN, COL_AF):
        if col not in fieldnames:
            fieldnames.append(col)

    # single streaming pass over the sorted sites VCF
    n_hit = 0
    opener = gzip.open if db.endswith(".gz") else open
    with opener(db, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            # fast field slice: CHROM POS ID REF ALT QUAL FILTER INFO
            F = line.rstrip("\n").split("\t", 8)
            if len(F) < 8:
                continue
            k = variant_key(F[0], F[1], F[3], F[4])
            idxs = index.get(k)
            if not idxs:
                continue
            info = F[7]
            ac = info_get(info, "INHOUSE_AC")
            an = info_get(info, "INHOUSE_AN")
            af = info_get(info, "INHOUSE_AF")
            for i in idxs:
                rows[i][COL_AC] = ac or ""
                rows[i][COL_AN] = an or ""
                rows[i][COL_AF] = af or ""
            n_hit += len(idxs)

    # atomic replace
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

    print(f"[inhouse-af] {len(rows)} variants, {n_hit} matched in-house AF DB")
    return 0


def selftest() -> int:
    assert norm_chrom("1") == "chr1"
    assert norm_chrom("chr1") == "chr1"
    assert norm_chrom("chrM") == "chrM"
    assert variant_key("1", "100", "a", "t") == ("chr1", "100", "A", "T")
    assert info_get("INHOUSE_AC=61;INHOUSE_AN=1352;INHOUSE_AF=0.045", "INHOUSE_AN") == "1352"
    assert info_get("X=1;INHOUSE_AF=0.5", "INHOUSE_AC") is None

    # end-to-end on a tiny TSV + tiny sites VCF
    import io
    d = tempfile.mkdtemp()
    tsv = Path(d) / "snv.tsv"
    tsv.write_text(
        "CHROM\tPOS\tREF\tALT\tGENE\n"
        "chr1\t100\tA\tG\tBRCA2\n"      # will match
        "chr1\t200\tAC\tA\tTP53\n"      # no match
        , encoding="utf-8")
    db = Path(d) / "inhouse.vcf.gz"
    with gzip.open(db, "wt", encoding="utf-8") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("chr1\t100\t.\tA\tG\t.\t.\tINHOUSE_AC=5;INHOUSE_AN=100;INHOUSE_AF=0.05\n")
    annotate(tsv, str(db))
    out = list(csv.DictReader(io.StringIO(tsv.read_text()), delimiter="\t"))
    assert out[0]["INHOUSE_AF"] == "0.05" and out[0]["INHOUSE_AN"] == "100", out[0]
    assert out[1]["INHOUSE_AF"] == "", out[1]
    # missing DB is a no-op (columns from the previous run stay, exit 0)
    assert annotate(tsv, str(Path(d) / "nope.vcf.gz")) == 0
    print("selftest OK")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", type=Path, help="snv_indel.annotated.tsv to annotate in place")
    ap.add_argument("--db", default=DEFAULT_DB,
                    help=f"in-house AF sites VCF (default {DEFAULT_DB})")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    if not args.tsv:
        ap.error("--tsv required (or --selftest)")
    if not args.tsv.is_file():
        raise SystemExit(f"--tsv not found: {args.tsv}")
    return annotate(args.tsv, args.db)


if __name__ == "__main__":
    raise SystemExit(main())
