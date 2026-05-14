#!/usr/bin/env python3
"""Backfill CLINVAR_* columns from a local ClinVar VCF.

The new tertiary pipeline (Phase 1) emits all CLINVAR_* columns as
'.' since the ClinVar annotation step isn't wired yet. This script
reads the same ClinVar VCF the pipeline already references
(clinvar_<date>.vcf.gz under <ref>/tertiary/clinvar/), joins by
(CHROM, POS, REF, ALT), and fills:
    CLINVAR_SIG       INFO/CLNSIG       (Pathogenic / Likely_pathogenic / ...)
    CLINVAR_STARS     review-status star count (0-4)
    CLINVAR_DN        INFO/CLNDN        disease names
    CLINVAR_SIGCONF   INFO/CLNSIGCONF

Pipeline-supplied non-empty cells are preserved (only '.' / blank are
filled). Stop running this once the pipeline's ClinVar step ships;
the join is a no-op when columns are already populated.

CLINVAR_STARS mapping (NCBI review-status hierarchy):
  practice_guideline                                                    → 4
  reviewed_by_expert_panel                                              → 3
  criteria_provided,_multiple_submitters,_no_conflicts                  → 2
  criteria_provided,_single_submitter / _conflicting_*                  → 1
  no_assertion_* / no_classification_provided                           → 0

Chromosome prefix mismatches ('chr1' vs '1') are handled both ways.
Assumes both VCFs are left-aligned + parsimonious (bcftools norm). If
they aren't, indel matches will silently miss.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
from pathlib import Path

_STAR_MAP = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "criteria_provided,_single_submitter": 1,
}


def _stars(revstat: str) -> int:
    return _STAR_MAP.get((revstat or "").strip(), 0)


def _info_field(info: str, key: str) -> str:
    """Extract INFO/<key>=value (URL-decoded)."""
    m = re.search(r"(?:^|;)" + re.escape(key) + r"=([^;]+)", info)
    if not m:
        return ""
    return m.group(1).replace("%2C", ",").replace("%3D", "=")


def _open_vcf(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") \
        if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def _norm_chr(c: str) -> str:
    """Strip leading 'chr' so we can match either prefix style."""
    c = c.strip()
    return c[3:] if c.startswith("chr") else c


def _has_value(row: dict, col: str) -> bool:
    v = (row.get(col) or "").strip()
    return v not in ("", ".", "NA", "N/A")


def collect_keys(tsv_path: Path) -> set:
    """First pass over TSV → set of (chr_no_prefix, pos, ref, alt)."""
    keys: set = set()
    with open(tsv_path, "r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            chrom = _norm_chr(row.get("CHROM", ""))
            pos   = (row.get("POS") or "").strip()
            ref   = (row.get("REF") or "").strip()
            alt   = (row.get("ALT") or "").strip()
            if not (chrom and pos and ref and alt):
                continue
            # ALT in the pipeline TSV may still carry multi-allelics like
            # "AC,C" — split so each sub-allele can match separately.
            for a in alt.split(","):
                keys.add((chrom, pos, ref, a))
    return keys


def index_clinvar(vcf_path: Path, want: set) -> dict:
    """Stream ClinVar VCF, keep only entries whose key is in `want`.

    Returns {(chr_no_prefix, pos, ref, alt): (sig, stars, dn, sigconf)}.
    """
    out: dict = {}
    n_seen = 0
    n_match = 0
    with _open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _, ref, alt, _, _, info = parts[:8]
            chrom_n = _norm_chr(chrom)
            alts = alt.split(",")
            matched = [k for k in
                       ((chrom_n, pos, ref, a) for a in alts)
                       if k in want]
            n_seen += 1
            if not matched:
                continue
            sig = _info_field(info, "CLNSIG")
            if not sig:
                continue
            stars   = _stars(_info_field(info, "CLNREVSTAT"))
            dn      = _info_field(info, "CLNDN")
            sigconf = _info_field(info, "CLNSIGCONF")
            for k in matched:
                out[k] = (sig, stars, dn, sigconf)
                n_match += 1
            if n_seen % 500000 == 0:
                print(f"[clinvar] scanned {n_seen} VCF rows; "
                      f"matched {n_match}", file=sys.stderr)
    print(f"[clinvar] scanned {n_seen} VCF rows total; "
          f"matched {n_match}", file=sys.stderr)
    return out


def merge_into_tsv(in_tsv: Path, out_tsv: Path, cv: dict) -> dict:
    overwriting = in_tsv.resolve() == out_tsv.resolve()
    target = Path(str(out_tsv) + ".tmp") if overwriting else out_tsv
    target.parent.mkdir(parents=True, exist_ok=True)

    n_in = 0
    n_filled = 0

    with open(in_tsv, "r", encoding="utf-8", newline="") as fi:
        reader = csv.DictReader(fi, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for col in ("CLINVAR_SIG", "CLINVAR_STARS",
                    "CLINVAR_DN", "CLINVAR_SIGCONF"):
            if col not in fieldnames:
                fieldnames.append(col)
        with open(target, "w", encoding="utf-8", newline="") as fo:
            writer = csv.DictWriter(fo, fieldnames=fieldnames, delimiter="\t",
                                    extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                n_in += 1
                chrom = _norm_chr(row.get("CHROM", ""))
                pos   = (row.get("POS") or "").strip()
                ref   = (row.get("REF") or "").strip()
                alts  = [a for a in (row.get("ALT") or "").split(",") if a]
                hit = None
                for a in alts:
                    h = cv.get((chrom, pos, ref, a))
                    if h is not None:
                        hit = h
                        break
                if hit:
                    sig, stars, dn, sigconf = hit
                    # CLINVAR_SIG drives the group: if the pipeline left
                    # SIG blank ('.' / ''), it hasn't done ClinVar
                    # annotation at all and the other three cells are
                    # placeholders (raw TSV defaults STARS to "0",
                    # SIGCONF to "." etc.). Treat the four columns as
                    # one unit and overwrite them all from ClinVar.
                    if not _has_value(row, "CLINVAR_SIG") and sig:
                        row["CLINVAR_SIG"]     = sig
                        row["CLINVAR_STARS"]   = str(stars)
                        row["CLINVAR_DN"]      = dn or row.get("CLINVAR_DN", "")
                        row["CLINVAR_SIGCONF"] = sigconf or row.get("CLINVAR_SIGCONF", "")
                        n_filled += 1
                writer.writerow(row)
    if overwriting:
        os.replace(target, out_tsv)
    return {"n_in": n_in, "n_filled": n_filled}


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True,
                    help="snv_indel.annotated.tsv (updated in place unless "
                         "--out-tsv given)")
    ap.add_argument("--out-tsv",
                    help="write merged TSV here instead of overwriting --tsv")
    ap.add_argument("--clinvar",
                    default=os.environ.get(
                        "CLINVAR_VCF",
                        "/home/pipeline/reference/hg38/tertiary/clinvar/"
                        "clinvar_20260418.vcf.gz"),
                    help="ClinVar VCF[.gz] (default: pipeline's reference; "
                         "or env CLINVAR_VCF)")
    args = ap.parse_args()

    in_tsv = Path(args.tsv).resolve()
    if not in_tsv.is_file():
        print(f"ERROR: --tsv 找不到：{in_tsv}", file=sys.stderr)
        return 2
    out_tsv = Path(args.out_tsv).resolve() if args.out_tsv else in_tsv
    clinvar = Path(args.clinvar)
    if not clinvar.is_file():
        print(f"ERROR: --clinvar 找不到：{clinvar}", file=sys.stderr)
        return 2

    print(f"[clinvar] tsv : {in_tsv}", file=sys.stderr)
    print(f"[clinvar] out : {out_tsv}", file=sys.stderr)
    print(f"[clinvar] db  : {clinvar}", file=sys.stderr)

    print(f"[clinvar] pass 1/2: collecting keys from TSV ...", file=sys.stderr)
    want = collect_keys(in_tsv)
    print(f"[clinvar]          {len(want)} unique sites in TSV",
          file=sys.stderr)

    print(f"[clinvar] pass 2/2: scanning ClinVar VCF ...", file=sys.stderr)
    cv = index_clinvar(clinvar, want)

    stats = merge_into_tsv(in_tsv, out_tsv, cv)
    print(f"[clinvar] backfilled CLINVAR_* for {stats['n_filled']}/"
          f"{stats['n_in']} TSV rows  (matched {len(cv)} sites)",
          file=sys.stderr)
    print(f"[clinvar] done → {out_tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
