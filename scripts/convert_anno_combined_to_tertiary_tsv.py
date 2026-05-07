#!/usr/bin/env python3
"""Convert anno_combined.txt.gz (legacy R annotation table) to
snv_indel.annotated.tsv per 三級輸出計畫.md spec.

The source has one row per (variant × transcript) so we dedupe by
(CHROM,POS,REF,ALT) keeping the row with the best transcript:

    MANE_SELECT  >  MANE_PLUS_CLINICAL  >  CANONICAL  >  any

Columns the spec expects but the source doesn't carry (TWB_AF,
PKNN_LLR, ESM2/Evo2, LOFTEE_*, PHASE_*, IN_ROH/IN_PANEL/IN_BLACKLIST,
…) are emitted blank or with a sensible placeholder so the adapter
keeps loading. This script is a stop-gap until the real Phase 4
tertiary pipeline lands; for now it lets reviewers eyeball the UI on
real-shape data.

Usage:
    python3 scripts/convert_anno_combined_to_tertiary_tsv.py \\
        --in  26WE0040_anno_combined.txt.gz \\
        --out tertiary_output/26WE0040/snv_indel.annotated.tsv
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path


OUTPUT_COLUMNS = [
    "CHROM", "POS", "REF", "ALT",
    "GENE", "TRANSCRIPT", "TRANSCRIPT_TYPE",
    "HGVS_C", "HGVS_P", "CONSEQUENCE",
    "MANE_ALL",
    "CALLERS",
    "ZYGOSITY", "GT_DV", "GT_HC",
    "EXON", "INTRON",
    "AD", "VAF",
    "CLINVAR_SIG", "CLINVAR_STARS", "CLINVAR_DN", "CLINVAR_CONF",
    "GNOMAD_G_AF", "GNOMAD_G_EAS_AF", "GNOMAD_E_AF", "GNOMAD_E_EAS_AF",
    "TWB_AF",
    "PKNN_LLR",
    "REVEL", "BAYESDEL", "ALPHAMISSENSE", "METARNN",
    "ESM2_SCORE", "EVO2_SCORE",
    "SPLICEAI_MAX", "CADD_PHRED",
    "LOFTEE_HC", "LOFTEE_FILTER", "LOFTEE_FLAGS",
    "ACMG_EVIDENCE", "ACMG_POINTS", "ACMG_CLASS",
    "PHASE_GROUP", "PHASE_RESULT",
    "IN_ROH", "IN_PANEL", "IN_BLACKLIST",
    "OMIM_LINK", "GNOMAD_LINK", "CLINVAR_LINK",
    "REPORT_CLASS",
]

# R writes "NA" for missing values; treat blank too. Keep the helper
# narrow so we never accidentally swallow a literal "NA" gene symbol.
def _na(v):
    if v is None:
        return ""
    s = str(v).strip()
    if s in ("", "NA", ".", "<NA>"):
        return ""
    return s


def _hgvs_strip(hgvs: str, prefix_re: str) -> str:
    s = _na(hgvs)
    if not s:
        return ""
    m = re.match(prefix_re, s)
    return m.group(1) if m else s


HGVS_C_RE = re.compile(r"^[^:]+:(c\..*)$")
HGVS_P_RE = re.compile(r"^[^:]+:(p\..*)$")


def _zyg_to_gt(z: str) -> str:
    z = (_na(z) or "").lower()
    return {"het": "0/1", "hom": "1/1", "hemi": "1"}.get(z, "./.")


# Lower number = better. Drives both TRANSCRIPT_TYPE and the per-variant
# dedupe across transcripts.
def _tx_priority(row: dict) -> int:
    mane = (row.get("MANE") or "").strip()
    if "MANE_Select" in mane:
        return 0
    if "MANE_Plus_Clinical" in mane:
        return 1
    if (row.get("CANONICAL") or "").upper() == "YES":
        return 2
    return 9


def _tx_type(row: dict) -> str:
    p = _tx_priority(row)
    return {
        0: "MANE_SELECT",
        1: "MANE_PLUS_CLINICAL",
        2: "CANONICAL",
    }.get(p, "")


# Stub LOFTEE_HC: HC for high-confidence LoF consequences. Replace with
# the actual VEP+LOFTEE output once Phase 4 wires that in.
LOF_TERMS = {
    "stop_gained", "frameshift_variant", "splice_acceptor_variant",
    "splice_donor_variant", "stop_lost", "start_lost",
    "transcript_ablation",
}


def _gnomad_link(chrom, pos, ref, alt, build="hg38") -> str:
    if not all([chrom, pos, ref, alt]):
        return ""
    c = str(chrom).replace("chr", "")
    dataset = "gnomad_r4" if build == "hg38" else "gnomad_r2_1"
    return (f"https://gnomad.broadinstitute.org/variant/"
            f"{c}-{pos}-{ref}-{alt}?dataset={dataset}")


def _clinvar_link(chrom, pos) -> str:
    if not chrom or not pos:
        return ""
    c = str(chrom).replace("chr", "")
    return (f"https://www.ncbi.nlm.nih.gov/clinvar/?term="
            f"{c}%5BCHR%5D+AND+{pos}%5BCHRPOS%5D")


def _omim_link(omim_id) -> str:
    s = _na(omim_id)
    if not s:
        return ""
    first = re.split(r"[,;\s]+", s)[0]
    return f"https://www.omim.org/entry/{first}"


def _acmg_evidence_pipe(criteria) -> str:
    s = _na(criteria)
    if not s:
        return ""
    parts = [p.strip() for p in re.split(r"[,;|]", s) if p.strip()]
    return "|".join(parts)


def convert_row(row: dict, build: str = "hg38") -> dict:
    chrom = _na(row.get("CHROM"))
    pos   = _na(row.get("POS"))
    ref   = _na(row.get("REF"))
    alt   = _na(row.get("ALT"))
    consequence = _na(row.get("Consequence"))
    transcript  = _na(row.get("Feature"))
    return {
        "CHROM": chrom,
        "POS":   pos,
        "REF":   ref,
        "ALT":   alt,
        "GENE":              _na(row.get("SYMBOL")) or _na(row.get("gene_symbol")),
        "TRANSCRIPT":        transcript,
        "TRANSCRIPT_TYPE":   _tx_type(row),
        "HGVS_C":            _hgvs_strip(row.get("HGVSc"), HGVS_C_RE.pattern),
        "HGVS_P":            _hgvs_strip(row.get("HGVSp"), HGVS_P_RE.pattern),
        "CONSEQUENCE":       consequence,
        "MANE_ALL":          "[]",
        "CALLERS":           "DV+HC",
        "ZYGOSITY":          _na(row.get("zygosity")),
        "GT_DV":             _zyg_to_gt(row.get("zygosity")),
        "GT_HC":             _zyg_to_gt(row.get("zygosity")),
        "EXON":              _na(row.get("EXON")),
        "INTRON":            _na(row.get("INTRON")),
        "AD":                _na(row.get("AD")),
        "VAF":               _na(row.get("alt_af")),
        "CLINVAR_SIG":       _na(row.get("CLNSIG_20260503"))
                              or _na(row.get("CLNSIG")),
        "CLINVAR_STARS":     _na(row.get("clinvar_stars")),
        "CLINVAR_DN":        "",
        "CLINVAR_CONF":      _na(row.get("CLNSIGCONF_20260503"))
                              or _na(row.get("CLNSIGCONF")),
        "GNOMAD_G_AF":       _na(row.get("gnomad41_genome_AF")),
        "GNOMAD_G_EAS_AF":   _na(row.get("gnomad41_genome_AF_eas")),
        "GNOMAD_E_AF":       _na(row.get("gnomad41_exome_AF")),
        "GNOMAD_E_EAS_AF":   _na(row.get("gnomad41_exome_AF_eas")),
        "TWB_AF":            "",
        "PKNN_LLR":          "",
        "REVEL":             _na(row.get("REVEL_score")),
        "BAYESDEL":          _na(row.get("BayesDel_noAF_score")),
        "ALPHAMISSENSE":     _na(row.get("AlphaMissense_score")),
        "METARNN":           _na(row.get("MetaRNN_score")),
        "ESM2_SCORE":        "",
        "EVO2_SCORE":        "",
        "SPLICEAI_MAX":      _na(row.get("SpliceAI_score")),
        "CADD_PHRED":        _na(row.get("CADD_phred")),
        "LOFTEE_HC":         "HC" if consequence in LOF_TERMS else "",
        "LOFTEE_FILTER":     "",
        "LOFTEE_FLAGS":      "",
        "ACMG_EVIDENCE":     _acmg_evidence_pipe(row.get("ACMG_criteria")),
        "ACMG_POINTS":       _na(row.get("ACMG_score")),
        "ACMG_CLASS":        _na(row.get("ACMG_classification")),
        "PHASE_GROUP":       "",
        "PHASE_RESULT":      "unphased",
        "IN_ROH":            "false",
        "IN_PANEL":          "false",
        "IN_BLACKLIST":      "false",
        "OMIM_LINK":         _omim_link(row.get("OMIM_id")),
        "GNOMAD_LINK":       _gnomad_link(chrom, pos, ref, alt, build=build),
        "CLINVAR_LINK":      _clinvar_link(chrom, pos),
        "REPORT_CLASS":      "",
    }


def _sort_key(row: dict):
    c = str(_na(row.get("CHROM"))).replace("chr", "").upper()
    chrom_num = {"X": 23, "Y": 24, "MT": 25, "M": 25}.get(c)
    if chrom_num is None:
        try:
            chrom_num = int(c)
        except ValueError:
            chrom_num = 99
    try:
        p = int(_na(row.get("POS")) or 0)
    except ValueError:
        p = 0
    return (chrom_num, p)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in",  dest="inp", required=True, help="anno_combined.txt[.gz]")
    ap.add_argument("--out", dest="out", required=True, help="snv_indel.annotated.tsv")
    ap.add_argument("--build", default="hg38", choices=["hg19", "hg38"])
    args = ap.parse_args()

    opener = gzip.open if str(args.inp).endswith(".gz") else open
    best: dict[tuple, tuple[int, dict]] = {}
    n_rows = 0
    with opener(args.inp, "rt", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            n_rows += 1
            key = (_na(row.get("CHROM")), _na(row.get("POS")),
                   _na(row.get("REF")),   _na(row.get("ALT")))
            if not all(key):
                continue
            prio = _tx_priority(row)
            existing = best.get(key)
            if existing is None or prio < existing[0]:
                best[key] = (prio, row)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted((r for _, r in best.values()), key=_sort_key)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(convert_row(r, build=args.build))

    print(f"Read {n_rows} (variant×transcript) rows from {args.inp}")
    print(f"Wrote {len(rows)} unique variants × {len(OUTPUT_COLUMNS)} cols → {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
