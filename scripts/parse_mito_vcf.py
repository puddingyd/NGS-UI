#!/usr/bin/env python3
"""parse_mito_vcf.py — GATK Mutect2-mito VCF + MITOMAP tables → mito.annotated.tsv

MITOMAP-only annotation — no VEP. For a 16.5 kb genome where the card
only ever shows MITOMAP-recorded (or MitoTIP-pathogenic) variants,
VEP's consequence/HGVS calls add little: the gene/locus comes from a
hard-coded rCRS map, the m. HGVS for SNVs is trivial, and the protein
change comes straight from MITOMAP's "Amino Acid Change" column. So
this script reads the Mutect2-mito VCF directly (gzip), needs only
the bundled-by-the-pipeline MITOMAP tables, and is pure Python.

Input:
  --vcf          <sample>.mito.vcf.gz  (GATK Mutect2 --mitochondria-mode,
                 FILTER applied by FilterMutectCalls; chrM / rCRS coords;
                 FORMAT carries AF=heteroplasmy, AD, DP)
  --mitomap_cc   mitomap_mutations_coding_control.tsv
  --mitomap_rna  mitomap_mutations_rna.tsv
  --sample_id    sample id (the FORMAT-column sample name; auto-detected if omitted)
  --output       mito.annotated.tsv

Output columns:
  CHROM POS REF ALT HGVS_M GENE LOCUS_TYPE CONSEQUENCE AA_CHANGE
  HETEROPLASMY AD DEPTH FILTER TLOD
  MITOMAP_DISEASE MITOMAP_STATUS MITOMAP_PLASMY MITOMAP_GB_FREQ MITOMAP_GB_SEQS
  MITOMAP_REFS MITOTIP_SCORE MITOMAP_ALLELE

Notes:
  * HGVS_M: SNV → m.{pos}{ref}>{alt}; simple del/ins/dup get a basic
    m.{...}del / m.{...}_{...}ins{seq} / m.{...}dup form (no HGVS
    3'-shifting — cosmetic for indels, which never join MITOMAP here).
  * GENE / LOCUS_TYPE from the rCRS gene-coordinate table (D-loop →
    MT-CR/control, the OriL gap → MT-OLR, other gaps → intergenic).
  * CONSEQUENCE: non-coding for tRNA/rRNA/control; for protein-coding
    it's derived from MITOMAP's "Amino Acid Change" (missense /
    synonymous / stop_gained / coding-other); AA_CHANGE = that column
    verbatim (e.g. "A52T").
  * MITOMAP joined by exact (POS, REF, ALT) only (cc keyed off
    "Nucleotide Change", rna off the <ref><pos><alt> "Allele"). No
    POS-only fallback — wrong-allele bleed is worse than a blank.
  * Multiallelic sites are split (one row per ALT). Records at the
    same (POS,REF,ALT) are de-duplicated, keeping the highest-TLOD copy.
  * The mito adapter on the server side then keeps FILTER=PASS only
    and only disease-relevant variants — this script emits everything.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path


# rCRS (NC_012920.1) gene map — start/end are 1-based inclusive.
# (name, start, end, locus_type). Control region split into two ranges.
_MT_GENES = [
    ("MT-CR",   16024, 16569, "control"),
    ("MT-CR",       1,   576, "control"),
    ("MT-TF",     577,   647, "tRNA"),
    ("MT-RNR1",   648,  1601, "rRNA"),
    ("MT-TV",    1602,  1670, "tRNA"),
    ("MT-RNR2",  1671,  3229, "rRNA"),
    ("MT-TL1",   3230,  3304, "tRNA"),
    ("MT-ND1",   3307,  4262, "protein"),
    ("MT-TI",    4263,  4331, "tRNA"),
    ("MT-TQ",    4329,  4400, "tRNA"),
    ("MT-TM",    4402,  4469, "tRNA"),
    ("MT-ND2",   4470,  5511, "protein"),
    ("MT-TW",    5512,  5579, "tRNA"),
    ("MT-TA",    5587,  5655, "tRNA"),
    ("MT-TN",    5657,  5729, "tRNA"),
    ("MT-OLR",   5730,  5760, "control"),   # origin of L-strand replication (OriL)
    ("MT-TC",    5761,  5826, "tRNA"),
    ("MT-TY",    5826,  5891, "tRNA"),
    ("MT-CO1",   5904,  7445, "protein"),
    ("MT-TS1",   7446,  7514, "tRNA"),
    ("MT-TD",    7518,  7585, "tRNA"),
    ("MT-CO2",   7586,  8269, "protein"),
    ("MT-TK",    8295,  8364, "tRNA"),
    ("MT-ATP8",  8366,  8572, "protein"),
    ("MT-ATP6",  8527,  9207, "protein"),
    ("MT-CO3",   9207,  9990, "protein"),
    ("MT-TG",    9991, 10058, "tRNA"),
    ("MT-ND3",  10059, 10404, "protein"),
    ("MT-TR",   10405, 10469, "tRNA"),
    ("MT-ND4L", 10470, 10766, "protein"),
    ("MT-ND4",  10760, 12137, "protein"),
    ("MT-TH",   12138, 12206, "tRNA"),
    ("MT-TS2",  12207, 12265, "tRNA"),
    ("MT-TL2",  12266, 12336, "tRNA"),
    ("MT-ND5",  12337, 14148, "protein"),
    ("MT-ND6",  14149, 14673, "protein"),
    ("MT-TE",   14674, 14742, "tRNA"),
    ("MT-CYB",  14747, 15887, "protein"),
    ("MT-TT",   15888, 15953, "tRNA"),
    ("MT-TP",   15956, 16023, "tRNA"),
]

_RNA_ALLELE_RE = re.compile(r"^([ACGTN])(\d+)([ACGTN])$", re.I)
# MITOMAP "Amino Acid Change" like "A52T", "G29S", "L329P".
_AA_RE = re.compile(r"^([A-Z\*])(\d+)([A-Z\*])$")
_STOP_AA = {"*", "X", "Term", "Stop"}


def _gene_at(pos: int) -> tuple[str, str]:
    """(gene, locus_type) for a position, from the rCRS map. D-loop →
    control region; the few tiny intergenic gaps → ('', 'intergenic')."""
    for name, start, end, lt in _MT_GENES:
        if start <= pos <= end:
            return name, lt
    return ("MT-CR", "control") if (pos <= 576 or pos >= 16024) else ("", "intergenic")


def _hgvs_m(pos: int, ref: str, alt: str) -> str:
    """m. HGVS for a chrM variant. SNV → m.{pos}{ref}>{alt}. Simple
    indels get a basic (non-3'-shifted) del / ins / dup form."""
    if len(ref) == 1 and len(alt) == 1:
        return f"m.{pos}{ref}>{alt}"
    # VCF anchors indels at the preceding base, so the changed bases
    # start at pos+1 for a deletion / are inserted after pos.
    if len(ref) > len(alt) and alt and ref.startswith(alt):
        del_seq = ref[len(alt):]
        s = pos + len(alt)
        e = s + len(del_seq) - 1
        return f"m.{s}del" if s == e else f"m.{s}_{e}del"
    if len(alt) > len(ref) and ref and alt.startswith(ref):
        ins_seq = alt[len(ref):]
        # naive dup detection: inserted run equals the preceding base(s)
        a = pos
        b = pos + len(ins_seq) - 1
        return f"m.{a}_{b}dup" if a != b else f"m.{a}dup"
    # mixed / complex — fall back to a delins-ish form
    return f"m.{pos}{ref}>{alt}"


def _to_float(s):
    if s is None:
        return None
    s = str(s).strip()
    if not s or s.upper() in ("NA", "N/A", "."):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _to_int(s):
    f = _to_float(s)
    if f is None:
        return None
    try:
        return int(f)
    except (TypeError, ValueError):
        return None


# ---- MITOMAP table loaders ------------------------------------------
# MITOMAP exports are Latin-1, not UTF-8 (stray 0xa0 etc.).

def _load_mitomap_cc(path: Path) -> dict[tuple[int, str, str], dict]:
    """coding/control table → {(pos, ref, alt): record}, keyed off the
    "Nucleotide Change" (ref-alt) column for SNVs. Also carries the
    "Amino Acid Change" column."""
    out: dict[tuple, dict] = {}
    if not path or not path.is_file():
        return out
    with path.open("r", encoding="latin-1", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                pos = int((row.get("Position") or "").strip())
            except ValueError:
                continue
            nuc = (row.get("Nucleotide Change") or "").strip()  # e.g. "T-C"
            rec = {
                "disease":  (row.get("Disease") or "").strip(),
                "status":   (row.get("Status") or "").strip(),
                "plasmy":   (row.get("Plasmy Reports (Homo/Hetero)") or "").strip(),
                "gb_freq":  (row.get("GB Freq FL(CR)") or "").strip(),
                "gb_seqs":  (row.get("GB Seqs FL(CR)") or "").strip(),
                "refs":     (row.get("References") or "").strip(),
                "mitotip":  "",
                "allele":   (row.get("Allele") or "").strip(),
                "aa_change": (row.get("Amino Acid Change") or "").strip(),
            }
            m = re.fullmatch(r"([ACGTN])-([ACGTN])", nuc, re.I)
            if m:
                out[(pos, m.group(1).upper(), m.group(2).upper())] = rec
    return out


def _load_mitomap_rna(path: Path) -> dict[tuple[int, str, str], dict]:
    """tRNA/rRNA table → {(pos, ref, alt): record}. Allele is
    <ref><pos><alt> for SNVs (e.g. A576G); no AA change for these."""
    out: dict[tuple, dict] = {}
    if not path or not path.is_file():
        return out
    with path.open("r", encoding="latin-1", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                pos = int((row.get("Position") or "").strip())
            except ValueError:
                continue
            allele = (row.get("Allele") or "").strip()
            mt = (row.get("MitoTIP") or "").strip()
            rec = {
                "disease":  (row.get("Disease") or "").strip(),
                "status":   (row.get("Status") or "").strip(),
                "plasmy":   f"{(row.get('Homoplasmy') or '').strip()}/{(row.get('Heteroplasmy') or '').strip()}",
                "gb_freq":  (row.get("GB Freq FL(CR)") or "").strip(),
                "gb_seqs":  (row.get("GB Seqs FL(CR)") or "").strip(),
                "refs":     (row.get("References") or "").strip(),
                "mitotip":  "" if mt.upper() in ("", "N/A") else mt,
                "allele":   allele,
                "aa_change": "",
            }
            m = _RNA_ALLELE_RE.fullmatch(allele)
            if m:
                out[(pos, m.group(1).upper(), m.group(3).upper())] = rec
    return out


_EMPTY_MM = {"disease": "", "status": "", "plasmy": "", "gb_freq": "",
             "gb_seqs": "", "refs": "", "mitotip": "", "allele": "", "aa_change": ""}


def _mitomap_lookup(cc: dict, rna: dict, pos: int, ref: str, alt: str, locus_type: str) -> dict:
    """Exact (pos, ref, alt) lookup. tRNA/rRNA prefer the rna table,
    everything else the coding/control table. No POS-only fallback."""
    tables = ([rna, cc] if locus_type in ("tRNA", "rRNA") else [cc, rna])
    for t in tables:
        rec = t.get((pos, ref.upper(), alt.upper()))
        if rec:
            return rec
    return dict(_EMPTY_MM)


def _consequence(locus_type: str, aa_change: str) -> str:
    """Coarse consequence label without VEP — non-coding for the RNA /
    control loci; for protein-coding, infer from the MITOMAP AA change."""
    if locus_type == "tRNA":     return "non-coding (tRNA)"
    if locus_type == "rRNA":     return "non-coding (rRNA)"
    if locus_type == "control":  return "control region"
    if locus_type == "intergenic": return "intergenic"
    # protein-coding:
    aa = (aa_change or "").strip()
    if not aa or aa.lower() in ("noncoding", "-", "na"):
        return "coding (other)"
    m = _AA_RE.fullmatch(aa)
    if m:
        a1, a2 = m.group(1), m.group(3)
        if a2 in _STOP_AA or a1 in _STOP_AA:
            return "stop_gained" if a2 in _STOP_AA else "stop_lost"
        return "synonymous" if a1 == a2 else "missense"
    if any(tok in aa for tok in _STOP_AA):
        return "stop_gained"
    return "coding (other)"


def _open(path: str):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--mitomap_cc", required=True)
    ap.add_argument("--mitomap_rna", required=True)
    ap.add_argument("--sample_id", default="")
    ap.add_argument("--output", required=True)
    args = ap.parse_args(argv)

    cc  = _load_mitomap_cc(Path(args.mitomap_cc))
    rna = _load_mitomap_rna(Path(args.mitomap_rna))

    out_cols = [
        "CHROM", "POS", "REF", "ALT", "HGVS_M", "GENE", "LOCUS_TYPE",
        "CONSEQUENCE", "AA_CHANGE",
        "HETEROPLASMY", "AD", "DEPTH", "FILTER", "TLOD",
        "MITOMAP_DISEASE", "MITOMAP_STATUS", "MITOMAP_PLASMY",
        "MITOMAP_GB_FREQ", "MITOMAP_GB_SEQS", "MITOMAP_REFS",
        "MITOTIP_SCORE", "MITOMAP_ALLELE",
    ]

    sample_idx = None
    rows: dict[tuple, dict] = {}     # (pos,ref,alt) → row (best TLOD wins)

    with _open(args.vcf) as f:
        for line in f:
            if line.startswith("#CHROM"):
                hdr = line.rstrip("\n").split("\t")
                if args.sample_id and args.sample_id in hdr:
                    sample_idx = hdr.index(args.sample_id)
                else:
                    sample_idx = 9 if len(hdr) > 9 else None
                continue
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos_s, _id, ref, alts_s = cols[0], cols[1], cols[2], cols[3], cols[4]
            # Skip nuclear rows when handed a whole-genome VCF (DRAGEN
            # emits all calls in one file). mtDNA coords used by
            # _gene_at / _mitomap_lookup are meaningless off chrM.
            if chrom not in ("chrM", "MT", "chrMT"):
                continue
            filt = cols[6]
            info = cols[7]
            try:
                pos = int(pos_s)
            except ValueError:
                continue

            # INFO → TLOD (per-alt list when multiallelic). Mutect2-mito
            # emits TLOD; DRAGEN's chrM caller doesn't, but ships a
            # roughly-equivalent FORMAT/SQ (somatic quality) per ALT.
            # Read TLOD if present, else fall back to SQ — same column
            # in the output TSV, so the downstream adapter doesn't have
            # to know which caller produced the file.
            info_d = {}
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_d[k] = v
                else:
                    info_d[kv] = True
            tlod_list = []
            tl = info_d.get("TLOD", "")
            if isinstance(tl, str) and tl:
                tlod_list = tl.split(",")

            # FORMAT / proband sample → AF (heteroplasmy, per-alt), AD, DP
            fmt_keys = cols[8].split(":") if len(cols) > 8 else []
            samp_vals = cols[sample_idx].split(":") if (sample_idx is not None and len(cols) > sample_idx) else []
            fd = {k: (samp_vals[i] if i < len(samp_vals) else "") for i, k in enumerate(fmt_keys)}
            af_list = (fd.get("AF") or "").split(",")
            dp = fd.get("DP", "")
            ad_all = (fd.get("AD") or "").split(",")   # [refDepth, alt1Depth, alt2Depth, ...]
            # DRAGEN fallback: FORMAT/SQ when INFO/TLOD is absent.
            if not tlod_list:
                sq_raw = (fd.get("SQ") or "").strip()
                if sq_raw and sq_raw not in (".", "NA"):
                    tlod_list = sq_raw.split(",")

            for ai, alt in enumerate(alts_s.split(",")):
                alt = alt.strip()
                if not alt or alt == "*":
                    continue
                het = af_list[ai].strip() if ai < len(af_list) else ""
                # per-alt AD = refDepth + this alt's depth (best effort)
                if ad_all and len(ad_all) > ai + 1:
                    ad = f"{ad_all[0]},{ad_all[ai + 1]}"
                else:
                    ad = fd.get("AD", "")
                tlod = tlod_list[ai].strip() if ai < len(tlod_list) else (tlod_list[0].strip() if tlod_list else "")
                try:
                    tlod_f = float(tlod) if tlod else float("-inf")
                except ValueError:
                    tlod_f = float("-inf")

                gene, locus_type = _gene_at(pos)
                hgvs_m = _hgvs_m(pos, ref, alt)
                mm = _mitomap_lookup(cc, rna, pos, ref, alt, locus_type)
                aa = mm["aa_change"]
                cons = _consequence(locus_type, aa)

                row = {
                    "CHROM": chrom or "chrM", "POS": pos, "REF": ref, "ALT": alt,
                    "HGVS_M": hgvs_m, "GENE": gene, "LOCUS_TYPE": locus_type,
                    "CONSEQUENCE": cons, "AA_CHANGE": aa,
                    "HETEROPLASMY": het, "AD": ad, "DEPTH": dp,
                    "FILTER": filt, "TLOD": tlod,
                    "MITOMAP_DISEASE": mm["disease"], "MITOMAP_STATUS": mm["status"],
                    "MITOMAP_PLASMY": mm["plasmy"], "MITOMAP_GB_FREQ": mm["gb_freq"],
                    "MITOMAP_GB_SEQS": mm["gb_seqs"], "MITOMAP_REFS": mm["refs"],
                    "MITOTIP_SCORE": mm["mitotip"], "MITOMAP_ALLELE": mm["allele"],
                    "_tlod": tlod_f,
                }
                key = (pos, ref, alt)
                prev = rows.get(key)
                if prev is None or row["_tlod"] > prev["_tlod"]:
                    rows[key] = row

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(out_cols)
        for key in sorted(rows.keys(), key=lambda k: (k[0], k[1], k[2])):
            r = rows[key]
            w.writerow([r.get(c, "") for c in out_cols])

    sys.stderr.write(f"[parse_mito_vcf] {len(rows)} variants → {out_path}\n")


if __name__ == "__main__":
    main(sys.argv[1:])
