#!/usr/bin/env python3
"""parse_mito_vep.py — VEP-annotated mtDNA VCF + MITOMAP tables → mito.annotated.tsv

Input:
  --vep_vcf      <sample>.mito.vep.vcf.gz  (from annotate_mito_vcf.sh)
  --mitomap_cc   mitomap_mutations_coding_control.tsv
  --mitomap_rna  mitomap_mutations_rna.tsv
  --sample_id    sample id (the FORMAT-column sample name; auto-detected if omitted)
  --output       mito.annotated.tsv

Output columns:
  CHROM POS REF ALT HGVS_M GENE LOCUS_TYPE CONSEQUENCE HGVS_C HGVS_P AA_CHANGE
  HETEROPLASMY AD DEPTH FILTER TLOD
  MITOMAP_DISEASE MITOMAP_STATUS MITOMAP_PLASMY MITOMAP_GB_FREQ MITOMAP_GB_SEQS
  MITOMAP_REFS MITOTIP_SCORE MITOMAP_ALLELE

Notes:
  * LOCUS_TYPE comes from VEP's BIOTYPE when the variant falls inside a
    gene (Mt_tRNA→tRNA, Mt_rRNA→rRNA, protein_coding→protein); D-loop /
    intergenic variants (only up/downstream consequences) → control,
    GENE=MT-CR.
  * Multiple records at the same (POS,REF,ALT) — Mutect2-mito re-emits
    from different assembly regions, or norm-split STR alts — are
    de-duplicated, keeping the highest-TLOD copy.
  * MITOMAP is joined by (POS, REF, ALT) only — no POS-only fallback,
    so a different allele at the same position never bleeds in.
    MITOMAP indels (with their own notation) just stay unmatched.
  * Heteroplasmy = FORMAT/AF, depth = FORMAT/DP, AD = FORMAT/AD.
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

_BIOTYPE_TO_LOCUS = {
    "Mt_tRNA": "tRNA",
    "Mt_rRNA": "rRNA",
    "protein_coding": "protein",
}

# VEP consequence severity — used to pick the "best" CSQ entry.
_CSQ_SEVERITY = {
    "transcript_ablation": 0, "stop_gained": 1, "frameshift_variant": 2,
    "stop_lost": 3, "start_lost": 4, "transcript_amplification": 5,
    "inframe_insertion": 6, "inframe_deletion": 7, "missense_variant": 8,
    "protein_altering_variant": 9, "splice_region_variant": 10,
    "incomplete_terminal_codon_variant": 11, "start_retained_variant": 12,
    "stop_retained_variant": 13, "synonymous_variant": 14,
    "coding_sequence_variant": 15, "mature_miRNA_variant": 16,
    "non_coding_transcript_exon_variant": 17, "intron_variant": 18,
    "NMD_transcript_variant": 19, "non_coding_transcript_variant": 20,
    "upstream_gene_variant": 21, "downstream_gene_variant": 22,
    "intergenic_variant": 23,
}
_NON_GENIC = {"upstream_gene_variant", "downstream_gene_variant", "intergenic_variant"}

_RNA_ALLELE_RE = re.compile(r"^([ACGTN])(\d+)([ACGTN])$", re.I)


def _gene_at(pos: int) -> tuple[str, str]:
    """(gene, locus_type) for a position, from the rCRS map. D-loop →
    control region; the few tiny intergenic gaps between genes →
    ('', 'intergenic')."""
    for name, start, end, lt in _MT_GENES:
        if start <= pos <= end:
            return name, lt
    return ("MT-CR", "control") if (pos <= 576 or pos >= 16024) else ("", "intergenic")


def _load_mitomap_cc(path: Path) -> dict[tuple[int, str, str], dict]:
    """coding/control table → {(pos, ref, alt): record}, keyed off the
    "Nucleotide Change" (ref-alt) column for SNVs.

    MITOMAP exports are Latin-1, not UTF-8 (stray 0xa0 etc.), so read
    them as latin-1 — it decodes any byte and the data is otherwise
    plain ASCII."""
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
            }
            m = re.fullmatch(r"([ACGTN])-([ACGTN])", nuc, re.I)
            if m:
                out[(pos, m.group(1).upper(), m.group(2).upper())] = rec
    return out


def _load_mitomap_rna(path: Path) -> dict[tuple[int, str, str], dict]:
    """tRNA/rRNA table → {(pos, ref, alt): record}.
    Allele here is `<ref><pos><alt>` for SNVs (e.g. A576G)."""
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
            }
            m = _RNA_ALLELE_RE.fullmatch(allele)
            if m:
                out[(pos, m.group(1).upper(), m.group(3).upper())] = rec
    return out


def _mitomap_lookup(cc: dict, rna: dict, pos: int, ref: str, alt: str, locus_type: str) -> dict:
    """Exact (pos, ref, alt) lookup only. tRNA/rRNA variants prefer the
    rna table, everything else the coding/control table. No POS-only
    fallback — matching a different allele at the same position would
    surface the wrong disease/freq data (e.g. m.114C>A vs m.114C>T)."""
    tables = ([rna, cc] if locus_type in ("tRNA", "rRNA") else [cc, rna])
    for t in tables:
        rec = t.get((pos, ref.upper(), alt.upper()))
        if rec:
            return rec
    return {"disease": "", "status": "", "plasmy": "", "gb_freq": "",
            "gb_seqs": "", "refs": "", "mitotip": "", "allele": ""}


def _open(path: str):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, "r", encoding="utf-8")


def _parse_csq_format(line: str) -> list[str]:
    m = re.search(r'Format:\s*([^"]+)"', line)
    return m.group(1).split("|") if m else []


def _pick_csq(csq_field: str, csq_cols: list[str]) -> dict:
    """Return the 'best' CSQ entry as {col: value}: prefer a genic
    consequence (not up/downstream), most-severe, then PICK=1."""
    best = None
    best_key = None
    for entry in csq_field.split(","):
        vals = entry.split("|")
        d = {c: (vals[i] if i < len(vals) else "") for i, c in enumerate(csq_cols)}
        cons = d.get("Consequence", "")
        first_cons = cons.split("&")[0]
        sev = min((_CSQ_SEVERITY.get(c, 99) for c in cons.split("&")), default=99)
        genic = first_cons not in _NON_GENIC
        pick = 1 if (d.get("PICK") or "").strip() == "1" else 0
        key = (0 if genic else 1, sev, -pick)   # genic first, then severity, then PICK
        if best_key is None or key < best_key:
            best_key, best = key, d
    return best or {}


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--vep_vcf", required=True)
    ap.add_argument("--mitomap_cc", required=True)
    ap.add_argument("--mitomap_rna", required=True)
    ap.add_argument("--sample_id", default="")
    ap.add_argument("--output", required=True)
    args = ap.parse_args(argv)

    cc = _load_mitomap_cc(Path(args.mitomap_cc))
    rna = _load_mitomap_rna(Path(args.mitomap_rna))

    out_cols = [
        "CHROM", "POS", "REF", "ALT", "HGVS_M", "GENE", "LOCUS_TYPE",
        "CONSEQUENCE", "HGVS_C", "HGVS_P", "AA_CHANGE",
        "HETEROPLASMY", "AD", "DEPTH", "FILTER", "TLOD",
        "MITOMAP_DISEASE", "MITOMAP_STATUS", "MITOMAP_PLASMY",
        "MITOMAP_GB_FREQ", "MITOMAP_GB_SEQS", "MITOMAP_REFS",
        "MITOTIP_SCORE", "MITOMAP_ALLELE",
    ]

    csq_cols: list[str] = []
    sample_idx = None
    rows: dict[tuple, dict] = {}   # (pos,ref,alt) → row dict (best TLOD wins)

    with _open(args.vep_vcf) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                csq_cols = _parse_csq_format(line)
                continue
            if line.startswith("#CHROM"):
                hdr = line.rstrip("\n").split("\t")
                if args.sample_id and args.sample_id in hdr:
                    sample_idx = hdr.index(args.sample_id)
                else:
                    sample_idx = 9 if len(hdr) > 9 else None   # first sample col
                continue
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos_s, _id, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            filt = cols[6]
            info = cols[7]
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            # multiallelic should already be split by bcftools norm; if a
            # comma somehow survives, take the first alt.
            alt = alt.split(",")[0]

            # INFO fields we want
            info_d = {}
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_d[k] = v
                else:
                    info_d[kv] = True
            tlod = info_d.get("TLOD", "")
            if isinstance(tlod, str) and "," in tlod:
                tlod = tlod.split(",")[0]
            try:
                tlod_f = float(tlod) if tlod not in ("", True) else float("-inf")
            except ValueError:
                tlod_f = float("-inf")

            # FORMAT / sample → heteroplasmy (AF), AD, DP
            het = ad = dp = ""
            if sample_idx is not None and len(cols) > sample_idx and len(cols) > 8:
                fmt = cols[8].split(":")
                vals = cols[sample_idx].split(":")
                fd = {k: (vals[i] if i < len(vals) else "") for i, k in enumerate(fmt)}
                af = fd.get("AF", "")
                het = af.split(",")[0] if af else ""
                ad = fd.get("AD", "")
                dp = fd.get("DP", "")

            # CSQ → gene / consequence / HGVS
            csq_raw = info_d.get("CSQ", "") if isinstance(info_d.get("CSQ"), str) else ""
            csq = _pick_csq(csq_raw, csq_cols) if csq_raw and csq_cols else {}
            cons = csq.get("Consequence", "")
            first_cons = cons.split("&")[0] if cons else ""
            biotype = csq.get("BIOTYPE", "")
            symbol = csq.get("SYMBOL", "")
            hgvsg = csq.get("HGVSg", "")             # e.g. chrM:m.73A>G
            hgvs_m = hgvsg.split(":", 1)[1] if ":" in hgvsg else hgvsg
            hgvs_c = csq.get("HGVSc", "")
            hgvs_p = csq.get("HGVSp", "")
            aa = csq.get("Amino_acids", "")

            # Gene / locus type: prefer VEP when the consequence is genic;
            # else fall back to the rCRS position map (control region etc.).
            if first_cons and first_cons not in _NON_GENIC and symbol:
                gene = symbol
                locus_type = _BIOTYPE_TO_LOCUS.get(biotype, "")
                if not locus_type:
                    _, locus_type = _gene_at(pos)
            else:
                gene, locus_type = _gene_at(pos)
                if locus_type == "control" and not cons:
                    cons = "non_coding_transcript_variant"

            mm = _mitomap_lookup(cc, rna, pos, ref, alt, locus_type)

            row = {
                "CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt,
                "HGVS_M": hgvs_m, "GENE": gene, "LOCUS_TYPE": locus_type,
                "CONSEQUENCE": cons, "HGVS_C": hgvs_c, "HGVS_P": hgvs_p,
                "AA_CHANGE": aa,
                "HETEROPLASMY": het, "AD": ad, "DEPTH": dp, "FILTER": filt,
                "TLOD": tlod if tlod not in (True,) else "",
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

    # Sort by position, write.
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(out_cols)
        for key in sorted(rows.keys(), key=lambda k: (k[0], k[1], k[2])):
            r = rows[key]
            w.writerow([r.get(c, "") for c in out_cols])

    sys.stderr.write(f"[parse_mito_vep] {len(rows)} variants → {out_path}\n")


if __name__ == "__main__":
    main(sys.argv[1:])
