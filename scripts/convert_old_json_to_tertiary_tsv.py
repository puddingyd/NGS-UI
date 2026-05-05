#!/usr/bin/env python3
"""
Convert old vcf-analysis-hg38-R webdata JSON (per-sample) into the
snv_indel.annotated.tsv schema defined in 三級輸出計畫.md.

Fields without a source in the old JSON are filled with sensible
placeholders (empty / inferred / hardcoded defaults) and clearly
marked in the docstring of fill_missing() below. Replace these with
real pipeline values once the upstream is wired up.

Usage:
    python scripts/convert_old_json_to_tertiary_tsv.py \\
        --in 26WE0064.json \\
        --out tertiary_output/26WE0064/snv_indel.annotated.tsv
"""
import argparse
import csv
import json
import re
from pathlib import Path

COLUMNS = [
    "CHROM", "POS", "REF", "ALT",
    "GENE", "TRANSCRIPT", "TRANSCRIPT_TYPE",
    "HGVS_C", "HGVS_P", "CONSEQUENCE",
    "MANE_ALL",
    "CALLERS",
    "ZYGOSITY", "GT_DV", "GT_HC",
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


HGVS_RE = re.compile(
    r"^(?P<gene>[^:]+):(?P<tx>[^:]+)(?::(?P<c>c\.[^:]+))?(?::(?P<p>p\.[^:]+))?$"
)


def parse_hgvs(hgvs: str):
    """Split 'PAX7:NM_002584.3:c.335C>T:p.Pro112Leu' into pieces."""
    if not hgvs:
        return {"gene": "", "tx": "", "c": "", "p": ""}
    m = HGVS_RE.match(hgvs)
    if not m:
        return {"gene": "", "tx": "", "c": "", "p": ""}
    return {k: (m.group(k) or "") for k in ("gene", "tx", "c", "p")}


def zygosity_to_gt(zyg: str) -> str:
    """het -> 0/1, hom -> 1/1, hemi -> 1, anything else -> ./."""
    if not zyg:
        return "./."
    z = zyg.lower()
    if z == "het":
        return "0/1"
    if z == "hom":
        return "1/1"
    if z == "hemi":
        return "1"
    return "./."


def empty_to_blank(v):
    if v is None:
        return ""
    return v


def fmt_num(v, digits=None):
    if v is None or v == "":
        return ""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if digits is None:
        return str(f)
    return f"{f:.{digits}f}".rstrip("0").rstrip(".")


def acmg_evidence_pipe(criteria: str) -> str:
    """ACMG_criteria 'PP3,BP4_Moderate,BS2' -> 'PP3|BP4_Moderate|BS2'."""
    if not criteria:
        return ""
    parts = [p.strip() for p in re.split(r"[,;|]", str(criteria)) if p.strip()]
    return "|".join(parts)


def derive_consequence_loftee_hc(cons: str) -> str:
    """Stub LOFTEE_HC: HC for high-confidence LoF consequences, else blank.

    NOTE: In the real pipeline this comes from LOFTEE plugin output. Here we
    just guess from the consequence string so the new TSV column has SOMETHING
    plausible — replace with real value when wiring up VEP+LOFTEE.
    """
    if not cons:
        return ""
    lof_terms = {
        "stop_gained", "frameshift_variant", "splice_acceptor_variant",
        "splice_donor_variant", "stop_lost", "start_lost",
        "transcript_ablation",
    }
    return "HC" if cons in lof_terms else ""


def to_omim_link(omim_id) -> str:
    if omim_id in (None, ""):
        return ""
    return f"https://www.omim.org/entry/{omim_id}"


def to_gnomad_link(chrom, pos, ref, alt, build="hg38") -> str:
    if not all([chrom, pos, ref, alt]):
        return ""
    chrom_clean = str(chrom).replace("chr", "")
    dataset = "gnomad_r4" if build == "hg38" else "gnomad_r2_1"
    return (
        f"https://gnomad.broadinstitute.org/variant/"
        f"{chrom_clean}-{pos}-{ref}-{alt}?dataset={dataset}"
    )


def to_clinvar_link(chrom, pos, ref, alt) -> str:
    if not all([chrom, pos, ref, alt]):
        return ""
    chrom_clean = str(chrom).replace("chr", "")
    return (
        f"https://www.ncbi.nlm.nih.gov/clinvar/?term="
        f"{chrom_clean}%5BCHR%5D+AND+{pos}%5BCHRPOS%5D"
    )


def convert_variant(v, pheno_genes_set, build="hg38"):
    """Map one old-style variant dict to the new TSV row dict."""
    parts = parse_hgvs(v.get("HGVS", ""))
    gene = v.get("gene_symbol") or parts["gene"]
    transcript = parts["tx"]
    hgvs_c = parts["c"]
    hgvs_p = parts["p"]
    cons = v.get("Consequence", "") or ""

    gt = zygosity_to_gt(v.get("zygosity", ""))

    return {
        "CHROM": v.get("CHROM", ""),
        "POS": v.get("POS", ""),
        "REF": v.get("REF", ""),
        "ALT": v.get("ALT", ""),
        "GENE": gene,
        "TRANSCRIPT": transcript,
        # Placeholder: old data has no MANE_SELECT flag. Default everything
        # to MANE_SELECT for now since the old pipeline reported one
        # canonical transcript per variant.
        "TRANSCRIPT_TYPE": "MANE_SELECT" if transcript else "",
        "HGVS_C": hgvs_c,
        "HGVS_P": hgvs_p,
        "CONSEQUENCE": cons,
        # Placeholder: empty array, real pipeline emits all MANE transcripts.
        "MANE_ALL": "[]",
        # Placeholder: assume both callers agreed (DV+HC) since old pipeline
        # only reported consensus calls.
        "CALLERS": "DV+HC",
        "ZYGOSITY": v.get("zygosity", ""),
        # Placeholder: derive both from zygosity — real pipeline will give
        # per-caller GT.
        "GT_DV": gt,
        "GT_HC": gt,
        "AD":  empty_to_blank(v.get("AD")),
        "VAF": empty_to_blank(v.get("alt_af")),
        "CLINVAR_SIG": empty_to_blank(v.get("CLNSIG")),
        "CLINVAR_STARS": empty_to_blank(v.get("clinvar_stars")),
        # Placeholder: old data has no CLNDN.
        "CLINVAR_DN": "",
        "CLINVAR_CONF": empty_to_blank(v.get("CLNSIGCONF")),
        "GNOMAD_G_AF": empty_to_blank(v.get("AF")),
        "GNOMAD_G_EAS_AF": empty_to_blank(v.get("AF_eas")),
        # Placeholder: old data only has gnomAD genome AF, not exome split.
        "GNOMAD_E_AF": "",
        "GNOMAD_E_EAS_AF": "",
        "TWB_AF": empty_to_blank(v.get("TaiwanBioBank")),
        # Placeholder: P-KNN not in old data.
        "PKNN_LLR": "",
        # Placeholder: REVEL / BayesDel not in old data.
        "REVEL": "",
        "BAYESDEL": "",
        "ALPHAMISSENSE": empty_to_blank(v.get("AlphaMissense_score")),
        "METARNN": empty_to_blank(v.get("MetaRNN_score")),
        # Placeholder: ESM2 / Evo2 not in old data.
        "ESM2_SCORE": "",
        "EVO2_SCORE": "",
        "SPLICEAI_MAX": empty_to_blank(v.get("SpliceAI_score")),
        # Placeholder: CADD not in old data.
        "CADD_PHRED": "",
        # Placeholder: LOFTEE not in old data — guess HC from consequence.
        "LOFTEE_HC": derive_consequence_loftee_hc(cons),
        "LOFTEE_FILTER": "",
        "LOFTEE_FLAGS": "",
        "ACMG_EVIDENCE": acmg_evidence_pipe(v.get("ACMG_criteria", "")),
        "ACMG_POINTS": empty_to_blank(v.get("ACMG_score")),
        "ACMG_CLASS": empty_to_blank(v.get("ACMG_classification")),
        # Placeholder: WhatsHap phasing not in old data.
        "PHASE_GROUP": "",
        "PHASE_RESULT": "unphased",
        # Placeholder: ROH overlap not computed in old data.
        "IN_ROH": "false",
        # Real: in panel iff gene is in pheno_genes (HPO-derived gene list).
        "IN_PANEL": "true" if gene and gene in pheno_genes_set else "false",
        # Placeholder: blacklist not in old data.
        "IN_BLACKLIST": "false",
        "OMIM_LINK": to_omim_link(v.get("OMIM_id")),
        "GNOMAD_LINK": to_gnomad_link(
            v.get("CHROM"), v.get("POS"), v.get("REF"), v.get("ALT"), build=build
        ),
        "CLINVAR_LINK": to_clinvar_link(
            v.get("CHROM"), v.get("POS"), v.get("REF"), v.get("ALT")
        ),
        # User-filled in UI.
        "REPORT_CLASS": "",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    args = ap.parse_args()

    data = json.loads(Path(args.inp).read_text(encoding="utf-8"))
    pheno_genes = set(data.get("pheno_genes") or [])
    build_in = data.get("genome_build", "hg38")
    build = "hg19" if build_in == "hg19" else "hg38"

    rows = []
    for vid, v in data.get("variants", {}).items():
        rows.append(convert_variant(v, pheno_genes, build=build))

    rows.sort(key=lambda r: (str(r["CHROM"]), int(r["POS"]) if r["POS"] else 0))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"Wrote {len(rows)} rows × {len(COLUMNS)} columns to {out_path}")


if __name__ == "__main__":
    main()
