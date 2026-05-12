"""mito.annotated.tsv → variant payload adapter.

Reads the per-sample mitochondrial TSV produced by
scripts/annotate_mito_vcf.sh (VEP + the local MITOMAP tables) and
shapes it for the frontend's Mitochondria card.

Tiers (exclusive — a variant lands in the highest one it qualifies for):
    MITO-1  Pathogenic   — MITOMAP status confirmed-ish ("Cfrm" /
                           "Confirmed" / "[P]" / "[LP]") or a MitoTIP
                           "(likely) pathogenic" call
    MITO-2  Clinical     — GENE is in the patient's pheno_score gene
                           set (HPO/panel match) and not already in MITO-1
    MITO-3  Other        — everything else (polymorphisms, haplogroup
                           variants, the poly-C/A tract artefacts)

Heteroplasmy (FORMAT/AF, 0–1) is the headline metric — every tier
sorts by it descending (then by position) so high-load calls float up.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

MITO_TIERS = ["MITO-1", "MITO-2", "MITO-3"]

_PATHO_STATUS_RE = re.compile(r"\bCfrm\b|\bConfirmed\b|\[L?P\]", re.I)
_PATHO_MITOTIP   = {"pathogenic", "likely pathogenic"}


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


def _is_pathogenic(status: str, mitotip: str) -> bool:
    if status and _PATHO_STATUS_RE.search(status):
        return True
    if mitotip and mitotip.strip().lower() in _PATHO_MITOTIP:
        return True
    return False


def load_mito_tsv(
    tsv_path: Path,
    *,
    pheno_by_gene: dict[str, float] | None = None,
) -> tuple[dict[str, dict], dict[str, list[str]]]:
    """Read mito.annotated.tsv → ({id: variant}, {tier: [ids]}).

    `id` is chrM-{pos}-{ref}-{alt} (distinct from SNV/CNV ids, so the
    one flat state.reports.{status,edits} namespace stays collision-free).
    """
    pheno_by_gene = pheno_by_gene or {}
    variants: dict[str, dict] = {}
    categories: dict[str, list[str]] = {t: [] for t in MITO_TIERS}

    if not tsv_path or not tsv_path.exists():
        return variants, categories

    with tsv_path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pos = _to_int(row.get("POS"))
            ref = (row.get("REF") or "").strip()
            alt = (row.get("ALT") or "").strip()
            if pos is None or not ref or not alt:
                continue
            vid = f"chrM-{pos}-{ref}-{alt}"
            gene = (row.get("GENE") or "").strip()
            locus_type = (row.get("LOCUS_TYPE") or "").strip() or "unknown"
            status = (row.get("MITOMAP_STATUS") or "").strip()
            mitotip = (row.get("MITOTIP_SCORE") or "").strip()
            het = _to_float(row.get("HETEROPLASMY"))
            pheno = pheno_by_gene.get(gene)
            in_panel = bool(pheno and pheno > 0)
            pathogenic = _is_pathogenic(status, mitotip)

            v = {
                "id":            vid,
                "CHROM":         (row.get("CHROM") or "chrM").strip(),
                "POS":           pos,
                "REF":           ref,
                "ALT":           alt,
                "HGVS_M":        (row.get("HGVS_M") or "").strip(),
                "gene_symbol":   gene,
                "locus_type":    locus_type,
                "consequence":   (row.get("CONSEQUENCE") or "").strip(),
                "HGVS_C":        (row.get("HGVS_C") or "").strip(),
                "HGVS_P":        (row.get("HGVS_P") or "").strip(),
                "aa_change":     (row.get("AA_CHANGE") or "").strip(),
                "heteroplasmy":  het,                       # 0-1 fraction
                "AD":            (row.get("AD") or "").strip(),
                "depth":         _to_int(row.get("DEPTH")),
                "filter":        (row.get("FILTER") or "").strip(),
                "TLOD":          _to_float(row.get("TLOD")),
                "mitomap_disease": (row.get("MITOMAP_DISEASE") or "").strip(),
                "mitomap_status":  status,
                "mitomap_plasmy":  (row.get("MITOMAP_PLASMY") or "").strip(),
                "mitomap_gb_freq": (row.get("MITOMAP_GB_FREQ") or "").strip(),
                "mitomap_gb_seqs": (row.get("MITOMAP_GB_SEQS") or "").strip(),
                "mitomap_refs":    (row.get("MITOMAP_REFS") or "").strip(),
                "mitotip_score":   mitotip,
                "mitomap_allele":  (row.get("MITOMAP_ALLELE") or "").strip(),
                "pheno_score":     round(pheno, 2) if pheno is not None else None,
                "in_panel":        in_panel,
                "pathogenic":      pathogenic,
            }
            variants[vid] = v
            tier = "MITO-1" if pathogenic else ("MITO-2" if in_panel else "MITO-3")
            categories[tier].append(vid)

    def _key(vid: str) -> tuple:
        v = variants[vid]
        het = v.get("heteroplasmy")
        return (-(het if het is not None else -1.0), v.get("POS") or 0)
    for t in categories:
        categories[t].sort(key=_key)

    return variants, categories
