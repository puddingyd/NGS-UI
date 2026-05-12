"""mito.annotated.tsv → variant payload adapter.

Reads the per-sample mitochondrial TSV produced by
scripts/annotate_mito_vcf.sh (VEP + the local MITOMAP tables) and
shapes it for the frontend's Mitochondria card.

Only disease-relevant variants are surfaced — the card lists variants
that are either (a) pathogenic per MITOMAP/MitoTIP, or (b) carry some
MITOMAP disease association. Polymorphisms / haplogroup variants with
no MITOMAP record are dropped entirely (the raw mito VCF has ~150
variants per sample, almost all of which are exactly that).

Tiers:
    MITO-1  Pathogenic   — MITOMAP status confirmed-ish ("Cfrm" /
                           "Confirmed" / "[P]" / "[LP]") or a MitoTIP
                           "(likely) pathogenic" call
    MITO-2  Clinical     — has a non-empty MITOMAP_DISEASE (and isn't
                           already in MITO-1)

Both tiers sort by a "disease-relevance" key, most-relevant first:
    (status_rank, mitotip_rank, in_panel_rank, -refs, -heteroplasmy, pos)
where status_rank   = Cfrm/Confirmed 0 · Reported 1 · Conflicting 2 · else 3
      mitotip_rank  = Pathogenic 0 · Likely pathogenic 1 · Possibly… 2 · else 3
      in_panel_rank = 0 if GENE is in the patient's pheno_score gene set else 1
heteroplasmy still tie-breaks (a *disease-associated* variant at high load
is more likely clinically significant) but is no longer the headline sort.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

MITO_TIERS = ["MITO-1", "MITO-2"]

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


def _status_rank(status: str) -> int:
    s = (status or "").lower()
    if "cfrm" in s or "confirmed" in s:
        return 0
    if "reported" in s:
        return 1
    if "conflicting" in s or "disputed" in s:
        return 2
    return 3


def _mitotip_rank(mitotip: str) -> int:
    m = (mitotip or "").strip().lower()
    if m == "pathogenic":
        return 0
    if m == "likely pathogenic":
        return 1
    if m.startswith("possibly"):
        return 2
    return 3


def _refs_count(refs: str) -> int:
    """MITOMAP "References" is usually an integer count; be lenient."""
    n = _to_int(refs)
    if n is not None:
        return n
    # fall back to counting non-empty tokens
    return len([t for t in re.split(r"[;, ]+", refs or "") if t])


def load_mito_tsv(
    tsv_path: Path,
    *,
    pheno_by_gene: dict[str, float] | None = None,
) -> tuple[dict[str, dict], dict[str, list[str]]]:
    """Read mito.annotated.tsv → ({id: variant}, {tier: [ids]}).

    `id` is chrM-{pos}-{ref}-{alt} (distinct from SNV/CNV ids, so the
    one flat state.reports.{status,edits} namespace stays collision-free).
    Only disease-relevant variants make it into the returned dicts.
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
            gene = (row.get("GENE") or "").strip()
            locus_type = (row.get("LOCUS_TYPE") or "").strip() or "unknown"
            status = (row.get("MITOMAP_STATUS") or "").strip()
            mitotip = (row.get("MITOTIP_SCORE") or "").strip()
            disease = (row.get("MITOMAP_DISEASE") or "").strip()
            pathogenic = _is_pathogenic(status, mitotip)
            # Only keep disease-relevant variants: pathogenic, or with
            # any MITOMAP disease association. Everything else (the bulk
            # — polymorphisms / haplogroup variants) is dropped.
            if not pathogenic and not disease:
                continue

            vid = f"chrM-{pos}-{ref}-{alt}"
            het = _to_float(row.get("HETEROPLASMY"))
            pheno = pheno_by_gene.get(gene)
            in_panel = bool(pheno and pheno > 0)
            refs_raw = (row.get("MITOMAP_REFS") or "").strip()

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
                "mitomap_disease": disease,
                "mitomap_status":  status,
                "mitomap_plasmy":  (row.get("MITOMAP_PLASMY") or "").strip(),
                "mitomap_gb_freq": (row.get("MITOMAP_GB_FREQ") or "").strip(),
                "mitomap_gb_seqs": (row.get("MITOMAP_GB_SEQS") or "").strip(),
                "mitomap_refs":    refs_raw,
                "mitotip_score":   mitotip,
                "mitomap_allele":  (row.get("MITOMAP_ALLELE") or "").strip(),
                "pheno_score":     round(pheno, 2) if pheno is not None else None,
                "in_panel":        in_panel,
                "pathogenic":      pathogenic,
            }
            variants[vid] = v
            categories["MITO-1" if pathogenic else "MITO-2"].append(vid)

    def _relevance_key(vid: str) -> tuple:
        v = variants[vid]
        het = v.get("heteroplasmy")
        return (
            _status_rank(v.get("mitomap_status", "")),
            _mitotip_rank(v.get("mitotip_score", "")),
            0 if v.get("in_panel") else 1,
            -_refs_count(v.get("mitomap_refs", "")),
            -(het if het is not None else -1.0),
            v.get("POS") or 0,
        )
    for t in categories:
        categories[t].sort(key=_relevance_key)

    return variants, categories
