"""AnnotSV TSV → variant payload adapter.

AnnotSV emits one TSV per SV-caller stream (CNV vs SV). Each TSV has
multiple rows per SV: one `Annotation_mode=full` row carrying the
SV-level scoring + pathogenicity overlap, plus N
`Annotation_mode=split` rows (one per gene the SV intersects)
carrying per-gene transcript / OMIM / Location info.

This adapter aggregates by AnnotSV_ID into one variant record per SV,
classifies into the 4 reviewer-facing tiers, and joins per-gene
pheno_score for the Clinical-block sort.

Tier rules (independent slices — a variant can live in both):
    *-Clinical    : at least one overlapped gene has score > 0 in
                    the sample's pheno_score table (panel/HPO match)
    *-Pathogenic  : ACMG_class ∈ {4, 5} (AnnotSV's own scoring)

Pheno score on each split-row's gene comes from the caller-provided
pheno_by_gene dict (keys are gene_symbols, values 0–100). The variant's
`max_pheno_score` is the highest among matched genes; the gene that
provided it is `trigger_gene` and gets ⭐ in the UI.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

# Tier names mirror the frontend's CNV_SV_TIER_ORDER.
CNV_TIERS = ["CNV-1A", "CNV-1B"]
SV_TIERS  = ["SV-2A",  "SV-2B"]


def _to_float(s: str | None) -> float | None:
    if s is None:
        return None
    s = s.strip()
    if not s or s.upper() in ("NA", "."):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _to_int(s: str | None) -> int | None:
    f = _to_float(s)
    if f is None:
        return None
    try:
        return int(f)
    except (TypeError, ValueError):
        return None


def _split_semi(s: str | None) -> list[str]:
    """AnnotSV concatenates parallel sources/coords/AFs with ';'."""
    if not s:
        return []
    return [p for p in s.split(";") if p]


def _parse_b_block(row: dict, prefix: str) -> dict:
    """B_loss / B_gain / B_ins / B_inv: parallel ;-separated source +
    coord + AFmax lists, one entry per overlapping benign-region
    record. The scalar `af_max` is the max of the per-source AFs
    (used by the frontend for the AF summary line)."""
    sources = _split_semi(row.get(f"{prefix}_source") or "")
    coords  = _split_semi(row.get(f"{prefix}_coord") or "")
    afs_raw = _split_semi(row.get(f"{prefix}_AFmax") or "")
    af_nums = [_to_float(s) for s in afs_raw]
    af_nums = [f for f in af_nums if f is not None]
    return {
        "sources": sources,
        "coords":  coords,
        "afs":     afs_raw,
        "af_max":  max(af_nums) if af_nums else None,
    }


def _parse_format_sample(fmt: str, sample: str) -> dict[str, str]:
    """VCF FORMAT/sample → dict (e.g. GT:CN, 0/1:1 → {GT:0/1, CN:1})."""
    if not fmt or not sample:
        return {}
    keys = fmt.split(":")
    vals = sample.split(":")
    return {k: v for k, v in zip(keys, vals)}


def _zygosity_from_gt(gt: str) -> str:
    """0/1, 1/0, 0|1 → 'het'; 1/1 → 'hom'; 1 (haploid X) → 'hemi'."""
    if not gt:
        return ""
    g = gt.replace("|", "/").strip()
    if g in ("0/1", "1/0"):
        return "het"
    if g == "1/1":
        return "hom"
    if g == "1":
        return "hemi"
    return ""


def _parse_acmg_class(s: str) -> int | None:
    """`ACMG_class` is mostly '1'..'5' on full rows but split rows
    sometimes carry strings like 'full=4'. Return the leading int."""
    if not s:
        return None
    s = s.strip()
    if "=" in s:
        s = s.split("=", 1)[1].strip()
    try:
        return int(s)
    except ValueError:
        return None


def _split_genes(gene_name: str) -> list[str]:
    """Gene_name on full rows is a `;`-separated list of overlapped genes."""
    return [g.strip() for g in (gene_name or "").split(";") if g.strip()]


def _split_row_to_gene(
    row: dict,
    pheno_by_gene: dict[str, float],
    pheno_matched: dict[str, float],
    pheno_total: float,
) -> dict:
    """Build a per-gene record from one AnnotSV split row.

    `pheno_score` (0-100, used for sorting) and `pheno_matched` /
    `pheno_total` (raw matched-weight / total-input-weight pair, used
    by the UI to render the fraction) all describe the same per-gene
    value at different stages of normalisation.
    """
    gene = (row.get("Gene_name") or "").strip()
    score = pheno_by_gene.get(gene)
    matched = pheno_matched.get(gene, 0.0) if pheno_matched else 0.0
    return {
        "gene":             gene,
        "tx":               (row.get("Tx") or "").strip(),
        "tx_version":       (row.get("Tx_version") or "").strip(),
        "location":         (row.get("Location") or "").strip(),
        "location2":        (row.get("Location2") or "").strip(),
        "exon_count":       _to_int(row.get("Exon_count")),
        "frameshift":       (row.get("Frameshift") or "").strip(),
        "overlap_cds_pct":  _to_float(row.get("Overlapped_CDS_percent")),
        "overlap_cds_len":  _to_int(row.get("Overlapped_CDS_length")),
        "omim_id":          (row.get("OMIM_ID") or "").strip(),
        "omim_phenotype":   (row.get("OMIM_phenotype") or "").strip(),
        "omim_inheritance": (row.get("OMIM_inheritance") or "").strip(),
        "loeuf_bin":        _to_float(row.get("LOEUF_bin")),
        "pli":              _to_float(row.get("GnomAD_pLI") or row.get("ExAC_pLI")),
        "pheno_score":      round(score, 2) if score is not None else None,
        "pheno_matched":    matched,
        "pheno_total":      pheno_total,
        "in_panel":         bool(matched and matched > 0),
    }


def _full_row_to_variant(
    full_row: dict,
    *,
    sample_col_idx: int,
    fieldnames: list[str],
    raw_full_values: list[str],
    source: str,
) -> dict:
    """Build the per-SV record from the `full` row."""
    sv_chrom  = (full_row.get("SV_chrom") or "").strip()
    sv_start  = _to_int(full_row.get("SV_start"))
    sv_end    = _to_int(full_row.get("SV_end"))
    sv_type   = (full_row.get("SV_type") or "").strip()
    annotsv_id = (full_row.get("AnnotSV_ID") or "").strip()

    fmt = (full_row.get("FORMAT") or "").strip()
    proband = ""
    if 0 <= sample_col_idx < len(raw_full_values):
        proband = raw_full_values[sample_col_idx]
    proband_dict = _parse_format_sample(fmt, proband)
    gt = proband_dict.get("GT", "")
    cn = _to_int(proband_dict.get("CN"))

    return {
        "id":                annotsv_id,
        "source":            source,                          # "cnv" | "sv"
        "CHROM":             sv_chrom,
        "POS":               sv_start,
        "END":               sv_end,
        "length":            _to_int(full_row.get("SV_length")),
        "sv_type":           sv_type,
        "cytoband":          (full_row.get("CytoBand") or "").strip(),
        "gene_count":        _to_int(full_row.get("Gene_count")),
        "gene_symbol":       (_split_genes(full_row.get("Gene_name") or "") or [""])[0],
        "gene_list":         _split_genes(full_row.get("Gene_name") or ""),
        "genes":             [],   # filled from split rows
        "acmg_class":        _parse_acmg_class(full_row.get("ACMG_class") or ""),
        "ranking_score":     _to_float(full_row.get("AnnotSV_ranking_score")),
        "ranking_criteria":  (full_row.get("AnnotSV_ranking_criteria") or "").strip(),
        "p_loss": {
            "phens":    (full_row.get("P_loss_phen") or "").strip(),
            "sources":  _split_semi(full_row.get("P_loss_source") or ""),
            "coords":   _split_semi(full_row.get("P_loss_coord") or ""),
        },
        "p_gain": {
            "phens":    (full_row.get("P_gain_phen") or "").strip(),
            "sources":  _split_semi(full_row.get("P_gain_source") or ""),
            "coords":   _split_semi(full_row.get("P_gain_coord") or ""),
        },
        "p_ins": {
            "phens":    (full_row.get("P_ins_phen") or "").strip(),
            "sources":  _split_semi(full_row.get("P_ins_source") or ""),
            "coords":   _split_semi(full_row.get("P_ins_coord") or ""),
        },
        "b_loss":            _parse_b_block(full_row, "B_loss"),
        "b_gain":            _parse_b_block(full_row, "B_gain"),
        "b_ins":             _parse_b_block(full_row, "B_ins"),
        "b_inv":             _parse_b_block(full_row, "B_inv"),
        "qual":              _to_float(full_row.get("QUAL")),
        "filter":            (full_row.get("FILTER") or "").strip(),
        "format":            fmt,
        "sample_call":       proband,
        "GT":                gt,
        "zygosity":          _zygosity_from_gt(gt),
        "copy_number":       cn,
        # Filled in classify():
        "in_panel":          False,
        "trigger_gene":      "",
        "max_pheno_score":   None,
    }


def _classify(v: dict, source: str) -> list[str]:
    """Return the tier list a variant belongs to (empty = none)."""
    out: list[str] = []
    if source == "cnv":
        if v.get("in_panel"):
            out.append("CNV-1A")
        ac = v.get("acmg_class")
        if isinstance(ac, int) and ac >= 4:
            out.append("CNV-1B")
    else:  # sv
        if v.get("in_panel"):
            out.append("SV-2A")
        ac = v.get("acmg_class")
        if isinstance(ac, int) and ac >= 4:
            out.append("SV-2B")
    return out


def load_annotsv_tsv(
    tsv_path: Path,
    *,
    source: str,                              # "cnv" | "sv"
    pheno_by_gene: dict[str, float] | None = None,
    pheno_matched: dict[str, float] | None = None,
    pheno_total: float = 0.0,
) -> tuple[dict[str, dict], dict[str, list[str]]]:
    """Read AnnotSV output → ({annotsv_id: variant}, {tier: [ids]}).

    Tier ordering inside each list: Clinical sorts by max pheno_score
    desc (then ranking_score, then id); Pathogenic sorts by
    ranking_score desc (then id). Tie-break is the AnnotSV_ID for
    stable ordering across reloads.
    """
    pheno_by_gene = pheno_by_gene or {}
    pheno_matched = pheno_matched or {}
    tiers = list(CNV_TIERS) if source == "cnv" else list(SV_TIERS)
    variants: dict[str, dict] = {}
    categories: dict[str, list[str]] = {t: [] for t in tiers}

    if not tsv_path.exists():
        return variants, categories

    # We need both fieldnames (for sample-column index) and the raw
    # row values (so we can fish out the proband sample column even
    # when its header is the LIS_ID, not a fixed name).
    with tsv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames: list[str] = list(reader.fieldnames or [])
        # First column after FORMAT is the proband sample. SV TSVs may
        # carry multiple sample columns; we take the first.
        try:
            fmt_idx = fieldnames.index("FORMAT")
            sample_idx = fmt_idx + 1
        except ValueError:
            sample_idx = -1

        # AnnotSV_ID + per-row mode lets us group full+split. We also
        # need to remember the RAW row tuple (positional) for the
        # sample column whose header is dynamic.
        # csv.DictReader doesn't expose ordered values directly, so we
        # keep the raw line via a parallel reader.
        f.seek(0)
        next(f)  # discard header line we already consumed
        raw_reader = csv.reader(f, delimiter="\t")
        for raw_values in raw_reader:
            if not raw_values:
                continue
            row = {k: (raw_values[i] if i < len(raw_values) else "")
                   for i, k in enumerate(fieldnames)}
            mode = (row.get("Annotation_mode") or "").strip()
            aid  = (row.get("AnnotSV_ID") or "").strip()
            if not aid:
                continue
            if mode == "full":
                variants[aid] = _full_row_to_variant(
                    row,
                    sample_col_idx=sample_idx,
                    fieldnames=fieldnames,
                    raw_full_values=raw_values,
                    source=source,
                )
            elif mode == "split":
                # Split rows ride alongside the full row; if we haven't
                # seen the full one yet (file not strictly ordered),
                # stash the gene record and merge later.
                gene_rec = _split_row_to_gene(row, pheno_by_gene, pheno_matched, pheno_total)
                if aid in variants:
                    variants[aid]["genes"].append(gene_rec)
                else:
                    variants.setdefault(aid, {"genes": []})["genes"].append(gene_rec)
            # Other modes (none expected) are ignored.

    # Drop any orphan AnnotSV_IDs that only had split rows (no full)
    # — without the full row we have no SV-level info to render, so
    # they'd be meaningless cards. AnnotSV always emits the full row
    # for SVs that pass its filters; missing one means the file was
    # truncated or pre-filtered.
    for aid in list(variants.keys()):
        v = variants[aid]
        if "id" not in v:
            del variants[aid]

    # Compute Clinical-trigger fields + classify.
    for aid, v in variants.items():
        max_score = None
        trigger = ""
        for g in v["genes"]:
            s = g.get("pheno_score")
            if isinstance(s, (int, float)) and s > 0:
                if max_score is None or s > max_score:
                    max_score = s
                    trigger = g.get("gene", "")
        v["in_panel"] = bool(trigger)
        v["trigger_gene"] = trigger
        v["max_pheno_score"] = max_score
        for tier in _classify(v, source):
            categories[tier].append(aid)

    # Sort each tier.
    def _clinical_key(aid: str) -> tuple:
        v = variants[aid]
        score = v.get("max_pheno_score") or -1
        rank  = v.get("ranking_score") or -999
        return (-float(score), -float(rank), aid)

    def _pathogenic_key(aid: str) -> tuple:
        v = variants[aid]
        rank = v.get("ranking_score") or -999
        return (-float(rank), aid)

    for tier in categories:
        if tier.endswith("-1A") or tier.endswith("-2A"):
            categories[tier].sort(key=_clinical_key)
        else:
            categories[tier].sort(key=_pathogenic_key)

    return variants, categories
