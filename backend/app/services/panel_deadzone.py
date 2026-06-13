"""Panel scope, VEP→HGNC symbol translation, and dead-zone lookups.

The data comes from the self-contained `ngs_panel_deadzone/` delivery
package. All helpers degrade to empty lookups when the package is absent so
sample loading and report generation remain available in development.
"""
from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

from ..config import NGS_PANEL_DEADZONE_DIR


def _read_tsv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


@lru_cache(maxsize=1)
def hgnc_id_to_symbol() -> dict[str, str]:
    path = NGS_PANEL_DEADZONE_DIR / "panel" / "hgnc_id_to_symbol.tsv"
    out: dict[str, str] = {}
    try:
        for row in _read_tsv(path):
            hid = (row.get("hgnc_id") or "").strip()
            sym = (row.get("symbol") or "").strip()
            if hid and sym:
                out[hid] = sym
    except OSError:
        return {}
    return out


@lru_cache(maxsize=1)
def symbol_to_hgnc_id() -> dict[str, str]:
    return {sym: hid for hid, sym in hgnc_id_to_symbol().items()}


@lru_cache(maxsize=1)
def _current_symbol_by_upper() -> dict[str, str]:
    return {sym.upper(): sym for sym in hgnc_id_to_symbol().values()}


@lru_cache(maxsize=1)
def disease_associated_genes() -> frozenset[str]:
    paths = [
        NGS_PANEL_DEADZONE_DIR / "panel" / "panel_loose_plus_clinical.hgnc_canonical.txt",
        NGS_PANEL_DEADZONE_DIR / "panel" / "panel_loose.hgnc_canonical.txt",
    ]
    for path in paths:
        try:
            with path.open("r", encoding="utf-8") as fh:
                return frozenset(line.strip() for line in fh if line.strip())
        except OSError:
            continue
    return frozenset()


@lru_cache(maxsize=1)
def disease_associated_hgnc_ids() -> frozenset[str]:
    by_symbol = symbol_to_hgnc_id()
    return frozenset(
        by_symbol[sym]
        for sym in disease_associated_genes()
        if sym in by_symbol
    )


@lru_cache(maxsize=1)
def _vep_missing_alias_by_hgnc_id() -> dict[str, tuple[str, str]]:
    """VEP HGNC_ID → (panel HGNC_ID, panel HGNC-current symbol)."""
    path = NGS_PANEL_DEADZONE_DIR / "panel" / "vep_alias_map.tsv"
    out: dict[str, tuple[str, str]] = {}
    try:
        for row in _read_tsv(path):
            if (row.get("kind") or "").strip() != "vep_missing_gene":
                continue
            vep_hid = (row.get("vep_hgnc_id") or "").strip()
            panel_hid = (row.get("panel_hgnc_id") or "").strip()
            panel_sym = (row.get("panel_hgnc_symbol") or "").strip()
            if vep_hid and panel_hid and panel_sym:
                out[vep_hid] = (panel_hid, panel_sym)
    except OSError:
        return {}
    return out


@lru_cache(maxsize=1)
def _legacy_symbol_alias() -> dict[str, str]:
    """VEP emitted symbol → HGNC-current symbol, fallback when HGNC_ID is absent."""
    path = NGS_PANEL_DEADZONE_DIR / "panel" / "vep_alias_map.tsv"
    out: dict[str, str] = {}
    try:
        for row in _read_tsv(path):
            vep_sym = (row.get("vep_symbol") or "").strip()
            panel_sym = (row.get("panel_hgnc_symbol") or "").strip()
            if vep_sym and panel_sym:
                out[vep_sym] = panel_sym
    except OSError:
        return {}
    return out


@lru_cache(maxsize=1)
def _panel_gene_alias() -> dict[str, str]:
    """Panel/custom-panel input alias → HGNC-current symbol.

    This intentionally lives beside, but separate from, the VEP alias map:
    VEP-specific `vep_missing_gene` rows can describe positional annotation
    drift and must not be applied blindly to user-supplied panel gene lists.
    """
    path = NGS_PANEL_DEADZONE_DIR / "panel" / "panel_gene_aliases.tsv"
    out: dict[str, str] = {}
    try:
        for row in _read_tsv(path):
            alias = (row.get("alias") or "").strip()
            symbol = (row.get("hgnc_symbol") or "").strip()
            if alias and symbol:
                out[alias] = symbol
                out.setdefault(alias.upper(), symbol)
    except OSError:
        return {}
    return out


def canonical_gene_symbol(symbol: str, hgnc_id: str = "") -> tuple[str, str]:
    """Return (HGNC-current symbol, HGNC_ID) for a VEP/AnnotSV gene.

    HGNC_ID is preferred because VEP may emit old symbols while keeping the
    stable HGNC ID. For known VEP-missing split genes, use the package's
    cross-reference back to the panel gene.
    """
    sym = (symbol or "").strip()
    hid = (hgnc_id or "").strip()

    if hid:
        alias = _vep_missing_alias_by_hgnc_id().get(hid)
        if alias:
            return alias[1], alias[0]
        mapped = hgnc_id_to_symbol().get(hid)
        if mapped:
            return mapped, hid

    mapped = _legacy_symbol_alias().get(sym)
    if mapped:
        return mapped, symbol_to_hgnc_id().get(mapped, hid)
    return sym, hid


def canonical_panel_gene_symbol(symbol: str) -> str:
    """Return HGNC-current symbol for HPO/panel/custom-panel input.

    Unlike `canonical_gene_symbol()`, this avoids VEP positional alias rows
    that are only safe when a variant's HGNC_ID or transcript context is
    available. It is therefore the right helper for pheno_score gene sets.
    """
    sym = (symbol or "").strip()
    if not sym:
        return ""
    by_upper = _current_symbol_by_upper()
    current = hgnc_id_to_symbol()
    if sym in current.values():
        return sym
    cased = by_upper.get(sym.upper())
    if cased:
        return cased
    mapped = _panel_gene_alias().get(sym) or _panel_gene_alias().get(sym.upper())
    if mapped:
        return by_upper.get(mapped.upper(), mapped)
    return sym


def is_disease_associated_gene(symbol: str, hgnc_id: str = "") -> bool:
    sym, hid = canonical_gene_symbol(symbol, hgnc_id)
    genes = disease_associated_genes()
    ids = disease_associated_hgnc_ids()
    return (bool(sym) and sym in genes) or (bool(hid) and hid in ids)


def _compact_int_runs(values: list[int]) -> str:
    if not values:
        return ""
    xs = sorted(set(v for v in values if isinstance(v, int) and v > 0))
    if not xs:
        return ""
    runs: list[str] = []
    start = end = xs[0]
    for val in xs[1:]:
        if val == end + 1:
            end = val
            continue
        runs.append(str(start) if start == end else f"{start}-{end}")
        start = end = val
    runs.append(str(start) if start == end else f"{start}-{end}")
    return ", ".join(runs)


def _parse_exon_list(raw: str) -> list[int]:
    out: list[int] = []
    for part in (raw or "").split(","):
        part = part.strip()
        if not part:
            continue
        try:
            out.append(int(part))
        except ValueError:
            continue
    return sorted(set(out))


def _dead_zone_file_for(test_type: str) -> tuple[Path, int, str, str, str]:
    if (test_type or "").upper() == "WES":
        return (
            NGS_PANEL_DEADZONE_DIR / "dead_zone" / "wes" / "wes_panel_dead_exon_summary.tsv",
            20,
            "dead_20x_exons",
            "dead_20x_pct",
            "dead_20x_len",
        )
    return (
        NGS_PANEL_DEADZONE_DIR / "dead_zone" / "dragen_wgs" / "wgs_dragen_panel_dead_exon_summary.tsv",
        15,
        "dead_15x_exons",
        "dead_15x_pct",
        "dead_15x_len",
    )


def _parse_float(raw: str) -> float:
    try:
        return float(str(raw or "").strip())
    except ValueError:
        return 0.0


def _parse_int(raw: str) -> int:
    try:
        return int(float(str(raw or "").strip()))
    except ValueError:
        return 0


@lru_cache(maxsize=4)
def _dead_exons_by_gene(test_type_key: str) -> tuple[dict[str, dict], int]:
    path, threshold, exon_col, pct_col, len_col = _dead_zone_file_for(test_type_key)
    out: dict[str, dict] = {}
    try:
        for row in _read_tsv(path):
            gene, _ = canonical_gene_symbol(row.get("gene", ""))
            exons = _parse_exon_list(row.get(exon_col, ""))
            if not gene or not exons:
                continue
            cds_dead_pct = _parse_float(row.get(pct_col, ""))
            cds_dead_len = _parse_int(row.get(len_col, ""))
            cds_len = _parse_int(row.get("cds_len", ""))
            entry = out.setdefault(
                gene,
                {
                    "gene": gene,
                    "threshold": threshold,
                    "exons": set(),
                    "transcripts": set(),
                    "cds_len": 0,
                    "cds_dead_len": 0,
                    "cds_dead_pct": 0.0,
                },
            )
            entry["exons"].update(exons)
            entry["cds_len"] = max(entry["cds_len"], cds_len)
            entry["cds_dead_len"] = max(entry["cds_dead_len"], cds_dead_len)
            entry["cds_dead_pct"] = max(entry["cds_dead_pct"], cds_dead_pct)
            tx = (row.get("transcript") or "").strip()
            if tx:
                entry["transcripts"].add(tx)
    except OSError:
        return {}, threshold

    normalized: dict[str, dict] = {}
    for gene, entry in out.items():
        exons = sorted(entry["exons"])
        normalized[gene] = {
            "gene": gene,
            "threshold": threshold,
            "exons": exons,
            "exons_label": _compact_int_runs(exons),
            "transcripts": sorted(entry["transcripts"]),
            "cds_len": entry["cds_len"],
            "cds_dead_len": entry["cds_dead_len"],
            "cds_dead_pct": round(entry["cds_dead_pct"], 1),
        }
    return normalized, threshold


def dead_zone_threshold(test_type: str) -> int:
    return _dead_zone_file_for(test_type)[1]


def dead_zone_for_genes(test_type: str, genes: set[str] | list[str] | tuple[str, ...]) -> dict[str, dict]:
    key = "WES" if (test_type or "").upper() == "WES" else "WGS"
    all_dead, _ = _dead_exons_by_gene(key)
    out: dict[str, dict] = {}
    for raw in genes:
        gene, _ = canonical_gene_symbol(str(raw or ""))
        if gene in all_dead:
            out[gene] = all_dead[gene]
    return dict(sorted(
        out.items(),
        key=lambda item: (-float(item[1].get("cds_dead_pct") or 0), item[0]),
    ))


def dead_zone_suffix(gene: str, test_type: str) -> str:
    hit = dead_zone_for_genes(test_type, {gene}).get(canonical_gene_symbol(gene)[0])
    if not hit:
        return ""
    return f"exon {hit['exons_label']} <{hit['threshold']}X"
