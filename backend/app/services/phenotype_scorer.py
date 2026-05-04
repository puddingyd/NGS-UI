"""Python port of the R script's in-house pheno_score (vcf-analysis-hg19.R §1542).

Given the patient's HPO list (with weights) and any selected custom
panels, compute a per-gene score:

    pheno_score(gene) = 100 * Σ matching_weight / total_weight

where each (HPO term ∪ panel name) the patient carries that maps to
`gene` contributes its weight to the numerator. Panels are folded into
the same lookup table by treating the panel name as a synthetic
"hpo_id" (matching the R script's behaviour).

Loaded once at startup:
  - phenotype_data/phenotype_to_genes.txt  (~1M rows, 65 MB)
  - phenotype_data/gene_panels/*.txt       (69 panels, ~27 k gene rows)

A reload happens on demand via reload_db() (no auto-watch).
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from ..config import REPO_ROOT

PHENO_DATA_DIR = REPO_ROOT / "phenotype_data"
PHENO_TO_GENES_PATH = PHENO_DATA_DIR / "phenotype_to_genes.txt"
PANELS_DIR = PHENO_DATA_DIR / "gene_panels"


# hpo_id (or panel_name) → set[gene_symbol]
_HPO_TO_GENES: dict[str, set[str]] = defaultdict(set)
# panel_name → set[gene_symbol]
_PANEL_TO_GENES: dict[str, set[str]] = {}
_LOADED = False


def _load_phenotype_to_genes(path: Path = PHENO_TO_GENES_PATH) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    if not path.exists():
        return out
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            hid = (row.get("hpo_id") or "").strip()
            gene = (row.get("gene_symbol") or "").strip()
            if hid and gene and gene != "-":
                out[hid].add(gene)
    return out


def _load_panels(panel_dir: Path = PANELS_DIR) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    if not panel_dir.exists():
        return out
    for fp in sorted(panel_dir.glob("*.txt")):
        name = fp.stem
        genes: set[str] = set()
        # A few legacy panel files were saved as Latin-1 (Windows export);
        # fall back so a single bad byte doesn't kill the whole loader.
        for enc in ("utf-8", "latin-1"):
            try:
                with fp.open("r", encoding=enc) as f:
                    for line in f:
                        g = line.strip().split("\t")[0].strip()
                        if g and not g.startswith("#"):
                            genes.add(g)
                break
            except UnicodeDecodeError:
                genes.clear()
                continue
        if genes:
            out[name] = genes
    return out


def load() -> tuple[int, int]:
    """Idempotent. Returns (n_hpo_terms_loaded, n_panels_loaded)."""
    global _LOADED
    if _LOADED:
        return len(_HPO_TO_GENES), len(_PANEL_TO_GENES)
    hpo_map = _load_phenotype_to_genes()
    _HPO_TO_GENES.clear()
    _HPO_TO_GENES.update(hpo_map)
    panels = _load_panels()
    _PANEL_TO_GENES.clear()
    _PANEL_TO_GENES.update(panels)
    # Fold panels into the same lookup table — panel name acts as a
    # synthetic hpo_id, mirroring the R script's `rbind(custom_panels_df, hp_db)`.
    for panel_name, genes in panels.items():
        _HPO_TO_GENES[panel_name] |= genes
    _LOADED = True
    return len(_HPO_TO_GENES), len(_PANEL_TO_GENES)


def reload_db() -> tuple[int, int]:
    global _LOADED
    _LOADED = False
    return load()


def list_panels() -> list[dict]:
    if not _LOADED:
        load()
    return [
        {"name": name, "gene_count": len(genes)}
        for name, genes in sorted(_PANEL_TO_GENES.items())
    ]


def gene_count(hpo_id: str) -> int:
    """Number of distinct genes annotated to this HPO term (or panel name)."""
    if not _LOADED:
        load()
    return len(_HPO_TO_GENES.get(hpo_id, ()))


def compute_pheno_score(
    hpo_terms: list[dict] | list[tuple[str, float]],
    panels: Iterable = (),
) -> dict[str, float]:
    """Return {gene_symbol: pheno_score} for genes with score > 0.

    `hpo_terms` accepts either:
        [{"phenotype": "HP:0001250", "weight": 2}, ...]   (dict form)
        [("HP:0001250", 2), ...]                          (tuple form)

    `panels` accepts mixed:
        ["HIE", "Marfan_panel"]                              (legacy strings → weight 1)
        [{"name": "HIE", "weight": 2}, ...]                  (preferred)
        [("HIE", 2), ...]                                    (tuple)
    """
    if not _LOADED:
        load()
    pairs: list[tuple[str, float]] = []
    for entry in hpo_terms or []:
        if isinstance(entry, dict):
            hid = (entry.get("phenotype") or entry.get("hpo_id") or "").strip()
            try:
                w = float(entry.get("weight", 1) or 1)
            except (TypeError, ValueError):
                w = 1.0
        else:
            hid, w = entry[0], float(entry[1])
        if hid:
            pairs.append((hid, w))
    for entry in panels or []:
        if isinstance(entry, dict):
            name = (entry.get("name") or "").strip()
            try:
                w = float(entry.get("weight", 1) or 1)
            except (TypeError, ValueError):
                w = 1.0
        elif isinstance(entry, str):
            name, w = entry.strip(), 1.0
        else:
            name, w = entry[0], float(entry[1])
        if name and name in _PANEL_TO_GENES:
            pairs.append((name, w))

    total_weight = sum(w for _, w in pairs)
    if total_weight <= 0 or not pairs:
        return {}

    accum: dict[str, float] = defaultdict(float)
    for hid, w in pairs:
        for gene in _HPO_TO_GENES.get(hid, ()):
            accum[gene] += w
    return {
        g: 100.0 * v / total_weight
        for g, v in accum.items()
        if v > 0
    }


def write_pheno_table(sample_id: str, pheno_score: dict[str, float]) -> Path:
    """Persist gene → score as `<sample_dir>/pheno_score.tsv` (sorted desc).

    This is the equivalent of the R script's `<ID>.pheno.txt` and lets
    downstream code (and an operator who SSH's in) inspect what got
    scored without re-running the math.
    """
    from ..config import TERTIARY_OUTPUT_ROOT
    sub = TERTIARY_OUTPUT_ROOT / sample_id
    sub.mkdir(parents=True, exist_ok=True)
    out = sub / "pheno_score.tsv"
    rows = sorted(pheno_score.items(), key=lambda kv: -kv[1])
    with out.open("w", encoding="utf-8", newline="") as f:
        f.write("gene_symbol\tpheno_score\n")
        for g, s in rows:
            f.write(f"{g}\t{s:.4f}\n")
    return out


def update_in_panel_column(sample_id: str, in_panel_genes: set[str]) -> int:
    """Rewrite IN_PANEL column of snv_indel.annotated.tsv. Returns rows updated."""
    from ..config import TERTIARY_OUTPUT_ROOT
    tsv = TERTIARY_OUTPUT_ROOT / sample_id / "snv_indel.annotated.tsv"
    if not tsv.exists():
        return 0
    rows: list[list[str]] = []
    with tsv.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        try:
            gene_idx = header.index("GENE")
            inpanel_idx = header.index("IN_PANEL")
        except ValueError:
            return 0
        rows.append(header)
        n_updated = 0
        for r in reader:
            if len(r) <= max(gene_idx, inpanel_idx):
                rows.append(r); continue
            new_val = "true" if r[gene_idx] in in_panel_genes else "false"
            if r[inpanel_idx] != new_val:
                r[inpanel_idx] = new_val
                n_updated += 1
            rows.append(r)
    # Atomic-ish: write to .tmp then rename
    tmp = tsv.with_suffix(".tsv.tmp")
    with tmp.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        for r in rows:
            w.writerow(r)
    tmp.replace(tsv)
    return n_updated
