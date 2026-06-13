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
  - repo phenotype_data/gene_panels/*.txt  (fixed panels, versioned in git)
  - repo phenotype_data/custom_panels/*.txt  (custom panels, versioned in git)

A reload happens on demand via reload_db() (no auto-watch).
"""
from __future__ import annotations

import csv
import re
import threading
from functools import lru_cache
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from ..config import CUSTOM_GENE_PANELS_DIR, GENE_PANELS_DIR, PHENO_DATA_DIR
from . import panel_deadzone

PHENO_TO_GENES_PATH = PHENO_DATA_DIR / "phenotype_to_genes.txt"
PANELS_DIR = GENE_PANELS_DIR
CUSTOM_PANELS_DIR = CUSTOM_GENE_PANELS_DIR


# hpo_id (or panel_name) → set[gene_symbol]
_HPO_TO_GENES: dict[str, set[str]] = defaultdict(set)
# panel_name → set[gene_symbol]
_PANEL_TO_GENES: dict[str, set[str]] = {}
# panel_name → small metadata parsed from comment headers in *.txt panel files
_PANEL_META: dict[str, dict[str, str]] = {}
# canonical_gene → {"hpo": set[hpo_id], "panel": set[panel_name]}
_GENE_TO_KEYS: dict[str, dict[str, set[str]]] = {}
_LOADED = False
_LOAD_LOCK = threading.Lock()


@lru_cache(maxsize=200_000)
def _canonical_gene(gene: str) -> str:
    return panel_deadzone.canonical_panel_gene_symbol(gene) or (gene or "").strip()


def _load_phenotype_to_genes(path: Path = PHENO_TO_GENES_PATH) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    if not path.exists():
        return out
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            hid = (row.get("hpo_id") or "").strip()
            gene = _canonical_gene(row.get("gene_symbol") or "")
            if hid and gene and gene != "-":
                out[hid].add(gene)
    return out


@lru_cache(maxsize=4096)
def _load_hpo_genes_for_key(hpo_id: str, path_str: str = str(PHENO_TO_GENES_PATH)) -> tuple[str, ...]:
    """Fast path for the phenotype tool's one-HPO gene-list drawer.

    During startup the full phenotype score cache may still be warming.
    A single HPO drawer lookup should not wait for the whole 1M-row map,
    so scan only the requested HPO term and cache that small result.
    """
    key = (hpo_id or "").strip()
    path = Path(path_str)
    if not key or not path.exists():
        return ()
    genes: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if (row.get("hpo_id") or "").strip() != key:
                continue
            gene = _canonical_gene(row.get("gene_symbol") or "")
            if gene and gene != "-":
                genes.add(gene)
    return tuple(sorted(genes))


def _load_panels_from_dir(panel_dir: Path) -> tuple[dict[str, set[str]], dict[str, dict[str, str]]]:
    out: dict[str, set[str]] = {}
    meta: dict[str, dict[str, str]] = {}
    if not panel_dir.exists():
        return out, meta
    for fp in sorted(panel_dir.glob("*.txt")):
        name = fp.stem
        genes: set[str] = set()
        panel_meta: dict[str, str] = {}
        # A few legacy panel files were saved as Latin-1 (Windows export);
        # fall back so a single bad byte doesn't kill the whole loader.
        for enc in ("utf-8", "latin-1"):
            try:
                with fp.open("r", encoding=enc) as f:
                    for line in f:
                        raw = line.strip()
                        if raw.startswith("#"):
                            key, sep, value = raw[1:].partition(":")
                            if sep:
                                panel_meta[key.strip().lower()] = value.strip()
                            continue
                        g = _canonical_gene(raw.split("\t")[0])
                        if g and not g.startswith("#"):
                            genes.add(g)
                break
            except UnicodeDecodeError:
                genes.clear()
                panel_meta.clear()
                continue
        if genes:
            out[name] = genes
            if panel_meta:
                meta[name] = panel_meta
    return out, meta


def _load_panels(panel_dirs: Iterable[Path] = (PANELS_DIR, CUSTOM_PANELS_DIR)) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    meta: dict[str, dict[str, str]] = {}
    for panel_dir in panel_dirs:
        panels, panel_meta = _load_panels_from_dir(panel_dir)
        for name, genes in panels.items():
            out.setdefault(name, set()).update(genes)
        meta.update(panel_meta)
    _PANEL_META.clear()
    _PANEL_META.update(meta)
    return out


def load() -> tuple[int, int]:
    """Idempotent. Returns (n_hpo_terms_loaded, n_panels_loaded)."""
    global _LOADED
    if _LOADED:
        return len(_HPO_TO_GENES), len(_PANEL_TO_GENES)
    with _LOAD_LOCK:
        if _LOADED:
            return len(_HPO_TO_GENES), len(_PANEL_TO_GENES)
        hpo_map = _load_phenotype_to_genes()
        panels = _load_panels()
        # Fold panels into the same lookup table — panel name acts as a
        # synthetic hpo_id, mirroring the R script's `rbind(custom_panels_df, hp_db)`.
        for panel_name, genes in panels.items():
            hpo_map.setdefault(panel_name, set()).update(genes)
        _HPO_TO_GENES.clear()
        _HPO_TO_GENES.update(hpo_map)
        _PANEL_TO_GENES.clear()
        _PANEL_TO_GENES.update(panels)
        gene_to_keys: dict[str, dict[str, set[str]]] = defaultdict(lambda: {"hpo": set(), "panel": set()})
        for key, genes in hpo_map.items():
            bucket = "panel" if key in panels else "hpo"
            for gene in genes:
                gene_to_keys[gene][bucket].add(key)
        _GENE_TO_KEYS.clear()
        _GENE_TO_KEYS.update(gene_to_keys)
        _LOADED = True
        return len(_HPO_TO_GENES), len(_PANEL_TO_GENES)


def reload_db() -> tuple[int, int]:
    global _LOADED
    _LOADED = False
    _canonical_gene.cache_clear()
    _load_hpo_genes_for_key.cache_clear()
    return load()


_PANEL_NAME_RE = re.compile(r"[^A-Za-z0-9_-]+")


def sanitize_panel_name(name: str) -> str:
    """Map an arbitrary user-typed panel name to a filename-safe id:
    non [A-Za-z0-9_-] runs (incl. spaces, dots) → '_', collapse
    repeats, trim leading/trailing '_'. Empty after cleanup → ''."""
    cleaned = _PANEL_NAME_RE.sub("_", (name or "").strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    return cleaned[:64]


def register_custom_panel(name: str, genes: Iterable[str], source: str = "") -> dict:
    """Create a reusable gene panel from a user-supplied gene list.

    The name is sanitised (see sanitize_panel_name); collisions with an
    existing panel are refused. Genes are canonicalised to HGNC-current
    symbols where a safe panel alias exists, then de-duplicated in
    first-seen order. Writes {name}.txt into the runtime custom panels
    directory with a `#source:` header and updates the in-memory tables
    so the panel is usable immediately — no reload / restart needed.

    Returns {"name": <sanitised>, "n_genes": int}.
    Raises ValueError on an empty/invalid name, empty gene list, or a
    name collision.
    """
    if not _LOADED:
        load()
    clean = sanitize_panel_name(name)
    if not clean:
        raise ValueError("panel 名稱清理後為空，請改用含英數的名稱")
    if clean in _PANEL_TO_GENES:
        raise ValueError(f"已存在名為 {clean} 的 panel，請改名")

    seen: set[str] = set()
    ordered: list[str] = []
    for g in genes or []:
        g = _canonical_gene(g or "")
        if g and g not in seen:
            seen.add(g)
            ordered.append(g)
    if not ordered:
        raise ValueError("沒有可用的基因（清空後為空）")

    CUSTOM_PANELS_DIR.mkdir(parents=True, exist_ok=True)
    out = CUSTOM_PANELS_DIR / f"{clean}.txt"
    # Defence in depth: clean is [A-Za-z0-9_-]{1,64} so this can't
    # escape, but resolve-and-check anyway.
    if out.resolve().parent != CUSTOM_PANELS_DIR.resolve():
        raise ValueError("panel 檔名不合法")
    safe_source = " ".join(str(source or "").strip().split())
    header = f"#source: {safe_source}\n" if safe_source else "#source:\n"
    out.write_text(header + "\n".join(ordered) + "\n", encoding="utf-8")

    gene_set = set(ordered)
    _PANEL_TO_GENES[clean] = gene_set
    _HPO_TO_GENES[clean] |= gene_set
    _PANEL_META[clean] = {"source": safe_source} if safe_source else {"source": ""}
    return {"name": clean, "n_genes": len(gene_set), "source": safe_source}


def list_panels() -> list[dict]:
    if not _LOADED:
        load()
    return [
        {
            "name": name,
            "gene_count": len(genes),
            "source": (_PANEL_META.get(name) or {}).get("source", ""),
        }
        for name, genes in sorted(_PANEL_TO_GENES.items())
    ]


def genes_for_key(key: str, kind: str = "") -> dict:
    """Return canonical gene symbols for an HPO term or panel key."""
    clean_key = (key or "").strip()
    clean_kind = (kind or "").strip().lower()
    if not clean_key:
        return {"kind": clean_kind, "key": "", "gene_count": 0, "genes": []}

    if clean_kind == "panel" or (not clean_kind and clean_key in _PANEL_TO_GENES):
        if not _LOADED:
            load()
        genes = sorted(_PANEL_TO_GENES.get(clean_key, ()))
        meta = _PANEL_META.get(clean_key) or {}
        return {
            "kind": "panel",
            "key": clean_key,
            "gene_count": len(genes),
            "genes": genes,
            "source": meta.get("source", ""),
        }

    if _LOADED:
        genes = sorted(_HPO_TO_GENES.get(clean_key, ()))
    else:
        genes = list(_load_hpo_genes_for_key(clean_key))
    return {
        "kind": "hpo",
        "key": clean_key,
        "gene_count": len(genes),
        "genes": genes,
        "source": "phenotype_to_genes.txt",
    }


def memberships_for_gene(gene: str, *, limit_hpo: int = 200, limit_panels: int = 200) -> dict:
    """Return HPO terms and panels whose canonical gene set contains gene."""
    if not _LOADED:
        load()
    query = _canonical_gene(gene or "")
    if not query:
        return {"query": gene or "", "canonical_gene": "", "hpo": [], "panels": []}
    hits = _GENE_TO_KEYS.get(query) or {"hpo": set(), "panel": set()}
    hpo_ids = sorted(hits.get("hpo", ()))
    panel_names = sorted(hits.get("panel", ()))
    return {
        "query": gene,
        "canonical_gene": query,
        "hpo": [{"id": hid} for hid in hpo_ids[:limit_hpo]],
        "hpo_total": len(hpo_ids),
        "panels": [
            {
                "name": name,
                "gene_count": len(_PANEL_TO_GENES.get(name, ())),
                "source": (_PANEL_META.get(name) or {}).get("source", ""),
            }
            for name in panel_names[:limit_panels]
        ],
        "panel_total": len(panel_names),
    }


def gene_count(hpo_id: str) -> int:
    """Number of distinct genes annotated to this HPO term (or panel name)."""
    if not _LOADED:
        load()
    return len(_HPO_TO_GENES.get(hpo_id, ()))


def compute_pheno_match(
    hpo_terms: list[dict] | list[tuple[str, float]],
    panels: Iterable = (),
) -> tuple[dict[str, float], float]:
    """Pre-multiplication state of compute_pheno_score.

    Returns ({gene_symbol: matched_weight}, total_input_weight) where
    matched_weight is the sum of weights of HPO terms / panels that
    contain the gene, and total_input_weight is the sum of all input
    weights. CNV/SV cards render `matched/total` as a fraction so the
    reviewer can see "this gene was implicated by 2 of the 3 panels".

    Input shapes: same as compute_pheno_score.
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
    accum: dict[str, float] = defaultdict(float)
    for hid, w in pairs:
        for gene in _HPO_TO_GENES.get(hid, ()):
            accum[gene] += w
    return dict(accum), total_weight


def compute_pheno_score(
    hpo_terms: list[dict] | list[tuple[str, float]],
    panels: Iterable = (),
) -> dict[str, float]:
    """Return {gene_symbol: pheno_score} for genes with score > 0.

    Score = 100 × matched_weight / total_input_weight. Identical to
    compute_pheno_match() followed by the per-gene normalisation —
    kept as a thin wrapper so all existing callers (SNV pheno join,
    pheno_score.tsv writer, response stats) keep working unchanged.
    """
    matched, total = compute_pheno_match(hpo_terms, panels)
    if total <= 0 or not matched:
        return {}
    return {
        g: 100.0 * w / total
        for g, w in matched.items()
        if w > 0
    }


def write_pheno_table(
    sample_id: str,
    pheno_score: dict[str, float],
    *,
    target_dir: Path | None = None,
) -> Path:
    """Persist gene → score as `pheno_score.tsv` (sorted desc).

    `target_dir` lets callers (e.g. analyses_store.write_version) write
    into a specific version's directory regardless of which version is
    currently active. When omitted, the file lands in the sample's
    active analysis dir — same behaviour as before.
    """
    if target_dir is None:
        from . import analyses_store
        target_dir = analyses_store.active_version_dir(sample_id)
    target_dir.mkdir(parents=True, exist_ok=True)
    out = target_dir / "pheno_score.tsv"
    rows = sorted(pheno_score.items(), key=lambda kv: -kv[1])
    with out.open("w", encoding="utf-8", newline="") as f:
        f.write("gene_symbol\tpheno_score\n")
        for g, s in rows:
            f.write(f"{g}\t{s:.4f}\n")
    return out
