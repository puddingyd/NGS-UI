"""Read tertiary_output/{SAMPLE_ID}/* and shape it for the frontend.

For Phase 2 we only handle SNV/indel + sample_metadata. CNV / SV / STR /
SF / PGx adapters land in later phases and will be merged into the same
returned dict.

Layout (post-migration):
    tertiary_output/{sid}/
      sample_metadata.json     (patient-level: meta + reviewer fields +
                                active_analysis)
      snv_indel.annotated.tsv  (variant calls; same TSV across versions)
      analyses/{ver}/
        analysis.json          (hpo + selected_panels + note)
        pheno_score.tsv
        exomiser_results.tsv
        lirical_results.tsv

Pre-migration layout (sidecars + hpo at the sample root) is still
recognized as a fallback so the loader keeps working between deploy
and migration.
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from ..adapters.snv_tsv import TIERS, load_snv_tsv
from ..config import INDEX_PATH, TERTIARY_OUTPUT_ROOT
from . import analyses_store


def _read_json_or(path: Path, default):
    if not path.exists():
        return default
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return default


def _read_tsv_dict(path: Path, key_col: str = "VARIANT_ID") -> dict[str, dict]:
    """Read a tiny sidecar TSV into {key: row_dict}. Returns {} if absent."""
    if not path.exists():
        return {}
    import csv as _csv
    out: dict[str, dict] = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = _csv.DictReader(f, delimiter="\t")
        for row in reader:
            k = (row.get(key_col) or "").strip()
            if k:
                out[k] = row
    return out


def _to_num(s):
    if s in (None, ""):
        return None
    try:
        f = float(s)
        return int(f) if f.is_integer() else f
    except (TypeError, ValueError):
        return s


def list_index() -> list[dict]:
    """Return the sample list for the top-bar combobox.

    Prefer `_index.json` if present; otherwise scan subdirectories and
    build a minimal entry from each `sample_metadata.json` (or just the
    directory name).
    """
    # Prefer the canonical NGS_UI/_index.json. Tolerate the legacy path
    # inside tertiary_output/ so deployments mid-migration keep working.
    for idx_path in (INDEX_PATH, TERTIARY_OUTPUT_ROOT / "_index.json"):
        if idx_path.exists():
            data = _read_json_or(idx_path, [])
            return data if isinstance(data, list) else []

    out: list[dict] = []
    if not TERTIARY_OUTPUT_ROOT.exists():
        return out
    for sub in sorted(TERTIARY_OUTPUT_ROOT.iterdir()):
        if not sub.is_dir() or sub.name.startswith("_"):
            continue
        meta = _read_json_or(sub / "sample_metadata.json", {}) or {}
        out.append({
            "sample_id":     meta.get("sample_id") or sub.name,
            "lis_id":        meta.get("lis_id") or sub.name,
            "name":          meta.get("name", ""),
            "mrn":           meta.get("mrn", ""),
            "test_type":     meta.get("test_type", ""),
            "category":      meta.get("category", ""),
            "run_date":      meta.get("run_date", ""),
            "tags":          meta.get("tags", []),
            "has_completed": (sub / "snv_indel.annotated.tsv").exists(),
            "tertiary_dir":  sub.name,
        })
    return out


def _normalize_phenotype(entries: list) -> list[dict]:
    out = []
    for entry in entries or []:
        if isinstance(entry, str):
            out.append({"phenotype": entry, "label": entry})
        elif isinstance(entry, dict):
            out.append({
                "phenotype": entry.get("phenotype") or entry.get("hpo_id") or "",
                "label":     entry.get("label")     or entry.get("hpo_name") or "",
                "weight":    entry.get("weight"),
            })
    return out


def _resolve_version(sample_id: str, requested: str | None,
                     meta_active: str | None) -> str | None:
    """Pick the analysis version to load.

    Priority: explicit query param → sample_metadata.active_analysis →
    'default' if it exists → first available version → None (legacy
    pre-migration sample, sidecars at sample root).
    """
    versions = analyses_store.list_versions(sample_id)
    names = {v["name"] for v in versions}
    if requested and requested in names:
        return requested
    if meta_active and meta_active in names:
        return meta_active
    if "default" in names:
        return "default"
    if versions:
        return versions[0]["name"]
    return None


def load_sample(sample_id: str, version: str | None = None) -> dict | None:
    """Build the per-sample webdata payload the frontend renders.

    `version` selects which analysis sidecar set to join in. When None,
    falls back to the sample's `active_analysis`, then `default`, then
    the legacy flat layout (sidecars at the sample root).
    """
    sub = TERTIARY_OUTPUT_ROOT / sample_id
    if not sub.is_dir():
        return None

    snv_tsv = sub / "snv_indel.annotated.tsv"
    if snv_tsv.exists():
        variants, categories = load_snv_tsv(snv_tsv)
    else:
        variants, categories = {}, {t: [] for t in TIERS}

    meta = _read_json_or(sub / "sample_metadata.json", {}) or {}

    # Decide which directory holds the sidecar TSVs for this load.
    chosen_version = _resolve_version(
        sample_id,
        requested=version,
        meta_active=meta.get("active_analysis"),
    )
    if chosen_version is not None:
        sidecar_dir = analyses_store.version_dir(sample_id, chosen_version)
    else:
        # Pre-migration fallback: sidecars used to live at the sample root.
        sidecar_dir = sub

    # HPO/panels: prefer the chosen analysis version; fall back to legacy
    # fields on sample_metadata.json for un-migrated samples.
    if chosen_version is not None:
        analysis = analyses_store.read_version(sample_id, chosen_version) or {}
        hpo_list      = analysis.get("hpo") or []
        panels_list   = analysis.get("selected_panels") or []
    else:
        hpo_list      = meta.get("hpo") or meta.get("patient_phenotype") or []
        panels_list   = meta.get("selected_panels") or []

    # Join per-variant Exomiser / LIRICAL scores from the sidecar TSVs that
    # the rerun worker writes. Either may be absent (not run yet); cards
    # silently omit those rows when the field is None.
    exo = _read_tsv_dict(sidecar_dir / "exomiser_results.tsv")
    lir = _read_tsv_dict(sidecar_dir / "lirical_results.tsv")
    pheno_by_gene: dict[str, float] = {}
    pheno_path = sidecar_dir / "pheno_score.tsv"
    if pheno_path.exists():
        import csv as _csv
        with pheno_path.open("r", encoding="utf-8", newline="") as f:
            for row in _csv.DictReader(f, delimiter="\t"):
                gene = (row.get("gene_symbol") or "").strip()
                try:
                    pheno_by_gene[gene] = float(row.get("pheno_score") or 0)
                except ValueError:
                    pass
    def _scale_to_100(s):
        n = _to_num(s)
        if not isinstance(n, (int, float)):
            return None
        return int(round(n * 100))

    for vid, v in variants.items():
        gene = v.get("gene_symbol", "")
        if gene and gene in pheno_by_gene:
            v["pheno_score"] = round(pheno_by_gene[gene], 2)
        # Total = variant + pheno; either missing → treat as 0 in the
        # sum but only emit a total when at least one component exists.
        gs = v.get("geno_score")
        ps = v.get("pheno_score")
        if gs is not None or ps is not None:
            v["total_score"] = (gs or 0) + (ps or 0)
        e = exo.get(vid)
        if e:
            # Exomiser writes 0–1 floats; rescale to the 0–100 grid the
            # other scores live on so the card has one consistent unit.
            v["total_score_exomiser_variant"] = _scale_to_100(e.get("EXOMISER_GENE_COMBINED_SCORE"))
            v["pheno_score_exomiser"]         = _scale_to_100(e.get("EXOMISER_GENE_PHENO_SCORE"))
            v["exomiser_variant_score"]       = _scale_to_100(e.get("EXOMISER_VARIANT_SCORE"))
            v["rank_exomiser_variant"]        = _to_num(e.get("EXOMISER_RANK"))
        l = lir.get(vid)
        if l:
            v["lirical_variant_score"] = _to_num(l.get("LIRICAL_VARIANT_SCORE"))
            v["rank_lirical_variant"]  = _to_num(l.get("RANK_LIRICAL_VARIANT"))
            v["lirical_disease_name"]  = l.get("DISEASE_NAME") or ""
            v["lirical_disease_curie"] = l.get("DISEASE_CURIE") or ""

    # Re-sort each tier by total_score desc now that pheno_score is
    # joined. Adapter only had ACMG_POINTS available; here we have the
    # full composite total = variant + pheno. Tie-break by id so the
    # ordering is stable across reloads.
    def _ts(vid: str) -> float:
        ts = variants.get(vid, {}).get("total_score")
        return float(ts) if isinstance(ts, (int, float)) else float("-inf")
    for t, ids in categories.items():
        categories[t] = sorted(ids, key=lambda i: (-_ts(i), i))

    qc  = _read_json_or(sub / "qc_summary.json",  {}) or {}
    roh = _read_json_or(sub / "roh_summary.json", {}) or {}

    return {
        "meta": {
            "LIS_ID":   meta.get("lis_id") or meta.get("sample_id") or sample_id,
            "Name":     meta.get("name", ""),
            "MRN":      meta.get("mrn", ""),
            "Test":     meta.get("test_type", ""),
            "Category": meta.get("category", ""),
        },
        "sample_id":         sample_id,
        "genome_build":      meta.get("genome_build", "hg38"),
        "generated_at":      meta.get("run_date") or datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "patient_phenotype": _normalize_phenotype(hpo_list),
        "selected_panels":   panels_list,
        "vcf_path":          meta.get("vcf_path", ""),
        "qc_summary":        qc,
        "roh_summary":       roh,
        "variants":          variants,
        "categories":        categories,
        "tiers":             TIERS,
        "pharmcat":          {},
        # Active version metadata so the frontend can show a version
        # picker / detect when re-analysis should ask for a target.
        "active_analysis":   chosen_version,
        "analyses":          analyses_store.list_versions(sample_id),
    }
