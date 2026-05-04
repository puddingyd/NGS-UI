"""Read tertiary_output/{SAMPLE_ID}/* and shape it for the frontend.

For Phase 2 we only handle SNV/indel + sample_metadata. CNV / SV / STR /
SF / PGx adapters land in later phases and will be merged into the same
returned dict.
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from ..adapters.snv_tsv import TIERS, load_snv_tsv
from ..config import TERTIARY_OUTPUT_ROOT


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
    idx_path = TERTIARY_OUTPUT_ROOT / "_index.json"
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


def load_sample(sample_id: str) -> dict | None:
    """Build the per-sample webdata payload the frontend renders."""
    sub = TERTIARY_OUTPUT_ROOT / sample_id
    if not sub.is_dir():
        return None

    snv_tsv = sub / "snv_indel.annotated.tsv"
    if snv_tsv.exists():
        variants, categories = load_snv_tsv(snv_tsv)
    else:
        variants, categories = {}, {t: [] for t in TIERS}

    # Join per-variant Exomiser / LIRICAL scores from the sidecar TSVs that
    # the rerun worker writes. Either may be absent (not run yet); cards
    # silently omit those rows when the field is None.
    exo = _read_tsv_dict(sub / "exomiser_results.tsv")
    lir = _read_tsv_dict(sub / "lirical_results.tsv")
    # In-house pheno_score is gene-level (one row per gene_symbol), so
    # build a quick lookup and assign by GENE per variant.
    pheno_by_gene: dict[str, float] = {}
    pheno_path = sub / "pheno_score.tsv"
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

    meta = _read_json_or(sub / "sample_metadata.json", {}) or {}
    qc = _read_json_or(sub / "qc_summary.json", {}) or {}
    roh = _read_json_or(sub / "roh_summary.json", {}) or {}

    phenotype = meta.get("hpo") or meta.get("patient_phenotype") or []
    norm_phenotype = []
    for entry in phenotype:
        if isinstance(entry, str):
            norm_phenotype.append({"phenotype": entry, "label": entry})
        elif isinstance(entry, dict):
            norm_phenotype.append({
                "phenotype": entry.get("phenotype") or entry.get("hpo_id") or "",
                "label":     entry.get("label") or entry.get("hpo_name") or "",
                "weight":    entry.get("weight"),
            })

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
        "patient_phenotype": norm_phenotype,
        "selected_panels":   meta.get("selected_panels", []) or [],
        "vcf_path":          meta.get("vcf_path", ""),
        "qc_summary":        qc,
        "roh_summary":       roh,
        "variants":          variants,
        "categories":        categories,
        "tiers":             TIERS,
        "pharmcat":          {},
    }
