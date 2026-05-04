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
        "qc_summary":        qc,
        "roh_summary":       roh,
        "variants":          variants,
        "categories":        categories,
        "tiers":             TIERS,
        "pharmcat":          {},
    }
