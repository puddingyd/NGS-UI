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
from ..config import TERTIARY_OUTPUT_ROOT
from . import analyses_store, omim_store


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

    Source of truth is a directory scan: every subdirectory of
    tertiary_output/ that has sample_metadata.json shows up here,
    sorted by created_at descending so newly-registered samples land
    at the top. _index.json gets rewritten as a side-effect on every
    scan so the file on disk stays fresh for human inspection /
    pipeline-side enrichment, but the UI never trusts the cache for
    correctness.
    """
    out: list[dict] = []
    if not TERTIARY_OUTPUT_ROOT.exists():
        return out
    for sub in sorted(TERTIARY_OUTPUT_ROOT.iterdir()):
        if not sub.is_dir() or sub.name.startswith("_"):
            continue
        meta_path = sub / "sample_metadata.json"
        if not meta_path.exists():
            # Pipeline-dropped sample that hasn't been registered yet —
            # surfaces through /samples/unregistered, not the search bar.
            continue
        meta = _read_json_or(meta_path, {}) or {}
        out.append({
            "sample_id":     meta.get("sample_id") or sub.name,
            "lis_id":        meta.get("lis_id") or sub.name,
            "name":          meta.get("name", ""),
            "mrn":           meta.get("mrn", ""),
            "sex":           meta.get("sex", ""),
            "test_type":     meta.get("test_type", ""),
            "category":      meta.get("category", ""),
            "run_date":      meta.get("run_date", ""),
            "created_at":    meta.get("created_at", ""),
            "tags":          meta.get("tags", []),
            "has_completed": (sub / "snv_indel.annotated.tsv").exists(),
            "tertiary_dir":  sub.name,
        })
    # Sort newest-first by registration date; samples without a stored
    # created_at fall to the bottom (stable order by lis_id thereafter).
    out.sort(key=lambda r: (r.get("created_at") or "", r.get("lis_id") or ""), reverse=True)

    # Best-effort cache write so an operator browsing the filesystem
    # sees an up-to-date listing. Failures are non-fatal — the read
    # path doesn't depend on this file existing.
    try:
        cache_path = TERTIARY_OUTPUT_ROOT / "_index.json"
        cache_path.write_text(
            json.dumps(out, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
    except OSError:
        pass

    return out


def list_unregistered() -> list[dict]:
    """Sample dirs the tertiary pipeline left behind without a metadata file.

    A directory under tertiary_output/ counts as "unregistered" when it
    has snv_indel.annotated.tsv but NO sample_metadata.json yet. The UI
    surfaces these in the 載入新個案 dropdown so the reviewer can attach
    basic info + HPO without having to retype the LIS_ID.

    For each entry we also look up the matching
        NGS_UI/patient_phenotype/{lis_id}_{mrn}_phenotype.txt
    and surface the parsed contents inline so the modal can show a
    preview chip row + auto-fill the MRN field on selection.

    Sorted by directory mtime descending (newest first).
    """
    from ..config import PHENOTYPE_DIR
    from . import phenotype_io
    out: list[dict] = []
    if not TERTIARY_OUTPUT_ROOT.exists():
        return out
    for sub in TERTIARY_OUTPUT_ROOT.iterdir():
        if not sub.is_dir() or sub.name.startswith("_"):
            continue
        tsv = sub / "snv_indel.annotated.tsv"
        meta = sub / "sample_metadata.json"
        if not tsv.exists() or meta.exists():
            continue
        lis_id = sub.name

        # Resolve a matching phenotype file in the central phenotype dir.
        # Filename convention: {lis_id}_{mrn}_phenotype.txt → strip both
        # ends to recover the MRN. If multiple match (different MRNs)
        # take the lexicographically first; reviewer can re-aim later.
        pheno_payload = None
        if PHENOTYPE_DIR.is_dir():
            matches = sorted(PHENOTYPE_DIR.glob(f"{lis_id}_*_phenotype.txt"))
            if matches:
                pf = matches[0]
                stem = pf.stem  # "lis_mrn_phenotype"
                # stem looks like "{lis_id}_{mrn}_phenotype"
                if stem.startswith(lis_id + "_") and stem.endswith("_phenotype"):
                    mrn = stem[len(lis_id) + 1 : -len("_phenotype")]
                else:
                    mrn = ""
                try:
                    hpo, panels = phenotype_io.parse(pf.read_text(encoding="utf-8"))
                except OSError:
                    hpo, panels = [], []
                pheno_payload = {
                    "path":   str(pf),
                    "mrn":    mrn,
                    "hpo":    hpo,
                    "panels": panels,
                }

        try:
            mtime = sub.stat().st_mtime
        except OSError:
            mtime = 0.0
        out.append({
            "lis_id":     lis_id,
            "tsv_size":   tsv.stat().st_size if tsv.exists() else 0,
            "mtime":      mtime,
            "phenotype":  pheno_payload,
        })
    out.sort(key=lambda r: r["mtime"], reverse=True)
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

    # Lazy backfill: legacy samples + any pheno_score.tsv predating its
    # analysis.json (e.g. HPO/panels touched by a tool that bypassed
    # write_version) get recomputed inline so the Clinical/in-panel
    # consumers downstream always see a fresh table.
    analysis_path = sidecar_dir / "analysis.json"
    needs_backfill = (
        analysis_path.is_file() and (
            not pheno_path.exists()
            or pheno_path.stat().st_mtime < analysis_path.stat().st_mtime
        )
    )
    if needs_backfill:
        try:
            from . import phenotype_scorer
            scores = phenotype_scorer.compute_pheno_score(hpo_list, panels_list)
            if scores:
                phenotype_scorer.write_pheno_table(
                    sample_id, scores, target_dir=sidecar_dir
                )
        except Exception:
            # Backfill is best-effort; the loader still degrades to no
            # pheno column rather than 5xx the whole sample load.
            pass

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

        # OMIM annotation (Disease1..5 + OMIM_id + OMIM_disease +
        # Inheritance). Frontend already renders these when present;
        # missing-file / missing-gene rows just stay empty.
        omim_id = omim_store.parse_omim_id_from_link(v.get("OMIM_link", ""))
        rec = omim_store.lookup(omim_id=omim_id, gene=gene)
        if rec:
            v["OMIM_id"]      = rec.get("OMIM_id", "")
            v["OMIM_disease"] = rec.get("OMIM_disease", "")
            v["Inheritance"]  = rec.get("Inheritance", "")
            for f in ("Disease1", "Disease2", "Disease3", "Disease4", "Disease5"):
                v[f] = rec.get(f, "")

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
            "Sex":      meta.get("sex", ""),
            "DOB":      meta.get("date_of_birth", ""),
            "Test":     meta.get("test_type", ""),
            "Category": meta.get("category", ""),
        },
        "genetic_counseling": meta.get("genetic_counseling", ""),
        "emr_synced_at":      meta.get("emr_synced_at", ""),
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
