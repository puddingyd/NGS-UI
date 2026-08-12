"""Read tertiary_output/{SAMPLE_ID}/* and shape it for the frontend.

For Phase 2 we only handle SNV/indel + sample_metadata. CNV / SV / STR /
SF / PGx adapters land in later phases and will be merged into the same
returned dict.

Layout (post-migration):
    tertiary_output/{sid}/08_postprocessing/
      sample_metadata.json     (patient-level: meta + reviewer fields +
                                active_analysis)
      snv_annotations.sqlite   (sparse overlay for 03_acmg raw TSV)
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
import sqlite3
import threading
import time
from collections import OrderedDict
from datetime import datetime
from pathlib import Path

from ..adapters.snv_tsv import (
    EXCLUDED_GENES,
    TIERS,
    OldFormatError,
    _DP_HARD_FLOOR_WES,
    _DP_LOW_FLAG_WES,
    _DP_LOW_FLAG_WGS,
    _is_mito_chrom,
    _row_to_variant,
    load_snv_tsv,
    merge_snv_variant_row,
)
from ..config import SNV_CACHE_MAX, SNV_CACHE_MAX_RAW_MB
from . import (
    annotation_versions,
    analyses_store,
    clinvar_latest_store,
    gene_disease_store,
    hpo_ontology,
    litvar2_on_demand,
    manual_acmg,
    omim_store,
    panel_deadzone,
    phenotype_scorer,
    ploidy,
    sample_layout,
    snv_gene_index,
    snv_overlay,
    snv_review,
    test_types,
)
from .snv_rows import is_reportable_raw_row


SECONDARY_SNV_PANELS = {
    "acmg_sf": "ACMG_SF_v3.3",
    "stroke": "WGS__神經科__Stroke",
    "carrier": "carrier_mackenzie_1300+",
}

_snv_cache: OrderedDict[tuple, tuple[dict, dict, dict, str]] = OrderedDict()
_snv_cache_lock = threading.Lock()
_CASE_SUMMARY_CACHE_MAX = 128
_case_summary_cache: OrderedDict[tuple, dict[str, str]] = OrderedDict()
_case_summary_cache_lock = threading.Lock()
CASE_SUMMARY_CACHE_NAME = "case_summary.json"
CASE_TABLE_CACHE_NAME = "_case_table.json"
CASE_TABLE_VERSION = 2
_case_table_lock = threading.Lock()


def _effective_test_type(meta: dict, sample_id: str, *, default: str = "") -> str:
    identity = str(meta.get("lis_id") or meta.get("sample_id") or sample_id)
    return test_types.normalize_test_type(
        meta.get("test_type") or "",
        sample_id=identity,
        default=default,
    )


def _fmt_size(path: Path) -> str:
    try:
        size = path.stat().st_size
    except OSError:
        return "missing"
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.1f}{unit}" if unit != "B" else f"{size}B"
        size /= 1024
    return f"{size:.1f}GB"


def _file_size(path: Path) -> int:
    try:
        return path.stat().st_size
    except OSError:
        return 0


def _snv_cache_allowed(tsv_path: Path) -> bool:
    """Keep compact SNV payloads hot, but do not retain giant raw fallbacks."""
    if SNV_CACHE_MAX <= 0:
        return False
    if not tsv_path.name.endswith(snv_review.REVIEW_TSV_NAME):
        limit_bytes = SNV_CACHE_MAX_RAW_MB * 1024 * 1024
        if limit_bytes <= 0 or _file_size(tsv_path) > limit_bytes:
            return False
    return True


def _log_perf(event: str, started: float, **fields) -> None:
    parts = [f"[perf] {event}", f"elapsed={time.perf_counter() - started:.3f}s"]
    parts.extend(f"{key}={value}" for key, value in fields.items())
    print(" ".join(parts), flush=True)


def clear_snv_cache() -> None:
    """Drop parsed SNV payloads after a sample directory is removed."""
    with _snv_cache_lock:
        _snv_cache.clear()
    with _case_summary_cache_lock:
        _case_summary_cache.clear()


def invalidate_sample_cache(sample_dir: Path) -> None:
    """Drop cached payloads for one sample without cold-starting every case."""
    # A unified sample has sources in 03/04/06/07 and state in 08, so a
    # parent-directory comparison no longer identifies every cached entry.
    clear_snv_cache()


def _sample_id_from_state_dir(sample_dir: Path) -> str:
    sample_dir = Path(sample_dir)
    if sample_dir.name == sample_layout.POSTPROCESSING_DIRNAME:
        return sample_dir.parent.name
    return sample_dir.name


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


def _sidecar_file(
    sample_id: str,
    sidecar_dir: Path,
    name: str,
    *,
    for_write: bool = False,
) -> Path:
    return sample_layout.scoped_file(
        sidecar_dir,
        sample_id,
        name,
        for_write=for_write,
    )


def _to_num(s):
    if s in (None, ""):
        return None
    try:
        f = float(s)
        return int(f) if f.is_integer() else f
    except (TypeError, ValueError):
        return s


def _file_signature(path: Path) -> tuple[str, int, int]:
    """Cheap cache signature; absent files are represented consistently."""
    try:
        st = path.stat()
        return (str(path), st.st_mtime_ns, st.st_size)
    except OSError:
        return (str(path), 0, 0)


def _case_summary_signature(sample_dir: Path, omim_sig: tuple | None = None) -> list:
    sample_id = _sample_id_from_state_dir(sample_dir)
    return [
        list(_file_signature(sample_layout.state_file(sample_id, "sample_metadata.json"))),
        list(_file_signature(sample_layout.snv_raw_tsv(sample_id))),
        list(_file_signature(sample_layout.snv_overlay_path(sample_id))),
        list(_file_signature(sample_layout.clinvar_comparison_path(sample_id))),
        list(_file_signature(sample_layout.snv_gene_index_path(sample_id))),
        list(_file_signature(sample_layout.review_tsv(sample_id))),
        list(_file_signature(sample_layout.mito_tsv(sample_id))),
        list(_file_signature(sample_layout.cnv_tsv(sample_id))),
        list(_file_signature(sample_layout.sv_tsv(sample_id))),
        list(omim_sig if omim_sig is not None else omim_store.cache_signature()),
        list(gene_disease_store.cache_signature()),
    ]


def _case_summary_cache_path(sample_dir: Path) -> Path:
    sample_id = _sample_id_from_state_dir(sample_dir)
    return sample_layout.scoped_file(
        sample_dir,
        sample_id,
        CASE_SUMMARY_CACHE_NAME,
        for_write=True,
    )


def _read_case_summary_disk(sample_dir: Path, signature: list) -> dict[str, str] | None:
    data = _read_json_or(_case_summary_cache_path(sample_dir), {}) or {}
    if data.get("signature") != signature:
        return None
    summary = data.get("summary")
    return summary if isinstance(summary, dict) else None


def _write_case_summary_disk(sample_dir: Path, signature: list, summary: dict[str, str]) -> None:
    payload = {
        "signature": signature,
        "summary": summary,
        "updated_at": datetime.utcnow().isoformat(timespec="seconds") + "Z",
    }
    try:
        _case_summary_cache_path(sample_dir).write_text(
            json.dumps(payload, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
    except OSError:
        pass


def _case_table_path() -> Path:
    return sample_layout.global_cache_path(CASE_TABLE_CACHE_NAME)


def _case_table_row_id(row: dict) -> str:
    return str(row.get("tertiary_dir") or row.get("sample_id") or row.get("lis_id") or "").strip()


def _read_case_table_payload() -> dict:
    payload = _read_json_or(_case_table_path(), {}) or {}
    if not isinstance(payload, dict) or payload.get("version") != CASE_TABLE_VERSION:
        return {}
    rows = payload.get("rows")
    if not isinstance(rows, list):
        return {}
    return payload


def _write_case_table_rows(rows: list[dict]) -> None:
    payload = {
        "version": CASE_TABLE_VERSION,
        "updated_at": datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "rows": rows,
    }
    path = _case_table_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    tmp.replace(path)


def _registered_sample_ids() -> set[str]:
    return {
        sample_id for sample_id in sample_layout.iter_sample_ids()
        if sample_layout.state_file(sample_id, "sample_metadata.json").is_file()
    }


def _case_table_is_complete(rows: list[dict]) -> bool:
    ids = _registered_sample_ids()
    row_ids = {_case_table_row_id(row) for row in rows}
    row_ids.discard("")
    return row_ids == ids


def _analysis_phenotype_summary(sample_id: str, meta: dict) -> str:
    version = meta.get("active_analysis") or analyses_store.active_version(sample_id) or "default"
    analysis = analyses_store.read_version(sample_id, version) or {}
    hpo_rows = analysis.get("hpo") or meta.get("hpo") or meta.get("patient_phenotype") or []
    panel_rows = analysis.get("selected_panels") or meta.get("selected_panels") or []
    lines: list[str] = []
    for row in hpo_rows:
        if isinstance(row, dict):
            hpo_id = str(row.get("phenotype") or row.get("hpo_id") or "").strip()
            label = str(row.get("label") or row.get("name") or "").strip()
        else:
            hpo_id = str(row or "").strip()
            label = ""
        if not hpo_id:
            continue
        if not label:
            term = hpo_ontology.get(hpo_id)
            label = term.name if term else ""
        lines.append(f"{label} {hpo_id}".strip())
    for row in panel_rows:
        if isinstance(row, dict):
            name = str(row.get("name") or row.get("panel") or "").strip()
        else:
            name = str(row or "").strip()
        if name:
            lines.append(name)
    return "\n".join(dict.fromkeys(lines))


def _case_table_row_from_sample_dir(
    sample_dir: Path,
    *,
    omim_sig: tuple | None = None,
) -> dict | None:
    dir_sample_id = _sample_id_from_state_dir(sample_dir)
    meta_path = sample_layout.state_file(dir_sample_id, "sample_metadata.json")
    if not meta_path.exists():
        return None
    meta = _read_json_or(meta_path, {}) or {}
    sample_id = str(meta.get("sample_id") or dir_sample_id)
    summary = _case_management_summary(sample_dir, meta, omim_sig=omim_sig)
    return {
        "sample_id":     sample_id,
        "lis_id":        meta.get("lis_id") or dir_sample_id,
        "name":          meta.get("name", ""),
        "mrn":           meta.get("mrn", ""),
        "sex":           meta.get("sex", ""),
        "test_type":     _effective_test_type(meta, dir_sample_id),
        "category":      meta.get("category", ""),
        "run_date":      meta.get("run_date", ""),
        "created_at":    meta.get("created_at", ""),
        "tags":          meta.get("tags", []),
        "has_completed": sample_layout.snv_raw_tsv(dir_sample_id).is_file(),
        "tertiary_dir":  dir_sample_id,
        "phenotype_summary": _analysis_phenotype_summary(sample_id, meta),
        **summary,
    }


def rebuild_case_table() -> list[dict]:
    """Rebuild the denormalized case-list table from all registered samples."""
    rows: list[dict] = []
    omim_sig = omim_store.cache_signature()
    for sample_id in sorted(sample_layout.iter_sample_ids()):
        row = _case_table_row_from_sample_dir(
            sample_layout.state_dir(sample_id), omim_sig=omim_sig
        )
        if row is not None:
            rows.append(row)
    rows.sort(key=lambda r: (r.get("created_at") or "", r.get("lis_id") or ""), reverse=True)
    with _case_table_lock:
        _write_case_table_rows(rows)
    return rows


def update_case_table_row(sample_id: str) -> None:
    """Best-effort refresh of one row in the case-list table."""
    sample_dir = sample_layout.state_dir(sample_id)
    try:
        row = _case_table_row_from_sample_dir(sample_dir)
        if row is None:
            remove_case_table_row(sample_id)
            return
        with _case_table_lock:
            payload = _read_case_table_payload()
            rows = payload.get("rows") if payload else []
            rows = [r for r in rows if _case_table_row_id(r) != sample_id]
            rows.append(row)
            rows.sort(key=lambda r: (r.get("created_at") or "", r.get("lis_id") or ""), reverse=True)
            _write_case_table_rows(rows)
    except Exception as e:
        print(f"[case-table] refresh failed for {sample_id}: {e}", flush=True)


def remove_case_table_row(sample_id: str) -> None:
    """Best-effort removal of one row after a sample is deleted."""
    try:
        with _case_table_lock:
            payload = _read_case_table_payload()
            rows = payload.get("rows") if payload else []
            rows = [r for r in rows if _case_table_row_id(r) != sample_id]
            _write_case_table_rows(rows)
    except Exception as e:
        print(f"[case-table] remove failed for {sample_id}: {e}", flush=True)


def _read_pheno_scores(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    if not path.exists():
        return out
    import csv as _csv
    with path.open("r", encoding="utf-8", newline="") as f:
        for row in _csv.DictReader(f, delimiter="\t"):
            gene = panel_deadzone.canonical_panel_gene_symbol(row.get("gene_symbol") or "")
            if not gene:
                continue
            try:
                score = float(row.get("pheno_score") or 0)
                out[gene] = max(out.get(gene, 0.0), score)
            except ValueError:
                pass
    return out


_ACMG_SHORT = {
    "pathogenic": "P",
    "likely pathogenic": "LP",
    "pathogenic/likely pathogenic": "P/LP",
    "likely pathogenic/pathogenic": "P/LP",
    "uncertain significance": "VUS",
    "vus": "VUS",
    "likely benign": "LB",
    "benign": "B",
}


def _case_variant_label(variant: dict, edits: dict) -> str:
    """Compact SNV label for the case-management table."""
    manual_snapshot = edits.get("manual_acmg")
    acmg = (
        (
            manual_snapshot.get("classification")
            if isinstance(manual_snapshot, dict) else ""
        )
        or edits.get("ACMG_classification")
        or variant.get("effective_acmg_class")
        or variant.get("genebe_acmg_class")
        or variant.get("ACMG_classification")
        or ""
    ).strip()
    acmg_short = _ACMG_SHORT.get(acmg.lower().replace("_", " "), acmg)
    raw_zyg = (variant.get("zygosity") or "").strip()
    zyg_key = raw_zyg.lower().replace("_", " ")
    if zyg_key in ("heterozygous", "heterozygote", "het"):
        zyg = "het"
    elif zyg_key in ("homozygous", "homozygote", "hom"):
        zyg = "hom"
    elif zyg_key in ("hemizygous", "hemizygote", "hemi"):
        zyg = "hemi"
    else:
        zyg = raw_zyg
    hgvs = _selected_snv_hgvs(variant, edits)
    return ", ".join(x for x in (hgvs, acmg_short, zyg) if x)


def _selected_snv_hgvs(variant: dict, edits: dict) -> str:
    options = variant.get("transcript_options") or []
    selected = str(edits.get("selected_transcript_key") or "").strip()
    if selected:
        for opt in options:
            if str(opt.get("key") or "") == selected:
                return str(opt.get("HGVS") or "")
    default_key = str(variant.get("default_transcript_key") or "").strip()
    if default_key:
        for opt in options:
            if str(opt.get("key") or "") == default_key:
                return str(opt.get("HGVS") or "")
    return str(variant.get("HGVS") or "")


def _report_gene_list(hpo_rows: list, panel_entries: list) -> dict:
    phenotype_scorer.load()
    sections: list[dict] = []

    def genes_for(key: str, kind: str = "") -> list[str]:
        genes: list[str] = []
        entry = phenotype_scorer.genes_for_key(key, kind=kind)
        for raw in entry.get("genes") or []:
            gene = panel_deadzone.canonical_panel_gene_symbol(raw)
            if panel_deadzone.is_disease_associated_gene(gene):
                genes.append(gene)
        return sorted(set(genes))

    for row in hpo_rows or []:
        hid = row.get("phenotype") or row.get("hpo_id") or ""
        if not hid:
            continue
        label = row.get("label") or hid
        sections.append({"name": label, "genes": genes_for(hid, kind="hpo")})
    for entry in panel_entries or []:
        name = entry.get("name") if isinstance(entry, dict) else str(entry)
        if not name:
            continue
        sections.append({"name": name, "genes": genes_for(name, kind="panel")})

    merged: set[str] = set()
    for section in sections:
        merged.update(section.get("genes") or [])
    return {"grouped": sections, "merged": sorted(merged)}


def load_report_gene_list(sample_id: str, version: str | None = None) -> dict | None:
    ctx = _load_pheno_context(sample_id, version)
    if ctx is None:
        return None
    _sub, _sidecar_dir, hpo_list, panels_list, _pheno_by_gene = ctx
    return _report_gene_list(_normalize_phenotype(hpo_list), panels_list)


def _case_selected_diseases(variant: dict, edits: dict) -> list[str]:
    """Return only OMIM diseases explicitly selected by the reviewer."""
    picked = edits.get("report_diseases") or {}
    if not isinstance(picked, dict) or not any(bool(v) for v in picked.values()):
        return []
    omim_id = omim_store.parse_omim_id_from_link(variant.get("OMIM_link", ""))
    rec = omim_store.lookup_cached(
        omim_id=omim_id,
        gene=variant.get("gene_symbol", ""),
    ) or {}
    out = []
    for idx, field in enumerate(omim_store.DISEASE_FIELDS, start=1):
        if picked.get(str(idx)) or picked.get(idx):
            disease = omim_store.compact_disease_label(
                rec.get(field) or ""
            )
            if disease and disease not in out:
                out.append(disease)
    return out


def _case_mito_label(variant: dict, edits: dict) -> str:
    """Compact mitochondrial label for the case-management table."""
    acmg = (
        edits.get("ACMG_classification_mito")
        or variant.get("CLNSIG")
        or ""
    ).strip()
    acmg_short = _ACMG_SHORT.get(acmg.lower().replace("_", " "), acmg)
    gene = str(variant.get("gene_symbol") or "").strip()
    hgvs = str(variant.get("HGVS_M") or "").strip()
    variant_label = " ".join(x for x in (gene, hgvs) if x).strip()
    het = variant.get("heteroplasmy")
    het_label = ""
    if isinstance(het, (int, float)):
        het_label = f"{het * 100:.1f}% heteroplasmy"
    return ", ".join(x for x in (variant_label, acmg_short, het_label) if x)


def _case_selected_mito_diseases(variant: dict, edits: dict) -> list[str]:
    """Return ClinVar diseases explicitly selected on the mito card."""
    picked = edits.get("report_diseases_clinvar") or {}
    if not isinstance(picked, dict) or not any(bool(v) for v in picked.values()):
        return []
    diseases = variant.get("clinvar_diseases") or []
    out: list[str] = []
    for idx, disease in enumerate(diseases):
        key = str(idx)
        if picked.get(key) or picked.get(idx):
            label = str(disease or "").strip()
            if label and label not in out:
                out.append(label)
    return out


def _case_cnv_sv_label(variant: dict, edits: dict) -> str:
    """Compact CNV/SV label for the case-management table."""
    raw_acmg = edits.get("ACMG_class_sv")
    try:
        acmg_num = int(raw_acmg) if raw_acmg not in (None, "") else variant.get("acmg_class")
    except (TypeError, ValueError):
        acmg_num = None
    acmg = {
        1: "B", 2: "LB", 3: "VUS", 4: "LP", 5: "P",
    }.get(acmg_num, str(raw_acmg or ""))
    chrom = str(variant.get("CHROM") or "")
    if chrom and not chrom.startswith("chr"):
        chrom = "chr" + chrom
    start = variant.get("POS")
    end = variant.get("END")
    coords = f"{chrom}:{start}-{end}" if chrom and start is not None and end is not None else ""
    sv_type = str(variant.get("sv_type") or "").upper()
    return ", ".join(x for x in (f"{coords} {sv_type}".strip(), acmg) if x)


def _scan_snv_review_rows_by_ids(tsv_path: Path, wanted: set[str]) -> list[dict[str, str]]:
    if not wanted or not tsv_path.exists():
        return []
    import csv as _csv
    out = []
    with tsv_path.open("r", encoding="utf-8", newline="") as f:
        for row in _csv.DictReader(f, delimiter="\t"):
            try:
                variant = _row_to_variant(row)
            except (KeyError, TypeError, ValueError):
                continue
            if variant.get("id") in wanted:
                out.append(row)
                if len(out) >= len(wanted):
                    break
    return out


def _case_snv_variants_by_id(sample_dir: Path, wanted: set[str]) -> dict[str, dict]:
    sample_id = _sample_id_from_state_dir(sample_dir)
    raw_tsv = sample_layout.snv_raw_tsv(sample_id)
    if not wanted or not raw_tsv.exists():
        return {}
    rows = snv_gene_index.query_rows_by_ids(
        raw_tsv, wanted, sample_layout.snv_gene_index_path(sample_id)
    )
    if rows is None:
        rows = _scan_snv_review_rows_by_ids(
            sample_layout.review_tsv(sample_id),
            wanted,
        )
        # Small WES files are acceptable as a last-resort fallback. Avoid
        # multi-GB DRAGEN/WGS raw scans on the case-list path.
        try:
            raw_size = raw_tsv.stat().st_size
        except OSError:
            raw_size = 0
        if not rows and raw_size and raw_size < 100 * 1024 * 1024:
            rows = _scan_snv_review_rows_by_ids(raw_tsv, wanted)
    variants: dict[str, dict] = {}
    with snv_overlay.OverlayReader(
        raw_tsv, sample_layout.snv_overlay_path(sample_id)
    ) as overlay:
        for row in overlay.apply_many(rows or []):
            try:
                variant = _row_to_variant(row)
            except (KeyError, TypeError, ValueError):
                continue
            if variant.get("id") in wanted:
                variants[variant["id"]] = variant
    if variants:
        _enrich_snv_variants(variants, sample_id, sample_dir)
        meta = _read_json_or(
            sample_layout.state_file(sample_id, "sample_metadata.json"), {}
        ) or {}
        try:
            current = manual_acmg.bulk_current(
                meta.get("genome_build") or "hg38", variants.keys()
            )
        except (OSError, ValueError, sqlite3.Error):
            current = {}
        for variant_id, variant in variants.items():
            assertion = current.get(manual_acmg.normalize_variant_id(variant_id))
            if assertion and (
                assertion.get("reusable_criteria")
                or not assertion.get("criteria")
            ):
                variant["effective_acmg_class"] = assertion[
                    "reusable_classification"
                ]
                variant["effective_acmg_score"] = assertion["reusable_score"]
                variant["effective_acmg_criteria"] = assertion[
                    "reusable_criteria_text"
                ]
                variant["effective_acmg_source"] = "manual"
                variant["effective_acmg_vus_subclass"] = (
                    assertion.get("reusable_vus_subclass")
                    or manual_acmg.vus_subclass(
                        assertion["reusable_classification"],
                        assertion["reusable_score"],
                    )
                )
    return variants


def _case_mito_variants_by_id(sample_dir: Path, wanted: set[str]) -> dict[str, dict]:
    """Read marked mitochondrial variants from the mito adapter, not SNV TSV."""
    sample_id = _sample_id_from_state_dir(sample_dir)
    mito_tsv = sample_layout.mito_tsv(sample_id)
    if not wanted or not mito_tsv.exists():
        return {}
    try:
        from ..adapters.mito_tsv import load_mito_tsv
        sidecar_dir = analyses_store.active_version_dir(sample_id)
        pheno_by_gene = _read_pheno_scores(
            _sidecar_file(sample_id, sidecar_dir, "pheno_score.tsv")
        )
        variants, _categories = load_mito_tsv(mito_tsv, pheno_by_gene=pheno_by_gene)
    except Exception as e:
        print(f"[case-summary] mito lookup failed for {sample_id}: {e}", flush=True)
        return {}
    return {vid: variant for vid, variant in variants.items() if vid in wanted}


def _variant_from_merged_id(variant_id: str) -> dict | None:
    parts = str(variant_id or "").split("-")
    if len(parts) < 6 or parts[0] != "MERGED":
        return None
    source = parts[1].lower()
    if source not in ("cnv", "sv"):
        return None
    sv_type = parts[-1]
    try:
        start = int(parts[-3])
        end = int(parts[-2])
    except ValueError:
        return None
    chrom = "-".join(parts[2:-3])
    if not chrom:
        return None
    return {
        "id": variant_id,
        "source": source,
        "CHROM": chrom,
        "POS": start,
        "END": end,
        "sv_type": sv_type,
    }


def _case_cnv_sv_variants_by_id(sample_dir: Path, wanted: set[str]) -> dict[str, dict]:
    if not wanted:
        return {}
    variants: dict[str, dict] = {}
    unresolved: set[str] = set()
    for variant_id in wanted:
        merged = _variant_from_merged_id(variant_id)
        if merged:
            variants[variant_id] = merged
        else:
            unresolved.add(variant_id)
    if unresolved:
        from ..adapters.annotsv_tsv import load_annotsv_variants_by_ids
        sample_id = _sample_id_from_state_dir(sample_dir)
        for source, tsv_path in (
            ("cnv", sample_layout.cnv_tsv(sample_id)),
            ("sv", sample_layout.sv_tsv(sample_id)),
        ):
            hits = load_annotsv_variants_by_ids(
                tsv_path,
                source=source,
                ids=unresolved,
            )
            variants.update(hits)
            unresolved.difference_update(hits.keys())
            if not unresolved:
                break
    merged_by_source: dict[str, list[dict]] = {"cnv": [], "sv": []}
    for variant in variants.values():
        if str(variant.get("id") or "").startswith("MERGED-"):
            merged_by_source.setdefault(str(variant.get("source") or "cnv"), []).append(variant)
    for source, merged_variants in merged_by_source.items():
        if not merged_variants:
            continue
        sample_id = _sample_id_from_state_dir(sample_dir)
        _augment_merged_cnv_sv_acmg(
            sample_layout.cnv_tsv(sample_id)
            if source == "cnv" else sample_layout.sv_tsv(sample_id),
            merged_variants,
        )
    return variants


def _augment_merged_cnv_sv_acmg(tsv_path: Path, merged_variants: list[dict]) -> None:
    """Fill ACMG for merged parents by scanning matching full rows only."""
    if not tsv_path.exists() or not merged_variants:
        return
    import csv as _csv
    with tsv_path.open("r", encoding="utf-8", newline="") as f:
        reader = _csv.DictReader(f, delimiter="\t")
        for row in reader:
            if (row.get("Annotation_mode") or "").strip() != "full":
                continue
            chrom = str(row.get("SV_chrom") or "")
            sv_type = str(row.get("SV_type") or "").upper()
            try:
                start = int(float(row.get("SV_start") or 0))
                end = int(float(row.get("SV_end") or 0))
                acmg = int(float((row.get("ACMG_class") or "").split("=")[-1]))
            except (TypeError, ValueError):
                continue
            for variant in merged_variants:
                if str(variant.get("CHROM") or "") != chrom:
                    continue
                if str(variant.get("sv_type") or "").upper() != sv_type:
                    continue
                v_start = int(variant.get("POS") or 0)
                v_end = int(variant.get("END") or 0)
                if v_start <= start <= v_end and v_start <= end <= v_end:
                    current = variant.get("acmg_class")
                    if not isinstance(current, int) or acmg > current:
                        variant["acmg_class"] = acmg


def _case_management_summary(
    sample_dir: Path,
    meta: dict,
    *,
    omim_sig: tuple | None = None,
) -> dict[str, str]:
    """Summarize marked SNV/CNV/SV variants for the case-list modal."""
    signature = _case_summary_signature(sample_dir, omim_sig=omim_sig)
    key = tuple(tuple(item) for item in signature)
    with _case_summary_cache_lock:
        cached = _case_summary_cache.get(key)
        if cached is not None:
            _case_summary_cache.move_to_end(key)
            return cached
    disk_cached = _read_case_summary_disk(sample_dir, signature)
    if disk_cached is not None:
        with _case_summary_cache_lock:
            _case_summary_cache[key] = disk_cached
            _case_summary_cache.move_to_end(key)
            while len(_case_summary_cache) > _CASE_SUMMARY_CACHE_MAX:
                _case_summary_cache.popitem(last=False)
        return disk_cached

    statuses = meta.get("status") or {}
    wanted = {
        vid for vid, status in statuses.items()
        if str(status).strip() in ("1", "2")
    }
    edits_by_id = meta.get("edits") or {}
    causative: list[str] = []
    diseases: list[str] = []
    other: list[str] = []
    if wanted:
        omim_store.ensure_loaded()

    mito_variants = _case_mito_variants_by_id(sample_dir, wanted)
    for vid, variant in mito_variants.items():
        edits = edits_by_id.get(vid) or {}
        label = _case_mito_label(variant, edits)
        if str(statuses.get(vid, "")).strip() == "1":
            if label:
                causative.append(label)
            for disease in _case_selected_mito_diseases(variant, edits):
                if disease not in diseases:
                    diseases.append(disease)
        elif label:
            other.append(label)

    snv_wanted = wanted.difference(mito_variants.keys())
    snv_variants = _case_snv_variants_by_id(sample_dir, snv_wanted)
    for vid, variant in snv_variants.items():
        edits = edits_by_id.get(vid) or {}
        label = _case_variant_label(variant, edits)
        if str(statuses.get(vid, "")).strip() == "1":
            if label:
                causative.append(label)
            for disease in _case_selected_diseases(variant, edits):
                if disease not in diseases:
                    diseases.append(disease)
        elif label:
            other.append(label)

    cnv_sv_wanted = wanted.difference(mito_variants.keys()).difference(snv_variants.keys())
    cnv_sv_variants = _case_cnv_sv_variants_by_id(sample_dir, cnv_sv_wanted)
    for vid, variant in cnv_sv_variants.items():
        edits = edits_by_id.get(vid) or {}
        label = _case_cnv_sv_label(variant, edits)
        if str(statuses.get(vid, "")).strip() == "1":
            if label:
                causative.append(label)
            disease = str(edits.get("disease") or "").strip()
            if disease and disease not in diseases:
                diseases.append(disease)
        elif label:
            other.append(label)

    result = {
        "causative_variants": "\n".join(causative),
        "diseases": "\n".join(diseases),
        "other_variants": "\n".join(other),
        "comment": str(meta.get("comment") or ""),
        "sign_received_at": str(meta.get("sign_received_at") or ""),
    }
    with _case_summary_cache_lock:
        _case_summary_cache[key] = result
        _case_summary_cache.move_to_end(key)
        while len(_case_summary_cache) > _CASE_SUMMARY_CACHE_MAX:
            _case_summary_cache.popitem(last=False)
    _write_case_summary_disk(sample_dir, signature, result)
    return result


def _load_enriched_snv_cached(
    snv_tsv: Path,
    *,
    sample_id: str,
    sidecar_dir: Path,
    test_type: str,
) -> tuple[dict, dict, dict, str]:
    """Parse + enrich SNVs once per input revision, with a bounded LRU."""
    started = time.perf_counter()
    exo_path = _sidecar_file(sample_id, sidecar_dir, "exomiser_results.tsv")
    lir_path = _sidecar_file(sample_id, sidecar_dir, "lirical_results.tsv")
    pheno_path = _sidecar_file(sample_id, sidecar_dir, "pheno_score.tsv")
    analysis_path = _sidecar_file(sample_id, sidecar_dir, "analysis.json")
    omim_sig = omim_store.cache_signature()
    gene_disease_sig = gene_disease_store.cache_signature()
    key = (
        _file_signature(snv_tsv),
        _file_signature(analysis_path),
        _file_signature(pheno_path),
        _file_signature(exo_path),
        _file_signature(lir_path),
        (test_type or "WES").upper(),
        omim_sig,
        gene_disease_sig,
    )
    cache_allowed = _snv_cache_allowed(snv_tsv)
    if cache_allowed:
        with _snv_cache_lock:
            cached = _snv_cache.get(key)
            if cached is not None:
                _snv_cache.move_to_end(key)
                _log_perf(
                    "snv.enriched_cache",
                    started,
                    cache="hit",
                    tsv=snv_tsv.name,
                    size=_fmt_size(snv_tsv),
                    variants=len(cached[0]),
                    max_entries=SNV_CACHE_MAX,
                )
                return cached

    old_format_error = ""
    parse_started = time.perf_counter()
    if snv_tsv.exists():
        try:
            variants, categories = load_snv_tsv(snv_tsv, test_type=test_type)
        except OldFormatError as e:
            variants, categories = {}, {t: [] for t in TIERS}
            old_format_error = str(e)
    else:
        variants, categories = {}, {t: [] for t in TIERS}
    parse_elapsed = time.perf_counter() - parse_started

    enrich_started = time.perf_counter()
    pheno_by_gene = _enrich_snv_variants(variants, sample_id, sidecar_dir)
    for t, ids in categories.items():
        categories[t] = sorted(ids, key=lambda i: (-_variant_total_score(variants, i), i))

    result = (variants, categories, pheno_by_gene, old_format_error)
    if cache_allowed:
        with _snv_cache_lock:
            # Callers treat the enriched maps as read-only; keeping the same
            # objects avoids a full deep-copy cost on cached review payloads.
            _snv_cache[key] = result
            _snv_cache.move_to_end(key)
            while len(_snv_cache) > SNV_CACHE_MAX:
                _snv_cache.popitem(last=False)
    _log_perf(
        "snv.enriched_cache",
        started,
        cache="miss" if cache_allowed else "skip",
        tsv=snv_tsv.name,
        size=_fmt_size(snv_tsv),
        variants=len(variants),
        parse=f"{parse_elapsed:.3f}s",
        enrich=f"{time.perf_counter() - enrich_started:.3f}s",
        max_entries=SNV_CACHE_MAX,
        raw_cache_limit_mb=SNV_CACHE_MAX_RAW_MB,
    )
    return result


def _variant_total_score(variants: dict, vid: str) -> float:
    ts = variants.get(vid, {}).get("total_score")
    return float(ts) if isinstance(ts, (int, float)) else float("-inf")


def _is_clinvar_plp(variant: dict) -> bool:
    sig = (variant.get("CLNSIG") or "").strip().lower().replace("_", " ")
    return sig in {
        "pathogenic",
        "likely pathogenic",
        "pathogenic/likely pathogenic",
        "likely pathogenic/pathogenic",
    }


def _secondary_panel_genes() -> dict[str, set[str]]:
    from . import phenotype_scorer
    out: dict[str, set[str]] = {}
    for category, panel_name in SECONDARY_SNV_PANELS.items():
        payload = phenotype_scorer.genes_for_key(panel_name, kind="panel")
        out[category] = {
            panel_deadzone.canonical_panel_gene_symbol(g)
            for g in payload.get("genes", [])
            if g
        }
        out[category].discard("")
    return out


def _secondary_variant_gene(variant: dict) -> str:
    gene, _hgnc_id = panel_deadzone.canonical_gene_symbol(
        variant.get("gene_symbol") or variant.get("GENE") or "",
        variant.get("HGNC_ID") or "",
    )
    return gene


def _secondary_vaf_pass(variant: dict) -> bool:
    vaf = _to_num(variant.get("alt_af"))
    return isinstance(vaf, (int, float)) and vaf >= 0.2


def _secondary_zygosity_pass(variant: dict) -> bool:
    zygosity = str(
        variant.get("zygosity")
        or variant.get("ZYGOSITY")
        or ""
    ).strip().lower().replace("_", " ").replace("-", " ")
    if not zygosity:
        return True
    return zygosity not in {
        "ref",
        "reference",
        "hom ref",
        "homozygous reference",
        "0/0",
        "0|0",
    }


def _is_secondary_snv_candidate(variant: dict) -> bool:
    # Keep every existing ClinVar P/LP candidate, then broaden retrieval with
    # the same 1A/1B/1C buckets used by the main SNV/Indel card.  Tier 1B is
    # LOFTEE HC; tier 1C covers ACMG points >=4, P-KNN LLR >=1, and the other
    # configured predictor triggers. Selection into the report remains
    # ClinVar-only by default.
    tier = str(variant.get("tier") or "").strip().upper()
    return (
        (_is_clinvar_plp(variant) or tier in {"1A", "1B", "1C"})
        and _secondary_vaf_pass(variant)
        and _secondary_zygosity_pass(variant)
    )


def _build_secondary_snv_categories(variants: dict[str, dict]) -> dict[str, list[str]]:
    panel_genes = _secondary_panel_genes()
    categories: dict[str, list[str]] = {}
    for category, genes in panel_genes.items():
        ids = [
            vid for vid, variant in variants.items()
            if _secondary_variant_gene(variant) in genes
            and _is_secondary_snv_candidate(variant)
        ]
        categories[category] = sorted(
            set(ids),
            key=lambda vid: (
                0 if _is_clinvar_plp(variants.get(vid, {})) else 1,
                -_variant_total_score(variants, vid),
                vid,
            ),
        )
    return categories


def _pipeline_clinvar_secondary_variants(
    variants: dict[str, dict],
) -> dict[str, dict]:
    """Restore report-baseline ClinVar before health candidate selection."""
    restored = clinvar_latest_store.restore_pipeline_variants(variants)
    for variant in restored.values():
        # A weekly upgrade may have assigned 1A. Once the 2026-07-20
        # classification is restored, retain only independent 1B/1C evidence.
        if str(variant.get("tier") or "").upper() != "1A" or _is_clinvar_plp(variant):
            continue
        if variant.get("loftee_hc"):
            variant["tier"] = "1B"
            continue
        points = _to_num(
            variant.get("effective_acmg_score")
            if variant.get("effective_acmg_score") is not None
            else variant.get("ACMG_score")
        )
        variant["tier"] = (
            "1C"
            if (
                isinstance(points, (int, float)) and points >= 4
            ) or variant.get("predicted_suspect_non_acmg")
            else "2"
        )
    return restored


def _enrich_snv_variants(
    variants: dict[str, dict],
    sample_id: str,
    sidecar_dir: Path,
) -> dict[str, float]:
    """Join pheno / Exomiser / LIRICAL / OMIM data into variant payloads."""
    exo_path = _sidecar_file(sample_id, sidecar_dir, "exomiser_results.tsv")
    lir_path = _sidecar_file(sample_id, sidecar_dir, "lirical_results.tsv")
    pheno_path = _sidecar_file(sample_id, sidecar_dir, "pheno_score.tsv")
    exo = _read_tsv_dict(exo_path)
    lir = _read_tsv_dict(lir_path)
    pheno_by_gene = _read_pheno_scores(pheno_path)

    def _scale_to_100(s):
        n = _to_num(s)
        if not isinstance(n, (int, float)):
            return None
        return int(round(n * 100))

    # Refresh OMIM and supplemental disease data once per batch; individual
    # variant joins use only in-memory dict lookups.
    omim_store.ensure_loaded()
    gene_disease_store.ensure_loaded()
    for vid, v in variants.items():
        gene = v.get("gene_symbol", "")
        pheno = pheno_by_gene.get(gene) if gene else None
        if pheno and pheno > 0:
            v["pheno_score"] = round(pheno, 2)
        v["in_panel"] = bool(pheno and pheno > 0)
        gs = v.get("geno_score")
        ps = v.get("pheno_score")
        if gs is not None or ps is not None:
            v["total_score"] = (gs or 0) + (ps or 0)
        e = exo.get(vid)
        if e:
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
        omim_id = omim_store.parse_omim_id_from_link(v.get("OMIM_link", ""))
        rec = omim_store.lookup_cached(omim_id=omim_id, gene=gene)
        if rec:
            v["OMIM_id"]      = rec.get("OMIM_id", "")
            v["OMIM_disease"] = rec.get("OMIM_disease", "")
            v["Inheritance"]  = rec.get("Inheritance", "")
            for f in omim_store.DISEASE_FIELDS:
                v[f] = rec.get(f, "")
        v["disease_associations"] = gene_disease_store.merged_associations(gene, rec, refresh=False)
    litvar2_on_demand.apply_cached(variants, sample_id)
    return pheno_by_gene


def _legacy_manual_snapshot(edit: dict) -> dict | None:
    """Read pre-modal free-text ACMG edits as a per-sample snapshot."""
    if not isinstance(edit, dict):
        return None
    if isinstance(edit.get("manual_acmg"), dict):
        return edit["manual_acmg"]
    values = (
        edit.get("ACMG_classification"),
        edit.get("ACMG_score"),
        edit.get("ACMG_criteria"),
    )
    if not any(value not in (None, "") for value in values):
        return None
    parsed, unknown = manual_acmg.parse_criteria_text(edit.get("ACMG_criteria"))
    return {
        "criteria": parsed,
        "criteria_text": str(edit.get("ACMG_criteria") or ""),
        "score": edit.get("ACMG_score"),
        "classification": str(edit.get("ACMG_classification") or ""),
        "reviewer_username": "",
        "source_sample_id": "",
        "created_at": None,
        "legacy": True,
        "unrecognized_criteria": unknown,
    }


def _apply_effective_acmg(
    variants: dict[str, dict],
    categories: dict[str, list[str]],
    *,
    sample_id: str,
    genome_build: str,
    meta: dict,
) -> None:
    """Bulk overlay per-sample/global manual ACMG and observation counts."""
    if not variants:
        return
    started = time.perf_counter()
    try:
        current = manual_acmg.bulk_current(genome_build, variants.keys())
        observed = manual_acmg.bulk_observed_counts(
            genome_build,
            variants.keys(),
            exclude_sample_id=sample_id,
        )
    except (OSError, ValueError, sqlite3.Error):
        current = {}
        observed = {}
    edits = meta.get("edits") if isinstance(meta.get("edits"), dict) else {}
    for variant_id, original in list(variants.items()):
        variant = dict(original)
        variants[variant_id] = variant
        sample_snapshot = _legacy_manual_snapshot(edits.get(variant_id) or {})
        global_assertion = current.get(manual_acmg.normalize_variant_id(variant_id))

        if sample_snapshot:
            source = "manual"
            selected = sample_snapshot
            source_scope = "sample"
        elif global_assertion and (
            global_assertion.get("reusable_criteria")
            or not global_assertion.get("criteria")
        ):
            source = "manual"
            selected = {
                "classification": global_assertion.get(
                    "reusable_classification", ""
                ),
                "score": global_assertion.get("reusable_score"),
                "criteria_text": global_assertion.get(
                    "reusable_criteria_text", ""
                ),
            }
            source_scope = "global"
        elif variant.get("genebe_acmg_class") not in (None, ""):
            source = "GeneBe"
            selected = {
                "classification": variant.get("genebe_acmg_class") or "",
                "score": variant.get("genebe_acmg_score"),
                "criteria_text": variant.get("genebe_acmg_criteria") or "",
            }
            source_scope = "genebe"
        else:
            source = "in-house"
            selected = {
                "classification": variant.get("ACMG_classification") or "",
                "score": variant.get("ACMG_score"),
                "criteria_text": variant.get("ACMG_criteria") or "",
            }
            source_scope = "in_house"

        variant["sample_acmg_snapshot"] = sample_snapshot
        variant["manual_acmg_current"] = global_assertion
        variant["effective_acmg_source"] = source
        variant["effective_acmg_scope"] = source_scope
        variant["effective_acmg_class"] = selected.get("classification") or ""
        variant["effective_acmg_score"] = selected.get("score")
        variant["effective_acmg_criteria"] = selected.get("criteria_text") or ""
        variant["effective_acmg_vus_subclass"] = manual_acmg.vus_subclass(
            selected.get("classification"), selected.get("score")
        )
        variant["observed_count"] = int(
            observed.get(manual_acmg.normalize_variant_id(variant_id), 0)
        )

        geno_score = manual_acmg.acmg_to_variant_score(selected.get("score"))
        variant["geno_score"] = geno_score
        pheno_score = variant.get("pheno_score")
        if geno_score is not None or pheno_score is not None:
            variant["total_score"] = (geno_score or 0) + (pheno_score or 0)
        else:
            variant.pop("total_score", None)

        original_tier = str(variant.get("tier") or "2").upper()
        if original_tier in {"1A", "1B"}:
            tier = original_tier
        else:
            try:
                points_trigger = float(selected.get("score")) >= 4
            except (TypeError, ValueError):
                points_trigger = False
            tier = (
                "1C"
                if points_trigger or variant.get("predicted_suspect_non_acmg")
                else "2"
            )
        variant["tier"] = tier

    for tier in ("1A", "1B", "1C", "2"):
        categories[tier] = sorted(
            (
                variant_id for variant_id, variant in variants.items()
                if variant.get("tier") == tier
            ),
            key=lambda variant_id: (
                -_variant_total_score(variants, variant_id),
                variant_id,
            ),
        )
    _log_perf(
        "sample.manual_acmg_overlay",
        started,
        sample=sample_id,
        variants=len(variants),
        manual_matches=len(current),
        observed_matches=sum(1 for count in observed.values() if count),
    )


def _variants_from_rows(rows: list[dict[str, str]], *, test_type: str) -> dict[str, dict]:
    """Shape raw TSV rows into variant payloads without parsing the whole TSV."""
    out: dict[str, dict] = {}
    is_wes = (test_type or "").upper() == "WES"
    low_flag = _DP_LOW_FLAG_WES if is_wes else _DP_LOW_FLAG_WGS
    for row in rows:
        if _is_mito_chrom(row.get("CHROM") or ""):
            continue
        canonical_gene, _ = panel_deadzone.canonical_gene_symbol(
            row.get("GENE", ""),
            row.get("HGNC_ID", ""),
        )
        if canonical_gene in EXCLUDED_GENES:
            continue
        v = _row_to_variant(row)
        dp = v.get("depth")
        if is_wes and (dp is None or dp < _DP_HARD_FLOOR_WES):
            continue
        v["low_depth"] = bool(dp is not None and dp < low_flag)
        merge_snv_variant_row(out, v)
    return out


def _variant_id_from_row(row: dict[str, str]) -> str:
    chrom = (row.get("CHROM") or "").strip()
    pos = (row.get("POS") or "").strip()
    ref = (row.get("REF") or "").strip()
    alt = (row.get("ALT") or "").strip()
    if not (chrom and pos and ref and alt):
        return ""
    return "-".join((chrom, pos, ref, alt))


def _scan_raw_snv_rows_by_ids(raw_tsv: Path, wanted_ids: set[str]) -> list[dict[str, str]]:
    """Fallback exact-id lookup for reviewer-marked rows when the index is absent.

    Gene search can still find a variant through the raw TSV fallback, so the
    report view must be able to rehydrate the same reviewer-marked variant even
    on old samples that predate snv_gene_index.sqlite. The wanted set is tiny
    (status=1/2/C ids), and we keep scanning after the first hit so multi-
    transcript rows for the same genomic variant merge into one complete card.
    """
    import csv as _csv

    wanted = {str(vid).strip() for vid in wanted_ids if str(vid).strip()}
    if not wanted or not raw_tsv.is_file():
        return []
    out: list[dict[str, str]] = []
    with raw_tsv.open("r", encoding="utf-8", newline="") as f:
        for row in _csv.DictReader(f, delimiter="\t"):
            if is_reportable_raw_row(row) and _variant_id_from_row(row) in wanted:
                out.append(row)
    return out


def _scan_raw_snv_rows_by_genes(raw_tsv: Path, genes: set[str]) -> list[dict[str, str]]:
    """Streaming fallback for gene search when snv_gene_index.sqlite is absent.

    This is slower than the SQLite index, but it keeps memory bounded by the
    number of matching rows instead of materializing a multi-GB WGS TSV.
    """
    import csv as _csv

    wanted = {str(g).strip().upper() for g in genes if str(g).strip()}
    if not wanted or not raw_tsv.is_file():
        return []
    out: list[dict[str, str]] = []
    with raw_tsv.open("r", encoding="utf-8", newline="") as f:
        reader = _csv.DictReader(f, delimiter="\t")
        for row in reader:
            if not is_reportable_raw_row(row):
                continue
            gene, _hgnc_id = panel_deadzone.canonical_gene_symbol(
                row.get("GENE") or "",
                row.get("HGNC_ID") or "",
            )
            if gene.upper() in wanted:
                out.append(row)
    return out


def _supplement_marked_snv_variants(
    variants: dict[str, dict],
    categories: dict[str, list[str]],
    raw_tsv: Path,
    *,
    sample_id: str,
    keep_ids: set[str],
    sidecar_dir: Path,
    test_type: str,
    index_path: Path | None = None,
    overlay_path: Path | None = None,
) -> None:
    """Add reviewer-marked SNVs absent from review.tsv via the gene index.

    This avoids rebuilding the compact review TSV every time status/edit
    metadata changes. If the index is absent/stale, skip the supplement;
    the main review TSV still renders immediately.
    """
    missing = {vid for vid in keep_ids if vid and vid not in variants}
    if not missing:
        return
    started = time.perf_counter()
    rows = snv_gene_index.query_rows_by_ids(raw_tsv, missing, index_path)
    if rows is None:
        rows = _scan_raw_snv_rows_by_ids(raw_tsv, missing)
        _log_perf(
            "sample.marked_snv_supplement",
            started,
            status="raw_scan_fallback",
            missing=len(missing),
            rows=len(rows),
        )
        if not rows:
            return
    with snv_overlay.OverlayReader(raw_tsv, overlay_path) as overlay:
        rows = overlay.apply_many(rows)
    extra = _variants_from_rows(rows, test_type=test_type)
    _enrich_snv_variants(extra, sample_id, sidecar_dir)
    for vid, variant in extra.items():
        variants[vid] = variant
        tier = variant.get("tier")
        if tier in categories and vid not in categories[tier]:
            categories[tier].append(vid)
    for tier, ids in categories.items():
        categories[tier] = sorted(set(ids), key=lambda i: (-_variant_total_score(variants, i), i))
    _log_perf(
        "sample.marked_snv_supplement",
        started,
        status="ok",
        missing=len(missing),
        added=len(extra),
    )


def list_index() -> list[dict]:
    """Return the sample list for the top-bar combobox.

    Source of truth is a directory scan: every subdirectory of
    tertiary_output/ that has sample_metadata.json shows up here,
    sorted by created_at descending so newly-registered samples land
    at the top. This endpoint is intentionally lightweight for browser
    boot: case-management summaries are computed by list_case_summaries().
    """
    out: list[dict] = []
    for sample_id in sorted(sample_layout.iter_sample_ids()):
        sub = sample_layout.state_dir(sample_id)
        meta_path = sample_layout.state_file(sample_id, "sample_metadata.json")
        if not meta_path.exists():
            # Pipeline-dropped sample that hasn't been registered yet —
            # surfaces through /samples/unregistered, not the search bar.
            continue
        meta = _read_json_or(meta_path, {}) or {}
        out.append({
            "sample_id":     meta.get("sample_id") or sample_id,
            "lis_id":        meta.get("lis_id") or sample_id,
            "name":          meta.get("name", ""),
            "mrn":           meta.get("mrn", ""),
            "sex":           meta.get("sex", ""),
            "test_type":     _effective_test_type(meta, sample_id),
            "category":      meta.get("category", ""),
            "run_date":      meta.get("run_date", ""),
            "created_at":    meta.get("created_at", ""),
            "tags":          meta.get("tags", []),
            "has_completed": sample_layout.snv_raw_tsv(sample_id).is_file(),
            "tertiary_dir":  sample_id,
        })
    # Sort newest-first by registration date; samples without a stored
    # created_at fall to the bottom (stable order by lis_id thereafter).
    out.sort(key=lambda r: (r.get("created_at") or "", r.get("lis_id") or ""), reverse=True)

    # Best-effort cache write so an operator browsing the filesystem
    # sees an up-to-date listing. Failures are non-fatal — the read
    # path doesn't depend on this file existing.
    try:
        cache_path = sample_layout.global_cache_path("_index.json")
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(
            json.dumps(out, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
    except OSError:
        pass

    return out


def list_case_summaries() -> list[dict]:
    """Return the denormalized case-list table rows.

    The modal should not scan variant TSVs on every open. Normal app writes
    refresh one row at a time; this path rebuilds only when the table is absent
    or the registered sample set changed outside the app.
    """
    payload = _read_case_table_payload()
    rows = payload.get("rows") if payload else []
    if rows and _case_table_is_complete(rows):
        rows = [
            {
                **row,
                "test_type": test_types.normalize_test_type(
                    row.get("test_type") or "",
                    sample_id=row.get("lis_id") or row.get("sample_id") or "",
                    default="",
                ),
            }
            for row in rows
        ]
        rows.sort(key=lambda r: (r.get("created_at") or "", r.get("lis_id") or ""), reverse=True)
        return rows
    return rebuild_case_table()


def list_unregistered() -> list[dict]:
    """Sample dirs the tertiary pipeline left behind without a metadata file.

    A unified sample counts as "unregistered" after layout.json activation
    when it has a 03_acmg raw TSV but no sample_metadata.json. Legacy UI
    copies remain supported during migration. The UI
    surfaces these in the 載入新個案 dropdown so the reviewer can attach
    basic info + HPO without having to retype the LIS_ID.

    Per-entry enrichment:
      * `roster` — {mrn, name, test_type, department} from the uploaded
        未完成報告清單 (patient_list_store). This is the preferred
        source for MRN + 姓名 + Test type in the modal.
      * `phenotype` — parsed HPO/panels from
        NGS_UI/patient_phenotype/{lis_id}_{mrn}_phenotype.txt, used for
        the preview chip row. (Its filename-derived MRN is kept as a
        fallback for samples not in the roster.)

    Sorted by directory mtime descending (newest first).
    """
    from ..config import PHENOTYPE_DIR
    from . import dragen_jobs, patient_list_store, phenotype_io
    roster = patient_list_store.load_roster()
    active_tertiary_samples = dragen_jobs.active_sample_ids()
    out: list[dict] = []
    for lis_id in sample_layout.iter_sample_ids():
        if not sample_layout.is_ui_ready(lis_id):
            continue
        sub = sample_layout.state_dir(lis_id)
        pipeline_sub = sample_layout.pipeline_sample_dir(lis_id)
        tsv = sample_layout.snv_raw_tsv(lis_id)
        meta = sample_layout.state_file(lis_id, "sample_metadata.json")
        if not tsv.exists() or meta.exists():
            continue
        if lis_id in active_tertiary_samples:
            continue

        source_info = {}
        source_vcf_path = ""
        source_vcf_size = 0
        pipeline_source = sample_layout.state_file(lis_id, "pipeline_source.json")
        if pipeline_source.is_file():
            try:
                source_info = json.loads(pipeline_source.read_text(encoding="utf-8")) or {}
            except (OSError, json.JSONDecodeError):
                source_info = {}
            source_vcf_path = str(source_info.get("source_vcf_path") or "")
            if source_vcf_path:
                try:
                    source_vcf_size = Path(source_vcf_path).stat().st_size
                except OSError:
                    source_vcf_size = 0

        source_sample_id = str(source_info.get("source_sample_id") or "")
        roster_entry, roster_lis_id = patient_list_store.lookup_with_key(
            lis_id,
            source_sample_id,
            roster=roster,
        )
        roster_mrn = (roster_entry or {}).get("mrn") or ""
        pheno_lis_ids = patient_list_store.lookup_candidates(lis_id, roster_lis_id, source_sample_id)

        # Resolve a matching phenotype file in the central phenotype dir.
        # Lookup order:
        #   1. {lis_id}_*_phenotype.txt, also trying the roster/source
        #      LIS_ID when the UI directory carries a caller suffix.
        #      If the roster gives an MRN, prefer exact {lis_id}_{mrn}.
        #   2. {mrn}_phenotype.txt       — the standalone HPO tool's
        #      MRN-only output (no LIS_ID).
        # Either way we still parse the file for the HPO/panel preview.
        # Lookup priority:
        #   1. {candidate_lis_id}_{roster_mrn}_phenotype.txt
        #   2. {candidate_lis_id}_phenotype.txt
        #   3. {candidate_lis_id}_*_phenotype.txt
        #   4. {roster_mrn}_phenotype.txt           — MRN-only file
        pheno_payload = None
        if PHENOTYPE_DIR.is_dir():
            pf = None
            if roster_mrn:
                for pheno_lis_id in pheno_lis_ids:
                    exact = PHENOTYPE_DIR / f"{pheno_lis_id}_{roster_mrn}_phenotype.txt"
                    if exact.is_file():
                        pf = exact
                        break
            if pf is None:
                for pheno_lis_id in pheno_lis_ids:
                    lis_only = PHENOTYPE_DIR / f"{pheno_lis_id}_phenotype.txt"
                    if lis_only.is_file():
                        pf = lis_only
                        break
            if pf is None:
                for pheno_lis_id in pheno_lis_ids:
                    lis_matches = sorted(PHENOTYPE_DIR.glob(f"{pheno_lis_id}_*_phenotype.txt"))
                    if lis_matches:
                        pf = lis_matches[0]
                        break
            if pf is None and roster_mrn:
                mrn_only = PHENOTYPE_DIR / f"{roster_mrn}_phenotype.txt"
                if mrn_only.is_file():
                    pf = mrn_only
            if pf is not None:
                # Recover the MRN from the filename for the modal's
                # fallback. {lis_id}_phenotype → no MRN; {lis_id}_{mrn}
                # → the second segment; {mrn}_phenotype → the whole core.
                stem = pf.stem
                core = stem[:-len("_phenotype")] if stem.endswith("_phenotype") else stem
                matched_pheno_lis = next(
                    (pheno_lis_id for pheno_lis_id in pheno_lis_ids if core == pheno_lis_id),
                    "",
                )
                matched_pheno_prefix = next(
                    (
                        pheno_lis_id
                        for pheno_lis_id in pheno_lis_ids
                        if core.startswith(pheno_lis_id + "_")
                    ),
                    "",
                )
                if matched_pheno_lis:
                    file_mrn = ""
                elif matched_pheno_prefix:
                    file_mrn = core[len(matched_pheno_prefix) + 1:]
                elif core == roster_mrn:
                    file_mrn = roster_mrn
                else:
                    file_mrn = core.split("_")[-1] if "_" in core else core
                try:
                    hpo, panels = phenotype_io.parse(pf.read_text(encoding="utf-8"))
                except OSError:
                    hpo, panels = [], []
                pheno_payload = {
                    "path":   str(pf),
                    "mrn":    file_mrn,
                    "hpo":    hpo,
                    "panels": panels,
                }

        try:
            mtime = pipeline_sub.stat().st_mtime
        except OSError:
            mtime = 0.0
        out.append({
            "lis_id":     lis_id,
            "tsv_size":   tsv.stat().st_size if tsv.exists() else 0,
            "mtime":      mtime,
            "pipeline_type": str(source_info.get("pipeline_type") or ""),
            "source_sample_id": source_sample_id,
            "roster_lis_id": roster_lis_id,
            "source_vcf_path": source_vcf_path,
            "source_vcf_size": source_vcf_size,
            "phenotype":  pheno_payload,
            # Roster-sourced identifiers (None when the LIS_ID isn't on
            # any uploaded clinic list yet).
            "roster":     {
                "mrn":              (roster_entry or {}).get("mrn", ""),
                "name":             (roster_entry or {}).get("name", ""),
                "test_type":        test_types.normalize_test_type(
                    (roster_entry or {}).get("test_type", ""),
                    sample_id=lis_id,
                    default="",
                ),
                "department":       (roster_entry or {}).get("department", ""),
                "physician":        (roster_entry or {}).get("physician", ""),
                "sign_received_at": (roster_entry or {}).get("sign_received_at", ""),
            } if roster_entry else None,
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


def _read_pheno_by_gene(sample_id: str, sidecar_dir: Path) -> dict:
    """Read pheno_score.tsv → {gene_symbol: 0-100 score}. Empty if absent."""
    out: dict[str, float] = {}
    p = _sidecar_file(sample_id, sidecar_dir, "pheno_score.tsv")
    if not p.exists():
        return out
    import csv as _csv
    with p.open("r", encoding="utf-8", newline="") as f:
        for row in _csv.DictReader(f, delimiter="\t"):
            g = panel_deadzone.canonical_panel_gene_symbol(row.get("gene_symbol") or "")
            if not g:
                continue
            try:
                score = float(row.get("pheno_score") or 0)
                out[g] = max(out.get(g, 0.0), score)
            except ValueError:
                pass
    return out


def _load_pheno_context(sample_id: str, version: str | None):
    """Resolve the active analysis version + its pheno gene set for a
    sample, without parsing any variant TSV. Used by the staged
    CNV/SV and Mito loaders so they don't have to re-do the full
    load_sample. Returns (sub, sidecar_dir, hpo_list, panels_list,
    pheno_by_gene) or None when the sample dir is missing."""
    sub = sample_layout.state_dir(sample_id)
    if not sub.is_dir():
        return None
    meta = _read_json_or(sample_layout.state_file(sample_id, "sample_metadata.json"), {}) or {}
    chosen_version = _resolve_version(
        sample_id, requested=version, meta_active=meta.get("active_analysis"))
    sidecar_dir = (analyses_store.version_dir(sample_id, chosen_version)
                   if chosen_version is not None else sub)
    if chosen_version is not None:
        analysis = analyses_store.read_version(sample_id, chosen_version) or {}
        hpo_list    = analysis.get("hpo") or []
        panels_list = analysis.get("selected_panels") or []
    else:
        hpo_list    = meta.get("hpo") or meta.get("patient_phenotype") or []
        panels_list = meta.get("selected_panels") or []
    return sub, sidecar_dir, hpo_list, panels_list, _read_pheno_by_gene(sample_id, sidecar_dir)


def load_sample_cnv_sv(sample_id: str, version: str | None = None) -> dict | None:
    """Staged loader: just the CNV/SV side-channels for a sample.
    {cnv_variants, cnv_categories, sv_variants, sv_categories} or None."""
    started = time.perf_counter()
    ctx = _load_pheno_context(sample_id, version)
    from ..adapters.annotsv_tsv import load_annotsv_tsv, CNV_TIERS, SV_TIERS
    if ctx is None:
        return None
    sub, _sd, hpo_list, panels_list, pheno_by_gene = ctx
    from . import phenotype_scorer
    pheno_matched, pheno_total = phenotype_scorer.compute_pheno_match(hpo_list, panels_list)
    cnv_path = sample_layout.cnv_tsv(sample_id)
    sv_path  = sample_layout.sv_tsv(sample_id)
    cnv_variants, cnv_categories = (
        load_annotsv_tsv(cnv_path, source="cnv", pheno_by_gene=pheno_by_gene,
                         pheno_matched=pheno_matched, pheno_total=pheno_total)
        if cnv_path.exists() else ({}, {t: [] for t in CNV_TIERS})
    )
    sv_variants, sv_categories = (
        load_annotsv_tsv(sv_path, source="sv", pheno_by_gene=pheno_by_gene,
                         pheno_matched=pheno_matched, pheno_total=pheno_total)
        if sv_path.exists() else ({}, {t: [] for t in SV_TIERS})
    )
    _log_perf(
        "sample.cnv_sv",
        started,
        sample=sample_id,
        cnv_size=_fmt_size(cnv_path),
        cnv_variants=len(cnv_variants),
        sv_size=_fmt_size(sv_path),
        sv_variants=len(sv_variants),
    )
    return {
        "cnv_variants": cnv_variants, "cnv_categories": cnv_categories,
        "sv_variants": sv_variants,   "sv_categories": sv_categories,
    }


def load_sample_cnv(sample_id: str, version: str | None = None) -> dict | None:
    """Staged loader: CNV side-channel only."""
    started = time.perf_counter()
    ctx = _load_pheno_context(sample_id, version)
    from ..adapters.annotsv_tsv import load_annotsv_tsv, CNV_TIERS
    if ctx is None:
        return None
    sub, _sd, hpo_list, panels_list, pheno_by_gene = ctx
    from . import phenotype_scorer
    pheno_matched, pheno_total = phenotype_scorer.compute_pheno_match(hpo_list, panels_list)
    cnv_path = sample_layout.cnv_tsv(sample_id)
    cnv_variants, cnv_categories = (
        load_annotsv_tsv(cnv_path, source="cnv", pheno_by_gene=pheno_by_gene,
                         pheno_matched=pheno_matched, pheno_total=pheno_total)
        if cnv_path.exists() else ({}, {t: [] for t in CNV_TIERS})
    )
    _log_perf(
        "sample.cnv",
        started,
        sample=sample_id,
        size=_fmt_size(cnv_path),
        variants=len(cnv_variants),
    )
    return {"cnv_variants": cnv_variants, "cnv_categories": cnv_categories}


def load_sample_sv(sample_id: str, version: str | None = None) -> dict | None:
    """Staged loader: SV side-channel only."""
    started = time.perf_counter()
    ctx = _load_pheno_context(sample_id, version)
    from ..adapters.annotsv_tsv import load_annotsv_tsv, SV_TIERS
    if ctx is None:
        return None
    sub, _sd, hpo_list, panels_list, pheno_by_gene = ctx
    from . import phenotype_scorer
    pheno_matched, pheno_total = phenotype_scorer.compute_pheno_match(hpo_list, panels_list)
    sv_path = sample_layout.sv_tsv(sample_id)
    sv_variants, sv_categories = (
        load_annotsv_tsv(sv_path, source="sv", pheno_by_gene=pheno_by_gene,
                         pheno_matched=pheno_matched, pheno_total=pheno_total)
        if sv_path.exists() else ({}, {t: [] for t in SV_TIERS})
    )
    _log_perf(
        "sample.sv",
        started,
        sample=sample_id,
        size=_fmt_size(sv_path),
        variants=len(sv_variants),
    )
    return {"sv_variants": sv_variants, "sv_categories": sv_categories}


def load_sample_mito(sample_id: str, version: str | None = None) -> dict | None:
    """Staged loader: just the Mitochondria side-channel for a sample.
    {mito_variants, mito_categories} or None."""
    started = time.perf_counter()
    ctx = _load_pheno_context(sample_id, version)
    from ..adapters.mito_tsv import load_mito_tsv, MITO_TIERS
    if ctx is None:
        return None
    sub, _sd, _h, _p, pheno_by_gene = ctx
    mito_path = sample_layout.mito_tsv(sample_id)
    mv, mc = (
        load_mito_tsv(mito_path, pheno_by_gene=pheno_by_gene)
        if mito_path.exists() else ({}, {t: [] for t in MITO_TIERS})
    )
    _log_perf(
        "sample.mito",
        started,
        sample=sample_id,
        size=_fmt_size(mito_path),
        variants=len(mv),
    )
    return {"mito_variants": mv, "mito_categories": mc}


def load_sample_str(sample_id: str, version: str | None = None) -> dict | None:
    """Staged loader: STR side-channel for a sample."""
    started = time.perf_counter()
    sub = sample_layout.state_dir(sample_id)
    if not sub.is_dir():
        return None
    from ..adapters.str_tsv import STR_TIERS, load_str_tsv
    str_path = sample_layout.str_tsv(sample_id)
    variants, categories = (
        load_str_tsv(str_path)
        if str_path.exists() else ({}, {t: [] for t in STR_TIERS})
    )
    _log_perf(
        "sample.str",
        started,
        sample=sample_id,
        size=_fmt_size(str_path),
        variants=len(variants),
    )
    return {
        "str_variants": variants,
        "str_categories": categories,
        "str_pending": False,
    }


def load_sample_pgx(sample_id: str, version: str | None = None) -> dict | None:
    """Staged loader: PGx / PharmCAT side-channel for a sample."""
    started = time.perf_counter()
    sub = sample_layout.state_dir(sample_id)
    if not sub.is_dir():
        return None
    from ..adapters.pgx_tsv import load_pgx
    pgx_path = sample_layout.pgx_tsv(sample_id)
    pharmcat_path = sample_layout.pharmcat_json(sample_id)
    pgx = load_pgx(pgx_path, pharmcat_path if pharmcat_path.exists() else None)
    # Lazily import the DOCX presenter to avoid a module-load cycle:
    # docx_export itself uses this staged loader when building health reports.
    # The resulting view is the single shared projection for DOCX + browser.
    from .docx_export import build_pgx_report_view
    pgx["report_view"] = build_pgx_report_view(pgx)
    _log_perf(
        "sample.pgx",
        started,
        sample=sample_id,
        pgx_size=_fmt_size(pgx_path),
        pharmcat_size=_fmt_size(pharmcat_path),
        genes=len(pgx.get("gene_order") or []),
    )
    return {
        "pgx": pgx,
        "pharmcat": pgx,
        "pgx_pending": False,
    }


def _sample_snv_sidecar_context(sample_id: str, version: str | None = None):
    sub = sample_layout.state_dir(sample_id)
    if not sub.is_dir():
        return None
    meta = _read_json_or(sample_layout.state_file(sample_id, "sample_metadata.json"), {}) or {}
    chosen_version = _resolve_version(
        sample_id,
        requested=version,
        meta_active=meta.get("active_analysis"),
    )
    sidecar_dir = (
        analyses_store.version_dir(sample_id, chosen_version)
        if chosen_version is not None else sub
    )
    return sub, meta, sidecar_dir, _effective_test_type(meta, sample_id, default="WES")


def load_sample_secondary_snv(
    sample_id: str,
    version: str | None = None,
    *,
    clinvar_baseline: bool = False,
) -> dict | None:
    """Staged loader for secondary-finding SNV panel categories.

    ``clinvar_baseline`` is reserved for fixed-release DOCX export. The API/UI
    default continues to classify candidates from the weekly comparison.
    """
    started = time.perf_counter()
    ctx = _sample_snv_sidecar_context(sample_id, version)
    if ctx is None:
        return None
    sub, meta, sidecar_dir, test_type = ctx
    raw_snv_tsv = sample_layout.snv_raw_tsv(sample_id)
    review_tsv = sample_layout.review_tsv(sample_id)
    try:
        snv_tsv = snv_review.ensure_review_tsv(
            raw_snv_tsv,
            keep_ids=set(),
            test_type=test_type,
            output_dir=sub,
            output_path=sample_layout.review_tsv(sample_id, for_write=True),
            manifest_path=sample_layout.review_manifest(sample_id, for_write=True),
            overlay_path=sample_layout.snv_overlay_path(sample_id),
        )
    except OSError:
        if not review_tsv.is_file():
            categories = {category: [] for category in SECONDARY_SNV_PANELS}
            _log_perf(
                "sample.secondary_snv",
                started,
                sample=sample_id,
                selected="missing_review",
                raw_size=_fmt_size(raw_snv_tsv),
                review_size="missing",
                variants=0,
                **{category: 0 for category in SECONDARY_SNV_PANELS},
            )
            return {
                "variants": {},
                "categories": categories,
                "pheno_scores": {},
                "secondary_pending": False,
            }
        snv_tsv = review_tsv
    all_variants, tiers, _pheno, _old_format_error = _load_enriched_snv_cached(
        snv_tsv,
        sample_id=sample_id,
        sidecar_dir=sidecar_dir,
        test_type=test_type,
    )
    all_variants = {
        variant_id: dict(variant) for variant_id, variant in all_variants.items()
    }
    tiers = {tier: list(ids) for tier, ids in tiers.items()}
    _apply_effective_acmg(
        all_variants,
        tiers,
        sample_id=sample_id,
        genome_build=meta.get("genome_build") or "hg38",
        meta=meta,
    )
    if clinvar_baseline:
        all_variants = _pipeline_clinvar_secondary_variants(all_variants)
    categories = _build_secondary_snv_categories(all_variants)
    wanted_ids = {vid for ids in categories.values() for vid in ids}
    variants = {vid: all_variants[vid] for vid in wanted_ids if vid in all_variants}
    _log_perf(
        "sample.secondary_snv",
        started,
        sample=sample_id,
        selected=snv_tsv.name,
        raw_size=_fmt_size(raw_snv_tsv),
        review_size=_fmt_size(snv_tsv),
        variants=len(variants),
        **{category: len(categories.get(category) or []) for category in SECONDARY_SNV_PANELS},
    )
    return {
        "variants": variants,
        "categories": categories,
        "secondary_pending": False,
    }


def load_sample(sample_id: str, version: str | None = None,
                include_aux: bool = True) -> dict | None:
    """Build the per-sample webdata payload the frontend renders.

    `version` selects which analysis sidecar set to join in. When None,
    falls back to the sample's `active_analysis`, then `default`, then
    the legacy flat layout (sidecars at the sample root).

    `include_aux=False` skips the CNV/SV and Mito TSV parsing (returns
    empty dicts for them + `aux_pending: True`); the frontend then
    pulls those from /samples/{id}/cnv, /samples/{id}/sv, and /samples/{id}/mito so
    the SNV/Indel view appears without waiting on the rest.
    """
    started = time.perf_counter()
    sub = sample_layout.state_dir(sample_id)
    if not sub.is_dir():
        return None

    raw_snv_tsv = sample_layout.snv_raw_tsv(sample_id)
    _meta_early = _read_json_or(
        sample_layout.state_file(sample_id, "sample_metadata.json"), {}
    ) or {}
    _meta_early["test_type"] = _effective_test_type(_meta_early, sample_id, default="WES")
    reported_ids = set((_meta_early.get("status") or {}).keys())
    snv_tsv = raw_snv_tsv
    if raw_snv_tsv.exists():
        try:
            review_started = time.perf_counter()
            snv_tsv = snv_review.ensure_review_tsv(
                raw_snv_tsv,
                keep_ids=reported_ids,
                test_type=_meta_early["test_type"],
                output_dir=sub,
                output_path=sample_layout.review_tsv(sample_id, for_write=True),
                manifest_path=sample_layout.review_manifest(sample_id, for_write=True),
                overlay_path=sample_layout.snv_overlay_path(sample_id),
            )
            _log_perf(
                "sample.review_tsv",
                review_started,
                sample=sample_id,
                selected=snv_tsv.name,
                raw_size=_fmt_size(raw_snv_tsv),
                review_size=_fmt_size(snv_tsv),
            )
        except OSError:
            # Read-only / temporarily full disks should not break the
            # reviewer page. Prefer the last compact derived TSV (which
            # already contains overlay fields), then fall back to raw.
            existing_review = sample_layout.review_tsv(sample_id)
            snv_tsv = existing_review if existing_review.is_file() else raw_snv_tsv
            _log_perf(
                "sample.review_tsv",
                review_started,
                sample=sample_id,
                selected=snv_tsv.name,
                error="fallback_raw",
                raw_size=_fmt_size(raw_snv_tsv),
            )
    # Read meta early so the WES/WGS depth gate in load_snv_tsv kicks
    # in correctly. meta is re-read below for the response payload —
    # this duplicate is cheap (one small JSON) and keeps the gating
    # logic out of the adapter's API.
    _test_type = _meta_early["test_type"]
    meta = _read_json_or(sample_layout.state_file(sample_id, "sample_metadata.json"), {}) or {}
    meta["test_type"] = _effective_test_type(meta, sample_id, default="WES")

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

    pheno_path = _sidecar_file(sample_id, sidecar_dir, "pheno_score.tsv")

    # Lazy backfill: legacy samples + any pheno_score.tsv predating its
    # analysis.json (e.g. HPO/panels touched by a tool that bypassed
    # write_version) get recomputed inline so the Clinical/in-panel
    # consumers downstream always see a fresh table.
    analysis_path = _sidecar_file(sample_id, sidecar_dir, "analysis.json")
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

    variants, categories, pheno_by_gene, old_format_error = (
        _load_enriched_snv_cached(
            snv_tsv,
            sample_id=sample_id,
            sidecar_dir=sidecar_dir,
            test_type=_test_type,
        )
    )
    # The cache stores shared read-only maps. Per-load supplements
    # and ACMG overlays depend on reviewer/global state, so copy before edits.
    variants = {variant_id: dict(variant) for variant_id, variant in variants.items()}
    categories = {tier: list(ids) for tier, ids in categories.items()}
    if raw_snv_tsv.exists() and snv_tsv.name.endswith(snv_review.REVIEW_TSV_NAME):
        _supplement_marked_snv_variants(
            variants,
            categories,
            raw_snv_tsv,
            sample_id=sample_id,
            keep_ids=reported_ids,
            sidecar_dir=sidecar_dir,
            test_type=_test_type,
            index_path=sample_layout.snv_gene_index_path(sample_id),
            overlay_path=sample_layout.snv_overlay_path(sample_id),
        )
    _apply_effective_acmg(
        variants,
        categories,
        sample_id=sample_id,
        genome_build=meta.get("genome_build") or "hg38",
        meta=meta,
    )
    review_variants_for_secondary = variants

    # CNV / SV: load only when the AnnotSV outputs are present beside
    # the SNV TSV (pipeline drops them per-sample). Empty dicts when
    # absent → frontend just shows "（無資料）" placeholders. When
    # include_aux is False these are deferred to /samples/{id}/cnv,
    # /samples/{id}/sv, and /samples/{id}/mito (staged loading).
    from ..adapters.annotsv_tsv import load_annotsv_tsv, CNV_TIERS, SV_TIERS
    from ..adapters.mito_tsv import load_mito_tsv, MITO_TIERS
    if include_aux:
        cnv_path = sample_layout.cnv_tsv(sample_id)
        sv_path  = sample_layout.sv_tsv(sample_id)
        # CNV/SV gene tables render the pheno column as `matched/total`
        # (e.g. "2/3" — gene was implicated by 2 of 3 input HPO/panel
        # weights). Recompute the raw matched-weight + total-weight pair
        # here from the active analysis version's HPO/panels so the
        # numerator on each gene matches what compute_pheno_score
        # multiplied by 100 to write pheno_score.tsv.
        from . import phenotype_scorer
        pheno_matched, pheno_total = phenotype_scorer.compute_pheno_match(
            hpo_list, panels_list
        )
        cnv_variants, cnv_categories = (
            load_annotsv_tsv(
                cnv_path, source="cnv",
                pheno_by_gene=pheno_by_gene,
                pheno_matched=pheno_matched, pheno_total=pheno_total,
            )
            if cnv_path.exists() else ({}, {t: [] for t in CNV_TIERS})
        )
        sv_variants, sv_categories = (
            load_annotsv_tsv(
                sv_path, source="sv",
                pheno_by_gene=pheno_by_gene,
                pheno_matched=pheno_matched, pheno_total=pheno_total,
            )
            if sv_path.exists() else ({}, {t: [] for t in SV_TIERS})
        )
        # Mitochondria: per-sample mito.annotated.tsv (VEP + local MITOMAP),
        # produced by scripts/annotate_mito_vcf.sh. Pheno gene set drives the
        # "Clinical (in panel)" tier just like the SNV/CNV/SV sides.
        mito_path = sample_layout.mito_tsv(sample_id)
        mito_variants, mito_categories = (
            load_mito_tsv(mito_path, pheno_by_gene=pheno_by_gene)
            if mito_path.exists() else ({}, {t: [] for t in MITO_TIERS})
        )
    else:
        cnv_variants, cnv_categories = {}, {t: [] for t in CNV_TIERS}
        sv_variants,  sv_categories  = {}, {t: [] for t in SV_TIERS}
        mito_variants, mito_categories = {}, {t: [] for t in MITO_TIERS}

    qc = _read_json_or(sample_layout.state_file(sample_id, "qc_summary.json"), {}) or {}
    roh = _read_json_or(sample_layout.state_file(sample_id, "roh_summary.json"), {}) or {}
    ploidy_result = ploidy.load_sample_ploidy(sample_id)
    dead_zone_hits = panel_deadzone.dead_zone_for_genes(_test_type, set(pheno_by_gene.keys()))
    dead_zone_entries = []
    for gene, hit in dead_zone_hits.items():
        entry = dict(hit)
        entry["pheno_score"] = round(float(pheno_by_gene.get(gene) or 0), 2)
        dead_zone_entries.append(entry)
    secondary_categories = (
        _build_secondary_snv_categories(review_variants_for_secondary)
        if include_aux else {category: [] for category in SECONDARY_SNV_PANELS}
    )
    for category, ids in secondary_categories.items():
        categories[category] = ids

    annotation_metadata = annotation_versions.load_annotation_versions(
        raw_snv_tsv, sub
    )
    annotation_payload = {
        key: value
        for key, value in annotation_metadata.items()
        if key != "metadata_path"
    }
    clinvar_comparison = _read_json_or(
        sample_layout.clinvar_comparison_path(sample_id), {}
    ) or {}
    if clinvar_comparison.get("status") != "complete":
        clinvar_comparison = {}
    pipeline_clinvar_date = str(annotation_payload.get("clinvar_date") or "")
    latest_clinvar_date = str(clinvar_comparison.get("latest_release") or "")
    payload = {
        "meta": {
            "LIS_ID":         meta.get("lis_id") or meta.get("sample_id") or sample_id,
            "Name":           meta.get("name", ""),
            "MRN":            meta.get("mrn", ""),
            "Sex":            meta.get("sex", ""),
            "DOB":            meta.get("date_of_birth", ""),
            "Test":           _effective_test_type(meta, sample_id),
            "Category":       meta.get("category", ""),
            "Department":     meta.get("department", ""),
            "Physician":      meta.get("physician", ""),
            "SignReceivedAt": meta.get("sign_received_at", ""),
        },
        "genetic_counseling": meta.get("genetic_counseling", ""),
        "emr_synced_at":      meta.get("emr_synced_at", ""),
        "sample_id":         sample_id,
        "genome_build":      meta.get("genome_build", "hg38"),
        "generated_at":      meta.get("run_date") or datetime.utcnow().isoformat(timespec="seconds") + "Z",
        # Database release metadata comes from a compact pipeline sidecar;
        # never substitute today's date or a hard-coded report constant.
        "annotation_versions": annotation_payload,
        "clinvar_comparison": clinvar_comparison,
        "litvar2_lookup":     litvar2_on_demand.sample_status(sample_id),
        "patient_phenotype": _normalize_phenotype(hpo_list),
        "selected_panels":   panels_list,
        "vcf_path":          meta.get("vcf_path", ""),
        "qc_summary":        qc,
        "roh_summary":       roh,
        "ploidy":            ploidy_result,
        "dead_zone": {
            "threshold": panel_deadzone.dead_zone_threshold(_test_type),
            "entries": dead_zone_entries,
        },
        "variants":          variants,
        # When non-empty, the SNV/Indel TSV is in the pre-2026-05
        # layout and load_snv_tsv refused to parse it. Frontend uses
        # this to render a "請重跑新版 pipeline" banner instead of an
        # empty SNV card.
        "snv_tsv_error":     old_format_error,
        "categories":        categories,
        "tiers":             TIERS,
        # CNV / SV side-channels (independent variant maps + tier
        # categories so the frontend can render them in their own
        # tab group without colliding with SNV ids).
        "cnv_variants":      cnv_variants,
        "cnv_categories":    cnv_categories,
        "sv_variants":       sv_variants,
        "sv_categories":     sv_categories,
        "mito_variants":     mito_variants,
        "mito_categories":   mito_categories,
        "str_variants":      {},
        "str_categories":    {},
        "pgx":               {},
        # When True the CNV/SV + Mito side-channels above are empty
        # placeholders; the frontend fetches them from the dedicated
        # /samples/{id}/cnv, /samples/{id}/sv, and /samples/{id}/mito endpoints.
        "aux_pending":       not include_aux,
        "secondary_pending": not include_aux,
        # Whether the sample has phenotype configured at all — the
        # frontend uses this to show a "Clinical 區塊空白是因為沒有
        # phenotype" hint instead of leaving the panel silently empty.
        "has_phenotype":     bool(hpo_list or panels_list),
        "pharmcat":          {},
        # Active version metadata so the frontend can show a version
        # picker / detect when re-analysis should ask for a target.
        "active_analysis":   chosen_version,
        "analyses":          analyses_store.list_versions(sample_id),
    }
    payload["clinvar_pipeline_date"] = pipeline_clinvar_date
    payload["clinvar_latest_date"] = latest_clinvar_date
    # Every ordinary ClinVar label/report remains pinned to the Nextflow
    # baseline. The weekly date is exposed separately for change-arrow hover.
    payload["clinvar_date"] = pipeline_clinvar_date
    _log_perf(
        "sample.load",
        started,
        sample=sample_id,
        include_aux=include_aux,
        snv_tsv=snv_tsv.name,
        snv_variants=len(variants),
        cnv_variants=len(cnv_variants),
        sv_variants=len(sv_variants),
        mito_variants=len(mito_variants),
    )
    return payload


def search_snv_by_genes(
    sample_id: str,
    genes: list[str],
    *,
    version: str | None = None,
) -> dict | None:
    """Search the complete raw SNV TSV while the main UI uses review.tsv."""
    started = time.perf_counter()
    sub = sample_layout.state_dir(sample_id)
    if not sub.is_dir():
        return None
    raw_tsv = sample_layout.snv_raw_tsv(sample_id)
    if not raw_tsv.is_file():
        return {"variants": {}, "snv_tsv_error": ""}
    meta = _read_json_or(sample_layout.state_file(sample_id, "sample_metadata.json"), {}) or {}
    chosen_version = _resolve_version(
        sample_id, requested=version, meta_active=meta.get("active_analysis"),
    )
    sidecar_dir = (
        analyses_store.version_dir(sample_id, chosen_version)
        if chosen_version is not None else sub
    )
    wanted = {g.strip().upper() for g in genes if g.strip()}
    test_type = _effective_test_type(meta, sample_id, default="WES")
    indexed_rows = snv_gene_index.query_rows(
        raw_tsv, genes, sample_layout.snv_gene_index_path(sample_id)
    )
    old_format_error = ""
    if indexed_rows is not None:
        with snv_overlay.OverlayReader(
            raw_tsv, sample_layout.snv_overlay_path(sample_id)
        ) as overlay:
            indexed_rows = overlay.apply_many(indexed_rows)
        variants = _variants_from_rows(indexed_rows, test_type=test_type)
        _enrich_snv_variants(variants, sample_id, sidecar_dir)
        matches = variants
        search_source = "gene_index"
        raw_variant_count = "indexed"
    else:
        indexed_rows = _scan_raw_snv_rows_by_genes(raw_tsv, wanted)
        with snv_overlay.OverlayReader(
            raw_tsv, sample_layout.snv_overlay_path(sample_id)
        ) as overlay:
            indexed_rows = overlay.apply_many(indexed_rows)
        variants = _variants_from_rows(indexed_rows, test_type=test_type)
        _enrich_snv_variants(variants, sample_id, sidecar_dir)
        matches = variants
        search_source = "raw_stream_fallback"
        raw_variant_count = "streamed"
    _log_perf(
        "sample.snv_search",
        started,
        sample=sample_id,
        genes=len(wanted),
        source=search_source,
        raw_size=_fmt_size(raw_tsv),
        raw_variants=raw_variant_count,
        matches=len(matches),
    )
    return {
        "variants": matches,
        "snv_tsv_error": old_format_error,
        "litvar2_lookup": litvar2_on_demand.sample_status(sample_id),
    }


def warm_raw_snv_cache(sample_id: str, version: str | None = None) -> None:
    """Best-effort background prewarm for the complete-TSV gene search."""
    started = time.perf_counter()
    try:
        search_snv_by_genes(sample_id, [], version=version)
        _log_perf("sample.raw_snv_warm", started, sample=sample_id, status="ok")
    except Exception:
        # Registration already succeeded. Search falls back to a
        # foreground parse later if this opportunistic warm-up fails.
        _log_perf("sample.raw_snv_warm", started, sample=sample_id, status="failed")
        pass
