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
    _row_to_variant,
    load_snv_tsv,
    merge_snv_variant_row,
)
from ..config import TERTIARY_OUTPUT_ROOT
from . import analyses_store, omim_store, panel_deadzone, phenotype_scorer, snv_gene_index, snv_review


SECONDARY_SNV_PANELS = {
    "acmg_sf": "ACMG_SF_v3.3",
    "proactive": "proactive",
    "carrier": "carrier_mackenzie_1300+",
}

_SNV_CACHE_MAX = 8
_snv_cache: OrderedDict[tuple, tuple[dict, dict, dict, str]] = OrderedDict()
_snv_cache_lock = threading.Lock()
_CASE_SUMMARY_CACHE_MAX = 128
_case_summary_cache: OrderedDict[tuple, dict[str, str]] = OrderedDict()
_case_summary_cache_lock = threading.Lock()
CASE_SUMMARY_CACHE_NAME = "case_summary.json"


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
    sample_dir = Path(sample_dir)
    with _snv_cache_lock:
        for key in list(_snv_cache):
            if Path(key[0][0]).parent == sample_dir:
                del _snv_cache[key]
    with _case_summary_cache_lock:
        for key in list(_case_summary_cache):
            if Path(key[0][0]).parent == sample_dir:
                del _case_summary_cache[key]


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


def _file_signature(path: Path) -> tuple[str, int, int]:
    """Cheap cache signature; absent files are represented consistently."""
    try:
        st = path.stat()
        return (str(path), st.st_mtime_ns, st.st_size)
    except OSError:
        return (str(path), 0, 0)


def _case_summary_signature(sample_dir: Path, omim_sig: tuple | None = None) -> list:
    return [
        list(_file_signature(sample_dir / "sample_metadata.json")),
        list(_file_signature(sample_dir / "snv_indel.annotated.tsv")),
        list(_file_signature(sample_dir / "snv_gene_index.sqlite")),
        list(_file_signature(sample_dir / "snv_indel.review.tsv")),
        list(_file_signature(sample_dir / "cnv.annotated.tsv")),
        list(_file_signature(sample_dir / "sv.annotated.tsv")),
        list(omim_sig if omim_sig is not None else omim_store.cache_signature()),
    ]


def _case_summary_cache_path(sample_dir: Path) -> Path:
    return sample_dir / CASE_SUMMARY_CACHE_NAME


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
    acmg = (
        edits.get("ACMG_classification")
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

    def genes_for(key: str) -> list[str]:
        genes: list[str] = []
        for raw in phenotype_scorer._HPO_TO_GENES.get(key, set()):
            gene, hgnc_id = panel_deadzone.canonical_gene_symbol(raw)
            if panel_deadzone.is_disease_associated_gene(gene, hgnc_id):
                genes.append(gene)
        return sorted(set(genes))

    for row in hpo_rows or []:
        hid = row.get("phenotype") or row.get("hpo_id") or ""
        if not hid:
            continue
        label = row.get("label") or hid
        sections.append({"name": label, "genes": genes_for(hid)})
    for entry in panel_entries or []:
        name = entry.get("name") if isinstance(entry, dict) else str(entry)
        if not name:
            continue
        sections.append({"name": name, "genes": genes_for(name)})

    merged: set[str] = set()
    for section in sections:
        merged.update(section.get("genes") or [])
    return {"grouped": sections, "merged": sorted(merged)}


def load_report_gene_list(sample_id: str, version: str | None = None) -> dict | None:
    sub = TERTIARY_OUTPUT_ROOT / sample_id
    if not sub.is_dir():
        return None
    meta = _read_json_or(sub / "sample_metadata.json", {}) or {}
    if meta is None:
        return None
    ctx = _analysis_context(sample_id, meta, version)
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
    for idx in range(1, 6):
        if picked.get(str(idx)) or picked.get(idx):
            disease = omim_store.compact_disease_label(
                rec.get(f"Disease{idx}") or ""
            )
            if disease and disease not in out:
                out.append(disease)
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
    raw_tsv = sample_dir / "snv_indel.annotated.tsv"
    if not wanted or not raw_tsv.exists():
        return {}
    rows = snv_gene_index.query_rows_by_ids(raw_tsv, wanted)
    if rows is None:
        rows = _scan_snv_review_rows_by_ids(
            sample_dir / "snv_indel.review.tsv",
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
    for row in rows or []:
        try:
            variant = _row_to_variant(row)
        except (KeyError, TypeError, ValueError):
            continue
        if variant.get("id") in wanted:
            variants[variant["id"]] = variant
    if variants:
        _enrich_snv_variants(variants, sample_dir)
    return variants


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
        for source, filename in (("cnv", "cnv.annotated.tsv"), ("sv", "sv.annotated.tsv")):
            hits = load_annotsv_variants_by_ids(
                sample_dir / filename,
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
        _augment_merged_cnv_sv_acmg(
            sample_dir / ("cnv.annotated.tsv" if source == "cnv" else "sv.annotated.tsv"),
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
    snv_variants = _case_snv_variants_by_id(sample_dir, wanted)
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

    cnv_sv_wanted = wanted.difference(snv_variants.keys())
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
    sidecar_dir: Path,
    test_type: str,
) -> tuple[dict, dict, dict, str]:
    """Parse + enrich SNVs once per input revision, with a bounded LRU."""
    started = time.perf_counter()
    exo_path = sidecar_dir / "exomiser_results.tsv"
    lir_path = sidecar_dir / "lirical_results.tsv"
    pheno_path = sidecar_dir / "pheno_score.tsv"
    analysis_path = sidecar_dir / "analysis.json"
    omim_sig = omim_store.cache_signature()
    key = (
        _file_signature(snv_tsv),
        _file_signature(analysis_path),
        _file_signature(pheno_path),
        _file_signature(exo_path),
        _file_signature(lir_path),
        (test_type or "WES").upper(),
        omim_sig,
    )
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
    pheno_by_gene = _enrich_snv_variants(variants, sidecar_dir)
    for t, ids in categories.items():
        categories[t] = sorted(ids, key=lambda i: (-_variant_total_score(variants, i), i))

    result = (variants, categories, pheno_by_gene, old_format_error)
    with _snv_cache_lock:
        # Callers treat the enriched maps as read-only; keeping the same
        # objects avoids a full deep-copy cost on large cached payloads.
        _snv_cache[key] = result
        _snv_cache.move_to_end(key)
        while len(_snv_cache) > _SNV_CACHE_MAX:
            _snv_cache.popitem(last=False)
    _log_perf(
        "snv.enriched_cache",
        started,
        cache="miss",
        tsv=snv_tsv.name,
        size=_fmt_size(snv_tsv),
        variants=len(variants),
        parse=f"{parse_elapsed:.3f}s",
        enrich=f"{time.perf_counter() - enrich_started:.3f}s",
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


def _is_plp_text(raw: str) -> bool:
    key = str(raw or "").strip().lower().replace("_", " ")
    return key in {
        "pathogenic",
        "likely pathogenic",
        "pathogenic/likely pathogenic",
        "likely pathogenic/pathogenic",
    }


def _is_acmg_plp(variant: dict) -> bool:
    cls = (
        variant.get("genebe_acmg_class")
        or variant.get("ACMG_classification")
        or ""
    )
    return _is_plp_text(cls)


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
    return (
        (_is_clinvar_plp(variant) or _is_acmg_plp(variant))
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


def _enrich_snv_variants(variants: dict[str, dict], sidecar_dir: Path) -> dict[str, float]:
    """Join pheno / Exomiser / LIRICAL / OMIM data into variant payloads."""
    exo_path = sidecar_dir / "exomiser_results.tsv"
    lir_path = sidecar_dir / "lirical_results.tsv"
    pheno_path = sidecar_dir / "pheno_score.tsv"
    exo = _read_tsv_dict(exo_path)
    lir = _read_tsv_dict(lir_path)
    pheno_by_gene = _read_pheno_scores(pheno_path)

    def _scale_to_100(s):
        n = _to_num(s)
        if not isinstance(n, (int, float)):
            return None
        return int(round(n * 100))

    # Refresh OMIM once per batch; individual variant joins use only
    # in-memory dict lookups and no longer stat OMIM.xlsx repeatedly.
    omim_store.ensure_loaded()
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
            for f in ("Disease1", "Disease2", "Disease3", "Disease4", "Disease5"):
                v[f] = rec.get(f, "")
    return pheno_by_gene


def _variants_from_rows(rows: list[dict[str, str]], *, test_type: str) -> dict[str, dict]:
    """Shape raw TSV rows into variant payloads without parsing the whole TSV."""
    out: dict[str, dict] = {}
    is_wes = (test_type or "").upper() == "WES"
    low_flag = _DP_LOW_FLAG_WES if is_wes else _DP_LOW_FLAG_WGS
    for row in rows:
        canonical_gene, _ = panel_deadzone.canonical_gene_symbol(
            row.get("GENE", ""),
            row.get("HGNC_ID", ""),
        )
        if canonical_gene in EXCLUDED_GENES:
            continue
        v = _row_to_variant(row)
        dp = v.get("depth") or 0
        if is_wes and dp < _DP_HARD_FLOOR_WES:
            continue
        v["low_depth"] = bool(dp and dp < low_flag)
        merge_snv_variant_row(out, v)
    return out


def _supplement_marked_snv_variants(
    variants: dict[str, dict],
    categories: dict[str, list[str]],
    raw_tsv: Path,
    *,
    keep_ids: set[str],
    sidecar_dir: Path,
    test_type: str,
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
    rows = snv_gene_index.query_rows_by_ids(raw_tsv, missing)
    if rows is None:
        _log_perf(
            "sample.marked_snv_supplement",
            started,
            status="skipped_no_index",
            missing=len(missing),
        )
        return
    extra = _variants_from_rows(rows, test_type=test_type)
    _enrich_snv_variants(extra, sidecar_dir)
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


def list_case_summaries() -> list[dict]:
    """Return the case-list table rows, including cached variant summaries."""
    out: list[dict] = []
    if not TERTIARY_OUTPUT_ROOT.exists():
        return out
    omim_sig = omim_store.cache_signature()
    for sub in sorted(TERTIARY_OUTPUT_ROOT.iterdir()):
        if not sub.is_dir() or sub.name.startswith("_"):
            continue
        meta_path = sub / "sample_metadata.json"
        if not meta_path.exists():
            continue
        meta = _read_json_or(meta_path, {}) or {}
        summary = _case_management_summary(sub, meta, omim_sig=omim_sig)
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
            **summary,
        })
    out.sort(key=lambda r: (r.get("created_at") or "", r.get("lis_id") or ""), reverse=True)
    return out


def list_unregistered() -> list[dict]:
    """Sample dirs the tertiary pipeline left behind without a metadata file.

    A directory under tertiary_output/ counts as "unregistered" when it
    has snv_indel.annotated.tsv but NO sample_metadata.json yet. The UI
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
        if lis_id in active_tertiary_samples:
            continue

        source_info = {}
        source_vcf_path = ""
        source_vcf_size = 0
        pipeline_source = sub / "pipeline_source.json"
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
            mtime = sub.stat().st_mtime
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
                "test_type":        (roster_entry or {}).get("test_type", ""),
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


def _read_pheno_by_gene(sidecar_dir) -> dict:
    """Read pheno_score.tsv → {gene_symbol: 0-100 score}. Empty if absent."""
    out: dict[str, float] = {}
    p = sidecar_dir / "pheno_score.tsv"
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
    sub = TERTIARY_OUTPUT_ROOT / sample_id
    if not sub.is_dir():
        return None
    meta = _read_json_or(sub / "sample_metadata.json", {}) or {}
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
    return sub, sidecar_dir, hpo_list, panels_list, _read_pheno_by_gene(sidecar_dir)


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
    cnv_path = sub / "cnv.annotated.tsv"
    sv_path  = sub / "sv.annotated.tsv"
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
    cnv_path = sub / "cnv.annotated.tsv"
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
    sv_path = sub / "sv.annotated.tsv"
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
    mito_path = sub / "mito.annotated.tsv"
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
    sub = TERTIARY_OUTPUT_ROOT / sample_id
    if not sub.is_dir():
        return None
    from ..adapters.str_tsv import STR_TIERS, load_str_tsv
    str_path = sub / "str.tsv"
    if not str_path.exists():
        str_path = sub / "str.annotated.tsv"
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
    sub = TERTIARY_OUTPUT_ROOT / sample_id
    if not sub.is_dir():
        return None
    from ..adapters.pgx_tsv import load_pgx
    pgx_path = sub / "pgx.tsv"
    if not pgx_path.exists():
        hits = sorted(sub.glob("*.pgx.tsv"))
        pgx_path = hits[0] if hits else pgx_path
    pharmcat_path = sub / "pharmcat.report.json"
    if not pharmcat_path.exists():
        hits = sorted(sub.glob("*pharmcat*.json"))
        pharmcat_path = hits[0] if hits else pharmcat_path
    pgx = load_pgx(pgx_path, pharmcat_path if pharmcat_path.exists() else None)
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
    sub = TERTIARY_OUTPUT_ROOT / sample_id
    if not sub.is_dir():
        return None
    meta = _read_json_or(sub / "sample_metadata.json", {}) or {}
    chosen_version = _resolve_version(
        sample_id,
        requested=version,
        meta_active=meta.get("active_analysis"),
    )
    sidecar_dir = (
        analyses_store.version_dir(sample_id, chosen_version)
        if chosen_version is not None else sub
    )
    return sub, meta, sidecar_dir, (meta.get("test_type") or "WES").upper()


def load_sample_secondary_snv(sample_id: str, version: str | None = None) -> dict | None:
    """Staged loader: ACMG SF / Proactive / Carrier SNV side-channel."""
    started = time.perf_counter()
    ctx = _sample_snv_sidecar_context(sample_id, version)
    if ctx is None:
        return None
    sub, _meta, sidecar_dir, test_type = ctx
    raw_snv_tsv = sub / "snv_indel.annotated.tsv"
    review_tsv = raw_snv_tsv.with_name(snv_review.REVIEW_TSV_NAME)
    try:
        snv_tsv = snv_review.ensure_review_tsv(
            raw_snv_tsv,
            keep_ids=set(),
            test_type=test_type,
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
                acmg_sf=0,
                proactive=0,
                carrier=0,
            )
            return {
                "variants": {},
                "categories": categories,
                "pheno_scores": {},
                "secondary_pending": False,
            }
        snv_tsv = review_tsv
    all_variants, _tiers, _pheno, _old_format_error = _load_enriched_snv_cached(
        snv_tsv,
        sidecar_dir=sidecar_dir,
        test_type=test_type,
    )
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
        acmg_sf=len(categories.get("acmg_sf") or []),
        proactive=len(categories.get("proactive") or []),
        carrier=len(categories.get("carrier") or []),
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
    sub = TERTIARY_OUTPUT_ROOT / sample_id
    if not sub.is_dir():
        return None

    raw_snv_tsv = sub / "snv_indel.annotated.tsv"
    _meta_early = _read_json_or(sub / "sample_metadata.json", {}) or {}
    reported_ids = set((_meta_early.get("status") or {}).keys())
    snv_tsv = raw_snv_tsv
    if raw_snv_tsv.exists():
        try:
            review_started = time.perf_counter()
            snv_tsv = snv_review.ensure_review_tsv(
                raw_snv_tsv,
                keep_ids=reported_ids,
                test_type=(_meta_early.get("test_type") or "WES").upper(),
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
            # reviewer page. Fall back to the complete source TSV.
            snv_tsv = raw_snv_tsv
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
    _test_type  = (_meta_early.get("test_type") or "WES").upper()
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

    variants, categories, pheno_by_gene, old_format_error = (
        _load_enriched_snv_cached(
            snv_tsv, sidecar_dir=sidecar_dir, test_type=_test_type,
        )
    )
    review_variants_for_secondary = variants
    # The cache stores shared read-only maps. Per-load supplements
    # depend on reviewer status, so copy before adding marked variants.
    variants = dict(variants)
    categories = {tier: list(ids) for tier, ids in categories.items()}
    if raw_snv_tsv.exists() and snv_tsv.name == snv_review.REVIEW_TSV_NAME:
        _supplement_marked_snv_variants(
            variants,
            categories,
            raw_snv_tsv,
            keep_ids=reported_ids,
            sidecar_dir=sidecar_dir,
            test_type=_test_type,
        )

    # CNV / SV: load only when the AnnotSV outputs are present beside
    # the SNV TSV (pipeline drops them per-sample). Empty dicts when
    # absent → frontend just shows "（無資料）" placeholders. When
    # include_aux is False these are deferred to /samples/{id}/cnv,
    # /samples/{id}/sv, and /samples/{id}/mito (staged loading).
    from ..adapters.annotsv_tsv import load_annotsv_tsv, CNV_TIERS, SV_TIERS
    from ..adapters.mito_tsv import load_mito_tsv, MITO_TIERS
    if include_aux:
        cnv_path = sub / "cnv.annotated.tsv"
        sv_path  = sub / "sv.annotated.tsv"
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
        mito_path = sub / "mito.annotated.tsv"
        mito_variants, mito_categories = (
            load_mito_tsv(mito_path, pheno_by_gene=pheno_by_gene)
            if mito_path.exists() else ({}, {t: [] for t in MITO_TIERS})
        )
    else:
        cnv_variants, cnv_categories = {}, {t: [] for t in CNV_TIERS}
        sv_variants,  sv_categories  = {}, {t: [] for t in SV_TIERS}
        mito_variants, mito_categories = {}, {t: [] for t in MITO_TIERS}

    qc  = _read_json_or(sub / "qc_summary.json",  {}) or {}
    roh = _read_json_or(sub / "roh_summary.json", {}) or {}
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

    payload = {
        "meta": {
            "LIS_ID":         meta.get("lis_id") or meta.get("sample_id") or sample_id,
            "Name":           meta.get("name", ""),
            "MRN":            meta.get("mrn", ""),
            "Sex":            meta.get("sex", ""),
            "DOB":            meta.get("date_of_birth", ""),
            "Test":           meta.get("test_type", ""),
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
        "patient_phenotype": _normalize_phenotype(hpo_list),
        "selected_panels":   panels_list,
        "vcf_path":          meta.get("vcf_path", ""),
        "qc_summary":        qc,
        "roh_summary":       roh,
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
    sub = TERTIARY_OUTPUT_ROOT / sample_id
    if not sub.is_dir():
        return None
    raw_tsv = sub / "snv_indel.annotated.tsv"
    if not raw_tsv.is_file():
        return {"variants": {}, "snv_tsv_error": ""}
    meta = _read_json_or(sub / "sample_metadata.json", {}) or {}
    chosen_version = _resolve_version(
        sample_id, requested=version, meta_active=meta.get("active_analysis"),
    )
    sidecar_dir = (
        analyses_store.version_dir(sample_id, chosen_version)
        if chosen_version is not None else sub
    )
    wanted = {g.strip().upper() for g in genes if g.strip()}
    test_type = (meta.get("test_type") or "WES").upper()
    indexed_rows = snv_gene_index.query_rows(raw_tsv, genes)
    old_format_error = ""
    if indexed_rows is not None:
        variants = _variants_from_rows(indexed_rows, test_type=test_type)
        _enrich_snv_variants(variants, sidecar_dir)
        matches = variants
        search_source = "gene_index"
        raw_variant_count = "indexed"
    else:
        variants, _categories, _pheno, old_format_error = _load_enriched_snv_cached(
            raw_tsv,
            sidecar_dir=sidecar_dir,
            test_type=test_type,
        )
        matches = {
            vid: v for vid, v in variants.items()
            if (v.get("gene_symbol") or "").upper() in wanted
        }
        search_source = "raw_tsv_fallback"
        raw_variant_count = len(variants)
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
    return {"variants": matches, "snv_tsv_error": old_format_error}


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
