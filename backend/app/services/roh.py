"""ROH source selection, normalization, loading, and SNV interval joins.

All outputs owned by this module live in ``08_postprocessing``.  The upstream
Nextflow 00-07 tree is treated as read-only input.
"""
from __future__ import annotations

import csv
import json
import os
import shutil
import tempfile
from bisect import bisect_right
from pathlib import Path
from typing import Callable

from . import sample_layout


REGIONS_NAME = "roh_regions.tsv"
SUMMARY_NAME = "roh_summary.json"
AUTOMAP_TSV_NAME = "roh.source.automap.tsv"
AUTOMAP_PDF_NAME = "roh.source.automap.pdf"
BCFTOOLS_NAME = "roh.source.bcftools.txt"
DRAGEN_BED_NAME = "roh.source.dragen.bed"
DRAGEN_METRICS_NAME = "roh.source.dragen_metrics.csv"

MANAGED_NAMES = {
    REGIONS_NAME,
    SUMMARY_NAME,
    AUTOMAP_TSV_NAME,
    AUTOMAP_PDF_NAME,
    BCFTOOLS_NAME,
    DRAGEN_BED_NAME,
    DRAGEN_METRICS_NAME,
}

REGION_COLUMNS = (
    "region_id", "chrom", "start0", "end0", "display_start", "display_end",
    "length_bp", "length_mb", "source", "passes_default_filter",
    "n_markers", "n_hom", "n_het", "homozygosity_pct", "quality", "score",
)


def _chrom(value: object) -> str:
    raw = str(value or "").strip()
    if raw.lower().startswith("chr"):
        raw = raw[3:]
    raw = {"23": "X", "24": "Y"}.get(raw.upper(), raw.upper())
    return f"chr{raw}" if raw in {str(i) for i in range(1, 23)} | {"X", "Y"} else ""


def _number(value: object, kind=float):
    try:
        return kind(str(value).strip())
    except (TypeError, ValueError):
        return None


def _region(chrom: object, start0: int, end0: int, source: str, **fields) -> dict | None:
    canonical = _chrom(chrom)
    if not canonical or start0 < 0 or end0 <= start0:
        return None
    length = end0 - start0
    row = {
        "chrom": canonical,
        "start0": start0,
        "end0": end0,
        "display_start": start0 + 1,
        "display_end": end0,
        "length_bp": length,
        "length_mb": round(length / 1_000_000, 6),
        "source": source,
        "passes_default_filter": True,
        "n_markers": "",
        "n_hom": "",
        "n_het": "",
        "homozygosity_pct": "",
        "quality": "",
        "score": "",
    }
    row.update(fields)
    return row


def parse_automap(path: Path) -> tuple[list[dict], dict]:
    regions: list[dict] = []
    version = ""
    parameters: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if "AutoMap" in line and not version:
            version = line.lstrip("# ")
        if line.startswith("#"):
            for token in line.lstrip("# ").split():
                if "=" in token:
                    key, value = token.split("=", 1)
                    parameters[key] = value.rstrip(",")
            continue
        parts = line.split("\t")
        if len(parts) < 6 or parts[0].lower().lstrip("#") in {"chr", "chrom"}:
            continue
        begin, end = _number(parts[1], int), _number(parts[2], int)
        if begin is None or end is None:
            continue
        item = _region(
            parts[0], begin - 1, end, "automap",
            n_markers=_number(parts[4], int) or "",
            homozygosity_pct=_number(parts[5], float) or "",
        )
        if item:
            regions.append(item)
    return regions, {"tool_version": version or "AutoMap", "parameters": parameters}


def parse_bcftools(path: Path) -> tuple[list[dict], dict]:
    regions: list[dict] = []
    version = ""
    command = ""
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            if "bcftools" in line.lower() and not version:
                version = line.lstrip("# ")
            if "command" in line.lower() and "bcftools" in line.lower():
                command = line.lstrip("# ")
            continue
        parts = line.split("\t")
        if len(parts) < 8 or parts[0] != "RG":
            continue
        start, end = _number(parts[3], int), _number(parts[4], int)
        markers, quality = _number(parts[6], int), _number(parts[7], float)
        if start is None or end is None:
            continue
        item = _region(
            parts[2], start - 1, end, "bcftools",
            n_markers=markers if markers is not None else "",
            quality=quality if quality is not None else "",
        )
        if item:
            item["passes_default_filter"] = bool(
                item["length_bp"] >= 1_000_000
                and (markers or 0) >= 25
                and (quality or 0) >= 20
            )
            regions.append(item)
    return regions, {"tool_version": version or "BCFtools roh", "command": command}


def parse_dragen(path: Path) -> tuple[list[dict], dict]:
    regions: list[dict] = []
    version = "DRAGEN"
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                parts = line.split()
            if len(parts) < 6:
                continue
            start0, end0 = _number(parts[1], int), _number(parts[2], int)
            if start0 is None or end0 is None:
                continue
            item = _region(
                parts[0], start0, end0, "dragen",
                score=_number(parts[3], float) if _number(parts[3], float) is not None else "",
                n_hom=_number(parts[4], int) if _number(parts[4], int) is not None else "",
                n_het=_number(parts[5], int) if _number(parts[5], int) is not None else "",
            )
            if item:
                item["passes_default_filter"] = item["length_bp"] >= 3_000_000
                regions.append(item)
    return regions, {"tool_version": version}


def parse_dragen_metrics(path: Path | None) -> dict:
    if path is None or not path.is_file():
        return {}
    metrics: dict[str, float | int | str] = {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.reader(handle):
            if len(row) < 2:
                continue
            key = (row[2] if len(row) >= 4 and row[2].strip() else row[0]).strip()
            value = row[-1].strip()
            parsed = _number(value, float)
            metrics[key] = parsed if parsed is not None else value
    return metrics


def _merged_total_bp(regions: list[dict]) -> int:
    intervals: dict[str, list[tuple[int, int]]] = {}
    for row in regions:
        intervals.setdefault(row["chrom"], []).append((row["start0"], row["end0"]))
    total = 0
    for spans in intervals.values():
        current_start = current_end = None
        for start, end in sorted(spans):
            if current_end is None or start > current_end:
                if current_end is not None:
                    total += current_end - current_start
                current_start, current_end = start, end
            else:
                current_end = max(current_end, end)
        if current_end is not None:
            total += current_end - current_start
    return total


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        os.replace(tmp_name, path)
    finally:
        Path(tmp_name).unlink(missing_ok=True)


def _atomic_regions(path: Path, regions: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=REGION_COLUMNS)
            writer.writeheader()
            for index, row in enumerate(
                sorted(regions, key=lambda r: (int(r["chrom"][3:]) if r["chrom"][3:].isdigit() else 23, r["start0"])),
                1,
            ):
                record = dict(row)
                record["region_id"] = f"roh-{index}"
                record["passes_default_filter"] = "1" if record["passes_default_filter"] else "0"
                writer.writerow({key: record.get(key, "") for key in REGION_COLUMNS})
        os.replace(tmp_name, path)
    finally:
        Path(tmp_name).unlink(missing_ok=True)


def _copy(source: Path | None, destination: Path) -> bool:
    if source is None or not source.is_file():
        return False
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    return True


def _inhouse_root(source_vcf: Path) -> Path:
    return source_vcf.parent.parent if source_vcf.parent.name == "04_snv_indel" else source_vcf.parent


def _dragen_root(source_vcf: Path) -> Path:
    return source_vcf.parent.parent if source_vcf.parent.name == "vcf.gz" else source_vcf.parent


def _scoped(post_dir: Path, sample_id: str, name: str) -> Path:
    return sample_layout.scoped_file(
        post_dir, sample_id, name, for_write=True, force_prefixed=True,
    )


def prepare_roh_outputs(
    *,
    mode: str,
    source_vcf: str | Path,
    source_sample_id: str,
    sample_id: str,
    post_dir: Path,
    research_only: bool,
    automap_runner: Callable[[Path, str], tuple[Path | None, Path | None]] | None = None,
) -> dict:
    """Select and normalize one ROH source into a staged post directory."""
    source_vcf = Path(source_vcf)
    warnings: list[str] = []
    generated_automap = False
    source_path: Path | None = None
    aux_path: Path | None = None
    source = "missing"
    reason = ""
    parser_meta: dict = {}
    regions: list[dict] = []

    if mode == "dragen":
        root = _dragen_root(source_vcf)
        germline = root / "other" / source_sample_id / "germline_seq"
        bed = germline / f"{source_sample_id}.roh.bed"
        metrics = germline / f"{source_sample_id}.roh_metrics.csv"
        if bed.is_file():
            source, source_path, aux_path = "dragen", bed, metrics if metrics.is_file() else None
            regions, parser_meta = parse_dragen(bed)
            parser_meta["metrics"] = parse_dragen_metrics(aux_path)
            reason = "DRAGEN case：固定使用原生 DRAGEN ROH"
        else:
            reason = "找不到同一 DRAGEN sample 的原生 roh.bed"
    elif mode == "inhouse":
        root = _inhouse_root(source_vcf)
        roh_dir = root / "08_roh"
        automap_tsv = roh_dir / f"{source_sample_id}.HomRegions.tsv"
        automap_pdf = roh_dir / f"{source_sample_id}.HomRegions.pdf"
        bcftools = roh_dir / f"{source_sample_id}.roh.txt"
        if not automap_tsv.is_file() and research_only:
            hc_vcf = root / "04_snv_indel" / f"{source_sample_id}.haplotypecaller.vcf.gz"
            if automap_runner is None:
                warnings.append("Research-only 已勾選，但 AutoMap runner 未設定；改用 BCFtools")
            elif not hc_vcf.is_file():
                warnings.append("Research-only 已勾選，但找不到 HaplotypeCaller VCF；改用 BCFtools")
            else:
                try:
                    generated_tsv, generated_pdf = automap_runner(hc_vcf, source_sample_id)
                    if generated_tsv and generated_tsv.is_file():
                        automap_tsv, automap_pdf = generated_tsv, generated_pdf or Path("")
                        generated_automap = True
                    else:
                        warnings.append("AutoMap 未產生 HomRegions.tsv；改用 BCFtools")
                except Exception as exc:
                    warnings.append(f"AutoMap 執行失敗（{exc}）；改用 BCFtools")
        if automap_tsv.is_file():
            source, source_path = "automap", automap_tsv
            aux_path = automap_pdf if automap_pdf.is_file() else None
            regions, parser_meta = parse_automap(automap_tsv)
            reason = "使用既有 AutoMap" if not generated_automap else "Research-only：三級 postprocessing 補跑 AutoMap"
        elif bcftools.is_file():
            source, source_path = "bcftools", bcftools
            regions, parser_meta = parse_bcftools(bcftools)
            reason = "缺少 AutoMap；使用二級 BCFtools roh"
        else:
            reason = "找不到 AutoMap HomRegions.tsv 或 BCFtools roh.txt"
    else:
        raise ValueError(f"unsupported ROH mode: {mode}")

    artifact_names: dict[str, str] = {}
    if source == "automap":
        dst = _scoped(post_dir, sample_id, AUTOMAP_TSV_NAME)
        if _copy(source_path, dst):
            artifact_names["automap_tsv"] = dst.name
        pdf_dst = _scoped(post_dir, sample_id, AUTOMAP_PDF_NAME)
        if _copy(aux_path, pdf_dst):
            artifact_names["automap_pdf"] = pdf_dst.name
        root = _inhouse_root(source_vcf)
        bcf = root / "08_roh" / f"{source_sample_id}.roh.txt"
        bcf_dst = _scoped(post_dir, sample_id, BCFTOOLS_NAME)
        if _copy(bcf, bcf_dst):
            artifact_names["bcftools_txt"] = bcf_dst.name
    elif source == "bcftools":
        dst = _scoped(post_dir, sample_id, BCFTOOLS_NAME)
        if _copy(source_path, dst):
            artifact_names["bcftools_txt"] = dst.name
    elif source == "dragen":
        dst = _scoped(post_dir, sample_id, DRAGEN_BED_NAME)
        if _copy(source_path, dst):
            artifact_names["dragen_bed"] = dst.name
        metrics_dst = _scoped(post_dir, sample_id, DRAGEN_METRICS_NAME)
        if _copy(aux_path, metrics_dst):
            artifact_names["dragen_metrics"] = metrics_dst.name

    selected = [row for row in regions if row["passes_default_filter"]]
    autosomal = [row for row in selected if row["chrom"][3:].isdigit()]
    all_autosomal = [row for row in regions if row["chrom"][3:].isdigit()]
    default_filter = {
        "automap": {"min_length_mb": 1.0, "description": "AutoMap 原生門檻（≥1 Mb）"},
        "bcftools": {"min_length_mb": 1.0, "min_markers": 25, "min_quality": 20, "description": "≥1 Mb、markers ≥25、quality ≥20"},
        "dragen": {"min_length_mb": 3.0, "description": "DRAGEN large ROH（≥3 Mb）"},
    }.get(source, {})
    status = "source_missing" if source == "missing" else ("complete" if regions else "complete_no_regions")
    summary = {
        "schema_version": 1,
        "status": status,
        "source": source,
        "source_label": {"automap": "AutoMap", "bcftools": "BCFtools roh", "dragen": "DRAGEN ROH"}.get(source, "No ROH source"),
        "selection_reason": reason,
        "pipeline_type": mode,
        "research_only": bool(research_only),
        "generated_automap": generated_automap,
        "genome_build": "GRCh38",
        "coordinate_system": "normalized 0-based half-open; UI 1-based inclusive",
        "default_filter": default_filter,
        "region_count_all": len(regions),
        "region_count_default": len(selected),
        "autosomal_region_count_all": len(all_autosomal),
        "autosomal_region_count_default": len(autosomal),
        "autosomal_total_bp_default": _merged_total_bp(autosomal),
        "autosomal_total_mb_default": round(_merged_total_bp(autosomal) / 1_000_000, 3),
        "warnings": warnings,
        "source_artifacts": artifact_names,
        "source_input_path": str(source_path) if source_path is not None else "",
        "source_aux_path": str(aux_path) if aux_path is not None else "",
        "source_details": parser_meta,
    }
    _atomic_regions(_scoped(post_dir, sample_id, REGIONS_NAME), regions)
    _atomic_json(_scoped(post_dir, sample_id, SUMMARY_NAME), summary)
    return summary


def load_regions(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    rows: list[dict] = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            start0, end0 = _number(row.get("start0"), int), _number(row.get("end0"), int)
            if start0 is None or end0 is None:
                continue
            item = dict(row)
            item.update({
                "start0": start0,
                "end0": end0,
                "display_start": _number(row.get("display_start"), int),
                "display_end": _number(row.get("display_end"), int),
                "length_bp": _number(row.get("length_bp"), int),
                "length_mb": _number(row.get("length_mb"), float),
                "passes_default_filter": str(row.get("passes_default_filter") or "").lower() in {"1", "true", "yes"},
            })
            for key in ("n_markers", "n_hom", "n_het"):
                item[key] = _number(row.get(key), int)
            for key in ("homozygosity_pct", "quality", "score"):
                item[key] = _number(row.get(key), float)
            rows.append(item)
    return rows


def load_sample_roh(sample_id: str) -> dict | None:
    sub = sample_layout.state_dir(sample_id)
    if not sub.is_dir():
        return None
    summary_path = sample_layout.state_file(sample_id, SUMMARY_NAME)
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        summary = {"status": "source_missing", "source": "missing", "warnings": ["ROH 尚未完成 postprocessing"]}
    regions = load_regions(sample_layout.state_file(sample_id, REGIONS_NAME))
    return {"roh_summary": summary, "roh_regions": regions, "roh_pending": False}


def annotate_variants(variants: dict[str, dict], sample_id: str) -> None:
    """Set ``in_roh`` from the selected source's default-filter intervals."""
    regions_path = sample_layout.state_file(sample_id, REGIONS_NAME)
    # Preserve legacy TSV-provided IN_ROH values until that case has gone
    # through the new normalized ROH postprocessing.
    if not regions_path.is_file():
        return
    regions = [
        row for row in load_regions(regions_path)
        if row.get("passes_default_filter")
    ]
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for row in regions:
        by_chrom.setdefault(row["chrom"], []).append((row["start0"], row["end0"]))
    spans_by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, spans in by_chrom.items():
        merged: list[list[int]] = []
        for start, end in sorted(spans):
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        spans_by_chrom[chrom] = [(start, end) for start, end in merged]
    starts = {chrom: [span[0] for span in spans] for chrom, spans in spans_by_chrom.items()}
    for variant in variants.values():
        chrom = _chrom(variant.get("CHROM") or variant.get("chrom"))
        pos = _number(variant.get("POS") or variant.get("pos"), int)
        inside = False
        if chrom and pos is not None and chrom in starts:
            point = pos - 1
            index = bisect_right(starts[chrom], point) - 1
            if index >= 0:
                start, end = spans_by_chrom[chrom][index]
                inside = start <= point < end
        variant["in_roh"] = inside
