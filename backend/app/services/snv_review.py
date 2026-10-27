"""Build the compact SNV/Indel TSV used by the main review screen.

The pipeline-owned 03_acmg TSV remains the complete source of truth. This
derived file is intentionally disposable: when the raw TSV or overlay
changes, the next sample load rebuilds it atomically.
"""
from __future__ import annotations

import csv
import bisect
import json
import os
import re
import time
from pathlib import Path

from . import gpn_msa, gpn_star
from .snv_overlay import OverlayReader, overlay_signature, row_key
from .snv_rows import is_reportable_raw_row

REVIEW_TSV_NAME = "snv_indel.review.tsv"
MAX_GNOMAD_G_AF = 0.01
WES_DP_HARD_FLOOR = 20
DEFAULT_CANDIDATE_BED = Path.home() / "NGS_UI" / "biotools" / "cds_combined.bed"
BedIndex = dict[str, tuple[list[int], list[tuple[int, int]]]]

_PATHOGENIC_RE = re.compile(r"(?:^|[/|,; ])(?:likely[_ ]?)?pathogenic(?:$|[/|,; ])", re.I)
_CLINVAR_CHANGE_RESCUES = {"UP_TO_PLP", "DOWN_FROM_PLP"}


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


def _write_manifest(path: Path, payload: dict[str, object]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    os.replace(tmp, path)


def _to_float(raw: str) -> float | None:
    value = (raw or "").strip()
    if not value or value in {".", "NA", "N/A"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _to_int(raw: str) -> int:
    value = (raw or "").strip()
    if not value or value in {".", "NA", "N/A"}:
        return 0
    try:
        return int(float(value))
    except ValueError:
        return 0


def _depth_from_ad(ad: str) -> int:
    total = 0
    for part in str(ad or "").split(","):
        value = part.strip()
        if not value or value == "." or value.upper() in {"NA", "N/A"}:
            continue
        try:
            total += int(float(value))
        except ValueError:
            continue
    return total


def _coalesce(*vals: str) -> str:
    for value in vals:
        text = (value or "").strip()
        if text and text not in {".", "NA", "N/A"}:
            return text
    return ""


def _row_depth(row: dict[str, str]) -> int:
    return (
        _to_int(row.get("DP") or "")
        or _to_int(row.get("DP_DV") or "")
        or _to_int(row.get("DP_HC") or "")
        or _depth_from_ad(_coalesce(row.get("AD"), row.get("AD_DV"), row.get("AD_HC")))
    )


def _candidate_bed_path() -> Path:
    return Path(os.environ.get("NGS_UI_CDS_CANDIDATE_BED") or DEFAULT_CANDIDATE_BED)


def _norm_chrom(chrom: str) -> str:
    s = chrom.strip()
    if s.lower().startswith("chr"):
        s = s[3:]
    if s in {"M", "MT", "m", "mt"}:
        return "MT"
    return s.upper() if s.upper() in {"X", "Y"} else s


def _is_mito_chrom(chrom: str) -> bool:
    return _norm_chrom(chrom or "") == "MT"


def _load_bed(path: Path) -> BedIndex | None:
    if not path.is_file():
        return None
    raw: dict[str, list[tuple[int, int]]] = {}
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if end <= start:
                continue
            raw.setdefault(_norm_chrom(parts[0]), []).append((start, end))

    out: BedIndex = {}
    for chrom, intervals in raw.items():
        intervals.sort()
        merged: list[tuple[int, int]] = []
        for start, end in intervals:
            if not merged or start > merged[-1][1]:
                merged.append((start, end))
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        out[chrom] = ([start for start, _ in merged], merged)
    return out


def _bed_signature(path: Path) -> dict[str, object]:
    try:
        st = path.stat()
    except OSError:
        return {"path": str(path), "exists": False}
    return {
        "path": str(path),
        "exists": True,
        "mtime_ns": st.st_mtime_ns,
        "size": st.st_size,
    }


def _variant_span_0based(row: dict[str, str]) -> tuple[int, int] | None:
    try:
        start = int((row.get("POS") or "").strip()) - 1
    except ValueError:
        return None
    if start < 0:
        return None
    ref = (row.get("REF") or "").strip()
    return start, start + max(1, len(ref))


def _overlaps_bed(row: dict[str, str], bed: BedIndex | None) -> bool:
    if bed is None:
        return True
    span = _variant_span_0based(row)
    if span is None:
        return False
    start, end = span
    chrom = _norm_chrom(row.get("CHROM") or "")
    item = bed.get(chrom)
    if not item:
        return False
    starts, intervals = item
    i = bisect.bisect_right(starts, start) - 1
    if i >= 0 and intervals[i][1] > start:
        return True
    j = i + 1
    return j < len(intervals) and intervals[j][0] < end


def _passes_basic_review_filter(row: dict[str, str], *, is_wes: bool) -> bool:
    """Apply invariant row eligibility before any clinical rescue rule."""
    if not is_reportable_raw_row(row):
        return False
    if _is_mito_chrom(row.get("CHROM") or ""):
        return False
    if is_wes and _row_depth(row) < WES_DP_HARD_FLOOR:
        return False
    return True


def _keep_row(row: dict[str, str], bed: BedIndex | None, *, is_wes: bool) -> bool:
    """Keep ClinVar P/LP rescue rows, then rare/unknown-AF BED rows."""
    if not _passes_basic_review_filter(row, is_wes=is_wes):
        return False
    sig = (row.get("CLINVAR_SIG") or "").strip()
    if _PATHOGENIC_RE.search(sig):
        return True
    af = _to_float(row.get("GNOMAD_G_AF") or "")
    return (af is None or af < MAX_GNOMAD_G_AF) and _overlaps_bed(row, bed)


def load_candidate_bed(path: Path | None = None) -> BedIndex | None:
    """Load the configured review candidate BED.

    GeneBe live-API fallback uses this same public helper so its network
    candidate scope cannot drift from the rows retained in review.tsv.
    A missing BED intentionally preserves the existing review behaviour:
    rare/unknown-AF rows are not BED-restricted.
    """
    return _load_bed(Path(path) if path is not None else _candidate_bed_path())


def is_review_candidate(
    row: dict[str, str],
    *,
    test_type: str,
    bed: BedIndex | None,
) -> bool:
    """Return whether a raw row belongs to the compact review candidate set."""
    return _keep_row(row, bed, is_wes=(test_type or "WES").upper() == "WES")


def is_review_retained(
    row: dict[str, str],
    *,
    test_type: str,
    bed: BedIndex | None,
    clinically_changed: bool | None = None,
) -> bool:
    """Return the authoritative final review-TSV retention decision.

    SpliceAI, GeneBe live fallback, LitVar2, and review TSV materialisation all
    call this function.  ``clinically_changed`` is only supplied while reading
    the immutable raw TSV, where CLINVAR_CHANGE still lives in the overlay.
    """
    is_wes = (test_type or "WES").upper() == "WES"
    if _keep_row(row, bed, is_wes=is_wes):
        return True
    if clinically_changed is None:
        clinically_changed = (
            str(row.get("CLINVAR_CHANGE") or "").strip().upper()
            in _CLINVAR_CHANGE_RESCUES
        )
    return bool(clinically_changed) and _passes_basic_review_filter(
        row, is_wes=is_wes
    )


def _refresh_gpn_annotations(
    tsv: Path, source: dict, previous: dict, *,
    msa_db: Path, star_dir: Path, research_only: bool,
    require_msa: bool, require_star: bool, force_msa: bool, force_star: bool,
) -> tuple[dict, bool]:
    """Refresh each independent predictor only when its own signature changes."""
    payload = {**previous, **source}
    refreshed = False
    for name, module, database, required, force, options in (
        ("gpn_msa", gpn_msa, msa_db, require_msa, force_msa, {}),
        ("gpn_star", gpn_star, star_dir, require_star, force_star, {}),
        ("gpn_star_v", gpn_star, star_dir, require_star, force_star, {"model": gpn_star.V_MODEL}),
    ):
        if name.startswith("gpn_star") and not research_only:
            # This controls new annotation, not visibility of existing scores.
            # Keep previous *_values provenance for scores already on disk.
            payload[f"{name}_annotation"] = {"status": "skipped_not_research_only"}
            continue
        status = previous.get(f"{name}_annotation") or {}
        if (previous.get(name) == source[name] and status.get("status")
                and not force and (not required or status["status"] == "complete")):
            continue
        stats = module.annotate_review_tsv(tsv, database, required=required, **options)
        payload[f"{name}_annotation"] = stats
        if stats["status"] == "complete":
            payload[f"{name}_values"] = {
                **source[name],
                **({"revision": stats["revision"]} if name.startswith("gpn_star") else {}),
            }
        refreshed = True
    return payload, refreshed or payload != previous


def ensure_review_tsv(
    raw_tsv: Path,
    *,
    keep_ids: set[str] | None = None,
    test_type: str = "WES",
    output_dir: Path | None = None,
    output_path: Path | None = None,
    manifest_path: Path | None = None,
    overlay_path: Path | None = None,
    gpn_msa_db: Path | None = None,
    require_gpn_msa: bool = False,
    force_gpn_msa: bool = False,
    research_only: bool | None = None,
    gpn_star_dir: Path | None = None,
    require_gpn_star: bool = False,
    force_gpn_star: bool = False,
) -> Path:
    """Return a compact review TSV with MSA and opt-in Research-only Star.

    Loaders omit research_only to inherit the persisted explicit opt-in.
    Legacy manifests without that flag never enable new Star annotation.
    Pipeline calls always pass a bool, so a new run records its own choice.
    """
    started = time.perf_counter()
    keep_ids = keep_ids or set()
    test_type_key = (test_type or "WES").upper()
    is_wes = test_type_key == "WES"
    output_dir = (
        Path(output_path).parent
        if output_path is not None
        else Path(output_dir) if output_dir else raw_tsv.parent
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    review_tsv = Path(output_path) if output_path is not None else output_dir / REVIEW_TSV_NAME
    manifest = (
        Path(manifest_path)
        if manifest_path is not None
        else review_tsv.with_suffix(review_tsv.suffix + ".source.json")
    )
    try:
        current = json.loads(manifest.read_text(encoding="utf-8"))
        if not isinstance(current, dict):
            current = {}
    except (OSError, json.JSONDecodeError):
        current = {}
    star_enabled = (research_only if research_only is not None
                    else current.get("gpn_star_research_only") is True)
    if require_gpn_star and not star_enabled:
        raise ValueError("GPN-Star requires explicit Research-only opt-in")
    candidate_bed_path = _candidate_bed_path()
    selected_gpn_msa_db = Path(gpn_msa_db or gpn_msa.GPN_MSA_DB)
    selected_gpn_star_dir = Path(gpn_star_dir or gpn_star.GPN_STAR_DIR)
    if require_gpn_msa:
        gpn_msa.validate_database(selected_gpn_msa_db)
    if require_gpn_star:
        for model in gpn_star.MODELS:
            gpn_star.validate_database(selected_gpn_star_dir, model=model)
    source = {
        "raw_mtime_ns": raw_tsv.stat().st_mtime_ns,
        "raw_size": raw_tsv.stat().st_size,
        "max_gnomad_g_af": MAX_GNOMAD_G_AF,
        "test_type": test_type_key,
        "wes_dp_hard_floor": WES_DP_HARD_FLOOR if is_wes else None,
        "candidate_bed": _bed_signature(candidate_bed_path),
        "overlay": overlay_signature(raw_tsv, overlay_path),
        "gpn_msa": gpn_msa.database_signature(selected_gpn_msa_db),
        "gpn_star_research_only": star_enabled,
        "gpn_star_v": (gpn_star.database_signature(selected_gpn_star_dir, model=gpn_star.V_MODEL)
                       if star_enabled else {"enabled": False}),
        "gpn_star": (gpn_star.database_signature(selected_gpn_star_dir)
                     if star_enabled else {"enabled": False}),
    }
    annotation_options = {
        "msa_db": selected_gpn_msa_db, "star_dir": selected_gpn_star_dir,
        "research_only": star_enabled,
        "require_msa": require_gpn_msa, "require_star": require_gpn_star,
        "force_msa": force_gpn_msa, "force_star": force_gpn_star,
    }
    non_gpn_matches = all(
        current.get(key) == value for key, value in source.items()
        if key not in {"gpn_msa", "gpn_star", "gpn_star_v", "gpn_star_research_only"}
    )
    if review_tsv.is_file() and non_gpn_matches:
        payload, refreshed = _refresh_gpn_annotations(
            review_tsv, source, current, **annotation_options,
        )
        if refreshed:
            _write_manifest(manifest, payload)
        _log_perf(
            "snv_review.ensure", started,
            action="gpn-refresh" if refreshed else "hit",
            raw_size=_fmt_size(raw_tsv), review_size=_fmt_size(review_tsv),
            keep_ids=len(keep_ids),
            gpn_msa_status=payload["gpn_msa_annotation"]["status"],
            gpn_star_status=payload["gpn_star_annotation"]["status"],
            gpn_star_v_status=payload["gpn_star_v_annotation"]["status"],
        )
        return review_tsv

    candidate_bed = _load_bed(candidate_bed_path)
    tmp = review_tsv.with_suffix(review_tsv.suffix + ".tmp")
    kept = 0
    scanned = 0
    with raw_tsv.open("r", encoding="utf-8", newline="") as src, \
            OverlayReader(raw_tsv, overlay_path) as overlay:
        clinvar_change_keys = overlay.keys_with_field_values(
            "CLINVAR_CHANGE", {"UP_TO_PLP", "DOWN_FROM_PLP"}
        )
        reader = csv.DictReader(src, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for field in overlay.fields:
            if field not in fieldnames:
                fieldnames.append(field)
        with tmp.open("w", encoding="utf-8", newline="") as dst:
            writer = csv.DictWriter(
                dst, fieldnames=fieldnames, delimiter="\t",
                extrasaction="ignore", lineterminator="\n",
            )
            writer.writeheader()
            batch: list[dict[str, str]] = []
            for row in reader:
                scanned += 1
                clinically_changed = row_key(row) in clinvar_change_keys
                keep = is_review_retained(
                    row,
                    test_type=test_type_key,
                    bed=candidate_bed,
                    clinically_changed=clinically_changed,
                )
                if keep:
                    batch.append(row)
                    kept += 1
                    if len(batch) >= 1000:
                        writer.writerows(overlay.apply_many(batch))
                        batch.clear()
            if batch:
                writer.writerows(overlay.apply_many(batch))
    try:
        payload, _ = _refresh_gpn_annotations(tmp, source, {}, **annotation_options)
        os.replace(tmp, review_tsv)
    except Exception:
        try:
            tmp.unlink()
        except OSError:
            pass
        raise
    _write_manifest(manifest, payload)
    _log_perf(
        "snv_review.ensure",
        started,
        action="rebuild",
        raw_size=_fmt_size(raw_tsv),
        review_size=_fmt_size(review_tsv),
        scanned=scanned,
        kept=kept,
        keep_ids=len(keep_ids),
        bed_exists=bool(candidate_bed),
        test_type=test_type_key,
        gpn_msa_status=payload["gpn_msa_annotation"]["status"],
        gpn_star_status=payload["gpn_star_annotation"]["status"],
        gpn_star_v_status=payload["gpn_star_v_annotation"]["status"],
    )
    return review_tsv
