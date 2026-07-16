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

from .snv_overlay import OverlayReader, overlay_signature
from .snv_rows import is_reportable_raw_row

REVIEW_TSV_NAME = "snv_indel.review.tsv"
MAX_GNOMAD_G_AF = 0.01
WES_DP_HARD_FLOOR = 20
DEFAULT_CANDIDATE_BED = Path.home() / "NGS_UI" / "biotools" / "cds_combined.bed"
BedIndex = dict[str, tuple[list[int], list[tuple[int, int]]]]

_PATHOGENIC_RE = re.compile(r"(?:^|[/|,; ])(?:likely[_ ]?)?pathogenic(?:$|[/|,; ])", re.I)


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


def _keep_row(row: dict[str, str], bed: BedIndex | None, *, is_wes: bool) -> bool:
    """Keep ClinVar P/LP rescue rows, then rare/unknown-AF BED rows."""
    if not is_reportable_raw_row(row):
        return False
    if _is_mito_chrom(row.get("CHROM") or ""):
        return False
    if is_wes and _row_depth(row) < WES_DP_HARD_FLOOR:
        return False
    sig = (row.get("CLINVAR_SIG") or "").strip()
    if _PATHOGENIC_RE.search(sig):
        return True
    af = _to_float(row.get("GNOMAD_G_AF") or "")
    return (af is None or af < MAX_GNOMAD_G_AF) and _overlaps_bed(row, bed)


def ensure_review_tsv(
    raw_tsv: Path,
    *,
    keep_ids: set[str] | None = None,
    test_type: str = "WES",
    output_dir: Path | None = None,
    overlay_path: Path | None = None,
) -> Path:
    """Return an up-to-date compact review TSV derived from *raw_tsv*."""
    started = time.perf_counter()
    keep_ids = keep_ids or set()
    test_type_key = (test_type or "WES").upper()
    is_wes = test_type_key == "WES"
    output_dir = Path(output_dir) if output_dir else raw_tsv.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    review_tsv = output_dir / REVIEW_TSV_NAME
    manifest = review_tsv.with_suffix(review_tsv.suffix + ".source.json")
    candidate_bed_path = _candidate_bed_path()
    source = {
        "raw_mtime_ns": raw_tsv.stat().st_mtime_ns,
        "raw_size": raw_tsv.stat().st_size,
        "max_gnomad_g_af": MAX_GNOMAD_G_AF,
        "test_type": test_type_key,
        "wes_dp_hard_floor": WES_DP_HARD_FLOOR if is_wes else None,
        "candidate_bed": _bed_signature(candidate_bed_path),
        "overlay": overlay_signature(raw_tsv, overlay_path),
    }
    if review_tsv.is_file() and manifest.is_file():
        try:
            current = json.loads(manifest.read_text(encoding="utf-8"))
            if all(current.get(key) == value for key, value in source.items()):
                _log_perf(
                    "snv_review.ensure",
                    started,
                    action="hit",
                    raw_size=_fmt_size(raw_tsv),
                    review_size=_fmt_size(review_tsv),
                    keep_ids=len(keep_ids),
                )
                return review_tsv
        except (OSError, json.JSONDecodeError):
            pass

    candidate_bed = _load_bed(candidate_bed_path)
    tmp = review_tsv.with_suffix(review_tsv.suffix + ".tmp")
    manifest_tmp = manifest.with_suffix(manifest.suffix + ".tmp")
    kept = 0
    scanned = 0
    with raw_tsv.open("r", encoding="utf-8", newline="") as src, \
            OverlayReader(raw_tsv, overlay_path) as overlay:
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
                if _keep_row(row, candidate_bed, is_wes=is_wes):
                    batch.append(row)
                    kept += 1
                    if len(batch) >= 1000:
                        writer.writerows(overlay.apply_many(batch))
                        batch.clear()
            if batch:
                writer.writerows(overlay.apply_many(batch))
    os.replace(tmp, review_tsv)
    manifest_tmp.write_text(
        json.dumps(source, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    os.replace(manifest_tmp, manifest)
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
    )
    return review_tsv
