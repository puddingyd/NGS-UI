"""Patient-level phenotype.txt persistence and legacy resolution.

The canonical file is ``patient_phenotype/{MRN}_phenotype.txt`` so a later
specimen for the same patient can reuse the reviewer-curated HPO/panel set.
Older LIS-specific files remain readable as fallbacks but are never allowed to
override the MRN-only file.
"""
from __future__ import annotations

import json
import logging
import os
import re
import threading
from datetime import datetime, timezone
from pathlib import Path

from ..config import PHENOTYPE_DIR
from . import phenotype_io, phenotype_scorer


_TOKEN_RE = re.compile(r"^[A-Za-z0-9_-]{1,32}$")
MAX_CONTENT_BYTES = 64 * 1024
_WRITE_LOCK = threading.RLock()
logger = logging.getLogger(__name__)


def check_token(name: str, value: str, *, required: bool) -> str:
    token = str(value or "").strip()
    if not token:
        if required:
            raise ValueError(f"{name} 為必填")
        return ""
    if not _TOKEN_RE.match(token):
        raise ValueError(f"{name} 只能是英數 / - / _（最多 32 字）")
    return token


def filename_for(*, mrn: str = "", code: str = "") -> str:
    """Return the canonical filename, with LIS_ID as a no-MRN fallback."""
    mrn = check_token("MRN", mrn, required=False)
    code = check_token("LIS_ID", code, required=False)
    if mrn:
        return f"{mrn}_phenotype.txt"
    if code:
        return f"{code}_phenotype.txt"
    raise ValueError("請至少提供 MRN 或 LIS_ID")


def _path_for(*, mrn: str = "", code: str = "") -> Path:
    return Path(PHENOTYPE_DIR) / filename_for(mrn=mrn, code=code)


def _atomic_write(hpo: list, panels: list, out: Path, *, overwrite: bool = True) -> bool:
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.resolve().parent != Path(PHENOTYPE_DIR).resolve():
        raise ValueError("檔名不合法")
    tmp = out.with_name(f".{out.name}.{os.getpid()}.tmp")
    with _WRITE_LOCK:
        try:
            phenotype_io.write(hpo or [], panels or [], tmp)
            if overwrite:
                os.replace(tmp, out)
            else:
                # Publish only if still absent, including across web workers.
                # A concurrent manual save must always remain authoritative.
                try:
                    os.link(tmp, out)
                except FileExistsError:
                    return False
            return True
        finally:
            try:
                tmp.unlink(missing_ok=True)
            except OSError:
                pass


def save(
    *,
    mrn: str = "",
    code: str = "",
    hpo: list | None = None,
    panels: list | None = None,
) -> dict:
    """Atomically save a canonical patient phenotype snapshot.

    An empty HPO/panel set intentionally writes a header-only file. This is an
    explicit empty marker and prevents a stale legacy file from resurfacing.
    """
    mrn = check_token("MRN", mrn, required=False)
    code = check_token("LIS_ID", code, required=False)
    out = _path_for(mrn=mrn, code=code)
    _atomic_write(hpo or [], panels or [], out)
    return {
        "path": str(out),
        "filename": out.name,
        "mrn": mrn,
        "code": code,
        "n_hpo": len(hpo or []),
        "n_panels": len(panels or []),
    }


def save_text(*, mrn: str = "", code: str = "", content: str = "") -> dict:
    """Parse and save a phenotype.txt body through the canonical writer."""
    if not isinstance(content, str) or not content.strip():
        raise ValueError("內容為空")
    if len(content.encode("utf-8")) > MAX_CONTENT_BYTES:
        raise ValueError("內容過大")
    hpo, panels = phenotype_io.parse(content)
    return save(mrn=mrn, code=code, hpo=hpo, panels=panels)


def find(
    *,
    mrn: str = "",
    code: str = "",
    code_candidates: list[str] | tuple[str, ...] | None = None,
) -> Path | None:
    """Find phenotype data, always preferring the MRN-only canonical file."""
    mrn = check_token("MRN", mrn, required=False)
    code = check_token("LIS_ID", code, required=False)
    codes: list[str] = []
    for raw in [code, *(code_candidates or [])]:
        token = check_token("LIS_ID", raw or "", required=False)
        if token and token not in codes:
            codes.append(token)
    if not mrn and not codes:
        raise ValueError("請至少提供 MRN 或 LIS_ID")

    root = Path(PHENOTYPE_DIR)
    if not root.is_dir():
        return None

    if mrn:
        canonical = root / f"{mrn}_phenotype.txt"
        if canonical.is_file():
            return canonical

    if mrn:
        for candidate_code in codes:
            legacy = root / f"{candidate_code}_{mrn}_phenotype.txt"
            if legacy.is_file():
                return legacy

    for candidate_code in codes:
        legacy = root / f"{candidate_code}_phenotype.txt"
        if legacy.is_file():
            return legacy

    for candidate_code in codes:
        matches = sorted(root.glob(f"{candidate_code}_*_phenotype.txt"))
        if matches:
            return matches[0]

    if mrn:
        matches = sorted(root.glob(f"*_{mrn}_phenotype.txt"))
        if matches:
            return matches[0]
    return None


def _parse_filename(path: Path, *, mrn: str = "", codes: list[str] | None = None) -> dict:
    core = path.stem
    if core.endswith("_phenotype"):
        core = core[:-len("_phenotype")]
    if mrn and core == mrn:
        return {"code": "", "mrn": mrn}
    for code in codes or []:
        if core == code:
            return {"code": code, "mrn": ""}
        if core.startswith(code + "_"):
            return {"code": code, "mrn": core[len(code) + 1:]}
    if mrn and core.endswith("_" + mrn):
        return {"code": core[:-(len(mrn) + 1)], "mrn": mrn}
    return {"code": "", "mrn": core}


def _snapshot_time(value: str | None) -> float:
    try:
        parsed = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
        return parsed.replace(tzinfo=parsed.tzinfo or timezone.utc).timestamp()
    except (TypeError, ValueError, OverflowError):
        return 0.0


def backfill_registered(*, mrn: str = "") -> dict:
    """Create missing MRN snapshots from saved default analyses, without scoring.

    Existing patient files (including explicit empty/legacy files) win. For
    multiple specimens, use the most recently saved default, without merging
    experimental analysis versions or different patients' selections.
    """
    from . import category_store, sample_layout

    requested = check_token("MRN", mrn, required=False)
    candidates = {}
    blocked = set()
    created = errors = 0
    # Coordinate with metadata/MRN edits, unregister and analysis publication.
    # Only small JSON files are read; no variant files or gene scoring here.
    with category_store.locked():
        for sample_id in sorted(sample_layout.iter_sample_ids()):
            path = sample_layout.state_file(sample_id, "sample_metadata.json")
            try:
                meta = json.loads(path.read_text(encoding="utf-8"))
            except FileNotFoundError:
                continue  # Pipeline output without a saved registration.
            except (OSError, ValueError):
                errors += 1
                continue
            if not isinstance(meta, dict):
                errors += 1
                continue
            try:
                patient_mrn = check_token("MRN", meta.get("mrn"), required=True)
            except ValueError:
                continue
            if requested and patient_mrn != requested:
                continue
            if find(mrn=patient_mrn) is not None:
                continue
            default_path = sample_layout.analysis_file(sample_id, "default", "analysis.json")
            try:
                if default_path.is_file():
                    snapshot = json.loads(default_path.read_text(encoding="utf-8"))
                else:
                    snapshot = meta  # Pre-versioning registrations stored selections here.
                if not isinstance(snapshot, dict):
                    raise ValueError("invalid phenotype snapshot")
                if not any(key in snapshot for key in ("hpo", "patient_phenotype", "selected_panels")):
                    continue  # No saved default; never infer from an experimental version.
                hpo = snapshot.get("hpo", snapshot.get("patient_phenotype", [])) or []
                panels = snapshot.get("selected_panels") or []
                if not isinstance(hpo, list) or not isinstance(panels, list):
                    raise ValueError("invalid phenotype lists")
                if any(not isinstance(row, (dict, str)) for row in hpo + panels):
                    raise ValueError("invalid phenotype entry")
                hpo = [row if isinstance(row, dict) else {"phenotype": row} for row in hpo]
                panels = phenotype_scorer.normalize_panel_entries(panels)
            except (OSError, ValueError):
                # Do not substitute an older specimen when a saved default
                # for this same MRN cannot be read reliably.
                blocked.add(patient_mrn)
                errors += 1
                continue
            saved_at = (snapshot.get("updated_at") or meta.get("phenotype_updated_at")
                        or meta.get("created_at") or meta.get("run_date"))
            rank = (_snapshot_time(saved_at), sample_id)
            if patient_mrn not in candidates or rank > candidates[patient_mrn][0]:
                candidates[patient_mrn] = (rank, hpo, panels)

        for patient_mrn, (_rank, hpo, panels) in candidates.items():
            if patient_mrn in blocked or find(mrn=patient_mrn) is not None:
                continue
            if _atomic_write(hpo, panels, _path_for(mrn=patient_mrn), overwrite=False):
                created += 1
    if errors:
        logger.warning("patient phenotype backfill skipped %d unreadable/invalid saved records", errors)
    return {"created": created, "errors": errors}


def load(
    *,
    mrn: str = "",
    code: str = "",
    code_candidates: list[str] | tuple[str, ...] | None = None,
) -> dict:
    """Return the resolved file, parsed HPO/panels, and original content."""
    safe_mrn = check_token("MRN", mrn, required=False)
    safe_code = check_token("LIS_ID", code, required=False)
    codes: list[str] = []
    for raw in [safe_code, *(code_candidates or [])]:
        token = check_token("LIS_ID", raw or "", required=False)
        if token and token not in codes:
            codes.append(token)
    path = find(mrn=safe_mrn, code=safe_code, code_candidates=codes)
    if path is None and safe_mrn:
        backfill_registered(mrn=safe_mrn)
        path = find(mrn=safe_mrn, code=safe_code, code_candidates=codes)
    if path is None:
        return {}
    try:
        content = path.read_text(encoding="utf-8")
    except OSError:
        return {}
    hpo, panels = phenotype_io.parse(content)
    panels = phenotype_scorer.normalize_panel_entries(panels)
    return {
        "filename": path.name,
        "path": str(path),
        "content": content,
        "hpo": hpo,
        "panels": panels,
        **_parse_filename(path, mrn=safe_mrn, codes=codes),
    }
