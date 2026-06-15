"""MRN/LIS_ID keyed Clinical presentation sidecar files.

The standalone /phenotype/ tool can be used before a sample is registered,
so this free text lives next to phenotype.txt and is pulled into
sample_metadata.json when the case is opened in the main reviewer UI.
"""
from __future__ import annotations

import re
from pathlib import Path

from ..config import PHENOTYPE_DIR

_TOKEN_RE = re.compile(r"^[A-Za-z0-9_-]{1,32}$")
_SUFFIX = "_clinical_presentation.txt"
MAX_CONTENT_BYTES = 64 * 1024


def check_token(name: str, value: str, *, required: bool) -> str:
    v = (value or "").strip()
    if not v:
        if required:
            raise ValueError(f"{name} 為必填")
        return ""
    if not _TOKEN_RE.match(v):
        raise ValueError(f"{name} 只能是英數 / - / _（最多 32 字）")
    return v


def filename_for(*, code: str = "", mrn: str = "") -> str:
    code = check_token("LIS_ID", code, required=False)
    mrn = check_token("MRN", mrn, required=False)
    if code and mrn:
        return f"{code}_{mrn}{_SUFFIX}"
    if code:
        return f"{code}{_SUFFIX}"
    if mrn:
        return f"{mrn}{_SUFFIX}"
    raise ValueError("請至少提供 MRN 或 LIS_ID")


def _path_for(*, code: str = "", mrn: str = "") -> Path:
    return PHENOTYPE_DIR / filename_for(code=code, mrn=mrn)


def _read(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except OSError:
        return ""


def find(*, code: str = "", mrn: str = "", code_candidates: list[str] | tuple[str, ...] | None = None) -> Path | None:
    """Find the best clinical presentation sidecar.

    Precedence mirrors phenotype.txt lookup: exact LIS+MRN, LIS-only, any
    LIS+MRN for the LIS candidate, then MRN-only / any LIS for MRN.
    """
    code = check_token("LIS_ID", code, required=False)
    mrn = check_token("MRN", mrn, required=False)
    candidates = []
    seen = set()
    for item in [code, *(code_candidates or [])]:
        item = check_token("LIS_ID", item or "", required=False)
        if item and item not in seen:
            candidates.append(item)
            seen.add(item)

    if not PHENOTYPE_DIR.is_dir():
        return None

    for c in candidates:
        if mrn:
            p = PHENOTYPE_DIR / f"{c}_{mrn}{_SUFFIX}"
            if p.is_file():
                return p
        p = PHENOTYPE_DIR / f"{c}{_SUFFIX}"
        if p.is_file():
            return p
        matches = sorted(PHENOTYPE_DIR.glob(f"{c}_*{_SUFFIX}"))
        if matches:
            return matches[0]

    if mrn:
        p = PHENOTYPE_DIR / f"{mrn}{_SUFFIX}"
        if p.is_file():
            return p
        matches = sorted(PHENOTYPE_DIR.glob(f"*_{mrn}{_SUFFIX}"))
        if matches:
            return matches[0]
    return None


def load(*, code: str = "", mrn: str = "", code_candidates: list[str] | tuple[str, ...] | None = None) -> dict:
    path = find(code=code, mrn=mrn, code_candidates=code_candidates)
    if path is None:
        return {}
    return {
        "filename": path.name,
        "path": str(path),
        "content": _read(path),
        **parse_filename(path.name, code=code, mrn=mrn),
    }


def save(*, code: str = "", mrn: str = "", content: str = "") -> dict:
    code = check_token("LIS_ID", code, required=False)
    mrn = check_token("MRN", mrn, required=False)
    if not code and not mrn:
        raise ValueError("請至少提供 MRN 或 LIS_ID")
    if not isinstance(content, str):
        raise ValueError("Clinical presentation 內容格式不合法")
    if len(content.encode("utf-8")) > MAX_CONTENT_BYTES:
        raise ValueError("Clinical presentation 內容過大")

    PHENOTYPE_DIR.mkdir(parents=True, exist_ok=True)
    out = _path_for(code=code, mrn=mrn)
    if out.resolve().parent != PHENOTYPE_DIR.resolve():
        raise ValueError("檔名不合法")
    out.write_text(content if not content or content.endswith("\n") else content + "\n", encoding="utf-8")
    return {"path": str(out), "filename": out.name, "mrn": mrn, "code": code}


def parse_filename(filename: str, *, code: str = "", mrn: str = "") -> dict:
    stem = filename[:-len(_SUFFIX)] if filename.endswith(_SUFFIX) else Path(filename).stem
    if code and stem == code:
        return {"code": code, "mrn": ""}
    if mrn and stem == mrn:
        return {"code": "", "mrn": mrn}
    if code and stem.startswith(code + "_"):
        return {"code": code, "mrn": stem[len(code) + 1:]}
    parts = stem.split("_")
    if len(parts) >= 2:
        return {"code": parts[0], "mrn": parts[1]}
    return {"code": "", "mrn": stem}
