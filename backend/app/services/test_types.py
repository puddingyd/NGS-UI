"""Canonical sample test-type labels and TITAN-WGS inference."""
from __future__ import annotations

import re


WES = "WES"
WGS = "WGS"
TITAN_WGS = "TITAN-WGS"
VALID_TEST_TYPES = frozenset({WES, WGS, TITAN_WGS})

_TITAN_SAMPLE_ID_RE = re.compile(r"^\d{2}T", re.IGNORECASE)


def is_titan_sample_id(sample_id: str) -> bool:
    """Return True for year-plus-T LIS IDs such as 25T..., 26T..., 27T...."""
    return bool(_TITAN_SAMPLE_ID_RE.match(str(sample_id or "").strip()))


def normalize_test_type(
    value: str = "",
    *,
    sample_id: str = "",
    default: str = WES,
) -> str:
    """Normalize labels and force year-plus-T sample IDs into TITAN-WGS."""
    if is_titan_sample_id(sample_id):
        return TITAN_WGS
    normalized = str(value or "").strip().upper()
    return normalized or default


def is_wgs_type(value: str) -> bool:
    """TITAN-WGS uses the same analysis/report behavior as ordinary WGS."""
    return str(value or "").strip().upper() in {WGS, TITAN_WGS}
