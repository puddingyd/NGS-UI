"""Read reproducibility metadata emitted beside tertiary annotation outputs.

The SNV TSV intentionally contains only row-level annotations.  Database
release dates belong in a compact per-sample sidecar so they are recorded once
and can be extended without changing the TSV schema.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any


_DATE_RE = re.compile(r"(?<!\d)(20\d{2})[-_]?([01]\d)[-_]?([0-3]\d)(?!\d)")
_CLINVAR_KEYS = ("clinvar", "clin_var", "clinvar_version")


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _normalise_date(value: Any) -> str:
    """Return an ISO release date, or an empty string for unknown values."""
    if isinstance(value, dict):
        for key in ("release_date", "version_date", "date", "version"):
            date = _normalise_date(value.get(key))
            if date:
                return date
        return ""
    if value is None:
        return ""
    match = _DATE_RE.search(str(value).strip())
    if not match:
        return ""
    year, month, day = match.groups()
    try:
        # Reject impossible dates such as 2026-19-99 while avoiding a heavier
        # dependency for one metadata field.
        from datetime import date

        return date(int(year), int(month), int(day)).isoformat()
    except ValueError:
        return ""


def _clinvar_date(payload: dict[str, Any]) -> str:
    direct = _normalise_date(payload.get("clinvar_date"))
    if direct:
        return direct

    for container in (
        payload,
        payload.get("annotation_versions"),
        payload.get("databases"),
        payload.get("annotations"),
    ):
        if not isinstance(container, dict):
            continue
        for key in _CLINVAR_KEYS:
            date = _normalise_date(container.get(key))
            if date:
                return date
    return ""


def sidecar_candidates(raw_tsv: Path, state_dir: Path | None = None) -> list[Path]:
    """Return supported metadata paths in precedence order.

    The preferred tertiary-pipeline filename is
    ``{source}.annotation_versions.json`` beside
    ``{source}.snv_indel.acmg.tsv``.  ``*.tsv.meta.json`` and the shared
    directory filename are accepted for pipelines that already have a generic
    output-metadata convention.  State-directory files are migration fallbacks.
    """
    raw_tsv = Path(raw_tsv)
    suffix = ".snv_indel.acmg.tsv"
    source_name = (
        raw_tsv.name[: -len(suffix)]
        if raw_tsv.name.endswith(suffix)
        else raw_tsv.stem
    )
    out = [
        raw_tsv.with_name(f"{source_name}.annotation_versions.json"),
        raw_tsv.with_name(f"{raw_tsv.name}.meta.json"),
        raw_tsv.parent / "annotation_versions.json",
    ]
    if state_dir is not None:
        state_dir = Path(state_dir)
        out.extend(
            [
                state_dir / "annotation_versions.json",
                state_dir / "pipeline_source.json",
            ]
        )
    deduped: list[Path] = []
    for path in out:
        if path not in deduped:
            deduped.append(path)
    return deduped


def load_annotation_versions(
    raw_tsv: Path, state_dir: Path | None = None
) -> dict[str, Any]:
    """Load known annotation versions without inventing a release date."""
    for path in sidecar_candidates(raw_tsv, state_dir):
        payload = _read_json(path)
        release_date = _clinvar_date(payload)
        if not release_date:
            continue
        return {
            "clinvar": {"release_date": release_date},
            # Kept as a flat compatibility field for the current frontend and
            # report exporters.  New consumers should prefer the nested value.
            "clinvar_date": release_date,
            "metadata_path": str(path),
        }
    return {}
