"""Resolve unified and legacy tertiary-analysis sample paths.

New samples use this layout::

    <PIPELINE_OUT_ROOT>/<sample>/00-07...
    <PIPELINE_OUT_ROOT>/<sample>/08_postprocessing/...

The sample-prefixed ``08_postprocessing/{sample_id}.layout.json`` marker is
written last for layout v3.  Layout v2 used the unprefixed ``layout.json``.
Both remain readable so existing cases do not need an eager rename.
"""
from __future__ import annotations

import json
import os
import shutil
import time
from pathlib import Path
from typing import Iterable

from .. import config


POSTPROCESSING_DIRNAME = "08_postprocessing"
LAYOUT_MARKER_NAME = "layout.json"
LAYOUT_VERSION = 3
SUPPORTED_LAYOUT_VERSIONS = frozenset({2, 3})


def unified_root() -> Path:
    return Path(config.PIPELINE_OUT_ROOT)


def legacy_ui_root() -> Path:
    return Path(config.LEGACY_TERTIARY_OUTPUT_ROOT)


def legacy_pipeline_root() -> Path:
    return Path(config.LEGACY_PIPELINE_OUT_ROOT)


def unified_sample_dir(sample_id: str) -> Path:
    return unified_root() / sample_id


def unified_postprocessing_dir(sample_id: str) -> Path:
    return unified_sample_dir(sample_id) / POSTPROCESSING_DIRNAME


def prefixed_filename(sample_id: str, name: str) -> str:
    """Return the layout-v3 filename for one sample-owned artifact."""
    return f"{sample_id}.{name}"


def _scoped_candidates(directory: Path, sample_id: str, name: str) -> tuple[Path, Path]:
    return (
        directory / prefixed_filename(sample_id, name),
        directory / name,
    )


def scoped_file(
    directory: Path,
    sample_id: str,
    name: str,
    *,
    for_write: bool = False,
    force_prefixed: bool = False,
) -> Path:
    """Resolve a sample-owned file with layout-v3 → legacy precedence.

    Reads always prefer the sample-prefixed filename.  Ordinary writes keep a
    layout-v2 case internally consistent until it is explicitly reprocessed;
    new/v3 cases write the prefixed form.  Post-processing can pass
    ``force_prefixed=True`` while constructing a complete v3 state tree before
    atomically publishing the v3 marker.
    """
    prefixed, legacy = _scoped_candidates(Path(directory), sample_id, name)
    if not for_write:
        if prefixed.is_file():
            return prefixed
        if legacy.is_file():
            return legacy
        return prefixed
    if force_prefixed:
        return prefixed
    post = unified_postprocessing_dir(sample_id)
    try:
        in_unified_post = Path(directory).resolve().is_relative_to(post.resolve())
    except (OSError, ValueError):
        in_unified_post = False
    if not in_unified_post:
        return legacy
    prefixed_marker, legacy_marker = _scoped_candidates(
        post, sample_id, LAYOUT_MARKER_NAME
    )
    if prefixed_marker.is_file() or not legacy_marker.is_file():
        return prefixed
    return legacy


def postprocessing_file(
    sample_id: str,
    name: str,
    *,
    for_write: bool = False,
    force_prefixed: bool = False,
) -> Path:
    return scoped_file(
        unified_postprocessing_dir(sample_id),
        sample_id,
        name,
        for_write=for_write,
        force_prefixed=force_prefixed,
    )


def layout_marker_path(sample_id: str) -> Path:
    return postprocessing_file(sample_id, LAYOUT_MARKER_NAME)


def uses_unified_layout(sample_id: str) -> bool:
    path = layout_marker_path(sample_id)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return (
        isinstance(payload, dict)
        and payload.get("layout_version") in SUPPORTED_LAYOUT_VERSIONS
    )


def is_ui_ready(sample_id: str) -> bool:
    """True only after post-processing completed or for a legacy UI copy."""
    return uses_unified_layout(sample_id) or (
        legacy_ui_root() / sample_id / "snv_indel.annotated.tsv"
    ).is_file()


def state_dir(sample_id: str, *, for_write: bool = False) -> Path:
    """Return the directory containing UI-owned per-sample state.

    Writes for a new sample go to 08_postprocessing.  An old sample keeps
    writing to its legacy directory until migration activates the marker.
    """
    post = unified_postprocessing_dir(sample_id)
    if uses_unified_layout(sample_id):
        return post
    old = legacy_ui_root() / sample_id
    if old.is_dir():
        return old
    return post


def pipeline_sample_dir(sample_id: str) -> Path:
    """Return the best available 00-07 pipeline directory."""
    new = unified_sample_dir(sample_id)
    if new.is_dir():
        return new
    old = legacy_pipeline_root() / sample_id
    if old.is_dir():
        return old
    return new


def iter_sample_ids() -> Iterable[str]:
    """Yield the union of new pipeline, legacy pipeline, and legacy UI IDs."""
    seen: set[str] = set()
    for root in (unified_root(), legacy_ui_root(), legacy_pipeline_root()):
        if not root.is_dir():
            continue
        try:
            children = root.iterdir()
        except OSError:
            continue
        for child in children:
            if not child.is_dir() or child.name.startswith(("_", ".")):
                continue
            if child.name in seen:
                continue
            seen.add(child.name)
            yield child.name


def _first_file(directory: Path, exact_name: str, pattern: str) -> Path:
    exact = directory / exact_name
    if exact.is_file():
        return exact
    if directory.is_dir():
        try:
            matches = sorted(p for p in directory.glob(pattern) if p.is_file())
        except OSError:
            matches = []
        if matches:
            return matches[0]
    return exact


def _pipeline_source(sample_id: str) -> dict:
    path = state_file(sample_id, "pipeline_source.json")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _layout_marker(sample_id: str) -> dict:
    try:
        value = json.loads(layout_marker_path(sample_id).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def source_sample_id(sample_id: str) -> str:
    value = str(_pipeline_source(sample_id).get("source_sample_id") or "").strip()
    return value or sample_id


def _pipeline_candidates(sample_id: str) -> list[Path]:
    source_id = source_sample_id(sample_id)
    out: list[Path] = []
    for root in (unified_root(), legacy_pipeline_root()):
        for sid in (sample_id, source_id):
            candidate = root / sid
            if candidate not in out:
                out.append(candidate)
    return out


def snv_raw_tsv(sample_id: str) -> Path:
    """Return the immutable pipeline SNV source (or legacy full copy)."""
    if uses_unified_layout(sample_id):
        sample = unified_sample_dir(sample_id)
        marker = _layout_marker(sample_id)
        source = _pipeline_source(sample_id)
        raw_value = str(marker.get("raw_tsv") or source.get("source_path") or "")
        if raw_value:
            raw_path = Path(raw_value)
            if not raw_path.is_absolute():
                raw_path = sample / raw_path
            if raw_path.is_file():
                return raw_path
        source_id = source_sample_id(sample_id)
        return _first_file(
            sample / "03_acmg",
            f"{source_id}.snv_indel.acmg.tsv",
            "*.snv_indel.acmg.tsv",
        )

    # Old UI copies may contain post-processing fields that are absent from
    # the old 03_acmg source, so they stay authoritative until migration.
    legacy_copy = legacy_ui_root() / sample_id / "snv_indel.annotated.tsv"
    if legacy_copy.is_file():
        return legacy_copy
    for sample in _pipeline_candidates(sample_id):
        source_id = source_sample_id(sample_id)
        found = _first_file(
            sample / "03_acmg",
            f"{source_id}.snv_indel.acmg.tsv",
            "*.snv_indel.acmg.tsv",
        )
        if found.is_file():
            return found
    return unified_sample_dir(sample_id) / "03_acmg" / f"{sample_id}.snv_indel.acmg.tsv"


def _pipeline_artifact(
    sample_id: str,
    subdir: str,
    exact_suffix: str,
    pattern: str,
) -> Path:
    source_id = source_sample_id(sample_id)
    for sample in _pipeline_candidates(sample_id):
        found = _first_file(
            sample / subdir,
            f"{source_id}{exact_suffix}",
            pattern,
        )
        if found.is_file():
            return found
    return unified_sample_dir(sample_id) / subdir / f"{source_id}{exact_suffix}"


def review_tsv(sample_id: str, *, for_write: bool = False) -> Path:
    return state_file(sample_id, "snv_indel.review.tsv", for_write=for_write)


def review_manifest(sample_id: str, *, for_write: bool = False) -> Path:
    return state_file(
        sample_id,
        "snv_indel.review.tsv.source.json",
        for_write=for_write,
    )


def snv_overlay_path(sample_id: str, *, for_write: bool = False) -> Path:
    return state_file(sample_id, "snv_annotations.sqlite", for_write=for_write)


def snv_gene_index_path(sample_id: str, *, for_write: bool = False) -> Path:
    return state_file(sample_id, "snv_gene_index.sqlite", for_write=for_write)


def cnv_tsv(sample_id: str) -> Path:
    legacy = state_file(sample_id, "cnv.annotated.tsv")
    if not uses_unified_layout(sample_id) and legacy.is_file():
        return legacy
    source = _pipeline_artifact(
        sample_id, "06_cnv_sv", ".cnv.annotated.tsv", "*.cnv.annotated.tsv"
    )
    if source.is_file():
        return source
    return legacy


def sv_tsv(sample_id: str) -> Path:
    legacy = state_file(sample_id, "sv.annotated.tsv")
    if not uses_unified_layout(sample_id) and legacy.is_file():
        return legacy
    source = _pipeline_artifact(
        sample_id, "06_cnv_sv", ".sv.annotated.tsv", "*.sv.annotated.tsv"
    )
    if source.is_file():
        return source
    return legacy


def mito_tsv(sample_id: str) -> Path:
    # A locally enriched mito file is intentionally preferred; it is small and
    # may contain MITOMAP columns not present in the pipeline source.
    derived = state_file(sample_id, "mito.annotated.tsv")
    if derived.is_file():
        return derived
    source = _pipeline_artifact(sample_id, "04_mito", ".mito.tsv", "*.mito.tsv")
    if source.is_file():
        return source
    return derived


def str_tsv(sample_id: str) -> Path:
    if not uses_unified_layout(sample_id):
        for name in ("str.tsv", "str.annotated.tsv"):
            path = state_file(sample_id, name)
            if path.is_file():
                return path
    source = _pipeline_artifact(sample_id, "05_str", ".str.tsv", "*.str.tsv")
    if source.is_file():
        return source
    for name in ("str.tsv", "str.annotated.tsv"):
        path = state_file(sample_id, name)
        if path.is_file():
            return path
    return state_file(sample_id, "str.tsv")


def pgx_tsv(sample_id: str) -> Path:
    legacy = state_file(sample_id, "pgx.tsv")
    if not uses_unified_layout(sample_id) and legacy.is_file():
        return legacy
    source = _pipeline_artifact(sample_id, "07_pgx", ".pgx.tsv", "*.pgx.tsv")
    if source.is_file():
        return source
    return legacy


def pharmcat_json(sample_id: str) -> Path:
    legacy = state_file(sample_id, "pharmcat.report.json")
    if not uses_unified_layout(sample_id) and legacy.is_file():
        return legacy
    source_id = source_sample_id(sample_id)
    for sample in _pipeline_candidates(sample_id):
        directory = sample / "07_pgx"
        for name in (
            f"{source_id}.pharmcat.report.json",
            "pharmcat.report.json",
        ):
            path = directory / name
            if path.is_file():
                return path
        found = _first_file(directory, "pharmcat.report.json", "*.pharmcat.report.json")
        if found.is_file():
            return found
    return legacy


def state_file(sample_id: str, name: str, *, for_write: bool = False) -> Path:
    return scoped_file(
        state_dir(sample_id),
        sample_id,
        name,
        for_write=for_write,
    )


def state_file_candidates(sample_id: str, name: str) -> tuple[Path, Path]:
    """Return (sample-prefixed, legacy-unprefixed) candidates."""
    return _scoped_candidates(state_dir(sample_id), sample_id, name)


def analysis_file(
    sample_id: str,
    version: str,
    name: str,
    *,
    for_write: bool = False,
) -> Path:
    directory = state_dir(sample_id) / "analyses" / version
    return scoped_file(
        directory,
        sample_id,
        name,
        for_write=for_write,
    )


def analysis_file_candidates(
    sample_id: str,
    version: str,
    name: str,
) -> tuple[Path, Path]:
    directory = state_dir(sample_id) / "analyses" / version
    return _scoped_candidates(directory, sample_id, name)


def promote_state_tree_to_v3(sample_id: str) -> list[tuple[Path, Path]]:
    """Non-destructively copy remaining v2 files to sample-prefixed names.

    A full worker reprocess calls this immediately before publishing the v3
    marker. Existing files are never overwritten or deleted; freshly-built
    prefixed artifacts therefore win, while legacy reviewer/analysis state is
    copied forward for exact-name v3 reads.
    """
    post = unified_postprocessing_dir(sample_id)
    if not post.is_dir():
        return []
    copied: list[tuple[Path, Path]] = []
    for source in sorted(path for path in post.rglob("*") if path.is_file()):
        if source.name == LAYOUT_MARKER_NAME:
            continue
        if source.name.startswith((".", f"{sample_id}.")):
            continue
        destination = source.with_name(prefixed_filename(sample_id, source.name))
        if destination.exists():
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination, follow_symlinks=False)
        copied.append((source, destination))
    return copied


def write_layout_marker(
    sample_id: str,
    *,
    source_id: str | None = None,
    raw_tsv: Path | None = None,
    migration: bool = False,
) -> Path:
    """Atomically activate the unified layout for one sample."""
    post = unified_postprocessing_dir(sample_id)
    post.mkdir(parents=True, exist_ok=True)
    raw = Path(raw_tsv) if raw_tsv else snv_raw_tsv(sample_id)
    try:
        raw_value = str(raw.relative_to(unified_sample_dir(sample_id)))
    except ValueError:
        raw_value = str(raw)
    payload = {
        "layout_version": LAYOUT_VERSION,
        "sample_id": sample_id,
        "source_sample_id": source_id or sample_id,
        "raw_tsv": raw_value,
        "activated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "migration": bool(migration),
    }
    path = postprocessing_file(
        sample_id,
        LAYOUT_MARKER_NAME,
        for_write=True,
        force_prefixed=True,
    )
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    os.replace(tmp, path)
    return path


def global_cache_path(name: str) -> Path:
    return unified_root() / name
