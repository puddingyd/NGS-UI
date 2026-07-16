#!/usr/bin/env python3
"""Canary/bulk migration into the unified tertiary-output layout.

The default is a read-only dry run. ``--apply`` copies 00-07 and reviewer
state, rebuilds the SNV sparse overlay/review/index, verifies the result, and
only then writes ``08_postprocessing/layout.json``. Old data is never deleted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app import config  # noqa: E402
from app.services import sample_layout, snv_gene_index, snv_overlay, snv_review  # noqa: E402


ALWAYS_REBUILD_OR_SKIP = {
    "snv_indel.annotated.tsv",
    "snv_indel.review.tsv",
    "snv_indel.review.tsv.source.json",
    "snv_gene_index.sqlite",
    "snv_annotations.sqlite",
    "case_summary.json",
}

PIPELINE_COPY_PATTERNS = {
    "cnv.annotated.tsv": ("06_cnv_sv", "*.cnv.annotated.tsv"),
    "sv.annotated.tsv": ("06_cnv_sv", "*.sv.annotated.tsv"),
    "mito.annotated.tsv": ("04_mito", "*.mito.tsv"),
    "str.tsv": ("05_str", "*.str.tsv"),
    "str.annotated.tsv": ("05_str", "*.str.tsv"),
    "pgx.tsv": ("07_pgx", "*.pgx.tsv"),
    "pharmcat.report.json": ("07_pgx", "*.pharmcat.report.json"),
    "outside_calls.tsv": ("07_pgx", "*.outside_calls.tsv"),
    "stellarpgx.tsv": ("07_pgx", "*.stellarpgx.tsv"),
    "optitype.tsv": ("07_pgx", "*.optitype.tsv"),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _copy_entry(source: Path, destination: Path) -> None:
    if source.is_dir():
        shutil.copytree(source, destination, dirs_exist_ok=True, symlinks=True)
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination, follow_symlinks=False)


def _same_content(left: Path, right: Path) -> bool:
    try:
        if left.stat().st_size != right.stat().st_size:
            return False
    except OSError:
        return False
    return _sha256(left) == _sha256(right)


def _matching_pipeline_copy(entry: Path, target_sample: Path) -> Path | None:
    spec = PIPELINE_COPY_PATTERNS.get(entry.name)
    if not spec or not entry.is_file():
        return None
    subdir, pattern = spec
    hits = sorted(path for path in (target_sample / subdir).glob(pattern) if path.is_file())
    if len(hits) != 1:
        return None
    return hits[0] if _same_content(entry, hits[0]) else None


def _read_json(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def _source_id(sample_id: str, legacy_state: Path) -> str:
    info = _read_json(legacy_state / "pipeline_source.json")
    return str(info.get("source_sample_id") or sample_id).strip() or sample_id


def _find_raw(sample_root: Path, source_id: str) -> Path:
    exact = sample_root / "03_acmg" / f"{source_id}.snv_indel.acmg.tsv"
    if exact.is_file():
        return exact
    hits = sorted((sample_root / "03_acmg").glob("*.snv_indel.acmg.tsv"))
    if hits:
        return hits[0]
    return exact


def _test_type(legacy_state: Path) -> str:
    value = str(_read_json(legacy_state / "sample_metadata.json").get("test_type") or "WES").upper()
    return value if value in {"WES", "WGS"} else "WES"


def _rollback(sample_id: str, target_root: Path, *, apply: bool) -> bool:
    marker = target_root / sample_id / sample_layout.POSTPROCESSING_DIRNAME / sample_layout.LAYOUT_MARKER_NAME
    if not marker.is_file():
        print(f"[{sample_id}] no active unified marker: {marker}")
        return True
    disabled = marker.with_name(
        f"{marker.name}.disabled.{time.strftime('%Y%m%dT%H%M%S')}"
    )
    print(f"[{sample_id}] rollback marker: {marker} -> {disabled}")
    if apply:
        os.replace(marker, disabled)
    return True


def migrate_one(
    sample_id: str,
    *,
    source_pipeline_root: Path,
    source_ui_root: Path,
    target_root: Path,
    apply: bool,
) -> bool:
    legacy_state = source_ui_root / sample_id
    source_id = _source_id(sample_id, legacy_state)
    source_pipeline = source_pipeline_root / sample_id
    if not source_pipeline.is_dir() and source_id != sample_id:
        source_pipeline = source_pipeline_root / source_id
    target_sample = target_root / sample_id
    post_dir = target_sample / sample_layout.POSTPROCESSING_DIRNAME
    marker_path = post_dir / sample_layout.LAYOUT_MARKER_NAME
    old_annotated = legacy_state / "snv_indel.annotated.tsv"

    print(f"[{sample_id}] source pipeline: {source_pipeline}")
    print(f"[{sample_id}] source UI state: {legacy_state}")
    print(f"[{sample_id}] target: {target_sample}")
    if marker_path.is_file():
        print(f"[{sample_id}] already active; skip (rollback first if re-migration is intended)")
        return True
    if not source_pipeline.is_dir() and not (target_sample / "03_acmg").is_dir():
        print(f"[{sample_id}] ERROR: pipeline source missing", file=sys.stderr)
        return False
    if not old_annotated.is_file():
        print(f"[{sample_id}] ERROR: legacy enriched SNV TSV missing: {old_annotated}", file=sys.stderr)
        return False
    if not apply:
        print(f"[{sample_id}] DRY-RUN: copy 00-07, copy state except pure duplicates, rebuild overlay/review/index, activate marker")
        return True

    target_sample.mkdir(parents=True, exist_ok=True)
    if source_pipeline.is_dir() and source_pipeline.resolve() != target_sample.resolve():
        for entry in sorted(source_pipeline.iterdir()):
            if entry.name == sample_layout.POSTPROCESSING_DIRNAME:
                continue
            _copy_entry(entry, target_sample / entry.name)

    post_dir.mkdir(parents=True, exist_ok=True)
    omitted_exact_copies: list[str] = []
    if legacy_state.is_dir() and legacy_state.resolve() != post_dir.resolve():
        for entry in sorted(legacy_state.iterdir()):
            if entry.name in ALWAYS_REBUILD_OR_SKIP:
                continue
            pipeline_copy = _matching_pipeline_copy(entry, target_sample)
            if pipeline_copy is not None:
                omitted_exact_copies.append(entry.name)
                print(f"[{sample_id}] omit exact copy: {entry.name} == {pipeline_copy}")
                continue
            _copy_entry(entry, post_dir / entry.name)

    raw_tsv = _find_raw(target_sample, source_id)
    if not raw_tsv.is_file():
        print(f"[{sample_id}] ERROR: target raw TSV missing: {raw_tsv}", file=sys.stderr)
        return False

    overlay_path = post_dir / snv_overlay.OVERLAY_NAME
    snv_overlay.build_overlay(raw_tsv, old_annotated, overlay_path)
    review_path = snv_review.ensure_review_tsv(
        raw_tsv,
        test_type=_test_type(legacy_state),
        output_dir=post_dir,
        overlay_path=overlay_path,
    )
    index_path = snv_gene_index.build_index(
        raw_tsv, post_dir / snv_gene_index.INDEX_NAME
    )

    checks = {
        "raw_exists": raw_tsv.is_file(),
        "overlay_current": snv_overlay.is_current(raw_tsv, overlay_path),
        "review_exists": review_path.is_file(),
        "index_current": snv_gene_index.is_current(raw_tsv, index_path),
    }
    old_meta = legacy_state / "sample_metadata.json"
    new_meta = post_dir / "sample_metadata.json"
    if old_meta.is_file():
        checks["metadata_checksum"] = new_meta.is_file() and _sha256(old_meta) == _sha256(new_meta)
    if not all(checks.values()):
        print(f"[{sample_id}] ERROR: verification failed: {checks}", file=sys.stderr)
        return False

    manifest = {
        "sample_id": sample_id,
        "source_sample_id": source_id,
        "source_pipeline_root": str(source_pipeline_root),
        "source_ui_root": str(source_ui_root),
        "target_root": str(target_root),
        "raw_tsv": str(raw_tsv),
        "raw_size": raw_tsv.stat().st_size,
        "overlay": str(overlay_path),
        "review": str(review_path),
        "gene_index": str(index_path),
        "verification": checks,
        "omitted_exact_pipeline_copies": omitted_exact_copies,
        "old_data_deleted": False,
        "migrated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    manifest_path = post_dir / "migration_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    # Configure the resolver used in this process, then activate last.
    config.PIPELINE_OUT_ROOT = target_root
    sample_layout.write_layout_marker(
        sample_id,
        source_id=source_id,
        raw_tsv=raw_tsv,
        migration=True,
    )
    print(f"[{sample_id}] OK: unified layout activated; legacy data retained")
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", action="append", default=[], help="can be repeated")
    parser.add_argument("--all", action="store_true", help="migrate every legacy sample")
    parser.add_argument("--apply", action="store_true", help="perform writes; default is dry-run")
    parser.add_argument("--rollback", action="store_true", help="disable marker only; keep both trees")
    parser.add_argument("--source-pipeline-root", type=Path, default=config.LEGACY_PIPELINE_OUT_ROOT)
    parser.add_argument("--source-ui-root", type=Path, default=config.LEGACY_TERTIARY_OUTPUT_ROOT)
    parser.add_argument("--target-root", type=Path, default=config.PIPELINE_OUT_ROOT)
    args = parser.parse_args()

    sample_ids = list(dict.fromkeys(str(value).strip() for value in args.sample if str(value).strip()))
    if args.all:
        if args.source_ui_root.is_dir():
            sample_ids.extend(
                child.name
                for child in sorted(args.source_ui_root.iterdir())
                if child.is_dir()
                and not child.name.startswith(("_", "."))
                and (child / "snv_indel.annotated.tsv").is_file()
            )
        sample_ids = list(dict.fromkeys(sample_ids))
    if not sample_ids:
        parser.error("provide --sample ID (recommended canary) or --all")

    print("APPLY" if args.apply else "DRY-RUN (no files will be changed)")
    ok = True
    for sample_id in sample_ids:
        if args.rollback:
            ok = _rollback(sample_id, args.target_root, apply=args.apply) and ok
        else:
            ok = migrate_one(
                sample_id,
                source_pipeline_root=args.source_pipeline_root,
                source_ui_root=args.source_ui_root,
                target_root=args.target_root,
                apply=args.apply,
            ) and ok
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
