"""Per-sample analysis-version management.

Each sample owns one or more analysis versions, stored as
    tertiary_output/{sample_id}/analyses/{version}/analysis.json
plus the sidecar TSVs (pheno_score.tsv / exomiser_results.tsv /
lirical_results.tsv) and the analysis_files/ run directory used by the
Exomiser/LIRICAL worker.

The active version (the one auto-loaded next time the sample opens) is
recorded in sample_metadata.json under 'active_analysis'.

'default' is reserved as the first version created for any sample and
cannot be deleted or renamed; reviewers can still edit its HPO/panel
list and re-run analysis on top of it.
"""
from __future__ import annotations

import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

from ..config import TERTIARY_OUTPUT_ROOT

VERSION_NAME_RE = re.compile(r"^[A-Za-z0-9_\-]{1,32}$")
RESERVED_NAMES = {"default"}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def validate_name(name: str) -> None:
    if not VERSION_NAME_RE.match(name or ""):
        raise ValueError(
            "version name must match [A-Za-z0-9_-]{1,32}"
        )


def sample_dir(sample_id: str) -> Path:
    return TERTIARY_OUTPUT_ROOT / sample_id


def analyses_dir(sample_id: str) -> Path:
    return sample_dir(sample_id) / "analyses"


def version_dir(sample_id: str, version: str) -> Path:
    return analyses_dir(sample_id) / version


def active_version(sample_id: str) -> str | None:
    """Return the active analysis name from sample_metadata.json, or None.

    Falls back to 'default' if it exists, then the first available
    version, else None (un-migrated sample with no analyses dir).
    """
    meta_path = sample_dir(sample_id) / "sample_metadata.json"
    name: str | None = None
    if meta_path.exists():
        try:
            data = json.loads(meta_path.read_text(encoding="utf-8"))
            name = data.get("active_analysis") if isinstance(data, dict) else None
        except (json.JSONDecodeError, OSError):
            pass
    if name and (version_dir(sample_id, name) / "analysis.json").exists():
        return name
    if (version_dir(sample_id, "default") / "analysis.json").exists():
        return "default"
    versions = list_versions(sample_id)
    return versions[0]["name"] if versions else None


def active_version_dir(sample_id: str) -> Path:
    """Where sidecar TSVs (pheno/exomiser/lirical) for the active version go.

    Falls back to the sample root for pre-migration samples (no
    analyses/ dir yet) so the workers keep landing files somewhere
    valid. The loader's matching fallback path picks them up.
    """
    name = active_version(sample_id)
    if name is None:
        return sample_dir(sample_id)
    return version_dir(sample_id, name)


def set_active(sample_id: str, version: str) -> None:
    """Persist `active_analysis` on sample_metadata.json. Validates name."""
    if version not in {v["name"] for v in list_versions(sample_id)}:
        raise ValueError(f"unknown analysis version: {version}")
    meta_path = sample_dir(sample_id) / "sample_metadata.json"
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if not isinstance(meta, dict):
            meta = {}
    except (json.JSONDecodeError, OSError):
        meta = {}
    meta["active_analysis"] = version
    meta["updated_at"] = _now()
    meta_path.write_text(
        json.dumps(meta, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def list_versions(sample_id: str) -> list[dict]:
    """Return [{name, created_at, updated_at, n_hpo, n_panels, note}].

    Sorted by name with 'default' first when present.
    """
    root = analyses_dir(sample_id)
    if not root.is_dir():
        return []
    out: list[dict] = []
    for sub in sorted(root.iterdir()):
        if not sub.is_dir():
            continue
        ap = sub / "analysis.json"
        if not ap.exists():
            continue
        try:
            data = json.loads(ap.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            data = {}
        out.append({
            "name": sub.name,
            "created_at": data.get("created_at"),
            "updated_at": data.get("updated_at"),
            "n_hpo":    len(data.get("hpo") or []),
            "n_panels": len(data.get("selected_panels") or []),
            "note":     data.get("note", ""),
        })
    # 'default' first; then alphabetical.
    out.sort(key=lambda r: (0 if r["name"] == "default" else 1, r["name"]))
    return out


def read_version(sample_id: str, version: str) -> dict | None:
    p = version_dir(sample_id, version) / "analysis.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def write_version(
    sample_id: str,
    version: str,
    hpo: list,
    panels: list,
    note: str = "",
) -> dict:
    """Create or overwrite the analysis.json for a version.

    For overwrites, 'created_at' is preserved; 'updated_at' is set to now.
    """
    validate_name(version)
    vdir = version_dir(sample_id, version)
    vdir.mkdir(parents=True, exist_ok=True)
    p = vdir / "analysis.json"
    existing: dict = {}
    if p.exists():
        try:
            existing = json.loads(p.read_text(encoding="utf-8")) or {}
        except json.JSONDecodeError:
            existing = {}
    now = _now()
    payload = {
        "hpo": hpo or [],
        "selected_panels": panels or [],
        "note": note or "",
        "created_at": existing.get("created_at") or now,
        "updated_at": now,
    }
    p.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    # Side effect: keep pheno_score.tsv in sync with this version's
    # HPO/panels. Computing it here means every caller that touches a
    # version (register, edit, copy/rename) gets a consistent
    # gene-score sidecar without having to remember the second step.
    # Empty HPO+panels → wipe any stale pheno_score.tsv from a
    # previous edit instead of leaving misleading content behind.
    from . import phenotype_scorer
    pheno_tsv = vdir / "pheno_score.tsv"
    if (hpo or panels):
        scores = phenotype_scorer.compute_pheno_score(hpo or [], panels or [])
        phenotype_scorer.write_pheno_table(sample_id, scores, target_dir=vdir)
    elif pheno_tsv.exists():
        try:
            pheno_tsv.unlink()
        except OSError:
            pass

    return payload


def delete_version(sample_id: str, version: str) -> bool:
    """Remove a non-reserved version directory entirely.

    Returns True if removed, False if missing. Raises ValueError on
    reserved names.
    """
    if version in RESERVED_NAMES:
        raise ValueError(f"cannot delete reserved version '{version}'")
    vdir = version_dir(sample_id, version)
    if not vdir.is_dir():
        return False
    shutil.rmtree(vdir)
    return True


def clear_sidecars(sample_id: str, version: str) -> None:
    """Wipe pheno/exomiser/lirical sidecars + analysis_files/ in a version.

    Used before a re-run overwrites the version, so a partial failure
    can't leave stale results blended with new ones.
    """
    vdir = version_dir(sample_id, version)
    for fname in ("pheno_score.tsv", "exomiser_results.tsv", "lirical_results.tsv"):
        p = vdir / fname
        if p.exists():
            p.unlink()
    af = vdir / "analysis_files"
    if af.exists():
        shutil.rmtree(af)
