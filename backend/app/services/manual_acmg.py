"""Manual ACMG evidence, reusable coordinate assertions, and case observations.

Manual assertions are matched only by genome build + normalized
``chr-pos-ref-alt``.  The reusable latest assertion is deliberately separate
from per-sample snapshots in ``sample_metadata.json`` and from the observation
registry of currently active Causative/Other cases.
"""
from __future__ import annotations

import json
import math
import re
import sqlite3
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from ..config import MANUAL_ACMG_DB


ACMG_2015_URL = "https://pubmed.ncbi.nlm.nih.gov/25741868/"
TAVTIGIAN_URL = "https://pubmed.ncbi.nlm.nih.gov/32720330/"
VUS_SUBCLASS_URL = "https://doi.org/10.1016/j.gim.2024.101276"
CRITERIA_ORDER = (
    "PVS1",
    "PS1", "PS2", "PS3", "PS4",
    "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5",
    "BA1",
    "BS1", "BS2", "BS3", "BS4",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
)
EVIDENCE_GROUPS = (
    {
        "value": "population",
        "label": "Population data",
        "criteria": ("BA1", "BS1", "BS2", "PM2", "PS4"),
    },
    {
        "value": "computational_predictive",
        "label": "Computational and predictive data",
        "criteria": ("BP4", "BP1", "BP7", "BP3", "PP3", "PM5", "PM4", "PS1", "PVS1"),
    },
    {
        "value": "functional",
        "label": "Functional data",
        "criteria": ("BS3", "PP2", "PM1", "PS3"),
    },
    {
        "value": "segregation",
        "label": "Segregation data",
        "criteria": ("BS4", "PP1"),
    },
    {
        "value": "de_novo",
        "label": "De novo data",
        "criteria": ("PM6", "PS2"),
    },
    {
        "value": "allelic",
        "label": "Allelic data",
        "criteria": ("BP2", "PM3"),
    },
    {
        "value": "other_database",
        "label": "Other database",
        "criteria": ("BP6", "PP5"),
    },
    {
        "value": "other_data",
        "label": "Other data",
        "criteria": ("BP5", "PP4"),
    },
)
_CRITERION_EVIDENCE_GROUP = {
    code: (group["value"], group["label"], group_order, criterion_order)
    for group_order, group in enumerate(EVIDENCE_GROUPS)
    for criterion_order, code in enumerate(group["criteria"])
}
CASE_SCOPED_CRITERIA = frozenset({
    "PS2", "PM3", "PM6", "PP1", "PP4", "BS4", "BP2", "BP5",
})

_DESCRIPTIONS = {
    "PVS1": "Null variant (nonsense, frameshift, canonical ±1/2 splice, initiation codon, or exon/gene deletion) in a gene where loss of function is an established disease mechanism.",
    "PS1": "Same amino-acid change as an established pathogenic variant, regardless of the nucleotide substitution.",
    "PS2": "De novo variant in the patient, with maternity and paternity confirmed, in a compatible disorder.",
    "PS3": "Well-established functional studies show a damaging effect on the gene or gene product.",
    "PS4": "Variant prevalence is significantly increased in affected individuals compared with controls.",
    "PM1": "Located in a mutational hot spot or critical, well-established functional domain without benign variation.",
    "PM2": "Absent from controls, or extremely rare when appropriate for a recessive disorder, in population databases.",
    "PM3": "For a recessive disorder, detected in trans with a pathogenic variant.",
    "PM4": "Protein length changes caused by an in-frame insertion/deletion outside a repetitive region, or by stop loss.",
    "PM5": "Novel missense change at a residue where a different missense change is established as pathogenic.",
    "PM6": "Assumed de novo, but maternity and paternity have not both been confirmed.",
    "PP1": "Cosegregation with disease in multiple affected relatives in a gene definitively known to cause the disorder.",
    "PP2": "Missense variant in a gene with few benign missense variants and where missense is a common disease mechanism.",
    "PP3": "Multiple computational lines support a deleterious effect on the gene or gene product.",
    "PP4": "The phenotype or family history is highly specific for a disorder with a single genetic cause.",
    "PP5": "A reputable source reports pathogenicity, but the underlying evidence is unavailable for independent evaluation.",
    "BA1": "Allele frequency is greater than 5% in population databases (stand-alone benign evidence).",
    "BS1": "Allele frequency is greater than expected for the disorder.",
    "BS2": "Observed in a healthy adult despite an expectation of full penetrance at an early age.",
    "BS3": "Well-established functional studies show no damaging effect on protein function or splicing.",
    "BS4": "Lack of segregation in affected members of a family.",
    "BP1": "Missense variant in a gene where disease is caused primarily by truncating variants.",
    "BP2": "Observed in trans with a pathogenic variant for a fully penetrant dominant disorder, or in cis with a pathogenic variant.",
    "BP3": "In-frame insertion/deletion in a repetitive region without a known function.",
    "BP4": "Multiple computational lines suggest no impact on the gene or gene product.",
    "BP5": "Variant found in a case with an alternate molecular basis for disease.",
    "BP6": "A reputable source reports benignity, but the underlying evidence is unavailable for independent evaluation.",
    "BP7": "Synonymous variant at a nonconserved nucleotide with no predicted splice impact.",
}

_DEFAULT_STRENGTH = {
    "PVS1": "very_strong",
    **{code: "strong" for code in ("PS1", "PS2", "PS3", "PS4")},
    **{code: "moderate" for code in ("PM1", "PM2", "PM3", "PM4", "PM5", "PM6")},
    **{code: "supporting" for code in ("PP1", "PP2", "PP3", "PP4", "PP5")},
    "BA1": "stand_alone",
    **{code: "strong" for code in ("BS1", "BS2", "BS3", "BS4")},
    **{code: "supporting" for code in ("BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7")},
}

_GUIDANCE = {
    "PVS1": [
        ("ClinGen PVS1 recommendations", "https://pubmed.ncbi.nlm.nih.gov/30192042/"),
        ("ClinGen SVI splicing recommendations", "https://pubmed.ncbi.nlm.nih.gov/37352859/"),
    ],
    "PS1": [("ClinGen SVI splicing recommendations", "https://pubmed.ncbi.nlm.nih.gov/37352859/")],
    "PS2": [("ClinGen SVI de novo recommendation", "https://www.clinicalgenome.org/site/assets/files/3461/svi_proposal_for_de_novo_criteria_v1_1.pdf")],
    "PM6": [("ClinGen SVI de novo recommendation", "https://www.clinicalgenome.org/site/assets/files/3461/svi_proposal_for_de_novo_criteria_v1_1.pdf")],
    "PS3": [("ClinGen SVI functional evidence recommendations", "https://pubmed.ncbi.nlm.nih.gov/31892348/")],
    "BS3": [("ClinGen SVI functional evidence recommendations", "https://pubmed.ncbi.nlm.nih.gov/31892348/")],
    "PM2": [("ClinGen SVI PM2 recommendation", "https://clinicalgenome.org/site/assets/files/5182/pm2_-_svi_recommendation_-_approved_sept2020.pdf")],
    "PM3": [("ClinGen SVI PM3 recommendation", "https://www.clinicalgenome.org/docs/pm3-recommendation-for-in-trans-criterion-pm3-version-1.0/")],
    "PP1": [("ClinGen SVI segregation recommendations", "https://pubmed.ncbi.nlm.nih.gov/38103548/")],
    "BS4": [("ClinGen SVI segregation recommendations", "https://pubmed.ncbi.nlm.nih.gov/38103548/")],
    "PP3": [("ClinGen SVI computational evidence recommendations", "https://pubmed.ncbi.nlm.nih.gov/36413997/")],
    "BP4": [("ClinGen SVI computational evidence recommendations", "https://pubmed.ncbi.nlm.nih.gov/36413997/")],
    "BP7": [("ClinGen SVI splicing recommendations", "https://pubmed.ncbi.nlm.nih.gov/37352859/")],
    "PP4": [("ClinGen SVI phenotype specificity recommendations", "https://pubmed.ncbi.nlm.nih.gov/38103548/")],
    "BA1": [("ClinGen SVI BA1 recommendations", "https://pubmed.ncbi.nlm.nih.gov/30311383/")],
    "PP5": [("ClinGen recommendation to retire PP5/BP6", "https://pubmed.ncbi.nlm.nih.gov/29543229/")],
    "BP6": [("ClinGen recommendation to retire PP5/BP6", "https://pubmed.ncbi.nlm.nih.gov/29543229/")],
}

_POINTS = {
    "supporting": 1,
    "moderate": 2,
    "strong": 4,
    "very_strong": 8,
    "stand_alone": 8,
}
_STRENGTH_LABELS = {
    "supporting": "Supporting",
    "moderate": "Moderate",
    "strong": "Strong",
    "very_strong": "Very strong",
    "stand_alone": "Stand-alone",
}
_DB_LOCK = threading.RLock()
_VARIANT_RE = re.compile(r"^(?:chr)?([^-]+)-(\d+)-([^-]+)-([^-]+)$", re.IGNORECASE)


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def normalize_build(value: Any) -> str:
    raw = str(value or "hg38").strip().lower().replace("_", "")
    if raw in {"38", "grch38", "hg38"}:
        return "hg38"
    return raw or "hg38"


def normalize_variant_id(value: Any) -> str:
    match = _VARIANT_RE.match(str(value or "").strip())
    if not match:
        raise ValueError("variant 必須是 chr-pos-ref-alt")
    chrom, pos, ref, alt = match.groups()
    chrom = chrom.upper()
    if chrom == "MT":
        chrom = "M"
    elif chrom not in {"X", "Y", "M"}:
        try:
            chrom = str(int(chrom))
        except ValueError:
            pass
    return f"chr{chrom}-{int(pos)}-{ref.upper()}-{alt.upper()}"


def is_snv_variant_id(value: Any) -> bool:
    try:
        return not normalize_variant_id(value).startswith("chrM-")
    except ValueError:
        return False


def catalog() -> dict[str, Any]:
    criteria = []
    for code in CRITERIA_ORDER:
        direction = "benign" if code.startswith(("BA", "BS", "BP")) else "pathogenic"
        group_value, group_label, group_order, criterion_order = (
            _CRITERION_EVIDENCE_GROUP[code]
        )
        refs = [{"title": "ACMG/AMP 2015", "url": ACMG_2015_URL}]
        refs.extend({"title": title, "url": url} for title, url in _GUIDANCE.get(code, []))
        criteria.append({
            "code": code,
            "direction": direction,
            "evidence_group": group_value,
            "evidence_group_label": group_label,
            "evidence_group_order": group_order,
            "evidence_order": criterion_order,
            "default_strength": _DEFAULT_STRENGTH[code],
            "description": _DESCRIPTIONS[code],
            "scope": "case" if code in CASE_SCOPED_CRITERIA else "global",
            "scope_note": (
                "Case/family-specific evidence: saved for this sample and revision, but not automatically applied to other samples."
                if code in CASE_SCOPED_CRITERIA else ""
            ),
            "deprecated_warning": (
                "ClinGen SVI recommends retiring this criterion; it remains usable here for compatibility with GeneBe."
                if code in {"PP5", "BP6"} else ""
            ),
            "references": refs,
        })
    return {
        "criteria": criteria,
        "evidence_groups": [
            {"value": group["value"], "label": group["label"]}
            for group in EVIDENCE_GROUPS
        ],
        "strengths": [
            {"value": key, "label": label, "points": points}
            for key, label in _STRENGTH_LABELS.items()
            for points in [_POINTS[key]]
        ],
        "scoring_reference": {"title": "Tavtigian et al. point system", "url": TAVTIGIAN_URL},
        "vus_subclasses": [
            {"value": "VUS-low", "min_points": 0, "max_points": 1},
            {"value": "VUS-mid", "min_points": 2, "max_points": 3},
            {"value": "VUS-high", "min_points": 4, "max_points": 5},
        ],
        "vus_subclass_reference": {
            "title": "Eldomery et al. VUS sub-tiering",
            "url": VUS_SUBCLASS_URL,
        },
    }


def parse_criteria_text(value: Any) -> tuple[dict[str, dict[str, Any]], list[str]]:
    """Convert GeneBe/pipeline criterion tokens into structured evidence."""
    parsed: dict[str, dict[str, Any]] = {}
    unknown: list[str] = []
    for token in re.split(r"[,|;\s]+", str(value or "").strip()):
        if not token:
            continue
        match = re.match(
            r"^(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])"
            r"(?:[_-](supporting|moderate|strong|very[_-]?strong|stand[_-]?alone))?$",
            token,
            re.IGNORECASE,
        )
        if not match:
            unknown.append(token)
            continue
        code = match.group(1).upper()
        strength = (match.group(2) or _DEFAULT_STRENGTH[code]).lower().replace("-", "_")
        strength = strength.replace("verystrong", "very_strong").replace("standalone", "stand_alone")
        parsed[code] = {"enabled": True, "strength": strength}
    return parsed, unknown


def normalize_criteria(value: Any) -> dict[str, dict[str, Any]]:
    if isinstance(value, str):
        return parse_criteria_text(value)[0]
    if not isinstance(value, dict):
        raise ValueError("criteria 必須是物件")
    out: dict[str, dict[str, Any]] = {}
    for code, raw in value.items():
        code = str(code or "").upper()
        if code not in _DESCRIPTIONS:
            raise ValueError(f"未知 ACMG criterion: {code}")
        if isinstance(raw, bool):
            raw = {"enabled": raw}
        if not isinstance(raw, dict) or not raw.get("enabled", False):
            continue
        strength = str(raw.get("strength") or _DEFAULT_STRENGTH[code]).lower().replace("-", "_")
        if strength not in _POINTS:
            raise ValueError(f"{code} strength 無效: {strength}")
        if strength == "stand_alone" and code != "BA1":
            raise ValueError(f"{code} 不支援 Stand-alone strength")
        out[code] = {"enabled": True, "strength": strength}
    return out


def criteria_text(criteria: dict[str, dict[str, Any]]) -> str:
    tokens: list[str] = []
    for code in CRITERIA_ORDER:
        evidence = criteria.get(code)
        if not evidence or not evidence.get("enabled"):
            continue
        strength = str(evidence.get("strength") or _DEFAULT_STRENGTH[code])
        if strength == _DEFAULT_STRENGTH[code]:
            tokens.append(code)
        else:
            tokens.append(f"{code}_{_STRENGTH_LABELS[strength].replace(' ', '_')}")
    return ",".join(tokens)


def vus_subclass(classification: Any, points: Any) -> str:
    """Return the display-only VUS point sub-tier without changing ACMG class."""
    normalized = re.sub(
        r"[_-]+", " ", str(classification or "").strip().lower()
    )
    explicit = {
        "vus low": "VUS-low",
        "vus mid": "VUS-mid",
        "vus high": "VUS-high",
    }
    if normalized in explicit:
        return explicit[normalized]
    if normalized not in {
        "vus",
        "uncertain significance",
        "variant of uncertain significance",
    }:
        return ""
    try:
        score = float(points)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(score) or not score.is_integer():
        return ""
    point = int(score)
    if 0 <= point <= 1:
        return "VUS-low"
    if 2 <= point <= 3:
        return "VUS-mid"
    if 4 <= point <= 5:
        return "VUS-high"
    return ""


def calculate(criteria_value: Any) -> dict[str, Any]:
    criteria = normalize_criteria(criteria_value)
    score = 0
    has_pathogenic = False
    has_benign = False
    for code, evidence in criteria.items():
        points = _POINTS[evidence["strength"]]
        if code.startswith(("BA", "BS", "BP")):
            score -= points
            has_benign = True
        else:
            score += points
            has_pathogenic = True
    if score >= 10:
        classification = "Pathogenic"
    elif score >= 6:
        classification = "Likely pathogenic"
    elif score >= 0:
        classification = "Uncertain significance"
    elif score >= -6:
        classification = "Likely benign"
    else:
        classification = "Benign"
    return {
        "criteria": criteria,
        "criteria_text": criteria_text(criteria),
        "score": score,
        "classification": classification,
        "vus_subclass": vus_subclass(classification, score),
        "conflicting_evidence": has_pathogenic and has_benign,
    }


def reusable_criteria(criteria_value: Any) -> dict[str, dict[str, Any]]:
    criteria = normalize_criteria(criteria_value)
    return {
        code: evidence
        for code, evidence in criteria.items()
        if code not in CASE_SCOPED_CRITERIA
    }


def acmg_to_variant_score(points: Any) -> int | None:
    try:
        value = float(points)
    except (TypeError, ValueError):
        return None
    value = max(-10.0, min(10.0, value))
    return int(round((value + 10.0) / 20.0 * 100.0))


def _connect(path: Path | None = None) -> sqlite3.Connection:
    db_path = Path(path or MANUAL_ACMG_DB)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path, timeout=5.0)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=5000")
    conn.execute("PRAGMA foreign_keys=ON")
    _ensure_schema(conn)
    return conn


def _ensure_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS manual_revisions (
          id INTEGER PRIMARY KEY AUTOINCREMENT,
          genome_build TEXT NOT NULL,
          variant_id TEXT NOT NULL,
          criteria_json TEXT NOT NULL,
          criteria_text TEXT NOT NULL,
          acmg_score REAL NOT NULL,
          acmg_classification TEXT NOT NULL,
          reusable_criteria_json TEXT,
          reusable_criteria_text TEXT,
          reusable_acmg_score REAL,
          reusable_acmg_classification TEXT,
          reviewer_user_id INTEGER,
          reviewer_username TEXT NOT NULL,
          source_sample_id TEXT NOT NULL,
          created_at TEXT NOT NULL
        );
        CREATE INDEX IF NOT EXISTS idx_manual_revision_key
          ON manual_revisions(genome_build, variant_id, id DESC);
        CREATE TABLE IF NOT EXISTS manual_current (
          genome_build TEXT NOT NULL,
          variant_id TEXT NOT NULL,
          revision_id INTEGER NOT NULL REFERENCES manual_revisions(id),
          PRIMARY KEY (genome_build, variant_id)
        );
        CREATE TABLE IF NOT EXISTS observations (
          genome_build TEXT NOT NULL,
          variant_id TEXT NOT NULL,
          sample_id TEXT NOT NULL,
          status TEXT NOT NULL CHECK(status IN ('1', '2')),
          reviewer_user_id INTEGER,
          reviewer_username TEXT NOT NULL,
          updated_at TEXT NOT NULL,
          PRIMARY KEY (genome_build, variant_id, sample_id)
        );
        CREATE INDEX IF NOT EXISTS idx_observation_lookup
          ON observations(genome_build, variant_id, sample_id);
        """
    )
    columns = {
        row["name"] if isinstance(row, sqlite3.Row) else row[1]
        for row in conn.execute("PRAGMA table_info(manual_revisions)").fetchall()
    }
    additions = {
        "reusable_criteria_json": "TEXT",
        "reusable_criteria_text": "TEXT",
        "reusable_acmg_score": "REAL",
        "reusable_acmg_classification": "TEXT",
    }
    for name, kind in additions.items():
        if name not in columns:
            conn.execute(f"ALTER TABLE manual_revisions ADD COLUMN {name} {kind}")


def _row_to_assertion(row: sqlite3.Row | None) -> dict[str, Any] | None:
    if row is None:
        return None
    full_criteria = json.loads(row["criteria_json"])
    reusable_json = row["reusable_criteria_json"]
    if reusable_json:
        reusable = json.loads(reusable_json)
        reusable_result = {
            "criteria": reusable,
            "criteria_text": row["reusable_criteria_text"] or "",
            "score": row["reusable_acmg_score"],
            "classification": row["reusable_acmg_classification"] or "",
        }
    else:
        reusable_result = calculate(reusable_criteria(full_criteria))
    score = row["acmg_score"]
    reusable_score = reusable_result["score"]
    return {
        "revision_id": row["id"],
        "genome_build": row["genome_build"],
        "variant_id": row["variant_id"],
        "criteria": full_criteria,
        "criteria_text": row["criteria_text"],
        "score": int(score) if float(score).is_integer() else score,
        "classification": row["acmg_classification"],
        "vus_subclass": vus_subclass(row["acmg_classification"], score),
        "reusable_criteria": reusable_result["criteria"],
        "reusable_criteria_text": reusable_result["criteria_text"],
        "reusable_score": (
            int(reusable_score)
            if float(reusable_score).is_integer() else reusable_score
        ),
        "reusable_classification": reusable_result["classification"],
        "reusable_vus_subclass": vus_subclass(
            reusable_result["classification"], reusable_score
        ),
        "reviewer_user_id": row["reviewer_user_id"],
        "reviewer_username": row["reviewer_username"],
        "source_sample_id": row["source_sample_id"],
        "created_at": row["created_at"],
    }


def save_assertion(
    genome_build: Any,
    variant_id: Any,
    criteria_value: Any,
    *,
    reviewer_user_id: int | None,
    reviewer_username: str,
    source_sample_id: str,
    path: Path | None = None,
) -> dict[str, Any]:
    build = normalize_build(genome_build)
    vid = normalize_variant_id(variant_id)
    result = calculate(criteria_value)
    reusable_result = calculate(reusable_criteria(result["criteria"]))
    created_at = now_iso()
    with _DB_LOCK, _connect(path) as conn:
        cursor = conn.execute(
            """
            INSERT INTO manual_revisions (
              genome_build, variant_id, criteria_json, criteria_text,
              acmg_score, acmg_classification,
              reusable_criteria_json, reusable_criteria_text,
              reusable_acmg_score, reusable_acmg_classification,
              reviewer_user_id,
              reviewer_username, source_sample_id, created_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                build, vid,
                json.dumps(result["criteria"], ensure_ascii=False, sort_keys=True),
                result["criteria_text"], result["score"], result["classification"],
                json.dumps(
                    reusable_result["criteria"],
                    ensure_ascii=False,
                    sort_keys=True,
                ),
                reusable_result["criteria_text"],
                reusable_result["score"],
                reusable_result["classification"],
                reviewer_user_id, str(reviewer_username or ""),
                str(source_sample_id or ""), created_at,
            ),
        )
        revision_id = int(cursor.lastrowid)
        conn.execute(
            """
            INSERT INTO manual_current(genome_build, variant_id, revision_id)
            VALUES (?, ?, ?)
            ON CONFLICT(genome_build, variant_id)
            DO UPDATE SET revision_id=excluded.revision_id
            """,
            (build, vid, revision_id),
        )
        row = conn.execute(
            "SELECT * FROM manual_revisions WHERE id=?", (revision_id,)
        ).fetchone()
    assertion = _row_to_assertion(row)
    assert assertion is not None
    assertion["conflicting_evidence"] = result["conflicting_evidence"]
    return assertion


def current_assertion(
    genome_build: Any, variant_id: Any, *, path: Path | None = None
) -> dict[str, Any] | None:
    build = normalize_build(genome_build)
    vid = normalize_variant_id(variant_id)
    with _connect(path) as conn:
        row = conn.execute(
            """
            SELECT r.* FROM manual_current c
            JOIN manual_revisions r ON r.id=c.revision_id
            WHERE c.genome_build=? AND c.variant_id=?
            """,
            (build, vid),
        ).fetchone()
    return _row_to_assertion(row)


def bulk_current(
    genome_build: Any, variant_ids: Iterable[Any], *, path: Path | None = None
) -> dict[str, dict[str, Any]]:
    build = normalize_build(genome_build)
    normalized = []
    for value in variant_ids:
        try:
            normalized.append(normalize_variant_id(value))
        except ValueError:
            continue
    if not normalized:
        return {}
    out: dict[str, dict[str, Any]] = {}
    with _connect(path) as conn:
        for start in range(0, len(normalized), 800):
            chunk = normalized[start:start + 800]
            placeholders = ",".join("?" for _ in chunk)
            rows = conn.execute(
                f"""
                SELECT r.* FROM manual_current c
                JOIN manual_revisions r ON r.id=c.revision_id
                WHERE c.genome_build=? AND c.variant_id IN ({placeholders})
                """,
                (build, *chunk),
            ).fetchall()
            for row in rows:
                assertion = _row_to_assertion(row)
                if assertion:
                    out[assertion["variant_id"]] = assertion
    return out


def sync_observations(
    genome_build: Any,
    sample_id: str,
    statuses: dict[str, Any],
    *,
    reviewer_user_id: int | None,
    reviewer_username: str,
    updated_at: str | None = None,
    status_audit: dict[str, Any] | None = None,
    path: Path | None = None,
) -> None:
    """Replace one sample's active SNV Causative/Other observations."""
    build = normalize_build(genome_build)
    active: dict[str, str] = {}
    for variant_id, raw_status in (statuses or {}).items():
        status = str(raw_status or "").strip()
        if status not in {"1", "2"} or not is_snv_variant_id(variant_id):
            continue
        active[normalize_variant_id(variant_id)] = status
    normalized_audit: dict[str, Any] = {}
    for variant_id, audit in (status_audit or {}).items():
        try:
            normalized_audit[normalize_variant_id(variant_id)] = audit
        except ValueError:
            continue
    timestamp = updated_at or now_iso()
    with _DB_LOCK, _connect(path) as conn:
        conn.execute(
            "DELETE FROM observations WHERE genome_build=? AND sample_id=?",
            (build, str(sample_id)),
        )
        rows = []
        for vid, status in active.items():
            audit = normalized_audit.get(vid) or {}
            rows.append((
                build,
                vid,
                str(sample_id),
                status,
                audit.get("reviewer_user_id", reviewer_user_id),
                str(audit.get("reviewer_username", reviewer_username) or ""),
                str(audit.get("updated_at", timestamp) or timestamp),
            ))
        conn.executemany(
            """
            INSERT INTO observations (
              genome_build, variant_id, sample_id, status,
              reviewer_user_id, reviewer_username, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            rows,
        )


def backfill_observations_from_samples() -> int:
    """Rebuild active observations from existing sample metadata.

    This reads only small JSON sidecars.  It intentionally does not retain
    cancelled/Candidate statuses and is safe to rerun at startup.
    """
    from . import sample_layout

    count = 0
    for sample_id in sample_layout.iter_sample_ids():
        path = sample_layout.state_file(sample_id, "sample_metadata.json")
        try:
            meta = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if not isinstance(meta, dict):
            continue
        sync_observations(
            meta.get("genome_build") or "hg38",
            sample_id,
            meta.get("status") or {},
            reviewer_user_id=meta.get("last_reviewer_user_id"),
            reviewer_username=meta.get("last_reviewer_username") or "",
            updated_at=meta.get("updated_at") or now_iso(),
            status_audit=meta.get("status_audit") or {},
        )
        count += 1
    return count


def bulk_observed_counts(
    genome_build: Any,
    variant_ids: Iterable[Any],
    *,
    exclude_sample_id: str = "",
    path: Path | None = None,
) -> dict[str, int]:
    build = normalize_build(genome_build)
    normalized = []
    for value in variant_ids:
        try:
            normalized.append(normalize_variant_id(value))
        except ValueError:
            continue
    if not normalized:
        return {}
    out: dict[str, int] = {}
    with _connect(path) as conn:
        for start in range(0, len(normalized), 800):
            chunk = normalized[start:start + 800]
            placeholders = ",".join("?" for _ in chunk)
            rows = conn.execute(
                f"""
                SELECT variant_id, COUNT(*) AS n
                FROM observations
                WHERE genome_build=? AND sample_id<>?
                  AND variant_id IN ({placeholders})
                GROUP BY variant_id
                """,
                (build, str(exclude_sample_id or ""), *chunk),
            ).fetchall()
            out.update({row["variant_id"]: int(row["n"]) for row in rows})
    return out


def observed_cases(
    genome_build: Any,
    variant_id: Any,
    *,
    exclude_sample_id: str = "",
    path: Path | None = None,
) -> list[dict[str, Any]]:
    build = normalize_build(genome_build)
    vid = normalize_variant_id(variant_id)
    with _connect(path) as conn:
        rows = conn.execute(
            """
            SELECT sample_id, status, reviewer_user_id, reviewer_username, updated_at
            FROM observations
            WHERE genome_build=? AND variant_id=? AND sample_id<>?
            ORDER BY updated_at DESC, sample_id
            """,
            (build, vid, str(exclude_sample_id or "")),
        ).fetchall()
    return [
        {
            "sample_id": row["sample_id"],
            "status": row["status"],
            "status_label": "Causative" if row["status"] == "1" else "Other",
            "reviewer_user_id": row["reviewer_user_id"],
            "reviewer_username": row["reviewer_username"],
            "updated_at": row["updated_at"],
        }
        for row in rows
    ]


def remove_sample_observations(
    sample_id: str, *, path: Path | None = None
) -> None:
    with _DB_LOCK, _connect(path) as conn:
        conn.execute(
            "DELETE FROM observations WHERE sample_id=?",
            (str(sample_id),),
        )
