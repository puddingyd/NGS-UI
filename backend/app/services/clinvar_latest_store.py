"""Compact local ClinVar snapshot used for weekly UI-only comparisons."""
from __future__ import annotations

import gzip
import json
import os
import re
import sqlite3
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable
from urllib.parse import unquote


SCHEMA_VERSION = "2"
_FILE_DATE_RE = re.compile(r"^##fileDate=(\d{4})-?(\d{2})-?(\d{2})")
_BASELINE_VARIANT_FIELDS = {
    "CLNSIG": "CLNSIG_old",
    "CLNSIGCONF": "CLNSIGCONF_old",
    "clinvar_stars": "clinvar_stars_old",
    "clinvar_dn": "clinvar_dn_old",
    "clinvar_variation_id": "clinvar_variation_id_old",
}


def normalize_chrom(value: str) -> str:
    chrom = str(value or "").strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if chrom.upper() in {"M", "MT"}:
        return "MT"
    return chrom.upper() if chrom.upper() in {"X", "Y"} else chrom


def normalize_variant(chrom: str, pos: str | int, ref: str, alt: str) -> tuple[str, int, str, str]:
    """Minimal representation shared by VCF ingestion and TSV lookup."""
    position = int(pos)
    reference = str(ref or "").strip().upper()
    alternate = str(alt or "").strip().upper()
    if not reference or not alternate:
        raise ValueError("REF and ALT are required")
    while len(reference) > 1 and len(alternate) > 1 and reference[-1] == alternate[-1]:
        reference = reference[:-1]
        alternate = alternate[:-1]
    while len(reference) > 1 and len(alternate) > 1 and reference[0] == alternate[0]:
        reference = reference[1:]
        alternate = alternate[1:]
        position += 1
    return normalize_chrom(chrom), position, reference, alternate


def variant_key(chrom: str, pos: str | int, ref: str, alt: str) -> str:
    normalized = normalize_variant(chrom, pos, ref, alt)
    return ":".join(str(value) for value in normalized)


def parse_info(raw: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for field in str(raw or "").split(";"):
        if not field:
            continue
        key, separator, value = field.partition("=")
        out[key] = unquote(value) if separator else "1"
    return out


def review_stars(review_status: str) -> int:
    value = str(review_status or "").lower().replace(" ", "_")
    if "practice_guideline" in value:
        return 4
    if "reviewed_by_expert_panel" in value:
        return 3
    if "criteria_provided,_multiple_submitters,_no_conflicts" in value:
        return 2
    if "multiple_submitters" in value and "no_conflicts" in value:
        return 2
    if "criteria_provided" in value or "conflicting" in value:
        return 1
    return 0


def normalize_significance(value: str) -> str:
    return str(value or "").replace("_", " ").replace("|", "; ").strip()


def normalize_variation_id(value: str) -> str:
    identifier = str(value or "").strip().strip(".")
    if identifier.upper().startswith("VCV"):
        identifier = identifier[3:].split(".", 1)[0]
    return identifier.lstrip("0") or ("0" if identifier else "")


def significance_bucket(value: str) -> str:
    """Return a clinically meaningful bucket for change-arrow decisions."""
    text = normalize_significance(value).lower()
    if not text or text in {".", "no record", "not provided"}:
        return "no_record"
    if "conflict" in text:
        return "conflicting"
    if "pathogenic" in text and "benign" not in text:
        return "plp"
    if "benign" in text and "pathogenic" not in text:
        return "blb"
    if "uncertain" in text or "vus" in text:
        return "vus"
    return "other"


def meaningful_change(old_significance: str, new_significance: str) -> str:
    old_bucket = significance_bucket(old_significance)
    new_bucket = significance_bucket(new_significance)
    if old_bucket != "plp" and new_bucket == "plp":
        return "UP_TO_PLP"
    if old_bucket == "plp" and new_bucket != "plp":
        return "DOWN_FROM_PLP"
    return ""


def restore_pipeline_variant(variant: dict) -> dict:
    """Copy one UI variant with its fixed Nextflow ClinVar values restored."""
    if not variant.get("clinvar_latest_applied"):
        return dict(variant)
    restored = dict(variant)
    for current_field, baseline_field in _BASELINE_VARIANT_FIELDS.items():
        value = variant.get(baseline_field)
        restored[current_field] = (
            value
            if current_field == "clinvar_stars"
            else value or ""
        )
    return restored


def restore_pipeline_variants(variants: dict) -> dict:
    """Restore baseline ClinVar without mutating sample-loader caches."""
    return {
        variant_id: restore_pipeline_variant(variant)
        for variant_id, variant in (variants or {}).items()
    }


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open(
        "r", encoding="utf-8", errors="replace"
    )


def _atomic_database_target(path: Path) -> tuple[Path, Path]:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    return Path(name), path


def build_database(vcf_path: Path, db_path: Path, *, release_date: str = "") -> dict[str, object]:
    """Build a replacement SQLite snapshot from the official GRCh38 VCF."""
    vcf_path = Path(vcf_path)
    temp_db, final_db = _atomic_database_target(Path(db_path))
    parsed_release = release_date
    record_count = 0
    try:
        with sqlite3.connect(temp_db) as connection:
            connection.executescript(
                """
                PRAGMA journal_mode=OFF;
                PRAGMA synchronous=OFF;
                CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
                CREATE TABLE records (
                    variant_key TEXT PRIMARY KEY,
                    variation_id TEXT NOT NULL DEFAULT '',
                    chrom TEXT NOT NULL,
                    pos INTEGER NOT NULL,
                    ref TEXT NOT NULL,
                    alt TEXT NOT NULL,
                    significance TEXT NOT NULL DEFAULT '',
                    stars INTEGER NOT NULL DEFAULT 0,
                    disease TEXT NOT NULL DEFAULT '',
                    significance_conflict TEXT NOT NULL DEFAULT '',
                    review_status TEXT NOT NULL DEFAULT ''
                );
                CREATE INDEX records_variation_id_idx ON records(variation_id);
                """
            )
            batch: list[tuple[object, ...]] = []
            with _open_text(vcf_path) as handle:
                for line in handle:
                    if line.startswith("##"):
                        match = _FILE_DATE_RE.match(line)
                        if match and not parsed_release:
                            parsed_release = "-".join(match.groups())
                        continue
                    if line.startswith("#") or not line.strip():
                        continue
                    columns = line.rstrip("\r\n").split("\t")
                    if len(columns) < 8:
                        continue
                    info = parse_info(columns[7])
                    for alt in columns[4].split(","):
                        if alt in {"", ".", "*"} or alt.startswith("<"):
                            continue
                        try:
                            chrom, pos, ref, normalized_alt = normalize_variant(
                                columns[0], columns[1], columns[3], alt
                            )
                        except (TypeError, ValueError):
                            continue
                        key = f"{chrom}:{pos}:{ref}:{normalized_alt}"
                        batch.append(
                            (
                                key,
                                normalize_variation_id(
                                    columns[2] or info.get("CLNVID") or info.get("ALLELEID") or ""
                                ),
                                chrom,
                                pos,
                                ref,
                                normalized_alt,
                                str(info.get("CLNSIG", "")).strip(),
                                review_stars(info.get("CLNREVSTAT", "")),
                                str(info.get("CLNDN", "")).strip(),
                                str(info.get("CLNSIGCONF", "")).strip(),
                                str(info.get("CLNREVSTAT", "")).strip(),
                            )
                        )
                        if len(batch) >= 20_000:
                            connection.executemany(
                                "INSERT OR REPLACE INTO records VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                                batch,
                            )
                            record_count += len(batch)
                            batch.clear()
                if batch:
                    connection.executemany(
                        "INSERT OR REPLACE INTO records VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                        batch,
                    )
                    record_count += len(batch)
            record_count = int(connection.execute("SELECT count(*) FROM records").fetchone()[0])
            metadata = {
                "schema_version": SCHEMA_VERSION,
                "release_date": parsed_release,
                "source_path": str(vcf_path),
                "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "record_count": str(record_count),
            }
            connection.executemany("INSERT INTO meta(key, value) VALUES (?, ?)", metadata.items())
            connection.commit()
        os.replace(temp_db, final_db)
    except BaseException:
        temp_db.unlink(missing_ok=True)
        raise
    return {
        "release_date": parsed_release,
        "record_count": record_count,
        "database": str(final_db),
    }


def metadata(db_path: Path) -> dict[str, str]:
    try:
        with sqlite3.connect(db_path) as connection:
            return dict(connection.execute("SELECT key, value FROM meta"))
    except (OSError, sqlite3.Error):
        return {}


def _chunks(values: list[str], size: int = 400) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield values[start:start + size]


def lookup_records(
    db_path: Path,
    keys: Iterable[str],
    variation_ids: Iterable[str] = (),
) -> tuple[dict[str, dict[str, object]], dict[str, dict[str, object]]]:
    """Batch lookup by normalized allele key and, secondarily, ClinVar ID."""
    unique_keys = sorted({str(value) for value in keys if value})
    unique_ids = sorted({str(value) for value in variation_ids if value})
    by_key: dict[str, dict[str, object]] = {}
    by_id: dict[str, dict[str, object]] = {}
    columns = (
        "variant_key, variation_id, chrom, pos, ref, alt, significance, stars, "
        "disease, significance_conflict, review_status"
    )
    with sqlite3.connect(db_path) as connection:
        connection.row_factory = sqlite3.Row
        for chunk in _chunks(unique_keys):
            placeholders = ",".join("?" for _ in chunk)
            for row in connection.execute(
                f"SELECT {columns} FROM records WHERE variant_key IN ({placeholders})", chunk
            ):
                record = dict(row)
                by_key[str(row["variant_key"])] = record
        for chunk in _chunks(unique_ids):
            placeholders = ",".join("?" for _ in chunk)
            candidates: dict[str, list[dict[str, object]]] = {}
            for row in connection.execute(
                f"SELECT {columns} FROM records WHERE variation_id IN ({placeholders})", chunk
            ):
                record = dict(row)
                by_key.setdefault(str(row["variant_key"]), record)
                candidates.setdefault(str(row["variation_id"]), []).append(record)
            # A Variation ID can describe a haplotype/set with multiple
            # alleles. Use it as a representation fallback only when unique;
            # exact normalized allele matches above remain authoritative.
            for variation_id, records in candidates.items():
                if len(records) == 1:
                    by_id[variation_id] = records[0]
    return by_key, by_id


def write_manifest(path: Path, payload: dict[str, object]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    target = path.with_suffix(path.suffix + ".tmp")
    target.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(target, path)
