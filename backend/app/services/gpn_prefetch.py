"""Batch prefetch of GPN-MSA / GPN-Star scores for several review TSVs.

Building each sample's review TSV queries GPN-MSA (tabix) and GPN-Star V/M
(one Parquet shard scan per chromosome).  For a batch of samples the shard
scans dominate and are repeated once per sample.  Post-processing therefore
queries the *union* of every sample's review keys once and stores the outcome
in a job-private SQLite file; ``gpn_msa.query_scores`` / ``gpn_star.query_scores``
consult it (via ``NGS_UI_GPN_PREFETCH_DB``) before touching the databases.

Scores are allele-exact database lookups, so a key's value is the same
whether it was queried alone or with other keys.  Keys missing from the
prefetch (e.g. a reviewer-marked rescue row) still go to the database, and
the review manifest / status bookkeeping is unchanged.
"""
from __future__ import annotations

import os
import sqlite3
from pathlib import Path

ENV_VAR = "NGS_UI_GPN_PREFETCH_DB"

Key = tuple[str, int, str, str]


def provider_id(kind: str, db_path: Path | str, model: str = "") -> str:
    return f"{kind}|{Path(db_path).resolve()}|{model}"


def _connect(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(path, timeout=60)
    conn.execute(
        "CREATE TABLE IF NOT EXISTS queried ("
        "provider TEXT NOT NULL, chrom TEXT NOT NULL, pos INTEGER NOT NULL, "
        "ref TEXT NOT NULL, alt TEXT NOT NULL, score TEXT, "
        "PRIMARY KEY (provider, chrom, pos, ref, alt)) WITHOUT ROWID"
    )
    return conn


def store(path: Path, provider: str, keys: set[Key], scores: dict[Key, str]) -> None:
    """Record that ``keys`` were queried for ``provider`` and what was found."""
    with _connect(Path(path)) as conn:
        conn.executemany(
            "INSERT OR REPLACE INTO queried VALUES (?, ?, ?, ?, ?, ?)",
            ((provider, chrom, int(pos), ref, alt, scores.get((chrom, pos, ref, alt)))
             for chrom, pos, ref, alt in keys),
        )


def split(provider: str, keys: set[Key]) -> tuple[dict[Key, str], set[Key]]:
    """Return ``(prefetched scores, keys still to query)``.

    Without an active prefetch file every key is returned as missing.
    """
    path = os.environ.get(ENV_VAR, "").strip()
    keys = set(keys)
    if not path or not keys or not Path(path).is_file():
        return {}, keys
    found: dict[Key, str] = {}
    known: set[Key] = set()
    try:
        conn = sqlite3.connect(f"file:{Path(path).resolve()}?mode=ro", uri=True, timeout=60)
    except sqlite3.Error:
        return {}, keys
    try:
        conn.execute("CREATE TEMP TABLE wanted (chrom TEXT, pos INTEGER, ref TEXT, alt TEXT)")
        conn.executemany("INSERT INTO temp.wanted VALUES (?, ?, ?, ?)", keys)
        for chrom, pos, ref, alt, score in conn.execute(
            "SELECT q.chrom, q.pos, q.ref, q.alt, q.score FROM temp.wanted w "
            "JOIN queried q ON q.provider = ? AND q.chrom = w.chrom AND q.pos = w.pos "
            "AND q.ref = w.ref AND q.alt = w.alt",
            (provider,),
        ):
            key = (str(chrom), int(pos), str(ref), str(alt))
            known.add(key)
            if score is not None:
                found[key] = str(score)
    except sqlite3.Error:
        return {}, keys
    finally:
        conn.close()
    return found, keys - known


def build_for_tsvs(
    path: Path,
    items: list[tuple[Path, str]],
    *,
    msa_db: Path | None,
    star_dir: Path | None,
    research_only: bool,
    log=print,
    workers: int = 1,
) -> dict[str, dict]:
    """Prefetch every provider for the review-candidate keys of ``items``.

    ``items`` are ``(working_tsv, test_type)``. Best-effort: a provider that
    is unavailable is simply not prefetched, and the per-sample review build
    then reports its own status exactly as before.
    """
    import csv

    from . import gpn_msa, gpn_star, snv_review
    from .parallel_tsv import fork_map

    bed = snv_review.load_candidate_bed()

    def keys_of(item) -> tuple[set[Key], set[Key]]:
        tsv, test_type = item
        msa_part: set[Key] = set()
        star_part: set[Key] = set()
        with Path(tsv).open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if not snv_review.is_review_retained(row, test_type=test_type, bed=bed):
                    continue
                key = gpn_msa._normalise_key(row)
                if key is not None:
                    msa_part.add(key)
                key = gpn_star._normalise_key(row)
                if key is not None:
                    star_part.add(key)
        return msa_part, star_part

    msa_keys: set[Key] = set()
    star_keys: set[Key] = set()
    for msa_part, star_part in fork_map(keys_of, items, workers):
        msa_keys |= msa_part
        star_keys |= star_part
    summary: dict[str, dict] = {}
    providers = [("gpn_msa", None)]
    if research_only:
        providers += [("gpn_star", gpn_star.MODEL), ("gpn_star_v", gpn_star.V_MODEL)]
    for name, model in providers:
        try:
            if name == "gpn_msa":
                db, tabix = gpn_msa.validate_database(Path(msa_db or gpn_msa.GPN_MSA_DB))
                keys = msa_keys
                scores = gpn_msa.query_scores_db(db, keys, tabix_bin=tabix)
                provider = provider_id("gpn_msa", db)
            else:
                root = Path(star_dir or gpn_star.GPN_STAR_DIR)
                gpn_star.validate_database(root, model=model)
                keys = star_keys
                scores = gpn_star.query_scores_db(root, keys, model)
                provider = provider_id("gpn_star", root, model)
            store(path, provider, keys, scores)
            summary[name] = {"keys": len(keys), "found": len(scores)}
            log(f"[gpn-prefetch] {name}: {len(keys)} unique keys, {len(scores)} scored")
        except Exception as exc:  # noqa: BLE001 — per-sample build reports status
            summary[name] = {"error": str(exc)}
            log(f"[gpn-prefetch] {name}: not prefetched ({exc})")
    return summary
