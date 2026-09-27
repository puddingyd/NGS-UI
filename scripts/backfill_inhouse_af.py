#!/usr/bin/env python3
"""Backfill INHOUSE_AC/AN/AF into already-analysed samples.

Run this for samples analysed before the in-house AF step existed, and after
every refresh of the in-house AF DB (the cohort grows, so AN and every AF
change). New analyses don't need it — run_stopgaps.sh does the same work
inline.

It mirrors what the tertiary worker does. Which of the two modes applies is
decided by THE RESOLVED FILE, not by sample_layout.uses_unified_layout(): a
"legacy" sample can still resolve to a pipeline 03_acmg file, because
snv_raw_tsv() falls back to <legacy pipeline root>/<sample>/03_acmg/… when no
UI copy exists. Only a file literally named snv_indel.annotated.tsv outside
03_acmg is the UI's own copy and safe to rewrite; do_inplace() additionally
refuses any path containing 03_acmg.

OVERLAY mode — the resolved TSV is a pipeline source (03_acmg, under either
the unified or the legacy pipeline root):
    03_acmg/<source>.snv_indel.acmg.tsv   IMMUTABLE pipeline source of truth;
                                          exported reports read it. NEVER written.
    08_postprocessing/<sid>.snv_annotations.sqlite
                                          sparse overlay holding every field
                                          post-processing adds (GeneBe,
                                          SpliceAI, MANE, LitVar2, INHOUSE_*)
    08_postprocessing/<sid>.snv_indel.review.tsv
                                          main-screen rows, built from raw+overlay
  The overlay is updated IN PLACE, only for the INHOUSE_* fields: read the raw
  TSV (read-only), join it against the DB, strip every old INHOUSE_* value from
  the overlay, merge the new ones into each row's existing payload, all in one
  SQLite transaction. Other annotations are never touched, and the overlay's
  source signature stays valid because the raw TSV does not change. Then only
  the INHOUSE_* columns of the review TSV are patched (its row set does not
  depend on them); it is rebuilt in full only when it wasn't current or was
  built for another test type. The test type follows the backend loader's rule
  (see sample_test_type), so the UI doesn't rebuild it again on open. The gene index is built from the raw TSV, which never
  changes, so it is NOT rebuilt here.

  The earlier approach (materialise raw+overlay into a full working TSV,
  annotate it, re-diff it with build_snv_annotation_overlay.py) produced the
  same overlay but moved several times the size of a WGS TSV over NFS per
  sample — 15+ minutes each. --selftest checks the direct update against that
  reference path.

IN-PLACE mode — the resolved TSV is the legacy UI copy
(<NGS_UI_HOME>/tertiary_output/<sample>/snv_indel.annotated.tsv):
    that file is itself the annotated copy, so it is rewritten in place and the
    gene index MUST be rebuilt with it — the atomic replace shifts every byte
    offset and the index is only auto-rebuilt when missing.

Derived artifacts (overlay, review TSV, gene index) always live where
state_dir() says, which for most samples is the unified 08_postprocessing even
when the source TSV still sits under the legacy pipeline root.

Paths come from backend/app/services/sample_layout.py, so this follows the
backend exactly.

Usage:
    scripts/backfill_inhouse_af.py                  # every sample
    scripts/backfill_inhouse_af.py SID1 SID2 ...    # only these
    scripts/backfill_inhouse_af.py --dry-run        # show plan, run nothing
    scripts/backfill_inhouse_af.py --selftest       # synthetic check, no real data
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import signal
import sqlite3
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "backend"))

from app import config  # noqa: E402
from app.services import sample_layout as layout  # noqa: E402
from app.services import snv_overlay, snv_review, test_types  # noqa: E402
from app.services.snv_overlay import KEY_FIELDS, OverlayReader  # noqa: E402

SCRIPTS = REPO / "scripts"
sys.path.insert(0, str(SCRIPTS))
import annotate_inhouse_af as iaf  # noqa: E402

INH = (iaf.COL_AC, iaf.COL_AN, iaf.COL_AF)
csv.field_size_limit(4 * 1024 * 1024)
_T0 = time.monotonic()


def log(msg: str) -> None:
    print(f"      [{time.monotonic() - _T0:7.1f}s] {msg}", flush=True)


def derived_paths(sample_id: str):
    """Same names run_stopgaps.sh uses: sample-prefixed under 08_postprocessing,
    bare in the legacy tree."""
    post = layout.state_dir(sample_id, for_write=True)
    pre = f"{sample_id}." if post.name == layout.POSTPROCESSING_DIRNAME else ""
    return {
        "post": post,
        "review": post / f"{pre}snv_indel.review.tsv",
        "manifest": post / f"{pre}snv_indel.review.tsv.source.json",
        "gene_index": post / f"{pre}snv_gene_index.sqlite",
        "overlay": post / f"{pre}snv_annotations.sqlite",
    }


def materialise(raw: Path, overlay: Path, out: Path) -> None:
    """raw + overlay -> full annotated TSV. Only used by --selftest as the
    reference the direct overlay update must reproduce; far too much I/O for
    real WGS samples."""
    have = overlay.is_file()
    with OverlayReader(raw, overlay if have else None) as rd, \
            open(raw, "r", encoding="utf-8", newline="") as fi, \
            open(out, "w", encoding="utf-8", newline="") as fo:
        rdr = csv.DictReader(fi, delimiter="\t")
        fields = list(rdr.fieldnames or [])
        fields += [f for f in (rd.fields if rd.active else []) if f not in fields]
        wtr = csv.DictWriter(fo, fieldnames=fields, delimiter="\t",
                             extrasaction="ignore", lineterminator="\n")
        wtr.writeheader()
        for row in rdr:
            wtr.writerow(rd.apply(row) if rd.active else row)


def _cell(row, i) -> str:
    return row[i] if i is not None and i < len(row) else ""


def _row_values(i, nalts, hits):
    """INHOUSE_AC/AN/AF for raw row i, exactly as annotate_inhouse_af.write_out
    fills them: one slot per ALT token ('.' for misses), '' when no ALT hit."""
    n = nalts[i] if i < len(nalts) else 0
    base = i << 8
    ac, an, af, any_hit = [], [], [], False
    for j in range(min(n, iaf._MAX_ALT + 1)):
        v = hits.get(base | j)
        if v:
            any_hit = True
            ac.append(v[0]); an.append(v[1]); af.append(v[2])
        else:
            ac.append("."); an.append("."); af.append(".")
    if not any_hit:
        return ("", "", "")
    return (",".join(ac), ",".join(an), ",".join(af))


_JSON_RM = "json_remove(payload_json, " + ", ".join(f"'$.{c}'" for c in INH) + ")"
_SQL_UPSERT = (
    "INSERT INTO annotations(row_key, payload_json) VALUES (?, ?) "
    "ON CONFLICT(row_key) DO UPDATE SET payload_json = "
    f"json_patch({_JSON_RM}, excluded.payload_json)")
_SQL_STRIP = (f"UPDATE annotations SET payload_json = {_JSON_RM} "
              "WHERE row_key = ? AND payload_json LIKE '%\"INHOUSE_%'")


def _has_sql_json(conn) -> bool:
    """UPSERT (3.24+) and JSON1 let SQLite merge payloads in C, with no
    SELECT round trip and no Python json per row."""
    if sqlite3.sqlite_version_info < (3, 24, 0):
        return False
    try:
        conn.execute("SELECT json_patch(json_remove('{}', '$.a'), '{}')").fetchone()
        return True
    except sqlite3.OperationalError:
        return False


class _PyMerge:
    """Fallback for an SQLite without UPSERT/JSON1: same result, in Python."""

    def __init__(self, conn):
        self.conn = conn

    def _load(self, keys):
        out = {}
        for start in range(0, len(keys), 800):
            part = keys[start:start + 800]
            q = ",".join("?" for _ in part)
            for key, raw_payload in self.conn.execute(
                    f"SELECT row_key, payload_json FROM annotations WHERE row_key IN ({q})",
                    part):
                try:
                    payload = json.loads(raw_payload)
                except (TypeError, json.JSONDecodeError):
                    payload = {}
                out[key] = payload if isinstance(payload, dict) else {}
        return out

    def upsert(self, pairs):
        existing = self._load([k for k, _ in pairs])
        out = []
        for key, add in pairs:
            merged = {k: v for k, v in existing.get(key, {}).items() if k not in INH}
            merged.update(json.loads(add))
            existing[key] = merged
            out.append((key, json.dumps(merged, ensure_ascii=False)))
        self.conn.executemany(
            "INSERT OR REPLACE INTO annotations(row_key, payload_json) VALUES (?, ?)", out)

    def strip(self, keys):
        existing = self._load([k for (k,) in keys])
        for key, payload in existing.items():
            if not any(c in payload for c in INH):
                continue
            left = {k: v for k, v in payload.items() if k not in INH}
            self.conn.execute("UPDATE annotations SET payload_json = ? WHERE row_key = ?",
                              (json.dumps(left, ensure_ascii=False), key))


def update_overlay(raw: Path, overlay: Path, db: Path) -> int:
    """Put fresh INHOUSE_* into the sparse overlay without rewriting anything
    else. Reads the raw TSV twice (read-only) and the DB once. Returns the
    number of raw rows with an in-house hit."""
    with open(raw, "r", encoding="utf-8", newline="") as f:
        header = next(csv.reader(f, delimiter="\t"), None) or []
    cols = {name: i for i, name in enumerate(header)}
    missing = [c for c in ("CHROM", "POS", "REF", "ALT") if c not in cols]
    if missing:
        raise RuntimeError(f"{raw} has no {','.join(missing)} column")

    have = overlay.is_file()
    if have and not snv_overlay.is_current(raw, overlay):
        # The UI already ignores a stale overlay; refilling it would still
        # leave GeneBe/SpliceAI/... describing an older raw TSV.
        raise RuntimeError(
            f"overlay is stale (raw TSV changed after post-processing): {overlay}\n"
            f"         re-run the tertiary post-processing for this sample first")

    log(f"scan raw TSV ({raw.stat().st_size / 1e9:.2f} GB)")
    index, nalts = iaf.scan_keys(raw, cols["CHROM"], cols["POS"], cols["REF"], cols["ALT"])
    log(f"{len(nalts):,} rows, {len(index):,} allele keys; join in-house DB")
    hits = iaf.join_scan(index, str(db))
    del index
    hit_row = bytearray(len(nalts))
    for packed in hits:
        hit_row[packed >> 8] = 1
    n_hit = sum(hit_row)
    log(f"{n_hit:,} rows matched; update overlay")

    raw_inh = [cols.get(c) for c in INH]
    # If the raw TSV already carries INHOUSE_* columns (not the case for the
    # pipeline 03_acmg today), a miss must still override them with ''.
    every_row = any(i is not None for i in raw_inh)
    key_idx = [cols.get(f) for f in KEY_FIELDS]

    if have:
        conn = sqlite3.connect(overlay, isolation_level=None)
        target = None
    else:
        fd, tmp_name = tempfile.mkstemp(dir=str(overlay.parent),
                                        prefix=overlay.name + ".", suffix=".tmp")
        os.close(fd)
        target = Path(tmp_name)
        target.unlink()
        conn = sqlite3.connect(target, isolation_level=None)
        conn.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
        conn.execute("CREATE TABLE annotations ("
                     "row_key TEXT PRIMARY KEY, payload_json TEXT NOT NULL)")
    conn.execute("PRAGMA cache_size=-262144")   # 256 MB: big WGS overlays
    sql_json = _has_sql_json(conn)
    py = None if sql_json else _PyMerge(conn)
    try:
        conn.execute("BEGIN IMMEDIATE")
        # One pass over raw. A row with a hit gets its old INHOUSE_* replaced
        # by the new values (other fields kept); a row without one only has
        # stale INHOUSE_* removed. Payloads left empty are deleted at the end.
        ups, strips = [], []

        def flush():
            if ups:
                conn.executemany(_SQL_UPSERT, ups) if sql_json else py.upsert(ups)
                ups.clear()
            if strips:
                conn.executemany(_SQL_STRIP, strips) if sql_json else py.strip(strips)
                strips.clear()

        with open(raw, "r", encoding="utf-8", newline="") as f:
            rdr = csv.reader(f, delimiter="\t")
            next(rdr, None)
            for i, row in enumerate(rdr):
                key = json.dumps([_cell(row, k) for k in key_idx],
                                 ensure_ascii=False, separators=(",", ":"))
                add = None
                if every_row or (i < len(hit_row) and hit_row[i]):
                    vals = _row_values(i, nalts, hits)
                    add = {c: v for c, v, ri in zip(INH, vals, raw_inh)
                           if v != _cell(row, ri)}
                if add:
                    ups.append((key, json.dumps(add, ensure_ascii=False)))
                else:
                    strips.append((key,))
                if len(ups) + len(strips) >= 20000:
                    flush()
        flush()
        conn.execute("DELETE FROM annotations WHERE payload_json = '{}'")

        meta = dict(conn.execute("SELECT key, value FROM meta").fetchall())
        raw_fields = json.loads(meta.get("raw_fields_json") or "null") or header
        fields = json.loads(meta.get("annotated_fields_json") or "null") or list(header)
        fields += [c for c in INH if c not in fields]
        meta.update({
            "annotated_fields_json": json.dumps(fields, ensure_ascii=False),
            "overlay_fields_json": json.dumps(
                sorted(set(fields) - set(raw_fields)), ensure_ascii=False),
            "overlay_rows": str(conn.execute(
                "SELECT count(*) FROM annotations").fetchone()[0]),
        })
        if not have:
            meta.update({
                "schema_version": snv_overlay.SCHEMA_VERSION,
                **snv_overlay.source_signature(raw),
                "raw_fields_json": json.dumps(header, ensure_ascii=False),
                "raw_rows": str(len(nalts)),
                "annotated_rows": str(len(nalts)),
                "skipped_raw_rows": "0",
            })
        conn.executemany("INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)",
                         sorted(meta.items()))
        conn.execute("COMMIT")
    except BaseException:
        if conn.in_transaction:
            conn.execute("ROLLBACK")
        conn.close()
        if target is not None:
            target.unlink(missing_ok=True)
        raise
    conn.close()
    if target is not None:
        os.replace(target, overlay)
    return n_hit


def is_ui_copy(raw: Path) -> bool:
    """True only for the legacy UI *annotated copy*, which is ours to rewrite.

    Do NOT decide this from sample_layout.uses_unified_layout(): a sample can
    be 'legacy' yet still resolve to a pipeline 03_acmg file, because
    snv_raw_tsv() falls back to <legacy pipeline root>/<sample>/03_acmg/… when
    no UI copy exists. Anything under 03_acmg is a pipeline source of truth
    (exported reports read it) and must never be modified."""
    return raw.name == "snv_indel.annotated.tsv" and "03_acmg" not in raw.parts


def run(cmd) -> None:
    subprocess.run([str(c) for c in cmd], check=True)


def _read_json(path: Path) -> dict:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return value if isinstance(value, dict) else {}


def sample_test_type(sid: str, old_manifest: dict) -> str:
    """The test type the backend loader will ask the review TSV for, so the
    review built here is not thrown away and rebuilt on first open.

    Same rule as sample_loader._effective_test_type (year+T IDs are always
    TITAN-WGS); with no metadata yet (sample not registered), fall back to what
    the tertiary worker used (old manifest), then DRAGEN -> WGS, then WES."""
    meta = _read_json(layout.state_file(sid, "sample_metadata.json"))
    src = _read_json(layout.state_file(sid, "pipeline_source.json"))
    fallback = (str(old_manifest.get("test_type") or "").upper()
                or ("WGS" if str(src.get("pipeline_type") or "").lower() == "dragen" else "")
                or "WES")
    identity = str(meta.get("lis_id") or meta.get("sample_id") or sid)
    return test_types.normalize_test_type(
        meta.get("test_type") or "", sample_id=identity, default=fallback)


def build_review(raw: Path, p: dict, test_type: str) -> None:
    """Full rebuild (one pass over the whole raw TSV)."""
    snv_review.ensure_review_tsv(
        raw, test_type=test_type, output_dir=p["post"], output_path=p["review"],
        manifest_path=p["manifest"], overlay_path=p["overlay"])


def patch_review(raw: Path, p: dict, old_manifest: dict, old_overlay_sig: dict,
                 test_type: str) -> bool:
    """Refresh only INHOUSE_* in the existing review TSV.

    The review row set never depends on INHOUSE_* (selection uses raw fields,
    the BED, gnomAD AF and CLINVAR_CHANGE), so if the review TSV was current
    for (raw, overlay before this update, test type), then after replacing
    those three columns it is current for the new overlay. Anything else it
    was stale for (BED, GPN-MSA) stays exactly as stale as before. Returns
    False — caller does a full rebuild — whenever that can't be shown."""
    review, manifest = p["review"], p["manifest"]
    st = raw.stat()
    if not (review.is_file() and old_manifest
            and old_manifest.get("raw_mtime_ns") == st.st_mtime_ns
            and old_manifest.get("raw_size") == st.st_size
            and str(old_manifest.get("test_type") or "").upper() == test_type.upper()
            and old_manifest.get("overlay") == old_overlay_sig):
        return False
    with open(raw, "r", encoding="utf-8", newline="") as f:
        raw_header = next(csv.reader(f, delimiter="\t"), None) or []
    if any(c in raw_header for c in INH):
        return False          # values would need raw fallbacks; just rebuild

    with open(review, "r", encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows:
        return False
    header = rows[0]
    if any(f not in header for f in KEY_FIELDS):
        return False
    for c in INH:
        if c not in header:
            header.append(c)
    width = len(header)
    key_idx = [header.index(f) for f in KEY_FIELDS]
    inh_idx = [header.index(c) for c in INH]
    body = rows[1:]
    keys = [json.dumps([_cell(r, k) for k in key_idx], ensure_ascii=False,
                       separators=(",", ":")) for r in body]
    payloads = {}
    with sqlite3.connect(p["overlay"]) as conn:
        uniq = list(dict.fromkeys(keys))
        for start in range(0, len(uniq), 800):
            part = uniq[start:start + 800]
            q = ",".join("?" for _ in part)
            for key, raw_payload in conn.execute(
                    f"SELECT row_key, payload_json FROM annotations WHERE row_key IN ({q})",
                    part):
                try:
                    payloads[key] = json.loads(raw_payload)
                except (TypeError, json.JSONDecodeError):
                    pass
    for r, key in zip(body, keys):
        if len(r) < width:
            r.extend([""] * (width - len(r)))
        pl = payloads.get(key) or {}
        for c, i in zip(INH, inh_idx):
            r[i] = str(pl.get(c) or "")

    tmp = review.with_suffix(review.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(body)
    os.replace(tmp, review)
    snv_review._write_manifest(manifest, {
        **old_manifest,
        "overlay": snv_overlay.overlay_signature(raw, p["overlay"]),
    })
    return True


def do_overlay(sid: str, raw: Path, p: dict, db: Path) -> None:
    old_manifest = _read_json(p["manifest"])
    test_type = sample_test_type(sid, old_manifest)
    old_sig = snv_overlay.overlay_signature(raw, p["overlay"])
    n_hit = update_overlay(raw, p["overlay"], db)
    log(f"overlay updated ({n_hit:,} rows with INHOUSE_AF)")
    if patch_review(raw, p, old_manifest, old_sig, test_type):
        log(f"review TSV patched in place (test type {test_type})")
    else:
        why = ("no review TSV yet" if not p["review"].is_file() else
               f"test type {old_manifest.get('test_type')} -> {test_type}"
               if str(old_manifest.get("test_type") or "").upper() != test_type.upper()
               else "review TSV was not current")
        log(f"rebuild review TSV ({why}; test type {test_type})")
        build_review(raw, p, test_type)
    log("done")


def do_inplace(sid: str, raw: Path, p: dict, db: Path) -> None:
    # belt and braces: never rewrite a pipeline source of truth
    if "03_acmg" in raw.parts:
        raise RuntimeError(f"refusing to rewrite pipeline source in place: {raw}")
    log("annotate UI copy in place")
    run([SCRIPTS / "annotate_inhouse_af.py", "--tsv", raw, "--db", db])
    test_type = sample_test_type(sid, _read_json(p["manifest"]))
    log(f"rebuild review TSV (test type {test_type})")
    build_review(raw, p, test_type)
    log("rebuild gene index")
    # raw was rewritten in place -> every byte offset moved
    run([SCRIPTS / "build_snv_gene_index.py", "--tsv", raw, "--out", p["gene_index"]])


def _dump_overlay(path: Path):
    with sqlite3.connect(path) as conn:
        ann = {k: json.loads(v) for k, v in conn.execute(
            "SELECT row_key, payload_json FROM annotations")}
        meta = dict(conn.execute("SELECT key, value FROM meta"))
    return ann, meta


def _selftest_overlay() -> None:
    """Direct overlay update == reference path (materialise -> annotate ->
    build_overlay), raw untouched, other annotations kept, stale INHOUSE_*
    removed, idempotent."""
    import hashlib
    import shutil
    d = Path(tempfile.mkdtemp())
    raw = d / "raw.acmg.tsv"
    raw.write_text(
        "CHROM\tPOS\tREF\tALT\tGENE\tTRANSCRIPT\tHGVS_C\tHGVS_P\tCONSEQUENCE\tNOTE\n"
        "chr1\t45330228\tCAA\tC,CA\tMUTYH\tNM_1\tc.1\tp.1\tfs\t\"[{\"\"tx\"\": 1}]\"\n"
        "chr1\t45330228\tCAA\tC,CA\tMUTYH\tNM_2\tc.1\tp.1\tfs\tsecond transcript\n"
        "chr1\t200\tA\tG\tTP53\tNM_3\tc.2\tp.2\tmis\twas a hit, now a miss\n"
        "chr1\t300\tA\tT\tTP53\tNM_3\tc.3\tp.3\tmis\tgenebe only\n"
        "chr7\t500\tAT\t*,A\tSPAN\tNM_4\tc.4\tp.4\tdel\t\"a,b\"\n"
        "chr2\t9\tA\t*\tX\tNM_5\t\t\t\tsymbolic only\n", encoding="utf-8")
    raw_md5 = hashlib.md5(raw.read_bytes()).hexdigest()

    # existing post-processing result: GeneBe on some rows, OLD in-house AF
    ann = d / "annotated.tsv"
    rows = list(csv.reader(raw.open(encoding="utf-8", newline=""), delimiter="\t"))
    # row 1: GeneBe + old AF; row 3: old AF only, now a DB miss (row must go);
    # row 4: GeneBe only, now a DB hit (GeneBe must survive the merge)
    extra = {1: ("P", "0.1,0.2", "10,20", "2000,2000"),
             3: ("", "0.5", "5", "10"), 4: ("B", "", "", "")}
    with ann.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(rows[0] + ["GENEBE_ACMG_CLASS", "INHOUSE_AF", "INHOUSE_AC", "INHOUSE_AN"])
        for i, r in enumerate(rows[1:], 1):
            g, af, ac, an = extra.get(i, ("", "", "", ""))
            w.writerow(r + [g, af, ac, an])
    base = d / "base.sqlite"
    snv_overlay.build_overlay(raw, ann, base)

    db = d / "inhouse.vcf.gz"
    iaf._mk_db(str(db), [
        ("chr1", 45330228, "CAA", "C", 1235, 2794, "0.442019"),
        ("chr1", 45330228, "CA", "C", 1136, 2794, "0.406586"),
        ("chr7", 500, "AT", "A", 10, 1000, "0.01"),
        ("chr1", 300, "A", "T", 1, 2794, "0.000358"),
    ])

    def reference(ov: Path) -> None:
        work = d / "work.tsv"
        materialise(raw, ov, work)
        iaf.annotate(work, str(db))
        snv_overlay.build_overlay(raw, work, ov)
        work.unlink()

    def full(ov: Path):
        out = d / "full.tsv"
        materialise(raw, ov, out)
        text = out.read_text(encoding="utf-8")
        out.unlink()
        return text

    for label, start in (("existing overlay", base), ("no overlay", None)):
        ref, new = d / "ref.sqlite", d / "new.sqlite"
        for x in (ref, new):
            x.unlink(missing_ok=True)
            if start is not None:
                shutil.copy(start, x)
        reference(ref)
        update_overlay(raw, new, db)
        assert full(new) == full(ref), f"{label}: materialised rows differ"
        a_new, m_new = _dump_overlay(new)
        a_ref, m_ref = _dump_overlay(ref)
        assert a_new == a_ref, f"{label}: overlay rows differ\n{a_new}\n{a_ref}"
        for k in ("annotated_fields_json", "overlay_fields_json", "overlay_rows",
                  "source_path", "source_mtime_ns", "source_size", "schema_version"):
            assert m_new[k] == m_ref[k], f"{label}: meta {k}: {m_new[k]} != {m_ref[k]}"
        assert snv_overlay.is_current(raw, new), f"{label}: overlay no longer current"
        update_overlay(raw, new, db)                         # idempotent
        assert _dump_overlay(new)[0] == a_ref, f"{label}: re-run not idempotent"

    shutil.copy(base, new)
    update_overlay(raw, new, db)
    a, _ = _dump_overlay(new)
    by = {(json.loads(k)[5], json.loads(k)[6]): v for k, v in a.items()}
    assert by[("NM_1", "c.1")] == {"GENEBE_ACMG_CLASS": "P", "INHOUSE_AF": "0.442019,0.406586",
                                   "INHOUSE_AC": "1235,1136", "INHOUSE_AN": "2794,2794"}, by
    assert ("NM_3", "c.2") not in by, by                     # stale-only row removed
    assert by[("NM_3", "c.3")] == {"GENEBE_ACMG_CLASS": "B", "INHOUSE_AF": "0.000358",
                                   "INHOUSE_AC": "1", "INHOUSE_AN": "2794"}, by
    assert by[("NM_4", "c.4")]["INHOUSE_AF"] == ".,0.01", by  # '*' keeps its slot
    assert raw_md5 == hashlib.md5(raw.read_bytes()).hexdigest(), "raw TSV modified"

    # a stale overlay (raw changed after post-processing) is refused, untouched
    stale = d / "stale.sqlite"
    shutil.copy(base, stale)
    before = stale.read_bytes()
    os.utime(raw, ns=(1, 1))
    try:
        update_overlay(raw, stale, db)
        raise AssertionError("stale overlay was not refused")
    except RuntimeError as e:
        assert "stale" in str(e)
    assert stale.read_bytes() == before, "stale overlay was modified"
    shutil.rmtree(d, ignore_errors=True)


def _selftest_review() -> None:
    """patch_review() on the old review TSV == a full rebuild after the update."""
    import shutil
    d = Path(tempfile.mkdtemp())
    raw = d / "raw.acmg.tsv"
    raw.write_text(
        "CHROM\tPOS\tREF\tALT\tGENE\tTRANSCRIPT\tHGVS_C\tHGVS_P\tCONSEQUENCE\tIMPACT\tDP_DV\tNOTE\n"
        "chr1\t1000\tA\tG\tAAA\tNM_1\tc.1\tp.1\tmis\tMODERATE\t30\t\"x,\"\"y\"\"\"\n"
        "chr1\t1100\tC\tT\tAAA\tNM_1\tc.2\tp.2\tmis\tMODERATE\t30\tmiss now\n"
        "chr1\t1200\tG\tA\tAAA\tNM_1\tc.3\tp.3\tmis\tMODERATE\t30\tnew hit\n"
        "chr1\t900000\tG\tA\tBBB\tNM_2\tc.4\tp.4\tmis\tMODERATE\t30\toff-BED\n",
        encoding="utf-8")
    bed = d / "cds.bed"
    bed.write_text("chr1\t900\t2000\n", encoding="utf-8")
    ann = d / "ann.tsv"
    rows = list(csv.reader(raw.open(encoding="utf-8", newline=""), delimiter="\t"))
    old = {1: ("P", "0.1"), 2: ("", "0.2"), 3: ("B", "")}
    with ann.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(rows[0] + ["GENEBE_ACMG_CLASS", "INHOUSE_AF"])
        for i, r in enumerate(rows[1:], 1):
            w.writerow(r + list(old.get(i, ("", ""))))
    db = d / "inhouse.vcf.gz"
    iaf._mk_db(str(db), [("chr1", 1000, "A", "G", 3, 2794, "0.00107"),
                         ("chr1", 1200, "G", "A", 9, 2794, "0.00322")])

    def paths(tag):
        post = d / tag
        post.mkdir()
        return {"post": post, "review": post / "review.tsv",
                "manifest": post / "review.tsv.source.json",
                "overlay": post / "ov.sqlite"}

    def as_rows(path):
        with open(path, encoding="utf-8", newline="") as f:
            return [dict(r) for r in csv.DictReader(f, delimiter="\t")]

    env_old = os.environ.get("NGS_UI_CDS_CANDIDATE_BED")
    os.environ["NGS_UI_CDS_CANDIDATE_BED"] = str(bed)
    try:
        pa, pb = paths("patch"), paths("full")
        for x in (pa, pb):
            snv_overlay.build_overlay(raw, ann, x["overlay"])
            build_review(raw, x, "WGS")
        before = as_rows(pa["review"])
        assert len(before) == 3 and before[1]["INHOUSE_AF"] == "0.2", before

        m = _read_json(pa["manifest"])
        sig = snv_overlay.overlay_signature(raw, pa["overlay"])
        update_overlay(raw, pa["overlay"], db)
        assert patch_review(raw, pa, m, sig, "WGS"), "patch refused"
        update_overlay(raw, pb["overlay"], db)
        build_review(raw, pb, "WGS")
        got, want = as_rows(pa["review"]), as_rows(pb["review"])
        assert got == want, f"patched review != rebuilt review\n{got}\n{want}"
        assert [r["INHOUSE_AF"] for r in got] == ["0.00107", "", "0.00322"], got
        assert got[0]["NOTE"] == 'x,"y"' and got[0]["GENEBE_ACMG_CLASS"] == "P", got
        # the patched manifest is what the backend expects -> no rebuild on open
        mtime = pa["review"].stat().st_mtime_ns
        build_review(raw, pa, "WGS")
        assert pa["review"].stat().st_mtime_ns == mtime, "loader would rebuild"

        # wrong test type / overlay changed since the review -> refuse (full rebuild)
        m = _read_json(pa["manifest"])
        sig = snv_overlay.overlay_signature(raw, pa["overlay"])
        assert not patch_review(raw, pa, m, sig, "WES"), "test-type change not detected"
        # overlay changed by something else after the review was built: the
        # next backfill's pre-update signature no longer matches the manifest
        st = pa["overlay"].stat()      # (explicit bump: coarse fs clocks)
        os.utime(pa["overlay"], ns=(st.st_atime_ns, st.st_mtime_ns + 10**9))
        sig = snv_overlay.overlay_signature(raw, pa["overlay"])
        update_overlay(raw, pa["overlay"], db)
        assert not patch_review(raw, pa, m, sig, "WGS"), "stale review not detected"
    finally:
        if env_old is None:
            os.environ.pop("NGS_UI_CDS_CANDIDATE_BED", None)
        else:
            os.environ["NGS_UI_CDS_CANDIDATE_BED"] = env_old
        shutil.rmtree(d, ignore_errors=True)

    assert sample_test_type("26T00028-dragen", {"test_type": "WES"}) == "TITAN-WGS"
    assert sample_test_type("NOPE-X", {"test_type": "wgs"}) == "WGS"
    assert sample_test_type("NOPE-X", {}) == "WES"


def selftest() -> int:
    global _has_sql_json
    import contextlib
    import io
    real = _has_sql_json
    for label, force_py in (("SQLite JSON1 upsert", False), ("Python fallback", True)):
        if force_py:
            _has_sql_json = lambda _conn: False   # noqa: E731
        elif not real(sqlite3.connect(":memory:")):
            print(f"  ({label}: not available in this SQLite, skipped)")
            continue
        try:
            with contextlib.redirect_stdout(io.StringIO()), \
                    contextlib.redirect_stderr(io.StringIO()):
                _selftest_overlay()
                _selftest_review()
        finally:
            _has_sql_json = real
        print(f"  {label}: OK")
    print(f"selftest OK (SQLite {sqlite3.sqlite_version}) — overlay update == "
          "materialise+annotate+build_overlay, raw untouched, idempotent, stale "
          "overlay refused; patched review TSV == full rebuild; test type rule")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples", nargs="*", help="sample IDs (default: all)")
    ap.add_argument("--db", default=str(config.INHOUSE_AF_DB),
                    help=f"in-house AF sites VCF (default {config.INHOUSE_AF_DB})")
    ap.add_argument("--dry-run", action="store_true",
                    help="list each sample, its layout and resolved paths; run nothing")
    ap.add_argument("--continue-on-error", action="store_true",
                    help="keep going when one sample fails (default: stop)")
    ap.add_argument("--selftest", action="store_true",
                    help="synthetic check of the overlay update; touches no real data")
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    try:               # behave like a normal tool under `| head`
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except (AttributeError, ValueError):
        pass

    db = Path(args.db)
    if not db.is_file():
        print(f"ERROR: in-house AF DB not found: {db}\n"
              f"       deploy it first with scripts/inhouse_af/deploy_inhouse_af_db.sh",
              file=sys.stderr)
        return 1

    sids = args.samples or sorted(layout.iter_sample_ids())
    if not sids:
        print(f"no samples under {layout.unified_root()}")
        return 0

    print(f"backfilling INHOUSE_AF for {len(sids)} sample(s)")
    print(f"  DB   : {db}")
    print(f"  root : {layout.unified_root()}")

    n_ok = n_skip = n_fail = 0
    for sid in sids:
        raw = layout.snv_raw_tsv(sid)
        if not raw.is_file():
            print(f"  - {sid}: no SNV TSV ({raw}), skip", file=sys.stderr)
            n_skip += 1
            continue
        p = derived_paths(sid)
        inplace = is_ui_copy(raw)
        if args.dry_run:
            kind = "legacy UI copy" if inplace else "pipeline source + overlay"
            fate = "REWRITTEN IN PLACE" if inplace else "read-only"
            print(f"  • {sid}  [{kind}]")
            print(f"      raw      {raw}   ({fate})")
            print(f"      overlay  {p['overlay']}")
            print(f"      review   {p['review']}")
            if inplace:
                print(f"      index    {p['gene_index']}   (rebuilt)")
            n_ok += 1
            continue
        p["post"].mkdir(parents=True, exist_ok=True)
        print(f"  • {sid}  [{'in-place' if inplace else 'overlay'}]")
        try:
            (do_inplace if inplace else do_overlay)(sid, raw, p, db)
            n_ok += 1
        except subprocess.CalledProcessError as e:
            n_fail += 1
            print(f"  ! {sid}: {Path(e.cmd[0]).name} failed (exit {e.returncode})",
                  file=sys.stderr)
            if not args.continue_on_error:
                print("stopping (use --continue-on-error to keep going)", file=sys.stderr)
                break
        except Exception as e:  # noqa: BLE001
            n_fail += 1
            print(f"  ! {sid}: {e}", file=sys.stderr)
            if not args.continue_on_error:
                break

    print(f"done. {n_ok} annotated, {n_skip} skipped, {n_fail} failed.")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
