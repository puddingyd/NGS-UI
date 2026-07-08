#!/usr/bin/env python3
"""Build the NGS-UI supplemental gene-disease SQLite index."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sqlite3
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "backend"))

from app.config import GENE_DISEASE_DB, GENE_DISEASE_RAW_DIR  # noqa: E402
from app.services import panel_deadzone  # noqa: E402

ALLOWED_CLASSIFICATIONS = {
    "definitive": 50,
    "strong": 40,
    "moderate": 30,
    "limited": 20,
}
EXCLUDED_CLASSIFICATIONS = {
    "supportive",
    "animal model only",
    "disputed",
    "refuted",
    "no known disease relationship",
    "no reported evidence",
}
MOI_MAP = {
    "Autosomal dominant": "AD",
    "Autosomal recessive": "AR",
    "X-linked": "XL",
    "X-linked dominant": "XLD",
    "X-linked recessive": "XLR",
    "Mitochondrial": "MT",
    "Semidominant": "SD",
    "Digenic": "DG",
}
CURIE_RE = re.compile(r"^([A-Za-z0-9_.-]+):(.+)$")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _file_meta(path: Path, source: str, url: str = "") -> dict:
    st = path.stat()
    return {
        "source": source,
        "path": str(path),
        "url": url,
        "size": st.st_size,
        "mtime_ns": st.st_mtime_ns,
        "sha256": _sha256(path),
    }


def _canonical_gene(symbol: str) -> str:
    return panel_deadzone.canonical_panel_gene_symbol(symbol or "") or (symbol or "").strip()


def _norm_mondo(curie_or_iri: str) -> str:
    value = (curie_or_iri or "").strip()
    if not value:
        return ""
    if value.startswith("http://purl.obolibrary.org/obo/MONDO_"):
        return "MONDO:" + value.rsplit("MONDO_", 1)[1].replace("_", ":")
    if value.startswith("MONDO_"):
        return "MONDO:" + value.split("MONDO_", 1)[1].replace("_", ":")
    if value.upper().startswith("MONDO:"):
        return "MONDO:" + value.split(":", 1)[1]
    return value


def _curie_parts(value: str) -> tuple[str, str]:
    match = CURIE_RE.match((value or "").strip())
    if not match:
        return "", ""
    return match.group(1).upper(), match.group(2)


def _omim_from_curie(value: str) -> str:
    db, ident = _curie_parts(value)
    return ident if db in {"OMIM", "MIM"} and re.fullmatch(r"\d{6}", ident or "") else ""


def _inheritance(value: str) -> str:
    value = (value or "").strip()
    return MOI_MAP.get(value, value)


def _classification_rank(value: str) -> int:
    return ALLOWED_CLASSIFICATIONS.get((value or "").strip().lower(), 0)


def _include_classification(value: str) -> bool:
    key = (value or "").strip().lower()
    if key in EXCLUDED_CLASSIFICATIONS:
        return False
    return key in ALLOWED_CLASSIFICATIONS


def _read_mondo(path: Path) -> tuple[dict[str, dict], dict[str, list[str]], dict[str, list[str]]]:
    with path.open(encoding="utf-8") as fh:
        data = json.load(fh)
    graph = (data.get("graphs") or [{}])[0]
    terms: dict[str, dict] = {}
    xrefs_by_mondo: dict[str, list[str]] = {}
    mondo_by_xref: dict[str, list[str]] = {}
    for node in graph.get("nodes") or []:
        mondo_id = _norm_mondo(node.get("id") or "")
        if not mondo_id.startswith("MONDO:"):
            continue
        meta = node.get("meta") or {}
        if meta.get("deprecated"):
            continue
        terms[mondo_id] = {
            "mondo_id": mondo_id,
            "label": node.get("lbl") or "",
            "definition": ((meta.get("definition") or {}).get("val") or ""),
            "obsolete": 0,
        }
        xrefs = []
        for xref in meta.get("xrefs") or []:
            val = xref.get("val") if isinstance(xref, dict) else str(xref)
            if not val:
                continue
            xrefs.append(val)
            mondo_by_xref.setdefault(val.upper(), []).append(mondo_id)
        xrefs_by_mondo[mondo_id] = xrefs
    return terms, xrefs_by_mondo, mondo_by_xref


def _omim_for_mondo(mondo_id: str, xrefs_by_mondo: dict[str, list[str]]) -> str:
    omims = []
    for xref in xrefs_by_mondo.get(mondo_id, []):
        omim = _omim_from_curie(xref)
        if omim and omim not in omims:
            omims.append(omim)
    return omims[0] if len(omims) == 1 else ""


def _parse_gencc(path: Path, xrefs_by_mondo: dict[str, list[str]]) -> list[dict]:
    out = []
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            classification = (row.get("classification_title") or "").strip()
            if not _include_classification(classification):
                continue
            gene = _canonical_gene(row.get("gene_symbol") or row.get("submitted_as_hgnc_symbol") or "")
            disease_name = (row.get("disease_title") or row.get("submitted_as_disease_name") or "").strip()
            mondo_id = _norm_mondo(row.get("disease_curie") or "")
            phenotype_mim = _omim_from_curie(row.get("disease_original_curie") or "")
            if not phenotype_mim and mondo_id:
                phenotype_mim = _omim_for_mondo(mondo_id, xrefs_by_mondo)
            if not gene or not disease_name:
                continue
            out.append({
                "gene": gene,
                "disease_name": disease_name,
                "source": "GenCC",
                "classification": classification,
                "mondo_id": mondo_id,
                "phenotype_mim": phenotype_mim,
                "source_disease_id": row.get("uuid") or row.get("disease_original_curie") or mondo_id,
                "inheritance": _inheritance(row.get("moi_title") or row.get("submitted_as_moi_name") or ""),
                "evidence_rank": _classification_rank(classification),
                "include_in_ui": 1,
            })
    return out


def _clingen_reader(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh))
    header_idx = None
    for idx, row in enumerate(rows):
        normalized = [c.strip().upper() for c in row]
        if "GENE SYMBOL" in normalized and "CLASSIFICATION" in normalized:
            header_idx = idx
            break
    if header_idx is None:
        raise RuntimeError(f"could not find ClinGen CSV header in {path}")
    headers = rows[header_idx]
    for row in rows[header_idx + 1:]:
        if not row or row[0].startswith("+"):
            continue
        if len(row) < len(headers):
            row += [""] * (len(headers) - len(row))
        yield dict(zip(headers, row))


def _parse_clingen(path: Path, xrefs_by_mondo: dict[str, list[str]]) -> list[dict]:
    out = []
    for row in _clingen_reader(path):
        classification = (row.get("CLASSIFICATION") or "").strip()
        if not _include_classification(classification):
            continue
        gene = _canonical_gene(row.get("GENE SYMBOL") or "")
        disease_name = (row.get("DISEASE LABEL") or "").strip()
        mondo_id = _norm_mondo(row.get("DISEASE ID (MONDO)") or "")
        phenotype_mim = _omim_for_mondo(mondo_id, xrefs_by_mondo) if mondo_id else ""
        if not gene or not disease_name:
            continue
        out.append({
            "gene": gene,
            "disease_name": disease_name,
            "source": "ClinGen",
            "classification": classification,
            "mondo_id": mondo_id,
            "phenotype_mim": phenotype_mim,
            "source_disease_id": row.get("ONLINE REPORT") or mondo_id,
            "inheritance": _inheritance(row.get("MOI") or ""),
            "evidence_rank": _classification_rank(classification),
            "include_in_ui": 1,
        })
    return out


def _create_schema(con: sqlite3.Connection) -> None:
    con.executescript(
        """
        DROP TABLE IF EXISTS source_files;
        DROP TABLE IF EXISTS mondo_terms;
        DROP TABLE IF EXISTS mondo_xrefs;
        DROP TABLE IF EXISTS gene_disease_associations;

        CREATE TABLE source_files (
            source TEXT PRIMARY KEY,
            path TEXT NOT NULL,
            url TEXT,
            size INTEGER NOT NULL,
            mtime_ns INTEGER NOT NULL,
            sha256 TEXT NOT NULL
        );
        CREATE TABLE mondo_terms (
            mondo_id TEXT PRIMARY KEY,
            label TEXT,
            definition TEXT,
            obsolete INTEGER NOT NULL DEFAULT 0
        );
        CREATE TABLE mondo_xrefs (
            mondo_id TEXT NOT NULL,
            db TEXT NOT NULL,
            db_id TEXT NOT NULL,
            xref TEXT NOT NULL,
            PRIMARY KEY (mondo_id, xref)
        );
        CREATE TABLE gene_disease_associations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene TEXT NOT NULL,
            disease_name TEXT NOT NULL,
            source TEXT NOT NULL,
            classification TEXT NOT NULL,
            mondo_id TEXT,
            phenotype_mim TEXT,
            source_disease_id TEXT,
            inheritance TEXT,
            evidence_rank INTEGER NOT NULL,
            include_in_ui INTEGER NOT NULL DEFAULT 1
        );
        CREATE INDEX idx_gene_disease_gene ON gene_disease_associations(gene);
        CREATE INDEX idx_gene_disease_mondo ON gene_disease_associations(mondo_id);
        CREATE INDEX idx_gene_disease_mim ON gene_disease_associations(phenotype_mim);
        CREATE INDEX idx_mondo_xrefs_xref ON mondo_xrefs(xref);
        """
    )


def _write_db(
    target: Path,
    *,
    raw_dir: Path,
    terms: dict[str, dict],
    xrefs_by_mondo: dict[str, list[str]],
    associations: list[dict],
    source_files: list[dict],
) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(delete=False, suffix=".sqlite", dir=str(target.parent)) as tmp:
        tmp_path = Path(tmp.name)
    con = sqlite3.connect(str(tmp_path))
    try:
        _create_schema(con)
        con.executemany(
            "INSERT INTO source_files(source,path,url,size,mtime_ns,sha256) VALUES (:source,:path,:url,:size,:mtime_ns,:sha256)",
            source_files,
        )
        con.executemany(
            "INSERT INTO mondo_terms(mondo_id,label,definition,obsolete) VALUES (:mondo_id,:label,:definition,:obsolete)",
            terms.values(),
        )
        xref_rows = []
        for mondo_id, xrefs in xrefs_by_mondo.items():
            for xref in xrefs:
                db, db_id = _curie_parts(xref)
                if not db:
                    continue
                xref_rows.append({
                    "mondo_id": mondo_id,
                    "db": db,
                    "db_id": db_id,
                    "xref": xref,
                })
        con.executemany(
            "INSERT OR IGNORE INTO mondo_xrefs(mondo_id,db,db_id,xref) VALUES (:mondo_id,:db,:db_id,:xref)",
            xref_rows,
        )
        con.executemany(
            """
            INSERT INTO gene_disease_associations(
                gene,disease_name,source,classification,mondo_id,phenotype_mim,
                source_disease_id,inheritance,evidence_rank,include_in_ui
            ) VALUES (
                :gene,:disease_name,:source,:classification,:mondo_id,:phenotype_mim,
                :source_disease_id,:inheritance,:evidence_rank,:include_in_ui
            )
            """,
            associations,
        )
        con.commit()
    finally:
        con.close()
    tmp_path.replace(target)


def _write_flat_tsv(path: Path, associations: list[dict]) -> None:
    fields = [
        "gene", "disease_name", "source", "classification", "mondo_id",
        "phenotype_mim", "source_disease_id", "inheritance", "evidence_rank",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in associations:
            writer.writerow({f: row.get(f, "") for f in fields})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-dir", type=Path, default=GENE_DISEASE_RAW_DIR)
    ap.add_argument("--db", type=Path, default=GENE_DISEASE_DB)
    ap.add_argument("--flat-tsv", type=Path, default=None)
    args = ap.parse_args()

    raw_dir = args.raw_dir
    gencc = raw_dir / "gencc_submissions.csv"
    clingen = raw_dir / "clingen_gene_validity.csv"
    mondo = raw_dir / "mondo.json"
    for path in (gencc, clingen, mondo):
        if not path.is_file():
            raise SystemExit(f"missing source file: {path}")

    print("[build] reading MONDO", flush=True)
    terms, xrefs_by_mondo, _mondo_by_xref = _read_mondo(mondo)
    print(f"[build] MONDO terms={len(terms)}", flush=True)
    associations = []
    associations.extend(_parse_gencc(gencc, xrefs_by_mondo))
    associations.extend(_parse_clingen(clingen, xrefs_by_mondo))
    associations.sort(key=lambda r: (r["gene"], -r["evidence_rank"], r["source"], r["disease_name"]))
    source_files = [
        _file_meta(gencc, "GenCC", "https://thegencc.org/download/action/submissions-export-csv"),
        _file_meta(clingen, "ClinGen", "https://search.clinicalgenome.org/kb/gene-validity/download"),
        _file_meta(mondo, "MONDO", "https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.json"),
    ]

    print(f"[build] associations={len(associations)} genes={len({r['gene'] for r in associations})}", flush=True)
    _write_db(
        args.db,
        raw_dir=raw_dir,
        terms=terms,
        xrefs_by_mondo=xrefs_by_mondo,
        associations=associations,
        source_files=source_files,
    )
    flat = args.flat_tsv or (args.db.parent / "gene_disease.tsv")
    _write_flat_tsv(flat, associations)
    counts = Counter(r["source"] for r in associations)
    manifest = {
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "db": str(args.db),
        "flat_tsv": str(flat),
        "source_files": source_files,
        "association_count": len(associations),
        "gene_count": len({r["gene"] for r in associations}),
        "source_counts": dict(counts),
        "mondo_term_count": len(terms),
    }
    manifest_path = args.db.parent / "source_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    report_path = args.db.parent / "build_report.tsv"
    with report_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key in ("association_count", "gene_count", "mondo_term_count"):
            writer.writerow([key, manifest[key]])
        for source, count in counts.items():
            writer.writerow([f"source_count.{source}", count])
    print(f"[build] wrote {args.db}")
    print(f"[build] wrote {flat}")
    print(f"[build] wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
