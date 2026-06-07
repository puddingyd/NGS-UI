"""MITOMAP annotation helpers for pipeline mito TSVs.

The tertiary pipeline's v3.2 mito TSV is VEP/ClinVar/gnomAD based.
For backward compatibility with the old local MITOMAP-only workflow we
silently add MITOMAP columns after copying the TSV into NGS-UI's
tertiary_output directory. The GUI does not render these columns, but
the mito adapter can use confirmed/pathogenic MITOMAP statuses for tier
classification.
"""
from __future__ import annotations

import csv
import os
import re
from pathlib import Path

_RNA_ALLELE_RE = re.compile(r"^([ACGTN])(\d+)([ACGTN])$", re.I)

MITOMAP_COLUMNS = [
    "MITOMAP_DISEASE",
    "MITOMAP_STATUS",
    "MITOMAP_PLASMY",
    "MITOMAP_GB_FREQ",
    "MITOMAP_GB_SEQS",
    "MITOMAP_REFS",
    "MITOTIP_SCORE",
    "MITOMAP_ALLELE",
]

_EMPTY_MM = {
    "disease": "",
    "status": "",
    "plasmy": "",
    "gb_freq": "",
    "gb_seqs": "",
    "refs": "",
    "mitotip": "",
    "allele": "",
}


def _mitomap_dir() -> Path:
    ref_dir = os.environ.get("REF_DIR", "/home/pipeline/reference/hg38")
    return Path(os.environ.get("MITOMAP_DIR", f"{ref_dir}/tertiary/mitomap"))


def _load_mitomap_cc(path: Path) -> dict[tuple[int, str, str], dict]:
    out: dict[tuple[int, str, str], dict] = {}
    if not path.is_file():
        return out
    with path.open("r", encoding="latin-1", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                pos = int((row.get("Position") or "").strip())
            except ValueError:
                continue
            nuc = (row.get("Nucleotide Change") or "").strip()
            match = re.fullmatch(r"([ACGTN])-([ACGTN])", nuc, re.I)
            if not match:
                continue
            out[(pos, match.group(1).upper(), match.group(2).upper())] = {
                "disease": (row.get("Disease") or "").strip(),
                "status": (row.get("Status") or "").strip(),
                "plasmy": (row.get("Plasmy Reports (Homo/Hetero)") or "").strip(),
                "gb_freq": (row.get("GB Freq FL(CR)") or "").strip(),
                "gb_seqs": (row.get("GB Seqs FL(CR)") or "").strip(),
                "refs": (row.get("References") or "").strip(),
                "mitotip": "",
                "allele": (row.get("Allele") or "").strip(),
            }
    return out


def _load_mitomap_rna(path: Path) -> dict[tuple[int, str, str], dict]:
    out: dict[tuple[int, str, str], dict] = {}
    if not path.is_file():
        return out
    with path.open("r", encoding="latin-1", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                pos = int((row.get("Position") or "").strip())
            except ValueError:
                continue
            allele = (row.get("Allele") or "").strip()
            match = _RNA_ALLELE_RE.fullmatch(allele)
            if not match:
                continue
            mitotip = (row.get("MitoTIP") or "").strip()
            out[(pos, match.group(1).upper(), match.group(3).upper())] = {
                "disease": (row.get("Disease") or "").strip(),
                "status": (row.get("Status") or "").strip(),
                "plasmy": (
                    f"{(row.get('Homoplasmy') or '').strip()}/"
                    f"{(row.get('Heteroplasmy') or '').strip()}"
                ),
                "gb_freq": (row.get("GB Freq FL(CR)") or "").strip(),
                "gb_seqs": (row.get("GB Seqs FL(CR)") or "").strip(),
                "refs": (row.get("References") or "").strip(),
                "mitotip": "" if mitotip.upper() in ("", "N/A") else mitotip,
                "allele": allele,
            }
    return out


def _lookup(
    cc: dict[tuple[int, str, str], dict],
    rna: dict[tuple[int, str, str], dict],
    pos: int,
    ref: str,
    alt: str,
    gene: str,
    biotype: str,
) -> dict:
    key = (pos, ref.upper(), alt.upper())
    g = (gene or "").upper()
    b = (biotype or "").lower()
    prefer_rna = "trna" in b or "rrna" in b or re.match(r"^MT-T[A-Z0-9]+$", g) or g in {"MT-RNR1", "MT-RNR2"}
    tables = (rna, cc) if prefer_rna else (cc, rna)
    for table in tables:
        rec = table.get(key)
        if rec:
            return rec
    return dict(_EMPTY_MM)


def annotate_mito_tsv(tsv_path: Path) -> bool:
    """Append/fill MITOMAP columns in-place.

    Returns True when the file was rewritten. Missing MITOMAP resources
    are treated as a no-op so the tertiary job can continue silently.
    """
    if not tsv_path.is_file():
        return False
    mm_dir = _mitomap_dir()
    cc = _load_mitomap_cc(mm_dir / "mitomap_mutations_coding_control.tsv")
    rna = _load_mitomap_rna(mm_dir / "mitomap_mutations_rna.tsv")
    if not cc and not rna:
        return False

    with tsv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    if not fieldnames:
        return False

    added = False
    for col in MITOMAP_COLUMNS:
        if col not in fieldnames:
            fieldnames.append(col)
            added = True

    changed = added
    for row in rows:
        try:
            pos = int((row.get("POS") or "").strip())
        except ValueError:
            continue
        ref = (row.get("REF") or "").strip()
        alt = (row.get("ALT") or "").strip()
        if not ref or not alt:
            continue
        rec = _lookup(
            cc,
            rna,
            pos,
            ref,
            alt,
            row.get("GENE") or row.get("MITOMAP_LOCUS") or "",
            row.get("BIOTYPE") or "",
        )
        mapping = {
            "MITOMAP_DISEASE": rec["disease"],
            "MITOMAP_STATUS": rec["status"],
            "MITOMAP_PLASMY": rec["plasmy"],
            "MITOMAP_GB_FREQ": rec["gb_freq"],
            "MITOMAP_GB_SEQS": rec["gb_seqs"],
            "MITOMAP_REFS": rec["refs"],
            "MITOTIP_SCORE": rec["mitotip"],
            "MITOMAP_ALLELE": rec["allele"],
        }
        for col, val in mapping.items():
            if val and not (row.get(col) or "").strip():
                row[col] = val
                changed = True
            else:
                row.setdefault(col, row.get(col, ""))

    if not changed:
        return False
    tmp = tsv_path.with_suffix(tsv_path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(tsv_path)
    return True
