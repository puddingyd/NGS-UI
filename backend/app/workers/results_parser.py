"""Parse Exomiser / LIRICAL CLI outputs into spec-compliant sidecars.

Output format (one row per chr-pos-ref-alt that the tool scored):

    exomiser_results.tsv
        VARIANT_ID, GENE_SYMBOL, EXOMISER_GENE_PHENO_SCORE,
        EXOMISER_GENE_COMBINED_SCORE, EXOMISER_VARIANT_SCORE,
        EXOMISER_RANK, MOI

    lirical_results.tsv
        VARIANT_ID, RANK_LIRICAL_VARIANT, LIRICAL_VARIANT_SCORE,
        LIRICAL_PATHOGENICITY, DISEASE_NAME, DISEASE_CURIE
"""
from __future__ import annotations

import csv
import re
from pathlib import Path


def _norm_chrom(c: str) -> str:
    c = (c or "").strip()
    return c if c.startswith("chr") else f"chr{c}"


def parse_exomiser_variants_tsv(in_path: Path, out_path: Path) -> int:
    """Pull the columns we want from EXOMISER's variant TSV.

    The header line in 14.x outputs starts with `#RANK`; rows follow.
    Return the number of rows written (excluding header).
    """
    if not in_path.exists():
        # Still write a header-only sidecar so downstream join code doesn't blow up.
        out_path.write_text(
            "VARIANT_ID\tGENE_SYMBOL\tEXOMISER_GENE_PHENO_SCORE\tEXOMISER_GENE_COMBINED_SCORE"
            "\tEXOMISER_VARIANT_SCORE\tEXOMISER_RANK\tMOI\n",
            encoding="utf-8",
        )
        return 0

    written = 0
    seen: set[str] = set()
    with in_path.open("r", encoding="utf-8", newline="") as f, \
         out_path.open("w", encoding="utf-8", newline="") as fo:
        # Find the header row (line starting with #RANK)
        header: list[str] | None = None
        rest_lines = []
        for line in f:
            if line.startswith("#RANK"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                break
            # Skip any preamble lines
        if header is None:
            fo.write("VARIANT_ID\tGENE_SYMBOL\tEXOMISER_GENE_PHENO_SCORE\tEXOMISER_GENE_COMBINED_SCORE"
                     "\tEXOMISER_VARIANT_SCORE\tEXOMISER_RANK\tMOI\n")
            return 0

        rest = csv.DictReader(f, fieldnames=header, delimiter="\t")
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow([
            "VARIANT_ID", "GENE_SYMBOL", "EXOMISER_GENE_PHENO_SCORE",
            "EXOMISER_GENE_COMBINED_SCORE", "EXOMISER_VARIANT_SCORE",
            "EXOMISER_RANK", "MOI",
        ])
        for row in rest:
            chrom = _norm_chrom(row.get("CONTIG") or "")
            pos = (row.get("START") or "").strip()
            ref = (row.get("REF") or "").strip()
            alt = (row.get("ALT") or "").strip()
            if not (chrom and pos and ref and alt):
                continue
            vid = f"{chrom}-{pos}-{ref}-{alt}"
            # Take the first (best) row per variant; later duplicates skipped.
            if vid in seen:
                continue
            seen.add(vid)
            w.writerow([
                vid,
                (row.get("GENE_SYMBOL") or "").strip(),
                (row.get("EXOMISER_GENE_PHENO_SCORE") or "").strip(),
                (row.get("EXOMISER_GENE_COMBINED_SCORE") or "").strip(),
                (row.get("EXOMISER_VARIANT_SCORE") or "").strip(),
                (row.get("RANK") or row.get("#RANK") or "").strip(),
                (row.get("MOI") or "").strip(),
            ])
            written += 1
    return written


_LIR_PATHO_RE = re.compile(r"pathogenicity\s*:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")


def _composite_to_score(composite_lr: str | float | None) -> int | str:
    """Clamp compositeLR to [-10, 10] and rescale to [0, 100] integer.

    Mirrors vcf-analysis-hg19.R lines 1496-1499.
    """
    if composite_lr in (None, ""):
        return ""
    try:
        x = float(composite_lr)
    except (TypeError, ValueError):
        return ""
    x = max(-10.0, min(10.0, x))
    return int(round((x + 10.0) / 20.0 * 100.0))


def parse_lirical_variant_tsv(in_path: Path, out_path: Path) -> int:
    """Explode LIRICAL's per-disease rows into one row per variant.

    LIRICAL writes one row per disease, with all the variants that
    voted for it (semicolon-separated) crammed into a `variants` cell.
    We expand that cell back out, parse `pathogenicity:` from each
    sub-token, and keep the row with the highest LIRICAL score per
    chr-pos-ref-alt.
    """
    out_path.write_text(
        "VARIANT_ID\tRANK_LIRICAL_VARIANT\tLIRICAL_VARIANT_SCORE\t"
        "LIRICAL_PATHOGENICITY\tDISEASE_NAME\tDISEASE_CURIE\n",
        encoding="utf-8",
    )
    if not in_path.exists():
        return 0

    # Find the header line; LIRICAL prefixes it with comments starting with `!`.
    with in_path.open("r", encoding="utf-8") as f:
        lines = f.readlines()
    header_idx = next(
        (i for i, ln in enumerate(lines) if ln.startswith("rank\t")),
        -1,
    )
    if header_idx < 0:
        return 0
    header = lines[header_idx].rstrip("\n").split("\t")
    body_lines = lines[header_idx + 1:]

    # Expand into per-variant rows
    best: dict[str, dict] = {}
    for ln in body_lines:
        if not ln.strip() or ln.startswith("!"):
            continue
        cells = ln.rstrip("\n").split("\t")
        row = dict(zip(header, cells))
        rank = (row.get("rank") or "").strip()
        disease_name  = (row.get("diseaseName") or "").strip()
        disease_curie = (row.get("diseaseCurie") or "").strip()
        composite_lr  = row.get("compositeLR")
        score = _composite_to_score(composite_lr)

        for token in (row.get("variants") or "").split(";"):
            t = token.strip()
            if not t:
                continue
            # LIRICAL token format: "2:73985020T>C ... pathogenicity:0.0 [1/1]"
            head = t.split(" ", 1)[0]
            m_chr = re.match(r"(\d+|[XYM]|MT):(\d+)([ACGT*]+)>([ACGT*]+)", head)
            if not m_chr:
                continue
            chrom_raw, pos, ref, alt = m_chr.groups()
            chrom = _norm_chrom(chrom_raw)
            vid = f"{chrom}-{pos}-{ref}-{alt}"
            patho_m = _LIR_PATHO_RE.search(t)
            patho = patho_m.group(1) if patho_m else ""
            current = best.get(vid)
            cur_score = current.get("score") if current else -1
            new_score = score if isinstance(score, int) else -1
            if current is None or new_score > cur_score:
                best[vid] = {
                    "rank":  rank,
                    "score": new_score if isinstance(score, int) else "",
                    "patho": patho,
                    "disease_name":  disease_name,
                    "disease_curie": disease_curie,
                }

    written = 0
    with out_path.open("a", encoding="utf-8", newline="") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        for vid, e in best.items():
            w.writerow([vid, e["rank"], e["score"], e["patho"],
                        e["disease_name"], e["disease_curie"]])
            written += 1
    return written
