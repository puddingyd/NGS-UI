"""mito.annotated.tsv → variant payload adapter.

Reads the per-sample mitochondrial TSV copied into
tertiary_output/{sample}/mito.annotated.tsv and shapes it for the
frontend's Mitochondria card. The file can be either the older
NGS-UI MITOMAP-only TSV or the newer tertiary-pipeline v3.2
04_mito/{sample}.mito.tsv; field access below intentionally accepts
both schemas.

Two screens are applied before a variant reaches the card:
  1. If a FILTER column exists, it must be PASS — non-PASS Mutect2-mito flags (weak_evidence,
     base_qual, blacklisted_site, possible_numt, contamination, …) are
     all "don't trust this call" reasons; the GATK Mitochondria
     best-practices keep PASS only. (TLOD is already baked into the
     FILTER decision, so no separate TLOD threshold is needed.) Pipeline
     v3.2 mito TSVs currently omit FILTER, so this screen is skipped there.
  2. Disease-relevance — the variant is either (a) pathogenic per
     MITOMAP/MitoTIP, or (b) carries some MITOMAP disease association.
     Polymorphisms / haplogroup variants with no MITOMAP record are
     dropped (the raw mito VCF has ~150 variants per sample, almost
     all of which are exactly that).

Tiers:
    MITO-1  Pathogenic   — MITOMAP status confirmed-ish ("Cfrm" /
                           "Confirmed" / "[P]" / "[LP]") or a MitoTIP
                           "(likely) pathogenic" call
    MITO-2  Clinical     — has a non-empty MITOMAP_DISEASE (and isn't
                           already in MITO-1)

Both tiers sort by a "disease-relevance" key, most-relevant first:
    (status_rank, mitotip_rank, in_panel_rank, -refs, -heteroplasmy, pos)
where status_rank   = Cfrm/Confirmed 0 · Reported 1 · Conflicting 2 · else 3
      mitotip_rank  = Pathogenic 0 · Likely pathogenic 1 · Possibly… 2 · else 3
      in_panel_rank = 0 if GENE is in the patient's pheno_score gene set else 1
heteroplasmy still tie-breaks (a *disease-associated* variant at high load
is more likely clinically significant) but is no longer the headline sort.
"""
from __future__ import annotations

import csv
import re
from urllib.parse import unquote
from pathlib import Path

from ..services import clinvar_mito

MITO_TIERS = ["MITO-1", "MITO-2"]

_PATHO_STATUS_RE = re.compile(r"\bCfrm\b|\bConfirmed\b|\[L?P\]", re.I)
_PATHO_MITOTIP   = {"pathogenic", "likely pathogenic"}


def _to_float(s):
    if s is None:
        return None
    s = str(s).strip()
    if not s or s.upper() in ("NA", "N/A", "."):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _to_int(s):
    f = _to_float(s)
    if f is None:
        return None
    try:
        return int(f)
    except (TypeError, ValueError):
        return None


def _is_pathogenic(status: str, mitotip: str) -> bool:
    if status and _PATHO_STATUS_RE.search(status):
        return True
    if mitotip and mitotip.strip().lower() in _PATHO_MITOTIP:
        return True
    return False


def _status_rank(status: str) -> int:
    s = (status or "").lower()
    if "cfrm" in s or "confirmed" in s:
        return 0
    if "reported" in s:
        return 1
    if "conflicting" in s or "disputed" in s:
        return 2
    return 3


def _mitotip_rank(mitotip: str) -> int:
    m = (mitotip or "").strip().lower()
    if m == "pathogenic":
        return 0
    if m == "likely pathogenic":
        return 1
    if m.startswith("possibly"):
        return 2
    return 3


def _refs_count(refs: str) -> int:
    """MITOMAP "References" is usually an integer count; be lenient."""
    n = _to_int(refs)
    if n is not None:
        return n
    # fall back to counting non-empty tokens
    return len([t for t in re.split(r"[;, ]+", refs or "") if t])


def _first(row: dict, *names: str) -> str:
    for name in names:
        val = row.get(name)
        if val is None:
            continue
        s = str(val).strip()
        if s and s.upper() not in ("NA", "N/A", "."):
            return s
    return ""


def _derive_locus_type(row: dict, gene: str) -> str:
    explicit = _first(row, "LOCUS_TYPE")
    if explicit:
        return explicit
    biotype = _first(row, "BIOTYPE").lower()
    consequence = _first(row, "CONSEQUENCE").lower()
    g = (gene or "").upper()
    if "trna" in biotype or re.match(r"^MT-T[A-Z0-9]+$", g):
        return "tRNA"
    if "rrna" in biotype or g in {"MT-RNR1", "MT-RNR2"}:
        return "rRNA"
    if "protein" in biotype or g in {
        "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
        "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ATP6", "MT-ATP8",
    }:
        return "protein"
    if "upstream" in consequence or "control" in consequence or g in {"MT-CR", "MT-OLR"}:
        return "control"
    if g:
        return "unknown"
    return "intergenic"


def _hgvs_m(row: dict, pos: int, ref: str, alt: str) -> str:
    hgvs = _first(row, "HGVS_M")
    if hgvs and hgvs.startswith("m."):
        return hgvs
    if len(ref) == 1 and len(alt) == 1:
        return f"m.{pos}{ref}>{alt}"
    # VCF anchors indels at the preceding base, so the changed bases
    # start at pos+1 for a deletion / are inserted after pos.
    if len(ref) > len(alt) and alt and ref.startswith(alt):
        del_seq = ref[len(alt):]
        s = pos + len(alt)
        e = s + len(del_seq) - 1
        return f"m.{s}del" if s == e else f"m.{s}_{e}del"
    if len(alt) > len(ref) and ref and alt.startswith(ref):
        ins_seq = alt[len(ref):]
        a = pos
        b = pos + len(ins_seq) - 1
        return f"m.{a}_{b}dup" if a != b else f"m.{a}dup"
    return f"m.{pos}{ref}>{alt}"


def _clean_hgvs_p(s: str) -> str:
    s = _first({"v": s}, "v")
    if not s:
        return ""
    s = unquote(s)
    if ":" in s:
        s = s.rsplit(":", 1)[1]
    return s


def load_mito_tsv(
    tsv_path: Path,
    *,
    pheno_by_gene: dict[str, float] | None = None,
) -> tuple[dict[str, dict], dict[str, list[str]]]:
    """Read mito.annotated.tsv → ({id: variant}, {tier: [ids]}).

    `id` is chrM-{pos}-{ref}-{alt} (distinct from SNV/CNV ids, so the
    one flat state.reports.{status,edits} namespace stays collision-free).
    Only disease-relevant variants make it into the returned dicts.
    """
    pheno_by_gene = pheno_by_gene or {}
    variants: dict[str, dict] = {}
    categories: dict[str, list[str]] = {t: [] for t in MITO_TIERS}

    if not tsv_path or not tsv_path.exists():
        return variants, categories

    with tsv_path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pos = _to_int(row.get("POS"))
            ref = (row.get("REF") or "").strip()
            alt = (row.get("ALT") or "").strip()
            if pos is None or not ref or not alt:
                continue
            filt = _first(row, "FILTER")
            # Keep PASS only — every non-PASS Mutect2-mito flag is a
            # reason to distrust the call (NUMT / blacklist artefacts,
            # weak/low-base-qual noise, contamination, strand bias, …).
            if filt and filt not in ("PASS", "."):
                continue
            gene = _first(row, "GENE", "MITOMAP_LOCUS")
            locus_type = _derive_locus_type(row, gene)
            status = _first(row, "MITOMAP_STATUS")
            mitotip = _first(row, "MITOTIP_SCORE", "MITOMAP_MITOTIP")
            disease = _first(row, "MITOMAP_DISEASE")
            pathogenic = _is_pathogenic(status, mitotip)
            # Only keep disease-relevant variants: pathogenic, or with
            # any MITOMAP disease association. Everything else (the bulk
            # — polymorphisms / haplogroup variants) is dropped.
            if not pathogenic and not disease:
                continue

            vid = f"chrM-{pos}-{ref}-{alt}"
            het = _to_float(_first(row, "HETEROPLASMY", "AF_SAMPLE"))
            pheno = pheno_by_gene.get(gene)
            in_panel = bool(pheno and pheno > 0)
            refs_raw = _first(row, "MITOMAP_REFS")
            plasmy = _first(row, "MITOMAP_PLASMY")
            if not plasmy:
                homo = _first(row, "MITOMAP_HOMO")
                hetero = _first(row, "MITOMAP_HETERO")
                if homo or hetero:
                    plasmy = f"{homo or '-'}/{hetero or '-'}"

            v = {
                "id":            vid,
                "CHROM":         (row.get("CHROM") or "chrM").strip(),
                "POS":           pos,
                "REF":           ref,
                "ALT":           alt,
                "HGVS_M":        _hgvs_m(row, pos, ref, alt),
                "HGVS_P":        _clean_hgvs_p(row.get("HGVS_P") or ""),
                "gene_symbol":   gene,
                "locus_type":    locus_type,
                "consequence":   _first(row, "CONSEQUENCE"),
                "aa_change":     _first(row, "AA_CHANGE") or _clean_hgvs_p(row.get("HGVS_P") or ""),
                "heteroplasmy":  het,                       # 0-1 fraction
                "AD":            _first(row, "AD"),
                "depth":         _to_int(_first(row, "DEPTH", "DP")),
                "filter":        filt,
                "has_filter":    bool(filt),
                "TLOD":          _to_float(row.get("TLOD")),
                "mitomap_disease": disease,
                "mitomap_status":  status,
                "mitomap_plasmy":  plasmy,
                "mitomap_gb_freq": _first(row, "MITOMAP_GB_FREQ"),
                "mitomap_gb_seqs": _first(row, "MITOMAP_GB_SEQS"),
                "gnomad_mito_af":  _first(row, "GNOMAD_MITO_AF"),
                "mitomap_refs":    refs_raw,
                "mitotip_score":   mitotip,
                "mitomap_allele":  _first(row, "MITOMAP_ALLELE"),
                "pheno_score":     round(pheno, 2) if pheno is not None else None,
                "in_panel":        in_panel,
                "pathogenic":      pathogenic,
            }
            # ClinVar — runtime lookup against the chrM-only ClinVar
            # VCF (decoupled from the 三級分析 pipeline). Field names
            # match the SNV adapter so the UI / docx_export pick them
            # up uniformly; absent variants get blanks.
            cv = clinvar_mito.lookup(pos, ref, alt)
            v["CLNSIG"]        = _first(row, "CLINVAR_SIG") or cv.get("CLNSIG", "")
            v["clinvar_stars"] = cv.get("stars", 0)
            v["clinvar_dn"]    = _first(row, "CLINVAR_DN") or cv.get("CLNDN", "")
            v["CLNSIGCONF"]    = _first(row, "CLINVAR_SIGCONF") or cv.get("CLNSIGCONF", "")
            variants[vid] = v
            categories["MITO-1" if pathogenic else "MITO-2"].append(vid)

    def _relevance_key(vid: str) -> tuple:
        v = variants[vid]
        het = v.get("heteroplasmy")
        return (
            _status_rank(v.get("mitomap_status", "")),
            _mitotip_rank(v.get("mitotip_score", "")),
            0 if v.get("in_panel") else 1,
            -_refs_count(v.get("mitomap_refs", "")),
            -(het if het is not None else -1.0),
            v.get("POS") or 0,
        )
    for t in categories:
        categories[t].sort(key=_relevance_key)

    return variants, categories
