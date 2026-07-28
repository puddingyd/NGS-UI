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
  2. Only primary chrM/MT records with POS/REF/ALT are emitted. New v3.2
     pipeline mito annotation uses ClinVar/gnomAD-mito, with MITOMAP
     columns optionally appended by the NGS-UI worker after copy.

Tiers:
    MITO-1  Pathogenic   — ClinVar P/LP or MITOMAP confirmed/pathogenic
    MITO-2  Rare/reported — not ClinVar B/LB, and either not observed /
                            rare in gnomAD-mito or reported in
                            ClinVar/MITOMAP
    MITO-3  Other        — every other PASS mtDNA variant

Within tiers, sort by in-panel, heteroplasmy, depth, and position so the
most clinically useful rows float upward while preserving all variants.
"""
from __future__ import annotations

import csv
import re
from urllib.parse import unquote
from pathlib import Path

from ..services import clinvar_mito, inhouse_af_mito, panel_deadzone

MITO_TIERS = ["MITO-1", "MITO-2", "MITO-3"]

_PATHO_STATUS_RE = re.compile(
    r"\bCfrm\b|\bConfirmed\b|\[(?:L?P)\]|\bLikely[_ ]pathogenic\b|\bPathogenic\b",
    re.I,
)
_PATHO_MITOTIP   = {"pathogenic", "likely pathogenic"}
_RARE_AF_CUTOFF = 0.01
_CLINVAR_DISEASE_SKIP = {
    "not_specified",
    "not provided",
    "not_provided",
    "see cases",
    "see_cases",
}


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


def _is_clinvar_pathogenic(sig: str) -> bool:
    s = (sig or "").strip().lower().replace("_", " ")
    if not s or s == ".":
        return False
    if "conflicting" in s:
        return False
    return (
        "pathogenic/likely pathogenic" in s
        or "likely pathogenic/pathogenic" in s
        or s in {"pathogenic", "likely pathogenic"}
    )


def _is_clinvar_benign(sig: str) -> bool:
    s = (sig or "").strip().lower().replace("_", " ")
    if not s or s == ".":
        return False
    if "conflicting" in s:
        return False
    return (
        "benign/likely benign" in s
        or "likely benign/benign" in s
        or s in {"benign", "likely benign"}
    )


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


def _clinvar_diseases(raw: str) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for part in re.split(r"[&|]", raw or ""):
        name = unquote(part).strip().strip(";:,")
        if not name:
            continue
        name = re.sub(r"\s+", " ", name.replace("_", " ")).strip()
        key = name.lower()
        if key in _CLINVAR_DISEASE_SKIP or key in seen:
            continue
        seen.add(key)
        out.append(name)
    return out


def _gnomad_af_values(row: dict) -> list[float]:
    vals: list[float] = []
    for key in ("GNOMAD_MITO_AF_HOM", "GNOMAD_MITO_AF_HET", "GNOMAD_MITO_AF"):
        val = _to_float(_first(row, key))
        if val is not None:
            vals.append(val)
    return vals


def load_mito_tsv(
    tsv_path: Path,
    *,
    pheno_by_gene: dict[str, float] | None = None,
) -> tuple[dict[str, dict], dict[str, list[str]]]:
    """Read mito.annotated.tsv → ({id: variant}, {tier: [ids]}).

    `id` is chrM-{pos}-{ref}-{alt} (distinct from SNV/CNV ids, so the
    one flat state.reports.{status,edits} namespace stays collision-free).
    All PASS mtDNA variants make it into the returned dicts.
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
            gene = panel_deadzone.canonical_panel_gene_symbol(_first(row, "GENE", "MITOMAP_LOCUS"))
            locus_type = _derive_locus_type(row, gene)
            status = _first(row, "MITOMAP_STATUS")
            mitotip = _first(row, "MITOTIP_SCORE", "MITOMAP_MITOTIP")
            disease = _first(row, "MITOMAP_DISEASE")
            pathogenic = _is_pathogenic(status, mitotip)

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

            cv = clinvar_mito.lookup(pos, ref, alt)
            inhouse = inhouse_af_mito.lookup(pos, ref, alt)
            clnsig = _first(row, "CLINVAR_SIG") or cv.get("CLNSIG", "")
            clinvar_dn = _first(row, "CLINVAR_DN") or cv.get("CLNDN", "")
            gnomad_afs = _gnomad_af_values(row)
            gnomad_max = max(gnomad_afs) if gnomad_afs else None
            clinvar_pathogenic = _is_clinvar_pathogenic(clnsig)
            clinvar_benign = _is_clinvar_benign(clnsig)
            reported = bool(_clinvar_diseases(clinvar_dn)) or bool(disease)
            rare_or_unobserved = gnomad_max is None or gnomad_max < _RARE_AF_CUTOFF
            rare = (not clinvar_pathogenic) and (not pathogenic) and (not clinvar_benign) and (
                rare_or_unobserved or reported
            )

            v = {
                "id":            vid,
                "CHROM":         (row.get("CHROM") or "chrM").strip(),
                "POS":           pos,
                "REF":           ref,
                "ALT":           alt,
                "HGVS_M":        _hgvs_m(row, pos, ref, alt),
                "HGVS_C":        _first(row, "HGVS_C"),
                "HGVS_P":        _clean_hgvs_p(row.get("HGVS_P") or ""),
                "gene_symbol":   gene,
                "locus_type":    locus_type,
                "impact":        _first(row, "IMPACT"),
                "biotype":       _first(row, "BIOTYPE"),
                "consequence":   _first(row, "CONSEQUENCE"),
                "aa_change":     _first(row, "AA_CHANGE") or _clean_hgvs_p(row.get("HGVS_P") or ""),
                "heteroplasmy":  het,                       # 0-1 fraction
                "genotype":      _first(row, "GENOTYPE"),
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
                "gnomad_mito_af_hom": _first(row, "GNOMAD_MITO_AF_HOM"),
                "gnomad_mito_af_het": _first(row, "GNOMAD_MITO_AF_HET"),
                "gnomad_mito_an":  _first(row, "GNOMAD_MITO_AN"),
                "gnomad_mito_af_max": gnomad_max,
                # In-house chrM values are carrier counts/frequency, not
                # diploid nuclear-allele arithmetic. They are context only:
                # tiering and ACMG classification remain unchanged.
                "inhouse_af":      inhouse.get("inhouse_af"),
                "inhouse_ac":      inhouse.get("inhouse_ac"),
                "inhouse_an":      inhouse.get("inhouse_an"),
                "inhouse_nhom":    inhouse.get("inhouse_nhom"),
                "inhouse_het_mt":  inhouse.get("inhouse_het_mt"),
                "mitomap_refs":    refs_raw,
                "mitotip_score":   mitotip,
                "mitomap_allele":  _first(row, "MITOMAP_ALLELE"),
                "clinvar_variation_id": _first(row, "CLINVAR_VARIATION_ID") or cv.get("VariationID", ""),
                "clinvar_diseases": _clinvar_diseases(clinvar_dn),
                "omim_ids":       _first(row, "OMIM_IDS"),
                "pipeline":       _first(row, "PIPELINE"),
                "pheno_score":     round(pheno, 2) if pheno is not None else None,
                "in_panel":        in_panel,
                "pathogenic":      clinvar_pathogenic or pathogenic,
                "clinvar_pathogenic": clinvar_pathogenic,
                "clinvar_benign":  clinvar_benign,
                "mitomap_pathogenic": pathogenic,
                "mitomap_reported": bool(disease),
                "rare_mito":       rare,
            }
            # ClinVar — runtime lookup against the chrM-only ClinVar
            # VCF (decoupled from the 三級分析 pipeline). Field names
            # match the SNV adapter so the UI / docx_export pick them
            # up uniformly; absent variants get blanks.
            v["CLNSIG"]        = clnsig
            v["clinvar_stars"] = cv.get("stars", 0)
            v["clinvar_dn"]    = clinvar_dn
            v["CLNSIGCONF"]    = _first(row, "CLINVAR_SIGCONF") or cv.get("CLNSIGCONF", "")
            variants[vid] = v
            if clinvar_pathogenic or pathogenic:
                categories["MITO-1"].append(vid)
            elif rare:
                categories["MITO-2"].append(vid)
            else:
                categories["MITO-3"].append(vid)

    def _relevance_key(vid: str) -> tuple:
        v = variants[vid]
        het = v.get("heteroplasmy")
        return (
            0 if v.get("in_panel") else 1,
            -(het if het is not None else -1.0),
            -(v.get("depth") or -1),
            v.get("gnomad_mito_af_max") if v.get("gnomad_mito_af_max") is not None else 9,
            v.get("POS") or 0,
        )
    for t in categories:
        categories[t].sort(key=_relevance_key)

    return variants, categories
