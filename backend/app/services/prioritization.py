"""In-house prioritization: variant + phenotype + inheritance (plan §5.3, §6, §7; phase 3).

For every SNV/Indel candidate of a case::

    Total = V′ + P(g, d) + C

* V′ — the ranking points the loader already uses (In-house incl. the ClinVar
  ≥2★ / VCEP adopt and the P-KNN LLR → GeneBe, minus the 院內常見 penalty),
  plus the variant-type prior (``type_prior``: synonymous, UTR and intronic
  variants with a negative splice prediction start lower than a missense,
  estimated from ClinVar); a recessive pair uses the combined points of plan §5.3
  (the weaker variant is the second allele, or benign with the second allele missed).
* P(g, d) — phenotype points of disease d in gene g from the phenotype engine
  (plan §4.8.7).  Genes without a ranked disease get the lowest value −6.1.  With no
  usable HPO term (panel only, no knowledge base) P is not counted and the
  case keeps the old ordering (``mode == "old"``).
* C — inheritance compatibility (plan §6.3, provisional): dominant / X-linked
  dominant / unknown 0; recessive homozygous 0, two different variants 0
  (phase unknown; user decision 2026-11-24, no deduction), a single heterozygote −3 (WES) / −4 (WGS); X-linked
  recessive male 0, female heterozygote −3 (two female hets pair like AR).

Each variant takes the best candidate it takes part in (plan §6.5): over the
gene's diseases and their inheritance modes, the highest total, then higher P,
then higher V′.  The disease inheritance comes from the platform OMIM data
(per phenotype MIM; Orphanet entries borrow the OMIM entries of their MONDO
group, otherwise the gene's OMIM inheritance).

Everything here is read-only and deterministic; nothing is written to disk.
"""
from __future__ import annotations

import json
import math
import re
import threading
from collections import OrderedDict
from dataclasses import dataclass, field

L0 = math.log(0.1 / 0.9)
K = 0.733                       # ln 2.08 — one ACMG point
SINGLE_HET = {"WES": -3.0, "WGS": -4.0}
PAIR_PHASE_UNKNOWN = 0.0       # two variants, phase unknown: no deduction (2026-11-24)
PAIR_MIN_GAIN = 0.05           # a pair must beat the single het by half a shown digit
XLR_FEMALE_HET = -3.0
P_FLOOR = -6.1                 # = pheno_engine.engine.P_MIN (lowest phenotype points)
VERSION = "prio-v4"
CANDIDATE_LIMIT = 200          # candidate rows sent to the browser (20 shown, "顯示更多")

# Variant-type prior (2026-11-25, user decision): ranking points added to V′ for
# variant types that are far less often pathogenic than a missense, applied only
# when the splice prediction is negative (max(SpliceAI, |Pangolin|) ≤ 0.1, the
# BP4 cut-off); a missing prediction, or any protein-coding transcript where the
# variant changes the protein, adds nothing.  Values = log_2.08 of the odds of
# P/LP vs (VUS + B/LB) relative to missense, ClinVar 2026-11-22, rare (AF < 0.001)
# germline ≥1★ small variants in 4,321 genes with a P/LP (gene bootstrap 95% CI in
# the plan §5.2.12).  5′UTR, 3′UTR and deep intronic variants (away from the coding
# region and its splice junctions) are pooled into one class (v10.76, user decision:
# the UTR classes alone have few P/LP).  It is a prior, so BP4 / BP7 still count on top of it.
TYPE_PRIOR = {"synonymous": -5.1, "intron_near": -1.6, "noncoding": -3.3}
TYPE_PRIOR_LABEL = {"synonymous": "同義變異", "intron_near": "近剪接位內含子（+3～+6／−3～−20）",
                    "noncoding": "非編碼區（UTR、內含子 ≥+7／≤−21）"}
_TYPE_POOL = {"5utr": "noncoding", "3utr": "noncoding", "intron_deep": "noncoding"}
TYPE_PRIOR_SPLICE_MAX = 0.1
TYPE_PRIOR_SOURCE = "ClinVar 2026-11-22"
_PROTEIN_TERMS = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained",
    "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "feature_elongation",
    "feature_truncation", "inframe_insertion", "inframe_deletion", "missense_variant",
    "protein_altering_variant", "incomplete_terminal_codon_variant", "coding_sequence_variant",
}
_SYNONYMOUS_TERMS = {"synonymous_variant", "stop_retained_variant", "start_retained_variant"}
_INTRONIC_TERMS = {"intron_variant", "splice_region_variant", "splice_polypyrimidine_tract_variant",
                   "splice_donor_region_variant", "splice_donor_5th_base_variant"}
_INTRON_OFFSET = re.compile(r"c\.[*-]?\d+([+-])(\d+)")

_MIM = re.compile(r"\((\d{6})\)")
_INH = re.compile(r"\(([^()]*)\)\s*$")


# ── small helpers ──────────────────────────────────────────────────────────
def _prob(v: float) -> float:
    """ACMG-scale points → probability (0 points ≈ 10%, one point = odds ×2.08)."""
    x = L0 + K * v
    return 1.0 / (1.0 + math.exp(-x)) if x >= 0 else math.exp(x) / (1.0 + math.exp(x))


def _points(p: float) -> float:
    p = min(max(p, 1e-300), 1 - 1e-15)
    return (math.log(p / (1.0 - p)) - L0) / K


def pair_points(v1: float, v2: float, *, single_c: float) -> float:
    """Combined points of two recessive heterozygotes, phase unknown (plan §5.3).

    Call the stronger variant a and the weaker b.  Either b is the second
    pathogenic allele — both alleles found, the points of a — or b is benign
    and the second allele was missed — a alone with the single-het deduction
    ``single_c``.  The two are weighted by the probability that b is
    pathogenic (user decision 2026-11-27), so a second het never lowers the
    gene below the single het, and a benign one adds almost nothing.
    """
    a, b = max(v1, v2), min(v1, v2)
    pb = _prob(b)
    return _points(pb * _prob(a) + (1.0 - pb) * _prob(a + single_c))


def _type_class(consequence: str, hgvs_c: str) -> str:
    """'coding' / 'synonymous' / '5utr' / '3utr' / 'intron_near' / 'intron_deep' / '' for one transcript."""
    terms = {t.strip() for t in re.split(r"[&,|]", str(consequence or "")) if t.strip()}
    if not terms:
        return ""
    if terms & _PROTEIN_TERMS:
        return "coding"
    if terms & _SYNONYMOUS_TERMS:
        return "synonymous"
    if "5_prime_UTR_variant" in terms:
        return "5utr"
    if "3_prime_UTR_variant" in terms:
        return "3utr"
    if terms & _INTRONIC_TERMS:
        m = _INTRON_OFFSET.search(str(hgvs_c or ""))
        if not m:
            return ""
        k = int(m.group(2))
        if m.group(1) == "+":
            return "intron_deep" if k >= 7 else "intron_near" if k >= 3 else "coding"
        return "intron_deep" if k >= 21 else "intron_near" if k >= 3 else "coding"
    return ""


def _number(value):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def type_prior(variant: dict) -> dict | None:
    """Variant-type prior for one loaded SNV/Indel, or None when it does not apply.

    Every protein-coding transcript is looked at (other transcripts only when
    there is none); the class with the smallest deduction wins, so a variant that
    changes the protein in any transcript gets nothing."""
    options = [o for o in (variant.get("transcript_options") or []) if isinstance(o, dict)]
    if not options:
        options = [variant]
    coding = [o for o in options if str(o.get("transcript_type") or "").lower() == "protein_coding"]
    classes = [_TYPE_POOL.get(c, c) for c in (
        _type_class(o.get("Consequence") or o.get("consequence") or "", o.get("HGVS_C") or o.get("hgvs_c") or "")
        for o in (coding or options))]
    if not classes or any(c not in TYPE_PRIOR for c in classes):
        return None
    cls = max(classes, key=lambda c: TYPE_PRIOR[c])
    scores = [x for x in (_number(variant.get("SpliceAI_score")), _number(variant.get("Pangolin_score")))
              if x is not None]
    if not scores:
        return None
    splice = max(abs(x) for x in scores)
    if splice > TYPE_PRIOR_SPLICE_MAX:
        return None
    return {"class": cls, "label": TYPE_PRIOR_LABEL[cls], "points": TYPE_PRIOR[cls],
            "splice": round(splice, 4), "source": TYPE_PRIOR_SOURCE}


def norm_zygosity(value) -> str:
    text = str(value or "").strip().lower()
    if text in {"ref", "hom_ref", "homref", "hom-ref", "reference", "0/0", "0|0"}:
        return "ref"
    if text.startswith("hemi"):
        return "hemi"
    if text.startswith("hom"):
        return "hom"
    if text.startswith("het"):
        return "het"
    return ""


def norm_sex(karyotype: str = "", reported: str = "") -> str:
    """XY / XX from the ploidy karyotype, else the reported sex; '' when unknown."""
    k = str(karyotype or "").strip().upper()
    if k in {"XY", "XX"}:
        return k
    s = str(reported or "").strip().upper()
    if s in {"M", "MALE", "男"}:
        return "XY"
    if s in {"F", "FEMALE", "女"}:
        return "XX"
    return ""


def norm_test_type(value) -> str:
    return "WES" if str(value or "").strip().upper() == "WES" else "WGS"


def moi_classes(codes) -> set[str]:
    """OMIM inheritance codes → {'dominant','ar','xlr','xld','other'}; empty when unknown."""
    out = set()
    for code in codes or ():
        c = str(code).strip().upper()
        if c in {"AD", "DD", "ISOL", "IC"}:
            out.add("dominant")
        elif c in {"AR", "DR"}:
            out.add("ar")
        elif c == "XLR":
            out.add("xlr")
        elif c == "XLD":
            out.add("xld")
        elif c == "XL":
            out.update({"xlr", "xld"})
        elif c in {"YL", "MT", "MI"}:
            out.add("other")
    return out


def _codes(text: str) -> set[str]:
    return {c for c in re.split(r"[,;/&\s]+", str(text or "").upper()) if c}


# ── OMIM inheritance per disease ───────────────────────────────────────────
_OMIM_LOCK = threading.Lock()
_OMIM_CACHE: dict = {"sig": None, "genes": {}}


def omim_disease_moi(gene: str) -> tuple[dict[str, set[str]], set[str]]:
    """({'OMIM:xxxxxx': classes}, gene-level classes) from the platform OMIM data."""
    try:
        from . import omim_store
        sig = omim_store.cache_signature()
    except Exception:  # noqa: BLE001
        return {}, set()
    with _OMIM_LOCK:
        if _OMIM_CACHE["sig"] != sig:
            _OMIM_CACHE["sig"], _OMIM_CACHE["genes"] = sig, {}
        hit = _OMIM_CACHE["genes"].get(gene)
    if hit is not None:
        return hit
    per: dict[str, set[str]] = {}
    union: set[str] = set()
    try:
        row = omim_store.lookup(gene=gene) or {}
    except Exception:  # noqa: BLE001
        row = {}
    for line in str(row.get("OMIM_disease") or "").splitlines():
        text = line.strip()
        m = _MIM.search(text)
        inh = _INH.search(text)
        if not text or not m or not inh or inh.group(1).isdigit():
            continue
        cls = moi_classes(_codes(inh.group(1)))
        if cls:
            per.setdefault(f"OMIM:{m.group(1)}", set()).update(cls)
            union |= cls
    union |= moi_classes(_codes(row.get("Inheritance") or ""))
    result = (per, union)
    with _OMIM_LOCK:
        _OMIM_CACHE["genes"][gene] = result
    return result


def canonical_gene(gene: str) -> str:
    gene = str(gene or "").strip()
    if not gene:
        return ""
    try:
        from . import panel_deadzone
        return panel_deadzone.canonical_panel_gene_symbol(gene) or gene
    except Exception:  # noqa: BLE001
        return gene


# ── phenotype result cache ─────────────────────────────────────────────────
_PHENO_LOCK = threading.Lock()
_PHENO_CACHE: "OrderedDict[tuple, object]" = OrderedDict()
_PHENO_MAX = 16


def hpo_items(hpo_list) -> list[dict]:
    """Only HPO terms (panels never enter the phenotype engine)."""
    out = []
    for entry in hpo_list or []:
        if isinstance(entry, str):
            hid, w = entry, 1
        elif isinstance(entry, dict):
            hid = entry.get("phenotype") or entry.get("hpo_id") or ""
            w = entry.get("weight", 1)
        else:
            continue
        hid = str(hid).strip()
        if hid.upper().startswith("HP:"):
            out.append({"phenotype": hid, "weight": 1 if w in (None, "") else w})
    return out


def pheno_result(hpo_list, engine=None):
    """(engine, Result) for the case's HPO terms, cached; (None, None) without a knowledge base."""
    items = hpo_items(hpo_list)
    if not items:
        return None, None
    if engine is None:
        try:
            from .pheno_engine.engine import load_engine
            engine = load_engine()
        except Exception:  # noqa: BLE001  (numpy missing, unreadable snapshot)
            return None, None
    if engine is None:
        return None, None
    key = (id(engine), engine.manifest_sha256, json.dumps(items, sort_keys=True, default=str))
    with _PHENO_LOCK:
        hit = _PHENO_CACHE.get(key)
        if hit is not None:
            _PHENO_CACHE.move_to_end(key)
            return engine, hit
    result = engine.score(items)
    with _PHENO_LOCK:
        _PHENO_CACHE[key] = result
        while len(_PHENO_CACHE) > _PHENO_MAX:
            _PHENO_CACHE.popitem(last=False)
    return engine, result


# ── candidate composition (plan §6.5) ──────────────────────────────────────
@dataclass
class Var:
    id: str
    v: float
    zygosity: str = ""
    somatic: bool = False
    v_missing: bool = False
    type_prior: dict | None = None


@dataclass
class Disease:
    id: str
    name: str = ""
    P: float | None = None
    moi: set = field(default_factory=set)
    moi_source: str = ""
    rank: float | None = None


def _options(var: Var, cls: str, *, sex: str, test_type: str, partner: Var | None):
    """[(v_combo, C, mode, partner_id)] for one variant under one inheritance class."""
    if var.somatic:
        return [(var.v, 0.0, "somatic", "")]
    if cls in {"dominant", "xld", "other", "unknown"}:
        return [(var.v, 0.0, cls, "")]
    if cls == "xlr" and sex == "XY":
        return [(var.v, 0.0, "hemi", "")]
    # ar, or xlr in a female / unknown sex
    if var.zygosity in {"hom", "hemi"}:
        return [(var.v, 0.0, "hemi" if cls == "xlr" and var.zygosity == "hemi" else "hom", "")]
    single_c = SINGLE_HET[test_type] if cls == "ar" else XLR_FEMALE_HET
    out = [(var.v, single_c, "single" if cls == "ar" else "xlr_female_single", "")]
    if partner is not None:
        v_pair = pair_points(var.v, partner.v, single_c=single_c)
        # a benign partner barely moves the score: keep the single het (the
        # partner still shows under the candidate row as 同基因另一個變異)
        if v_pair + PAIR_PHASE_UNKNOWN - (max(var.v, partner.v) + single_c) >= PAIR_MIN_GAIN:
            out.append((v_pair, PAIR_PHASE_UNKNOWN, "pair", partner.id))
    return out


def _hets(variants: list[Var]) -> list[Var]:
    """Possible second alleles, strongest first (somatic, homozygous and ref calls never pair)."""
    return sorted((x for x in variants if not x.somatic and x.zygosity not in {"hom", "hemi", "ref"}),
                  key=lambda x: (-x.v, x.id))


def _enumerate(variants: list[Var], diseases: list[Disease], *, sex: str, test_type: str):
    """Every (variant, candidate) pair: each disease × inheritance class × genotype option."""
    test_type = norm_test_type(test_type)
    hets = _hets(variants)
    by_id = {x.id: x for x in variants}
    for var in variants:
        partner = None if var.zygosity == "ref" else next((x for x in hets if x.id != var.id), None)
        for d_order, dis in enumerate(diseases):
            for cls in sorted(dis.moi) or ["unknown"]:
                for v_combo, c, mode, partner_id in _options(var, cls, sex=sex, test_type=test_type, partner=partner):
                    p = dis.P
                    total = v_combo + (p or 0.0) + c
                    key = (round(total, 6), -1e9 if p is None else round(p, 6), round(v_combo, 6), -d_order)
                    cand = {"total": total, "v": v_combo, "p": p, "c": c, "mode": mode, "moi_class": cls,
                            "disease_id": dis.id, "disease_name": dis.name, "disease_rank": dis.rank,
                            "moi": sorted(dis.moi), "moi_source": dis.moi_source,
                            "partner": partner_id,
                            "partner_v": by_id[partner_id].v if partner_id else None,
                            "v_single": var.v, "v_missing": var.v_missing, "type_prior": var.type_prior}
                    yield var, key, cand


def evaluate_gene(variants: list[Var], diseases: list[Disease], *, sex: str, test_type: str) -> dict[str, dict]:
    """Best candidate per variant of one gene."""
    best: dict[str, tuple] = {}
    for var, key, cand in _enumerate(variants, diseases, sex=sex, test_type=test_type):
        if var.id not in best or key > best[var.id][0]:
            best[var.id] = (key, cand)
    return {vid: item[1] for vid, item in best.items()}


def gene_candidates(variants: list[Var], diseases: list[Disease], *, sex: str, test_type: str,
                    others: int = 5) -> dict | None:
    """The gene's candidate diagnosis (plan §6.5): its best combination, plus the next few.

    A recessive pair is one candidate (both variants), not two.  Reference
    calls never form a candidate.  None when the gene has no usable variant."""
    seen: dict[tuple, tuple] = {}
    for var, key, cand in _enumerate(variants, diseases, sex=sex, test_type=test_type):
        if var.zygosity == "ref":
            continue
        ids = tuple(sorted(filter(None, (var.id, cand["partner"]))))
        ident = (cand["disease_id"], cand["moi_class"], cand["mode"], ids)
        if ident not in seen:
            item = dict(cand)
            item.pop("partner", None)
            item["v_each"] = {var.id: var.v, **({cand["partner"]: cand["partner_v"]} if cand["partner"] else {})}
            # stronger variant (higher V′) first, so a pair reads LP + VUS
            item["variants"] = sorted(item["v_each"], key=lambda i: (-item["v_each"][i], i))
            seen[ident] = (key, item)
    if not seen:
        return None
    # highest total → P → V′ → disease order, then variant IDs for a fixed order
    ranked = [item for _, item in sorted(seen.values(),
                                         key=lambda x: (tuple(-k for k in x[0]), ",".join(sorted(x[1]["variants"]))))]
    best = dict(ranked[0])
    best["others"] = ranked[1:1 + others]
    best["n_combinations"] = len(ranked)
    best["same_gene"] = same_gene_hets(variants, best)
    return best


SECOND_ALLELE_MODES = {"single", "pair", "xlr_female_single"}
SAME_GENE_LIMIT = 5


def same_gene_hets(variants: list[Var], cand: dict, limit: int = SAME_GENE_LIMIT) -> list[str]:
    """Other heterozygous variants of the gene, strongest first, for a recessive row.

    Shown under the candidate so a possible second allele is never hidden,
    whatever its score (a benign one included).  Empty for dominant / hom rows."""
    if cand.get("mode") not in SECOND_ALLELE_MODES:
        return []
    used = set(cand.get("variants") or [])
    return [x.id for x in _hets(variants) if x.id not in used][:limit]


# ── case level ─────────────────────────────────────────────────────────────
def ranking_points(variant: dict) -> tuple[float, bool]:
    """V′ for a loaded variant: the loader's ranking points (missing → 0, flagged)."""
    for key in ("prio_v_points", "effective_rank_score", "inhouse_rank_score", "ACMG_score", "genebe_acmg_score"):
        value = variant.get(key)
        try:
            if value not in (None, ""):
                return float(value), False
        except (TypeError, ValueError):
            continue
    return 0.0, True


def _gene_diseases(gene: str, engine, result) -> tuple[list[Disease], dict]:
    per, union = omim_disease_moi(gene) if gene else ({}, set())
    info: dict = {}
    if result is not None and gene:
        gi = engine.gid.get(gene)
        if gi is not None:
            g = result.gene(gene) or {}
            info = {"gene_rank": g.get("rank"), "gene_rank_mid": g.get("rank_mid"), "gene_points": g.get("P"),
                    "lead": g.get("lead"), "near_tie": g.get("near_tie")}
            out = []
            groups = engine.groups
            for row in result.gene_disease_points(gene):
                did = row["disease_id"]
                moi, src = set(per.get(did, set())), "OMIM" if did in per else ""
                if not moi:
                    grp = groups.get(did)
                    for other, other_moi in per.items():
                        if grp and groups.get(other) == grp:
                            moi |= other_moi
                    src = "OMIM（MONDO 同組）" if moi else ""
                if not moi and union:
                    moi, src = set(union), "OMIM（基因）"
                out.append(Disease(id=did, name=engine.kb.disease_names[engine.did[did]], P=float(row["P"]),
                                   moi=moi, moi_source=src, rank=row.get("group_rank")))
            return out, info
        # no ranked disease: floor
        info = {"gene_rank": None, "gene_points": P_FLOOR, "no_disease": True}
        return [Disease(id="", name="", P=P_FLOOR, moi=set(union), moi_source="OMIM（基因）" if union else "")], info
    if result is not None:
        return [Disease(id="", P=P_FLOOR)], {"gene_points": P_FLOOR, "no_disease": True}
    # no phenotype: inheritance only
    if per:
        return [Disease(id=did, moi=set(moi), moi_source="OMIM") for did, moi in sorted(per.items())], info
    return [Disease(id="", moi=set(union), moi_source="OMIM（基因）" if union else "")], info


def _context(hpo, engine=None, old_label: str = "Score (old)"):
    """(case info, engine, phenotype result or None when P is not counted)."""
    info = {"version": VERSION, "mode": "old", "reason": "", "snapshot": "", "terms": 0}
    try:
        engine, result = pheno_result(hpo, engine=engine)
    except Exception as exc:  # noqa: BLE001
        engine, result = None, None
        info["reason"] = f"表型引擎無法計算（{type(exc).__name__}）"
    if result is not None and not result.empty:
        info.update(mode="prio", snapshot=engine.snapshot, terms=len(result.q))
        return info, engine, result
    if not info["reason"]:
        if not hpo_items(hpo):
            info["reason"] = f"沒有 HPO（只有 panel 或未設定），表型分數不計，排序沿用 {old_label}"
        elif engine is None:
            info["reason"] = f"尚未啟用表型知識庫，排序沿用 {old_label}"
        else:
            info["reason"] = f"HPO 都被移除（非表型或無效），排序沿用 {old_label}"
    return info, engine, None


def annotate(variants: dict, *, hpo, sex: str = "", test_type: str = "WES", rank: bool = True,
             engine=None) -> dict:
    """Add ``variant['prio']`` (and ``old_rank``) to every SNV/Indel variant; return case info.

    ``mode`` is "prio" when phenotype points are available, else "old" (sort by
    Score (old)).  Ranks are competition ranks over all variants passed in.
    """
    info = {"version": VERSION, "mode": "old", "reason": "", "snapshot": "", "terms": 0}
    if not variants:
        return info
    info, engine, result = _context(hpo, engine)
    sex = norm_sex(sex)
    test_type = norm_test_type(test_type)
    # the browser pairs CNV deletions with SNVs (plan §5.4) using the same deductions
    info["single_het"], info["xlr_female_het"] = SINGLE_HET[test_type], XLR_FEMALE_HET
    by_gene: dict[str, list[Var]] = {}
    for vid, variant in variants.items():
        if not isinstance(variant, dict):
            continue
        v, missing = ranking_points(variant)
        tp = type_prior(variant)
        if tp is not None:
            v += tp["points"]
        gene = canonical_gene(variant.get("gene_symbol") or "")
        by_gene.setdefault(gene, []).append(
            Var(id=vid, v=v, zygosity=norm_zygosity(variant.get("zygosity")),
                somatic=bool(variant.get("somatic")), v_missing=missing, type_prior=tp))
    genes = []
    for gene in sorted(by_gene):
        diseases, ginfo = _gene_diseases(gene, engine, result)
        best = evaluate_gene(by_gene[gene], diseases, sex=sex, test_type=test_type)
        for vid, cand in best.items():
            cand.update(ginfo)
            cand["gene"] = gene
            cand["pheno_counted"] = result is not None
            if result is None:
                cand["p"] = None
            variants[vid]["prio"] = cand
            variants[vid]["prio_total"] = round(cand["total"], 3)
        if rank and result is not None and gene:
            cand = gene_candidates(by_gene[gene], diseases, sex=sex, test_type=test_type)
            if cand is not None:
                cand.update(ginfo, gene=gene)
                genes.append(cand)
    if rank:
        _rank(variants)
        if result is not None:
            info["candidates"], info["candidate_count"] = rank_candidates(genes)
    return info


def rank_candidates(genes: list[dict], limit: int = CANDIDATE_LIMIT) -> tuple[list[dict], int]:
    """One row per gene, ordered as plan §7.4 (total, P, V′, gene); competition ranks on the shown total."""
    genes = sorted(genes, key=lambda g: (-round(g["total"], 3), -(g["p"] if g["p"] is not None else -1e9),
                                         -round(g["v"], 3), g["gene"]))
    prev, pos = None, 0
    for i, g in enumerate(genes):
        shown = round(g["total"], 1)
        if prev is None or shown != prev:
            pos, prev = i + 1, shown
        g["rank"] = pos
    return genes[:limit], len(genes)


def sort_key(variant: dict, vid: str, mode: str):
    """Tier ordering: plan §7.4 for "prio", old total_score otherwise."""
    if mode == "prio":
        prio = variant.get("prio") or {}
        p = prio.get("p")
        return (-round(float(prio.get("total", -1e9)), 3), -(p if p is not None else -1e9),
                -round(float(prio.get("v", -1e9)), 3), str(prio.get("gene") or ""), vid)
    try:
        t = float(variant.get("total_score"))
    except (TypeError, ValueError):
        t = -math.inf
    return (-t, vid)


def _competition(values: list[tuple[float, str]]) -> dict[str, int]:
    ranks, prev, pos = {}, None, 0
    for i, (val, vid) in enumerate(sorted(values, key=lambda x: (-x[0], x[1]))):
        if prev is None or val != prev:
            pos, prev = i + 1, val
        ranks[vid] = pos
    return ranks


def _rank(variants: dict) -> None:
    prio_vals, old_vals = [], []
    for vid, variant in variants.items():
        prio = variant.get("prio") if isinstance(variant, dict) else None
        if prio:
            prio_vals.append((round(prio["total"], 1), vid))
        try:
            old_vals.append((round(float(variant.get("total_score")), 1), vid))
        except (TypeError, ValueError, AttributeError):
            pass
    pr = _competition(prio_vals)
    orank = _competition(old_vals)
    for vid, variant in variants.items():
        if not isinstance(variant, dict):
            continue
        if vid in pr:
            variant["prio"]["rank"] = pr[vid]
        variant["old_rank"] = orank.get(vid)


def resort(categories: dict, variants: dict, mode: str, tiers=("1A", "1B", "1C", "2")) -> None:
    for tier in tiers:
        ids = categories.get(tier)
        if ids:
            categories[tier] = sorted(set(ids), key=lambda i: sort_key(variants.get(i, {}), i, mode))


# ── CNV / SV events (plan §5.4, phase 5b; 2026-11-24 user decision) ─────────
# Scores are shown on the cards and order the CNV/SV tiers.  The candidate
# list stays SNV-based: a CNV deletion enters it only as the second allele of
# a recessive gene with an SNV (``ar_genes`` below; paired in the browser).
#
# V′ from the AnnotSV ranking score (ACMG CNV points) — piecewise linear:
#   ≥0.99 → +10, 0.90–0.98 → +6…+9, 0…0.89 → 0…+5,
#   −0.89…0 → −1…0, −0.90…−0.98 → −1…−6, ≤−0.99 → −7.
CNV_KNOTS = ((-0.99, -7.0), (-0.98, -6.0), (-0.90, -1.0), (0.0, 0.0),
             (0.89, 5.0), (0.90, 6.0), (0.98, 9.0), (0.99, 10.0))
CNV_AR_GENE_LIMIT = 200         # recessive genes kept per deletion for SNV pairing
X_PAR = ((10_001, 2_781_479), (155_701_383, 156_030_895))


def annotsv_points(score) -> float | None:
    """V′ of a CNV/SV from its AnnotSV ranking score (None when missing)."""
    try:
        s = float(score)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(s):
        return None
    if s <= CNV_KNOTS[0][0]:
        return CNV_KNOTS[0][1]
    if s >= CNV_KNOTS[-1][0]:
        return CNV_KNOTS[-1][1]
    for (x0, y0), (x1, y1) in zip(CNV_KNOTS, CNV_KNOTS[1:]):
        if x0 <= s <= x1:
            return round(y0 + (y1 - y0) * (s - x0) / (x1 - x0), 3)
    return 0.0


def _x_non_par(chrom, start, end) -> bool:
    c = str(chrom or "").lower().removeprefix("chr")
    if c != "x":
        return False
    try:
        a, b = int(start), int(end if end not in (None, "") else start)
    except (TypeError, ValueError):
        return True
    return not any(a <= pe and b >= ps for ps, pe in X_PAR)


def event_zygosity(event: dict, sex: str) -> str:
    """hom (CN 0 / homozygous deletion), hemi (male chrX outside PAR) or het."""
    sv_type = str(event.get("sv_type") or "").upper()
    cn = event.get("copy_number")
    if sv_type == "DEL" and (cn == 0 or str(event.get("zygosity") or "").lower().startswith("hom")):
        return "hom"
    if sex == "XY" and _x_non_par(event.get("CHROM"), event.get("POS"), event.get("END")):
        return "hemi"
    return "het"


def _better(key_a, key_b) -> bool:
    return key_b is None or key_a > key_b


def annotate_events(events: dict, *, hpo, sex: str = "", test_type: str = "WES", source: str = "cnv",
                    engine=None) -> dict:
    """``event['prio']`` for every CNV/SV event; ranks are competition ranks per source.

    The genes scored are those the event affects (``_prio_genes`` from the
    adapter: coding / splice / unclear impact; pure UTR or intronic genes are
    left out).  An event without such a gene gets the lowest phenotype points.
    """
    info, engine, result = _context(hpo, engine, old_label="原本 CNV/SV 排序")
    info["source"] = source
    sex = norm_sex(sex)
    test_type = norm_test_type(test_type)
    gene_cache: dict[str, tuple] = {}
    totals = []
    for eid, ev in events.items():
        if not isinstance(ev, dict):
            continue
        genes = ev.pop("_prio_genes", None) or []
        score = ev.get("ranking_score")
        v = annotsv_points(score)
        missing = v is None
        v = 0.0 if v is None else v
        zyg = event_zygosity(ev, sex)
        best, best_key, ar_genes = None, None, []
        for order, (gene, impact) in enumerate(genes):
            if gene not in gene_cache:
                gene_cache[gene] = _gene_diseases(gene, engine, result)
            diseases, ginfo = gene_cache[gene]
            cand = evaluate_gene([Var(id=eid, v=v, zygosity=zyg, v_missing=missing)], diseases,
                                 sex=sex, test_type=test_type)[eid]
            p = cand["p"]
            key = (round(cand["total"], 6), -1e9 if p is None else round(p, 6), -order)
            if _better(key, best_key):
                best, best_key = dict(cand, gene=gene, impact=impact, **ginfo), key
            # a heterozygous deletion can be the second allele of a recessive SNV
            if (source == "cnv" and result is not None and zyg == "het" and impact == "functional"
                    and str(ev.get("sv_type") or "").upper() == "DEL" and len(ar_genes) < CNV_AR_GENE_LIMIT):
                rec = [{"disease_id": d.id, "disease_name": d.name, "p": d.P, "disease_rank": d.rank,
                        "moi_class": cls, "moi": sorted(d.moi), "moi_source": d.moi_source}
                       for d in diseases for cls in ("ar", "xlr") if cls in d.moi
                       and (cls == "ar" or sex != "XY") and d.id]
                if rec:
                    ar_genes.append({"gene": gene, "diseases": rec, **ginfo})
        if best is None:
            p = P_FLOOR if result is not None else None
            best = {"total": v + (p or 0.0), "v": v, "p": p, "c": 0.0, "mode": "no_gene", "moi_class": "",
                    "disease_id": "", "disease_name": "", "disease_rank": None, "moi": [], "moi_source": "",
                    "partner": "", "partner_v": None, "v_single": v, "v_missing": missing,
                    "gene": "", "impact": "", "no_disease": True, "gene_points": p}
        best.update(annotsv_score=score, zygosity=zyg, n_genes=len(genes), pheno_counted=result is not None)
        if result is None:
            best["p"] = None
        if ar_genes:
            best["ar_genes"] = ar_genes
        ev["prio"] = best
        ev["prio_total"] = round(best["total"], 3)
        totals.append((round(best["total"], 1), eid))
    ranks = _competition(totals)
    for eid, r in ranks.items():
        events[eid]["prio"]["rank"] = r
    return info


def resort_events(categories: dict, events: dict, mode: str) -> None:
    if mode != "prio":
        return
    for tier, ids in categories.items():
        categories[tier] = sorted(ids, key=lambda i: sort_key(events.get(i, {}), i, mode))


# ── mitochondrial variants (plan §5.5, phase 5c) ───────────────────────────
MITO_COMMON_AF = 0.01


def mito_points(variant: dict) -> tuple[float, str]:
    """V′ of an mtDNA variant (ClinVar + MITOMAP + gnomAD-mito); heteroplasmy is display only."""
    clinvar_plp = bool(variant.get("clinvar_pathogenic"))
    mitomap_p = bool(variant.get("mitomap_pathogenic"))
    if variant.get("clinvar_benign"):
        return -7.0, "ClinVar B/LB"
    if clinvar_plp and mitomap_p:
        return 10.0, "ClinVar P/LP＋MITOMAP 確認致病"
    if clinvar_plp:
        return 7.0, "ClinVar P/LP"
    if mitomap_p:
        return 7.0, "MITOMAP 確認致病"
    af = variant.get("gnomad_mito_af_max")
    if isinstance(af, (int, float)) and af >= MITO_COMMON_AF:
        return -4.0, "gnomAD-mito ≥1%"
    sig = str(variant.get("CLNSIG") or "").lower()
    if "uncertain" in sig or "conflicting" in sig or variant.get("mitomap_reported"):
        return 2.0, "ClinVar VUS／MITOMAP 有報告"
    return 1.0, "罕見（gnomAD-mito <1%），無報告"


def annotate_mito(variants: dict, *, hpo, engine=None) -> dict:
    """``variant['prio']`` for every mtDNA variant: V′ + P(g, d), inheritance C = 0."""
    info, engine, result = _context(hpo, engine, old_label="原本粒線體排序")
    info["source"] = "mito"
    gene_cache: dict[str, tuple] = {}
    totals = []
    for vid, variant in variants.items():
        if not isinstance(variant, dict):
            continue
        v, basis = mito_points(variant)
        gene = canonical_gene(variant.get("gene_symbol") or "")
        if gene not in gene_cache:
            gene_cache[gene] = _gene_diseases(gene, engine, result)
        diseases, ginfo = gene_cache[gene]
        best = max(enumerate(diseases),
                   key=lambda x: (-1e9 if x[1].P is None else x[1].P, -x[0]))[1] if diseases else Disease(id="")
        p = best.P if result is not None else None
        cand = {"total": v + (p or 0.0), "v": v, "p": p, "c": 0.0, "mode": "mito", "moi_class": "other",
                "disease_id": best.id, "disease_name": best.name, "disease_rank": best.rank,
                "moi": sorted(best.moi), "moi_source": best.moi_source, "partner": "", "partner_v": None,
                "v_single": v, "v_missing": False, "v_basis": basis, "gene": gene,
                "pheno_counted": result is not None, **ginfo}
        variant["prio"] = cand
        variant["prio_total"] = round(cand["total"], 3)
        totals.append((round(cand["total"], 1), vid))
    for vid, r in _competition(totals).items():
        variants[vid]["prio"]["rank"] = r
    return info


# ── which patient HPO terms a candidate disease explains ───────────────────
MATCH_RELATIONS = ("same", "patient_specific", "patient_broader")
EXPLAIN_LIMIT = 300


def match_class(relation: str, net: float) -> str:
    """match = the disease is annotated with the term, a broader or a more specific one;
    partial = only a related (collateral) term that still scores above the disease baseline."""
    if relation in MATCH_RELATIONS:
        return "match"
    if relation == "collateral" and net > 0:
        return "partial"
    return "none"


def explain_diseases(hpo, disease_ids, *, engine=None) -> dict:
    """Per patient HPO term: the disease annotation it pairs with and how (plan §4.7)."""
    from .pheno_engine import engine as E
    engine, result = pheno_result(hpo, engine=engine)
    if result is None or result.empty:
        return {"terms": [], "diseases": {}, "snapshot": getattr(engine, "snapshot", "")}
    kb = engine.kb
    terms = [{"hpo_id": kb.terms[int(q)], "name": kb.names[int(q)], "weight": float(w)}
             for q, w in zip(result.q, result.w)]
    from .pheno_engine import lineage as L
    lin = L.lineage(engine)
    tid = engine.tid
    out = {}
    seen = set()
    for did in disease_ids:
        if did in seen or did not in engine.did or len(seen) >= EXPLAIN_LIMIT:
            continue
        seen.add(did)
        di = engine.did[did]
        rows = []
        for r in result.explain(di):
            steps = L.pair_steps(lin, tid[r["hpo_id"]], tid[r["annotation"]], r["relation"], tid[r["mica"]])
            none = r["relation"] == "none"
            rows.append({"hpo_id": r["hpo_id"], "annotation": r["annotation"], "annotation_name": r["annotation_name"],
                         "relation": r["relation"], "relation_zh": E.RELATION_ZH.get(r["relation"], r["relation"]),
                         "mica": "" if none else r["mica"], "mica_name": "" if none else r["mica_name"],
                         "ic_mica": round(float(r["ic_mica"]), 3), "s": round(float(r["s"]), 3),
                         "net": round(float(r["net"]), 3),
                         "steps_patient": steps["steps_patient"], "steps_annotation": steps["steps_annotation"],
                         "match": match_class(r["relation"], float(r["net"]))})
        out[did] = {"name": kb.disease_names[di], "n_annotations": int(kb.dir_ptr[di + 1] - kb.dir_ptr[di]),
                    "matched": sum(r["match"] == "match" for r in rows),
                    "partial": sum(r["match"] == "partial" for r in rows), "rows": rows}
    return {"terms": terms, "diseases": out, "snapshot": engine.snapshot}


def disease_graph(hpo, disease_id: str, *, all_annotations: bool = False, engine=None) -> dict | None:
    """HPO sub-graph linking the patient terms to one disease (plan: candidate list picture)."""
    from .pheno_engine import engine as E
    from .pheno_engine import lineage as L
    engine, result = pheno_result(hpo, engine=engine)
    if result is None or result.empty or disease_id not in engine.did:
        return None
    kb = engine.kb
    di = engine.did[disease_id]
    graph = L.disease_graph(engine, result, di, all_annotations=all_annotations)
    rows = {r["hpo_id"]: r for r in explain_diseases(hpo, [disease_id], engine=engine)["diseases"][disease_id]["rows"]}
    for pair in graph["pairs"]:
        r = rows.get(pair["hpo_id"], {})
        pair.update(relation_zh=E.RELATION_ZH.get(pair["relation"], pair["relation"]), match=r.get("match", "none"),
                    ic_mica=r.get("ic_mica"), s=r.get("s"), net=r.get("net"),
                    annotation_name=r.get("annotation_name", ""), mica_name=r.get("mica_name", ""))
    graph["terms"] = [{"hpo_id": kb.terms[int(q)], "name": kb.names[int(q)], "weight": float(w)}
                      for q, w in zip(result.q, result.w)]
    graph["snapshot"] = engine.snapshot
    return graph
