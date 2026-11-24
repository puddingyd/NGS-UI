"""In-house prioritization: variant + phenotype + inheritance (plan §5.3, §6, §7; phase 3).

For every SNV/Indel candidate of a case::

    Total = V′ + P(g, d) + C

* V′ — the ranking points the loader already uses (In-house incl. the ClinVar
  ≥2★ / VCEP adopt and the P-KNN LLR → GeneBe, minus the 院內常見 penalty);
  a recessive pair uses the combined points of plan §5.3.
* P(g, d) — phenotype points of disease d in gene g from the phenotype engine
  (plan §4.8.7).  Genes without a ranked disease get the floor −3.1.  With no
  usable HPO term (panel only, no knowledge base) P is not counted and the
  case keeps the old ordering (``mode == "old"``).
* C — inheritance compatibility (plan §6.3, provisional): dominant / X-linked
  dominant / unknown 0; recessive homozygous 0, two different variants −1
  (phase unknown), a single heterozygote −3 (WES) / −4 (WGS); X-linked
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
PAIR_PHASE_UNKNOWN = -1.0
XLR_FEMALE_HET = -3.0
P_FLOOR = -3.1
VERSION = "prio-v1"
CANDIDATE_LIMIT = 200          # candidate rows sent to the browser (20 shown, "顯示更多")

_MIM = re.compile(r"\((\d{6})\)")
_INH = re.compile(r"\(([^()]*)\)\s*$")


# ── small helpers ──────────────────────────────────────────────────────────
def pair_points(v1: float, v2: float) -> float:
    """Combined points of two recessive alleles (plan §5.3)."""
    p1 = 1.0 / (1.0 + math.exp(-(L0 + K * v1)))
    p2 = 1.0 / (1.0 + math.exp(-(L0 + K * v2)))
    p = min(max(p1 * p2, 1e-300), 1 - 1e-15)
    return (math.log(p / (1.0 - p)) - L0) / K


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
        out.append((pair_points(var.v, partner.v), PAIR_PHASE_UNKNOWN, "pair", partner.id))
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
                            "v_single": var.v, "v_missing": var.v_missing}
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
            item = dict(cand, variants=[var.id] + ([cand["partner"]] if cand["partner"] else []))
            item.pop("partner", None)
            item["v_each"] = {var.id: var.v, **({cand["partner"]: cand["partner_v"]} if cand["partner"] else {})}
            seen[ident] = (key, item)
    if not seen:
        return None
    # highest total → P → V′ → disease order, then variant IDs for a fixed order
    ranked = [item for _, item in sorted(seen.values(),
                                         key=lambda x: (tuple(-k for k in x[0]), ",".join(sorted(x[1]["variants"]))))]
    best = dict(ranked[0])
    best["others"] = ranked[1:1 + others]
    best["n_combinations"] = len(ranked)
    return best


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


def annotate(variants: dict, *, hpo, sex: str = "", test_type: str = "WES", rank: bool = True,
             engine=None) -> dict:
    """Add ``variant['prio']`` (and ``old_rank``) to every SNV/Indel variant; return case info.

    ``mode`` is "prio" when phenotype points are available, else "old" (sort by
    Score (old)).  Ranks are competition ranks over all variants passed in.
    """
    info = {"version": VERSION, "mode": "old", "reason": "", "snapshot": "", "terms": 0}
    if not variants:
        return info
    try:
        engine, result = pheno_result(hpo, engine=engine)
    except Exception as exc:  # noqa: BLE001
        engine, result = None, None
        info["reason"] = f"表型引擎無法計算（{type(exc).__name__}）"
    if result is not None and not result.empty:
        info.update(mode="prio", snapshot=engine.snapshot, terms=len(result.q))
    else:
        result = None if result is None or result.empty else result
        if not info["reason"]:
            if not hpo_items(hpo):
                info["reason"] = "沒有 HPO（只有 panel 或未設定），表型分數不計，排序沿用 Score (old)"
            elif engine is None:
                info["reason"] = "尚未啟用表型知識庫，排序沿用 Score (old)"
            else:
                info["reason"] = "HPO 都被移除（非表型或無效），排序沿用 Score (old)"
    sex = norm_sex(sex)
    test_type = norm_test_type(test_type)
    by_gene: dict[str, list[Var]] = {}
    for vid, variant in variants.items():
        if not isinstance(variant, dict):
            continue
        v, missing = ranking_points(variant)
        gene = canonical_gene(variant.get("gene_symbol") or "")
        by_gene.setdefault(gene, []).append(
            Var(id=vid, v=v, zygosity=norm_zygosity(variant.get("zygosity")),
                somatic=bool(variant.get("somatic")), v_missing=missing))
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
