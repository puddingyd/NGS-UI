"""NGS-UI ACMG rule engine (plan docs/research/NGS_Integrated_Prioritization_Plan_v3.md §5.2).

``evaluate_variant(rows, context)`` scores one genomic variant from its TSV
transcript rows (raw 03_acmg columns merged with the 08 overlay).  It is pure:
reference tables come from ``context.refs`` and thresholds from the versioned
rulebook, so the same inputs always give the same result.

Each criterion ends as ``met`` (with a strength), ``not_met`` or
``not_evaluated`` (data or reference table missing).  Only ``met`` criteria
carry points.  ClinVar never adds points: 3–4★ (or a ClinGen VCEP call) is
adopted, 2★ P/LP sets a floor of 6 points.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

from . import parse as P
from .refs import References

RULEBOOK_PATH = Path(__file__).with_name("rulebook_v1.json")

CRITERIA_ORDER = (
    "PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5", "BA1", "BS1", "BS2", "BS3", "BS4",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
)
DEFAULT_STRENGTH = {
    "PVS1": "very_strong", "PS1": "strong", "PM1": "moderate", "PM2": "moderate",
    "PM4": "moderate", "PM5": "moderate", "PP2": "supporting", "PP3": "supporting", "BA1": "stand_alone",
    "BS1": "strong", "BP3": "supporting", "BP4": "supporting", "BP7": "supporting",
}
STRENGTH_LABEL = {
    "supporting": "Supporting", "moderate": "Moderate", "strong": "Strong",
    "very_strong": "Very_strong", "stand_alone": "Stand-alone",
}
CLASS_RANGES = {
    "P": ("Pathogenic", 10, None),
    "LP": ("Likely pathogenic", 6, 9),
    "VUS": ("Uncertain significance", 0, 5),
    "LB": ("Likely benign", -6, -1),
    "B": ("Benign", None, -7),
}


@lru_cache(maxsize=4)
def load_rulebook(path: str | None = None) -> dict:
    with open(path or RULEBOOK_PATH, encoding="utf-8") as handle:
        return json.load(handle)


@dataclass
class Context:
    """Per-case inputs that do not live in the variant rows."""

    test_type: str = "WES"
    moi: set | None = None          # OMIM inheritance codes for the gene; None = unknown
    refs: References = field(default_factory=References)
    clinvar: dict | None = None     # override {"category", "stars", "conflict"} (time-split validation)


@dataclass
class Evidence:
    code: str
    status: str                     # met / not_met / not_evaluated
    strength: str | None = None
    reason: str = ""
    points: int = 0

    def as_dict(self) -> dict:
        out = {"code": self.code, "status": self.status, "reason": self.reason}
        if self.status == "met":
            out["strength"] = self.strength
            out["points"] = self.points
        return out


def _met(code: str, strength: str, reason: str, rb: dict) -> Evidence:
    pts = int(rb["points"][strength])
    sign = -1 if code.startswith(("BA", "BS", "BP")) else 1
    return Evidence(code, "met", strength, reason, sign * pts)


def _token(code: str, strength: str) -> str:
    if strength == DEFAULT_STRENGTH.get(code):
        return code
    return f"{code}_{STRENGTH_LABEL[strength]}"


def _class_for(points: int) -> str:
    if points >= 10:
        return "Pathogenic"
    if points >= 6:
        return "Likely pathogenic"
    if points >= 0:
        return "Uncertain significance"
    if points >= -6:
        return "Likely benign"
    return "Benign"


# ── row selection ──────────────────────────────────────────────────────────
def _primary_row(rows: list[dict]) -> dict:
    """MANE Select > MANE Plus Clinical > canonical > other; then most severe impact."""
    def rank(row):
        ttype = P.text(row.get("TRANSCRIPT_TYPE")).upper()
        return (P.TRANSCRIPT_RANK.get(ttype, 9), P.IMPACT_RANK.get(P.text(row.get("IMPACT")).upper(), 9))
    return min(rows, key=rank)


# ── ClinVar ────────────────────────────────────────────────────────────────
def _clinvar(row: dict, ctx: Context) -> dict:
    if ctx.clinvar is not None:
        info = dict(ctx.clinvar)
        info.setdefault("source", "override")
        return info
    if P.text(row.get("CLINVAR_LATEST_SIG")):
        sig, stars, conf, source = (row.get("CLINVAR_LATEST_SIG"), row.get("CLINVAR_LATEST_STARS"),
                                    row.get("CLINVAR_LATEST_SIGCONF"), "ClinVar latest")
    else:
        sig, stars, conf, source = (row.get("CLINVAR_SIG"), row.get("CLINVAR_STARS"),
                                    row.get("CLINVAR_SIGCONF"), "ClinVar")
    category = P.clinvar_category(sig)
    star_value = P.first_number(stars)
    return {
        "category": None if category == "CONFLICT" else category,
        "stars": int(star_value) if star_value is not None else 0,
        "conflict": category == "CONFLICT" or bool(P.text(conf)),
        "source": source,
        "significance": P.text(sig),
    }


def _plp_supported(clinvar: dict, expert: str | None) -> bool:
    if expert in ("P", "LP"):
        return True
    return (clinvar.get("category") in ("P", "LP") and not clinvar.get("conflict")
            and int(clinvar.get("stars") or 0) >= 2)


# ── frequency ──────────────────────────────────────────────────────────────
def _frequency(row: dict, rb: dict) -> dict:
    cfg = rb["frequency"]
    gnomad = {}
    for col in cfg["large_sources"]:
        vals = P.numbers(row.get(col))
        if vals:
            gnomad[col] = max(vals)
    inhouse_af = P.first_number(row.get("INHOUSE_AF"))
    inhouse_ac = P.first_number(row.get("INHOUSE_AC"))
    inhouse_an = P.first_number(row.get("INHOUSE_AN"))
    inhouse_ok = (inhouse_af is not None and inhouse_an is not None and inhouse_ac is not None
                  and inhouse_an >= cfg["inhouse_min_an"] and inhouse_ac >= cfg["inhouse_min_ac"])
    values = list(gnomad.values()) + ([inhouse_af] if inhouse_ok else [])
    return {
        "max_af": max(values) if values else None,
        "gnomad_max": max(gnomad.values()) if gnomad else None,
        "inhouse_af": inhouse_af,
        "inhouse_an": inhouse_an,
        "inhouse_used": inhouse_ok,
        "sources": gnomad,
    }


def _ad_like(ctx: Context, rb: dict) -> bool | None:
    if not ctx.moi:
        return None
    codes = P.moi_codes(ctx.moi)
    return codes.issubset(set(rb["frequency"]["ad_like_moi"]))


def _xlinked_recessive(ctx: Context, rb: dict) -> bool:
    """OMIM lists only X-linked recessive (or unspecified X-linked) inheritance.

    Female carriers are unaffected, so the variant may sit in population data
    at a low frequency; PM2 then uses ``pm2_xlinked_max`` instead of 0.  The
    annotation has no hemizygote (male) counts, so the threshold is on the
    overall AF."""
    codes = P.moi_codes(ctx.moi) if ctx.moi else set()
    return bool(codes) and codes <= set(rb["frequency"].get("pm2_xlinked_moi") or ())


def _recessive_only(ctx: Context) -> bool:
    """OMIM lists only autosomal recessive inheritance for the gene."""
    return bool(ctx.moi) and P.moi_codes(ctx.moi) <= {"AR"}


def _vcep_label(spec: dict) -> str:
    if not spec:
        return ""
    name = (spec.get("vcep") or "VCEP").replace("Variant Curation Expert Panel", "VCEP").strip()
    ident = " ".join(x for x in (spec.get("svi"), f"v{spec['version']}" if spec.get("version") else "") if x)
    return f"ClinGen {name}（{ident}）" if ident else f"ClinGen {name}"


def _vcep_allows(spec: dict, code: str) -> bool:
    """False only when the VCEP marks every strength of ``code`` not applicable."""
    criteria = (spec or {}).get("criteria") or {}
    return code not in criteria or bool(criteria[code])


def _vcep_threshold(spec: dict, name: str, ad_like, recessive: bool):
    thr = (spec or {}).get("thresholds") or {}
    if ad_like and f"{name}_ad" in thr:
        return thr[f"{name}_ad"]
    if recessive and f"{name}_ar" in thr:
        return thr[f"{name}_ar"]
    if name in thr:
        return thr[name]
    if name == "bs1":
        return (spec or {}).get("bs1")          # older spec format
    return None


def _apply_vcep(evidence: list, spec: dict, rb: dict) -> None:
    """Drop criteria the VCEP does not use for this gene; cap strengths above its maximum."""
    criteria = (spec or {}).get("criteria") or {}
    if not criteria:
        return
    label = _vcep_label(spec)
    rank = {"supporting": 1, "moderate": 2, "strong": 3, "very_strong": 4, "stand_alone": 5}
    for e in evidence:
        if e.status != "met" or e.code in ("BA1",) or e.code not in criteria:
            continue
        allowed = criteria[e.code]
        if not allowed:
            e.status, e.points = "not_met", 0
            e.reason += f"（{label} 規格不使用 {e.code}）"
            continue
        top = max(allowed, key=lambda k: rank.get(k, 0))
        if rank.get(e.strength, 0) > rank.get(top, 0):
            sign = -1 if e.code.startswith(("BA", "BS", "BP")) else 1
            e.reason += f"（{label} 規格最高 {STRENGTH_LABEL[top]}，原為 {STRENGTH_LABEL[e.strength]}）"
            e.strength, e.points = top, sign * int(rb["points"][top])


def _pm2_with_pvs1(evidence: list, rb: dict) -> None:
    """PVS1 + PM2_Moderate would reach Pathogenic (10 points) on its own; ClinGen
    counts PM2 as Supporting next to PVS1 so that the pair stays Likely pathogenic."""
    target = (rb.get("pm2") or {}).get("with_pvs1_strength")
    if not target or not any(e.code == "PVS1" and e.status == "met" for e in evidence):
        return
    rank = {"supporting": 1, "moderate": 2, "strong": 3, "very_strong": 4}
    for e in evidence:
        if e.code == "PM2" and e.status == "met" and rank.get(e.strength, 0) > rank.get(target, 0):
            e.reason += f"（與 PVS1 同時成立，以 {STRENGTH_LABEL[target]} 計，避免只靠這兩條就達致病）"
            e.strength, e.points = target, int(rb["points"][target])


def _fmt_af(value) -> str:
    return "缺值" if value is None else f"{value:.5g}"


def _recurrence(row: dict, key: str, ctx: Context, rb: dict, freq: dict, protected: bool) -> tuple[str | None, str | None]:
    """(guard reason, not-evaluated note).  Guard only blocks PVS1/PM2."""
    cfg = rb["recurrence"]
    gnomad = freq["gnomad_max"]
    if ctx.test_type.upper().startswith("WES"):
        result = ctx.refs.wes_recurrence(key)
        if result is None:
            return None, "WES 出現次數表未載入"
        count, patients = result
        if patients < cfg["wes_min_patients"]:
            return None, f"WES 出現次數表只有 {patients} 位病人（滿 {cfg['wes_min_patients']} 位才啟用）"
        frac = count / patients if patients else 0.0
        if frac < cfg["wes_fraction"]:
            return None, None
        inhouse_af = frac / 2.0
        detail = f"院內 WES {count}/{patients} 位病人（{frac:.1%}）"
    else:
        af, an = freq["inhouse_af"], freq["inhouse_an"]
        if af is None or an is None or an < cfg["wgs_min_an"] or af < cfg["wgs_inhouse_af"]:
            return None, None
        inhouse_af = af
        detail = f"院內 AF {af:.3g}（AN {int(an)}）"
    if gnomad is not None and gnomad >= inhouse_af * cfg["gnomad_ratio"]:
        return None, None
    if protected:
        return None, None
    return f"{detail}，gnomAD {_fmt_af(gnomad)}：院內重複出現，疑似假象或東亞常見變異", None


# ── PVS1 ───────────────────────────────────────────────────────────────────
def _pvs1(rows: list[dict], row: dict, ctx: Context, rb: dict, recurrent: str | None) -> tuple[Evidence | None, bool]:
    """Returns (evidence or None if not a LoF variant, splice_fallback)."""
    csq = P.consequences(row.get("CONSEQUENCE"))
    lof = csq & P.LOF_CONSEQUENCES
    if not lof:
        others = [r for r in rows if r is not row and P.consequences(r.get("CONSEQUENCE")) & P.LOF_CONSEQUENCES]
        if others:
            tx = ", ".join(sorted({P.text(r.get("TRANSCRIPT")) for r in others}))
            return Evidence("PVS1", "not_met", reason=f"只在非報告轉錄本（{tx}）造成功能喪失"), False
        return None, False
    if recurrent:
        return Evidence("PVS1", "not_evaluated", reason=recurrent), False
    notes = []
    if P.text(row.get("TRANSCRIPT_TYPE")).upper() not in P.MANE_TYPES:
        notes.append("非 MANE 轉錄本")
    splice_only = bool(lof & P.CANONICAL_SPLICE) and not (lof & (P.NONSENSE | P.FRAMESHIFT | P.START_LOST))
    if splice_only:
        hit = P.canonical_splice_hit(row.get("HGVS_C"))
        if hit is False:
            return Evidence("PVS1", "not_met",
                            reason=f"HGVS 寫法（{P.text(row.get('HGVS_C'))}）不在 ±1、2，剪接位未改變；改用剪接預測"), True
        if hit is None:
            notes.append("HGVS 無法解析，依 consequence 判斷")
    loftee = P.text(row.get("LOFTEE")).upper()
    lfilter = P.text(row.get("LOFTEE_FILTER")).upper()
    end_trunc = False
    if loftee == "LC":
        if "ANC_ALLELE" in lfilter:
            return Evidence("PVS1", "not_met", reason="LOFTEE LC：參考序列為祖先型（ANC_ALLELE）"), False
        if any(f in lfilter for f in ("NON_CAN_SPLICE", "NAGNAG_SITE", "3UTR_SPLICE", "5UTR_SPLICE")):
            return Evidence("PVS1", "not_met", reason=f"LOFTEE LC（{lfilter}）：改用剪接預測"), True
        if any(f in lfilter for f in ("INCOMPLETE_CDS", "EXON_INTRON_UNDEF", "SMALL_INTRON")):
            return Evidence("PVS1", "not_evaluated", reason=f"LOFTEE LC（{lfilter}）：轉錄本結構不完整"), False
        end_trunc = "END_TRUNC" in lfilter
    gene = P.text(row.get("GENE"))
    mechanism = ctx.refs.lof_mechanism(gene)
    if mechanism is None:
        reason = ("功能喪失機轉基因表未載入" if "lof_genes" not in ctx.refs.available
                  else f"{gene}：功能喪失機轉未確立（ClinGen 未評估，ClinVar 致病截斷變異不足）")
        return Evidence("PVS1", "not_evaluated", reason=reason), False
    if mechanism[0] != "yes":
        return Evidence("PVS1", "not_met", reason=f"{gene}：{mechanism[1]}"), False
    basis = [mechanism[1]] + notes
    if lof & P.START_LOST and not (lof & (P.NONSENSE | P.FRAMESHIFT)):
        return _met("PVS1", "moderate", "；".join(["起始密碼子喪失"] + basis), rb), False
    pos, length = P.fraction(row.get("PROTEIN_POSITION"))
    cut = rb["pvs1"]["truncation_fraction"]
    if lof & (P.NONSENSE | P.FRAMESHIFT):
        nmd = _nmd_from_exons(row, ctx, rb)
        if nmd is not None:
            escape, why = nmd
        elif end_trunc:
            escape, why = True, "LOFTEE END_TRUNC（截在尾端）"
        elif "NMD" in row:
            escape = bool(P.text(row.get("NMD")))
            why = "預測逃過 NMD" if escape else "預測會被 NMD 降解"
        else:
            exon, total = P.fraction(row.get("EXON"))
            if not exon or not total:
                return Evidence("PVS1", "not_evaluated", reason="缺 NMD 與 exon 資訊"), False
            escape = exon == total
            why = ("缺 NMD 欄位，依 exon 位置推估：" + ("最後一個 exon，逃過 NMD" if escape else "會被 NMD 降解"))
        if not escape:
            return _met("PVS1", "very_strong", "；".join([why] + basis), rb), False
        return _truncation_strength(gene, pos, length, cut, [why] + basis, ctx, rb), False
    # canonical splice site (passed the HGVS check)
    intron, _total = P.fraction(row.get("INTRON"))
    skipped = None
    if intron:
        skipped = intron if "splice_donor_variant" in lof else intron + 1
    info = ctx.refs.exon_skip(P.text(row.get("TRANSCRIPT")), skipped) if skipped else None
    if info is None:
        return _met("PVS1", "strong", "；".join(["剪接位變異，未確認跳過 exon 後是否 frameshift，保守給 Strong"] + basis), rb), False
    if info.get("frameshift") and not info.get("last_exon"):
        return _met("PVS1", "very_strong", "；".join([f"跳過 exon {skipped} 造成 frameshift，預測 NMD"] + basis), rb), False
    fraction_removed = float(info.get("fraction") or 0.0)
    strength = "strong" if fraction_removed > cut else "moderate"
    return _met("PVS1", strength, "；".join([f"跳過 exon {skipped}（in-frame 或逃過 NMD），移除約 {fraction_removed:.0%} 蛋白質"] + basis), rb), False


def _nmd_from_exons(row: dict, ctx: Context, rb: dict) -> tuple[bool, str] | None:
    """(escapes NMD, reason) from the MANE exon structure (ClinGen SVI 2018 PVS1).

    The premature stop escapes NMD in the last exon or within the 3'-most
    ``nmd_window`` (50) nt of the penultimate exon; otherwise NMD is predicted.
    Counted from the first base of the new stop codon to the last base of the
    penultimate exon.  The VEP NMD plugin is not used first: it flags 52 nt
    (inclusive ±51), uses the variant rather than the new stop of a frameshift
    and treats the first 100 coding bases as escaping, which SVI does not.
    None when the transcript or the stop position is unknown (caller falls back)."""
    lengths = ctx.refs.cds_exons(P.text(row.get("TRANSCRIPT")))
    codon = P.ptc_codon(row.get("HGVS_P"))
    if not lengths or not codon:
        return None
    window = int(rb["pvs1"].get("nmd_window", 50))
    if len(lengths) == 1:
        return True, "單一 exon 轉錄本，預測逃過 NMD"
    ptc = 3 * codon - 2
    junction = sum(lengths[:-1])            # c. position of the last coding base before the last junction
    if ptc > junction:
        if not lengths[-1]:
            return None                     # beyond the CDS: inconsistent with this transcript
        return True, "新的終止密碼子在最後一個 exon，預測逃過 NMD"
    distance = junction - ptc + 1
    if distance > window:
        return False, f"預測會被 NMD 降解（新的終止密碼子距最後一個 exon 交界 {distance} nt，>{window} nt）"
    if not lengths[-1]:
        return None                         # last exon is all UTR: distance also includes UTR not stored
    return True, f"新的終止密碼子在倒數第二個 exon 最後 {distance} nt（≤{window} nt），預測逃過 NMD"


def _truncation_strength(gene, pos, length, cut, basis, ctx, rb) -> Evidence:
    if not pos or not length:
        return Evidence("PVS1", "not_evaluated", reason="逃過 NMD，但缺蛋白質位置，無法判斷截斷比例")
    downstream = ctx.refs.truncating_downstream(gene, pos)
    if downstream:
        return _met("PVS1", "strong", "；".join(basis + ["下游有 ClinVar 致病截斷變異（重要區域）"]), rb)
    removed = (length - pos) / length
    note = [] if downstream is not None else ["截斷變異位置表未載入"]
    if removed > cut:
        return _met("PVS1", "strong", "；".join(basis + [f"移除約 {removed:.0%} 蛋白質"] + note), rb)
    return _met("PVS1", "moderate", "；".join(basis + [f"只移除約 {removed:.0%} 蛋白質"] + note), rb)


# ── computational evidence ────────────────────────────────────────────────
def _splice_score(row: dict) -> tuple[float | None, str]:
    spliceai = P.max_abs(row.get("SPLICEAI_MAX"))
    if spliceai is not None:
        return spliceai, "SpliceAI"
    pangolin = P.max_abs(row.get("PANGOLIN_SCORE"))
    if pangolin is not None:
        return pangolin, "Pangolin"
    return None, ""


def _missense(row: dict, rb: dict, splice: float | None, splice_tool: str) -> Evidence | None:
    if not (P.consequences(row.get("CONSEQUENCE")) & P.MISSENSE):
        return None
    cfg, bp4_max = rb["missense"], rb["splicing"]["bp4_max"]
    pknn = P.text(row.get(cfg["primary"]))
    llr = P.first_number(row.get("PKNN_LLR"))
    code = strength = None
    if pknn:
        m = re.match(r"^(PP3|BP4)(?:[_-](supporting|moderate|strong|very[_-]?strong))?$", pknn, re.I)
        if not m:
            return Evidence("PP3", "not_met", reason=f"P-KNN 不確定（{pknn}），不換工具")
        code, strength = m.group(1).upper(), (m.group(2) or "supporting").lower().replace("-", "_")
        tool = f"P-KNN {pknn}"
    elif llr is not None:
        return Evidence("PP3", "not_met", reason=f"P-KNN 不確定（LLR {llr:g}），不換工具")
    else:
        values = P.numbers(row.get(cfg["fallback"]))
        if not values:
            return Evidence("PP3", "not_evaluated", reason="缺 missense 預測分數（P-KNN、BayesDel）")
        low, high = min(values), max(values)
        for threshold, label in cfg["bayesdel_pp3"]:
            if low >= threshold:
                code, strength = "PP3", label
                break
        if code is None:
            for threshold, label in cfg["bayesdel_bp4"]:
                if high <= threshold:
                    code, strength = "BP4", label
                    break
        if code is None:
            return Evidence("PP3", "not_met", reason=f"BayesDel_noAF 不確定（{low:g}～{high:g}）")
        tool = f"BayesDel_noAF {low:g}" if code == "PP3" else f"BayesDel_noAF {high:g}"
    if code == "BP4":
        if splice is None:
            return Evidence("BP4", "not_met", reason=f"{tool} 偏良性，但缺剪接預測，不給 BP4")
        if splice > bp4_max:
            return Evidence("BP4", "not_met", reason=f"{tool} 偏良性，但 {splice_tool} {splice:g} 預測影響剪接")
    return _met(code, strength, f"missense：{tool}", rb)


def _splicing(row: dict, rb: dict, splice: float | None, tool: str, splice_fallback: bool) -> Evidence | None:
    csq = P.consequences(row.get("CONSEQUENCE"))
    canonical = bool(csq & P.CANONICAL_SPLICE)
    if canonical and not splice_fallback:
        return None
    relevant = csq & (P.INTRONIC | P.SYNONYMOUS | P.MISSENSE | P.CANONICAL_SPLICE)
    if splice is None:
        return Evidence("PP3", "not_evaluated", reason="缺剪接預測（SpliceAI、Pangolin）") if relevant else None
    cfg = rb["splicing"]
    if splice >= cfg["pp3_min"]:
        return _met("PP3", "supporting", f"剪接：{tool} {splice:g} ≥ {cfg['pp3_min']}", rb)
    if splice <= cfg["bp4_max"] and (csq & (P.INTRONIC | P.SYNONYMOUS)) and not (csq & P.MISSENSE):
        return _met("BP4", "supporting", f"剪接：{tool} {splice:g} ≤ {cfg['bp4_max']}", rb)
    return None


def _bp7(row: dict, rb: dict, splice: float | None, tool: str) -> Evidence | None:
    csq = P.consequences(row.get("CONSEQUENCE"))
    cfg = rb["splicing"]
    synonymous = bool(csq & P.SYNONYMOUS)
    deep = bool(csq & P.INTRONIC) and P.deep_intronic(row.get("HGVS_C"), cfg["deep_intron_donor"], cfg["deep_intron_acceptor"])
    if not (synonymous or deep):
        return None
    where = "同義變異" if synonymous else "深層內含子"
    if splice is None:
        return Evidence("BP7", "not_evaluated", reason=f"{where}，缺剪接預測")
    if splice > cfg["bp4_max"]:
        return Evidence("BP7", "not_met", reason=f"{where}，{tool} {splice:g} 預測影響剪接")
    phylop = P.first_number(row.get("PHYLOP100"))
    limit = rb["conservation"]["phylop_not_conserved_below"]
    if phylop is None:
        if rb["splicing"].get("bp7_without_phylop"):
            return _met("BP7", "supporting", f"{where}；{tool} {splice:g}；缺 phyloP，只依剪接預測", rb)
        return Evidence("BP7", "not_evaluated", reason=f"{where}，缺 phyloP")
    if phylop >= limit:
        return Evidence("BP7", "not_met", reason=f"{where}，位置保守（phyloP {phylop:g}）")
    return _met("BP7", "supporting", f"{where}；{tool} {splice:g}；phyloP {phylop:g}", rb)


def _inframe(row: dict, rb: dict) -> Evidence | None:
    csq = P.consequences(row.get("CONSEQUENCE"))
    if csq & P.STOP_LOST:
        return _met("PM4", "moderate", "stop-loss，蛋白質延長", rb)
    if not (csq & P.INFRAME):
        return None
    if "GIAB_STRATA" not in row:
        return _met("PM4", "moderate", "in-frame indel（未確認是否位於重複區）", rb)
    repeat = P.giab_labels(row.get("GIAB_STRATA")) & {"tandem_repeat", "homopolymer"}
    if not repeat:
        return _met("PM4", "moderate", "in-frame indel，不在重複區", rb)
    phylop = P.first_number(row.get("PHYLOP100"))
    limit = rb["conservation"]["phylop_not_conserved_below"]
    if phylop is not None and phylop < limit:
        return _met("BP3", "supporting", f"in-frame indel 位於重複區（{','.join(sorted(repeat))}），位置不保守", rb)
    return Evidence("PM4", "not_met", reason=f"in-frame indel 位於重複區（{','.join(sorted(repeat))}）")


def _same_residue(row: dict, ctx: Context, key: str, rb: dict) -> Evidence | None:
    if not (P.consequences(row.get("CONSEQUENCE")) & P.MISSENSE):
        return None
    parsed = P.hgvs_p(row.get("HGVS_P"))
    if parsed is None:
        return Evidence("PS1", "not_evaluated", reason="無法解析 HGVS_P")
    ref_aa, pos, alt = parsed
    gene = P.text(row.get("GENE"))
    entries = ctx.refs.clinvar_missense(gene, pos)
    if entries is None:
        return Evidence("PS1", "not_evaluated", reason="ClinVar 蛋白質位置表未載入")
    entries = [e for e in entries if e.get("variant_key") != key and not e.get("splicing")]
    ps1_min = int((rb.get("ps1") or {}).get("min_stars", 2))
    pm5_min = int((rb.get("pm5") or {}).get("min_stars", 2))
    same = [e for e in entries if e.get("alt_aa") == alt and int(e.get("stars") or 0) >= ps1_min]
    if same:
        return _met("PS1", "strong", f"與 ClinVar ≥{ps1_min}★ 致病變異 {same[0].get('variant_key')} 造成相同胺基酸改變 p.{ref_aa}{pos}{alt}", rb)
    other = [e for e in entries if e.get("alt_aa") != alt and int(e.get("stars") or 0) >= pm5_min]
    if other:
        # Entries from schema-1 tables have no category: treat as pathogenic (previous behaviour).
        pathogenic = [e for e in other if e.get("category", "P") == "P"]
        best = (pathogenic or other)[0]
        strength = "moderate" if pathogenic else rb.get("pm5", {}).get("lp_strength", "moderate")
        grade = "致病" if pathogenic else "疑似致病"
        stars = int(best.get("stars") or 0)
        return _met("PM5", strength, f"同一位置 p.{ref_aa}{pos} 另有 ClinVar {grade} missense（{best.get('alt_aa')}，{stars}★）", rb)
    return None


def _pm1(row: dict, ctx: Context, rb: dict) -> Evidence | None:
    if not (P.consequences(row.get("CONSEQUENCE")) & P.MISSENSE):
        return None
    gene = P.text(row.get("GENE"))
    spec = ctx.refs.vcep(gene) or {}
    parsed = P.hgvs_p(row.get("HGVS_P"))
    if parsed is None:
        return None
    ref_aa, pos, alt = parsed
    label = _vcep_label(spec)
    if spec and not _vcep_allows(spec, "PM1"):
        return Evidence("PM1", "not_met", reason=f"{label} 規格不使用 PM1")
    if spec.get("pm1_other_gene"):
        return Evidence("PM1", "not_met", reason=f"{label} 的 PM1 熱點只定義給同規格的其他基因")
    rules = spec.get("pm1") or ([{"strength": "moderate", "ranges": spec["pm1_regions"]}]
                                if spec.get("pm1_regions") else [])
    if rules:                                     # VCEP hotspots win; outside them no PM1
        rank = {"supporting": 1, "moderate": 2, "strong": 3, "very_strong": 4}
        hits = [(r["strength"], f"p.{a}" if a == b else f"p.{a}–{b}")
                for r in rules for a, b in r.get("ranges") or () if a <= pos <= b]
        hits += [(r["strength"], f"p.{ref_aa}{pos}{alt}，失去或產生 Cys") for r in rules
                 if r.get("residue_change") == "Cys" and "Cys" in (ref_aa, alt) and alt not in ("=", ref_aa)]
        if hits:
            strength, where = max(hits, key=lambda h: rank.get(h[0], 0))
            return _met("PM1", strength, f"位於 {label} 定義的 PM1 熱點／重要區域（{where}）", rb)
        what = "（規格只列失去或產生 Cys 的變異，其他關鍵位置需人工判斷）" if any(
            r.get("residue_change") for r in rules) else ""
        return Evidence("PM1", "not_met", reason=f"不在 {label} 定義的 PM1 熱點／重要區域{what}")
    collagen = ctx.refs.collagen_gly(gene)
    cfg_col = rb.get("pm1_collagen")
    if collagen and cfg_col and ref_aa == "Gly" and alt not in ("=", "Gly") and pos in collagen["positions"]:
        _, length = P.fraction(row.get("PROTEIN_POSITION"))
        if not length or int(length) == int(collagen.get("length") or 0):
            return _met("PM1", cfg_col["strength"],
                        f"膠原蛋白三股螺旋 Gly-X-Y 的甘胺酸被取代（p.Gly{pos}；UniProt {collagen.get('accession')}）", rb)
    note = f"（{label} 的 PM1 熱點表無法自動讀取，改用 ClinVar 推算）" if spec.get("pm1_unparsed") else ""
    cfg = rb.get("pm1_clinvar")
    if not cfg:
        return None
    window = int(cfg["window"])
    entries = ctx.refs.clinvar_missense_window(gene, pos - window, pos + window)
    if entries is None:
        return None
    residues = {e["protein_pos"] for e in entries
                if int(e.get("stars") or 0) >= cfg["min_stars"] and e["protein_pos"] != pos}
    if len(residues) < cfg["supporting_residues"]:
        return None
    benign = [b for b in ctx.refs.benign_missense_window(gene, pos - window, pos + window) or []
              if int(b.get("stars") or 0) >= cfg["benign_min_stars"]]
    where = f"p.{pos} 前後 {window} 個胺基酸內有 {len(residues)} 個位置有 ClinVar ≥{cfg['min_stars']}★ 致病 missense"
    if benign:
        return Evidence("PM1", "not_met", reason=f"{where}，但也有 {len(benign)} 個良性 missense")
    moderate_at = cfg.get("moderate_residues")    # None (v1.2): ClinVar hotspots stay Supporting
    strength = "moderate" if moderate_at and len(residues) >= moderate_at else "supporting"
    return _met("PM1", strength, f"ClinVar 突變熱點：{where}，沒有良性 missense{note}", rb)


def _pp2(row: dict, ctx: Context, rb: dict) -> Evidence | None:
    if not (P.consequences(row.get("CONSEQUENCE")) & P.MISSENSE):
        return None
    cfg = rb.get("pp2")
    gene = P.text(row.get("GENE"))
    counts = ctx.refs.missense_gene_counts(gene) if cfg else None
    if counts is None:
        return None
    plp, blb = counts
    if plp >= cfg["min_plp"] and blb / (plp + blb) <= cfg["max_benign_fraction"]:
        return _met("PP2", "supporting", f"{gene} 的 ClinVar missense 以致病為主（致病 {plp}、良性 {blb}）", rb)
    return None


# ── main entry ─────────────────────────────────────────────────────────────
# ── PS4 (proband counting) ─────────────────────────────────────────────────
def _ps4(row: dict, key: str, ctx: Context, rb: dict, pm2_met: bool) -> tuple[Evidence | None, int | None]:
    """PS4 from reported affected probands instead of PP5 (ClinGen SVI 2018).

    Dominant disorders only (recessive cases are PM3) and only when PM2 is met:
    a variant seen in affected people means little if it is also common.
    A: other unrelated patients on this platform whose case marked the
    variant Causative.  B: ClinVar submitters reporting it in affected
    individuals, capped at Supporting.  The stronger of the two counts.
    Returns (evidence or None when PS4 does not apply, in-house count used)."""
    cfg = rb.get("ps4")
    if not cfg:
        return None, None
    if not ctx.moi:
        return Evidence("PS4", "not_evaluated", reason="遺傳模式未知"), None
    codes = P.moi_codes(ctx.moi)
    if not codes & set(cfg["dominant_moi"]):
        return Evidence("PS4", "not_met", reason="隱性疾病以 PM3（與致病變異 trans）評估，不用 PS4"), None
    if cfg.get("requires_pm2", True) and not pm2_met:
        return Evidence("PS4", "not_met", reason="PM2 不成立（族群資料庫中不夠罕見）"), None
    inhouse = ctx.refs.inhouse_probands(key)
    # without a Variation ID ClinVar has no record: {} when the case table is loaded, else None
    clinvar = ctx.refs.clinvar_affected(_variation_id(row))
    rank = {"supporting": 1, "moderate": 2, "strong": 3}
    options, parts = [], []
    if inhouse is None:
        parts.append("院內判讀紀錄未載入")
    elif inhouse:
        strength = next((s for n, s in cfg["inhouse_strength"] if inhouse >= n), None)
        parts.append(f"院內另有 {inhouse} 位無親緣病人判為致病（Causative）")
        if strength:
            options.append(strength)
    else:
        parts.append("院內沒有其他病人判為致病")
    if clinvar is None:
        parts.append("ClinVar 病例資料未載入")
    elif clinvar.get("submitters", 0) >= cfg["clinvar_min_submitters"]:
        n, m = clinvar.get("submitters", 0), clinvar.get("individuals") or 0
        extra = f"、至少 {m} 位個案" if m else ""
        parts.append(f"ClinVar {n} 個提交者在受影響個案中報告{extra}（最多 {STRENGTH_LABEL[cfg['clinvar_strength']]}）")
        options.append(cfg["clinvar_strength"])
    else:
        parts.append("ClinVar 沒有受影響個案的提交")
    reason = "；".join(parts)
    if options:
        return _met("PS4", max(options, key=lambda x: rank[x]), reason, rb), inhouse
    if inhouse is None and clinvar is None:
        return Evidence("PS4", "not_evaluated", reason=reason), None
    return Evidence("PS4", "not_met", reason=reason), inhouse


def _variation_id(row: dict) -> str:
    variation = P.text(row.get("CLINVAR_VARIATION_ID")) or P.text(row.get("CLINVAR_LATEST_VARIATION_ID"))
    return variation.split(",")[0].strip() if variation else ""


def _pm3(row: dict, ctx: Context, rb: dict, too_common: bool) -> Evidence | None:
    """PM3 from ClinVar submissions describing homozygous / compound-heterozygous patients.

    Recessive (AR) genes only.  The partner variant and the phase cannot be
    verified from ClinVar text, so the strength is capped (rulebook ``pm3``)."""
    cfg = rb.get("pm3")
    if not cfg:
        return None
    if not ctx.moi:
        return Evidence("PM3", "not_evaluated", reason="遺傳模式未知")
    if not P.moi_codes(ctx.moi) & set(cfg["recessive_moi"]):
        return None                                   # dominant-only genes: PM3 does not apply
    if too_common:
        return Evidence("PM3", "not_met", reason="族群頻率過高（BA1／BS1）")
    found = ctx.refs.clinvar_biallelic(_variation_id(row))
    if found is None:
        return Evidence("PM3", "not_evaluated", reason="ClinVar 病例資料未載入（或為舊版，未含同型合子／複合異型合子）")
    n = int(found.get("submitters") or 0)
    if not n:
        return Evidence("PM3", "not_met", reason="ClinVar 沒有提交者描述同型合子或複合異型合子病人")
    strength = next((st for k, st in cfg["clinvar_strength"] if n >= k), None)
    kinds = "、".join(x for x, on in (("同型合子", found.get("hom")), ("複合異型合子／in trans", found.get("comphet")))
                     if on)
    quotes = "；".join(f"{d.get('scv') or 'SCV'}：「{d.get('snippet')}」" for d in (found.get("detail") or [])[:2])
    reason = f"ClinVar {n} 個提交者描述{kinds}病人（另一個變異的分類與相位無法確認，請核對）"
    if quotes:
        reason += f"；{quotes}"
    if not strength:
        return Evidence("PM3", "not_met", reason=reason)
    return _met("PM3", strength, reason, rb)


def _rank_points(evidence: list, row: dict, rb: dict, engine_points: int) -> float:
    """Ranking points: the P-KNN PP3/BP4 strength is replaced by the P-KNN LLR itself.

    The LLR is calibrated on the same points scale (±1 Supporting, ±2
    Moderate, ±4 Strong), so 3.442 ranks above 2.1 although both are
    PP3_Moderate.  An indeterminate P-KNN call (−1 < LLR < 1) adds its LLR
    when no other PP3/BP4 counted and PVS1 is not met.  Classification is
    unchanged."""
    llr = P.first_number(row.get("PKNN_LLR"))
    if llr is None or not (P.consequences(row.get("CONSEQUENCE")) & P.MISSENSE):
        return float(engine_points)
    pknn = [e for e in evidence if e.code in ("PP3", "BP4") and e.reason.startswith("missense：P-KNN")]
    met = [e for e in pknn if e.status == "met"]
    if met:
        return float(engine_points - met[0].points + llr)
    indeterminate = any(e.code == "PP3" and e.status == "not_met" and e.reason.startswith("P-KNN 不確定")
                        for e in evidence)
    others = any(e.status == "met" and e.code in ("PP3", "BP4", "PVS1") for e in evidence)
    if indeterminate and not others:
        return float(engine_points + llr)
    return float(engine_points)


def evaluate_variant(rows, context: Context | None = None, rulebook: dict | None = None) -> dict:
    rows = [rows] if isinstance(rows, dict) else list(rows)
    if not rows:
        raise ValueError("evaluate_variant needs at least one row")
    ctx = context or Context()
    rb = rulebook or load_rulebook()
    row = _primary_row(rows)
    key = P.variant_key(row)
    clinvar = _clinvar(row, ctx)
    expert = P.clinvar_category(row.get("CLINGEN_VCEP_CLASS"))
    expert = None if expert == "CONFLICT" else expert
    freq = _frequency(row, rb)
    flags: list[dict] = []
    protected = _plp_supported(clinvar, expert)
    # ClinVar P/LP above the BA1 level is trusted only on the ClinGen BA1 exception
    # list or from a VCEP; otherwise BA1 applies and ClinVar is not adopted (rule v1.12)
    clinvar_overridden = False
    override_af = rb["clinvar"].get("override_max_af")
    if (override_af is not None and expert is None and freq["max_af"] is not None
            and freq["max_af"] > override_af and clinvar.get("category") in ("P", "LP")
            and not ctx.refs.ba1_exception(key)):
        protected, clinvar_overridden = False, True
        flags.append({"code": "clinvar_overridden",
                      "message": f"ClinVar {clinvar.get('category')}（{int(clinvar.get('stars') or 0)}★）但最大 AF "
                                 f"{_fmt_af(freq['max_af'])} > {override_af}，且不在 ClinGen BA1 例外清單："
                                 "不採用 ClinVar，依頻率判讀"})
    recurrent, recurrence_note = _recurrence(row, key, ctx, rb, freq, protected)
    if recurrent:
        flags.append({"code": "recurrent", "message": recurrent})
    evidence: list[Evidence] = []

    # frequency: BA1 / BS1 / PM2
    cfg = rb["frequency"]
    max_af = freq["max_af"]
    ad_like = _ad_like(ctx, rb)
    recessive = _recessive_only(ctx)
    exempt = protected or ctx.refs.ba1_exception(key)
    spec = ctx.refs.vcep(P.text(row.get("GENE"))) or {}
    vlabel = _vcep_label(spec)
    if spec:
        flags.append({"code": "vcep", "message": f"套用 {vlabel} 規格（{spec.get('url') or ''}）"})
    ba1 = bs1 = False
    ba1_threshold = _vcep_threshold(spec, "ba1", ad_like, recessive)
    ba1_from = f"（{vlabel}）" if ba1_threshold is not None else ""
    ba1_threshold = cfg["ba1"] if ba1_threshold is None else ba1_threshold
    if max_af is not None and max_af > ba1_threshold and _vcep_allows(spec, "BA1"):
        if exempt:
            evidence.append(Evidence("BA1", "not_met", reason=f"最大 AF {_fmt_af(max_af)} > {ba1_threshold}{ba1_from}，但為已知致病變異或例外清單"))
            flags.append({"code": "frequency_conflict", "message": "族群頻率高但 ClinVar ≥2★ 致病或列於例外清單"})
        else:
            evidence.append(_met("BA1", "stand_alone", f"最大 AF {_fmt_af(max_af)} > {ba1_threshold}{ba1_from}", rb))
            ba1 = True
    bs1_threshold = _vcep_threshold(spec, "bs1", ad_like, recessive)
    bs1_from = f"（{vlabel}）" if bs1_threshold is not None else ""
    if bs1_threshold is None:
        bs1_threshold = cfg["bs1_ad_like"] if ad_like else cfg["bs1_other"]
    if not ba1 and max_af is not None and max_af > bs1_threshold and _vcep_allows(spec, "BS1"):
        if exempt:
            evidence.append(Evidence("BS1", "not_met", reason=f"最大 AF {_fmt_af(max_af)} > {bs1_threshold}，但為已知致病變異或例外清單"))
        else:
            mode = "顯性／X 性聯" if ad_like else ("遺傳模式未知" if ad_like is None else "隱性或混合")
            evidence.append(_met("BS1", "strong", f"最大 AF {_fmt_af(max_af)} > {bs1_threshold}（{mode}）{bs1_from}", rb))
            bs1 = True
    if ba1 or bs1:
        evidence.append(Evidence("PM2", "not_met", reason="族群頻率過高"))
    elif recurrent:
        evidence.append(Evidence("PM2", "not_met", reason="院內重複出現"))
    elif not _vcep_allows(spec, "PM2"):
        evidence.append(Evidence("PM2", "not_met", reason=f"{vlabel} 規格不使用 PM2"))
    else:
        limit = _vcep_threshold(spec, "pm2", ad_like, recessive)
        vcep_limit = limit is not None
        pm2_from = f"（{vlabel}）" if vcep_limit else ""
        xlr = False
        if limit is None:
            xlr = _xlinked_recessive(ctx, rb) and "pm2_xlinked_max" in cfg
            limit = (cfg["pm2_xlinked_max"] if xlr
                     else cfg["pm2_ad_like_max"] if ad_like
                     else cfg.get("pm2_recessive_max", cfg["pm2_other_max"]) if recessive
                     else cfg["pm2_other_max"])
            pm2_from = ("（X 性聯隱性，女性可帶因）" if xlr
                        else "（隱性遺傳）" if recessive and not ad_like else "")
        pm2cfg = rb.get("pm2") or {}
        if not max_af:                                # missing or 0 in every large population database
            evidence.append(_met("PM2", pm2cfg.get("absent_strength", "supporting"),
                                 "未在族群資料庫出現（ACMG 2015 原文；也可能是覆蓋不足）", rb))
        elif max_af <= limit:
            # an expert-panel threshold: its own strength applies (_apply_vcep caps to the spec)
            low = (pm2cfg.get("absent_strength") if vcep_limit
                   else pm2cfg.get("recessive_low_strength") if xlr or (recessive and not ad_like) else None)
            strength = low or pm2cfg.get("rare_strength", "supporting")
            note = (("，X 性聯隱性頻率極低（ACMG 2015 原文）" if xlr
                     else "，隱性遺傳頻率極低（ACMG 2015 原文）") if low and not vcep_limit else "")
            evidence.append(_met("PM2", strength, f"最大 AF {_fmt_af(max_af)} ≤ {limit}{pm2_from}{note}", rb))
        else:
            evidence.append(Evidence("PM2", "not_met", reason=f"最大 AF {_fmt_af(max_af)} > {limit}{pm2_from}"))

    # variant effect
    pvs1, splice_fallback = _pvs1(rows, row, ctx, rb, recurrent)
    if pvs1 is not None:
        evidence.append(pvs1)
    pvs1_met = pvs1 is not None and pvs1.status == "met"
    splice, splice_tool = _splice_score(row)
    for item in (_missense(row, rb, splice, splice_tool),
                 _splicing(row, rb, splice, splice_tool, splice_fallback),
                 _bp7(row, rb, splice, splice_tool),
                 _inframe(row, rb),
                 _same_residue(row, ctx, key, rb),
                 _pm1(row, ctx, rb),
                 _pp2(row, ctx, rb)):
        if item is not None:
            evidence.append(item)

    _apply_vcep(evidence, spec, rb)
    _pm2_with_pvs1(evidence, rb)

    pm2_met = any(e.code == "PM2" and e.status == "met" for e in evidence)
    ps4, ps4_inhouse = _ps4(row, key, ctx, rb, pm2_met)
    if ps4 is not None:
        _apply_vcep([ps4], spec, rb)
        evidence.append(ps4)
    pm3 = _pm3(row, ctx, rb, ba1 or bs1)
    if pm3 is not None:
        _apply_vcep([pm3], spec, rb)
        evidence.append(pm3)
        if ps4 is not None and ps4.status == "met" and pm3.status == "met":
            weaker = min((ps4, pm3), key=lambda e: (e.points, e.code == "PS4"))
            weaker.reason += "（PS4 與 PM3 可能描述同一批病人，只計較強者）"
            weaker.status, weaker.points = "not_met", 0

    # double counting
    met = [e for e in evidence if e.status == "met"]
    if pvs1_met:
        for e in met:
            if e.code in ("PP3", "PM4"):
                e.status, e.reason = "not_met", e.reason + "（PVS1 成立，不重複計分）"
                e.points = 0
    for code in ("PP3", "BP4"):
        same = [e for e in evidence if e.status == "met" and e.code == code]
        if len(same) > 1:
            keep = max(same, key=lambda e: abs(e.points))
            for e in same:
                if e is not keep:
                    e.status, e.reason, e.points = "not_met", e.reason + f"（{code} 只計一次）", 0
    if any(e.status == "met" and e.code == "BP7" for e in evidence):
        for e in evidence:
            if e.status == "met" and e.code == "BP4" and e.reason.startswith("剪接"):
                e.status, e.reason, e.points = "not_met", e.reason + "（已給 BP7）", 0
    pp3 = next((e for e in evidence if e.status == "met" and e.code == "PP3"), None)
    pm1 = next((e for e in evidence if e.status == "met" and e.code == "PM1"), None)
    if pp3 and pm1 and pp3.points + pm1.points > rb["pm1_pp3_cap"]:
        pm1.status, pm1.reason, pm1.points = "not_met", pm1.reason + "（PM1 + PP3 合計上限 4 點）", 0

    met = [e for e in evidence if e.status == "met"]
    engine_points = sum(e.points for e in met)
    engine_class = "Benign" if ba1 else _class_for(engine_points)
    rank = _rank_points(evidence, row, rb, engine_points)
    if ba1:                                  # stand-alone benign ranks as Benign, not as the point sum (v1.12)
        rank = min(rank, float(CLASS_RANGES["B"][2]))

    # ClinVar: adopt (3–4★ / VCEP) or floor (2★ P/LP); never additive
    points, classification = engine_points, engine_class
    action = {"action": "none", "source": clinvar.get("source"),
              "category": clinvar.get("category"), "stars": clinvar.get("stars")}
    adopt = None
    if clinvar_overridden:
        action["action"] = "overridden"
    elif expert in CLASS_RANGES:
        adopt, action["source"], action["category"] = expert, "ClinGen VCEP", expert
    elif (clinvar.get("category") in CLASS_RANGES and not clinvar.get("conflict")
          and int(clinvar.get("stars") or 0) >= rb["clinvar"]["adopt_min_stars"]):
        adopt = clinvar["category"]
    if adopt:
        name, low, high = CLASS_RANGES[adopt]
        if low is not None:
            points = max(points, low)
            rank = max(rank, low)
        if high is not None:
            points = min(points, high)
            rank = min(rank, high)
        classification = name
        action["action"] = "adopt"
    elif (not clinvar_overridden and clinvar.get("category") in ("P", "LP") and not clinvar.get("conflict")
          and int(clinvar.get("stars") or 0) >= rb["clinvar"]["floor_min_stars"]
          and points < rb["clinvar"]["floor_points"]):
        points = rb["clinvar"]["floor_points"]
        rank = max(rank, float(points))
        classification = _class_for(points)
        action["action"] = "floor"

    criteria = {e.code: {"enabled": True, "strength": e.strength} for e in met}
    text = ",".join(_token(code, criteria[code]["strength"]) for code in CRITERIA_ORDER if code in criteria)
    not_evaluated = [{"code": e.code, "reason": e.reason} for e in evidence if e.status == "not_evaluated"]
    if recurrence_note:
        not_evaluated.append({"code": "RECURRENCE", "reason": recurrence_note})
    return {
        "version": rb["version"],
        "variant_key": key,
        "gene": P.text(row.get("GENE")),
        "transcript": P.text(row.get("TRANSCRIPT")),
        "transcript_type": P.text(row.get("TRANSCRIPT_TYPE")),
        "points": int(points),
        "classification": classification,
        "engine_points": int(engine_points),
        "engine_classification": engine_class,
        "criteria": criteria,
        "criteria_text": text,
        "evidence": [e.as_dict() for e in evidence],
        "not_evaluated": not_evaluated,
        "flags": flags,
        "clinvar_rule": action,
        "vcep": ({"svi": spec.get("svi"), "version": spec.get("version"), "vcep": spec.get("vcep"),
                  "url": spec.get("url")} if spec else None),
        "ps4_inhouse": ps4_inhouse,
        "rank_points": round(rank, 3),
        "frequency": {"max_af": max_af, "gnomad_max": freq["gnomad_max"],
                      "inhouse_af": freq["inhouse_af"], "inhouse_used": freq["inhouse_used"]},
    }
