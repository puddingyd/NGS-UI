"""Offline validation of the ACMG rule engine (plan §5.2.10, implementation step 5).

One pass over every published case reads the 08 review TSV and the stored
engine results (``{LIS_ID}.snv_acmg_engine.sqlite``) and answers:

1. ClinVar-labelled variants (≥2★, no conflict): how well do the engine's
   evidence-only points separate P/LP from B/LB?  The engine is re-run with
   ClinVar and the ClinGen VCEP class masked, so it cannot simply adopt the
   label.  Optional time-split: with ``t0_clinvar`` (an older
   ``variant_summary``) only variants first labelled after T0 are scored;
   build the reference tables from the same old file and pass them as refs.
2. ClinGen VCEP criteria (``CLINGEN_VCEP_CRITERIA``): per-criterion agreement.
3. Reviewer manual ACMG (reusable criteria only): per-criterion agreement.
4. Ranking of reviewer-marked causative SNVs (status 1) with each ACMG source.
5. Coverage: PVS1 on LOFTEE HC loss-of-function variants by inheritance;
   criteria that could not be evaluated.
6. Recurrence guard: flagged variants and WES threshold sensitivity.

The report holds aggregate numbers only.  Per-variant tables carry variant
keys, genes and HGVS but never sample IDs, patient keys or reviewer names.
"""
from __future__ import annotations

import csv
import json
import re
import sqlite3
import statistics
import time
from collections import Counter, defaultdict
from pathlib import Path

from . import parse as P
from . import runner
from .engine import Context, _clinvar, _primary_row, evaluate_variant, load_rulebook

ENGINE_CODES = ("PVS1", "PS1", "PM1", "PM2", "PM4", "PM5", "PP2", "PP3", "BA1", "BS1", "BP3", "BP4", "BP7")
CLASSES = ("P", "LP", "VUS", "LB", "B")
SOURCES = ("engine", "pipeline", "genebe")
RANK_SOURCES = ("engine", "engine_evidence", "pipeline", "genebe", "platform_current")
RECURRENCE_THRESHOLDS = (0.03, 0.05, 0.10)
MASKED_CLINVAR = {"category": None, "stars": 0, "conflict": False, "source": "masked"}
PREDICTOR_COLUMNS = ("SPLICEAI_MAX", "PANGOLIN_SCORE", "PKNN_EVIDENCE", "PKNN_LLR", "BAYESDEL_NOAF", "PHYLOP100")
_DIGITS = re.compile(r"\d+(?:\.\d+)?(?:e-?\d+)?%?")


def af_bin(value) -> str:
    if value is None:
        return "no AF"
    for limit, label in ((0.001, "<0.1%"), (0.005, "0.1–0.5%"), (0.05, "0.5–5%")):
        if value < limit:
            return label
    return ">5%"


def stars_label(info: dict) -> str:
    if info.get("conflict"):
        return "conflicting"
    category = info.get("category")
    if not category:
        return "none"
    return f"{category} {int(info.get('stars') or 0)}★"

_CODE_RE = re.compile(
    r"(?<![A-Z0-9])(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])"
    r"(?:[_ (:]*(very[_ ]?strong|stand[-_ ]?alone|strong|moderate|supporting))?", re.I)
_STRENGTHS = {"verystrong": "very_strong", "standalone": "stand_alone", "strong": "strong",
              "moderate": "moderate", "supporting": "supporting"}


def default_strength(code: str) -> str:
    return {"PVS": "very_strong", "PS": "strong", "PM": "moderate", "PP": "supporting",
            "BA": "stand_alone", "BS": "strong", "BP": "supporting"}[code[:2] if code[:2] != "PV" else "PVS"]


def parse_criteria(text) -> dict[str, str]:
    """Criteria text from any source -> {code: strength}."""
    out: dict[str, str] = {}
    for code, strength in _CODE_RE.findall(P.text(text)):
        code = code.upper()
        norm = re.sub(r"[^a-z]", "", (strength or "").lower())
        out[code] = _STRENGTHS.get(norm) or default_strength(code)
    return out


def class_code(value) -> str | None:
    category = P.clinvar_category(value)
    return category if category in CLASSES else None


def consequence_class(row: dict) -> str:
    cons = P.consequences(row.get("CONSEQUENCE"))
    if cons & (P.NONSENSE | P.FRAMESHIFT | P.START_LOST):
        return "truncating"
    if cons & P.CANONICAL_SPLICE:
        return "canonical_splice"
    if cons & P.MISSENSE:
        return "missense"
    if cons & (P.INFRAME | P.STOP_LOST | {"protein_altering_variant"}):
        return "inframe_other_coding"
    if cons & P.INTRONIC:
        return "splice_region_intronic"
    if cons & P.SYNONYMOUS:
        return "synonymous"
    return "other"


def clinvar_label(row: dict) -> dict:
    info = _clinvar(row, Context())
    category, stars = info.get("category"), int(info.get("stars") or 0)
    group = None
    if not info.get("conflict") and stars >= 2:
        group = "PLP" if category in ("P", "LP") else ("BLB" if category in ("LB", "B") else None)
    return {"group": group, "category": category, "stars": stars, "conflict": bool(info.get("conflict"))}


def first_value(rows: list[dict], column: str) -> str:
    for row in rows:
        value = P.text(row.get(column))
        if value:
            return value
    return ""


def geno_score(points) -> int:
    """Same 0–100 mapping as ``adapters.snv_tsv._acmg_to_geno_score`` (missing -> 0)."""
    if points is None:
        return 0
    x = max(-10.0, min(10.0, float(points)))
    return int(round((x + 10.0) / 20.0 * 100.0))


def auc(positives, negatives) -> float | None:
    """Mann–Whitney AUC (ties count half)."""
    pos, neg = list(positives), list(negatives)
    if not pos or not neg:
        return None
    combined = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    rank_sum, i = 0.0, 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        avg = (i + 1 + j) / 2.0
        rank_sum += avg * sum(1 for k in range(i, j) if combined[k][1] == 1)
        i = j
    return round((rank_sum - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg)), 4)


def _rate(num, den):
    return round(num / den, 4) if den else None


def _moi_group(moi: set | None) -> str:
    if not moi:
        return "unknown"
    dominant = moi & {"AD", "XL", "XLD", "YL"}
    recessive = moi & {"AR", "XLR"}
    if recessive and not dominant:
        return "recessive"
    if dominant and not recessive:
        return "dominant_or_xlinked"
    return "mixed" if dominant and recessive else "other"


def load_t0_labels(path: Path) -> dict[str, str]:
    """``variant_summary`` at T0 -> {variant_key: PLP|BLB} for ≥2★ non-conflicting GRCh38 records."""
    from .ref_store import _open, _variant_key, review_stars

    out: dict[str, str] = {}
    with _open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row = {k.lstrip("#"): v for k, v in row.items() if k}
            if (row.get("Assembly") or "") != "GRCh38":
                continue
            category = P.clinvar_category(row.get("ClinicalSignificance"))
            if category not in CLASSES or review_stars(row.get("ReviewStatus") or "") < 2:
                continue
            key = _variant_key(row)
            if key:
                out[key] = "PLP" if category in ("P", "LP") else ("BLB" if category in ("LB", "B") else "VUS")
    return out


def _wes_recurrence_candidates(db: Path | None) -> tuple[dict[str, float], int]:
    """{variant_key: carrier patient fraction} for fractions ≥ the lowest threshold."""
    if not db or not Path(db).is_file():
        return {}, 0
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        total = conn.execute("SELECT COUNT(DISTINCT patient_key) FROM samples").fetchone()[0] or 0
        if not total:
            return {}, 0
        minimum = max(2, int(min(RECURRENCE_THRESHOLDS) * total))
        rows = conn.execute(
            "SELECT v.variant_key, COUNT(DISTINCT s.patient_key) AS c FROM sample_variants v "
            "JOIN samples s USING (sample_id) GROUP BY v.variant_key HAVING c >= ?", (minimum,))
        return {key: count / total for key, count in rows}, total
    finally:
        conn.close()


def _manual_assertions(db: Path | None) -> dict[str, dict[str, dict]]:
    """{source_sample_id: {variant_id: assertion}} (current revisions, reusable criteria)."""
    out: dict[str, dict[str, dict]] = defaultdict(dict)
    if not db or not Path(db).is_file():
        return out
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        rows = conn.execute(
            "SELECT r.* FROM manual_current c JOIN manual_revisions r ON r.id=c.revision_id "
            "WHERE c.genome_build='hg38'").fetchall()
    except sqlite3.Error:
        rows = []
    finally:
        conn.close()
    for row in rows:
        keys = set(row.keys())
        raw = row["reusable_criteria_json"] if "reusable_criteria_json" in keys else None
        try:
            criteria = json.loads(raw or row["criteria_json"] or "{}")
        except ValueError:
            continue
        criteria = {c.upper(): (v.get("strength") if isinstance(v, dict) else None) or default_strength(c.upper())
                    for c, v in criteria.items() if (v.get("enabled", True) if isinstance(v, dict) else v)}
        cls = (row["reusable_acmg_classification"] if "reusable_acmg_classification" in keys else None) \
            or row["acmg_classification"]
        score = row["reusable_acmg_score"] if "reusable_acmg_score" in keys and row["reusable_acmg_score"] is not None \
            else row["acmg_score"]
        out[row["source_sample_id"]][row["variant_id"]] = {
            "criteria": criteria, "classification": class_code(cls), "score": score}
    return out


def _causative(db: Path | None) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    if not db or not Path(db).is_file():
        return out
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        for sid, vid in conn.execute(
                "SELECT sample_id, variant_id FROM observations WHERE genome_build='hg38' AND status='1'"):
            out[sid].add(vid)
    except sqlite3.Error:
        pass
    finally:
        conn.close()
    return out


def normalize_variant_id(value: str) -> str:
    m = re.match(r"^(?:chr)?([^-]+)-(\d+)-([^-]+)-([^-]+)$", str(value or "").strip(), re.I)
    if not m:
        return ""
    chrom, pos, ref, alt = m.groups()
    chrom = chrom.upper()
    chrom = "M" if chrom == "MT" else chrom
    return f"chr{chrom}-{int(pos)}-{ref.upper()}-{alt.upper()}"


def _default_pheno(sample_id: str):
    """{canonical gene: pheno score} of the active analysis; None when unavailable."""
    try:
        from .. import sample_loader
        ctx = sample_loader._load_pheno_context(sample_id, None)
    except Exception:
        return None
    return (ctx[-1] or {}) if ctx else {}


def _canonical(gene: str) -> str:
    try:
        from .. import panel_deadzone
        return panel_deadzone.canonical_panel_gene_symbol(gene) or gene
    except Exception:
        return gene


class _Confusion:
    def __init__(self):
        self.counts = Counter()

    def add(self, truth, call):
        self.counts[(truth or "none", call or "none")] += 1

    def table(self):
        rows = sorted({t for t, _ in self.counts}, key=lambda c: (CLASSES + ("none",)).index(c)
                      if c in CLASSES + ("none",) else 99)
        cols = sorted({c for _, c in self.counts}, key=lambda c: (CLASSES + ("none",)).index(c)
                      if c in CLASSES + ("none",) else 99)
        return {"rows": rows, "cols": cols, "counts": {f"{t}|{c}": n for (t, c), n in self.counts.items()}}


class _CriterionAgreement:
    def __init__(self):
        self.tp, self.fp, self.fn, self.tn = Counter(), Counter(), Counter(), Counter()
        self.same_strength = Counter()
        self.n = 0
        self.reference_only = Counter()

    def add(self, reference: dict[str, str], engine: dict[str, str]):
        self.n += 1
        for code in ENGINE_CODES:
            ref, eng = code in reference, code in engine
            if ref and eng:
                self.tp[code] += 1
                if reference[code] == engine[code]:
                    self.same_strength[code] += 1
            elif eng:
                self.fp[code] += 1
            elif ref:
                self.fn[code] += 1
            else:
                self.tn[code] += 1
        for code in reference:
            if code not in ENGINE_CODES:
                self.reference_only[code] += 1

    def table(self):
        out = {}
        for code in ENGINE_CODES:
            tp, fp, fn = self.tp[code], self.fp[code], self.fn[code]
            out[code] = {"tp": tp, "fp": fp, "fn": fn, "tn": self.tn[code],
                         "precision": _rate(tp, tp + fp), "recall": _rate(tp, tp + fn),
                         "same_strength": _rate(self.same_strength[code], tp)}
        return {"n": self.n, "criteria": out, "reference_only_criteria": dict(self.reference_only.most_common())}


def _rank(scores: dict[str, float], vid: str) -> float:
    mine = scores[vid]
    higher = sum(1 for v, s in scores.items() if v != vid and s > mine)
    ties = sum(1 for v, s in scores.items() if v != vid and s == mine)
    return 1 + higher + ties / 2.0


def _rank_summary(ranks: list[float]) -> dict:
    if not ranks:
        return {"n": 0}
    return {"n": len(ranks), "median": statistics.median(ranks),
            "top1": _rate(sum(r <= 1 for r in ranks), len(ranks)),
            "top5": _rate(sum(r <= 5 for r in ranks), len(ranks)),
            "top10": _rate(sum(r <= 10 for r in ranks), len(ranks)),
            "top20": _rate(sum(r <= 20 for r in ranks), len(ranks)),
            "mrr": round(sum(1.0 / r for r in ranks) / len(ranks), 4)}


# ── main pass ──────────────────────────────────────────────────────────────
def run(*, out_dir: Path, sample_ids=None, refs=None, recurrence=None, recurrence_db: Path | None = None,
        manual_db: Path | None = None, t0_clinvar: Path | None = None, moi_fn=runner.omim_moi,
        pheno_fn=_default_pheno, progress=None) -> dict:
    from .. import sample_layout

    started = time.time()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    refs = refs if refs is not None else runner.open_references()
    rb = load_rulebook()
    t0 = load_t0_labels(t0_clinvar) if t0_clinvar else None
    wes_candidates, wes_patients = _wes_recurrence_candidates(recurrence_db)
    manual = _manual_assertions(manual_db)
    causative = _causative(manual_db)
    sample_ids = sorted(set(sample_layout.iter_sample_ids())) if sample_ids is None else list(sample_ids)

    moi_cache: dict[str, set | None] = {}

    def moi(gene):
        if gene not in moi_cache:
            moi_cache[gene] = moi_fn(gene) if moi_fn else None
        return moi_cache[gene]

    cases = Counter()
    labelled: dict[str, dict] = {}
    vcep_seen: set[str] = set()
    vcep = _CriterionAgreement()
    vcep_class = _Confusion()
    manual_agree = _CriterionAgreement()
    manual_class, manual_class_final = _Confusion(), _Confusion()
    manual_diff: list[float] = []
    manual_counts = Counter()
    ranks: dict[str, dict[str, dict[str, list]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    rank_counts = Counter()
    pheno_missing = 0
    coverage = defaultdict(lambda: {"occurrences": 0, "classes": Counter(), "not_evaluated": Counter(),
                                    "reasons": Counter(), "flags": Counter(), "lof": defaultdict(Counter),
                                    "predictors": defaultdict(Counter), "columns_absent": Counter(),
                                    "vs_pipeline": _Confusion(), "vs_genebe": _Confusion()})
    vcep_frequency = Counter()
    causal_detail: list[dict] = []
    flagged: dict[str, dict] = {}
    wes_info: dict[str, dict] = {}

    for index, sid in enumerate(sample_ids, 1):
        if progress and (index == 1 or index % 25 == 0 or index == len(sample_ids)):
            progress(f"[{index}/{len(sample_ids)}] {sid}")
        case = runner.case_inputs(sid)
        if not case["review"].is_file() or not Path(case["existing"]).is_file():
            cases["skipped_no_results"] += 1
            continue
        manifest = runner.read_manifest(case["existing"]) or {}
        test_type = manifest.get("test_type") or case["test_type"]
        test_type = "WGS" if str(test_type).upper().startswith(("WGS", "TITAN")) else "WES"
        cases[test_type] += 1
        with open(case["review"], encoding="utf-8", errors="replace", newline="") as handle:
            groups = runner.group_variant_rows(csv.DictReader(handle, delimiter="\t"))
        case_refs = runner.CaseRecurrence(refs, recurrence if test_type == "WES" else None, sid,
                                          runner.case_patient(recurrence, sid, case["mrn"]))
        case_manual = {normalize_variant_id(v): a for v, a in manual.get(sid, {}).items()}
        case_causal = {normalize_variant_id(v) for v in causative.get(sid, ())}
        want_rank = bool(case_causal)
        case_causal_info: dict[str, dict] = {}
        pheno = pheno_fn(sid) if want_rank else {}
        if want_rank and pheno is None:
            pheno_missing += 1
            pheno = {}
        rank_scores = {s: {} for s in RANK_SOURCES}
        acmg_only = {s: {} for s in RANK_SOURCES}
        cov = coverage[test_type]

        def masked(rows, gene):
            rows = [dict(r, CLINGEN_VCEP_CLASS="") for r in rows]
            ctx = Context(test_type=test_type, moi=moi(gene), refs=case_refs, clinvar=MASKED_CLINVAR)
            return evaluate_variant(rows, ctx, rb)

        for vid, result in runner.iter_results(case["existing"], detail=True):
            rows = groups.get(vid)
            if not rows:
                continue
            row = _primary_row(rows)
            key = P.variant_key(row)
            gene = P.text(row.get("GENE"))
            conseq = consequence_class(row)
            pipe_pts = P.first_number(first_value(rows, "ACMG_SCORE"))
            gb_pts = P.first_number(first_value(rows, "GENEBE_ACMG_SCORE"))
            engine_codes = {c: v.get("strength") for c, v in (result.get("criteria") or {}).items()}

            # 5. coverage
            cov["occurrences"] += 1
            cov["classes"][class_code(result.get("classification")) or "none"] += 1
            for item in result.get("not_evaluated") or []:
                cov["not_evaluated"][item.get("code")] += 1
                cov["reasons"][f"{item.get('code')}: {_DIGITS.sub('#', item.get('reason') or '')}"] += 1
            pred = cov["predictors"][conseq]
            pred["n"] += 1
            for column in PREDICTOR_COLUMNS:
                if column not in row:
                    cov["columns_absent"][column] += 1
                elif P.text(row.get(column)):
                    pred[column] += 1
            for flag in result.get("flags") or []:
                cov["flags"][flag.get("code")] += 1
            cov["vs_pipeline"].add(class_code(first_value(rows, "ACMG_CLASS")), class_code(result.get("classification")))
            cov["vs_genebe"].add(class_code(first_value(rows, "GENEBE_ACMG_CLASS")),
                                 class_code(result.get("classification")))
            if (P.consequences(row.get("CONSEQUENCE")) & P.LOF_CONSEQUENCES
                    and P.text(row.get("LOFTEE")).upper() == "HC"):
                bucket = cov["lof"][_moi_group(moi(gene))]
                bucket["n"] += 1
                bucket["engine"] += "PVS1" in engine_codes
                bucket["pipeline"] += "PVS1" in parse_criteria(first_value(rows, "ACMG_CRITERIA"))
                bucket["genebe"] += "PVS1" in parse_criteria(first_value(rows, "GENEBE_ACMG_CRITERIA"))

            # 6. recurrence
            if any(f.get("code") == "recurrent" for f in result.get("flags") or []):
                item = flagged.setdefault(key, {
                    "variant": key, "gene": gene, "hgvs_c": P.text(row.get("HGVS_C")).split(":")[-1],
                    "consequence": conseq, "gnomad_max": (result.get("frequency") or {}).get("gnomad_max"),
                    "inhouse_af": (result.get("frequency") or {}).get("inhouse_af"),
                    "clinvar": P.text(row.get("CLINVAR_LATEST_SIG")) or P.text(row.get("CLINVAR_SIG")),
                    "cases_wes": 0, "cases_wgs": 0})
                item[f"cases_{test_type.lower()}"] += 1
            if test_type == "WES" and key in wes_candidates and key not in wes_info:
                info = clinvar_label(row)
                wes_info[key] = {"gene": gene, "consequence": conseq, "category": info["category"],
                                 "stars": info["stars"], "conflict": info["conflict"],
                                 "gnomad_max": (result.get("frequency") or {}).get("gnomad_max"),
                                 "hgvs_c": P.text(row.get("HGVS_C")).split(":")[-1]}

            # 1. ClinVar-labelled (unique variants)
            label = clinvar_label(row)
            masked_result = None
            if label["group"] and key not in labelled:
                masked_result = masked(rows, gene)
                labelled[key] = {
                    "variant": key, "gene": gene, "consequence": conseq, "label": label["group"],
                    "category": label["category"], "stars": label["stars"], "test_type": test_type,
                    "engine": masked_result["engine_points"], "pipeline": pipe_pts, "genebe": gb_pts,
                    "engine_criteria": masked_result["criteria_text"],
                    "pipeline_criteria": first_value(rows, "ACMG_CRITERIA"),
                    "genebe_criteria": first_value(rows, "GENEBE_ACMG_CRITERIA"),
                    "new_since_t0": (None if t0 is None else t0.get(key) != label["group"]),
                }

            # 2. ClinGen VCEP criteria
            vcep_text = first_value(rows, "CLINGEN_VCEP_CRITERIA")
            if vcep_text and key not in vcep_seen:
                vcep_seen.add(key)
                masked_result = masked_result or masked(rows, gene)
                expert_codes = parse_criteria(vcep_text)
                masked_codes = {c: v.get("strength") for c, v in masked_result["criteria"].items()}
                vcep.add(expert_codes, masked_codes)
                engine_freq = "BA1" if "BA1" in masked_codes else ("BS1" if "BS1" in masked_codes else "neither")
                for code in ("BA1", "BS1"):
                    if code in expert_codes:
                        vcep_frequency[(code, engine_freq, af_bin((masked_result.get("frequency") or {})
                                                                  .get("max_af")))] += 1
                vcep_class.add(class_code(first_value(rows, "CLINGEN_VCEP_CLASS")),
                               class_code(masked_result["engine_classification"]))

            # 3. reviewer manual ACMG
            nid = normalize_variant_id(vid)
            if nid in case_manual:
                assertion = case_manual.pop(nid)
                manual_counts["matched"] += 1
                manual_agree.add(assertion["criteria"], engine_codes)
                manual_class.add(assertion["classification"], class_code(result.get("engine_classification")))
                manual_class_final.add(assertion["classification"], class_code(result.get("classification")))
                if assertion["score"] is not None and result.get("engine_points") is not None:
                    manual_diff.append(float(result["engine_points"]) - float(assertion["score"]))

            # 4. ranking inputs
            if want_rank:
                pheno_score = float(pheno.get(_canonical(gene), 0) or 0)
                points = {"engine": result.get("points"), "engine_evidence": result.get("engine_points"),
                          "pipeline": pipe_pts, "genebe": gb_pts,
                          "platform_current": gb_pts if gb_pts is not None else pipe_pts}
                for source, pts in points.items():
                    rank_scores[source][nid] = geno_score(pts) + pheno_score
                    acmg_only[source][nid] = float(pts) if pts is not None else -99.0
                if nid in case_causal:
                    case_causal_info[nid] = {
                        "test_type": test_type, "clinvar": stars_label(clinvar_label(row)),
                        "consequence": conseq, "points": points,
                        "engine": set(engine_codes),
                        "genebe": set(parse_criteria(first_value(rows, "GENEBE_ACMG_CRITERIA"))),
                        "pipeline": set(parse_criteria(first_value(rows, "ACMG_CRITERIA")))}

        manual_counts["not_in_results"] += len(case_manual)
        if want_rank:
            for target in case_causal:
                if target not in rank_scores["engine"]:
                    rank_counts["causative_not_in_review"] += 1
                    continue
                rank_counts[f"causative_{test_type}"] += 1
                detail = dict(case_causal_info.get(target, {}), ranks={}, ties={})
                for source in RANK_SOURCES:
                    with_pheno = _rank(rank_scores[source], target)
                    ranks[test_type]["with_pheno"][source].append(with_pheno)
                    ranks[test_type]["acmg_only"][source].append(_rank(acmg_only[source], target))
                    detail["ranks"][source] = with_pheno
                    mine = acmg_only[source][target]
                    detail["ties"][source] = sum(1 for v, s in acmg_only[source].items() if v != target and s == mine)
                causal_detail.append(detail)
    # leftover manual assertions whose source case had no results
    manual_counts["source_case_missing"] = sum(
        len(v) for sid, v in manual.items() if sid not in set(sample_ids)) if manual else 0

    report = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "seconds": round(time.time() - started),
        "rulebook": rb.get("version"),
        "reference_tables": sorted(refs.available),
        "reference_versions": getattr(refs, "version", {}),
        "time_split": bool(t0),
        "cases": dict(cases),
        "clinvar": _clinvar_metrics(labelled, t0 is not None),
        "vcep": dict(vcep.table(), classes=vcep_class.table(),
                     frequency_crosstab={f"{a}|{b}|{c}": n for (a, b, c), n in sorted(vcep_frequency.items())}),
        "manual": dict(manual_agree.table(), counts=dict(manual_counts),
                       classes_evidence_only=manual_class.table(), classes_final=manual_class_final.table(),
                       points_difference=_describe(manual_diff)),
        "ranking": {"counts": dict(rank_counts), "pheno_unavailable_cases": pheno_missing,
                    "by_test_type": {tt: {mode: {s: _rank_summary(v) for s, v in by_src.items()}
                                          for mode, by_src in modes.items()} for tt, modes in ranks.items()}},
        "coverage": {tt: _coverage_summary(c) for tt, c in coverage.items()},
        "ranking_breakdown": _ranking_breakdown(causal_detail),
        "criteria_usage": _criteria_usage(labelled),
        "recurrence": _recurrence_summary(flagged, wes_candidates, wes_patients, wes_info, rb),
    }
    _write_tsv(out_dir / "clinvar_labelled.tsv", list(labelled.values()))
    _write_tsv(out_dir / "recurrent_flagged.tsv",
               sorted(flagged.values(), key=lambda r: -(r["cases_wes"] + r["cases_wgs"])))
    _write_tsv(out_dir / "wes_recurrence_candidates.tsv", sorted(
        ({"variant": k, "fraction": round(f, 4), **wes_info.get(k, {})} for k, f in wes_candidates.items()),
        key=lambda r: -r["fraction"]))
    (out_dir / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2, default=str),
                                          encoding="utf-8")
    (out_dir / "report.md").write_text(render_markdown(report), encoding="utf-8")
    return report


def _describe(values: list[float]) -> dict:
    if not values:
        return {"n": 0}
    return {"n": len(values), "mean": round(statistics.mean(values), 3), "median": statistics.median(values),
            "abs_le_1": _rate(sum(abs(v) <= 1 for v in values), len(values))}


def _clinvar_metrics(labelled: dict[str, dict], time_split: bool) -> dict:
    subsets = {"all": list(labelled.values())}
    if time_split:
        subsets["new_since_t0"] = [v for v in labelled.values() if v["new_since_t0"]]
    out = {}
    for name, items in subsets.items():
        by_class: dict[str, list] = defaultdict(list)
        for item in items:
            by_class["all"].append(item)
            by_class[item["consequence"]].append(item)
        out[name] = {cls: _label_block(group) for cls, group in sorted(by_class.items())}
    return out


def _label_block(items: list[dict]) -> dict:
    plp = [i for i in items if i["label"] == "PLP"]
    blb = [i for i in items if i["label"] == "BLB"]
    block = {"plp": len(plp), "blb": len(blb), "sources": {}}
    common = [i for i in items if all(i[s] is not None for s in SOURCES)]
    for source in SOURCES:
        ps = [i[source] for i in plp if i[source] is not None]
        bs = [i[source] for i in blb if i[source] is not None]
        block["sources"][source] = {
            "plp_scored": len(ps), "blb_scored": len(bs),
            "plp_ge6": _rate(sum(v >= 6 for v in ps), len(ps)),
            "plp_ge10": _rate(sum(v >= 10 for v in ps), len(ps)),
            "blb_le_minus1": _rate(sum(v <= -1 for v in bs), len(bs)),
            "blb_le_minus7": _rate(sum(v <= -7 for v in bs), len(bs)),
            "auc": auc(ps, bs),
            "auc_common": auc([i[source] for i in common if i["label"] == "PLP"],
                              [i[source] for i in common if i["label"] == "BLB"]),
        }
    block["common_n"] = len(common)
    return block


def _ranking_breakdown(details: list[dict]) -> dict:
    """Where the engine ranks the causative variant worse than the current platform source."""
    out = {}
    for test_type in sorted({d.get("test_type") for d in details if d.get("test_type")}):
        items = [d for d in details if d.get("test_type") == test_type and d.get("ranks")]
        worse = [d for d in items if d["ranks"]["engine"] > d["ranks"]["platform_current"]]
        better = [d for d in items if d["ranks"]["engine"] < d["ranks"]["platform_current"]]
        only_genebe, only_engine = Counter(), Counter()
        for d in worse:
            only_genebe.update(d.get("genebe", set()) - d.get("engine", set()))
            only_engine.update(d.get("engine", set()) - d.get("genebe", set()))
        out[test_type] = {
            "n": len(items), "engine_better": len(better), "engine_worse": len(worse),
            "same": len(items) - len(better) - len(worse),
            "worse_by_clinvar": dict(Counter(d.get("clinvar", "none") for d in worse).most_common()),
            "worse_by_consequence": dict(Counter(d.get("consequence", "") for d in worse).most_common()),
            "worse_points": _describe([float(d["points"]["engine"] or 0) for d in worse]),
            "worse_platform_points": _describe([float(d["points"]["platform_current"] or 0) for d in worse
                                                if d["points"].get("platform_current") is not None]),
            "worse_genebe_only_criteria": dict(only_genebe.most_common()),
            "worse_engine_only_criteria": dict(only_engine.most_common()),
            "acmg_only_ties": {s: _describe([float(d["ties"][s]) for d in items]) for s in ("engine", "platform_current")},
        }
    return out


def _criteria_usage(labelled: dict[str, dict]) -> dict:
    """How often each source uses each criterion on ClinVar P/LP vs B/LB variants (PP5/BP6 = ClinVar reuse)."""
    out = {}
    for source, column in (("engine", "engine_criteria"), ("pipeline", "pipeline_criteria"),
                           ("genebe", "genebe_criteria")):
        block = {}
        for label in ("PLP", "BLB"):
            items = [v for v in labelled.values() if v["label"] == label]
            counts = Counter()
            for item in items:
                counts.update(parse_criteria(item.get(column)).keys())
            block[label] = {"n": len(items), "criteria": {c: _rate(n, len(items)) for c, n in counts.most_common(12)}}
        out[source] = block
    return out


def _coverage_summary(c: dict) -> dict:
    occ = c["occurrences"]
    return {
        "occurrences": occ,
        "classes": dict(c["classes"]),
        "not_evaluated": {k: {"n": n, "rate": _rate(n, occ)} for k, n in c["not_evaluated"].most_common()},
        "not_evaluated_reasons": dict(c["reasons"].most_common(15)),
        "predictor_presence": {cls: {col: _rate(b[col], b["n"]) for col in PREDICTOR_COLUMNS} | {"n": b["n"]}
                               for cls, b in sorted(c["predictors"].items())},
        "columns_absent": dict(c["columns_absent"]),
        "flags": dict(c["flags"]),
        "lof_hc_pvs1": {group: {"n": b["n"], **{s: _rate(b[s], b["n"]) for s in SOURCES}}
                        for group, b in sorted(c["lof"].items())},
        "class_vs_pipeline": c["vs_pipeline"].table(),
        "class_vs_genebe": c["vs_genebe"].table(),
    }


def _recurrence_summary(flagged, candidates, patients, info, rb) -> dict:
    ratio = rb["recurrence"]["gnomad_ratio"]
    thresholds = {}
    for threshold in RECURRENCE_THRESHOLDS:
        keys = [k for k, f in candidates.items() if f >= threshold]
        blocked = []
        for k in keys:
            item = info.get(k, {})
            gnomad = item.get("gnomad_max")
            if gnomad is not None and gnomad >= candidates[k] / 2.0 * ratio:
                continue
            protected = (item.get("category") in ("P", "LP") and not item.get("conflict")
                         and int(item.get("stars") or 0) >= 2)
            if not protected:
                blocked.append(item)
        thresholds[f"{threshold:.0%}"] = {
            "variants": len(keys), "would_block": len(blocked),
            "would_block_clinvar_plp_0_1_star": sum(1 for i in blocked if i.get("category") in ("P", "LP")
                                                   and not i.get("conflict")),
            "would_block_lof": sum(1 for i in blocked if i.get("consequence") in ("truncating", "canonical_splice")),
        }
    return {"flagged_variants": len(flagged),
            "flagged_occurrences_wes": sum(v["cases_wes"] for v in flagged.values()),
            "flagged_occurrences_wgs": sum(v["cases_wgs"] for v in flagged.values()),
            "wes_patients": patients, "wes_thresholds": thresholds}


def _write_tsv(path: Path, rows: list[dict]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or ["empty"], delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if v is None else v) for k, v in row.items()})


# ── report ─────────────────────────────────────────────────────────────────
_SOURCE_NAMES = {"engine": "規則引擎", "engine_evidence": "規則引擎（不含 ClinVar 採用）",
                 "pipeline": "pipeline ACMG", "genebe": "GeneBe", "platform_current": "平台現行（GeneBe→pipeline）"}


def _pct(value) -> str:
    return "—" if value is None else f"{value * 100:.1f}%"


def _num(value) -> str:
    return "—" if value is None else f"{value}"


def render_markdown(report: dict) -> str:
    lines = [f"# ACMG 規則引擎離線驗證（{report['generated_at'][:10]}）", "",
             f"- 規則版本：{report['rulebook']}；參考表：{', '.join(report['reference_tables'])}",
             f"- 個案：{report['cases']}；耗時 {report['seconds']} 秒",
             f"- ClinVar 時間切分：{'有' if report['time_split'] else '沒有（使用目前的 ClinVar，見限制）'}", ""]

    lines += ["## 1. ClinVar 已有定論的變異（≥2★、無衝突；引擎遮住 ClinVar 重新評估）", ""]
    for subset, blocks in report["clinvar"].items():
        lines += [f"### {('全部' if subset == 'all' else 'T0 之後才有定論')}", "",
                  "| 變異類型 | P/LP | B/LB | 來源 | P/LP ≥6 | P/LP ≥10 | B/LB ≤−1 | AUC | AUC（三者皆有分數） |",
                  "|---|---|---|---|---|---|---|---|---|"]
        for cls, block in blocks.items():
            for source in SOURCES:
                s = block["sources"][source]
                lines.append(f"| {cls} | {block['plp']} | {block['blb']} | {_SOURCE_NAMES[source]} | "
                             f"{_pct(s['plp_ge6'])} | {_pct(s['plp_ge10'])} | {_pct(s['blb_le_minus1'])} | "
                             f"{_num(s['auc'])} | {_num(s['auc_common'])} |")
        lines.append("")

    v = report["vcep"]
    lines += [f"## 2. 對照 ClinGen 專家判讀（{v['n']} 個變異）", "",
              "| 條目 | 雙方都給 | 只有引擎 | 只有專家 | precision | recall | 強度一致 |", "|---|---|---|---|---|---|---|"]
    lines += [f"| {c} | {s['tp']} | {s['fp']} | {s['fn']} | {_pct(s['precision'])} | {_pct(s['recall'])} | "
              f"{_pct(s['same_strength'])} |" for c, s in v["criteria"].items()]
    lines += ["", f"專家用了但引擎不評的條目：{v['reference_only_criteria']}", "",
              "專家給 BA1／BS1 的變異，引擎給了什麼（依最大 AF 分組）：", "",
              "| 專家 | 引擎 | 最大 AF | n |", "|---|---|---|---|"]
    lines += [f"| {k.split('|')[0]} | {k.split('|')[1]} | {k.split('|')[2]} | {n} |"
              for k, n in v.get("frequency_crosstab", {}).items()]
    lines.append("")

    m = report["manual"]
    lines += [f"## 3. 對照判讀者手動 ACMG（{m['n']} 個變異，只比可重複使用的條目）", "",
              f"對應情形：{m['counts']}；引擎點數減判讀者點數：{m['points_difference']}", "",
              "| 條目 | 雙方都給 | 只有引擎 | 只有判讀者 | precision | recall |", "|---|---|---|---|---|---|"]
    lines += [f"| {c} | {s['tp']} | {s['fp']} | {s['fn']} | {_pct(s['precision'])} | {_pct(s['recall'])} |"
              for c, s in m["criteria"].items()]
    lines += ["", f"判讀者用了但引擎不評的條目：{m['reference_only_criteria']}", ""]

    r = report["ranking"]
    lines += ["## 4. 判讀者標為致病（狀態 1）的 SNV 名次", "", f"{r['counts']}；讀不到表型分數的個案：{r['pheno_unavailable_cases']}", ""]
    for tt, modes in r["by_test_type"].items():
        for mode, by_src in modes.items():
            title = "ACMG 點數＋表型分數（平台現行公式）" if mode == "with_pheno" else "只看 ACMG 點數"
            lines += [f"### {tt}：{title}", "", "| 來源 | n | 中位名次 | 第 1 名 | 前 5 | 前 10 | 前 20 | MRR |",
                      "|---|---|---|---|---|---|---|---|"]
            for source, s in by_src.items():
                if s.get("n"):
                    lines.append(f"| {_SOURCE_NAMES[source]} | {s['n']} | {s['median']} | {_pct(s['top1'])} | "
                                 f"{_pct(s['top5'])} | {_pct(s['top10'])} | {_pct(s['top20'])} | {s['mrr']} |")
            lines.append("")

    lines += ["### 引擎名次比平台現行差的致病變異", ""]
    for tt, b in report.get("ranking_breakdown", {}).items():
        lines += [f"- {tt}：共 {b['n']}，引擎較好 {b['engine_better']}、相同 {b['same']}、較差 {b['engine_worse']}",
                  f"  - 較差者的 ClinVar：{b['worse_by_clinvar']}",
                  f"  - 較差者的變異類型：{b['worse_by_consequence']}",
                  f"  - 較差者的點數（引擎）：{b['worse_points']}；（平台現行）：{b['worse_platform_points']}",
                  f"  - 較差者 GeneBe 有、引擎沒有的條目：{b['worse_genebe_only_criteria']}",
                  f"  - 較差者引擎有、GeneBe 沒有的條目：{b['worse_engine_only_criteria']}",
                  f"  - 只看 ACMG 點數時，與致病變異同分的變異數：{b['acmg_only_ties']}"]
    lines += ["", "### 各來源在 ClinVar P/LP 與 B/LB 變異上用了哪些條目（比例）", ""]
    for source, block in report.get("criteria_usage", {}).items():
        lines.append(f"- {_SOURCE_NAMES[source]}：P/LP {block['PLP']['criteria']}；B/LB {block['BLB']['criteria']}")
    lines.append("")

    lines += ["## 5. 覆蓋率", ""]
    for tt, c in report["coverage"].items():
        lines += [f"### {tt}（{c['occurrences']} 個變異次）", "", f"分類：{c['classes']}", "",
                  "LOFTEE HC 功能喪失變異拿到 PVS1 的比例：", "",
                  "| 遺傳模式 | n | 規則引擎 | pipeline | GeneBe |", "|---|---|---|---|---|"]
        lines += [f"| {g} | {b['n']} | {_pct(b['engine'])} | {_pct(b['pipeline'])} | {_pct(b['genebe'])} |"
                  for g, b in c["lof_hc_pvs1"].items()]
        lines += ["", f"無法評估：{ {k: v['n'] for k, v in c['not_evaluated'].items()} }", f"旗標：{c['flags']}", "",
                  "無法評估的原因（前 15）：", ""]
        lines += [f"- {reason}：{n}" for reason, n in c["not_evaluated_reasons"].items()]
        lines += ["", f"TSV 沒有這些欄位的變異次：{c['columns_absent']}", "",
                  "| 變異類型 | n | " + " | ".join(PREDICTOR_COLUMNS) + " |",
                  "|---|---|" + "---|" * len(PREDICTOR_COLUMNS)]
        lines += [f"| {cls} | {b['n']} | " + " | ".join(_pct(b[col]) for col in PREDICTOR_COLUMNS) + " |"
                  for cls, b in c["predictor_presence"].items()]
        lines.append("")

    rec = report["recurrence"]
    lines += ["## 6. 院內重複出現的檢查", "",
              f"被擋下的變異 {rec['flagged_variants']} 個（WES {rec['flagged_occurrences_wes']} 次、"
              f"WGS {rec['flagged_occurrences_wgs']} 次）；WES 出現次數表 {rec['wes_patients']} 位病人。", "",
              "| WES 門檻 | 達門檻變異 | 會被擋（gnomAD 不夠高且非 ClinVar ≥2★ 致病） | 其中 ClinVar 0–1★ P/LP | 其中功能喪失 |",
              "|---|---|---|---|---|"]
    lines += [f"| {t} | {s['variants']} | {s['would_block']} | {s['would_block_clinvar_plp_0_1_star']} | "
              f"{s['would_block_lof']} |" for t, s in rec["wes_thresholds"].items()]
    lines += ["", "## 限制", "",
              "- 沒有時間切分時，PS1／PM5 與功能喪失基因表仍用目前的 ClinVar（變異本身已排除），結果偏樂觀。",
              "- pipeline 與 GeneBe 的分數可能已參考 ClinVar，和引擎遮住 ClinVar 後的分數相比對它們有利。",
              "- 名次中致病變異可能由本實驗室提交到 ClinVar，規則引擎會直接採用，所以另列不含 ClinVar 採用的版本。",
              "- 只涵蓋 review TSV 內的變異；不在 review TSV 的判讀者標記不列入。", ""]
    return "\n".join(lines)
