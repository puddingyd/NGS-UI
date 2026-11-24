"""Literature validation of the phenotype engine (plan §11, §15 stage 6).

Read-only evaluation on phenopacket-store cases:

* **case sets** — locked ID lists (``research/pheno_benchmark/*.txt``); new sets
  are drawn once with :func:`select` (≤3 cases per gene cohort, fixed seed);
* **conditions** (plan §11.2) — *time-split*: a frozen old knowledge base and
  only cases whose publication the old annotations do not cite (the main
  condition); *LOPO*: the current knowledge base with the case publication's
  annotations removed from the target disease (a pessimistic stress test);
  *full*: the current knowledge base as is (an optimistic upper bound);
* **methods** — the engine (plan ch. 4), the engine without baseline
  correction, classic Resnik best-match, and the platform's previous
  presence/absence table (``hpo_gene.tsv``); external tools (LIRICAL,
  Exomiser) are read from per-case result files when given;
* **statistics** (plan §11.5–11.6) — top-1/5/10, MRR and median with mid
  ranks for ties, exact McNemar on discordant cases, paired bootstrap of the
  MRR difference resampling gene cohorts, coverage reported separately;
* **calibration** — rank-bin points log₂.₀₈(f₁/f₀) for genes and diseases
  (plan §4.8.2), with a cohort bootstrap interval.

Nothing here writes to the platform; :func:`write_report` writes only to the
folder it is given.  Case IDs are public phenopacket IDs.
"""
from __future__ import annotations

import csv
import glob
import json
import math
import os
import random
import re
import statistics
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import engine as E

ROUND = 9
GENE_BINS = ((1, 1), (2, 5), (6, 10), (11, 50), (51, 200), (201, 1000), (1001, None))
BASE = 2.08
METHOD_LABELS = {
    "engine": "表型引擎（本計畫）",
    "nobaseline": "表型引擎，不做基準分校正",
    "resnik": "Resnik best-match（經典）",
    "platform": "平台舊表型分數（有或無查表）",
    "lirical": "LIRICAL（表型）",
    "exomiser": "Exomiser（表型）",
}
_PMID = re.compile(r"PMID[:_](\d+)")


# ── cases ──────────────────────────────────────────────────────────────────
@dataclass
class Case:
    id: str
    obs: list
    dis: list
    genes: list
    pmids: frozenset = frozenset()
    cohort: str = ""


def load_cases(folder: Path) -> list[Case]:
    """phenopacket-store JSON files below ``folder`` (cohort = parent folder, usually the gene)."""
    out = []
    for path in sorted(glob.glob(os.path.join(str(folder), "**", "*.json"), recursive=True)):
        try:
            data = json.load(open(path, encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if not isinstance(data, dict) or "phenotypicFeatures" not in data:
            continue
        obs = [f.get("type", {}).get("id", "") for f in data.get("phenotypicFeatures") or [] if not f.get("excluded")]
        dis = [d.get("term", {}).get("id", "") for d in data.get("diseases") or []]
        genes = set()
        for it in data.get("interpretations") or []:
            diag = it.get("diagnosis") or {}
            for gi in diag.get("genomicInterpretations") or []:
                vd = (gi.get("variantInterpretation") or {}).get("variationDescriptor") or {}
                for sym in ((vd.get("geneContext") or {}).get("symbol"), (gi.get("gene") or {}).get("symbol")):
                    if sym:
                        genes.add(sym)
            if (diag.get("disease") or {}).get("id"):
                dis.append(diag["disease"]["id"])
        refs = json.dumps((data.get("metaData") or {}).get("externalReferences") or []) + str(data.get("id") or "")
        out.append(Case(id=str(data.get("id")), obs=[o for o in obs if o], dis=sorted(set(filter(None, dis))),
                        genes=sorted(genes), pmids=frozenset(_PMID.findall(refs)),
                        cohort=os.path.basename(os.path.dirname(path))))
    return out


def read_ids(path: Path) -> list[str]:
    text = Path(path).read_text(encoding="utf-8")
    try:
        data = json.loads(text)
    except ValueError:
        return [line.split("\t")[0].strip() for line in text.splitlines() if line.strip() and not line.startswith("#")]
    return [r["case"] if isinstance(r, dict) else str(r) for r in data]


def write_ids(path: Path, ids, header: str = "") -> None:
    lines = [f"# {line}" for line in header.splitlines() if line] + list(ids)
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def hpoa_citations(hpoa: Path) -> dict[str, set]:
    """disease → PMIDs cited by its phenotype (aspect P) annotations."""
    cited: dict[str, set] = defaultdict(set)
    header = None
    with open(hpoa, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            row = dict(zip(header, parts))
            if row.get("aspect", "P") != "P":
                continue
            cited[row.get("database_id", "")].update(_PMID.findall(row.get("reference", "")))
    return cited


def target_disease(case: Case, engine) -> str | None:
    return next((d for d in case.dis if d in engine.did), None)


def timesplit_eligible(cases, engine, cited, *, min_terms: int = 3) -> list[Case]:
    """Plan §11.2: target disease in the old knowledge base, a gene, ≥min_terms terms,
    and no case PMID among the old annotations of that disease."""
    out = []
    for c in cases:
        tgt = target_disease(c, engine)
        if tgt is None or len(c.obs) < min_terms or not c.genes:
            continue
        if c.pmids & cited.get(tgt, set()):
            continue
        out.append(c)
    return out


def select(cases, *, per_cohort: int = 3, n: int | None = None, seed: int = 0,
           exclude_ids=(), exclude_cohorts=()) -> list[Case]:
    """At most ``per_cohort`` cases per cohort, then ``n`` at random; deterministic for a seed."""
    rng = random.Random(seed)
    skip_ids, skip_coh = set(exclude_ids), set(exclude_cohorts)
    by = defaultdict(list)
    for c in sorted(cases, key=lambda x: x.id):
        if c.id in skip_ids or c.cohort in skip_coh:
            continue
        by[c.cohort].append(c)
    pool = []
    for coh in sorted(by):
        items = by[coh]
        pool.extend(rng.sample(items, min(per_cohort, len(items))))
    if n is not None and n < len(pool):
        pool = rng.sample(pool, n)
    return sorted(pool, key=lambda x: x.id)


# ── ranking helpers ────────────────────────────────────────────────────────
def mid_rank(values: np.ndarray, i: int) -> float:
    v = np.round(values, ROUND)
    x = v[i]
    return float((v > x).sum() + ((v == x).sum() + 1) / 2)


def gene_max(engine, dscores: np.ndarray) -> np.ndarray:
    G = np.full(len(engine.genes), -np.inf)
    np.maximum.at(G, engine.lg, dscores[engine.ld])
    return G


# ── LOPO (plan §11.2) ──────────────────────────────────────────────────────
@dataclass
class LopoSpec:
    case_id: str
    disease: int
    kept: np.ndarray                 # kept direct annotation term indices
    removed: int
    total: int
    mu: float = field(default=float("nan"))


def disease_annotation_refs(snapshot: Path) -> dict[tuple[str, str], set]:
    """(disease, HPO term) → PMIDs, from the snapshot's ``disease_terms.tsv``."""
    refs: dict[tuple[str, str], set] = defaultdict(set)
    path = Path(snapshot) / "disease_terms.tsv"
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        cols = reader.fieldnames or []
        ref_col = next((c for c in ("references", "reference", "refs") if c in cols), None)
        for row in reader:
            text = row.get(ref_col, "") if ref_col else ""
            did = row.get("database_id") or row.get("disease_id") or ""
            refs[(did, row.get("hpo_id", ""))].update(_PMID.findall(text))
    return refs


def lopo_specs(engine, cases, refs) -> dict[str, LopoSpec]:
    kb = engine.kb
    out = {}
    for c in cases:
        tgt = target_disease(c, engine)
        if tgt is None or not c.pmids:
            continue
        di = engine.did[tgt]
        terms = kb.dir_idx[kb.dir_ptr[di]:kb.dir_ptr[di + 1]]
        kept = np.array([t for t in terms if not (refs.get((tgt, kb.terms[int(t)]), set()) & c.pmids)], dtype=np.int64)
        out[c.id] = LopoSpec(c.id, di, kept, int(len(terms) - len(kept)), int(len(terms)))
    return out


def lopo_baselines(engine, specs: dict, progress=None) -> None:
    """μ for each case's reduced target disease, in one pass over the baseline pool."""
    kb, kern = engine.kb, engine.kernel
    todo = [s for s in specs.values() if s.removed and len(s.kept)]
    if not todo:
        for s in specs.values():
            s.mu = 0.0 if not len(s.kept) else float(engine.mu[s.disease])
        return
    pool = np.unique(kb.dir_idx)
    total = np.zeros(len(todo))
    for k, q in enumerate(pool):
        ev = kern.evec(int(q))
        for j, s in enumerate(todo):
            total[j] += ev[s.kept].max()
        if progress and (k % 2000 == 0 or k == len(pool) - 1):
            progress(k + 1, len(pool))
    for j, s in enumerate(todo):
        s.mu = float(total[j] / len(pool))
    for s in specs.values():
        if not s.removed:
            s.mu = float(engine.mu[s.disease])
        elif not len(s.kept):
            s.mu = 0.0


def apply_lopo(res, spec: LopoSpec) -> None:
    """Replace the target disease's scores in a Result and re-rank (other diseases unchanged)."""
    di = spec.disease
    if res.empty:
        return
    if len(spec.kept):
        col = np.array([ev[spec.kept].max() for ev in res.evecs])
    else:
        col = np.zeros(len(res.q))
    res.s[:, di] = col
    S = float(res.w @ col)
    res.S = res.S.copy(); res.S_star = res.S_star.copy(); res.X = res.X.copy()
    res.S[di] = S
    res.S_star[di] = S - spec.mu * res.wsum
    res.X[di] = S / res.info_total if res.info_total > 0 else 0.0
    res._rank_diseases()
    res._rank_genes()


# ── comparators ────────────────────────────────────────────────────────────
class Comparators:
    """Alternative phenotype scorers on the same knowledge base."""

    def __init__(self, engine, hpo_gene: Path | None = None):
        self.engine = engine
        self._mica_cache: dict[int, np.ndarray] = {}
        self.term_genes: dict[str, set] = {}
        path = Path(hpo_gene) if hpo_gene else Path(engine.folder) / "hpo_gene.tsv"
        if path.is_file():
            with open(path, encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    self.term_genes.setdefault(row.get("hpo_id", ""), set()).add(row.get("gene_symbol", ""))

    def mica(self, q: int) -> np.ndarray:
        hit = self._mica_cache.get(q)
        if hit is not None:
            return hit
        kb, kern = self.engine.kb, self.engine.kernel
        anc = kb.anc(q)
        m = np.zeros(kb.T)
        for a in anc[np.argsort(kern.icz[anc], kind="stable")]:
            if kb.n[a] > 0:
                d = kern.desc[a]
                m[d] = np.maximum(m[d], kern.icz[a])
        if len(self._mica_cache) > 4096:
            self._mica_cache.clear()
        self._mica_cache[q] = m
        return m

    def resnik(self, res) -> np.ndarray:
        kb, kern = self.engine.kb, self.engine.kernel
        S = np.zeros(kb.N)
        for q, w in zip(res.q, res.w):
            S += w * np.maximum.reduceat(self.mica(int(q))[kb.dir_idx], kern.dir_offs)
        return S

    def platform(self, res) -> np.ndarray:
        """Gene scores of the previous platform table: matched weight over all terms (§4.2.3, §10.7)."""
        eng = self.engine
        G = np.zeros(len(eng.genes))
        gid = eng.gid
        for q, w in zip(res.q, res.w):
            for g in self.term_genes.get(eng.kb.terms[int(q)], ()):
                i = gid.get(g)
                if i is not None:
                    G[i] += w
        return G


def external_ranks(folder: Path | None) -> dict[str, dict]:
    """Per-case ranks from an external tool: ``{case}.tsv`` with columns gene_symbol / disease_id in rank order."""
    out: dict[str, dict] = {}
    if folder is None or not Path(folder).is_dir():
        return out
    for path in Path(folder).glob("*.tsv"):
        genes, diseases = [], []
        with open(path, encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                g = row.get("gene_symbol") or row.get("geneSymbol") or ""
                d = row.get("disease_id") or row.get("diseaseCurie") or ""
                if g and g not in genes:
                    genes.append(g)
                if d and d not in diseases:
                    diseases.append(d)
        out[path.stem] = {"genes": genes, "diseases": diseases}
    return out


def _list_rank(items: list, target, universe: int) -> float:
    if target in items:
        return float(items.index(target) + 1)
    rest = max(universe - len(items), 1)
    return len(items) + (rest + 1) / 2


# ── evaluation ─────────────────────────────────────────────────────────────
def evaluate(engine, cases, *, methods=("engine", "nobaseline", "resnik", "platform"), lopo=None,
             comparators: Comparators | None = None, external: dict | None = None, progress=None) -> tuple[list, dict]:
    """Per-case ranks of the known disease and gene under each method.

    Returns (rows, coverage); a case whose disease is not in the knowledge base
    counts in coverage and gets no row (plan §11.2)."""
    comparators = comparators or Comparators(engine)
    external = external or {}
    rows, coverage = [], {"cases": 0, "disease_missing": 0, "gene_missing": 0, "no_terms": 0}
    for k, c in enumerate(cases):
        coverage["cases"] += 1
        tgt = target_disease(c, engine)
        if tgt is None:
            coverage["disease_missing"] += 1
            continue
        res = engine.score(c.obs)
        if res.empty:
            coverage["no_terms"] += 1
            continue
        spec = (lopo or {}).get(c.id)
        if spec is not None:
            apply_lopo(res, spec)
        di = engine.did[tgt]
        gene = next((g for g in c.genes if g in engine.gid), None)
        if gene is None:
            coverage["gene_missing"] += 1
        gi = engine.gid.get(gene) if gene else None
        row = {"case": c.id, "cohort": c.cohort, "disease": tgt, "gene": gene or "", "terms": len(res.q),
               "lopo_removed": spec.removed if spec else None, "lopo_total": spec.total if spec else None}
        if "engine" in methods:
            row["engine|d"] = float(res.d_rank_mid[di])
            row["engine|dg"] = float(res.d_group_rank[di])
            row["engine|g"] = float(res.g_rank_mid[gi]) if gi is not None else None
        if "nobaseline" in methods:
            row["nobaseline|d"] = mid_rank(res.S, di)
            row["nobaseline|g"] = mid_rank(gene_max(engine, res.S), gi) if gi is not None else None
        if "resnik" in methods:
            S = comparators.resnik(res)
            if spec is not None:
                S[di] = float(sum(w * (comparators.mica(int(q))[spec.kept].max() if len(spec.kept) else 0.0)
                                  for q, w in zip(res.q, res.w)))
            row["resnik|d"] = mid_rank(S, di)
            row["resnik|g"] = mid_rank(gene_max(engine, S), gi) if gi is not None else None
        if "platform" in methods and gi is not None:
            row["platform|g"] = mid_rank(comparators.platform(res), gi)
        for tool, data in external.items():
            hit = data.get(c.id)
            if hit is None:
                continue
            row[f"{tool}|g"] = _list_rank(hit["genes"], gene, len(engine.genes)) if gene else None
            row[f"{tool}|d"] = _list_rank(hit["diseases"], tgt, engine.kb.N)
        rows.append(row)
        if progress and (k % 50 == 0 or k == len(cases) - 1):
            progress(k + 1, len(cases))
    return rows, coverage


# ── statistics ─────────────────────────────────────────────────────────────
def metrics(ranks) -> dict:
    r = [x for x in ranks if x is not None]
    n = len(r)
    if not n:
        return {"n": 0}
    return {"n": n, "top1": sum(x <= 1 for x in r) / n, "top5": sum(x <= 5 for x in r) / n,
            "top10": sum(x <= 10 for x in r) / n, "mrr": sum(1 / x for x in r) / n,
            "median": float(statistics.median(r))}


def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact McNemar p (binomial on discordant pairs)."""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n
    return min(1.0, 2 * tail)


def paired(rows, a: str, b: str, *, k: int = 10, boot: int = 2000, seed: int = 1) -> dict:
    """Method a vs b on cases where both have ranks: ΔMRR with a cohort bootstrap, McNemar on top-k."""
    pairs = [(r["cohort"], r[a], r[b]) for r in rows if r.get(a) is not None and r.get(b) is not None]
    if not pairs:
        return {"n": 0}
    d = [1 / x - 1 / y for _, x, y in pairs]
    by = defaultdict(list)
    for (coh, _, _), v in zip(pairs, d):
        by[coh].append(v)
    groups = list(by.values())
    rng = random.Random(seed)
    stats = []
    for _ in range(boot):
        pick = [groups[rng.randrange(len(groups))] for _ in range(len(groups))]
        flat = [v for g in pick for v in g]
        stats.append(sum(flat) / len(flat))
    stats.sort()
    better = sum(x <= k < y for _, x, y in pairs)
    worse = sum(y <= k < x for _, x, y in pairs)
    return {"n": len(pairs), "delta_mrr": sum(d) / len(d), "ci": [stats[int(0.025 * boot)], stats[int(0.975 * boot) - 1]],
            "k": k, "better": better, "worse": worse, "p_mcnemar": mcnemar_exact(better, worse)}


def calibrate(ranks, universe: int, *, bins=GENE_BINS, cohorts=None, boot: int = 2000, seed: int = 1) -> list[dict]:
    """Rank-bin points log₂.₀₈(f₁/f₀) with f₀ = bin width / universe (plan §4.8.2)."""
    r = np.array([x for x in ranks if x is not None], dtype=float)
    coh = [c for x, c in zip(ranks, cohorts or [None] * len(ranks)) if x is not None]

    def points(arr):
        out = []
        for lo, hi in bins:
            hi_v = universe if hi is None else hi
            f1 = float(((arr >= lo - 0.5) & (arr < hi_v + 0.5)).sum() / max(len(arr), 1))
            f0 = (hi_v - lo + 1) / universe
            out.append((f1, f0, math.log(f1 / f0) / math.log(BASE) if f1 > 0 else None))
        return out

    base = points(r)
    rng = random.Random(seed)
    groups = defaultdict(list)
    for x, c in zip(r, coh):
        groups[c].append(x)
    keys = list(groups)
    draws = [[] for _ in bins]
    for _ in range(boot):
        arr = np.array([x for _ in keys for x in groups[keys[rng.randrange(len(keys))]]])
        for i, (_, _, p) in enumerate(points(arr)):
            draws[i].append(p if p is not None else -99.0)
    rows = []
    for (lo, hi), (f1, f0, p), d in zip(bins, base, draws):
        d.sort()
        rows.append({"bin": f"{lo}" if hi == lo else (f">{lo - 1}" if hi is None else f"{lo}–{hi}"),
                     "k": int(round(f1 * len(r))), "f1": f1, "f0": f0, "lr": (f1 / f0) if f0 else None,
                     "points": None if p is None else round(p, 1),
                     "ci": [round(d[int(0.025 * boot)], 1), round(d[int(0.975 * boot) - 1], 1)]})
    return rows


def subgroups(rows, key: str) -> list[dict]:
    out = []
    for label, test in (("HPO ≤3", lambda t: t <= 3), ("HPO 4–7", lambda t: 4 <= t <= 7), ("HPO ≥8", lambda t: t >= 8)):
        sel = [r[key] for r in rows if r.get(key) is not None and test(r["terms"])]
        out.append({"group": label, **metrics(sel)})
    return out


# ── report ─────────────────────────────────────────────────────────────────
def _fmt(x, digits=3):
    return "—" if x is None else f"{x:.{digits}f}"


def method_table(rows, methods, level: str) -> list[str]:
    lines = ["| 方法 | n | top-1 | top-5 | top-10 | MRR | 名次中位數 |", "|---|---:|---:|---:|---:|---:|---:|"]
    for m in methods:
        mt = metrics([r.get(f"{m}|{level}") for r in rows])
        if not mt.get("n"):
            continue
        lines.append(f"| {METHOD_LABELS.get(m, m)} | {mt['n']} | {_fmt(mt['top1'])} | {_fmt(mt['top5'])} | "
                     f"{_fmt(mt['top10'])} | {_fmt(mt['mrr'])} | {mt['median']:.1f} |")
    return lines


def write_report(out: Path, sections: list[dict], meta: dict) -> None:
    """sections: [{title, note, rows, coverage, methods}], meta: versions and calibration."""
    out.mkdir(parents=True, exist_ok=True)
    md = [f"# 表型引擎文獻驗證（階段 6）", "", meta.get("intro", ""), ""]
    summary = {"meta": {k: v for k, v in meta.items() if k != "intro"}, "sections": []}
    for sec in sections:
        rows, methods = sec["rows"], sec["methods"]
        cov = sec["coverage"]
        md += [f"## {sec['title']}", "", sec.get("note", ""), "",
               f"病例 {cov['cases']}：可評估 {len(rows)}；疾病不在知識庫 {cov['disease_missing']}"
               f"（覆蓋率 {(cov['cases'] - cov['disease_missing']) / max(cov['cases'], 1):.1%}）；"
               f"基因不在疾病–基因對應 {cov['gene_missing']}；沒有可用 HPO {cov['no_terms']}。", ""]
        md += ["**基因層級**", ""] + method_table(rows, methods, "g") + [""]
        md += ["**疾病層級**", ""] + method_table(rows, [m for m in methods if m != "platform"], "d") + [""]
        comps = []
        for m in methods:
            if m == "engine":
                continue
            for level in ("g", "d"):
                p = paired(rows, f"engine|{level}", f"{m}|{level}")
                if p.get("n"):
                    comps.append((m, level, p))
        if comps:
            md += ["**成對比較（本計畫 − 對照；ΔMRR 以基因群組 bootstrap 2,000 次的 95% CI；top-10 改善／變差與 exact McNemar）**", "",
                   "| 對照 | 層級 | n | ΔMRR [95% CI] | top-10 改善 vs 變差 | p |", "|---|---|---:|---|---|---:|"]
            for m, level, p in comps:
                md.append(f"| {METHOD_LABELS.get(m, m)} | {'基因' if level == 'g' else '疾病'} | {p['n']} | "
                          f"{p['delta_mrr']:+.3f} [{p['ci'][0]:+.3f}, {p['ci'][1]:+.3f}] | {p['better']} vs {p['worse']} | "
                          f"{p['p_mcnemar']:.2g} |")
            md.append("")
        sg = subgroups(rows, "engine|g")
        md += ["**次群組（表型引擎，基因層級）**", "", "| HPO 數 | n | top-1 | top-10 | MRR |", "|---|---:|---:|---:|---:|"]
        for s in sg:
            if s.get("n"):
                md.append(f"| {s['group']} | {s['n']} | {_fmt(s['top1'])} | {_fmt(s['top10'])} | {_fmt(s['mrr'])} |")
        md.append("")
        summary["sections"].append({"title": sec["title"], "coverage": cov,
                                    "metrics": {f"{m}|{lv}": metrics([r.get(f'{m}|{lv}') for r in rows])
                                                for m in methods for lv in ("g", "d")},
                                    "paired": [{"method": m, "level": lv, **p} for m, lv, p in comps]})
        with open(out / f"rows_{sec['key']}.tsv", "w", encoding="utf-8", newline="") as handle:
            keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("case", "cohort", "disease", "gene", "terms"), k))
            w = csv.DictWriter(handle, fieldnames=keys, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    for title, cal in meta.get("calibration", {}).items():
        md += [f"## 名次點數校準：{title}", "", cal.get("note", ""), "",
               "| 名次 | 答案數 | $f_1$ | $f_0$ | 似然比 | 點數 | 95% 區間 |", "|---|---:|---:|---:|---:|---:|---|"]
        for r in cal["rows"]:
            pts = "—" if r["points"] is None else f"{r['points']:+.1f}"
            md.append(f"| {r['bin']} | {r['k']} | {r['f1']:.1%} | {r['f0']:.3%} | {r['lr']:.3g} | {pts} | "
                      f"{r['ci'][0]:+.1f} 至 {r['ci'][1]:+.1f} |")
        md.append("")
    (out / "report.md").write_text("\n".join(md), encoding="utf-8")
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
