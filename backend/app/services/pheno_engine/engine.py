"""Phenotype engine (plan §4.3–4.8, phase 2).

Reads a phase-1 knowledge-base snapshot (``hpo_kb/{hp}_{hpoa}/``) and turns a
patient's HPO terms into disease and gene rankings, phenotype points P and a
per-term explanation.

* §4.3 preprocessing: alt_id / obsolete → current ID (no replacement → drop),
  keep phenotype-subtree terms only (HP:0000118 itself is uninformative and
  dropped too), remove an ancestor when a more specific patient term is
  present; the kept descendant takes the larger weight.
* §4.4–4.5 S(D) = Σ w·s(q, D), S*(D) = S(D) − μ(D)·Σw, X(D) = S / Σ w·IC*(q);
  order by S*, then X, then disease ID (ties use mid ranks).
* §4.6 G(g) = max S* over the gene's diseases (``disease_gene.tsv`` use=1);
  order by G, then X of the best disease, then gene symbol.  Genes whose
  diseases have no annotations are listed separately (NA).
* §4.8 P(g) from the gene's rank bin, clipped to [−6.1, 13];
  P(g, d) = max(−6.1, P(g) − [P_dis(d*) − P_dis(d)]) with disease ranks taken
  per MONDO group (OMIM and Orphanet entries of one disease share the better
  rank).  The bin values are the provisional development-set numbers.
* §4.7 output files ``pheno_engine_{diseases,genes,explain}.tsv`` and
  ``pheno_engine_run.json``.

The engine is read-only and deterministic: the same snapshot and input give
byte-identical files.  It does not change the current phenotype score.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import tempfile
import threading
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import kb as K

ENGINE_VERSION = "pheno-engine-v2"
LEAD_NATS = 3.0          # §4.8.5 "rank 1 and clearly ahead"
NEAR_TIE_NATS = 1.0      # §4.8.6 near-tie group
# 2026-11-25 (user decision, plan §4.8.3): rank → points is a smooth curve instead of
# 7 rank bins.  points(r) = a0 + c1·ln r + c2·(ln r)² + c3·[r = 1], fitted by maximum
# likelihood (an exponential-family pmf over ranks 1..N) on the development (396/400)
# and second (481) time-split case sets pooled, production engine, 2023-10-09 KB;
# genes over N = 4,888 genes, diseases by MONDO-group rank over 10,308 groups.  A
# shared mid rank (r > 1) gets no rank-1 term.  Both analyses (calibration bits and
# the candidate-ranking simulation) and the locked 152-case test set favoured it
# over the bins; the bins and an 11-bin table stay on the in-house validation list.
# Floors = the old >1000 estimates (genes −6.1, the lowest P; diseases −4.8).
P_MIN, P_MAX = -6.1, 13.0
DISEASE_MIN = -4.8
INF = float("inf")
GENE_CURVE = (9.087, -1.2461, -0.0819, 1.118)
DISEASE_CURVE = (9.921, -1.1761, -0.0767, 1.132)
CALIBRATION = "provisional (plan §4.8.3, phenopacket-store time-split development set)"
EXPLAIN_TOP = 50
ROUND = 9                # S*, G rounded for tie detection only

RELATIONS = ("same", "patient_broader", "patient_specific", "collateral", "none")
RELATION_ZH = {"same": "相同", "patient_broader": "病人較廣", "patient_specific": "病人較具體",
               "collateral": "旁系", "none": "無共同上層"}
REMOVED_ZH = {"unknown": "hp.obo 沒有這個 ID", "obsolete": "已廢除且沒有替代 term",
              "non_phenotype": "非表型 term（遺傳方式、發病年齡等）", "broken_chain": "到不了根節點",
              "root": "表型根節點，沒有資訊量", "redundant_ancestor": "已有更具體的病人 term",
              "duplicate": "重複輸入"}

FILES = ("pheno_engine_diseases.tsv", "pheno_engine_genes.tsv", "pheno_engine_explain.tsv",
         "pheno_engine_run.json")


def _curve(rank: float, curve) -> float:
    r = max(1.0, float(rank))
    a0, c1, c2, c3 = curve
    lr = math.log(r)
    return round(a0 + c1 * lr + c2 * lr * lr + (c3 if r < 1.25 else 0.0), 1)


def gene_points(rank: float) -> float:
    return min(P_MAX, max(P_MIN, _curve(rank, GENE_CURVE)))


def disease_points(rank: float) -> float:
    return min(P_MAX, max(DISEASE_MIN, _curve(rank, DISEASE_CURVE)))


def _mid_ranks(sorted_values: np.ndarray) -> np.ndarray:
    """Mid ranks for values already sorted best → worst (equal values share the mean rank)."""
    n = len(sorted_values)
    out = np.empty(n)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sorted_values[j + 1] == sorted_values[i]:
            j += 1
        out[i:j + 1] = (i + 1 + j + 1) / 2
        i = j + 1
    return out


def _weight(value) -> float:
    try:
        w = float(value)
    except (TypeError, ValueError):
        return 1.0
    return w if math.isfinite(w) and w > 0 else 1.0


def _sha(path: Path) -> str:
    try:
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    except OSError:
        return ""


# ── snapshot ───────────────────────────────────────────────────────────────
class Engine:
    """One loaded snapshot.  ``links`` overrides ``disease_gene.tsv`` (validation only)."""

    def __init__(self, folder: Path, links: list[dict] | None = None):
        self.folder = Path(folder)
        self.kb, self.mu = K.load_kb(self.folder)
        self.kernel = K.Kernel(self.kb)
        self.tid = {t: i for i, t in enumerate(self.kb.terms)}
        self.did = {d: i for i, d in enumerate(self.kb.diseases)}
        self.root = self.tid.get(K.PHENOTYPE_ROOT)
        self._anc = [frozenset(self.kb.anc(i).tolist()) for i in range(self.kb.T)]
        self.status, self.current = self._read_term_status()
        self.manifest = {}
        mpath = self.folder / "manifest.json"
        if mpath.is_file():
            try:
                self.manifest = json.loads(mpath.read_text(encoding="utf-8"))
            except ValueError:
                self.manifest = {}
        self.manifest_sha256 = _sha(mpath)
        self.groups, self.group_labels = self._read_groups()
        if links is None:
            links = self._read_links()
        self._set_links(links)

    # term_status.tsv: every hp.obo ID with its status and current ID
    def _read_term_status(self):
        status, current = {}, {}
        path = self.folder / "term_status.tsv"
        if not path.is_file():
            return status, current
        with open(path, encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                hid = row.get("hpo_id") or ""
                status[hid] = row.get("status") or ""
                current[hid] = row.get("current_id") or ""
        return status, current

    def _read_groups(self):
        groups, labels = {}, {}
        path = self.folder / "disease_groups.tsv"
        if path.is_file():
            with open(path, encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    groups[row["disease_id"]] = row["mondo_id"]
                    labels[row["mondo_id"]] = row.get("mondo_label") or ""
        return groups, labels

    def _read_links(self) -> list[dict]:
        path = self.folder / "disease_gene.tsv"
        if not path.is_file():
            return []
        with open(path, encoding="utf-8", newline="") as handle:
            return [row for row in csv.DictReader(handle, delimiter="\t") if str(row.get("use")) == "1"]

    def _set_links(self, links: list[dict]) -> None:
        by_pair = {}
        for row in links:
            key = (row["disease_id"], row["gene_symbol"])
            prev = by_pair.get(key)
            if prev is None:
                by_pair[key] = {"source": row.get("source", ""), "flags": set(filter(None, str(row.get("flags") or "").split(",")))}
            else:
                prev["source"] = ",".join(sorted((set(prev["source"].split(",")) | {row.get("source", "")}) - {""}))
                prev["flags"] |= set(filter(None, str(row.get("flags") or "").split(",")))
        self.link_info = by_pair
        genes = sorted({g for _, g in by_pair})
        ranked = sorted({g for d, g in by_pair if d in self.did})
        self.genes = ranked
        self.gid = {g: i for i, g in enumerate(ranked)}
        self.na_genes = sorted(set(genes) - set(ranked))
        pairs = sorted((self.gid[g], self.did[d]) for d, g in by_pair if d in self.did)
        self.lg = np.array([p[0] for p in pairs], dtype=np.int64)
        self.ld = np.array([p[1] for p in pairs], dtype=np.int64)
        self.genes_of = {}
        for gi, di in pairs:
            self.genes_of.setdefault(int(di), []).append(self.genes[gi])
        self.diseases_of = {}
        for gi, di in pairs:
            self.diseases_of.setdefault(self.genes[gi], []).append(int(di))
        self.na_diseases = {g: sorted(d for d, gg in by_pair if gg == g) for g in self.na_genes}

    @property
    def snapshot(self) -> str:
        return self.folder.name

    # ── §4.3 preprocessing ────────────────────────────────────────────────
    def normalize(self, hid: str) -> tuple[str | None, str]:
        """(current phenotype ID or None, how / reason)."""
        hid = (hid or "").strip()
        if self.status:
            st = self.status.get(hid)
            if st is None:
                return None, "unknown"
            how = "current"
            if st in ("alt_id", "obsolete"):
                cur = self.current.get(hid) or ""
                if not cur:
                    return None, "obsolete"
                how = "alt_id" if st == "alt_id" else "replaced_by"
                hid, st = cur, self.status.get(cur, "")
            if hid not in self.tid:
                return None, st if st in ("non_phenotype", "broken_chain") else "non_phenotype"
        elif hid not in self.tid:
            return None, "unknown"
        else:
            how = "current"
        if self.tid[hid] == self.root:
            return None, "root"
        return hid, how

    def preprocess(self, items) -> dict:
        """items: [{"phenotype"|"hpo_id": ID, "weight": w}] or [(ID, w)] or [ID]."""
        raw, removed, mapped = [], [], []
        for item in items or []:
            if isinstance(item, dict):
                hid = str(item.get("phenotype") or item.get("hpo_id") or "").strip()
                w = _weight(item.get("weight", 1))
            elif isinstance(item, (tuple, list)):
                hid, w = str(item[0]).strip(), _weight(item[1] if len(item) > 1 else 1)
            else:
                hid, w = str(item).strip(), 1.0
            if not hid.upper().startswith("HP:"):
                continue                                    # panels and other entries are not HPO terms
            raw.append({"hpo_id": hid, "weight": w})
            cur, how = self.normalize(hid)
            if cur is None:
                removed.append({"hpo_id": hid, "weight": w, "reason": how, "current_id": "", "kept": ""})
                continue
            mapped.append((hid, cur, how, w))
        best: dict[int, dict] = {}
        for hid, cur, how, w in mapped:
            i = self.tid[cur]
            if i in best:
                removed.append({"hpo_id": hid, "weight": w, "reason": "duplicate", "current_id": cur, "kept": cur})
                best[i]["weight"] = max(best[i]["weight"], w)
                continue
            best[i] = {"hpo_id": cur, "input_id": hid, "how": how, "weight": w}
        idx = sorted(best)
        keep = []
        for q in idx:
            more_specific = [u for u in idx if u != q and q in self._anc[u]]
            if more_specific:
                removed.append({"hpo_id": best[q]["input_id"], "weight": best[q]["weight"],
                                "reason": "redundant_ancestor", "current_id": self.kb.terms[q],
                                "kept": ";".join(self.kb.terms[u] for u in more_specific)})
            else:
                keep.append(q)
        terms = []
        for q in keep:
            w = max([best[q]["weight"]] + [best[a]["weight"] for a in idx if a != q and a in self._anc[q]])
            terms.append({"idx": q, "hpo_id": self.kb.terms[q], "name": self.kb.names[q], "input_id": best[q]["input_id"],
                          "how": best[q]["how"], "input_weight": best[q]["weight"], "weight": w,
                          "ic_star": float(self.kb.ic_star[q])})
        return {"raw": raw, "terms": terms, "removed": removed}

    # ── §4.4 pair detail (explanation) ────────────────────────────────────
    def pair(self, q: int, d: int) -> dict:
        kb, k = self.kb, self.kernel
        if d == q:
            return {"relation": "same", "mica": q, "ic_mica": float(kb.ic_star[q]), "kappa": 1.0, "lift": None,
                    "e": float(kb.ic_star[q])}
        if q in self._anc[d]:
            return {"relation": "patient_broader", "mica": q, "ic_mica": float(kb.ic_star[q]), "kappa": 1.0,
                    "lift": None, "e": float(kb.ic_star[q])}
        if d in self._anc[q]:
            ic = float(k.icz[d])
            return {"relation": "patient_specific", "mica": d, "ic_mica": ic, "kappa": 1.0, "lift": None, "e": ic}
        common = [a for a in self._anc[q] & self._anc[d] if kb.n[a] > 0]
        m = min(common, key=lambda a: (-float(k.icz[a]), a)) if common else (self.root if self.root is not None else q)
        ic = float(k.icz[m])
        if ic <= 0:
            return {"relation": "none", "mica": m, "ic_mica": 0.0, "kappa": 0.0, "lift": None, "e": 0.0}
        qa = k._count_anchor(q)
        reach_q = k.t_dis[k.t_ptr[qa]:k.t_ptr[qa + 1]]
        reach_d = k.t_dis[k.t_ptr[d]:k.t_ptr[d + 1]]
        n_qd = float(np.intersect1d(reach_q, reach_d, assume_unique=True).size)
        n_m = kb.N * math.exp(-ic)
        lift = (n_qd + 0.5) * (n_m + 0.5) / ((float(kb.n[qa]) + 0.5) * (float(kb.n[d]) + 0.5))
        kappa = K.LAMBDA * min(1.0, lift)
        return {"relation": "collateral", "mica": m, "ic_mica": ic, "kappa": kappa, "lift": lift, "e": kappa * ic}

    def best_annotation(self, q: int, ev: np.ndarray, disease: int) -> int:
        kb = self.kb
        ann = kb.dir_idx[kb.dir_ptr[disease]:kb.dir_ptr[disease + 1]]
        vals = ev[ann]
        top = ann[vals >= vals.max() - 1e-12]
        rel_rank = {r: i for i, r in enumerate(RELATIONS)}

        def key(d):
            d = int(d)
            if d == q:
                r = "same"
            elif q in self._anc[d]:
                r = "patient_broader"
            elif d in self._anc[q]:
                r = "patient_specific"
            else:
                r = "collateral"
            return (rel_rank[r], d)
        return int(min(top, key=key))

    # ── §4.5–4.8 scoring ──────────────────────────────────────────────────
    def score(self, items) -> "Result":
        pre = self.preprocess(items)
        return Result(self, pre)


_CACHE: dict = {}
_CACHE_LOCK = threading.Lock()


def load_engine(folder: Path | None = None) -> Engine | None:
    """Cached engine for ``folder`` (default: the active snapshot); None without a knowledge base."""
    if folder is None:
        from .. import hpo_gene_table as H
        table = H.active_table()
        if table.term_status is None:
            return None
        folder = table.gene_table.parent
    folder = Path(folder)
    if not (folder / "disease_closure.npz").is_file():
        return None
    try:
        st = (folder / "disease_closure.npz").stat()
        sig = (str(folder.resolve()), st.st_mtime_ns, st.st_size, _sha(folder / "manifest.json"))
    except OSError:
        return None
    with _CACHE_LOCK:
        hit = _CACHE.get(str(folder.resolve()))
        if hit and hit[0] == sig:
            return hit[1]
        engine = Engine(folder)
        _CACHE.clear()                      # one snapshot in memory at a time
        _CACHE[str(folder.resolve())] = (sig, engine)
        return engine


# ── result ─────────────────────────────────────────────────────────────────
@dataclass
class Result:
    engine: Engine
    pre: dict
    s: np.ndarray = field(init=False)          # |Q| × N best-match scores
    evecs: list = field(init=False)

    def __post_init__(self):
        eng, kb = self.engine, self.engine.kb
        terms = self.pre["terms"]
        self.q = np.array([t["idx"] for t in terms], dtype=np.int64)
        self.w = np.array([t["weight"] for t in terms], dtype=np.float64)
        self.evecs = [eng.kernel.evec(int(q)) for q in self.q]
        if len(self.q):
            self.s = np.stack([np.maximum.reduceat(ev[kb.dir_idx], eng.kernel.dir_offs) for ev in self.evecs])
        else:
            self.s = np.zeros((0, kb.N))
        self.wsum = float(self.w.sum())
        self.S = self.w @ self.s if len(self.q) else np.zeros(kb.N)
        self.S_star = self.S - eng.mu * self.wsum
        denom = float((self.w * kb.ic_star[self.q]).sum()) if len(self.q) else 0.0
        self.info_total = denom
        self.X = self.S / denom if denom > 0 else np.zeros(kb.N)
        self._rank_diseases()
        self._rank_genes()

    @property
    def empty(self) -> bool:
        return len(self.q) == 0

    def _rank_diseases(self):
        kb, eng = self.engine.kb, self.engine
        key_s = np.round(self.S_star, ROUND)
        key_x = np.round(self.X, ROUND)
        order = np.lexsort((np.arange(kb.N), -key_x, -key_s))
        self.d_order = order
        self.d_pos = np.empty(kb.N, dtype=np.int64)
        self.d_pos[order] = np.arange(kb.N)
        self.d_rank_mid = np.empty(kb.N)
        self.d_rank_mid[order] = _mid_ranks(key_s[order])
        # MONDO groups: OMIM and Orphanet entries of one disease share the better rank
        group_of = [eng.groups.get(d, d) for d in kb.diseases]
        first = {}
        for di in order:
            first.setdefault(group_of[di], int(di))
        g_best = list(first.values())                      # group representatives in rank order
        g_mid = _mid_ranks(key_s[np.array(g_best, dtype=np.int64)]) if g_best else np.zeros(0)
        g_rank = {group_of[di]: float(r) for di, r in zip(g_best, g_mid)}
        self.group_of = group_of
        self.d_group_rank = np.array([g_rank[group_of[i]] for i in range(kb.N)]) if kb.N else np.zeros(0)

    def _rank_genes(self):
        eng = self.engine
        G = len(eng.genes)
        if not G:
            self.g_order = np.zeros(0, dtype=np.int64)
            self.g_best = np.zeros(0, dtype=np.int64)
            self.G = np.zeros(0)
            self.g_rank_mid = np.zeros(0)
            return
        best_pos = np.full(G, np.iinfo(np.int64).max)
        np.minimum.at(best_pos, eng.lg, self.d_pos[eng.ld])
        self.g_best = self.d_order[best_pos]
        self.G = self.S_star[self.g_best]
        key_g = np.round(self.G, ROUND)
        key_x = np.round(self.X[self.g_best], ROUND)
        self.g_order = np.lexsort((np.arange(G), -key_x, -key_g))     # genes are sorted by symbol
        self.g_pos = np.empty(G, dtype=np.int64)
        self.g_pos[self.g_order] = np.arange(G)
        self.g_rank_mid = np.empty(G)
        self.g_rank_mid[self.g_order] = _mid_ranks(key_g[self.g_order])

    # ── per gene / per disease views ───────────────────────────────────────
    def gene_index(self, gene: str) -> int | None:
        return self.engine.gid.get(gene)

    def lead_nats(self) -> float | None:
        if len(self.g_order) < 2 or self.empty:
            return None
        return float(self.G[self.g_order[0]] - self.G[self.g_order[1]])

    def gene_disease_points(self, gene: str) -> list[dict]:
        """P(g, d) for every disease of the gene (plan §4.8.7)."""
        eng = self.engine
        gi = eng.gid.get(gene)
        if gi is None or self.empty:
            return []
        p_gene = gene_points(float(self.g_rank_mid[gi]))
        best = int(self.g_best[gi])
        p_best = disease_points(float(self.d_group_rank[best]))
        out = []
        for di in sorted(eng.diseases_of[gene], key=lambda x: int(self.d_pos[x])):
            p_dis = disease_points(float(self.d_group_rank[di]))
            out.append({"disease_id": eng.kb.diseases[di], "disease_rank": float(self.d_rank_mid[di]),
                        "group_rank": float(self.d_group_rank[di]), "S_star": float(self.S_star[di]),
                        "P": max(P_MIN, p_gene - (p_best - p_dis)), "best": di == best})
        return out

    def gene(self, gene: str) -> dict | None:
        eng = self.engine
        gi = eng.gid.get(gene)
        if gi is None:
            if gene in eng.na_diseases:
                return {"gene": gene, "G": None, "rank": None, "P": None, "status": "na",
                        "diseases": eng.na_diseases[gene]}
            return None
        if self.empty:
            return {"gene": gene, "G": None, "rank": None, "P": None, "status": "no_terms"}
        top = self.g_order[0]
        rank_mid = float(self.g_rank_mid[gi])
        best = int(self.g_best[gi])
        lead = self.lead_nats()
        return {"gene": gene, "G": float(self.G[gi]), "rank": int(self.g_pos[gi]) + 1, "rank_mid": rank_mid,
                "P": gene_points(rank_mid), "best_disease": eng.kb.diseases[best],
                "best_disease_name": eng.kb.disease_names[best], "explained": float(self.X[best]),
                "n_diseases": len(eng.diseases_of[gene]),
                "near_tie": bool(float(self.G[top]) - float(self.G[gi]) <= NEAR_TIE_NATS),
                "lead": bool(gi == top and lead is not None and lead >= LEAD_NATS),
                "status": "ranked"}

    def explain(self, disease: int) -> list[dict]:
        eng, kb = self.engine, self.engine.kb
        mu = float(eng.mu[disease])
        rows = []
        for k, q in enumerate(self.q):
            q = int(q)
            d = eng.best_annotation(q, self.evecs[k], disease)
            p = eng.pair(q, d)
            s = float(self.s[k, disease])
            w = float(self.w[k])
            rows.append({"hpo_id": kb.terms[q], "hpo_name": kb.names[q], "weight": w,
                         "annotation": kb.terms[d], "annotation_name": kb.names[d], "relation": p["relation"],
                         "mica": kb.terms[p["mica"]], "mica_name": kb.names[p["mica"]], "ic_mica": p["ic_mica"],
                         "kappa": p["kappa"], "lift": p["lift"], "s": s, "mu": mu, "net": w * (s - mu)})
        return rows

    # ── summaries ─────────────────────────────────────────────────────────
    def top_diseases(self, n: int = 20) -> list[dict]:
        eng, kb = self.engine, self.engine.kb
        out = []
        if self.empty:
            return out
        for di in self.d_order[:n]:
            di = int(di)
            out.append({"rank": int(self.d_pos[di]) + 1, "rank_mid": float(self.d_rank_mid[di]),
                        "disease_id": kb.diseases[di], "name": kb.disease_names[di],
                        "S_star": float(self.S_star[di]), "S": float(self.S[di]), "mu": float(eng.mu[di]),
                        "explained": float(self.X[di]), "genes": eng.genes_of.get(di, [])})
        return out

    def top_genes(self, n: int = 20) -> list[dict]:
        if self.empty:
            return []
        return [self.gene(self.engine.genes[int(gi)]) for gi in self.g_order[:n]]


# ── output files (§4.7) ────────────────────────────────────────────────────
def _f(x, digits: int = 6) -> str:
    if x is None:
        return ""
    x = float(x)
    if not math.isfinite(x):
        return ""
    return f"{x:.{digits}f}"


def _write_tsv(path: Path, header, rows) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def _annotation_notes(folder: Path, diseases: set) -> dict:
    """{(disease, hpo): (frequency, reference)} from disease_terms.tsv for the given diseases."""
    out = {}
    path = folder / "disease_terms.tsv"
    if not path.is_file() or not diseases:
        return out
    with open(path, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("database_id") in diseases:
                out[(row["database_id"], row["hpo_id"])] = (row.get("frequency_raw") or row.get("frequency") or "",
                                                            row.get("reference") or "")
    return out


def run_info(result: Result, extra: dict | None = None) -> dict:
    eng = result.engine
    return {
        "engine_version": ENGINE_VERSION,
        "snapshot": eng.snapshot,
        "hp_version": eng.manifest.get("hp_version", ""),
        "hpoa_version": eng.manifest.get("hpoa_version", ""),
        "kb_version": (eng.manifest.get("kb") or {}).get("kb_version", ""),
        "manifest_sha256": eng.manifest_sha256,
        "params": {"lambda": K.LAMBDA, "lead_nats": LEAD_NATS, "near_tie_nats": NEAR_TIE_NATS,
                   "p_min": P_MIN, "p_max": P_MAX, "disease_min": DISEASE_MIN,
                   "gene_curve": list(GENE_CURVE), "disease_curve": list(DISEASE_CURVE),
                   "calibration": CALIBRATION},
        "input": result.pre["raw"],
        "terms": [{k: t[k] for k in ("hpo_id", "name", "input_id", "how", "input_weight", "weight", "ic_star")}
                  for t in result.pre["terms"]],
        "removed": [dict(r, reason_zh=REMOVED_ZH.get(r["reason"], "")) for r in result.pre["removed"]],
        "counts": {"diseases": result.engine.kb.N, "genes_ranked": len(eng.genes), "genes_na": len(eng.na_genes),
                   "weight_total": result.wsum, "information_total": result.info_total},
        "lead_nats": result.lead_nats(),
        **(extra or {}),
    }


def write_outputs(result: Result, out_dir: Path, *, prefix: str = "", extra: dict | None = None) -> list[Path]:
    """Write the four §4.7 files into ``out_dir`` (atomic per file)."""
    eng, kb = result.engine, result.engine.kb
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix=".pheno-engine-", dir=out_dir))
    paths = []
    try:
        drows = []
        if not result.empty:
            for di in result.d_order:
                di = int(di)
                grp = eng.groups.get(kb.diseases[di], "")
                drows.append([int(result.d_pos[di]) + 1, _f(result.d_rank_mid[di], 1), kb.diseases[di],
                              kb.disease_names[di], grp, eng.group_labels.get(grp, ""), _f(result.S_star[di]),
                              _f(result.S[di]), _f(eng.mu[di]), _f(result.X[di], 4),
                              int(kb.dir_ptr[di + 1] - kb.dir_ptr[di]), ";".join(eng.genes_of.get(di, []))])
        _write_tsv(tmp / FILES[0], ["rank", "rank_mid", "disease_id", "disease_name", "mondo_id", "mondo_label",
                                    "S_star", "S", "mu", "explained", "n_annotations", "genes"], drows)
        grows = []
        if not result.empty:
            for gi in result.g_order:
                gene = eng.genes[int(gi)]
                g = result.gene(gene)
                best = g["best_disease"]
                info = eng.link_info.get((best, gene), {})
                dpts = ";".join(f"{x['disease_id']}={x['P']:+.1f}" for x in result.gene_disease_points(gene))
                grows.append([g["rank"], _f(g["rank_mid"], 1), gene, _f(g["G"]), _f(g["P"], 1), best,
                              g["best_disease_name"], _f(g["explained"], 4), g["n_diseases"],
                              info.get("source", ""), ",".join(sorted(info.get("flags", ()))),
                              "yes" if g["near_tie"] else "", "yes" if g["lead"] else "", dpts])
        for gene in eng.na_genes:
            grows.append(["NA", "", gene, "", "", "", "", "", len(eng.na_diseases.get(gene, [])), "", "",
                          "", "", ";".join(eng.na_diseases.get(gene, []))])
        _write_tsv(tmp / FILES[1], ["rank", "rank_mid", "gene", "G", "P", "best_disease", "best_disease_name",
                                    "explained", "n_diseases", "source", "flags", "near_tie", "lead",
                                    "disease_points"], grows)
        erows = []
        if not result.empty:
            top = [int(di) for di in result.d_order[:EXPLAIN_TOP]]
            notes = _annotation_notes(eng.folder, {kb.diseases[di] for di in top})
            for di in top:
                for r in result.explain(di):
                    freq, ref = notes.get((kb.diseases[di], r["annotation"]), ("", ""))
                    erows.append([int(result.d_pos[di]) + 1, kb.diseases[di], kb.disease_names[di], r["hpo_id"],
                                  r["hpo_name"], _f(r["weight"], 3), r["annotation"], r["annotation_name"],
                                  r["relation"], RELATION_ZH[r["relation"]], r["mica"], r["mica_name"],
                                  _f(r["ic_mica"], 4), _f(r["kappa"], 4), _f(r["lift"], 4), _f(r["s"], 4),
                                  _f(r["mu"], 4), _f(r["net"], 4), freq, ref])
        _write_tsv(tmp / FILES[2], ["disease_rank", "disease_id", "disease_name", "hpo_id", "hpo_name", "weight",
                                    "best_annotation", "annotation_name", "relation", "relation_zh", "mica",
                                    "mica_name", "ic_mica", "kappa", "lift", "match", "baseline", "net",
                                    "annotation_frequency", "annotation_reference"], erows)
        (tmp / FILES[3]).write_text(json.dumps(run_info(result, extra), ensure_ascii=False, indent=2, sort_keys=True)
                                    + "\n", encoding="utf-8")
        for name in FILES:
            dest = out_dir / f"{prefix}{name}"
            os.replace(tmp / name, dest)
            paths.append(dest)
    finally:
        for leftover in tmp.glob("*"):
            leftover.unlink()
        tmp.rmdir()
    return paths


# ── platform helper ────────────────────────────────────────────────────────
def run_for_analysis(sample_id: str, version: str | None = None, *, engine: Engine | None = None) -> dict:
    """Score a sample's analysis version and write the files next to its analysis.json."""
    from .. import analyses_store, sample_layout
    version = version or analyses_store.active_version(sample_id) or "default"
    data = analyses_store.read_version(sample_id, version)
    if data is None:
        raise FileNotFoundError(f"{sample_id} 沒有分析版本 {version}")
    engine = engine or load_engine()
    if engine is None:
        raise RuntimeError("尚未啟用表型知識庫（請先在首頁按「更新 HPO」）")
    result = engine.score(data.get("hpo") or [])
    target = sample_layout.analysis_file(sample_id, version, FILES[0], for_write=True)
    prefix = target.name[: -len(FILES[0])]
    paths = write_outputs(result, target.parent, prefix=prefix,
                          extra={"sample_id": sample_id, "analysis_version": version})
    return {"sample_id": sample_id, "version": version, "result": result, "paths": paths}
