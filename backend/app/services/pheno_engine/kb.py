"""Phenotype knowledge base: IC, the match kernel and baseline scores (plan §2.3, §4.2.5, §4.4, §4.5.4).

Built from the phase-0 cleaning result (``hpo_gene_table.build``):

* universe 𝒟 = OMIM + Orphanet diseases with at least one cleaned phenotype
  annotation (DECIPHER excluded), N = |𝒟|;
* p(t) = share of diseases whose ancestor-closed annotations contain t,
  IC(t) = −ln p(t); IC*(t) = IC of the most informative annotated ancestor
  for terms no disease reaches;
* e(q, d) = κ·IC(MICA): IC(q) when d is q or a descendant, IC(d) when d is an
  ancestor of q, else 0.5·min(1, lift)·IC(MICA) with
  lift = (n_qd+.5)(n_m+.5)/((n_q+.5)(n_d+.5));
* s(q, D) = max over D's direct annotations; μ(D) = mean of s(t, D) over the
  pool 𝒫 of all directly annotated phenotype terms.

Everything is a pure function of the inputs, so a rebuild of the same
versions gives identical files.
"""
from __future__ import annotations

import csv
import hashlib
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import npzio

KB_VERSION = "hpo-kb-v1"
LAMBDA = 0.5
PHENOTYPE_ROOT = "HP:0000118"


@dataclass
class KB:
    terms: list            # phenotype-subtree term ids, sorted
    names: list
    anc_ptr: np.ndarray    # CSR: ancestors (incl. self) inside the phenotype subtree
    anc_idx: np.ndarray
    diseases: list         # universe, sorted
    disease_names: list
    dir_ptr: np.ndarray    # CSR: direct annotations per disease
    dir_idx: np.ndarray
    cl_ptr: np.ndarray     # CSR: ancestor-closed annotations per disease
    cl_idx: np.ndarray
    n: np.ndarray          # diseases reaching each term
    ic: np.ndarray         # nan where n == 0
    ic_star: np.ndarray

    @property
    def T(self) -> int:
        return len(self.terms)

    @property
    def N(self) -> int:
        return len(self.diseases)

    def anc(self, i: int) -> np.ndarray:
        return self.anc_idx[self.anc_ptr[i]:self.anc_ptr[i + 1]]


# ── build ──────────────────────────────────────────────────────────────────
def _csr(lists: list) -> tuple[np.ndarray, np.ndarray]:
    ptr = np.zeros(len(lists) + 1, dtype=np.int64)
    ptr[1:] = np.cumsum([len(x) for x in lists])
    idx = np.fromiter((v for x in lists for v in x), dtype=np.int32, count=int(ptr[-1]))
    return ptr, idx


def build_kb(result: dict) -> KB:
    onto = result["ontology"]
    terms = sorted(t for t in onto.names if t not in onto.obsolete and onto.status(t) == "phenotype")
    tid = {t: i for i, t in enumerate(terms)}
    anc_lists = [sorted(tid[a] for a in onto.ancestors(t) if a in tid) for t in terms]
    direct: dict[str, set] = {}
    for (disease, term) in result["merged"]:
        if disease.startswith("DECIPHER:") or term not in tid:
            continue
        direct.setdefault(disease, set()).add(tid[term])
    diseases = sorted(direct)
    dir_lists = [sorted(direct[d]) for d in diseases]
    cl_lists = []
    for ts in dir_lists:
        cl = set()
        for t in ts:
            cl.update(anc_lists[t])
        cl_lists.append(sorted(cl))
    anc_ptr, anc_idx = _csr(anc_lists)
    dir_ptr, dir_idx = _csr(dir_lists)
    cl_ptr, cl_idx = _csr(cl_lists)
    n = np.bincount(cl_idx, minlength=len(terms)).astype(np.int64)
    with np.errstate(divide="ignore"):
        ic = np.where(n > 0, -np.log(np.maximum(n, 1) / max(len(diseases), 1)), np.nan) + 0.0  # no −0.0
    ic_star = ic.copy()
    for i in np.flatnonzero(n == 0):
        anc = np.asarray(anc_lists[i], dtype=np.int64)
        vals = ic[anc]
        ic_star[i] = float(np.nanmax(vals)) if np.isfinite(vals).any() else 0.0
    names = [onto.names.get(t, "") for t in terms]
    dnames = [result["disease_names"].get(d, "") for d in diseases]
    return KB(terms, names, anc_ptr, anc_idx, diseases, dnames, dir_ptr, dir_idx, cl_ptr, cl_idx, n, ic, ic_star)


# ── kernel ─────────────────────────────────────────────────────────────────
class Kernel:
    """Vectorised e(q, ·) over all terms and s(q, ·) over all diseases."""

    def __init__(self, kb: KB):
        self.kb = kb
        T = kb.T
        desc: list[list[int]] = [[] for _ in range(T)]
        for i in range(T):
            for a in kb.anc(i):
                desc[int(a)].append(i)
        self.desc = [np.asarray(d, dtype=np.int64) for d in desc]
        self.dis_of_entry = np.repeat(np.arange(kb.N), np.diff(kb.cl_ptr))
        order = np.argsort(kb.cl_idx, kind="stable")
        self.t_ptr = np.zeros(T + 1, dtype=np.int64)
        self.t_ptr[1:] = np.cumsum(np.bincount(kb.cl_idx, minlength=T))
        self.t_dis = self.dis_of_entry[order]              # diseases reaching each term
        self.icz = np.nan_to_num(kb.ic)
        self.nf = kb.n.astype(np.float64)
        self.dir_offs = kb.dir_ptr[:-1]

    def _count_anchor(self, q: int) -> int:
        """q itself, or its most informative annotated ancestor when no disease reaches q."""
        kb = self.kb
        if kb.n[q] > 0:
            return q
        anc = kb.anc(q)
        ok = anc[kb.n[anc] > 0]
        best = ok[np.lexsort((ok, -self.icz[ok]))]
        return int(best[0])

    def evec(self, q: int) -> np.ndarray:
        kb = self.kb
        T = kb.T
        anc = kb.anc(q)
        mica = np.zeros(T)
        for a in anc[np.argsort(self.icz[anc], kind="stable")]:
            if kb.n[a] > 0:
                d = self.desc[a]
                mica[d] = np.maximum(mica[d], self.icz[a])
        qa = self._count_anchor(q)
        mask = np.zeros(kb.N, dtype=bool)
        mask[self.t_dis[self.t_ptr[qa]:self.t_ptr[qa + 1]]] = True
        n_qd = np.bincount(kb.cl_idx[mask[self.dis_of_entry]], minlength=T).astype(np.float64)
        n_m = kb.N * np.exp(-mica)
        lift = (n_qd + 0.5) * (n_m + 0.5) / ((self.nf[qa] + 0.5) * (self.nf + 0.5))
        e = LAMBDA * np.minimum(1.0, lift) * mica
        e[anc] = self.icz[anc]
        e[self.desc[q]] = kb.ic_star[q]
        return e

    def best(self, q: int) -> np.ndarray:
        """s(q, D) for every disease D."""
        return np.maximum.reduceat(self.evec(q)[self.kb.dir_idx], self.dir_offs)


def baseline(kb: KB, kernel: Kernel | None = None, progress=None) -> tuple[np.ndarray, int]:
    kernel = kernel or Kernel(kb)
    pool = np.unique(kb.dir_idx)
    total = np.zeros(kb.N)
    for k, q in enumerate(pool):
        total += kernel.best(int(q))
        if progress and (k % 1000 == 0 or k == len(pool) - 1):
            progress(k + 1, len(pool))
    return total / max(len(pool), 1), int(len(pool))


# ── QC ─────────────────────────────────────────────────────────────────────
def qc(kb: KB, onto, mu: np.ndarray) -> dict:
    tid = {t: i for i, t in enumerate(kb.terms)}
    root = tid.get(PHENOTYPE_ROOT)
    unreachable = [kb.terms[i] for i in range(kb.T) if root is None or root not in set(kb.anc(i).tolist())]
    violations = []
    for t, i in tid.items():
        if kb.n[i] == 0:
            continue
        for p in onto.parents.get(t, ()):
            j = tid.get(p)
            if j is not None and kb.ic[j] > kb.ic[i] + 1e-12:
                violations.append(f"{p}>{t}")
    annotated = kb.n > 0
    checks = {
        "all_terms_reach_root": not unreachable,
        "ic_ancestor_le_descendant": not violations,
        "root_ic_zero": root is not None and kb.n[root] == kb.N and abs(float(kb.ic[root])) < 1e-12,
        "ic_finite": bool(np.isfinite(kb.ic[annotated]).all() and np.isfinite(kb.ic_star).all()),
        "baseline_finite": bool(np.isfinite(mu).all()),
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "N": kb.N,
        "terms": kb.T,
        "terms_with_ic": int(annotated.sum()),
        "terms_p0": int((~annotated).sum()),
        "unreachable": unreachable[:20],
        "ic_violations": violations[:20],
    }


# ── write / read ───────────────────────────────────────────────────────────
def _sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _fmt(x) -> str:
    return "" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.6f}"


def write_files(kb: KB, mu: np.ndarray, out: Path) -> list[str]:
    out = Path(out)
    direct = set(np.unique(kb.dir_idx).tolist())
    with open(out / "terms.tsv", "w", encoding="utf-8", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["idx", "hpo_id", "name", "n_diseases", "ic", "ic_star", "direct"])
        for i, t in enumerate(kb.terms):
            w.writerow([i, t, kb.names[i], int(kb.n[i]), _fmt(float(kb.ic[i])), _fmt(float(kb.ic_star[i])),
                        int(i in direct)])
    npzio.save(out / "ancestors.npz", terms=np.array(kb.terms), anc_ptr=kb.anc_ptr, anc_idx=kb.anc_idx,
               n=kb.n, ic=kb.ic, ic_star=kb.ic_star)
    npzio.save(out / "disease_closure.npz", diseases=np.array(kb.diseases), dir_ptr=kb.dir_ptr,
               dir_idx=kb.dir_idx, cl_ptr=kb.cl_ptr, cl_idx=kb.cl_idx, mu=mu)
    with open(out / "disease_baseline.tsv", "w", encoding="utf-8", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["disease_id", "disease_name", "n_direct", "n_closure", "mu"])
        for i, d in enumerate(kb.diseases):
            w.writerow([d, kb.disease_names[i], int(kb.dir_ptr[i + 1] - kb.dir_ptr[i]),
                        int(kb.cl_ptr[i + 1] - kb.cl_ptr[i]), _fmt(float(mu[i]))])
    return ["terms.tsv", "ancestors.npz", "disease_closure.npz", "disease_baseline.tsv"]


def content_hashes(folder: Path, names) -> dict:
    return {name: _sha(Path(folder) / name) for name in sorted(names)}


def load_kb(folder: Path) -> tuple[KB, np.ndarray]:
    """Read a built snapshot back (for the engine and tests)."""
    folder = Path(folder)
    a = npzio.load(folder / "ancestors.npz")
    c = npzio.load(folder / "disease_closure.npz")
    names, dnames = [], []
    with open(folder / "terms.tsv", encoding="utf-8", newline="") as handle:
        names = [r["name"] for r in csv.DictReader(handle, delimiter="\t")]
    with open(folder / "disease_baseline.tsv", encoding="utf-8", newline="") as handle:
        dnames = [r["disease_name"] for r in csv.DictReader(handle, delimiter="\t")]
    kb = KB([str(t) for t in a["terms"]], names, a["anc_ptr"], a["anc_idx"], [str(d) for d in c["diseases"]],
            dnames, c["dir_ptr"], c["dir_idx"], c["cl_ptr"], c["cl_idx"], a["n"], a["ic"], a["ic_star"])
    return kb, c["mu"]
