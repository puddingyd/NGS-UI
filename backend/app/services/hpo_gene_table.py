"""Cleaned HPO → gene lookup built from ``phenotype.hpoa`` (plan §4.2.3, §10.7).

The official ``phenotype_to_genes.txt`` propagates *every* annotation of a
disease to the term's ancestors: inheritance / onset / clinical-course terms
(aspect I/C/M/H), explicit exclusions (qualifier ``NOT``), the pre-2023
``Excluded`` frequency and ``0/m`` rows all add genes.  This module rebuilds the
table from the source files with the plan's rules:

* R1 keep ``aspect = P`` only;
* R2 ``qualifier = NOT`` → negative table;
* R3 frequency ``HP:0040285`` (Excluded) → negative table;
* R4 ``0/m`` rows: merged with the other rows of the same (disease, HPO) first;
  only when every positive row reports zero (summed n = 0) does the pair move
  to the negative table;
* R5 ``alt_id`` → current ID, obsolete → ``replaced_by`` (dropped and logged
  when there is none);
* R6 keep only terms inside the ``HP:0000118`` subtree that reach the root;
  broken chains are logged;
* R7 duplicate (disease, HPO) rows are merged, frequency Σn/Σm, conflicting
  frequency categories are flagged.

Disease → gene comes from HPO ``genes_to_disease.txt`` unchanged, i.e. the same
gene universe as the old table, so the diff isolates the annotation cleaning.

A snapshot lives in ``PHENO_DATA_DIR/hpo_kb/{hp}_{hpoa}/``; the scorer uses it
only after ``hpo_kb/active.json`` points at it (``build_hpo_gene_table.py
--activate``), otherwise it keeps reading ``phenotype_to_genes.txt``.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import tempfile
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from ..config import PHENO_DATA_DIR

BUILD_VERSION = "hpo-gene-clean-v1"
PHENOTYPE_ROOT = "HP:0000118"
ALL_ROOT = "HP:0000001"
EXCLUDED = "HP:0040285"
KB_ROOT = Path(os.environ.get("NGS_UI_HPO_KB_DIR", str(PHENO_DATA_DIR / "hpo_kb")))
ACTIVE_POINTER = "active.json"
GENE_TABLE = "hpo_gene.tsv"
TERM_STATUS = "term_status.tsv"
LEGACY_TABLE = PHENO_DATA_DIR / "phenotype_to_genes.txt"
LEGACY_TABLE_ID = "legacy:phenotype_to_genes.txt"

# HPO frequency categories (lower bound %, label) used to compare frequencies.
FREQ_CODES = {
    "HP:0040280": 100.0,   # Obligate
    "HP:0040281": 90.0,    # Very frequent 80–99
    "HP:0040282": 55.0,    # Frequent 30–79
    "HP:0040283": 17.0,    # Occasional 5–29
    "HP:0040284": 2.5,     # Very rare 1–4
    EXCLUDED: 0.0,         # Excluded
}
_FRACTION = re.compile(r"^\s*(\d+)\s*/\s*(\d+)\s*$")
_PERCENT = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*%\s*$")


# ── ontology ──────────────────────────────────────────────────────────────
@dataclass
class Ontology:
    version: str = ""
    names: dict[str, str] = field(default_factory=dict)
    parents: dict[str, set[str]] = field(default_factory=lambda: defaultdict(set))
    alt: dict[str, str] = field(default_factory=dict)          # alt_id → primary id
    obsolete: set[str] = field(default_factory=set)
    replaced_by: dict[str, str] = field(default_factory=dict)
    _anc: dict[str, frozenset] = field(default_factory=dict, repr=False)

    def ancestors(self, term: str) -> frozenset:
        """Term plus all is_a ancestors (iterative, cycle-safe)."""
        hit = self._anc.get(term)
        if hit is not None:
            return hit
        seen = {term}
        stack = [term]
        while stack:
            for parent in self.parents.get(stack.pop(), ()):
                if parent not in seen:
                    seen.add(parent)
                    stack.append(parent)
        out = frozenset(seen)
        self._anc[term] = out
        return out

    def normalize(self, term: str) -> tuple[str | None, str]:
        """(current id or None, how): how ∈ current / alt_id / replaced_by / obsolete / unknown."""
        term = (term or "").strip()
        if term in self.alt:
            term, how = self.alt[term], "alt_id"
        else:
            how = "current"
        for _ in range(10):                      # follow replaced_by chains
            if term not in self.obsolete:
                break
            nxt = self.replaced_by.get(term)
            if not nxt:
                return None, "obsolete"
            term, how = self.alt.get(nxt, nxt), "replaced_by"
        if term not in self.names:
            return None, "unknown"
        return term, how

    def status(self, term: str) -> str:
        """phenotype / non_phenotype / broken_chain for a current id."""
        anc = self.ancestors(term)
        if ALL_ROOT not in anc:
            return "broken_chain"
        return "phenotype" if PHENOTYPE_ROOT in anc else "non_phenotype"


def parse_obo(path: Path) -> Ontology:
    onto = Ontology()
    cur = None
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("data-version:"):
                onto.version = line.split(":", 1)[1].strip().rsplit("/", 1)[-1]
            elif line == "[Term]":
                cur = None
            elif line.startswith("["):
                cur = "-"
            elif cur == "-":
                continue
            elif line.startswith("id: "):
                cur = line[4:].strip()
            elif cur is None:
                continue
            elif line.startswith("name: "):
                onto.names[cur] = line[6:].strip()
            elif line.startswith("is_a: "):
                onto.parents[cur].add(line[6:].split("!")[0].strip())
            elif line.startswith("alt_id: "):
                onto.alt[line[8:].split("!")[0].strip()] = cur
            elif line.startswith("is_obsolete:") and "true" in line:
                onto.obsolete.add(cur)
            elif line.startswith("replaced_by: "):
                onto.replaced_by[cur] = line[13:].split("!")[0].strip()
    return onto


# ── frequency ─────────────────────────────────────────────────────────────
def parse_frequency(value: str):
    """('fraction', n, m) / ('percent', p) / ('code', code) / ('blank',) / ('other', raw)."""
    value = (value or "").strip()
    if not value:
        return ("blank",)
    m = _FRACTION.match(value)
    if m:
        return ("fraction", int(m.group(1)), int(m.group(2)))
    m = _PERCENT.match(value)
    if m:
        return ("percent", float(m.group(1)))
    if value in FREQ_CODES:
        return ("code", value)
    return ("other", value)


def frequency_category(parsed) -> str | None:
    """Collapse a frequency to an HPO frequency category id for conflict checks."""
    kind = parsed[0]
    if kind == "code":
        return parsed[1]
    if kind == "fraction":
        n, m = parsed[1], parsed[2]
        if m <= 0:
            return None
        pct = 100.0 * n / m
    elif kind == "percent":
        pct = parsed[1]
    else:
        return None
    if pct <= 0:
        return EXCLUDED
    if pct < 5:
        return "HP:0040284"
    if pct < 30:
        return "HP:0040283"
    if pct < 80:
        return "HP:0040282"
    if pct < 100:
        return "HP:0040281"
    return "HP:0040280"


def _reports_presence(parsed) -> bool:
    """A positive row is evidence of presence unless it reports zero."""
    kind = parsed[0]
    if kind == "fraction":
        return parsed[1] > 0
    if kind == "percent":
        return parsed[1] > 0
    if kind == "code":
        return parsed[1] != EXCLUDED
    return True                                   # blank / unrecognised


# ── build ─────────────────────────────────────────────────────────────────
def _read_hpoa(path: Path):
    version = ""
    with open(path, encoding="utf-8") as handle:
        header = None
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("#"):
                if line.lower().startswith(("#version:", "#date:")):
                    version = line.split(":", 1)[1].strip()
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                continue
            yield version, dict(zip(header, fields))


def read_disease_genes(path: Path) -> dict[str, dict[str, str]]:
    """disease_id → {gene_symbol: ncbi_gene_id} from genes_to_disease.txt."""
    out: dict[str, dict[str, str]] = defaultdict(dict)
    with open(path, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = (row.get("gene_symbol") or "").strip()
            disease = (row.get("disease_id") or "").strip()
            if gene and disease and gene != "-":
                out[disease][gene] = (row.get("ncbi_gene_id") or "").replace("NCBIGene:", "").strip()
    return out


def _date_version(text: str) -> str:
    m = re.search(r"(\d{4}-\d{2}-\d{2})", text or "")
    return m.group(1) if m else (text or "unknown")


def build(hp_obo: Path, hpoa: Path, genes_to_disease: Path) -> dict:
    """Apply R1–R7 and propagate.  Returns everything needed to write a snapshot."""
    onto = parse_obo(hp_obo)
    disease_genes = read_disease_genes(genes_to_disease)
    qc: Counter = Counter()
    remapped: Counter = Counter()
    dropped_terms: Counter = Counter()
    broken: Counter = Counter()
    non_pheno: Counter = Counter()
    disease_names: dict[str, str] = {}
    positive: dict[tuple, list] = defaultdict(list)    # (disease, hpo) → rows
    negative: dict[tuple, list] = defaultdict(list)
    hpoa_version = ""

    for version, row in _read_hpoa(hpoa):
        hpoa_version = hpoa_version or version
        qc["rows_total"] += 1
        disease = (row.get("database_id") or "").strip()
        disease_names.setdefault(disease, (row.get("disease_name") or "").strip())
        if (row.get("aspect") or "").strip() != "P":                                    # R1
            qc["R1_non_P_aspect"] += 1
            continue
        raw_id = (row.get("hpo_id") or "").strip()
        term, how = onto.normalize(raw_id)                                              # R5
        if term is None:
            qc[f"R5_dropped_{how}"] += 1
            dropped_terms[raw_id] += 1
            continue
        if how != "current":
            qc[f"R5_{how}"] += 1
            remapped[f"{raw_id}→{term}"] += 1
        status = onto.status(term)                                                     # R6
        if status == "broken_chain":
            qc["R6_broken_chain"] += 1
            broken[term] += 1
            continue
        if status == "non_phenotype":
            qc["R6_outside_phenotype_subtree"] += 1
            non_pheno[term] += 1
            continue
        parsed = parse_frequency(row.get("frequency"))
        item = {"row": row, "freq": parsed, "raw_id": raw_id}
        if (row.get("qualifier") or "").strip().upper() == "NOT":                      # R2
            qc["R2_not"] += 1
            negative[(disease, term)].append(dict(item, reason="NOT"))
        elif parsed == ("code", EXCLUDED):                                              # R3
            qc["R3_excluded"] += 1
            negative[(disease, term)].append(dict(item, reason="excluded"))
        else:
            positive[(disease, term)].append(item)

    # R7 merge + R4 zero-only pairs
    merged: dict[tuple, dict] = {}
    for key, rows in positive.items():
        if len(rows) > 1:
            qc["R7_duplicate_pairs"] += 1
            qc["R7_rows_merged"] += len(rows) - 1
        fractions = [r["freq"] for r in rows if r["freq"][0] == "fraction"]
        n_sum = sum(f[1] for f in fractions)
        m_sum = sum(f[2] for f in fractions)
        zero_rows = sum(1 for r in rows if not _reports_presence(r["freq"]))
        qc["R4_zero_rows"] += zero_rows
        if not any(_reports_presence(r["freq"]) for r in rows):
            qc["R4_zero_pairs_to_negative"] += 1
            negative[key].extend(dict(r, reason="zero_frequency") for r in rows)
            continue
        if zero_rows:
            qc["R4_zero_rows_kept_merged"] += zero_rows
        cats = {c for c in (frequency_category(r["freq"]) for r in rows) if c}
        conflict = len(cats) > 1
        if conflict:
            qc["R7_frequency_conflict_pairs"] += 1
        merged[key] = {
            "rows": rows,
            "frequency": f"{n_sum}/{m_sum}" if fractions else next(
                (r["row"].get("frequency", "") for r in rows if r["row"].get("frequency")), ""),
            "conflict": conflict,
            "also_negative": key in negative,
        }
    for key in negative:
        if key in merged:
            qc["pairs_positive_and_negative"] += 1

    # propagate to ancestors inside the phenotype subtree
    term_genes: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    no_gene_diseases = set()
    in_subtree_cache: dict[str, frozenset] = {}
    for (disease, term) in merged:
        genes = disease_genes.get(disease)
        if not genes:
            no_gene_diseases.add(disease)
            continue
        anc = in_subtree_cache.get(term)
        if anc is None:
            anc = frozenset(a for a in onto.ancestors(term) if onto.status(a) == "phenotype")
            in_subtree_cache[term] = anc
        for a in anc:
            bucket = term_genes[a]
            for gene in genes:
                bucket[gene].add(disease)
    qc["positive_pairs"] = len(merged)
    qc["negative_pairs"] = len(negative)
    qc["diseases_with_annotations"] = len({d for d, _ in merged})
    qc["diseases_without_gene"] = len(no_gene_diseases)
    qc["terms_with_genes"] = len(term_genes)
    qc["term_gene_pairs"] = sum(len(g) for g in term_genes.values())
    return {
        "ontology": onto, "hpoa_version": _date_version(hpoa_version), "disease_genes": disease_genes,
        "disease_names": disease_names, "merged": merged, "negative": negative,
        "term_genes": term_genes, "qc": qc, "remapped": remapped, "dropped_terms": dropped_terms,
        "broken": broken, "non_phenotype": non_pheno,
    }


# ── snapshot I/O ──────────────────────────────────────────────────────────
def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _refs(rows) -> str:
    return ";".join(sorted({(r["row"].get("reference") or "").strip() for r in rows} - {""}))


def snapshot_name(result: dict) -> str:
    return f"{_date_version(result['ontology'].version)}_{result['hpoa_version']}"


def write_snapshot(result: dict, out_root: Path, sources: dict[str, Path], extra=None) -> Path:
    """Write the snapshot directory atomically; returns its path.

    ``extra(tmp_dir, manifest) -> str`` (the phase-1 knowledge base) may add
    files and manifest keys before the swap; its return value is appended to
    ``build_qc.md``.  It may raise to abort the whole snapshot.
    """
    onto: Ontology = result["ontology"]
    out_root.mkdir(parents=True, exist_ok=True)
    final = out_root / snapshot_name(result)
    tmp = Path(tempfile.mkdtemp(prefix=".build-", dir=out_root))
    try:
        with open(tmp / GENE_TABLE, "w", encoding="utf-8", newline="") as handle:
            w = csv.writer(handle, delimiter="\t", lineterminator="\n")
            w.writerow(["hpo_id", "hpo_name", "ncbi_gene_id", "gene_symbol", "disease_id"])
            for term in sorted(result["term_genes"]):
                for gene, diseases in sorted(result["term_genes"][term].items()):
                    ncbi = next((result["disease_genes"][d].get(gene, "") for d in sorted(diseases)), "")
                    w.writerow([term, onto.names.get(term, ""), ncbi, gene, ";".join(sorted(diseases))])
        with open(tmp / TERM_STATUS, "w", encoding="utf-8", newline="") as handle:
            w = csv.writer(handle, delimiter="\t", lineterminator="\n")
            w.writerow(["hpo_id", "status", "current_id", "name"])
            for term in sorted(onto.names):
                if term in onto.obsolete:
                    cur, how = onto.normalize(term)
                    w.writerow([term, "obsolete", cur or "", onto.names.get(term, "")])
                else:
                    w.writerow([term, onto.status(term), term, onto.names.get(term, "")])
            for alt, primary in sorted(onto.alt.items()):
                cur, how = onto.normalize(alt)
                w.writerow([alt, "alt_id", cur or "", onto.names.get(primary, "")])
        cols = ["database_id", "disease_name", "hpo_id", "hpo_name", "frequency", "frequency_raw",
                "reference", "evidence", "n_rows", "conflict", "original_hpo_ids"]
        with open(tmp / "disease_terms.tsv", "w", encoding="utf-8", newline="") as handle:
            w = csv.writer(handle, delimiter="\t", lineterminator="\n")
            w.writerow(cols + ["also_negative"])
            for (disease, term), item in sorted(result["merged"].items()):
                rows = item["rows"]
                w.writerow([disease, result["disease_names"].get(disease, ""), term, onto.names.get(term, ""),
                            item["frequency"], ";".join(r["row"].get("frequency", "") for r in rows),
                            _refs(rows), ";".join(sorted({r["row"].get("evidence", "") for r in rows} - {""})),
                            len(rows), "conflict" if item["conflict"] else "",
                            ";".join(sorted({r["raw_id"] for r in rows if r["raw_id"] != term})),
                            "yes" if item["also_negative"] else ""])
        with open(tmp / "negative_terms.tsv", "w", encoding="utf-8", newline="") as handle:
            w = csv.writer(handle, delimiter="\t", lineterminator="\n")
            w.writerow(["database_id", "disease_name", "hpo_id", "hpo_name", "reason", "frequency_raw",
                        "reference", "n_rows"])
            for (disease, term), rows in sorted(result["negative"].items()):
                w.writerow([disease, result["disease_names"].get(disease, ""), term, onto.names.get(term, ""),
                            ";".join(sorted({r["reason"] for r in rows})),
                            ";".join(r["row"].get("frequency", "") for r in rows), _refs(rows), len(rows)])
        manifest = {
            "build_version": BUILD_VERSION,
            "built_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "hp_version": onto.version, "hpoa_version": result["hpoa_version"],
            "sources": {k: {"path": str(p), "sha256": _sha256(p), "size": p.stat().st_size}
                        for k, p in sources.items() if p and Path(p).is_file()},
            "qc": dict(result["qc"]),
            "disease_gene_source": "genes_to_disease.txt (all association types, unchanged from the old table)",
        }
        qc_extra = extra(tmp, manifest) if extra else ""
        (tmp / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=1) + "\n",
                                           encoding="utf-8")
        (tmp / "build_qc.md").write_text(render_qc(result, manifest) + (qc_extra or ""), encoding="utf-8")
        if final.exists():
            backup = final.with_name(f"{final.name}.replaced-{int(time.time())}")
            final.rename(backup)
        tmp.rename(final)
    except BaseException:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    return final


def render_qc(result: dict, manifest: dict) -> str:
    qc = result["qc"]
    onto: Ontology = result["ontology"]
    lines = [f"# HPO → 基因查表建置紀錄（{BUILD_VERSION}）", "",
             f"- hp.obo：{manifest['hp_version']}；phenotype.hpoa：{manifest['hpoa_version']}",
             f"- 疾病–基因：{manifest['disease_gene_source']}", "",
             "## 各規則影響列數", "", "| 項目 | 數量 |", "|---|---|"]
    labels = [
        ("rows_total", "phenotype.hpoa 總列數"),
        ("R1_non_P_aspect", "R1 非表型 aspect（I/C/M/H）捨棄"),
        ("R2_not", "R2 qualifier=NOT → negative"),
        ("R3_excluded", "R3 頻率 Excluded（HP:0040285）→ negative"),
        ("R4_zero_rows", "R4 頻率為 0 的列（0/m、0%）"),
        ("R4_zero_rows_kept_merged", "R4 其中與同一（疾病, HPO）其他 n>0 列合併後保留"),
        ("R4_zero_pairs_to_negative", "R4 合計為 0 的（疾病, HPO）→ negative"),
        ("R5_alt_id", "R5 alt_id 改為現行 ID"),
        ("R5_replaced_by", "R5 obsolete 依 replaced_by 替換"),
        ("R5_dropped_obsolete", "R5 obsolete 無替換，捨棄"),
        ("R5_dropped_unknown", "R5 不在 hp.obo，捨棄"),
        ("R6_outside_phenotype_subtree", "R6 不在 HP:0000118 子樹，捨棄"),
        ("R6_broken_chain", "R6 到不了根節點（斷鏈），捨棄"),
        ("R7_duplicate_pairs", "R7 重複的（疾病, HPO）組數"),
        ("R7_rows_merged", "R7 合併掉的列數"),
        ("R7_frequency_conflict_pairs", "R7 頻率類別互相矛盾的組數（標記 conflict）"),
        ("pairs_positive_and_negative", "同一（疾病, HPO）同時有支持與排除註記"),
        ("positive_pairs", "清理後（疾病, HPO）支持註記"),
        ("negative_pairs", "negative 表（疾病, HPO）"),
        ("diseases_with_annotations", "有支持註記的疾病"),
        ("diseases_without_gene", "其中沒有基因對應的疾病"),
        ("terms_with_genes", "有基因的 HPO term（含往上層傳遞）"),
        ("term_gene_pairs", "HPO term–基因配對"),
    ]
    for key, label in labels:
        lines.append(f"| {label} | {qc.get(key, 0):,} |")
    if result["broken"]:
        lines += ["", "## 斷鏈 term（R6）", ""]
        lines += [f"- {t} {onto.names.get(t, '')}（{n} 列）" for t, n in result["broken"].most_common()]
    if result["dropped_terms"]:
        lines += ["", "## 捨棄的 obsolete／未知 ID（R5，前 30）", ""]
        lines += [f"- {t}（{n} 列）" for t, n in result["dropped_terms"].most_common(30)]
    if result["non_phenotype"]:
        lines += ["", "## aspect=P 但不在表型子樹的 term（R6，前 30）", ""]
        lines += [f"- {t} {onto.names.get(t, '')}（{n} 列）" for t, n in result["non_phenotype"].most_common(30)]
    return "\n".join(lines) + "\n"


# ── activation / lookup for the scorer ────────────────────────────────────
def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        handle.write(text)
    os.replace(tmp, path)


def activate(snapshot: Path, kb_root: Path | None = None) -> dict:
    kb_root = Path(kb_root or KB_ROOT)
    snapshot = Path(snapshot).resolve()
    if not (snapshot / GENE_TABLE).is_file() or not (snapshot / TERM_STATUS).is_file():
        raise FileNotFoundError(f"{snapshot} is not a complete snapshot")
    manifest = json.loads((snapshot / "manifest.json").read_text(encoding="utf-8"))
    pointer = {"snapshot": snapshot.name, "build_version": manifest.get("build_version"),
               "activated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    _atomic_write(kb_root / ACTIVE_POINTER, json.dumps(pointer, ensure_ascii=False, indent=1) + "\n")
    return pointer


def deactivate(kb_root: Path | None = None) -> bool:
    kb_root = Path(kb_root or KB_ROOT)
    path = kb_root / ACTIVE_POINTER
    if path.exists():
        path.unlink()
        return True
    return False


@dataclass(frozen=True)
class ActiveTable:
    table_id: str                 # stamped next to pheno_score.tsv
    gene_table: Path
    term_status: Path | None      # None for the legacy table
    label: str                    # shown in the UI


def active_table(kb_root: Path | None = None, legacy: Path | None = None) -> ActiveTable:
    """The table the scorer should use: the activated snapshot, else the official file."""
    kb_root = Path(kb_root or KB_ROOT)
    legacy = Path(legacy or LEGACY_TABLE)
    pointer = kb_root / ACTIVE_POINTER
    try:
        data = json.loads(pointer.read_text(encoding="utf-8"))
        snap = kb_root / str(data.get("snapshot") or "")
        if data.get("snapshot") and (snap / GENE_TABLE).is_file() and (snap / TERM_STATUS).is_file():
            return ActiveTable(f"{data.get('build_version') or BUILD_VERSION}:{snap.name}",
                               snap / GENE_TABLE, snap / TERM_STATUS, f"清理後 HPO 查表 {snap.name}")
    except (OSError, ValueError):
        pass
    return ActiveTable(LEGACY_TABLE_ID, legacy, None, "phenotype_to_genes.txt")


def read_term_status(path: Path | None) -> tuple[dict[str, str], set[str]]:
    """({any id → current id}, {current ids that are not phenotype terms})."""
    aliases: dict[str, str] = {}
    non_phenotype: set[str] = set()
    if not path or not Path(path).is_file():
        return aliases, non_phenotype
    with open(path, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hid, status, cur = row.get("hpo_id") or "", row.get("status") or "", row.get("current_id") or ""
            if status in ("alt_id", "obsolete"):
                if cur:
                    aliases[hid] = cur
                else:
                    non_phenotype.add(hid)          # obsolete without replacement: cannot score
            elif status in ("non_phenotype", "broken_chain"):
                non_phenotype.add(hid)
    return aliases, non_phenotype


def read_gene_table(path: Path, canonical=lambda g: g) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    if not Path(path).is_file():
        return out
    with open(path, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hid = (row.get("hpo_id") or "").strip()
            gene = canonical(row.get("gene_symbol") or "")
            if hid and gene and gene != "-":
                out[hid].add(gene)
    return out
