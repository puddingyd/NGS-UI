#!/usr/bin/env python3
"""Build the phenotype knowledge base snapshot (plan §4.2, phase 1).

One snapshot ``phenotype_data/hpo_kb/{hp}_{hpoa}/`` holds the phase-0 cleaned
HPO → gene table plus:

    terms.tsv, ancestors.npz        IC / IC* per term, ancestor index
    disease_closure.npz             direct and ancestor-closed annotations, μ(D)
    disease_baseline.tsv            baseline score μ(D) per disease
    disease_gene.tsv                OMIM (platform DB) + Orphanet product6 links, use/flags
    disease_gene_discordant.tsv     links that disagree with HPO genes_to_disease.txt
    disease_groups.tsv              MONDO equivalence groups (OMIM ↔ Orphanet)
    build_qc.md, manifest.json      QC (§4.2.5) and versions / content hashes

Building does not change scoring; activation is the same as phase 0
(``build_hpo_gene_table.py --activate-only``).  A failed QC check aborts the
build and leaves the old snapshot in place.

    .venv/bin/python scripts/build_hpo_kb.py --release 2026-02-16 --download-orphanet --download-mondo
    .venv/bin/python scripts/build_hpo_kb.py --hpoa X --genes-to-disease Y --orphanet en_product6.xml --mondo mondo.obo

Needs numpy.  The baseline step takes a few minutes.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "backend"))

from app.config import PHENO_DATA_DIR  # noqa: E402
from app.services import hpo_gene_table as H  # noqa: E402

ORPHANET_URL = "https://www.orphadata.com/data/xml/en_product6.xml"
MONDO_URL = "https://purl.obolibrary.org/obo/mondo.obo"
_PHASE0 = REPO / "scripts/hpo/build_hpo_gene_table.py"


def _phase0_module():
    spec = importlib.util.spec_from_file_location("build_hpo_gene_table", _PHASE0)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _platform_omim_rows() -> tuple[list[dict], dict]:
    from app.services import omim_store
    omim_store.ensure_loaded()
    state = omim_store._state
    rows = list({id(r): r for r in list(state.get("by_omim_id", {}).values())
                 + list(state.get("by_gene", {}).values())}.values())
    return rows, {"source": omim_store.source(), "signature": str(omim_store.cache_signature()), "genes": len(rows)}


def _default_mondo() -> Path | None:
    from app import config
    for p in (PHENO_DATA_DIR / "mondo" / "mondo.obo", Path(getattr(config, "GENE_DISEASE_RAW_DIR", "")) / "mondo.json"):
        if p and p.is_file():
            return p
    return None


def _write_tsv(path: Path, header, rows) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def make_extra(result: dict, *, omim_rows, omim_info, orphanet: Path | None, mondo: Path | None,
               allow_missing: set, progress=print):
    """The phase-1 hook for ``hpo_gene_table.write_snapshot``."""
    from app.services.pheno_engine import disease_gene as DG, kb as KBM, mondo as MD

    def extra(tmp: Path, manifest: dict) -> str:
        started = time.time()
        kb = KBM.build_kb(result)
        progress(f"知識庫：{kb.N:,} 個疾病、{kb.T:,} 個表型 term；計算基準分 μ(D) …")
        kernel = KBM.Kernel(kb)
        mu, pool = KBM.baseline(kb, kernel, progress=lambda k, n: progress(f"  基準分 {k:,}/{n:,}"))
        qc = KBM.qc(kb, result["ontology"], mu)
        if not qc["passed"]:
            raise SystemExit(f"QC 未通過：{json.dumps(qc['checks'], ensure_ascii=False)}；"
                             f"斷鏈 {qc['unreachable']}；IC 反序 {qc['ic_violations']}")
        files = KBM.write_files(kb, mu, tmp)

        universe = set(kb.diseases)
        links, sources = [], {}
        if omim_rows:
            links += DG.omim_links(omim_rows)
            sources["omim"] = omim_info
        elif "omim" not in allow_missing:
            raise SystemExit("讀不到平台 OMIM 資料；若確定要略過，加 --allow-missing omim")
        if orphanet and Path(orphanet).is_file():
            orpha, info = DG.orphanet_links(Path(orphanet))
            links += orpha
            sources["orphanet"] = {"path": str(orphanet), "sha256": KBM._sha(orphanet), **info}
        elif "orphanet" not in allow_missing:
            raise SystemExit("找不到 Orphanet product6（--orphanet 或 --download-orphanet）；略過請加 --allow-missing orphanet")
        links.sort(key=lambda r: (r["disease_id"], r["gene_symbol"], r["source"]))
        _write_tsv(tmp / "disease_gene.tsv", DG.COLUMNS, [[r[c] for c in DG.COLUMNS] for r in links])
        discordant, cross = DG.cross_check(links, result["disease_genes"], universe)
        _write_tsv(tmp / "disease_gene_discordant.tsv", ["disease_id", "gene_symbol", "in_kb", "in_genes_to_disease"],
                   discordant)
        files += ["disease_gene.tsv", "disease_gene_discordant.tsv"]

        group_stats = {"grouped_diseases": 0, "groups": 0}
        group_rows = []
        if mondo and Path(mondo).is_file():
            groups, labels, minfo = MD.parse(Path(mondo))
            group_rows, group_stats = MD.disease_groups(groups, universe)
            sources["mondo"] = {"path": str(mondo), "sha256": KBM._sha(mondo), **minfo}
            group_rows = [(g, labels.get(g, ""), d) for g, d in group_rows]
        elif "mondo" not in allow_missing:
            raise SystemExit("找不到 MONDO（--mondo 或 --download-mondo）；略過請加 --allow-missing mondo")
        _write_tsv(tmp / "disease_groups.tsv", ["mondo_id", "mondo_label", "disease_id"], group_rows)
        files.append("disease_groups.tsv")

        used = [r for r in links if r["use"]]
        genes_with_disease = {r["gene_symbol"] for r in used if r["disease_id"] in universe}
        diseases_with_gene = {r["disease_id"] for r in used} & universe
        data_files = [p.name for p in tmp.iterdir() if p.is_file() and p.name not in ("manifest.json", "build_qc.md")]
        manifest["kb"] = {
            "kb_version": KBM.KB_VERSION,
            "lambda": KBM.LAMBDA,
            "baseline_pool": "all directly annotated phenotype terms, equal weight",
            "pool_size": pool,
            "qc": qc,
            "sources": sources,
            "disease_gene": {"links": len(links), "used": len(used), "genes_with_disease": len(genes_with_disease),
                             "diseases_with_gene": len(diseases_with_gene), "cross_check": cross},
            "groups": group_stats,
            "content_sha256": KBM.content_hashes(tmp, data_files),
            "build_seconds": round(time.time() - started, 1),
        }
        return render_kb_qc(manifest["kb"])

    return extra


def render_kb_qc(kb: dict) -> str:
    qc = kb["qc"]
    lines = ["", f"## 知識庫（{kb['kb_version']}）", "",
             f"- 疾病全集 N = {qc['N']:,}（OMIM + Orphanet，有清理後表型註記；不含 DECIPHER）",
             f"- 表型 term {qc['terms']:,}：有 IC {qc['terms_with_ic']:,}、沒有疾病註記（p=0，用 IC*）{qc['terms_p0']:,}",
             f"- 基準分症狀池 {kb['pool_size']:,} 個直接註記的表型 term；旁系係數上限 λ = {kb['lambda']}",
             "", "| QC 檢查（§4.2.5） | 結果 |", "|---|---|"]
    labels = {"all_terms_reach_root": "每個 term 都到得了根節點", "ic_ancestor_le_descendant": "IC(祖先) ≤ IC(後代)",
              "root_ic_zero": "根節點 HP:0000118 的 IC = 0", "ic_finite": "IC／IC* 沒有 NaN 或無窮大",
              "baseline_finite": "基準分 μ(D) 都有限"}
    for key, ok in qc["checks"].items():
        lines.append(f"| {labels.get(key, key)} | {'通過' if ok else '未通過'} |")
    dg = kb["disease_gene"]
    lines += ["", "### 疾病–基因對應", "",
              f"- 對應 {dg['links']:,} 筆，採用 {dg['used']:,} 筆；有基因的疾病 {dg['diseases_with_gene']:,}、有疾病的基因 {dg['genes_with_disease']:,}",
              "", "| 來源 | 本知識庫 | genes_to_disease | 兩邊都有 | 只在本知識庫 | 只在 genes_to_disease |",
              "|---|---|---|---|---|---|"]
    for prefix, s in dg["cross_check"].items():
        lines.append(f"| {prefix} | {s['kb']:,} | {s['genes_to_disease']:,} | {s['both']:,} | {s['kb_only']:,} | "
                     f"{s['genes_to_disease_only']:,} |")
    g = kb["groups"]
    lines += ["", "### MONDO 分組", "",
              f"- {g.get('grouped_diseases', 0):,} 個疾病分到 {g.get('groups', 0):,} 組；"
              f"其中同時含 OMIM 與 Orphanet 的組 {g.get('groups_omim_and_orpha', 0):,}；"
              f"對應到多個 MONDO 而略過的疾病 {g.get('ambiguous_diseases', 0):,}",
              "", "- 各來源版本：" + "；".join(f"{k} {v.get('version') or v.get('source') or ''}".strip()
                                         for k, v in kb["sources"].items())]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--hp-obo", type=Path, default=PHENO_DATA_DIR / "hp.obo")
    ap.add_argument("--hpoa", type=Path, default=PHENO_DATA_DIR / "phenotype.hpoa")
    ap.add_argument("--genes-to-disease", type=Path, default=PHENO_DATA_DIR / "genes_to_disease.txt")
    ap.add_argument("--release", default=None, help="use (and download if needed) this HPO release, e.g. 2026-02-16")
    ap.add_argument("--orphanet", type=Path, default=PHENO_DATA_DIR / "orphanet" / "en_product6.xml")
    ap.add_argument("--download-orphanet", action="store_true", help="download en_product6.xml first")
    ap.add_argument("--mondo", type=Path, default=None, help="mondo.obo (preferred) or mondo.json")
    ap.add_argument("--download-mondo", action="store_true", help="download mondo.obo first")
    ap.add_argument("--allow-missing", default="", help="comma list of omim,orphanet,mondo to build without")
    ap.add_argument("--kb-root", type=Path, default=H.KB_ROOT)
    args = ap.parse_args(argv)
    allow = {x.strip() for x in args.allow_missing.split(",") if x.strip()}

    phase0 = _phase0_module()
    if args.release:
        folder = phase0.fetch_release(args.release)
        args.hp_obo, args.hpoa, args.genes_to_disease = (folder / n for n in phase0.RELEASE_FILES)
    if args.download_orphanet:
        args.orphanet.parent.mkdir(parents=True, exist_ok=True)
        print(f"下載 {ORPHANET_URL} …")
        phase0._download(ORPHANET_URL, args.orphanet)
    if args.download_mondo:
        args.mondo = PHENO_DATA_DIR / "mondo" / "mondo.obo"
        args.mondo.parent.mkdir(parents=True, exist_ok=True)
        print(f"下載 {MONDO_URL} …")
        phase0._download(MONDO_URL, args.mondo)
    mondo = args.mondo or _default_mondo()
    missing = [str(p) for p in (args.hp_obo, args.hpoa, args.genes_to_disease) if not p.is_file()]
    if missing:
        print("找不到：" + "、".join(missing), file=sys.stderr)
        return 2
    try:
        omim_rows, omim_info = _platform_omim_rows()
    except Exception as exc:  # noqa: BLE001
        print(f"[warn] 讀平台 OMIM 失敗：{type(exc).__name__}: {exc}", file=sys.stderr)
        omim_rows, omim_info = [], {}
    result = H.build(args.hp_obo, args.hpoa, args.genes_to_disease)
    extra = make_extra(result, omim_rows=omim_rows, omim_info=omim_info,
                       orphanet=args.orphanet if args.orphanet.is_file() else None, mondo=mondo,
                       allow_missing=allow)
    sources = {"hp_obo": args.hp_obo, "hpoa": args.hpoa, "genes_to_disease": args.genes_to_disease}
    snap = H.write_snapshot(result, args.kb_root, sources, extra=extra)
    print(f"快照：{snap}\nQC 與統計：{snap / 'build_qc.md'}")
    print(f"啟用（評分改用此快照的 HPO → 基因查表）：.venv/bin/python scripts/hpo/build_hpo_gene_table.py "
          f"--activate-only {snap.name}，再重新啟動服務")
    return 0


if __name__ == "__main__":
    sys.exit(main())
