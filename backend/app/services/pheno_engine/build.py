"""Phase-1 knowledge-base build hook (plan §4.2) shared by the CLI and the HPO update job.

``make_extra(...)`` returns the hook ``hpo_gene_table.write_snapshot`` calls
to add the knowledge-base files to a snapshot.
"""
from __future__ import annotations

import csv
import json
import time
from pathlib import Path

from ...config import PHENO_DATA_DIR


class KBBuildError(RuntimeError):
    """The knowledge base cannot be built (missing source or failed QC)."""


def platform_omim_rows() -> tuple[list[dict], dict]:
    from .. import omim_store
    omim_store.ensure_loaded()
    state = omim_store._state
    rows = list({id(r): r for r in list(state.get("by_omim_id", {}).values())
                 + list(state.get("by_gene", {}).values())}.values())
    return rows, {"source": omim_store.source(), "signature": str(omim_store.cache_signature()), "genes": len(rows)}


def default_mondo() -> Path | None:
    from ... import config
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
    from . import disease_gene as DG, kb as KBM, mondo as MD

    def extra(tmp: Path, manifest: dict) -> str:
        started = time.time()
        kb = KBM.build_kb(result)
        progress(f"知識庫：{kb.N:,} 個疾病、{kb.T:,} 個表型 term；計算基準分 μ(D) …")
        kernel = KBM.Kernel(kb)
        mu, pool = KBM.baseline(kb, kernel, progress=lambda k, n: progress(f"  基準分 {k:,}/{n:,}"))
        qc = KBM.qc(kb, result["ontology"], mu)
        if not qc["passed"]:
            raise KBBuildError(f"QC 未通過：{json.dumps(qc['checks'], ensure_ascii=False)}；"
                             f"斷鏈 {qc['unreachable']}；IC 反序 {qc['ic_violations']}")
        files = KBM.write_files(kb, mu, tmp)

        universe = set(kb.diseases)
        links, sources = [], {}
        if omim_rows:
            links += DG.omim_links(omim_rows)
            sources["omim"] = omim_info
        elif "omim" not in allow_missing:
            raise KBBuildError("讀不到平台 OMIM 資料；若確定要略過，加 --allow-missing omim")
        if orphanet and Path(orphanet).is_file():
            orpha, info = DG.orphanet_links(Path(orphanet))
            # ORPHA diseases newer than this product6 (e.g. ORPHA:714652 PCDH19 clustering
            # epilepsy): HPO already annotates them, so take their genes from genes_to_disease.
            newer = DG.orphanet_links_from_g2d(result["disease_genes"], skip={r["disease_id"] for r in orpha})
            links += orpha + newer
            added = sorted({r["disease_id"] for r in newer} & universe)
            sources["orphanet"] = {"path": str(orphanet), "sha256": KBM._sha(orphanet), **info,
                                   "added_from_genes_to_disease": len(added)}
            if added:
                progress(f"Orphanet product6 沒有的 {len(added)} 個 ORPHA 疾病（HPO 已有表型註記），"
                         f"基因改用 HPO genes_to_disease.txt")
        elif "orphanet" not in allow_missing:
            # Server cannot reach orphadata.com: the ORPHA rows of HPO's
            # genes_to_disease.txt come from the same product6 file.
            links += DG.orphanet_links_from_g2d(result["disease_genes"])
            sources["orphanet"] = {"fallback": True, "source": "HPO genes_to_disease.txt",
                                   "version": result.get("hpoa_version", "")}
            progress("Orphanet product6 不在本機：改用 HPO genes_to_disease.txt 內的 Orphanet 對應")
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
            raise KBBuildError("找不到 MONDO（--mondo 或 --download-mondo）；略過請加 --allow-missing mondo")
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
    if (kb["sources"].get("orphanet") or {}).get("fallback"):
        lines.append("- Orphanet 對應取自 HPO genes_to_disease.txt（伺服器下載不到 product6，無法區分關聯類型，"
                     "全部視為致病關聯；ORPHA 交叉核對因此一致）")
    return "\n".join(lines) + "\n"
