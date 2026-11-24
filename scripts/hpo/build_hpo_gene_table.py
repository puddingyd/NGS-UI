#!/usr/bin/env python3
"""Build (and optionally activate) the cleaned HPO → gene table (plan §10.7, phase 0).

Reads hp.obo, phenotype.hpoa and genes_to_disease.txt (default: the files in
NGS_UI_PHENO_DATA_DIR), applies rules R1–R7 of plan §4.2.3, propagates to
ancestors and writes a snapshot under ``phenotype_data/hpo_kb/{hp}_{hpoa}/``
(hpo_gene.tsv, term_status.tsv, disease_terms.tsv, negative_terms.tsv,
manifest.json, build_qc.md).

Building does not change scoring.  The scorer switches only after
``--activate`` writes ``hpo_kb/active.json``; restart the service (or it is
picked up on the next start) and each case's pheno_score.tsv is recomputed the
next time it is opened.  ``--deactivate`` goes back to phenotype_to_genes.txt.

    .venv/bin/python scripts/hpo/build_hpo_gene_table.py                 # build only
    .venv/bin/python scripts/hpo/build_hpo_gene_table.py --activate      # build + activate
    .venv/bin/python scripts/hpo/build_hpo_gene_table.py --activate-only hpo_kb/2026-09-01_2026-09-01
    .venv/bin/python scripts/hpo/build_hpo_gene_table.py --deactivate
    .venv/bin/python scripts/hpo/build_hpo_gene_table.py --status
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))

from app.config import PHENO_DATA_DIR  # noqa: E402
from app.services import hpo_gene_table as H  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--hp-obo", type=Path, default=PHENO_DATA_DIR / "hp.obo")
    ap.add_argument("--hpoa", type=Path, default=PHENO_DATA_DIR / "phenotype.hpoa")
    ap.add_argument("--genes-to-disease", type=Path, default=PHENO_DATA_DIR / "genes_to_disease.txt")
    ap.add_argument("--kb-root", type=Path, default=H.KB_ROOT)
    ap.add_argument("--activate", action="store_true", help="activate the snapshot after building")
    ap.add_argument("--activate-only", type=Path, default=None, help="activate an existing snapshot directory")
    ap.add_argument("--deactivate", action="store_true", help="go back to phenotype_to_genes.txt")
    ap.add_argument("--status", action="store_true", help="show the table the scorer will use")
    args = ap.parse_args(argv)

    if args.status:
        table = H.active_table(args.kb_root)
        print(f"使用中：{table.label}（{table.table_id}）\n檔案：{table.gene_table}")
        return 0
    if args.deactivate:
        print("已停用清理後查表，改回 phenotype_to_genes.txt" if H.deactivate(args.kb_root) else "原本就沒有啟用")
        print("重新啟動服務後生效。")
        return 0
    if args.activate_only:
        snap = args.activate_only if args.activate_only.is_absolute() else (args.kb_root / args.activate_only.name)
        print(json.dumps(H.activate(snap, args.kb_root), ensure_ascii=False))
        print("重新啟動服務後生效。")
        return 0

    missing = [str(p) for p in (args.hp_obo, args.hpoa, args.genes_to_disease) if not p.is_file()]
    if missing:
        print("找不到：" + "、".join(missing), file=sys.stderr)
        print("phenotype.hpoa 與 genes_to_disease.txt 可從 HPO release 下載："
              "https://github.com/obophenotype/human-phenotype-ontology/releases", file=sys.stderr)
        return 2
    result = H.build(args.hp_obo, args.hpoa, args.genes_to_disease)
    snap = H.write_snapshot(result, args.kb_root, {
        "hp_obo": args.hp_obo, "hpoa": args.hpoa, "genes_to_disease": args.genes_to_disease})
    qc = result["qc"]
    print(f"快照：{snap}")
    print(f"清理後支持註記 {qc['positive_pairs']:,}、negative {qc['negative_pairs']:,}、"
          f"term–基因配對 {qc['term_gene_pairs']:,}；詳見 {snap / 'build_qc.md'}")
    if args.activate:
        print(json.dumps(H.activate(snap, args.kb_root), ensure_ascii=False))
        print("已啟用；重新啟動服務後生效。")
    else:
        print("尚未啟用（評分仍用 phenotype_to_genes.txt）。先跑 compare_pheno_tables.py 看差異。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
