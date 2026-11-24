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

Updating HPO: ``--release v2026-09-01`` (or ``latest``) downloads hp.obo,
phenotype.hpoa and genes_to_disease.txt of one HPO release into
``phenotype_data/hpo_release_{version}/`` and builds from them; with
``--activate`` it also installs that hp.obo as ``phenotype_data/hp.obo`` (the
HPO search reads it; the old file is kept as ``hp.obo.bak-{time}``).

    .venv/bin/python scripts/hpo/build_hpo_gene_table.py --release latest
    .venv/bin/python scripts/hpo/build_hpo_gene_table.py --release 2026-09-01 --activate
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))

from app.config import PHENO_DATA_DIR  # noqa: E402
from app.services import hpo_gene_table as H  # noqa: E402


RELEASE_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/{part}/{name}"
RELEASE_FILES = ("hp.obo", "phenotype.hpoa", "genes_to_disease.txt")


def _download(url: str, dest: Path) -> None:
    tmp = dest.with_name(dest.name + ".part")
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "ngs-ui-hpo/1.0"})
            with urllib.request.urlopen(request, timeout=120) as response, open(tmp, "wb") as out:
                shutil.copyfileobj(response, out)
            tmp.replace(dest)
            return
        except OSError:
            if attempt == 3:
                raise
            time.sleep(2 ** (attempt + 1))


def fetch_release(version: str, data_dir: Path = PHENO_DATA_DIR) -> Path:
    """Download one HPO release's three files; returns the directory holding them."""
    version = version.strip().lstrip("v")
    if version == "latest":
        staging = data_dir / "hpo_release_latest.part"
        staging.mkdir(parents=True, exist_ok=True)
        for name in RELEASE_FILES:
            print(f"下載 latest/{name} …")
            _download(RELEASE_URL.format(part="latest/download", name=name), staging / name)
        found = H.parse_obo(staging / "hp.obo").version
        version = H._date_version(found)
        final = data_dir / f"hpo_release_{version}"
        if final.exists():
            shutil.rmtree(staging)
        else:
            staging.rename(final)
        return final
    target = data_dir / f"hpo_release_{version}"
    target.mkdir(parents=True, exist_ok=True)
    for name in RELEASE_FILES:
        if not (target / name).is_file():
            print(f"下載 v{version}/{name} …")
            _download(RELEASE_URL.format(part=f"download/v{version}", name=name), target / name)
    return target


def install_obo(source: Path, data_dir: Path = PHENO_DATA_DIR) -> Path | None:
    live = data_dir / "hp.obo"
    if live.is_file() and live.read_bytes() == source.read_bytes():
        return None
    if live.is_file():
        backup = live.with_name(f"hp.obo.bak-{time.strftime('%Y%m%d%H%M%S')}")
        shutil.copy2(live, backup)
    tmp = live.with_name(".hp.obo.part")
    shutil.copy2(source, tmp)
    tmp.replace(live)
    return live


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--hp-obo", type=Path, default=PHENO_DATA_DIR / "hp.obo")
    ap.add_argument("--hpoa", type=Path, default=PHENO_DATA_DIR / "phenotype.hpoa")
    ap.add_argument("--genes-to-disease", type=Path, default=PHENO_DATA_DIR / "genes_to_disease.txt")
    ap.add_argument("--kb-root", type=Path, default=H.KB_ROOT)
    ap.add_argument("--release", default=None,
                    help="download this HPO release (e.g. 2026-09-01 or latest) and build from it")
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
        try:
            source = Path(json.loads((snap / "manifest.json").read_text(encoding="utf-8"))
                          ["sources"]["hp_obo"]["path"])
        except (OSError, ValueError, KeyError, TypeError):
            source = None
        live = PHENO_DATA_DIR / "hp.obo"
        if source and source.is_file() and source.resolve() != live.resolve() and install_obo(source):
            print(f"已更新 {live}（與此快照同版本；舊檔另存 hp.obo.bak-時間）")
        print("重新啟動服務後生效。")
        return 0

    release_dir = None
    if args.release:
        release_dir = fetch_release(args.release)
        args.hp_obo, args.hpoa, args.genes_to_disease = (release_dir / n for n in RELEASE_FILES)
        print(f"HPO 檔案：{release_dir}")
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
        if release_dir is not None and install_obo(args.hp_obo):
            print(f"已更新 {PHENO_DATA_DIR / 'hp.obo'}（舊檔另存 hp.obo.bak-時間）")
        print("已啟用；重新啟動服務後生效。")
    else:
        print("尚未啟用（評分仍用 phenotype_to_genes.txt）。先跑 compare_pheno_tables.py 看差異。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
