#!/usr/bin/env python3
"""Build the phenotype knowledge base snapshot (plan §4.2, phase 1).

The home-page "更新 HPO" button does the same thing (download the latest HPO,
build, activate); this command is for building a specific release by hand.

One snapshot ``phenotype_data/hpo_kb/{hp}_{hpoa}/`` holds the phase-0 cleaned
HPO → gene table plus:

    terms.tsv, ancestors.npz        IC / IC* per term, ancestor index
    disease_closure.npz             direct and ancestor-closed annotations, μ(D)
    disease_baseline.tsv            baseline score μ(D) per disease
    disease_gene.tsv                OMIM (platform DB) + Orphanet links, use/flags
    disease_gene_discordant.tsv     links that disagree with HPO genes_to_disease.txt
    disease_groups.tsv              MONDO equivalence groups (OMIM ↔ Orphanet)
    build_qc.md, manifest.json      QC (§4.2.5) and versions / content hashes

Orphanet product6 and MONDO are downloaded when ``--download-orphanet`` /
``--download-mondo`` is given; a failed download keeps the old local copy, and
without any Orphanet file the Orphanet rows of HPO genes_to_disease.txt are
used.  The QC summary is printed at the end.  ``--activate`` switches scoring
to the new snapshot (the running service picks it up within a minute).

    .venv/bin/python scripts/build_hpo_kb.py --release 2026-02-16 --download-orphanet --download-mondo

Needs numpy.  The baseline step takes about a minute.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "backend"))

from app.config import PHENO_DATA_DIR  # noqa: E402
from app.services import hpo_gene_table as H, hpo_update as U  # noqa: E402
from app.services.pheno_engine import build as B  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--hp-obo", type=Path, default=PHENO_DATA_DIR / "hp.obo")
    ap.add_argument("--hpoa", type=Path, default=PHENO_DATA_DIR / "phenotype.hpoa")
    ap.add_argument("--genes-to-disease", type=Path, default=PHENO_DATA_DIR / "genes_to_disease.txt")
    ap.add_argument("--release", default=None, help="use (and download if needed) this HPO release, e.g. 2026-02-16")
    ap.add_argument("--orphanet", type=Path, default=PHENO_DATA_DIR / "orphanet" / "en_product6.xml")
    ap.add_argument("--download-orphanet", action="store_true", help="try to download en_product6.xml first")
    ap.add_argument("--mondo", type=Path, default=None, help="mondo.obo (preferred) or mondo.json")
    ap.add_argument("--download-mondo", action="store_true", help="try to download mondo.obo first")
    ap.add_argument("--allow-missing", default="", help="comma list of omim,mondo to build without")
    ap.add_argument("--activate", action="store_true", help="also switch scoring to this snapshot")
    ap.add_argument("--kb-root", type=Path, default=H.KB_ROOT)
    args = ap.parse_args(argv)
    allow = {x.strip() for x in args.allow_missing.split(",") if x.strip()}

    if args.release:
        folder = U.fetch_release(args.release)
        args.hp_obo, args.hpoa, args.genes_to_disease = (folder / n for n in U.RELEASE_FILES)
    warnings = []
    if args.download_orphanet:
        path, warn = U.fetch_optional("Orphanet", U.ORPHANET_URLS, args.orphanet, max_age_days=0)
        warnings += [warn] if warn else []
    if args.download_mondo:
        args.mondo, warn = U.fetch_optional("MONDO", U.MONDO_URLS, PHENO_DATA_DIR / "mondo" / "mondo.obo",
                                            max_age_days=0)
        warnings += [warn] if warn else []
    mondo = args.mondo or B.default_mondo()
    missing = [str(p) for p in (args.hp_obo, args.hpoa, args.genes_to_disease) if not p.is_file()]
    if missing:
        print("找不到：" + "、".join(missing), file=sys.stderr)
        return 2
    try:
        snap, more, manifest = U.build_snapshot(
            args.hp_obo, args.hpoa, args.genes_to_disease, kb_root=args.kb_root,
            orphanet=args.orphanet if args.orphanet.is_file() else None, mondo=mondo, allow_missing=allow,
            omim=B.platform_omim_rows() if "omim" not in allow else ([], {}))
    except B.KBBuildError as exc:
        print(f"建置中止：{exc}", file=sys.stderr)
        return 1
    if args.activate:
        H.activate(snap, args.kb_root)
        if U.install_obo(args.hp_obo):
            print("已更新 phenotype_data/hp.obo（舊檔另存 hp.obo.bak-時間）")
        U.write_picker_data(PHENO_DATA_DIR / "hp.obo", snap / H.GENE_TABLE, PHENO_DATA_DIR / U.PICKER_DATA)
    print("\n" + "=" * 60)
    print(U.summary_text(snap, manifest, warnings + more))
    print("=" * 60)
    print("已啟用（服務約一分鐘內自動改用）" if args.activate else "尚未啟用；要啟用請加 --activate 再跑一次，或在首頁按「更新 HPO」")
    return 0


if __name__ == "__main__":
    sys.exit(main())
