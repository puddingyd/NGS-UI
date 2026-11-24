"""Versions of the reference databases the platform uses (home page panel).

Two groups:

* fixed with the tertiary Nextflow pipeline (VEP, gnomAD, dbNSFP, …) — these
  only change when the pipeline is upgraded, so they are listed as text;
* platform databases on this server — read from their own metadata, manifest
  or file date, so the panel shows what is really installed.

Every lookup is best-effort: a missing or unreadable source shows "未安裝"
instead of failing the page.
"""
from __future__ import annotations

import json
import sqlite3
import time
from pathlib import Path

from .. import config

PIPELINE = [
    ("VEP", "115（Ensembl／GENCODE 115）", "轉錄本、HGVS、consequence"),
    ("gnomAD", "v4.1（exome／genome）＋ v4.1 joint", "族群頻率；joint 欄位僅參考"),
    ("dbNSFP", "4.9c／5.3a（Research-only）", "AlphaMissense、BayesDel、ESM1b、REVEL、CADD 等預測分數"),
    ("ClinVar（報告版）", "2026-07-20", "報告與 tier 1A 固定使用此版"),
    ("1000 Genomes", "phase 3（popfreq_max_20150413）", "東亞 AF"),
    ("dbSNP", "v150", "rsID"),
    ("LOFTEE", "VEP plugin", "LoF 信心度"),
    ("Pangolin", "院內 GPU 版", "剪接預測（跨個案快取）"),
    ("P-KNN", "dbNSFP 5.3a 附帶", "missense 校準 LLR"),
    ("ClinGen ERepo", "隨 pipeline", "專家小組判讀對照"),
    ("gnomAD mito", "v3.1", "粒線體頻率"),
    ("STRchive", "隨 pipeline", "STR 致病門檻"),
    ("AnnotSV", "3.5.10", "CNV／SV 註解"),
    ("PharmCAT", "3.2.0（data 2026-02-09）", "藥物基因"),
]


def _date(ts: float | None) -> str:
    return time.strftime("%Y-%m-%d", time.localtime(ts)) if ts else ""


def _mtime(path) -> str:
    try:
        return _date(Path(path).stat().st_mtime)
    except (OSError, TypeError):
        return ""


def _json(path) -> dict:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return {}
    return value if isinstance(value, dict) else {}


def _sqlite_meta(path, table: str = "meta") -> dict:
    path = Path(path)
    if not path.is_file():
        return {}
    try:
        uri = f"file:{path.resolve()}?mode=ro"
        with sqlite3.connect(uri, uri=True, timeout=2) as conn:
            return {str(k): str(v) for k, v in conn.execute(f"SELECT key, value FROM {table}")}
    except sqlite3.Error:
        return {}


def _item(name: str, version: str, note: str = "", update: str = "", key: str = "") -> dict:
    return {"key": key or name, "name": name, "version": version or "未安裝", "note": note,
            "update": update, "installed": bool(version)}


def _hpo_items() -> list[dict]:
    from . import hpo_gene_table as H
    table = H.active_table()
    manifest = {}
    if table.term_status is not None:
        manifest = _json(table.gene_table.parent / "manifest.json")
    items = []
    if manifest:
        hp = manifest.get("hp_version", "")
        hpoa = manifest.get("hpoa_version", "")
        version = hp if hp == hpoa else f"hp.obo {hp}／註記 {hpoa}"
        items.append(_item("HPO", version, "清理後 HPO → 基因查表（表型分數）", "按「更新 HPO」", key="hpo"))
        kb = manifest.get("kb") or {}
        sources = kb.get("sources") or {}
        orpha = sources.get("orphanet") or {}
        if orpha:
            label = orpha.get("version") or orpha.get("date") or ""
            if orpha.get("fallback"):
                label = f"經 HPO genes_to_disease 轉載（{manifest.get('hpoa_version', '')}）"
            items.append(_item("Orphanet", label, "疾病–基因對應（表型知識庫）", "隨「更新 HPO」", key="orphanet"))
        mondo = sources.get("mondo") or {}
        if mondo:
            items.append(_item("MONDO", mondo.get("version") or "", "OMIM／Orphanet 疾病分組", "隨「更新 HPO」",
                               key="mondo"))
    else:
        obo = config.PHENO_DATA_DIR / "hp.obo"
        version = ""
        try:
            with open(obo, encoding="utf-8") as handle:
                for _ in range(10):
                    line = handle.readline()
                    if line.startswith("data-version:"):
                        version = line.split(":", 1)[1].strip().split("/")[-1]
                        break
        except OSError:
            pass
        items.append(_item("HPO", f"{version}（未清理查表）" if version else "", "phenotype_to_genes.txt",
                           "按「更新 HPO」", key="hpo"))
    return items


def _omim_item() -> dict:
    try:
        from . import omim_db, omim_store
        if omim_store.source() == "db":
            status = omim_db.status()
            counts = status.get("counts") or {}
            return _item("OMIM", f"平台資料庫（{counts.get('genes', 0):,} 個基因）",
                         f"資料庫更新 {_mtime(status.get('path'))}", "OMIM 擷取頁採用／每月更新清單", key="omim")
        return _item("OMIM", f"OMIM.xlsx（{_mtime(config.OMIM_XLSX)}）" if Path(config.OMIM_XLSX).is_file() else "",
                     "疾病說明與遺傳模式", "OMIM 擷取頁", key="omim")
    except Exception:  # noqa: BLE001
        return _item("OMIM", "", "疾病說明與遺傳模式", key="omim")


def _clinvar_latest_item() -> dict:
    meta = _sqlite_meta(config.CLINVAR_LATEST_DB)
    version = meta.get("release_date") or meta.get("source_date") or meta.get("version") or ""
    if not version and Path(config.CLINVAR_LATEST_DB).is_file():
        version = _mtime(config.CLINVAR_LATEST_DB)
    return _item("ClinVar（每週比對）", version, "卡片上的升／降 P/LP 箭頭與 In-house ACMG 參考", "每週三 04:00 自動",
                 key="clinvar_latest")


def _litvar2_item() -> dict:
    try:
        from . import litvar2_store
        meta = litvar2_store.database_metadata()
        version = str(meta.get("dataset_date") or "") if meta.get("available", True) else ""
    except Exception:  # noqa: BLE001
        version = ""
    return _item("LitVar2", version, "文獻（PMID）", "每月 1 日自動", key="litvar2")


def _acmg_item() -> dict:
    rulebook = _json(Path(__file__).parent / "acmg_engine" / "rulebook_v1.json")
    refs = _sqlite_meta(Path(config.ACMG_REFS_DIR) / "acmg_refs.sqlite")
    note = []
    if refs.get("clinvar_release"):
        note.append(f"參考表 ClinVar {refs['clinvar_release']}")
    if refs.get("built_at"):
        note.append(f"建置 {refs['built_at'][:10]}")
    return _item("In-house ACMG 規則", rulebook.get("version", ""), "；".join(note) or "平台自建 ACMG 規則",
                 "改規則後重算", key="acmg")


def _gene_disease_item() -> dict:
    manifest = _json(Path(config.GENE_DISEASE_DB).parent / "source_manifest.json")
    if not manifest:
        return _item("GenCC／ClinGen 基因–疾病", _mtime(config.GENE_DISEASE_DB), "卡片疾病清單補充來源", key="gencc")
    dates = []
    for source in manifest.get("source_files") or []:
        if isinstance(source, dict) and source.get("downloaded_at"):
            dates.append(str(source["downloaded_at"])[:10])
    version = max(dates) if dates else _mtime(Path(config.GENE_DISEASE_DB))
    return _item("GenCC／ClinGen 基因–疾病", f"下載 {version}" if version else "", "卡片疾病清單補充來源",
                 key="gencc")


def _gpn_items() -> list[dict]:
    items = []
    msa = Path(config.GPN_MSA_DB)
    items.append(_item("GPN-MSA", _mtime(msa) and f"安裝 {_mtime(msa)}", "保守性預測分數", key="gpn_msa"))
    revs = []
    for name, label in (("manifest.v.json", "V"), ("manifest.json", "M")):
        manifest = _json(Path(config.GPN_STAR_DIR) / name)
        if manifest.get("revision"):
            revs.append(f"{label} {str(manifest['revision'])[:8]}")
    items.append(_item("GPN-Star", "、".join(revs), "Research-only 預測分數", key="gpn_star"))
    return items


def _file_item(name: str, path, note: str, update: str = "", key: str = "", prefix: str = "更新") -> dict:
    date = _mtime(path)
    return _item(name, f"{prefix} {date}" if date else "", note, update, key=key)


def collect() -> dict:
    platform = []
    for getter in (_hpo_items,):
        try:
            platform += getter()
        except Exception:  # noqa: BLE001
            platform.append(_item("HPO", "", "表型分數", key="hpo"))
    platform += [
        _omim_item(),
        _clinvar_latest_item(),
        _file_item("GeneBe", config.GENEBE_DB, "第二意見 ACMG", key="genebe"),
        _acmg_item(),
        _litvar2_item(),
        _file_item("院內 AF", config.INHOUSE_AF_DB, "院內 WGS 族群頻率", key="inhouse_af"),
        _gene_disease_item(),
    ]
    mane = Path(config.BIOTOOLS_DIR) / "MANE.GRCh38.v1.5.summary.txt.gz"
    platform.append(_item("MANE", "v1.5" if mane.is_file() else "", "RefSeq 對應", key="mane"))
    platform.append(_file_item("GIAB 困難區", Path(config.GIAB_STRAT_DIR) / "strata_manifest.json",
                               "GIAB v3.6 stratification", key="giab", prefix="安裝"))
    platform += _gpn_items()
    spliceai = Path(config.BIOTOOLS_DIR) / "spliceai" / "spliceai_scores.raw.snv.hg38.vcf.gz"
    platform.append(_item("SpliceAI", "1.3 precomputed" if spliceai.is_file() else "", "Research-only 剪接預測",
                          key="spliceai"))
    return {
        "pipeline": [{"key": n, "name": n, "version": v, "note": note, "update": "隨三級 pipeline", "installed": True}
                     for n, v, note in PIPELINE],
        "platform": platform,
    }
