#!/usr/bin/env python3
"""
 * =========================================================
 * WGS/WES Germline Analysis Pipeline
 * =========================================================
 * Author   : Po-Yu Lin (林伯昱)
 * Institute: Department of Neurology and
 *            Department of Genomic Medicine,
 *            National Cheng Kung University Hospital
 * Contact  : p88124019@gs.ncku.edu.tw
 *
 * Copyright (c) 2026, Po-Yu Lin
 * Licensed under the MIT License
 *
 * This pipeline was developed for clinical germline variant
 * analysis. Please cite appropriately if used in research.
 *
 * DISCLAIMER: This pipeline is provided "as is" without
 * warranty of any kind. The authors and their institution
 * make no representations or warranties regarding the
 * accuracy, completeness, or suitability of the analysis
 * results for any clinical or research purpose. Users are
 * solely responsible for validating and interpreting all
 * results. This software shall not be held liable for any
 * direct, indirect, or consequential damages arising from
 * its use.
 * =========================================================
parse_vep_csq.py
================
將 VEP annotation VCF + Pangolin VCF 解析為結構化 TSV，
供後續 acmg_classifier.py 和 GUI 使用。

輸入：
  --vep_vcf      VEP annotation VCF（*.vep.vcf.gz）
  --pangolin_vcf Pangolin splice score VCF（*.pangolin.vcf.gz）
  --sample_id    樣本 ID（用於 GT 欄位解析）
  --output       輸出 TSV 路徑

輸出欄位（snv_indel.annotated.tsv）：
  CHROM, POS, REF, ALT,
  GENE, TRANSCRIPT, TRANSCRIPT_TYPE,
  HGVS_C, HGVS_P, CONSEQUENCE, IMPACT,
  EXON, INTRON,
  MANE_ALL (JSON),
  CALLERS, DP_DV, AD_DV, VAF_DV, DP_HC, AD_HC,
  ZYGOSITY, GT_DV, GT_HC,
  GNOMAD_G_AF, GNOMAD_G_EAS_AF, GNOMAD_E_AF, GNOMAD_E_EAS_AF,
  GNOMAD_E_AF_DBNSFP, GNOMAD_E_EAS_AF_DBNSFP,
  CLINVAR_SIG, CLINVAR_STARS, CLINVAR_DN, CLINVAR_SIGCONF,
  LOFTEE, LOFTEE_FILTER, LOFTEE_FLAGS,
  LOFTOOL,
  BAYESDEL_NOAF, BAYESDEL_NOAF_PRED,
  ALPHAMISSENSE, ALPHAMISSENSE_PRED,
  ESM1B, ESM1B_PRED,
  VARITY_R,
  SIFT, SIFT_PRED,
  DANN,
  PHACTBOOST,
  PHYLOP100,
  GERP,
  PANGOLIN_SCORE, PANGOLIN_DETAIL,
  DOMAINS, SWISSPROT

設計原則：
  - CSQ 欄位順序從 VCF header 動態解析，不硬編碼位置
  - PICK=1 的 transcript 作為主 row
  - MANE_ALL 收集所有 MANE Select + MANE Plus Clinical transcript
  - ClinVar stars 從 CLNREVSTAT 轉換（規則見 clnrevstat_to_stars()）
  - Pangolin 分數從獨立 VCF 對位回填（以 CHROM+POS+REF+ALT 為 key）
  - Zygosity 從 GT_DV 和 GT_HC 推導（het/hom/hemizygous）

作者：Po-Yu Lin（林伯昱）
機構：國立成功大學醫院基因醫學部
"""

import argparse
import gzip
import json
import re
import sys
from collections import defaultdict


# ──────────────────────────────────────────────────────────────
# ClinVar review status → stars 轉換
# 來源：https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
# ──────────────────────────────────────────────────────────────

CLNREVSTAT_STARS = {
    "practice_guideline":                                    4,
    "reviewed_by_expert_panel":                              3,
    "criteria_provided_multiple_submitters_no_conflicts":    2,
    "criteria_provided_conflicting_classifications":         1,
    "criteria_provided_single_submitter":                    1,
    "no_assertion_criteria_provided":                        0,
    "no_classification_provided":                           0,
    "no_classification_for_the_single_variant":              0,
}


def clnrevstat_to_stars(revstat: str) -> int:
    """
    將 ClinVar CLNREVSTAT 字串轉換為 0-4 星。
    VEP custom annotation 中空白被替換為底線，& 為分隔符。
    取多個值中最高的星數（多個 submitter 情況）。
    """
    if not revstat or revstat == ".":
        return 0
    # VEP 把空白換成底線，把 & 換成分隔符
    # 格式：criteria_provided&_multiple_submitters&_no_conflicts
    # 需要把 &_ 合併成單一字串再查表
    normalized = revstat.replace("&_", "_").replace("&", "_").lower()
    return CLNREVSTAT_STARS.get(normalized, 0)


# ──────────────────────────────────────────────────────────────
# Zygosity 推導
# ──────────────────────────────────────────────────────────────

def infer_zygosity(gt_dv: str, gt_hc: str, chrom: str) -> str:
    """
    從 DV 和 HC 的 GT 推導 zygosity。
    優先使用 DV（DeepVariant 對 het 的準確度較高）。
    hemizygous：chrX/chrY 的 hom alt（男性）
    """
    gt = gt_dv if gt_dv not in (".", "./.", ".|.") else gt_hc
    if gt in (".", "./.", ".|."):
        return "unknown"

    # 標準化 GT（phased → unphased）
    gt_norm = gt.replace("|", "/")
    alleles = gt_norm.split("/")

    if len(alleles) != 2:
        return "unknown"

    ref_count = alleles.count("0")
    alt_alleles = [a for a in alleles if a not in ("0", ".")]

    if ref_count == 2:
        return "ref"
    elif ref_count == 0 and len(alt_alleles) == 2:
        # chrX/Y 單條染色體（血液學上的 hemizygous）
        if chrom in ("chrX", "chrY", "X", "Y"):
            return "hemizygous"
        return "hom"
    elif ref_count == 1:
        return "het"
    else:
        return "unknown"


# ──────────────────────────────────────────────────────────────
# Pangolin VCF 解析
# ──────────────────────────────────────────────────────────────

def load_pangolin_scores(pangolin_vcf: str) -> dict:
    """
    讀取 Pangolin VCF，建立 (chrom, pos, ref, alt) → (score, detail) 的 dict。

    Pangolin INFO 格式：
      Pangolin=GENE_ID|pos:score_change|pos:score_change|Warnings:

    score 取所有 pos:score 中絕對值最大的（最強的 splice 影響）。
    detail 保留完整字串供 GUI 顯示。
    """
    scores = {}
    opener = gzip.open if pangolin_vcf.endswith(".gz") else open

    with opener(pangolin_vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            info = parts[7]

            # 從 INFO 欄位取 Pangolin=... 字串
            pangolin_str = ""
            for field in info.split(";"):
                if field.startswith("Pangolin="):
                    pangolin_str = field[len("Pangolin="):]
                    break

            if not pangolin_str:
                continue

            # 解析分數：格式 gene|pos:score|pos:score|Warnings:
            detail = pangolin_str
            max_score = 0.0
            segments = pangolin_str.split("|")
            for seg in segments[1:]:  # 跳過第一個（gene ID）
                if seg.startswith("Warnings"):
                    break
                if ":" in seg:
                    try:
                        score_val = float(seg.split(":")[1])
                        if abs(score_val) > abs(max_score):
                            max_score = score_val
                    except (ValueError, IndexError):
                        pass

            key = (chrom, pos, ref, alt)
            scores[key] = (max_score, detail)

    return scores


# ──────────────────────────────────────────────────────────────
# VCF header 解析：取得 CSQ 欄位順序
# ──────────────────────────────────────────────────────────────

def parse_csq_fields(vcf_path: str) -> dict:
    """
    從 VCF header 的 ##INFO=<ID=CSQ,...,Format="..."> 行
    動態解析 CSQ 每個欄位的名稱和索引。

    回傳 {欄位名稱: 索引} 的 dict（索引從 0 開始）。
    """
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if not line.startswith("##"):
                break
            if line.startswith("##INFO=<ID=CSQ"):
                # 取 Format: 後面到 "> 之前的字串
                m = re.search(r'Format: ([^"]+)"', line)
                if m:
                    fields = m.group(1).rstrip(">").split("|")
                    return {name: idx for idx, name in enumerate(fields)}
    raise ValueError("找不到 CSQ FORMAT 定義，請確認 VEP VCF header 格式")


# ──────────────────────────────────────────────────────────────
# 從 CSQ transcript dict 取值（空字串 → "."）
# ──────────────────────────────────────────────────────────────

def get(tx: dict, field: str) -> str:
    """從 transcript dict 安全取值，空字串回傳 '.'"""
    val = tx.get(field, "")
    return val if val else "."


# ──────────────────────────────────────────────────────────────
# MANE_ALL JSON 建立
# ──────────────────────────────────────────────────────────────

def build_mane_all(transcripts: list, csq_fields: dict) -> str:
    """
    從所有 transcript 中收集 MANE Select 和 MANE Plus Clinical，
    建立 MANE_ALL JSON 字串供 GUI 展開顯示。

    設計原則：
      - MANE Select 和 MANE Plus Clinical 都是臨床相關 transcript，兩者都收集
      - 同一個 ENST 不重複收錄（seen_tx 去重）
      - 只有在兩者都無的情況下，才不收錄任何 transcript（CANONICAL 不進 MANE_ALL）

    格式：
    [
      {"tx": "NM_007294.4", "type": "MANE_SELECT",
       "consequence": "missense_variant", "hgvsc": "c.5266dupC", "hgvsp": "p.Gln1756fs"},
      {"tx": "NM_007297.4", "type": "MANE_PLUS_CLINICAL",
       "consequence": "missense_variant", "hgvsc": "c.5266dupC", "hgvsp": "p.Gln1756fs"}
    ]
    """
    mane_entries = []
    seen_tx = set()

    for tx in transcripts:
        mane_select = tx.get("MANE_SELECT", "")
        mane_plus   = tx.get("MANE_PLUS_CLINICAL", "")
        feature     = tx.get("Feature", "")

        if feature in seen_tx:
            continue

        # MANE Select 和 MANE Plus Clinical 各自獨立判斷（不用 elif）
        # 同一個 variant 可能同時有兩個不同 ENST 的 MANE transcript
        if mane_select:
            seen_tx.add(feature)
            mane_entries.append({
                "tx":          mane_select,
                "enst":        feature,
                "type":        "MANE_SELECT",
                "consequence": tx.get("Consequence", ""),
                "hgvsc":       tx.get("HGVSc", ""),
                "hgvsp":       tx.get("HGVSp", ""),
                "impact":      tx.get("IMPACT", ""),
            })
        if mane_plus and feature not in seen_tx:
            # 注意：mane_plus 的 if 不是 elif，但要確認同一個 ENST 不重複
            # （理論上 MANE_SELECT 和 MANE_PLUS_CLINICAL 不會在同一個 ENST 上同時標記）
            seen_tx.add(feature)
            mane_entries.append({
                "tx":          mane_plus,
                "enst":        feature,
                "type":        "MANE_PLUS_CLINICAL",
                "consequence": tx.get("Consequence", ""),
                "hgvsc":       tx.get("HGVSc", ""),
                "hgvsp":       tx.get("HGVSp", ""),
                "impact":      tx.get("IMPACT", ""),
            })

    return json.dumps(mane_entries, ensure_ascii=False) if mane_entries else "[]"


# ──────────────────────────────────────────────────────────────
# TRANSCRIPT_TYPE 判斷
# ──────────────────────────────────────────────────────────────

def get_transcript_type(tx: dict) -> str:
    """
    判斷代表 transcript 的類型，供 GUI 顯示 badge。
    優先順序：MANE_SELECT > MANE_PLUS_CLINICAL > CANONICAL > OTHER
    """
    if tx.get("MANE_SELECT"):
        return "MANE_SELECT"
    elif tx.get("MANE_PLUS_CLINICAL"):
        return "MANE_PLUS_CLINICAL"
    elif tx.get("CANONICAL") == "YES":
        return "CANONICAL"
    else:
        return "OTHER"


# ──────────────────────────────────────────────────────────────
# GT 解析（從 FORMAT + sample column）
# ──────────────────────────────────────────────────────────────

def parse_gt_field(format_str: str, sample_str: str, field: str) -> str:
    """
    從 FORMAT 和 sample 欄位取出指定 field 的值。
    例如：FORMAT=GT:DP:AD，sample=0/1:30:15,15，field=GT → "0/1"
    """
    if not format_str or not sample_str or sample_str == ".":
        return "."
    fields = format_str.split(":")
    values = sample_str.split(":")
    if field not in fields:
        return "."
    idx = fields.index(field)
    return values[idx] if idx < len(values) else "."


# ──────────────────────────────────────────────────────────────
# 主解析流程
# ──────────────────────────────────────────────────────────────

def parse_vep_vcf(vep_vcf: str, pangolin_scores: dict,
                  sample_id: str, output: str):
    """
    逐行讀取 VEP VCF，解析每個 variant 的代表 transcript，
    整合 Pangolin 分數，輸出 TSV。
    """

    # 動態解析 CSQ 欄位順序
    csq_fields = parse_csq_fields(vep_vcf)

    # TSV 輸出欄位定義（與計畫書 Section 15 對應）
    output_columns = [
        "CHROM", "POS", "REF", "ALT",
        "GENE", "TRANSCRIPT", "TRANSCRIPT_TYPE",
        "HGVS_C", "HGVS_P", "CONSEQUENCE", "IMPACT",
        "EXON", "INTRON",
        "MANE_ALL",
        "CALLERS", "DP_DV", "AD_DV", "VAF_DV", "DP_HC", "AD_HC",
        "ZYGOSITY", "GT_DV", "GT_HC",
        "GNOMAD_G_AF", "GNOMAD_G_EAS_AF",
        "GNOMAD_E_AF", "GNOMAD_E_EAS_AF",
        "GNOMAD_E_AF_DBNSFP", "GNOMAD_E_EAS_AF_DBNSFP",
        "CLINVAR_SIG", "CLINVAR_STARS", "CLINVAR_DN", "CLINVAR_SIGCONF",
        "LOFTEE", "LOFTEE_FILTER", "LOFTEE_FLAGS",
        "LOFTOOL",
        "BAYESDEL_NOAF", "BAYESDEL_NOAF_PRED",
        "ALPHAMISSENSE", "ALPHAMISSENSE_PRED",
        "ESM1B", "ESM1B_PRED",
        "VARITY_R",
        "SIFT", "SIFT_PRED",
        "DANN",
        "PHACTBOOST",
        "PHYLOP100",
        "GERP",
        "PANGOLIN_SCORE", "PANGOLIN_DETAIL",
        "DOMAINS", "SWISSPROT",
    ]

    opener = gzip.open if vep_vcf.endswith(".gz") else open
    sample_dv = f"{sample_id}_DV"
    sample_hc = f"{sample_id}_HC"

    # sample column index（從 #CHROM 行解析）
    col_dv = None
    col_hc = None

    written = 0
    skipped = 0

    with opener(vep_vcf, "rt") as fin, open(output, "w") as fout:
        # 寫入 header
        fout.write("\t".join(output_columns) + "\n")

        for line in fin:
            line = line.rstrip("\n")

            # 從 #CHROM 行取 sample column 位置
            if line.startswith("#CHROM"):
                cols = line.split("\t")
                if sample_dv in cols:
                    col_dv = cols.index(sample_dv)
                if sample_hc in cols:
                    col_hc = cols.index(sample_hc)
                continue

            if line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 8:
                continue

            chrom   = parts[0]
            pos     = parts[1]
            ref     = parts[3]
            alt     = parts[4]
            info    = parts[7]
            fmt     = parts[8] if len(parts) > 8 else ""
            smp_dv  = parts[col_dv] if col_dv and col_dv < len(parts) else "."
            smp_hc  = parts[col_hc] if col_hc and col_hc < len(parts) else "."

            # ── INFO 欄位解析 ──────────────────────────────
            info_dict = {}
            for field in info.split(";"):
                if "=" in field:
                    k, v = field.split("=", 1)
                    info_dict[k] = v

            # CALLERS tag（來自 prepare_vcf）
            callers  = info_dict.get("CALLERS", ".")
            dp_dv    = info_dict.get("DP_DV", ".")
            ad_dv    = info_dict.get("AD_DV", ".")
            vaf_dv   = info_dict.get("VAF_DV", ".")
            dp_hc    = info_dict.get("DP_HC", ".")
            ad_hc    = info_dict.get("AD_HC", ".")

            # GT（從 FORMAT/sample 欄位解析）
            gt_dv = parse_gt_field(fmt, smp_dv, "GT")
            gt_hc = parse_gt_field(fmt, smp_hc, "GT")
            zygosity = infer_zygosity(gt_dv, gt_hc, chrom)

            # ClinVar（custom annotation，直接在 INFO 欄位）
            clinvar_sig      = info_dict.get("ClinVar_CLNSIG", ".")
            clinvar_revstat  = info_dict.get("ClinVar_CLNREVSTAT", ".")
            clinvar_dn       = info_dict.get("ClinVar_CLNDN", ".")
            clinvar_sigconf  = info_dict.get("ClinVar_CLNSIGCONF", ".")
            clinvar_stars    = clnrevstat_to_stars(clinvar_revstat)

            # ── CSQ 解析 ──────────────────────────────────
            csq_raw = info_dict.get("CSQ", "")
            if not csq_raw:
                skipped += 1
                continue

            # 每個 transcript 是 CSQ 字串中以 , 分隔的一段
            # 注意：CSQ 內部的值可能含 & 分隔的多個後果，不含 ,
            transcripts_raw = csq_raw.split(",")
            transcripts = []
            for tx_raw in transcripts_raw:
                vals = tx_raw.split("|")
                # 補齊欄位數量（有些欄位可能是空的）
                while len(vals) < len(csq_fields):
                    vals.append("")
                tx_dict = {name: vals[idx] for name, idx in csq_fields.items()
                           if idx < len(vals)}
                transcripts.append(tx_dict)

            # 找 PICK=1 的代表 transcript
            picked_tx = None
            for tx in transcripts:
                if tx.get("PICK") == "1":
                    picked_tx = tx
                    break

            # 若沒有 PICK=1（理論上不應發生），取第一個
            if picked_tx is None:
                picked_tx = transcripts[0] if transcripts else {}
                skipped_pick = True
            else:
                skipped_pick = False

            # ── 從代表 transcript 提取欄位 ────────────────
            gene            = get(picked_tx, "SYMBOL")
            transcript      = get(picked_tx, "Feature")
            transcript_type = get_transcript_type(picked_tx)
            hgvs_c          = get(picked_tx, "HGVSc")
            hgvs_p          = get(picked_tx, "HGVSp")
            consequence     = get(picked_tx, "Consequence")
            impact          = get(picked_tx, "IMPACT")
            exon            = get(picked_tx, "EXON")
            intron          = get(picked_tx, "INTRON")

            # gnomAD（VEP cache 內建）
            gnomad_g_af     = get(picked_tx, "gnomADg_AF")
            gnomad_g_eas_af = get(picked_tx, "gnomADg_EAS_AF")
            gnomad_e_af     = get(picked_tx, "gnomADe_AF")
            gnomad_e_eas_af = get(picked_tx, "gnomADe_EAS_AF")

            # gnomAD（dbNSFP 版本，exome only）
            gnomad_e_af_db      = get(picked_tx, "gnomAD_exomes_AF")
            gnomad_e_eas_af_db  = get(picked_tx, "gnomAD_exomes_EAS_AF")

            # LOFTEE
            loftee        = get(picked_tx, "LoF")
            loftee_filter = get(picked_tx, "LoF_filter")
            loftee_flags  = get(picked_tx, "LoF_flags")
            loftool       = get(picked_tx, "LoFtool")

            # In silico scores
            bayesdel_noaf      = get(picked_tx, "BayesDel_noAF_score")
            bayesdel_noaf_pred = get(picked_tx, "BayesDel_noAF_pred")
            alphamissense      = get(picked_tx, "AlphaMissense_score")
            alphamissense_pred = get(picked_tx, "AlphaMissense_pred")
            esm1b              = get(picked_tx, "ESM1b_score")
            esm1b_pred         = get(picked_tx, "ESM1b_pred")
            varity_r           = get(picked_tx, "VARITY_R_score")
            sift               = get(picked_tx, "SIFT_score")
            sift_pred          = get(picked_tx, "SIFT_pred")
            dann               = get(picked_tx, "DANN_score")
            phactboost         = get(picked_tx, "PHACTboost_score")
            phylop100          = get(picked_tx, "phyloP100way_vertebrate")
            gerp               = get(picked_tx, "GERP++_RS")

            # Protein domains 和 UniProt
            domains   = get(picked_tx, "DOMAINS")
            swissprot = get(picked_tx, "SWISSPROT")

            # ── MANE_ALL JSON ─────────────────────────────
            mane_all = build_mane_all(transcripts, csq_fields)

            # ── Pangolin 分數回填 ─────────────────────────
            # key：(chrom, pos, ref, alt)
            # alt 可能是多 allelic（只取第一個）
            alt_first = alt.split(",")[0]
            pang_key = (chrom, pos, ref, alt_first)
            if pang_key in pangolin_scores:
                pang_score, pang_detail = pangolin_scores[pang_key]
                pangolin_score  = f"{pang_score:.4f}"
                pangolin_detail = pang_detail
            else:
                pangolin_score  = "."
                pangolin_detail = "."

            # ── 組裝輸出 row ──────────────────────────────
            row = [
                chrom, pos, ref, alt,
                gene, transcript, transcript_type,
                hgvs_c, hgvs_p, consequence, impact,
                exon, intron,
                mane_all,
                callers, dp_dv, ad_dv, vaf_dv, dp_hc, ad_hc,
                zygosity, gt_dv, gt_hc,
                gnomad_g_af, gnomad_g_eas_af,
                gnomad_e_af, gnomad_e_eas_af,
                gnomad_e_af_db, gnomad_e_eas_af_db,
                clinvar_sig, str(clinvar_stars), clinvar_dn, clinvar_sigconf,
                loftee, loftee_filter, loftee_flags,
                loftool,
                bayesdel_noaf, bayesdel_noaf_pred,
                alphamissense, alphamissense_pred,
                esm1b, esm1b_pred,
                varity_r,
                sift, sift_pred,
                dann,
                phactboost,
                phylop100,
                gerp,
                pangolin_score, pangolin_detail,
                domains, swissprot,
            ]

            fout.write("\t".join(str(v) for v in row) + "\n")
            written += 1

    print(f"[parse_vep_csq] 完成：寫入 {written} variants，跳過 {skipped} variants（無 CSQ）",
          file=sys.stderr)


# ──────────────────────────────────────────────────────────────
# 主程式
# ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="解析 VEP CSQ 欄位 + Pangolin 分數，輸出結構化 TSV"
    )
    parser.add_argument("--vep_vcf",      required=True,  help="VEP annotation VCF（*.vep.vcf.gz）")
    parser.add_argument("--pangolin_vcf", required=True,  help="Pangolin splice score VCF（*.pangolin.vcf.gz）")
    parser.add_argument("--sample_id",    required=True,  help="樣本 ID（用於 GT 欄位解析）")
    parser.add_argument("--output",       required=True,  help="輸出 TSV 路徑")
    args = parser.parse_args()

    print(f"[parse_vep_csq] 載入 Pangolin 分數：{args.pangolin_vcf}", file=sys.stderr)
    pangolin_scores = load_pangolin_scores(args.pangolin_vcf)
    print(f"[parse_vep_csq] Pangolin 分數載入完成：{len(pangolin_scores)} variants", file=sys.stderr)

    print(f"[parse_vep_csq] 解析 VEP VCF：{args.vep_vcf}", file=sys.stderr)
    parse_vep_vcf(args.vep_vcf, pangolin_scores, args.sample_id, args.output)


if __name__ == "__main__":
    main()
