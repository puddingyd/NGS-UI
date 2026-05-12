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

add_callers_tag.py
==================
目的：
    將二級分析產生的 ensemble.fixed.vcf.gz（含兩個 sample column：
    {SAMPLE_ID}_DV 和 {SAMPLE_ID}_HC）處理為三級分析可用的單一樣本 VCF，
    並在 INFO 欄位新增 CALLERS tag，記錄每個 variant 是由哪些 caller 偵測到的。

輸出的 CALLERS 值：
    DV+HC  → DeepVariant 和 HaplotypeCaller 都有 call（最高信心度）
    DV     → 只有 DeepVariant call
    HC     → 只有 HaplotypeCaller call

使用方式（由 Nextflow module prepare_vcf.nf 呼叫）：
    bcftools view ensemble.fixed.vcf.gz | \\
    python3 add_callers_tag.py --sample NA12878_WES | \\
    bgzip -c > snv_for_annotation.vcf.gz
    tabix -p vcf snv_for_annotation.vcf.gz

    或直接指定輸入輸出檔：
    python3 add_callers_tag.py \\
        --input  NA12878_WES.ensemble.fixed.vcf.gz \\
        --sample NA12878_WES \\
        --output snv_for_annotation.vcf.gz

依賴：
    pip install cyvcf2   （比 pysam 更快，處理大型 VCF 效率更好）
    bgzip / tabix        （htslib，通常已在環境中）
"""

#!/usr/bin/env python3
"""
add_callers_tag.py
==================
目的：
    將二級分析產生的 ensemble.fixed.vcf.gz（含兩個 sample column：
    {SAMPLE_ID}_DV 和 {SAMPLE_ID}_HC）處理為三級分析可用的單一樣本 VCF。

    在 INFO 欄位新增以下 tag：
        CALLERS   → 哪些 caller 偵測到此 variant（DV+HC / DV / HC）
        DP_DV     → DeepVariant 的 read depth
        AD_DV     → DeepVariant 的 allelic depth（REF,ALT 逗號分隔）
        VAF_DV    → DeepVariant 的 variant allele fraction
        DP_HC     → HaplotypeCaller 的 read depth
        AD_HC     → HaplotypeCaller 的 allelic depth（REF,ALT 逗號分隔）

    為何保留兩個 caller 的 DP/AD：
        - 臨床審閱時可以比較兩個 caller 的支持度
        - DV 和 HC 的 DP 定義略有不同（DV 計算方式更保守）
        - AD 可以用來計算 VAF，確認 zygosity

使用方式（由 Nextflow module prepare_vcf.nf 呼叫）：
    bcftools view ensemble.fixed.vcf.gz | \\
    python3 add_callers_tag.py --sample NA12878_WES | \\
    bgzip -c > snv_for_annotation.vcf.gz

    或直接指定輸入輸出：
    python3 add_callers_tag.py \\
        --input  NA12878_WES.ensemble.fixed.vcf.gz \\
        --sample NA12878_WES \\
        --output snv_for_annotation.vcf.gz

依賴：
    pip install cyvcf2
"""

import argparse
import math
import sys

try:
    from cyvcf2 import VCF, Writer
except ImportError:
    print("[ERROR] 請先安裝 cyvcf2：pip install cyvcf2", file=sys.stderr)
    sys.exit(1)


# ──────────────────────────────────────────────
# GT 判斷
# ──────────────────────────────────────────────

def is_called(gt_tuple) -> bool:
    """
    判斷一個 sample 的 GT 是否為有效的 variant call。

    cyvcf2 genotypes 格式：[allele1, allele2, phased]
        allele = -1 → missing（.）
        allele =  0 → REF
        allele >= 1 → ALT
    """
    if gt_tuple is None:
        return False
    a1, a2 = gt_tuple[0], gt_tuple[1]
    if a1 == -1 or a2 == -1:
        return False
    if a1 == 0 and a2 == 0:
        return False
    return True


def determine_callers(variant, dv_idx: int, hc_idx: int) -> str:
    """DV 和 HC 的 call 狀態 → CALLERS 字串"""
    dv_called = is_called(variant.genotypes[dv_idx])
    hc_called = is_called(variant.genotypes[hc_idx])

    if dv_called and hc_called:
        return "DV+HC"
    elif dv_called:
        return "DV"
    else:
        return "HC"


# ──────────────────────────────────────────────
# FORMAT 欄位擷取輔助函式
# ──────────────────────────────────────────────

def get_dp(variant, sample_idx: int) -> str:
    """
    從 FORMAT/DP 取得 read depth。
    回傳字串，missing 時回傳 "."

    cyvcf2 integer missing value 為 -2147483648（INT_MIN）。
    """
    try:
        dp = variant.format("DP")
        if dp is None:
            return "."
        val = dp[sample_idx][0]
        if val < 0:
            return "."
        return str(val)
    except Exception:
        return "."


def get_ad(variant, sample_idx: int) -> str:
    """
    從 FORMAT/AD 取得 allelic depth。
    回傳 "REF,ALT" 格式字串，missing 時回傳 "."

    multiallelic site 的 AD 格式為 "REF,ALT1,ALT2"，完整保留。
    負數值（cyvcf2 的 missing 表示）替換為 0。
    """
    try:
        ad = variant.format("AD")
        if ad is None:
            return "."
        vals = ad[sample_idx]
        if all(v < 0 for v in vals):
            return "."
        cleaned = [str(max(int(v), 0)) for v in vals]
        return ",".join(cleaned)
    except Exception:
        return "."


def get_vaf(variant, sample_idx: int) -> str:
    """
    從 FORMAT/VAF 取得 variant allele fraction（DV 特有欄位）。
    HC 沒有 VAF 欄位，會在 except 中回傳 "."。
    四捨五入到小數點後 4 位。
    """
    try:
        vaf = variant.format("VAF")
        if vaf is None:
            return "."
        val = vaf[sample_idx][0]
        if math.isnan(val) or val < 0:
            return "."
        return f"{val:.4f}"
    except Exception:
        return "."


# ──────────────────────────────────────────────
# 主要處理函式
# ──────────────────────────────────────────────

def add_callers_tag(input_path: str, sample_id: str, output_path: str):
    """
    主要處理流程：
    1. 開啟輸入 VCF，確認 DV/HC sample column
    2. 在 header 新增 CALLERS、DP_DV、AD_DV、VAF_DV、DP_HC、AD_HC 定義
    3. 逐 variant 擷取值，寫出輸出 VCF
    """

    vcf_in = VCF(input_path)

    # ── 確認 sample column ────────────────────
    expected_dv = f"{sample_id}_DV"
    expected_hc = f"{sample_id}_HC"
    samples = vcf_in.samples

    print(f"[INFO] VCF 中的 sample columns：{samples}", file=sys.stderr)

    if expected_dv not in samples:
        print(f"[ERROR] 找不到 DV sample column：{expected_dv}", file=sys.stderr)
        sys.exit(1)
    if expected_hc not in samples:
        print(f"[ERROR] 找不到 HC sample column：{expected_hc}", file=sys.stderr)
        sys.exit(1)

    dv_idx = samples.index(expected_dv)
    hc_idx = samples.index(expected_hc)
    print(f"[INFO] DV index：{dv_idx}，HC index：{hc_idx}", file=sys.stderr)

    # ── 新增 INFO tag 定義到 header ───────────
    new_info_fields = [
        {
            'ID': 'CALLERS',
            'Number': '1',
            'Type': 'String',
            'Description': (
                'Variant callers that detected this variant: '
                'DV+HC (both), DV (DeepVariant only), HC (HaplotypeCaller only)'
            )
        },
        {
            'ID': 'DP_DV',
            'Number': '1',
            'Type': 'String',
            'Description': 'Read depth from DeepVariant (FORMAT/DP). Dot if missing.'
        },
        {
            'ID': 'AD_DV',
            'Number': '1',
            'Type': 'String',
            'Description': (
                'Allelic depths from DeepVariant (FORMAT/AD), REF,ALT comma-separated. '
                'Dot if missing.'
            )
        },
        {
            'ID': 'VAF_DV',
            'Number': '1',
            'Type': 'String',
            'Description': (
                'Variant allele fraction from DeepVariant (FORMAT/VAF). '
                'DV-specific field. Dot if missing or HC-only variant.'
            )
        },
        {
            'ID': 'DP_HC',
            'Number': '1',
            'Type': 'String',
            'Description': 'Read depth from HaplotypeCaller (FORMAT/DP). Dot if missing.'
        },
        {
            'ID': 'AD_HC',
            'Number': '1',
            'Type': 'String',
            'Description': (
                'Allelic depths from HaplotypeCaller (FORMAT/AD), REF,ALT comma-separated. '
                'Dot if missing or DV-only variant.'
            )
        },
    ]

    for field in new_info_fields:
        vcf_in.add_info_to_header(field)

    # ── 開啟輸出 ─────────────────────────────
    vcf_out = Writer(output_path, vcf_in, mode="w")

    # ── 逐 variant 處理 ──────────────────────
    n_total = 0
    n_dv_hc = 0
    n_dv_only = 0
    n_hc_only = 0

    for variant in vcf_in:
        n_total += 1

        # CALLERS tag
        callers = determine_callers(variant, dv_idx, hc_idx)
        variant.INFO["CALLERS"] = callers

        # DV 的 DP、AD、VAF
        variant.INFO["DP_DV"]  = get_dp(variant, dv_idx)
        variant.INFO["AD_DV"]  = get_ad(variant, dv_idx)
        variant.INFO["VAF_DV"] = get_vaf(variant, dv_idx)

        # HC 的 DP、AD（HC 無 VAF 欄位）
        variant.INFO["DP_HC"] = get_dp(variant, hc_idx)
        variant.INFO["AD_HC"] = get_ad(variant, hc_idx)

        # 統計
        if callers == "DV+HC":
            n_dv_hc += 1
        elif callers == "DV":
            n_dv_only += 1
        else:
            n_hc_only += 1

        vcf_out.write_record(variant)

    vcf_out.close()
    vcf_in.close()

    # ── 統計摘要 ─────────────────────────────
    print(f"[INFO] 處理完成", file=sys.stderr)
    print(f"[INFO]   總 variant 數：{n_total:,}", file=sys.stderr)
    print(f"[INFO]   DV+HC：{n_dv_hc:,} ({n_dv_hc/n_total*100:.1f}%)", file=sys.stderr)
    print(f"[INFO]   DV only：{n_dv_only:,} ({n_dv_only/n_total*100:.1f}%)", file=sys.stderr)
    print(f"[INFO]   HC only：{n_hc_only:,} ({n_hc_only/n_total*100:.1f}%)", file=sys.stderr)


# ──────────────────────────────────────────────
# 命令列介面
# ──────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "在 ensemble VCF 的 INFO 欄位新增 CALLERS tag 及 DP/AD/VAF 欄位"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用範例：
  # stdin → stdout（搭配 bgzip）
  bcftools view NA12878_WES.ensemble.fixed.vcf.gz | \\
  python3 add_callers_tag.py --sample NA12878_WES | \\
  bgzip -c > snv_for_annotation.vcf.gz && \\
  tabix -p vcf snv_for_annotation.vcf.gz

  # 直接指定輸入輸出
  python3 add_callers_tag.py \\
      --input  NA12878_WES.ensemble.fixed.vcf.gz \\
      --sample NA12878_WES \\
      --output snv_for_annotation.vcf
        """
    )
    parser.add_argument("--input",  "-i", default="-",
                        help="輸入 VCF 路徑（預設：stdin）")
    parser.add_argument("--sample", "-s", required=True,
                        help="樣本 ID（例如 NA12878_WES）")
    parser.add_argument("--output", "-o", default="-",
                        help="輸出 VCF 路徑（預設：stdout）")
    return parser.parse_args()


def main():
    args = parse_args()
    print(f"[INFO] 輸入：{args.input}", file=sys.stderr)
    print(f"[INFO] 樣本 ID：{args.sample}", file=sys.stderr)
    print(f"[INFO] 輸出：{args.output}", file=sys.stderr)
    add_callers_tag(args.input, args.sample, args.output)


if __name__ == "__main__":
    main()
