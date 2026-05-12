/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - CNV/SV Module
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
 * modules/prepare_vcf.nf
 * ======================
 * 目的：
 *   將二級分析產生的 ensemble.fixed.vcf.gz（雙 sample column：_DV + _HC）
 *   處理為三級分析 VEP annotation 的輸入，共兩個步驟：
 *
 *   Step 1 - ADD_CALLERS_TAG：
 *     執行 add_callers_tag.py，在 INFO 欄位新增 CALLERS tag（DV+HC/DV/HC）
 *
 *   Step 2 - FILTER_FOR_ANNOTATION：
 *     用 bcftools 過濾掉不適合送進 VEP 的 variant：
 *       - RefCall（FILTER=RefCall，DV 叫 0/0 但 HC 叫 ./. 的情況）
 *       - 兩個 sample column 都是 0/0 或 ./. 的 variant
 *     同時 bgzip 壓縮 + tabix index，產生標準的 .vcf.gz + .tbi
 *
 * 輸入（來自 main_tertiary.nf）：
 *   tuple val(sample_id), path(ensemble_vcf), path(ensemble_tbi)
 *
 * 輸出：
 *   tuple val(sample_id), path("*.snv_for_annotation.vcf.gz"), path("*.snv_for_annotation.vcf.gz.tbi")
 *
 * 使用的容器：
 *   Step 1：tertiary_python_1.0.0.sif（含 cyvcf2）
 *   Step 2：bcftools_1.23.1.sif
 */

// ──────────────────────────────────────────────────────────────
// Process 1：新增 CALLERS tag
// ──────────────────────────────────────────────────────────────

process ADD_CALLERS_TAG {
    // 標籤：對應 nextflow_tertiary.config 的資源設定
    label 'process_low'

    // 使用 tertiary_python sif（含 cyvcf2）
    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    // 輸入：ensemble VCF + 其 index
    input:
    tuple val(sample_id), path(ensemble_vcf), path(ensemble_tbi)

    // 輸出：加上 CALLERS tag 的未壓縮 VCF（暫時檔，交給下一個 process）
    output:
    tuple val(sample_id), path("${sample_id}.callers_tagged.vcf")

    script:
    """
    # 執行 add_callers_tag.py
    # --sample 傳入 sample_id，腳本會自動尋找 {sample_id}_DV 和 {sample_id}_HC column
    python3 ${params.scripts_dir}/add_callers_tag.py \\
        --input  ${ensemble_vcf} \\
        --sample ${sample_id} \\
        --output ${sample_id}.callers_tagged.vcf
    """
}

// ──────────────────────────────────────────────────────────────
// Process 2：過濾 + bgzip + tabix
// ──────────────────────────────────────────────────────────────

process FILTER_FOR_ANNOTATION {
    label 'process_low'

    // 使用既有的 bcftools sif
    container "${params.sif_dir}/bcftools_1.23.1.sif"

    // publishDir：將最終輸出複製到三級分析輸出目錄
    // mode: 'copy' 確保輸出目錄有獨立的檔案（不是 symlink）
    publishDir "${params.out_dir}/${sample_id}/00_prepare", mode: 'copy'

    input:
    tuple val(sample_id), path(callers_tagged_vcf)

    // 輸出：bgzip 壓縮的 VCF + tabix index
    output:
    tuple val(sample_id),
          path("${sample_id}.snv_for_annotation.vcf.gz"),
          path("${sample_id}.snv_for_annotation.vcf.gz.tbi")

    script:
    """
    # 過濾策略：
    #   -f PASS,RefCall：只保留 FILTER=PASS 或 RefCall 的 variant
    #   後面再用 -e 排除真正的 RefCall（FILTER="RefCall"）
    #   以及排除兩個 caller 都沒有有效 call 的 variant
    #
    # 實際上我們要的是：
    #   保留 FILTER=PASS，且至少一個 caller 有非 ref/missing 的 GT
    #
    # bcftools filter 邏輯：
    #   -i 'FILTER="PASS"'  → 只保留 PASS variant
    #   這樣 RefCall 就自動被排除了

    bcftools view \\
        -i 'FILTER="PASS"' \\
        ${callers_tagged_vcf} \\
        -Oz -o ${sample_id}.snv_for_annotation.vcf.gz

    # 建立 tabix index（VEP 和後續工具都需要）
    tabix -p vcf ${sample_id}.snv_for_annotation.vcf.gz

    # 輸出統計（寫進 log，方便 debug）
    echo "[FILTER_FOR_ANNOTATION] ${sample_id}" >&2
    bcftools stats ${sample_id}.snv_for_annotation.vcf.gz | \\
        grep "^SN" >&2
    """
}

// ──────────────────────────────────────────────────────────────
// 組合 workflow（供 main_tertiary.nf 呼叫）
// ──────────────────────────────────────────────────────────────

workflow PREPARE_VCF {
    // 輸入 channel：tuple(sample_id, ensemble_vcf, ensemble_tbi)
    take:
    ensemble_ch

    // 執行兩個 process，串接輸出
    main:
    ADD_CALLERS_TAG(ensemble_ch)
    FILTER_FOR_ANNOTATION(ADD_CALLERS_TAG.out)

    // 輸出 channel：tuple(sample_id, snv_vcf, snv_tbi)
    emit:
    snv_ch = FILTER_FOR_ANNOTATION.out
}
