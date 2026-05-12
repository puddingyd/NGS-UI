/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - Parse CSQ Module
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
 * modules/parse_csq.nf
 * ====================
 * 目的：
 *   執行 parse_vep_csq.py，將 VEP annotation VCF 和 Pangolin VCF
 *   解析為結構化 TSV（snv_indel.annotated.tsv），供後續
 *   acmg_classifier.py 和 GUI 使用。
 *
 * 輸入：
 *   vep_ch:      tuple val(sample_id), path(vep_vcf), path(vep_tbi)
 *   pangolin_ch: tuple val(sample_id), path(pangolin_vcf), path(pangolin_tbi)
 *
 * 輸出：
 *   tsv_ch: tuple val(sample_id), path("*.snv_indel.annotated.tsv")
 *
 * 使用的容器：
 *   tertiary_python_1.0.0.sif（只需要 Python 3 標準函式庫）
 */

process PARSE_CSQ {
    label 'process_medium'

    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    // publishDir 直接輸出到樣本根目錄（計畫書 Section 15 的輸出規格）
    publishDir "${params.out_dir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(vep_vcf), path(vep_tbi)
    tuple val(sample_id2), path(pangolin_vcf), path(pangolin_tbi)

    output:
    tuple val(sample_id),
          path("${sample_id}.snv_indel.annotated.tsv"),
          emit: tsv_ch

    script:
    """
    python3 ${params.scripts_dir}/parse_vep_csq.py \\
        --vep_vcf      ${vep_vcf} \\
        --pangolin_vcf ${pangolin_vcf} \\
        --sample_id    ${sample_id} \\
        --output       ${sample_id}.snv_indel.annotated.tsv

    # 輸出統計
    echo "[PARSE_CSQ] ${sample_id} 完成" >&2
    TOTAL=\$(wc -l < ${sample_id}.snv_indel.annotated.tsv)
    echo "[PARSE_CSQ] 輸出 \$(( TOTAL - 1 )) variants" >&2

    # 確認 Pangolin 有命中
    PANG_HITS=\$(awk -F'\\t' 'NR>1 && \$51 != "."' \\
        ${sample_id}.snv_indel.annotated.tsv | wc -l)
    echo "[PARSE_CSQ] Pangolin 命中：\${PANG_HITS} variants" >&2
    """
}

// ──────────────────────────────────────────────────────────────
// 組合 workflow（供 main_tertiary.nf 呼叫）
// ──────────────────────────────────────────────────────────────

workflow PARSE_VEP_CSQ {
    take:
    vep_ch      // tuple val(sample_id), path(vep_vcf), path(vep_tbi)
    pangolin_ch // tuple val(sample_id), path(pangolin_vcf), path(pangolin_tbi)

    main:
    PARSE_CSQ(vep_ch, pangolin_ch)

    emit:
    tsv_ch = PARSE_CSQ.out.tsv_ch
}
