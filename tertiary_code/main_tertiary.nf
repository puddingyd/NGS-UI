#!/usr/bin/env nextflow
/*
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
 * main_tertiary.nf
 * ================
 * 臨床 WGS/WES 三級分析主 workflow
 *
 * 使用方式：
 *   nextflow -c nextflow_tertiary.config run main_tertiary.nf \
 *       -profile local \
 *       --sample_id NA12878_WES \
 *       --input_dir /scratch/pylin1991/Pipeline_test/NA12878_WES_PON/NA12878_WES \
 *       --seq_type WES \
 *       --out_dir /scratch/pylin1991/tertiary_test \
 *       -resume
 *
 * 開發進度（Phase 1 進行中）：
 *   ✅ PREPARE_VCF   - ensemble VCF 前處理（CALLERS tag + 過濾）
 *   ✅ SNV_ANNOTATE  - VEP 115 annotation（Phase 1）
 *   ✅ PARSE_CSQ     - transcript 選取 + MANE_ALL JSON（Phase 1 下一步）
 *   🔲 ACMG_CLASSIFY - ACMG evidence 收集與分類（Phase 1）
 *   🔲 CNV_SV        - AnnotSV（Phase 2）
 *   🔲 MITO          - mtDNA annotation（Phase 3）
 *   🔲 STR           - STRchive rule-based（Phase 3）
 *   🔲 ROH           - consanguinity + UPD（Phase 2）
 *   🔲 PHENOTYPE     - Exomiser + LIRICAL（Phase 2）
 *   🔲 WHATS_HAP     - phasing（Phase 4，條件觸發）
 *   🔲 PGX           - Aldy（Phase 3）
 *   🔲 SECONDARY     - ACMG SF v3.2（Phase 3）
 */

nextflow.enable.dsl = 2

// ──────────────────────────────────────────────────────────────
// 匯入 modules
// ──────────────────────────────────────────────────────────────

include { PREPARE_VCF   } from './modules/prepare_vcf.nf'
include { SNV_ANNOTATE  } from './modules/snv_annotation.nf'
include { PARSE_VEP_CSQ } from './modules/parse_csq.nf'
// 後續 phase 逐步新增：
// include { ACMG_CLASSIFY } from './modules/acmg_classifier.nf'

// ──────────────────────────────────────────────────────────────
// 參數驗證
// ──────────────────────────────────────────────────────────────

def validate_params() {
    // 必填參數檢查
    if (!params.sample_id) {
        error "[ERROR] --sample_id 為必填參數，例如：--sample_id NA12878_WES"
    }
    if (!params.input_dir) {
        error "[ERROR] --input_dir 為必填參數（二級分析輸出目錄）"
    }
    if (!params.out_dir) {
        error "[ERROR] --out_dir 為必填參數（三級分析輸出目錄）"
    }

    // seq_type 只能是 WES 或 WGS
    if (!['WES', 'WGS'].contains(params.seq_type)) {
        error "[ERROR] --seq_type 必須是 WES 或 WGS，目前值：${params.seq_type}"
    }

    // 確認輸入目錄存在
    def input_dir = file(params.input_dir)
    if (!input_dir.exists()) {
        error "[ERROR] input_dir 不存在：${params.input_dir}"
    }

    // 確認必要的資料庫路徑存在（早期報錯，避免跑到一半才失敗）
    def vep_cache = file(params.vep_cache)
    if (!vep_cache.exists()) {
        error "[ERROR] VEP cache 目錄不存在：${params.vep_cache}"
    }

    def dbnsfp = file(params.dbnsfp)
    if (!dbnsfp.exists()) {
        error "[ERROR] dbNSFP 找不到：${params.dbnsfp}"
    }

    def loftee_dir = file(params.loftee_dir)
    if (!loftee_dir.exists()) {
        error "[ERROR] LOFTEE 資料目錄不存在：${params.loftee_dir}"
    }

    def clinvar = file(params.clinvar)
    if (!clinvar.exists()) {
        error "[ERROR] ClinVar VCF 找不到：${params.clinvar}"
    }
}

// ──────────────────────────────────────────────────────────────
// 主 workflow
// ──────────────────────────────────────────────────────────────

workflow {

    // 執行參數驗證
    validate_params()

    // 印出執行資訊
    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║         臨床三級分析 Pipeline  v1.0.0                ║
    ╚══════════════════════════════════════════════════════╝
    樣本 ID   : ${params.sample_id}
    定序類型  : ${params.seq_type}
    輸入目錄  : ${params.input_dir}
    輸出目錄  : ${params.out_dir}
    HPO 輸入  : ${params.hpo ?: '（未提供）'}
    Run Evo2  : ${params.run_evo2}
    Run Phase : ${params.run_phasing}
    容器目錄  : ${params.sif_dir}
    VEP cache : ${params.vep_cache}
    dbNSFP    : ${params.dbnsfp}
    ClinVar   : ${params.clinvar}
    """.stripIndent()

    // ── 建立輸入 channel ──────────────────────────────────
    // 從二級分析輸出目錄找到 ensemble VCF
    // 路徑規則：{input_dir}/04_snv_indel/{sample_id}.ensemble.fixed.vcf.gz

    def ensemble_vcf = file(
        "${params.input_dir}/04_snv_indel/${params.sample_id}.ensemble.fixed.vcf.gz"
    )
    def ensemble_tbi = file(
        "${params.input_dir}/04_snv_indel/${params.sample_id}.ensemble.fixed.vcf.gz.tbi"
    )

    // 確認輸入檔案存在
    if (!ensemble_vcf.exists()) {
        error "[ERROR] 找不到 ensemble VCF：${ensemble_vcf}"
    }
    if (!ensemble_tbi.exists()) {
        error "[ERROR] 找不到 ensemble VCF index：${ensemble_tbi}"
    }

    // 建立 channel：tuple(sample_id, vcf, tbi)
    ensemble_ch = Channel.of(
        tuple(params.sample_id, ensemble_vcf, ensemble_tbi)
    )

    // ── Phase 1：VCF 前處理 ───────────────────────────────
    // 輸入：ensemble VCF（雙 sample column）
    // 輸出：snv_for_annotation.vcf.gz（單 sample，CALLERS tag，PASS only）
    PREPARE_VCF(ensemble_ch)

    // ── Phase 1：VEP Annotation + Pangolin ───────────────
    // 輸入：PREPARE_VCF 輸出的 snv_ch（snv_for_annotation.vcf.gz）
    // 輸出：
    //   vep_ch      → *.vep.vcf.gz（含所有 annotation）
    //   pangolin_ch → *.pangolin.vcf.gz（splice candidate 的 Pangolin 分數）
    SNV_ANNOTATE(PREPARE_VCF.out.snv_ch)

    // ── Phase 1：CSQ 解析 + Pangolin 整合 ────────────────
    // 輸入：VEP VCF + Pangolin VCF
    // 輸出：snv_indel.annotated.tsv（直接輸出到樣本根目錄）
    PARSE_VEP_CSQ(
        SNV_ANNOTATE.out.vep_ch,
        SNV_ANNOTATE.out.pangolin_ch
    )

    // ── Phase 1 後續（尚未實作）──────────────────────────
    // ACMG_CLASSIFY(PARSE_VEP_CSQ.out.tsv_ch)
}
