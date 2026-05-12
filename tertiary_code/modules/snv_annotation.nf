/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - SNV Annotation Module
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
 * modules/snv_annotation.nf
 * =========================
 * 目的：
 *   對 prepare_vcf 輸出的 snv_for_annotation.vcf.gz 進行兩階段 annotation：
 *
 *   Process 1 - VEP_ANNOTATE：
 *     執行 Ensembl VEP 115，帶入以下 plugins 和 custom annotation：
 *       - dbNSFP 4.9c（BayesDel、AlphaMissense、ESM1b 等）
 *       - LOFTEE GRCh38（LoF HC/LC 判定，PVS1 基礎）
 *         * loftee 資料目錄透過 --bind 掛載至 /opt/vep/Plugins/loftee_data
 *         * gerp BigWig 透過 --bind 直接掛載至 /opt/vep/Plugins/（LOFTEE bug workaround）
 *       - LoFtool（基因 LoF 容忍度）
 *       - ClinVar 20260418（CLNSIG、CLNREVSTAT、CLNDN、CLNSIGCONF）
 *       - gnomAD（VEP cache 內建，--af_gnomadg + --af_gnomade）
 *
 *     annotation 旗標（CSQ 欄位內容）：
 *       - --hgvs    ：輸出 HGVSc / HGVSp（TSV 主表必要欄位，無此旗標欄位為空）
 *       - --symbol  ：輸出 HGNC gene symbol（GENE 欄位，無此旗標只有 Ensembl ID）
 *       - --numbers ：輸出 EXON / INTRON 編號（格式 2/15，PVS1 NMD escape 判斷需要）
 *       - --canonical：輸出 CANONICAL 旗標（transcript 備用選取依據）
 *       - --biotype ：輸出 BIOTYPE（PVS1 需確認 protein_coding）
 *       - --tsl     ：輸出 TSL（--pick_order 的 tsl 層才有資料可比較）
 *       - --appris  ：輸出 APPRIS（--pick_order 的 appris 層才有資料可比較）
 *       - --uniprot ：輸出 SWISSPROT / TREMBL（PM1 UniProt domain 判斷依據）
 *       - --domains ：輸出 DOMAINS（PM1 蛋白質 domain 資訊）
 *       - --safe    ：plugin 報錯時 VEP 正確退出（臨床 pipeline 必要）
 *
 *     transcript 策略：
 *       - --mane：標記 MANE Select / MANE Plus Clinical
 *       - --flag_pick：標記代表 transcript（PICK 欄位值為 1）
 *       - --pick_order：mane_select > mane_plus_clinical > canonical > appris > tsl > ...
 *       - 注意：不用 --pick（--pick 只保留一個 transcript，無法展開所有 MANE）
 *               只用 --flag_pick，保留所有 transcript，以 PICK=1 標記代表 transcript
 *     PICK 欄位說明：
 *       VCF 的 CSQ 欄位中，PICK 是第 22 個 pipe-separated 值（從 1 開始數）
 *       值為 "1" 表示此 transcript 是代表 transcript，值為空表示非代表
 *       加入 --hgvs/--symbol 等旗標後，PICK 的欄位位置會往後移，
 *       parse_vep_csq.py 應從 VCF header 的 CSQ 欄位定義動態解析位置，不要硬編碼
 *
 *   Process 2 - PANGOLIN_SCORE：
 *     對 VEP 輸出中帶有 splice consequence 的 variant 執行 Pangolin inference。
 *     取代 SpliceAI（SpliceAI 商用需付費授權，Pangolin 為 GPL-3.0）。
 *     GPU inference（--nv 由 config 的 process_gpu label containerOptions 設定）。
 *
 * 輸入：
 *   tuple val(sample_id), path(snv_vcf), path(snv_tbi)
 *
 * 輸出：
 *   vep_ch:      tuple val(sample_id), path("*.vep.vcf.gz"), path("*.vep.vcf.gz.tbi")
 *   pangolin_ch: tuple val(sample_id), path("*.pangolin.vcf.gz"), path("*.pangolin.vcf.gz.tbi")
 *
 * 已知踩雷：
 *   - --stats_file 和 --no_stats 不能同時使用，本 module 統一用 --no_stats（較快）
 *   - LOFTEE gerp_dist.pl bug：路徑拼接缺少 /，需要透過 --bind 直接掛載 gerp bw
 *   - tabix 必須在 VEP 輸出之後才能建 index（VEP 加 --compress_output bgzip 自動壓縮）
 *   - 加入 --hgvs/--symbol/--numbers 等旗標後，CSQ 欄位的 column 數增加，
 *     PICK 的欄位位置不再是固定第 22 欄，parse_vep_csq.py 必須從 VCF header
 *     動態解析 CSQ format（##INFO=<ID=CSQ,...,Format="..."> 那行）
 */

// ──────────────────────────────────────────────────────────────
// Process 1：VEP Annotation
// ──────────────────────────────────────────────────────────────

process VEP_ANNOTATE {
    label 'process_high'

    container "${params.sif_dir}/vep_115.sif"

    // LOFTEE bind mount 設定：
    //   loftee_data → 整個 loftee 資料目錄掛載至容器內 /opt/vep/Plugins/loftee_data
    //   gerp bw     → 直接掛載至 /opt/vep/Plugins/（LOFTEE gerp_dist.pl 硬編碼路徑的 workaround）
    //
    // 注意：apptainer_base_opts 提供 --bind /scratch,/data 基礎掛載
    //       這裡在基礎上追加 loftee 相關的兩個 bind
    containerOptions "${params.apptainer_base_opts} \
        --bind ${params.loftee_dir}:/opt/vep/Plugins/loftee_data \
        --bind ${params.loftee_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw:/opt/vep/Plugins/gerp_conservation_scores.homo_sapiens.GRCh38.bw"

    publishDir "${params.out_dir}/${sample_id}/01_vep", mode: 'copy'

    input:
    tuple val(sample_id), path(snv_vcf), path(snv_tbi)

    output:
    tuple val(sample_id),
          path("${sample_id}.vep.vcf.gz"),
          path("${sample_id}.vep.vcf.gz.tbi"),
          emit: vep_out

    script:
    """
    vep \\
        --input_file ${snv_vcf} \\
        --output_file ${sample_id}.vep.vcf.gz \\
        --vcf \\
        --compress_output bgzip \\
        \\
        --offline \\
        --cache \\
        --dir_cache ${params.vep_cache} \\
        --dir_plugins /opt/vep/Plugins \\
        --assembly GRCh38 \\
        --fasta ${params.ref_fasta} \\
        --fork ${task.cpus} \\
        \\
        --hgvs \\
        --symbol \\
        --numbers \\
        --canonical \\
        --biotype \\
        --tsl \\
        --appris \\
        --uniprot \\
        --domains \\
        \\
        --mane \\
        --flag_pick \\
        --pick_order mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length \\
        \\
        --plugin dbNSFP,${params.dbnsfp},\\
BayesDel_noAF_score,BayesDel_noAF_pred,\\
AlphaMissense_score,AlphaMissense_pred,\\
ESM1b_score,ESM1b_pred,\\
VARITY_R_score,\\
SIFT_score,SIFT_pred,\\
DANN_score,\\
PHACTboost_score,\\
phyloP100way_vertebrate,\\
GERP++_RS,\\
gnomAD_exomes_AF,gnomAD_exomes_EAS_AF \\
        \\
        --plugin LoF,\\
loftee_path:/opt/vep/Plugins/,\\
human_ancestor_fa:/opt/vep/Plugins/loftee_data/human_ancestor.fa.gz,\\
conservation_file:/opt/vep/Plugins/loftee_data/loftee.sql \\
        \\
        --plugin LoFtool,/opt/vep/Plugins/loftee_data/LoFtool_scores.txt \\
        \\
        --custom file=${params.clinvar},\\
short_name=ClinVar,\\
format=vcf,\\
type=exact,\\
coords=0,\\
fields=CLNSIG%CLNREVSTAT%CLNDN%CLNSIGCONF \\
        \\
        --af_gnomadg \\
        --af_gnomade \\
        --force_overwrite \\
        --no_stats \\
        --safe

    # 建立 tabix index（VEP 已用 --compress_output bgzip 壓縮，直接 tabix）
    tabix -p vcf ${sample_id}.vep.vcf.gz

    echo "[VEP_ANNOTATE] ${sample_id} 完成" >&2
    bcftools stats ${sample_id}.vep.vcf.gz | grep "^SN" >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 2：Pangolin Splice Scoring
// ──────────────────────────────────────────────────────────────

process PANGOLIN_SCORE {
    label 'process_gpu'

    container "${params.sif_dir}/pangolin_1.0.0.sif"

    publishDir "${params.out_dir}/${sample_id}/02_pangolin", mode: 'copy'

    input:
    tuple val(sample_id), path(vep_vcf), path(vep_tbi)

    output:
    tuple val(sample_id),
          path("${sample_id}.pangolin.vcf.gz"),
          path("${sample_id}.pangolin.vcf.gz.tbi"),
          emit: pangolin_out

    script:
    """
    # Step 1：從 VEP 輸出中篩選 splice candidate
    # bcftools view -h → 只取 header 行
    # bcftools view -H → 只取 variant 行，awk 篩 INFO 欄含 "splice" 字眼
    bcftools view -h ${vep_vcf} > splice_header.vcf
    bcftools view -H ${vep_vcf} | \\
        awk -F'\\t' '\$8 ~ /splice/' > splice_body.vcf

    SPLICE_COUNT=\$(wc -l < splice_body.vcf)
    echo "[PANGOLIN] ${sample_id}：splice candidates = \${SPLICE_COUNT}" >&2

    if [ "\${SPLICE_COUNT}" -eq 0 ]; then
        # 無 splice candidate，產生空輸出
        echo "[PANGOLIN] 無 splice candidate，產生空輸出" >&2
        bgzip -c splice_header.vcf > ${sample_id}.pangolin.vcf.gz

    else
        # Step 2：合併 header + splice body，送進 Pangolin
        cat splice_header.vcf splice_body.vcf > splice_for_pangolin.vcf

        # Step 3：執行 Pangolin（GPU）
        # output_file 是 prefix，Pangolin 自動加 .vcf
        # 不要給 ${sample_id}.pangolin.vcf，否則輸出變成 .vcf.vcf
        # IUPAC ambiguity code（Y/R/W 等）已在 sif build 時 patch 處理，
        # 碰到非 ACGT 字元一律當 N（全零 encoding），不會 crash
        pangolin \\
            splice_for_pangolin.vcf \\
            ${params.ref_fasta} \\
            ${params.pangolin_db} \\
            ${sample_id}.pangolin \\
            --distance 50 \\
            -m True

        bgzip -c ${sample_id}.pangolin.vcf > ${sample_id}.pangolin.vcf.gz
    fi

    tabix -p vcf ${sample_id}.pangolin.vcf.gz
    echo "[PANGOLIN] ${sample_id} 完成" >&2
    """
}

// ──────────────────────────────────────────────────────────────
// 組合 workflow（供 main_tertiary.nf 呼叫）
// ──────────────────────────────────────────────────────────────

workflow SNV_ANNOTATE {
    take:
    snv_ch   // tuple val(sample_id), path(snv_vcf), path(snv_tbi)

    main:
    VEP_ANNOTATE(snv_ch)
    PANGOLIN_SCORE(VEP_ANNOTATE.out.vep_out)

    emit:
    vep_ch      = VEP_ANNOTATE.out.vep_out
    pangolin_ch = PANGOLIN_SCORE.out.pangolin_out
}
