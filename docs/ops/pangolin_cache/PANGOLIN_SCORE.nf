// =====================================================================
// PANGOLIN_SCORE — drop-in replacement for modules/snv_annotation.nf
// =====================================================================
// Differences from the current production process (everything else is
// unchanged, including the output file names and contents read downstream):
//
//   1. params.containsKey(...) → plain params.xxx ?: default.
//      A method call on `params` makes Nextflow hash the WHOLE params map into
//      the task hash, and every NGS-UI job has its own out_dir / samplesheet,
//      so PANGOLIN_SCORE could never be -resume cached. Reading the two
//      properties directly hashes only their values. (All profiles declare
//      both params; if one ever does not, Nextflow only prints
//      "WARN: Access to undefined parameter" and the default is used.)
//
//   2. Persistent score cache (scripts/pangolin_cache.py), enabled when
//      params.pangolin_cache_db is a non-empty path. Cached variants are taken
//      from the cache and only never-seen variants go to the GPU. With
//      params.pangolin_cache_db = '' the process runs Pangolin exactly as before.
//
//   3. Pangolin options live in one place (pangolin_distance / pangolin_mask)
//      so the command line and the cache namespace can never drift apart.
//
// Replacing the process body changes its script hash once, so the first run
// after deployment re-scores every sample (with the cache filling up); later
// runs are cached by Nextflow and by the score cache.
// =====================================================================

process PANGOLIN_SCORE {
    label 'process_gpu'

    // ── 容器依機器的 GPU 世代切換（params.pangolin_sif，每個 profile 都要宣告）──
    //   （原註解保留：cu121 → dgm + dgx；cu130 → 開發機 local）
    container "${params.sif_dir}/${params.pangolin_sif}"

    publishDir "${params.out_dir}/${sample_id}/02_pangolin", mode: 'copy'

    input:
    tuple val(sample_id), path(vep_vcf), path(vep_tbi)

    output:
    tuple val(sample_id),
          path("${sample_id}.pangolin.vcf.gz"),
          path("${sample_id}.pangolin.vcf.gz.tbi"),
          emit: pangolin_out

    script:
    def use_gpu       = params.use_gpu_pangolin == null ? true
                        : params.use_gpu_pangolin.toString().toLowerCase() == 'true'
    // ⚠️ 不可用 params.containsKey(...) / params.get(...) / params[...]：
    //    對 params 呼叫方法會讓 Nextflow 把「整個 params」算進 task hash，
    //    out_dir / samplesheet 每個 job 都不同 → PANGOLIN_SCORE 永遠無法 -resume。
    //    三個 profile 都已宣告這兩個參數；萬一沒宣告，只會出現
    //    "WARN: Access to undefined parameter"，並使用 ?: 後面的預設值。
    def use_lock      = use_gpu && (params.use_gpu_lock ?: false)
    def lock_script   = params.gpu_lock_script
    def unlock_script = params.gpu_unlock_script
    def num_gpus      = (params.pangolin_num_gpus ?: 1)

    // Pangolin 參數只寫一次：命令列與 cache namespace 共用，避免兩邊不一致
    def pangolin_distance = 50
    def pangolin_mask     = 'True'

    // 分數快取：params.pangolin_cache_db 為空字串 → 停用，行為與舊版完全相同
    def cache_db   = (params.pangolin_cache_db ?: '').toString()
    def cache_args = cache_db ? """--cache ${cache_db} \\
            --sif-name ${params.pangolin_sif} \\
            --pangolin-db ${params.pangolin_db} \\
            --ref-fasta ${params.ref_fasta} \\
            --distance ${pangolin_distance} \\
            --mask ${pangolin_mask} \\
            --sample ${sample_id}""" : ''

    def gpu_block = !use_gpu ? """
    export CUDA_VISIBLE_DEVICES=""
    echo "[PANGOLIN] ${sample_id} CPU 模式（use_gpu_pangolin=false）" >&2
    """ : use_lock ? """
    eval \$(bash ${lock_script} ${num_gpus})
    export CUDA_VISIBLE_DEVICES=\${MY_GPUS}
    echo "[PANGOLIN] ${sample_id} 取得 GPU \${MY_GPUS}" >&2
    trap "bash ${unlock_script} \${MY_GPUS}; echo '[PANGOLIN] ${sample_id} 釋放 GPU \${MY_GPUS}' >&2" EXIT
    """ : """
    echo "[PANGOLIN] ${sample_id} 使用 config 指定的 GPU（未啟用 GPU lock）" >&2
    """

    def run_block = cache_db ? """
        # Step 3a：查分數快取 → 只把沒算過的變異送 Pangolin
        python3 ${params.scripts_dir}/pangolin_cache.py split \\
            ${cache_args} \\
            --in splice_for_pangolin.vcf \\
            --miss-out pangolin_miss.vcf \\
            --hits-out pangolin_hits.tsv

        MISS_COUNT=\$(grep -vc '^#' pangolin_miss.vcf || true)
        if [ "\${MISS_COUNT}" -gt 0 ]; then
            # Step 3b：執行 Pangolin（GPU）；output_file 是 prefix，Pangolin 自動加 .vcf
            pangolin \\
                pangolin_miss.vcf \\
                ${params.ref_fasta} \\
                ${params.pangolin_db} \\
                pangolin_miss_scored \\
                --distance ${pangolin_distance} \\
                -m ${pangolin_mask}
        fi

        # Step 3c：快取值 + 新算的分數按原順序組回，並把新分數寫入快取
        python3 ${params.scripts_dir}/pangolin_cache.py merge \\
            ${cache_args} \\
            --in splice_for_pangolin.vcf \\
            --hits pangolin_hits.tsv \\
            --pangolin-out pangolin_miss_scored.vcf \\
            --out ${sample_id}.pangolin.vcf
    """ : """
        # Step 3：執行 Pangolin（GPU）
        # output_file 是 prefix，Pangolin 自動加 .vcf
        # 不要給 ${sample_id}.pangolin.vcf，否則輸出變成 .vcf.vcf
        pangolin \\
            splice_for_pangolin.vcf \\
            ${params.ref_fasta} \\
            ${params.pangolin_db} \\
            ${sample_id}.pangolin \\
            --distance ${pangolin_distance} \\
            -m ${pangolin_mask}
    """
    """
    ${gpu_block}

    # Step 1：從 VEP 輸出中篩選 splice candidate
    bcftools view -h ${vep_vcf} > splice_header.vcf
    bcftools view -H ${vep_vcf} | \\
        awk -F'\\t' '\$8 ~ /splice/' > splice_body.vcf

    SPLICE_COUNT=\$(wc -l < splice_body.vcf)
    echo "[PANGOLIN] ${sample_id}：splice candidates = \${SPLICE_COUNT}" >&2

    if [ "\${SPLICE_COUNT}" -eq 0 ]; then
        echo "[PANGOLIN] 無 splice candidate，產生空輸出" >&2
        bgzip -c splice_header.vcf > ${sample_id}.pangolin.vcf.gz

    else
        # Step 2：移除 CSQ + 只保留標準染色體（避免 Pangolin segfault）
        cat splice_header.vcf splice_body.vcf \\
            | bcftools annotate -x INFO/CSQ \\
            | grep -E '^#|^chr([0-9]+|[XYM])\t' \\
            > splice_for_pangolin.vcf
${run_block}
        bgzip -c ${sample_id}.pangolin.vcf > ${sample_id}.pangolin.vcf.gz
    fi

    tabix -p vcf ${sample_id}.pangolin.vcf.gz
    echo "[PANGOLIN] ${sample_id} 完成" >&2
    """
}
