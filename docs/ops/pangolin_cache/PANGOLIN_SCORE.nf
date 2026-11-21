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
    //   Pangolin 用 PyTorch，而 PyTorch 的 prebuilt wheel 只包特定 compute
    //   capability 的 kernel，跨 GPU 世代不相容（cubin 只在同一個 major 版本內
    //   向前相容）。實測官方 wheel 的 arch list：
    //     cu121 → sm_50…sm_90   （含 sm_70 = V100，無 sm_120）
    //     cu128 → sm_75…sm_120  （含 sm_120 = Blackwell，**已砍 sm_70**；cu130 同理）
    //   → 沒有任何一顆 prebuilt wheel 同時涵蓋兩者，只能分成兩顆容器：
    //     pangolin_cu121_1.0.0.sif → production：dgm + dgx（V100 = sm_70）
    //     pangolin_cu130_1.0.0.sif → 開發機 local（RTX PRO 6000 Blackwell = sm_120）
    //   檔名標的是 CUDA 變體而非機器，因為差異的本質就是 wheel 的 CUDA 變體。
    //   建置方式與 arch 守門員見 DEVELOPMENT_NOTES.md「容器建立 → Pangolin」。
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
    // ── GPU lock（DGX-2 共用機器）─────────────────────────────────────────
    //   與二級分析同一套機制：gpu_lock.sh 搶下 N 張空閒 V100 並輸出可 eval 的環境變數，
    //   trap EXIT 時以 gpu_unlock.sh 歸還（即使 process 失敗也會釋放）。
    //   Apptainer 預設會把 host 環境變數帶進容器，故 CUDA_VISIBLE_DEVICES 會被 Pangolin
    //   看見；這裡再明確 export 一次，避免 lock script 只設了 MY_GPUS。
    //   非 DGX（local / dgm）use_gpu_lock=false → 走 config 既有的 GPU 設定，不插入此段。
    // ── GPU / CPU 切換（--use_gpu_pangolin）──────────────────────────────
    //   Pangolin 以 PyTorch 實作，會自動偵測 CUDA；把 CUDA_VISIBLE_DEVICES 設為空字串
    //   即可強制走 CPU（慢但結果相同），讓 pipeline 能部署到沒有 GPU 的環境。
    //   容器的 --nv 由 config 的 process_gpu label 依同一參數決定是否加上。
    def use_gpu       = params.use_gpu_pangolin == null ? true
                        : params.use_gpu_pangolin.toString().toLowerCase() == 'true'
    // ⚠️ use_gpu_lock / pangolin_num_gpus 只宣告在 profile 內（不能放全域 params 區塊，
    //    那個區塊在 profiles 之後會蓋掉 profile 的值）。三個 profile 都已宣告；萬一沒宣告，
    //    Nextflow 只會噴 "WARN: Access to undefined parameter"，並使用 ?: 後面的預設值。
    //    ★ 不可改用 params.containsKey(...) / params.get(...) / params[...]：
    //      對 params 呼叫方法會讓 Nextflow 把「整個 params」算進 task hash，
    //      out_dir / samplesheet 每個 job 都不同 → PANGOLIN_SCORE 永遠無法 -resume。
    // 不用 GPU 就不必搶卡
    def use_lock      = use_gpu && (params.use_gpu_lock ?: false)
    def lock_script   = params.gpu_lock_script
    def unlock_script = params.gpu_unlock_script
    def num_gpus      = (params.pangolin_num_gpus ?: 1)

    // ── Pangolin 參數（命令列與分數快取的 namespace 共用，避免兩邊不一致）──
    def pangolin_distance = 50
    def pangolin_mask     = 'True'

    // ── Pangolin 分數快取（scripts/pangolin_cache.py）────────────────────
    //   Pangolin 分數只取決於變異本身、參考序列、gene model DB、容器與參數，與樣本無關；
    //   cohort 內約 98% 的 splice candidate 先前都已算過。快取依上述條件分區，
    //   換容器 / DB / 參考序列 / 參數會自動使用新的空白分區。
    //   params.pangolin_cache_db 設為 '' → 停用，執行內容與舊版完全相同。
    def cache_db   = (params.pangolin_cache_db ?: '').toString()
    def cache_args = cache_db ? """--cache ${cache_db} \\
            --sif-name ${params.pangolin_sif} \\
            --pangolin-db ${params.pangolin_db} \\
            --ref-fasta ${params.ref_fasta} \\
            --distance ${pangolin_distance} \\
            --mask ${pangolin_mask} \\
            --sample ${sample_id}""" : ''
    def run_block = cache_db ? """
        # Step 3a：查分數快取，只把沒算過的變異寫進 pangolin_miss.vcf
        python3 ${params.scripts_dir}/pangolin_cache.py split \\
            ${cache_args} \\
            --in splice_for_pangolin.vcf \\
            --miss-out pangolin_miss.vcf \\
            --hits-out pangolin_hits.tsv

        # Step 3b：執行 Pangolin（GPU），只跑沒算過的變異；全部查得到就不啟動模型
        # output_file 是 prefix，Pangolin 自動加 .vcf
        # IUPAC ambiguity code（Y/R/W 等）已在 sif build 時 patch 處理，
        # 碰到非 ACGT 字元一律當 N（全零 encoding），不會 crash
        MISS_COUNT=\$(grep -vc '^#' pangolin_miss.vcf || true)
        if [ "\${MISS_COUNT}" -gt 0 ]; then
            pangolin \\
                pangolin_miss.vcf \\
                ${params.ref_fasta} \\
                ${params.pangolin_db} \\
                pangolin_miss_scored \\
                --distance ${pangolin_distance} \\
                -m ${pangolin_mask}
        fi

        # Step 3c：快取值 + 新算的分數依原順序組回 ${sample_id}.pangolin.vcf，新分數寫入快取
        #   （下游 parse_vep_csq.py 只讀 CHROM/POS/REF/ALT 與 INFO/Pangolin，兩者與直接跑相同）
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
        # IUPAC ambiguity code（Y/R/W 等）已在 sif build 時 patch 處理，
        # 碰到非 ACGT 字元一律當 N（全零 encoding），不會 crash
        pangolin \\
            splice_for_pangolin.vcf \\
            ${params.ref_fasta} \\
            ${params.pangolin_db} \\
            ${sample_id}.pangolin \\
            --distance ${pangolin_distance} \\
            -m ${pangolin_mask}
"""
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
    """
    ${gpu_block}

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
        # Step 2：移除 CSQ + 過濾 alt/random contig，再送 Pangolin
        # 原因 1：WGS CSQ 可能長達 270KB，Pangolin parse 時 segfault
        # 原因 2：Pangolin gencode DB 沒有 alt/random contig 的 gene model，
        #         碰到 chr*_alt / chr*_random / chrUn_* 也會 segfault
        # bcftools view --regions 需要 bgzip+tabix 輸入，無法直接接 pipe
        # 改用 grep + awk：先移除 CSQ，再只保留標準染色體（chr1-22,X,Y,M）
        # alt/random contig 的 Pangolin gencode DB 沒有 gene model，會 segfault
        cat splice_header.vcf splice_body.vcf \
            | bcftools annotate -x INFO/CSQ \
            | grep -E '^#|^chr([0-9]+|[XYM])\t' \
            > splice_for_pangolin.vcf

${run_block}
        bgzip -c ${sample_id}.pangolin.vcf > ${sample_id}.pangolin.vcf.gz
    fi

    tabix -p vcf ${sample_id}.pangolin.vcf.gz
    echo "[PANGOLIN] ${sample_id} 完成" >&2
    """
}
