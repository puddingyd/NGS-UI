# 檔案清單：路徑驗證與共用索引規劃

2026-11-05。**目前只有唯讀盤點工具；modal、正式索引與排程整合尚未實作。**

先由操作者在 server 驗證，再接右上角「三級分析 → 檔案清單 → 上傳個案清單」。

## 1. 搜尋位置

以下是現有程式的預設值；server 的環境設定可能覆寫。`S` 是原始 sequencing sample ID，`C` 是 UI 個案 ID，`R` 是定序 run，`B` 是 in-house 分析 batch。不把 `C` 的後綴任意剪掉當作 `S`。

| 類型 | 搜尋根與相對位置 | 依據 |
|---|---|---|
| WES FASTQ | `/home/datalake_Raw/NextSeq2000/**/{S}_S數字_R{1,2}_001.fastq.gz`；包含 run 下 `Analysis/...` | `secondary_analysis.list_wes_fastqs()` |
| 重分析 FASTQ | `/home/datalake_Raw/Other/Reanalysis/**/{資料夾sample}/{來源S}.R{1,2}.clean.fastq.gz`；也查 `/datalake_Raw/Other/Reanalysis/` | 同上；sample 採父資料夾，保留原始檔名的 source ID |
| WGS FASTQ | `/home/datalake_Raw/Novaseq/{R}/fastq.gz/{S}_S數字_L00{1..4}_R{1,2}_001.fastq.gz`；現有搜尋其實遞迴找 lane 格式 | `secondary_analysis.list_wgs_fastqs()`；四 lane 共八檔，不拿第一對代替全部，也不拿 merged 當 lane |
| In-house BAM | `/home/datalake_Intermediate/pipeline/nextflow_output/{B}/{S}/02_alignment/{S}.aligned.sorted.bam` | `routers/igv.py`；額外深層 batch 由實際目錄驗證 |
| DRAGEN BAM | `/home/datalake_Raw/Novaseq/{R}/bam/{S}.bam`；另查 `/home/datalake_Intermediate/n102968/{R}/bam/` | 同上；`*.repeats.bam` 分類為其他 BAM，不能誤當主 BAM |
| BAM index | 同目錄 `{name}.bam.bai` 或 `{name}.bai`；另列存在的 CSI | `igv._bam_index_for()` 支援兩種 BAI |
| In-house SNV/Indel VCF | `/home/datalake_Intermediate/pipeline/nextflow_output/**/{S}/04_snv_indel/{S}.ensemble.fixed.vcf.gz` | `dragen_jobs.list_inhouse_vcfs()` |
| DRAGEN SNV/Indel VCF | `/home/datalake_Raw/Novaseq/{R}/vcf.gz/{S}.hard-filtered.vcf.gz`；另查 `/home/datalake_Intermediate/n102968/{R}/vcf.gz/` | `dragen_jobs.list_dragen_vcfs()` |
| In-house CNV / SV | 同一 sample 的 `05_cnv_sv/{S}.gcnv.vcf.gz`、`{S}.delly.vcf.gz`；WGS CNV 另查 `{S}.call.cns` / `*.cnr`（不是 VCF，實際輸出待驗證） | 現有三級 siblings／`annotate_inhouse_cnv_sv_annotsv.sh`；CNVkit 路徑依三級 pipeline 文件補查 |
| DRAGEN CNV / SV | 同一 run 的 `vcf.gz/{S}.cnv.vcf.gz`、`{S}.sv.vcf.gz`、`{S}.cnv_sv.vcf.gz` | `annotate_dragen_cnv_annotsv.sh` 用獨立 CNV/SV；CNV rescue 用合併檔，三者不可混稱同一檔 |
| Mito VCF | In-house `07_mitochondria/{S}.mito.vcf.gz`；其他實際出現 VCF 保留待分類 | 三級輸入 siblings |
| 三級衍生 CNV / SV | `/home/datalake_Intermediate/pipeline/tertiary_output/{C}/06_cnv_sv/{S}.{cnv,sv}.annotated.tsv`；`08_postprocessing/{C}.cnv.review.tsv`、舊 `cnv.rescued.annotated.tsv` | 這些是 TSV；不當作原始 VCF。08 的 UI ID 用 exact `pipeline_source.json` 連回 source |
| 其他衍生 VCF | 同一個案 `08_postprocessing/{C}.ploidy.vcf.gz`、`{C}.vcf_from_tsv.vcf.gz` | 標示衍生用途，不當作三級輸入或 PRS gVCF |
| DRAGEN gVCF | 見下面五個 PRS 模板 | `prs_jobs.DEFAULT_TEMPLATES` |
| In-house gVCF | 無預設路徑；由 PRS config 的 `gvcf_templates.inhouse` 決定 | 不猜 caller、不把一般 VCF 當 gVCF |
| VCF / gVCF index | 對應檔案後的 `.tbi` / `.csi` | 只列存在；存在不代表索引有效 |

DRAGEN gVCF 五個現有模板（`anchor_dir` 為上述 SNV VCF 的 `vcf.gz`，`run_dir` 為其上層）：

```text
{anchor_dir}/{S}.hard-filtered.gvcf.gz
{run_dir}/gvcf.gz/{S}.hard-filtered.gvcf.gz
{run_dir}/other/{S}/germline_seq/{S}.hard-filtered.gvcf.gz
{anchor_dir}/{S}.gvcf.gz
{run_dir}/gvcf.gz/{S}.gvcf.gz
```

有自訂 `gvcf_templates` 時，PRS 實際候選只用自訂 mapping；盤點仍列出目錄內發現的全部 gVCF，**盤點存在不等於 PRS 允許使用**。若需查舊三級資料根，額外以 `--root` 指定既有 legacy root，不掃 `nf_work` 或容器暫存。

## 2. Server 驗證方式

開啟 `scripts/paste_probe_sample_files.sh`，全選內容貼到 server shell。它使用 Python 標準函式庫，不必先更新 repo，也不需要 sudo。程式本體由 `scripts/probe_sample_files.py` 產生，兩份內容須同步。

貼上後會建立一個新的 `/tmp/ngs-file-probe.XXXXXX/`，只將 JSON 報告及進度 Log 寫在這裡。讀來源目錄、stat、少量 source marker／PRS config，不解壓 FASTQ/BAM/VCF、不更新平台快取、不啟停服務。

工具使用**目前 shell** 的路徑環境變數，未設定則用上表預設；不假設 shell 與 `ngs-ui.service` 的設定相同。貼上程式會列出 service 設定檔位置供操作者核對，只需將相關 `NGS_UI_*ROOT*`／`NGS_UI_PRS_CONFIG` 路徑設定到 shell。不要貼回包含其他機密的完整 service 環境。

也可直接將 Python 檔複製到 server 執行：

```bash
python3 /tmp/probe_sample_files.py --sample 26T01345 > /tmp/sample-files.json
# --sample 只篩選報告內 sample，仍掃描全部 root；省略可列全部。
# --root /額外資料根 可重複；--prs-config /實際/config.json 指定 PRS mapping。
```

報告要確認：

- `samples` 以 FASTQ 建立；有 raw 但沒有任何分析結果者仍存在。四 lane 正常為八檔，但不把八當硬性篩選條件。
- `incomplete_pairs` 保留缺 R1/R2 的 sample；`unrecognized_fastqs` 保留未辨識檔名，不悄悄捨棄。正式上線前須補齊這些命名對應。
- `files` 保存每個匹配 source ID 的候選路徑、類型、大小與時間；相同 sample 不同 run/batch 全留。這是檔名候選配對，不能證明 BAM/VCF 是某個 FASTQ run 產物。
- 重分析保留資料夾 sample 與原始 source ID；08 的 prefixed source marker 優先於 legacy marker。
- `unmatched_files` 列出沒有 FASTQ 對應的結果；可能是 raw 已搬走／刪除，也可能是命名需要補對應。
- `scan_errors` 列 root 不存在、權限不足、壞 symlink、未跟隨的目錄 symlink。發生時 exit code 為 2，仍輸出已找到的報告；不能宣稱已完整掃描。
- 目錄 symlink 不自動跟隨；確認目標後以 `--root TARGET` 加入。來源只讀，檔案存在／大小不是完整性或單樣本 header 驗證。

建議先核對 `26T01345` 是否對到使用者提供的 `20261028_LH00873_0043_B25FKFHLT4/bam/26T01345.bam` 及同 run 八個 FASTQ，再看一個 in-house WES、一個 in-house WGS、一個重分析及一個尚無 BAM/VCF 的 sample。

## 3. 正式共用索引與每日更新（待實作）

現況：`ngs-ui-secondary-index-update.timer` 與 `ngs-ui-tertiary-index-update.timer` **各自每天 02:00 Asia/Taipei** 執行。各服務有各自 flock，但開窗遇過期仍可能重複排入完整掃描。只是把兩個 worker 放到一個 timer 不足以消除這個問題。

預計改為：

1. 一個 canonical 檔案索引（建議 `data/sample_file_index.sqlite`），一個掃描器。以正規化後的 root 聯集走訪一次，辨識 FASTQ/BAM/VCF/CNV/SV/gVCF 與其 index；記 source、run/batch、lane/read、檔案角色及掃描 generation。索引的 sample 清單來自 FASTQ，檔案表也保留找不到 FASTQ 的結果。
2. 二級讀相同 generation 的 FASTQ 配對視圖，維持完整 R1/R2、lane samplesheet 規則；「檔案清單」讀完整 sample 視圖，含不完整 raw。三級讀相同檔案表的合法 VCF anchor 視圖，**保留目前無 FASTQ 但可分析的 VCF 候選**，避免 raw 已移走的舊樣本失去三級功能。
3. 一個 `ngs-ui-file-index-update.timer`，每天 02:00 台北時間。手動更新、過期 fallback 與每日更新共用一把跨程序鎖及同一個更新協調器；同時請求合併到同一工作，鎖內重新檢查 generation/請求時間，避免等待前一輪後又重掃。
4. 新 generation 完整掃描成功後才在 transaction 中啟用；更新期間讀上一份完整 generation。必要 root 無法讀取時保留舊資料並標示失敗／過期，不拿部分掃描覆蓋造成大量「檔案消失」。掃描成功但舊檔確已刪除才移除該紀錄；已驗證為選用的不存在 root 明確記錄為選用。
5. modal 預設讀索引，提供 sample 搜尋、FASTQ/BAM/VCF/其他分類、各檔完整路徑、run/batch、大小、更新時間與缺配對狀態；四 lane 可展開八檔。同 source 多 run 保留所有候選，未驗證的產物來源關係不自動宣稱。
6. 由操作者安裝新排程並驗證一次完整 generation，再停用兩個舊 timer；現有 API/舊 worker 亦須改走同一協調器，避免舊入口偷偷繼續掃描。既有分析完成後可更新相關 sample 的輸出狀態；每日全掃處理外部新上機及移動／刪除。安裝與服務操作仍由使用者執行。

這份規劃不把 ClinVar/LitVar/HPO 等資料庫更新排程混入檔案索引；它們的資料來源與更新頻率不同。
