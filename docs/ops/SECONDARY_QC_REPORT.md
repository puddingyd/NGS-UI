# WES 二級分析 QC summary（DGX2）

## 執行方式與範圍

從 v9.14 起，NGS-UI 產生的 **WES** 執行指令會依序做：

1. 在 DGX2 檢查 QC 腳本、samplesheet、target BED 與 Samtools 容器。
2. 執行原有 Nextflow 二級分析。
3. Nextflow 成功後，在 DGX2 從各樣本 BAM 計算 QC，產生批次 `report_summary.csv`。

UI 仍只產生指令，需由使用者複製到 DGX2 執行；UI server 不讀 BAM 計算，也不遠端啟動 DGX2。這組 capture target／30M／50X 門檻目前套用 WES；WGS 保留原流程。

這是接在 Nextflow 後面的獨立步驟，不需修改 DGX2 的 `main.nf` 或 modules。直接執行舊版 Nextflow 指令不會自動產生本報表；請使用更新後 UI 重新產生的完整指令，或以下補算方式。

## 一次性部署

### 1. 將腳本放到 DGX2

來源為 NGS-UI repository 的 [`scripts/secondary_qc_report.py`](../../scripts/secondary_qc_report.py)。它是單一檔案，使用 Python 3.8+ 標準函式庫，不需在 DGX2 安裝 pip 套件。

建議先上傳至 DGX2 自己的家目錄 `~/secondary_qc_report.py`，再在 DGX2 貼上：

```bash
QC_SCRIPT_DIR=/datalake_Intermediate/pipeline/pipeline_code/scripts
mkdir -p "${QC_SCRIPT_DIR}"
install -m 664 "${HOME}/secondary_qc_report.py" "${QC_SCRIPT_DIR}/secondary_qc_report.py"
python3 "${QC_SCRIPT_DIR}/secondary_qc_report.py" --help
```

最終位置：

```text
/datalake_Intermediate/pipeline/pipeline_code/scripts/secondary_qc_report.py
```

若環境腳本設定了其他 `PIPELINE_CODE`，改放到 `${PIPELINE_CODE}/scripts/secondary_qc_report.py`。檔案及父目錄需讓執行分析的帳號可讀。DGX 的 pipeline code 與 NGS-UI 是不同部署位置；更新 NGS-UI 後，仍需同步這支腳本到 DGX。

Samtools 預設使用 `nextflow config -profile dgx -flat` 解析出的 `process.withName:SAMTOOLS…container`，目前 pipeline 文件設定為 `samtools_1.23.1.sif`。由 Apptainer 或 Singularity 執行，要求 Samtools ≥1.13；QC 使用 CPU，不配置 GPU。

### 2. 更新 UI server

先完成上面的 DGX 腳本部署，再於 DGM 的 NGS-UI checkout 更新 default branch，重新啟動 `ngs-ui` 服務並重新整理瀏覽器。之後重新建立 samplesheet、複製新指令。已經複製出去的舊指令不會自動更新。

## 輸入與輸出

只處理本次 samplesheet 的 `sample`，依第一次出現順序輸出；同一 sample 多個 lane 只產生一列。Sample ID 沿用 samplesheet，不從 FASTQ 名稱推測或附加 `_S1` 等序號。

必要輸入：

```text
<batch>/samplesheet.csv
<batch>/<sample>/02_alignment/<sample>.aligned.sorted.bam
<batch>/<sample>/02_alignment/<sample>.aligned.sorted.bam.bai  # 也支援 .csi 或 .bai
<batch>/<sample>/02_alignment/<sample>.duplicate_metrics.txt
```

target 取自該次解析設定的 `params.wes_targets`，不固定寫死試劑版本，也不使用 reviewer 的疾病 gene panel。BED 依 0-based、半開區間解析，先合併重疊或相鄰區間，避免重複計算；BAM 必須為 coordinate sorted，BED 染色體名稱與座標須符合 BAM header，若 BAM 有 RG/SM 則須與 sample 一致。

報表寫入 **batch 層級**，不是各 sample 的 `pipeline_info`：

```text
<batch>/pipeline_info/report_summary.csv
<batch>/pipeline_info/report_summary.details.json
<batch>/pipeline_info/report_summary.log
<batch>/<sample>/03_alignment_qc/<sample>.report_qc.json
```

CSV 固定 UTF-8、以下八欄、每樣本一列，與範例的欄名／順序／數值格式相同：

```csv
Sample ID,Total reads,Duplicated rate,Mapping rate,On target rate,Mean depth,Uniformity,QC
```

Total reads 為整數；四個比率為兩位小數並附 `%`；Mean depth 為兩位小數。新算法與先前報表的數值可能不同，保留相同的是輸出格式。

## 計算定義

| 欄位 | 算法 | 允收門檻 |
|---|---|---|
| Total reads | BAM primary read ends；R1、R2 各算一條，排除 secondary／supplementary，保留 unmapped、duplicate 與 QC-failed reads | ≥30,000,000 |
| Duplicated rate | Parabricks/Picard `duplicate_metrics.txt`：`(UNPAIRED_READ_DUPLICATES + 2 × READ_PAIR_DUPLICATES) / (UNPAIRED_READS_EXAMINED + 2 × READ_PAIRS_EXAMINED)`；多 library 先加總分子分母 | 僅列出，不參與 PASS/FAIL |
| Mapping rate | mapped primary reads / 全部 primary reads | ≥95% |
| On target rate | 至少一個實際 aligned base 落在 target 內的 primary mapped reads / 全部 primary mapped reads；每 read 只算一次，保留 duplicate 與 QC-failed reads，不套 MQ/BQ 篩選 | ≥40% |
| Mean depth | target 所有 bases 的深度總和 / target 合併後總長度，包含 0X bases | ≥50X |
| Uniformity | 深度 ≥ 未四捨五入平均深度 ×20% 的 target bases / target 合併後總長度 | ≥90% |

Mean depth／Uniformity 共用同一份逐鹼基深度：

- MQ ≥20、BQ ≥20。
- 排除 unmapped、secondary、supplementary、QC-failed、duplicate（排除 flags `3844`）。
- 使用 `samtools depth -s`，成對 reads 的重疊部分依 Samtools 規則只計一次。
- 不將 deletion／reference skip 當成覆蓋鹼基，不設定最大深度截斷。
- 完全沒有 alignment 的區間或染色體仍計入 target 分母，深度為 0。
- Uniformity 用精確整數計算 cutoff：`ceil(depth_sum / (5 × target_bases))`；所有 PASS/FAIL 使用未四捨五入數值。
- 整個 target 深度為 0 時 Uniformity 列 `NA`，QC 為 FAIL，避免把 0X 樣本顯示成 100% 均勻。

On target 計數先透過 Samtools 的 BED union iterator 取得 reads，再驗證 CIGAR 的 `M/=/X` 區段確實與 target 相交；只以 deletion／skip 跨過 target 的 read 不算命中，也不把一條 read 跨多個 target 算成多條。

實作依據：[Samtools depth](https://www.htslib.org/doc/samtools-depth.html)、[Samtools view](https://www.htslib.org/doc/samtools-view.html)、[Picard DuplicationMetrics](https://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics)。

## 完成、失敗與重跑

- **PASS**：五項門檻全部達標。
- **FAIL**：計算完成但至少一項未達標。程序仍成功結束；`DONE` 表示分析與報表完成，不能解讀為所有樣本 QC 通過。
- **ERROR**：缺檔、BAM/BED 不相容或工具執行失敗；該列數值為 `NA`。繼續處理其餘樣本，最後 exit 2，DGX runner 顯示 `FAILED: QC report`。
- preflight 失敗會在 Nextflow 啟動前停止；Nextflow 失敗不執行 QC。tmux pane 仍保留供查看，`report_summary.log` 保留計算進度，details JSON 保留失敗原因與方法。
- 同 batch 的 QC 以檔案鎖避免同時寫入。CSV 在全批處理結束後以原子替換寫入；中斷時舊 CSV 可能仍在，須同時看 details JSON 的 `state`／時間與 runner 狀態，不能將舊 CSV 當成新結果。
- QC 逐樣本順序執行，串流讀取 SAM／depth，不將大型逐鹼基文字存成永久檔。執行時間取決於 BAM 大小與 DGX 儲存 I/O，目前未用真實 DGX 批次測量。
- 每 sample 的小型 JSON 保存原始計數、失敗項目、方法與輸入簽章。相同 BAM/index/duplicate metrics 的路徑、大小、mtime/ctime，及 BED hash、腳本 hash、Samtools 版本／容器簽章未變時重用；`--force` 可強制重算。這是 metadata-based cache，不會每次對大型 BAM 算完整 checksum。

## 已完成 WES 批次補算

不需要重跑 Nextflow。以下命令在 DGX2 執行，請改 `QC_BATCH`。這裡明確指定已確認的 target 和 Samtools SIF；補算舊批次前須確認它們是**該批次實際使用的版本**。

```bash
bash <<'QC_BASH'
set -euo pipefail
umask 0002
QC_BASE=/datalake_Intermediate/pipeline
QC_BATCH=260907_WES
QC_OUT="${QC_BASE}/nextflow_output/${QC_BATCH}"
QC_SCRIPT="${QC_BASE}/pipeline_code/scripts/secondary_qc_report.py"
QC_TARGET="${QC_BASE}/reference/hg38/Illumina_Exome_TargetedRegions_v1.2.hg38.bed"
QC_SIF="${QC_BASE}/nextflow_containers/samtools_1.23.1.sif"

python3 "${QC_SCRIPT}" \
    --out-dir "${QC_OUT}" --samplesheet "${QC_OUT}/samplesheet.csv" \
    --target-bed "${QC_TARGET}" --samtools-sif "${QC_SIF}" --check-only

python3 "${QC_SCRIPT}" \
    --out-dir "${QC_OUT}" --samplesheet "${QC_OUT}/samplesheet.csv" \
    --target-bed "${QC_TARGET}" --samtools-sif "${QC_SIF}" \
    2>&1 | tee "${QC_OUT}/pipeline_info/report_summary.log"
QC_BASH
```

若該批次已由新版 UI 留下 `nextflow_launch/<batch>/secondary_qc.config.flat`，可用 `--nextflow-config /完整路徑/secondary_qc.config.flat` 取代上面的 `--target-bed` 與 `--samtools-sif`。若使用主機原生 Samtools，改成 `--samtools /完整路徑/samtools`。補算指令只使用 CPU，無需重新 source 會清除 GPU locks 的 pipeline 環境腳本；執行前 `apptainer` 或 `singularity` 須已可從 PATH 找到。

## 開發驗證

```bash
python3 -m pytest -q tests/test_secondary_qc_report.py tests/test_secondary_analysis.py
```

含已知答案的合成 BAM（需本機有 Samtools）及實際執行所生成 Bash 的 stub 測試，驗證各種 flags、BQ/MQ、paired overlap、CIGAR D/N、BED union／0X、精確門檻、快取失效、部分 ERROR，以及 preflight → Nextflow → QC 的成功／失敗順序。本機測試不能取代首次 DGX 真實批次的部署驗證。
