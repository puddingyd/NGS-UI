# PGS004879 integrated v2.1 — Execution orchestration add-on

**研究用途，非臨床風險確效。更新：2026-10-31。**

## 1. Scope／設計邊界

本 add-on **只增加 execution orchestration**：sample lock、immutable input fingerprint、受限制的 stage-aware resume、batch execution/summary、三案 regression。它**不會修改**已通過 2026-10-31 三案端對端驗證的 `pipeline.py`、`pipeline_v2.py`、`method2_score.py`、`qc_report.py` 或 vendor gVCF reconstruction。

- 原 v2 已實測：三案原始 gVCF reconstruction 3/3、Method 2 A/B 6/6；scorefile provenance 610,305/610,305；Model A 610,305、B 603,438；reference 3,330 unrelated／735 EAS；獨立 PLINK2 scoring 6/6。
- 本 add-on 在交付環境完成 **25 項 Python synthetic/mock tests**（含原 v2 的 12 項），但**尚未在使用者 Linux 伺服器真實 gVCF/Nextflow 執行 v2.1 的新入口**。所以現在是 v2.1 candidate／server acceptance pending，不要把已驗證的 v2 說成已驗證的 v2.1。
- 只為 *post-roundtrip* failure 提供真正的 Nextflow `-resume`。若 upstream gVCF reconstruction 或 chr-wise roundtrip 不完整，會明確拒絕續跑；使用**新 run ID** 重跑，**不會**刪除、覆蓋或清理舊輸出。
- `--work` 仍為舊工具/data root；vendor Nextflow 實際使用 `-w "$OUT/nextflow_work"` 且 `cd "$OUT"`，已確認各樣本隔離。

## 2. 安裝（新增檔案；不要覆蓋已驗證的 v2）

將 `PGS004879_integrated_v21_addon.zip` 上傳到伺服器任意暫存目錄。以 ZIP 在 `$HOME` 為例：

```bash
PROJECT=/home/n102968/2026_INTERVENE_PRS
PIPE="$PROJECT/pipeline/PGS004879_integrated_v2"
mkdir -p "$PROJECT/pipeline/v21_staging"
unzip -q ~/PGS004879_integrated_v21_addon.zip -d "$PROJECT/pipeline/v21_staging"
ADDON="$PROJECT/pipeline/v21_staging/PGS004879_integrated_v21_addon"

bash "$ADDON/install_addon.sh" "$PIPE"
cd "$PIPE"
python3 -m unittest discover -s tests -v
python3 pipeline_v21.py verify-baseline --config config.json
```

`install_addon.sh` 先用 SHA256 核對十個 v2 原檔，與此 add-on 所依據的**使用者先前上傳 v2 ZIP**一致，才新增 `pipeline_v21.py`、`v21_raw_resume.py`、測試、manifest 與 SOP。若 server 的 v2 程式不同會 **fail closed**，此時應先確認是否安裝了同一版，不能關閉 checksum guard。

**不需重新 `freeze-model`**：沿用原先 `config.json` 指定的 `method2_model_lock.json`；不要覆蓋既有 validated model 或跑 `nextflow clean`。

> 若日後更新 v2.1 add-on，應先封存原版並以新的版本／目錄發行。本 installer 遇到同名 add-on 檔案時會拒絕覆寫。

## 3. 第一個 server acceptance：不重新跑 gVCF

先重用已在原 v2 **fresh validated** 的三份 raw outputs，測試新的 v2.1 外層 orchestration、固定 Model A/B scoring 及 regression：

```bash
cd "$PIPE"
python3 pipeline_v21.py batch \
  --config config.json \
  --manifest manifest.reuse_3cases.example.tsv \
  --run-id method2_v21_reuse_validation_01 \
  --jobs 1

python3 pipeline_v21.py verify-pilot \
  --config config.json \
  --run-id method2_v21_reuse_validation_01
```

應看到 `V21 BATCH PASS: 3 samples` 與 `V21 PILOT REGRESSION: PASS (6/6)`。如 manifest 路徑與伺服器實際 v2 fresh outputs 不同，修改 manifest 的 `raw_source`，**不要修改原 v2 raw outputs**。這是 reuse-mode acceptance，**不等於**新入口的 raw gVCF end-to-end validation。

### 固定驗證 benchmark

| Sample | Model | SUM | EAS percentile | Observed | Missing |
|---|---|---:|---:|---:|---:|
| 26T00970 | A | 0.366424148909 | 35.1020408163 | 609032 | 1273 |
| 26T00970 | B | 0.363472617655 | 35.3741496599 | 602628 | 810 |
| 26T00969 | A | 0.474838562769 | 60.5442176871 | 609162 | 1143 |
| 26T00969 | B | 0.475804056111 | 61.4965986395 | 602768 | 670 |
| 26T00968 | A | 0.556804458045 | 79.4557823129 | 609272 | 1033 |
| 26T00968 | B | 0.552610418498 | 79.4557823129 | 602869 | 569 |

`verify-pilot` 檢查 6 組 score、percentile、observed/missing counts 及 v2.1 saved-output SHA256。`total_sum` tolerance ≤2e−6，percentile tolerance ≤1e−6 percentage points，counts 必須完全一致。

## 4. 第二個 server acceptance：真正 fresh raw gVCF

完成 reuse acceptance 後，若要驗證新入口從原始 gVCF 的完整 end-to-end execution：

```bash
cd "$PIPE"
python3 pipeline_v21.py batch \
  --config config.json \
  --manifest manifest.fresh_3cases.example.tsv \
  --run-id method2_v21_fresh_validation_01 \
  --jobs 1

python3 pipeline_v21.py verify-pilot \
  --config config.json \
  --run-id method2_v21_fresh_validation_01
```

上述命令**會重新執行 raw reconstruction、22-chromosome roundtrip、pgsc_calc 與 Method 2 scoring**，但輸出位於新的 run ID，原本 `pilot_fresh_02`／`method2_v2_fresh_validation_01` 均保持不變。三案 PASS 才能標記 v2.1 fresh end-to-end accepted。

`--jobs 1` 為預設與推薦起點；可選 `--jobs 2` 啟用兩案並行，**尚未完成 concurrency pressure testing**，要先確認機器 CPU、RAM、I/O 與 Nextflow 資源。最多接受 8。

## 5. 單一樣本與 resume

```bash
python3 pipeline_v21.py run --config config.json \
  --sample SAMPLE_ID \
  --run-id research_001 \
  --mode fresh \
  --gvcf /absolute/path/SAMPLE_ID.hard-filtered.gvcf.gz
```

當同一個 v2.1 run 執行失敗後，**請先檢查** `integrated_runs/research_001/SAMPLE_ID/orchestration/v21_status.json` 與 `logs/v21_v2_stage.log`。原始輸入與設定未更動時可嘗試：

```bash
python3 pipeline_v21.py run --config config.json \
  --sample SAMPLE_ID --run-id research_001 \
  --mode fresh --gvcf /absolute/path/SAMPLE_ID.hard-filtered.gvcf.gz \
  --resume
```

這條路徑的恢復策略：

1. **v2.1 input manifest 與 model lock 一律先核對。** 存在不同 gVCF metadata、QC parameters、scorefile、raw vendor script、v2/v2.1 code hash 或 reference metadata 時直接拒絕，不能假裝 resume。
2. **已完成 Method 2：** 需要原先 output SHA256 manifest 完全相同，驗證後原樣返回，不重新計分。若剛好在 v2 PASS 後、寫入 v2.1 SHA256 manifest 前中斷，明確指定 `--resume` 時會先在獨立暫存目錄重新計分，比對四份 deterministic outputs 完全相同後才補建立 checksum manifest；不覆寫原結果。
3. **raw 已完整、Method 2 未完成：** 原 v2 驗證 raw input/source 後重新做 Method 2，沿用已有 `inputs.lock.json`，不重新構建 genotype。
4. **raw partial，但完整 chr1–22 roundtrip 已完成：** 核對 gVCF 與原始 scorefile SHA256、vendor-generated verifier exact source、samplesheet、PGEN/VCF allele/dosage、`local_containers.config`；成功才執行原 verifier `reconcile`，必要時在**同一** `raw/nextflow_work` 上 `nextflow ... -resume`，再執行 `reconcile`。
5. **raw reconstruction／roundtrip 不完整：** **拒絕在原 output directory 恢復或覆寫**。保留失敗現場，使用新 run ID 重新 fresh。這是 v2.1 的明確限制，不是完整的全階段 checkpoint restart。

**禁止**為了通過 stage checks 手工編輯 `final_reconciliation.json`、`source_sha256.txt` 或 `inputs.lock.json`；也不要手動刪除 `raw/nextflow_work`、執行 `nextflow clean`、覆寫原模型或將 missing genotype 改為 0/0。

## 6. 批次執行及輸出

Manifest 是 **tab-delimited TSV**；欄位 `sample`, `mode`, `gvcf`, `raw_source`。每一個 sample 只可出現一次；`fresh` 提供 gVCF，`reuse` 提供已驗證 raw directory；混合模式亦可。

```bash
python3 pipeline_v21.py batch \
  --config config.json \
  --manifest YOUR_MANIFEST.tsv --run-id new_batch_run_01 --jobs 1

# 一次或部分執行失敗後，只在輸入不變的前提下：
python3 pipeline_v21.py batch \
  --config config.json \
  --manifest YOUR_MANIFEST.tsv --run-id new_batch_run_01 --jobs 1 --resume

# 不執行計分，只重新輸出執行狀態與結果彙總：
python3 pipeline_v21.py summarize \
  --config config.json --run-id new_batch_run_01
```

主要檔案：

```text
integrated_runs/<run-id>/
├── batch_summary_v21.tsv
├── pilot_regression_v21.tsv           # 僅 verify-pilot 後
└── <sample>/
    ├── orchestration/
    │   ├── v21_inputs.lock.json       # input fingerprint
    │   ├── v21_status.json            # STARTED/RAW_REPAIR/V2_RUNNING/PASS/FAILED
    │   └── v21_output_sha256.json     # completed output checksums
    ├── logs/v21_v2_stage.log
    ├── raw/                           # only fresh mode
    │   ├── nextflow_work/             # isolated Nextflow -resume cache
    │   └── logs/v21_nextflow_resume.log  # only if resumed
    └── method2/
        ├── summary.json
        ├── scores.tsv
        ├── missing_loci.tsv.gz
        ├── report.md
        ├── inputs.lock.json
        └── status.json
```

`batch_summary_v21.tsv` 輸出 A/B 的 `observed_loci`、`missing_loci`、SUM、EAS Z/percentile、分組和 status。**即使 batch 中某案失敗，其他樣本仍會完成並產生 summary；整個 batch 會回傳非零 exit code。**

## 7. 風險與研究限制

- 新版只改 execution layer，**不提升或改變** CHD prediction performance、PRS absolute-risk model 或 ancestry calibration。原 EAS 百分位數仍是 reference comparison，而非 CHD cumulative incidence。
- Sample locking：v2.1 另使用 `integrated_runs/.v21_locks/<sample>.lock`，並在原 v2 核心執行時沿用 v1 的 sample lock；多程序協作前仍須在伺服器驗證 lock contention、失敗恢復與 I/O contention。
- Input identity 包含小型 config/script SHA256、模型 lock、重要大型資源檔案的 path/size/mtime，以及 raw stage 已完成後 gVCF/source scorefile 的 SHA256。未完成階段的資源是否曾被外部替換，仍需依輸入檔案管理 SOP 防範。
- 欲升級大型批次正式 production 前，還需要 failure injection、parallel resource benchmark、audit schema versioning、job scheduling 和站點備份／回復策略。
- `--resume` **不能**從 upstream 未完成的 gVCF reconstruction 自動續跑；不能把此版本宣稱為全階段 resume 功能。

## 8. Golden v2 不可變更內容

部署前核對 `BASELINE_V2_SHA256SUMS.txt` 內的 10 個 validated v2 component；`config.json`、`method2_model_lock.json`、官方 scorefile 及 Method 2 reference files 沿用原伺服器，原路徑不改。若 check 失敗，停止部署，不要去除 checksum 驗證。
