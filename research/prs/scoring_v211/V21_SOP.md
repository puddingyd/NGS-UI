# PGS004879 integrated v2.1.1｜主操作手冊：gVCF → fixed Model A/B

**更新：2026-11-02。** 研究用途；非臨床 CHD absolute-risk 模型。此文件取代 2026-10-31 的 v2.1 初始 SOP；v2.0 科學演算法與 locked reference 未變。

## 1. 系統範圍與版本

- **v2 golden**：`pipeline.py`（v1 raw/ancestry）、`pipeline_v2.py`（固定 Method 2 入口）、`method2_score.py`、vendor 原始 gVCF reconstruction；`BASELINE_V2_SHA256SUMS.txt` 十項 **PASS**。
- **v2.1**：`pipeline_v21.py`、`v21_raw_resume.py`（sample `flock`、immutable fingerprints、stage-aware resume、batch/summarize、verify-pilot）。
- **v2.1.1 hotfix**：`batch()` 的 `expected={r.sample for r in jobs}`，修復原 `Namespace` indexing 錯誤；新增 `tests/test_v21_batch_success.py` 兩項回歸測試，**27/27 PASS**。
- **Model A/B 不變**：A 610,305、B 603,438；unrelated 3,330 的 `2*p_effect` 補缺值；EAS 735 的 fixed ECDF percentile。不得將 ECDF percentile 當成台灣 CHD absolute risk。

全套伺服器驗證紀錄見 [`VALIDATION_REPORT_v211_2026-11-02.md`](VALIDATION_REPORT_v211_2026-11-02.md)。目前已通過 3-case fresh、兩案 parallel fresh、sample lock、completed resume、post-roundtrip recovery；**mid-workflow partial-task interruption recovery 尚未測試**。

## 2. 重要路徑與啟動前檢查

```bash
PROJECT=/home/n102968/2026_INTERVENE_PRS
PIPE="$PROJECT/pipeline/PGS004879_integrated_v2"
cd "$PIPE"

python3 -m json.tool config.json >/dev/null
python3 -m unittest discover -s tests -v      # 27/27 PASS
python3 pipeline_v21.py verify-baseline --config config.json
# V2 GOLDEN PROGRAM BASELINE: PASS (10 files)
# METHOD 2 MODEL LOCK: PASS
```

伺服器 `config.json` 已將 `raw_work` 改為 `$PROJECT`，vendor command 也實測 `--work /home/n102968/2026_INTERVENE_PRS`；每個樣本獨立 `integrated_runs/<run-id>/<sample>/raw/nextflow_work`。`config.example.json` 同步修正，**正式設定不要再由 example 直接覆蓋**。

**注意發行形式：** 使用者的 server ZIP 是已整合 add-on 的完整程式目錄，無須重新執行過去 v2.1 的 `install_addon.sh`。新部署請於乾淨目錄解壓並確認所有 absolute paths、模型資料、containers 均存在；不要直接在原驗證目錄覆蓋程式。必要時先在 staging 實測。

## 3. Manifest、fresh／reuse 及批次入口

Manifest 必須是 TAB-separated，包含 `sample`, `mode`, `gvcf`, `raw_source`。每 sample 唯一；`fresh` 指定 gVCF，`reuse` 指定已通過 raw audit 的 `raw_source`。示例檔：

- `manifest.v211_fresh_3cases.tsv`：三案原始 gVCF；**會重算**。
- `manifest.reuse_3cases.example.tsv`：三案 reuse 的已驗證 raw audit；**不重算 raw**。
- `manifest.v211_parallel_2cases.tsv`：兩案 fresh（26T00970、26T00969），用於已完成的 parallel validation；不要重用既有 run ID 進行新批次。

```bash
# 新批次：請使用新的 run ID 和自己的 TSV，不能沿用 validated run ID
python3 pipeline_v21.py batch --config config.json \
  --manifest /absolute/path/NEW_BATCH.tsv \
  --run-id NEW_BATCH_ID --jobs 1

# 二案並行已驗證 correctness，但非 throughput/8-way 正式確效
python3 pipeline_v21.py batch --config config.json \
  --manifest /absolute/path/NEW_BATCH_2.tsv \
  --run-id NEW_PARALLEL_RUN_ID --jobs 2

# 單一樣本
python3 pipeline_v21.py run --config config.json \
  --sample SAMPLE_ID --run-id NEW_SINGLE_RUN_ID \
  --mode fresh --gvcf /absolute/path/SAMPLE_ID.hard-filtered.gvcf.gz
```

`--jobs 1` 預設；允許 1–8，**只有 1 和 2 有真實 sample-level 驗證**。大量資料前評估 CPU、RAM、I/O、NXF executors、container cache 與 disk quota。

## 4. 已驗證的 acceptance runs

| Stage | Run ID | 記錄結果 |
|---|---|---|
| V2 tool-root regression | `method2_v2_toolroot_check_01` | 6/6 PASS |
| V2.1.1 reuse | `method2_v21_reuse_validation_02` | batch 3/3、regression 6/6 PASS |
| V2.1.1 fresh | `method2_v211_fresh_validation_01` | gVCF 3/3、regression 6/6 PASS |
| Fresh run already-complete resume | 同上一行 | `V21_ALREADY_COMPLETE` 3/3；9/9 method2 checksum identical |
| Cross-process locking | isolated temp dirs | same-sample blocked、different sample allowed、reacquire PASS |
| Nextflow start-failure recovery | `method2_v211_recovery_test_01` | upstream unchanged，A/B 2/2 PASS |
| Completed Nextflow cache | 同 recovery run | cached tasks > 0（最大顯示 22）；3/3 SHA256 unchanged |
| **Parallel two-sample fresh** | **`method2_v211_parallel_fresh_01`** | **two starts at 2026-11-02 02:54:10 UTC；2/2 case PASS、4/4 benchmark PASS** |
| Mid-workflow task interruption | — | **NOT TESTED** |

### 4.1 Fixed golden A/B benchmarks

| Sample | A SUM | A EAS pctl | B SUM | B EAS pctl |
|---|---:|---:|---:|---:|
| 26T00970 | 0.366424148908984 | 35.1020408163 | 0.363472617655232 | 35.3741496599 |
| 26T00969 | 0.474838562769280 | 60.5442176871 | 0.475804056110600 | 61.4965986395 |
| 26T00968 | 0.556804458044592 | 79.4557823129 | 0.552610418497840 | 79.4557823129 |

`verify-pilot` 只適用含三個已知樣本的 run，檢查六組 SUM、percentile、observed/missing loci、output SHA256，不適合只含 2 個 samples 的並行 run。二案並行的 **4/4** 實測、expected observed/missing 及 log 證據見 [`VALIDATION_REPORT_v211_2026-11-02.md`](VALIDATION_REPORT_v211_2026-11-02.md)。

## 5. `--resume`：範圍、必要條件與拒絕邏輯

```bash
# 和首次相同 manifest / run ID / input / code / model
python3 pipeline_v21.py batch --config config.json \
  --manifest /absolute/path/NEW_BATCH.tsv \
  --run-id NEW_BATCH_ID --jobs 1 --resume
```

1. 已存在的 v2.1 run **必須**明確 `--resume`，並要求 `orchestration/v21_inputs.lock.json` 與目前 fingerprint 完全一致；包括 gVCF/source metadata、模型 lock、程式 hash、參數。程式更新（例如 v2.1→v2.1.1）後應使用**新 run ID**。
2. Completed-run：檢查 `v21_output_sha256.json` 內容對照當下 outputs，若結果未變回傳 `V21_ALREADY_COMPLETE`，**不重算**；完整 fresh validation 三案已實測 9 個核心檔案 SHA256 不變。
3. 原 raw 已完成但 Method 2 未完成：核對 `raw_integrity` 與 source fingerprint，然後執行 v2/Method 2，不改 raw genotype。
4. 原 raw 未完成但 chr1–22 roundtrip 完整：`v21_raw_resume.py` 逐一核對 gVCF/scorefile 來源、generated verifier、samplesheet 22chr、PGEN/VCF GT roundtrip、local Nextflow config，才以同一 `raw/nextflow_work` 執行 `nextflow ... -resume`，再檢查 reconciliation。
5. 原 raw reconstruction、GT QC 或 chr-wise roundtrip 尚未完成：**fail closed**；保留現場，使用全新 run ID fresh，不會對原資料做不安全的 in-place repair。
6. 若已有 saved output SHA256 而現在輸出變更，**拒絕覆寫**；不能為了通過檢查自行編輯 `v21_inputs.lock.json`、model lock、raw audit 或 summary。

Post-roundtrip start failure 已實測：在第一個 `nextflow run` 前刻意注入 error，`--resume` 後 raw / Method 2 正確，但這次第一次真正呼叫 Nextflow 時 `-resume` 被忽略是因 **沒有 session cache**；後續 completed-task cache probe 才真的看到 cached tasks。**這兩者都不能取代 partial-task failure-recovery test**。

## 6. 輸出、summary 與唯讀查看

```text
$PROJECT/integrated_runs/<run-id>/
├── batch_summary_v21.tsv
├── pilot_regression_v21.tsv              # only after verify-pilot
└── <sample>/
    ├── orchestration/
    │   ├── v21_inputs.lock.json
    │   ├── v21_status.json
    │   └── v21_output_sha256.json
    ├── logs/v21_v2_stage.log
    ├── raw/                              # only fresh mode
    │   ├── audit/final_reconciliation.json
    │   ├── roundtrip_22chr/
    │   ├── nextflow_work/
    │   └── logs/v21_nextflow_resume.log  # only if resumed via repair
    └── method2/
        ├── scores.tsv
        ├── summary.json
        ├── missing_loci.tsv.gz
        ├── inputs.lock.json
        ├── status.json
        └── report.md
```

`batch_summary_v21.tsv` 列 Model A/B 個別 score、missing/imputed、EAS Z/percentile、group、PASS/FAILED 與 model lock SHA。某個 sample 失敗後，其他 sample 仍會完成；整體 batch 回傳非零 code。

**唯讀查驗**（避免更新已封存檔案）：

```bash
PROJECT=/home/n102968/2026_INTERVENE_PRS
RUN="$PROJECT/integrated_runs/method2_v211_parallel_fresh_01"
column -s $'\t' -t "$RUN/batch_summary_v21.tsv"
sha256sum "$RUN"/{26T00970,26T00969}/method2/scores.tsv
```

`pipeline_v21.py summarize` **會重新寫出** `batch_summary_v21.tsv`；`verify-pilot` **會重新寫出** `pilot_regression_v21.tsv`。在 sealed evidence run 上，若僅需查看，優先直接讀現有檔案；不要不必要地重新生成。

## 7. Debug／常見錯誤

| 情況 | 解釋與安全做法 |
|---|---|
| `Existing v2.1 run requires explicit --resume` | 已存在 fingerprint；核對原 manifest 與輸入後，只有同版同 input 才能 `--resume` |
| `Immutable input fingerprint changed` | 程式／模型／資料／設定改變；**新 run ID**，不繞過 |
| `Namespace object is not subscriptable` | 原 v2.1 batch success-path bug，**已於 v2.1.1 hotfix 修正**；確認 `pipeline_v21.py` SHA 和 27 tests |
| `Cannot resume pre-roundtrip` | upstream/roundtrip checkpoint 未完整，不能 in-place resume；另開 run ID |
| `Previously completed outputs now invalid` | 已完成結果與 checksum 不符，停止稽核，勿直接重算覆蓋 |
| Nextflow `-resume` ignored | 若是第一次真正啟動 Nextflow，此警告正常；不表示完成 cache-reuse 驗證 |
| Nextflow pulling `pyyaml` image | 反映 container cache 未預先固定；大量批次前設定共享 `NXF_SINGULARITY_CACHEDIR`、預先驗證 images |
| `--jobs 2` 成功 | **correctness verified**，不等於確認 3..8 也正確或資源線性擴展 |

## 8. 未完成工作與研究限制

- **最近待辦：** controlled mid-workflow *partial-task* interruption / cached-task reuse evidence；另做 5–10 獨立 WGS 的資源和吞吐量測試。
- 儲存 output SHA256、source/model lock、Nextflow execution metadata 的受控封存與 immutable validation ledger。
- 獨立 target ancestry assignment、台灣族群 PRS calibration、與 INTERVENE CHD absolute-risk 轉換／outcome validation 分開進行；ECDF EAS percentile 不是醫療診斷與 CHD 事件機率。
- 後續正式規劃、可行的驗收條件見 [`NEXT_STEPS_v211.md`](NEXT_STEPS_v211.md)。
