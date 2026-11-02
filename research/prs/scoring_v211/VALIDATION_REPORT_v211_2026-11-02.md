# PGS004879 v2.1.1｜Server validation evidence（截至 2026-11-02）

**狀態：端對端科學重現、completed-run resume、sample lock、post-roundtrip orchestration recovery、completed-task cache reuse 與 2-sample parallel correctness 已通過。** 未進行 mid-workflow *partial task interruption* recovery。以下記錄根據使用者交付的伺服器輸出；**本 ZIP 沒有包含大型 `integrated_runs/` output 與 Nextflow `.nextflow` cache**，需在原伺服器留存作原始稽核依據。

## 1. 程式與模型基準

- Golden v2：`BASELINE_V2_SHA256SUMS.txt` **10 files PASS**；`pipeline.py`、`pipeline_v2.py`、`method2_score.py`、兩份 vendor script、v2 unit tests 等未變。
- Model 2：`METHOD 2 MODEL LOCK: PASS`；reference unrelated 3,330 / EAS 735；Model A 610,305 loci，Model B 603,438 loci，missing genotype 以 `2*p_effect_unrelated3330` 插補。
- 伺服器 Python unit tests：**27/27 PASS**（v2 12 + v2.1 13 + v2.1.1 2）；其中 batch success／error path 追加測試。
- Model lock 檔案 `method2_model_lock.json` 自身 SHA256：`a3b2bb7a52762d3232cd7e72040b26470c882e5afce66d8e574edc9faeda1083`（使用者實際 `batch_summary_v21.tsv` 報告）。
- 官方 harmonized scorefile 解壓 SHA256：`0ee9852635bc831bfbcbb144cb308157bed658e0e26fc4243737b3957aea4d3f`；Model A/B scorefile SHA256 分別 `358dc695afaee2387b14da0641f23a46b8d7cf33312f6ee14225eb9018d37f96` / `980541e5f1df2d5c5cd439f1f7020b68ce40a0af403bd6a6e44c02a62e675a42`。

## 2. 驗證矩陣（run ID + 判準 + 實際結果）

| Validation | Run / evidence | 實際結果 | 邊界 |
|---|---|---|---|
| Model / source lock | `verify-baseline --config config.json` | **10/10 + model lock PASS** | 不等於台灣 CHD outcome calibration |
| Tool-root 修正後 v2 regression | `method2_v2_toolroot_check_01` | **6/6 PASS** | 只使用既有 raw audits |
| v2.1.1 reuse | `method2_v21_reuse_validation_02` | **3/3; 6/6 PASS** | 未重新建立 raw gVCF |
| v2.1.1 fresh | `method2_v211_fresh_validation_01` | **3/3; 6/6 PASS** | 三案原始 raw gVCF 全流程 |
| Completed-run `--resume` | 同 `method2_v211_fresh_validation_01` | **3 `V21_ALREADY_COMPLETE`; 9/9 SHA256 identical** | 不包含不完整 task-restart 情境 |
| OS-process `flock` | 獨立 temporary test directory | **same-sample blocked; different-sample allowed; release PASS** | 非整套多樣本 workflow stress test |
| Post-roundtrip fault injection | `method2_v211_recovery_test_01` | **expected fail → `--resume` PASS; A/B 2/2; upstream SHA unchanged** | 故障在 Nextflow 啟動前，未產生先前完成的 task cache |
| Nextflow completed-run cache probe | 同 recovery test run | **cached tasks detected; 3/3 file SHA unchanged** | log `Maximum displayed cached tasks: 22` 非總 cached tasks |
| **2-case parallel fresh** | **`method2_v211_parallel_fresh_01`** | **2/2 samples; 4/4 A/B; batch PASS** | 驗證 `--jobs 2`，未驗證更高 concurrency |
| Mid-workflow task recovery | 尚無 run | **NOT TESTED** | 正式進一步 gate |

## 3. Golden Method 2 結果（全部由三案 fresh 6/6 重現）

| Sample | Model | n_loci | observed | missing | observed_sum | imputed_sum | total_sum | EAS Z | EAS percentile | Group |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 26T00970 | A | 610305 | 609032 | 1273 | 0.354661 | +0.011763148909 | **0.366424148909** | -0.346354 | **35.102041** | 20-40 |
| 26T00970 | B | 603438 | 602628 | 810 | 0.351645 | +0.011827617655 | **0.363472617655** | -0.333872 | **35.374150** | 20-40 |
| 26T00969 | A | 610305 | 609162 | 1143 | 0.475192 | -0.000353437231 | **0.474838562769** | +0.341565 | **60.544218** | 60-80 |
| 26T00969 | B | 603438 | 602768 | 670 | 0.476581 | -0.000776943889 | **0.475804056111** | +0.381683 | **61.496599** | 60-80 |
| 26T00968 | A | 610305 | 609272 | 1033 | 0.558069 | -0.001264541955 | **0.556804458045** | +0.861661 | **79.455782** | 60-80 |
| 26T00968 | B | 603438 | 602869 | 569 | 0.552992 | -0.000381581502 | **0.552610418498** | +0.870941 | **79.455782** | 60-80 |

其中 `26T00968` EAS percentile 接近 80th percentile，仍應列為分組邊界審查樣本；這是研究分組穩定性問題，並非 CHD clinical cutoff。

## 4. 2026-11-02：兩案平行 fresh validation（本次重點）

### 4.1 實際指令

```bash
cd /home/n102968/2026_INTERVENE_PRS/pipeline/PGS004879_integrated_v2
python3 pipeline_v21.py batch \
  --config config.json \
  --manifest manifest.v211_parallel_2cases.tsv \
  --run-id method2_v211_parallel_fresh_01 --jobs 2
```

### 4.2 並行證據

- `26T00970` 與 `26T00969` 兩條 vendor raw script 日誌都顯示 **`START 2026-11-02T02:54:10+00:00`**，代表同一 UTC 秒開始；兩者的 `chr` 處理、QC、pgsc_calc 及後續 stage 有交錯輸出。這比單純傳入 `--jobs 2` 更能支持真正的時間重疊。
- 兩者各自使用 `.../integrated_runs/method2_v211_parallel_fresh_01/<sample>/raw`，並從 vendor script 取得獨立的 `nextflow_work`；`--work` 實測為專案根目錄 `/home/n102968/2026_INTERVENE_PRS`，**不再指向 26T00970**。
- `26T00969`: raw matched **609,488**、dbSNP-confirmed flips **331**、unmatched **1,189**、22 chr round-trip **PASS**、raw `SUM=0.471545`，`final_reconciliation SUM error=1E-16`。
- `26T00970`: raw matched **609,346**、dbSNP-confirmed flips **330**、unmatched **1,331**、22 chr round-trip **PASS**、raw `SUM=0.351907`，`final_reconciliation SUM error=1E-16`。
- 兩者均 `METHOD2 SCORING: TECHNICAL PASS`、`V21 SAMPLE COMPLETE`，最後顯示 **`V21 BATCH PASS: 2 samples`**。

### 4.3 獨立比對四組 golden benchmarks

| Sample | Model | Fresh SUM | EAS percentile | Observed / Missing | 比對結果 |
|---|---|---:|---:|---|---|
| 26T00970 | A | 0.366424148909 | 35.102041 | 609032 / 1273 | PASS |
| 26T00970 | B | 0.363472617655 | 35.374150 | 602628 / 810 | PASS |
| 26T00969 | A | 0.474838562769 | 60.544218 | 609162 / 1143 | PASS |
| 26T00969 | B | 0.475804056111 | 61.496599 | 602768 / 670 | PASS |

使用者執行獨立 `csv.DictReader` assert，要求 `execution_status/raw_status/method2_status` 全部 PASS，`total_sum` 絕對誤差 `< 1e-9`，EAS percentile `< 1e-6`，observed/missing loci **完全相等**，實際印出 **`PARALLEL FRESH PRS REGRESSION: PASS (4/4)`**。這是**平行工作下的數值一致性**，不等於已量測可擴展性、效能加速比或資料中心級資源隔離。

### 4.4 本次輸出位置

```text
$PROJECT/integrated_runs/method2_v211_parallel_fresh_01/
├── batch_summary_v21.tsv
├── 26T00970/
│   ├── raw/{audit/, roundtrip_22chr/, nextflow_work/, results/}
│   ├── orchestration/{v21_inputs.lock.json, v21_status.json, v21_output_sha256.json}
│   └── method2/{scores.tsv, summary.json, report.md, missing_loci.tsv.gz, ...}
└── 26T00969/
    └── ... 同樣獨立的 raw/orchestration/method2 結構
```

## 5. Recovery 訊息的證據邊界

`method2_v211_recovery_test_01` 在 upstream 22-chromosome round-trip PASS 後，第一次呼叫 `nextflow run` 被 PATH shim 主動回傳 error 86。稍後正常 `--resume`：

1. 重新獨立核對 `roundtrip`；初次 `reconcile` 因尚無 aggregated scores 失敗（預期現象）。
2. 進入 `RAW_PGSC_CALC_RESUME`、同一個 `raw/nextflow_work`；Nextflow 顯示 **`WARN ... Option -resume is ignored`**，因第一次真正啟動時沒有先前 Nextflow session，**此警告在該測試設計下是合理的**。
3. 最終 `RAW POST-ROUNDTRIP RESUME: PASS`、Model A/B **2/2** 重現 golden；三個 upstream checkpoint SHA256 完全不變。
4. 後續對**已完成**的 Nextflow run 再做一次 `-resume`，log 顯示至少一個 process 的 **22 cached tasks**，且 VCF/reconciliation/aggregated score 三份檔案 SHA256 unchanged，故 completed-task cache reuse **PASS**。

**不能宣稱**「中斷時已完成部分 Nextflow tasks，而第二次使用它們繼續工作」已實測。這需要專門的 mid-workflow failure injection 與 cached-task audit；參閱 [`NEXT_STEPS_v211.md`](NEXT_STEPS_v211.md)。

## 6. 建議的伺服器端歸檔命令（唯讀，不覆寫）

```bash
PROJECT=/home/n102968/2026_INTERVENE_PRS
PIPE="$PROJECT/pipeline/PGS004879_integrated_v2"
cd "$PIPE"
python3 pipeline_v21.py verify-baseline --config config.json
python3 -m unittest discover -s tests -q
python3 pipeline_v21.py verify-pilot --config config.json \
  --run-id method2_v211_fresh_validation_01
python3 pipeline_v21.py summarize --config config.json \
  --run-id method2_v211_parallel_fresh_01
```

**注意**：最後 `summarize` 會重寫該 run 的 `batch_summary_v21.tsv`；若要完全唯讀歸檔，只需 `cat` 或 `sha256sum` 既有檔案，不要執行 `summarize`。`verify-pilot` 也會更新 `pilot_regression_v21.tsv`；封存後不應重跑這兩項而期待檔案 metadata 不變。
