# PGS004879 integrated v2.1.1｜使用說明與驗證索引

**文件更新：2026-11-02；來源：使用者伺服器封存 ZIP 與已提供的真實執行紀錄。** 研究用途（research use only）。本版是 **v2.0 不變的科學計分核心 + v2.1 orchestration + v2.1.1 batch hotfix**，不是新 PGS 模型版本。

## 本版已確認的狀態

| 項目 | 伺服器驗證狀態 |
|---|---|
| v2 原始程式／vendor baseline | `verify-baseline`: **10/10 PASS** |
| Model A/B fixed resources | `METHOD 2 MODEL LOCK: PASS` |
| Unit tests | **27/27 PASS**（v2 原有 12 + v2.1 原有 13 + v2.1.1 新增 2） |
| 3-case reuse batch | **3/3 PASS**；A/B **6/6 PASS** |
| 3-case fresh gVCF → PRS | **3/3 PASS**；A/B **6/6 PASS** |
| 已完成 run 的 `--resume` | **3/3 已完成直接沿用**；9 個核心檔案 SHA256 完全一致 |
| 跨程序相同 sample lock / 不同 sample lock | **PASS**（實際 OS-process test） |
| Post-roundtrip 故障注入後恢復 | **PASS**（Nextflow 啟動前注入失敗，復原後 A/B **2/2 PASS**） |
| Completed Nextflow task-cache reuse | **PASS**（log 出現 cached tasks，最大顯示數 22；不是整個 workflow 的 task 總數） |
| **2-case concurrent fresh** (`--jobs 2`) | **2/2 PASS**；A/B **4/4 PASS**；兩條 raw job 實際時間重疊 |
| Mid-workflow *partial task* interruption → resume | **尚未實測** |
| 大規模 throughput / clinical CHD outcome calibration | **尚未實測／未完成** |

**重要證據界線：** 本套件包含目前伺服器上的程式、設定與測試腳本，**不包含** `integrated_runs/` 的真實 WGS 產出、Nextflow caches、固定 Model A/B 的大型 reference 資源；上表的 server PASS 係依使用者在對話中提供的輸出紀錄歸納。完整 run ID 與 golden 數字參見 [`VALIDATION_REPORT_v211_2026-11-02.md`](VALIDATION_REPORT_v211_2026-11-02.md)。

## 建議閱讀順序

1. [`V21_SOP.md`](V21_SOP.md)：**目前 v2.1.1 主操作手冊**，安裝、`batch`、`--resume`、故障恢復、平行執行、異常處理。
2. [`VALIDATION_REPORT_v211_2026-11-02.md`](VALIDATION_REPORT_v211_2026-11-02.md)：2026-11-02 驗證證據、Model A/B 6 組 benchmark、並行測試 4/4。
3. [`NEXT_STEPS_v211.md`](NEXT_STEPS_v211.md)：mid-workflow task recovery、可觀測性、container cache、5–10 樣本負載試驗、正式 release gate。
4. [`METHOD2_INTEGRATION_SOP.md`](METHOD2_INTEGRATION_SOP.md)：固定 Model A/B 計分方法、缺失值 2p、EAS 735 ECDF、QC、source-to-canonical crosswalk。
5. [`RELEASE_NOTES_v211.md`](RELEASE_NOTES_v211.md)：程式與文件版本差異；v2.0 的 [`RELEASE_NOTES_v2.md`](RELEASE_NOTES_v2.md) 保留為歷史紀錄。

## 現行伺服器路徑與第一輪檢查

```bash
PROJECT=/home/n102968/2026_INTERVENE_PRS
PIPE="$PROJECT/pipeline/PGS004879_integrated_v2"
cd "$PIPE"

python3 -m json.tool config.json >/dev/null
python3 -m unittest discover -s tests -v   # 27/27
python3 pipeline_v21.py verify-baseline --config config.json
# V2 GOLDEN PROGRAM BASELINE: PASS (10 files)
# METHOD 2 MODEL LOCK: PASS
```

`config.json` 實際 `raw_work` 現為 `$PROJECT`，不再使用 `$PROJECT/26T00970`；每個樣本的 raw / Nextflow work 依 `--run-id` 與 sample 隔離。`config.example.json` 已與此一致，**不需要覆蓋已完成驗證的實際 `config.json`**。不要重新 `freeze-model` 或覆寫現有 model lock。

## 正式入口（新樣本）

單樣本：

```bash
python3 pipeline_v21.py run --config config.json \
  --sample SAMPLE_ID --run-id production_batch_001 \
  --mode fresh --gvcf /absolute/path/SAMPLE_ID.hard-filtered.gvcf.gz
```

批次 TSV 格式：`sample<TAB>mode<TAB>gvcf<TAB>raw_source`；`fresh` 用 gVCF，`reuse` 用完整 `raw_source`：

```bash
python3 pipeline_v21.py batch --config config.json \
  --manifest /absolute/path/batch_001.tsv \
  --run-id production_batch_001 --jobs 1

# 失敗後僅在所有原始輸入、程式和 model fingerprint 不變時續跑：
python3 pipeline_v21.py batch --config config.json \
  --manifest /absolute/path/batch_001.tsv \
  --run-id production_batch_001 --jobs 1 --resume

python3 pipeline_v21.py summarize --config config.json \
  --run-id production_batch_001
```

新 run ID 不可事先存在。`--jobs 2` 已通過兩案正確性驗證，但大規模部署前仍須資源負載測試；參數允許 1–8，**不是**宣稱 8-way 具備驗證。對已封存 run 做回歸時：

```bash
python3 pipeline_v21.py verify-pilot --config config.json \
  --run-id method2_v211_fresh_validation_01
```

## 不可混淆的事項

- Raw PGS SUM 與固定 Model A/B SUM 的 loci / imputation 方法不同，不應互相取代。
- EAS percentile 是 **固定 HGDP+1kGP EAS735 reference 的 ECDF 百分位**，並非 CHD 絕對風險。預設 raw-only **未**對新樣本確認 ancestry；必須有獨立 ancestry QC 才能做匹配解讀。
- 目前只支援 **post-roundtrip 的安全 resume**：upstream 重建或 22-chromosome roundtrip 尚未完整通過時拒絕在原 run 接續。
- `pipeline_v21.py` 的 batch hotfix 已改變程式 SHA256；**不要用此版對 hotfix 前的 run 強行 `--resume`**，應建立新 run ID。
- `SHA256SUMS_v2_ORIGINAL.txt` 是舊 v2.0 的歷史清單；當前文件更新後請使用此版 `SHA256SUMS.txt` 檢查檔案。`BASELINE_V2_SHA256SUMS.txt` 仍專門保護十個**未更動的 v2 golden 檔案**。
