# PGS004879 v2.1.1｜後續工作與建議驗收條件

**起點（2026-11-02）：** 已通過 27/27 unit tests、三案 fresh 6/6、completed-run resume 9/9 SHA256、cross-process locks、post-roundtrip start-failure recovery、completed-task cache reuse、兩案 parallel fresh 4/4。**本文件的後續步驟是待辦，不代表已完成。**

## P0：Mid-workflow partial-task interruption → resume（尚待實測）

**要驗證的不同情境：** real Nextflow 已實際完成部分 `pgsc_calc` tasks、其他 tasks 尚未完成時，指定的 failure injection 使工作失敗；第二次 `--resume` 正確重用之前已完成 tasks，且最終 PRS 與 golden 一致。

建議方法與安全要求：

1. 使用**全新** `run-id`，挑選已知樣本 26T00970；保持原 v2/v2.1.1 source / model lock 不變。測試目錄只供本次 injection，**不可使用** `method2_v211_fresh_validation_01` 或原 golden runs。
2. 事前備份 `run-id`、gVCF identity、程式 SHA、原始 model SHA、Nextflow 版本／config、磁碟空間；預設 `--jobs 1`。規劃**可控且可回滾**的 injection point（例如被測試隔離環境中受控失敗的後段 task），避免任意 `kill -9` 整個 process group、刪掉 `.nextflow/cache` 或修改 validated vendor source。
3. 在故障前確認 **completed tasks > 0 且 remaining tasks > 0**，保存 `.nextflow.log`/`trace`/工作目錄證據；故障後記錄非零 exit、下一階段未完成、`v21_status.json` FAILED、原始 raw roundtrip SHA。
4. 以**同一 run ID / 同一 gVCF / 相同程式與模型**執行 `pipeline_v21.py ... --resume`；復原日誌應有 `RAW_PGSC_CALC_RESUME`、**同一** `raw/nextflow_work`、並出現真正 cached tasks（不是只有接受 `-resume` 字串）。
5. 以 `method2/scores.tsv` 對照 26T00970 golden：A SUM `0.366424148908984`、B SUM `0.363472617655232`，A/B percentile `35.1020408163` / `35.3741496599`，observed/missing A=`609032/1273`、B=`602628/810`；確認 upstream input/roundtrip hash 未變。
6. 記錄本次的 process-level cache hit vs rerun 數目、failure cause、恢復時間、SHA256 manifest，**避免只憑最後結果 PASS 判定 cache recovery PASS**。

**通過標準：** A. 故障前存在已完成的 Nextflow tasks；B. 同一 run 的恢復顯示那些 tasks 被 cached；C. upstream GT/VCF/PGEN unchanged；D. final raw reconciliation 及 Model A/B 2/2 PASS。缺任何一項，標記 PARTIAL／NOT VERIFIED，不得寫作 PASS。

## P1：封存 v2.1.1 操作版本／證據

1. 原程式 SHA：`BASELINE_V2_SHA256SUMS.txt` 十項 PASS；v2.1.1 的 `pipeline_v21.py` 與 `v21_raw_resume.py` 另列 SHA256；模組與 model lock 作 immutable snapshot。
2. 將已完成 run 的 `batch_summary_v21.tsv`、`pilot_regression_v21.tsv`、重要 raw audit、`v21_inputs.lock.json`、`v21_output_sha256.json`、Nextflow `trace`/`.nextflow.log` 與 container image digests 歸入**院內的受控檔案庫**；ZIP 本身不含這些正式結果。
3. 建議維持 append-only validation ledger；新的 pipeline bugfix 要升版及更換新 run ID。**不可**以新程式對 input fingerprint 不符的舊 run 強行 resume。
4. 上線前固定 pgsc_calc `v2.3.0`、Nextflow `25.04.3`（此次 log 所用）、bcftools `1.19`、samtools `1.19.2`、PLINK `2.00a5.10`、Singularity/Apptainer 資料來源與 SHA，建立可重建環境。

## P2：5–10 樣本的資源負載／可靠性測試

- 先建立 5–10 個**不同樣本**（不能把同一 ID 複製當成獨立樣本並行），`--jobs 1` 作 correctness baseline，再在充分資源下測 `--jobs 2`；只有在 RAM、I/O、CPU 和 task executor 真實可承受時才考慮 3–8。CLI 上限 8 **不是**已驗證的最佳值。
- 每個 sample 記錄 wall-clock、CPU/RAM peak、Nextflow scheduled/peak tasks、disk writes/reads、`raw/nextflow_work` 高峰容量、總輸出容量、container cache 重複拉取、失敗率及重試次數。
- 明確為 `NXF_SINGULARITY_CACHEDIR` 設定可共享且權限正確的 image cache，預先準備並驗證本次 Nextflow log 顯示曾拉取的 `pgscatalog/pyyaml:6.0-singularity`，再評估離線執行；**不要**在 production run 中途改變 image config 或 task cache 路徑。
- 定義受限執行策略：全站 job queue、最多同時 Nextflow workflows、per-sample disk quota、失敗後不可自動強制清 cache、成功後受控分級保留 genotype/audit/report。
- 若結果與 golden 不符，停止擴量，審查 input identity / reference / missingness、sample lock 與資源競爭，而不是調低 QC 門檻。

## P3：研究方法與臨床外推（與 pipeline correctness 分離）

- 為新 WGS 建立**獨立的 ancestry inference**；EAS 735 ECDF 比較應附上 ancestry 一致性與 uncertainty，不能單憑 EAS percentile 宣稱 target 是 EAS。
- 以台灣 WGS population calibration 檢查 reference transferability、imputation/PRS rank/percentile stability；針對 26T00968（~79.46 percentile）等邊界個案做 sensitivity check。
- 若要呈現 INTERVENE 18 disease / CHD clinical risk，需找到並確效原研究 age/sex/cumulative incidence 轉換與適合台灣的 baseline incidence、疾病診斷、outcome association；現在 percentile **不是** absolute-risk。
- 建立 PRS/report versioning、研究倫理核准與資料存取控管，尤其是 50k WGS 的 PHI/基因資料傳輸、保留與稽核。

## 最新可用命令（只供新 run，不覆寫既有驗證）

```bash
cd /home/n102968/2026_INTERVENE_PRS/pipeline/PGS004879_integrated_v2
python3 pipeline_v21.py verify-baseline --config config.json
python3 -m unittest discover -s tests -q

# 由已有 manifest 複製並加入新 case；run-id 務必唯一
python3 pipeline_v21.py batch --config config.json \
  --manifest /absolute/path/next_batch.tsv \
  --run-id NEW_UNIQUE_RUN_ID --jobs 1

# 指紋完全相同的原 run 才能 resume
python3 pipeline_v21.py batch --config config.json \
  --manifest /absolute/path/next_batch.tsv \
  --run-id NEW_UNIQUE_RUN_ID --jobs 1 --resume
```
