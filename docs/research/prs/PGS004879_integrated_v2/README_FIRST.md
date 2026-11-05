# PGS004879 integrated v2 — raw gVCF → fixed PRS Model A/B

這是由使用者提供的 `PGS004879_integrated_v1(2).zip` 衍生的**新增整合版**。原版 `pipeline.py`、`vendor/PGS004879_FINAL_gVCF_resume.sh`、`vendor/PGS004879_FINAL_gVCF_ORIGINAL.sh`、既有 benchmark **保持不變**；V2 新增 `pipeline_v2.py` 與 `method2_score.py`。

## 核心用途

使用同一支入口程式，完成 **fresh raw gVCF → v1 已驗證的 genotype reconstruction／QC／raw PGS → 固定 Method 2 Model A/B 的 target SUM、EAS Z-score、percentile、INTERVENE strata**。預設**不執行**舊的 target-dependent ancestry stage；若需要原本 v1 的報告，可加上 `--with-legacy-ancestry`。

**本 ZIP 只包含程式和操作文件，沒有附上 610,305 個 loci 的模型大檔。** 預期會在原伺服器讀取已驗證的 `sensitivity/canonical_v0/`，包括 `reference_qc_master_v0.tsv.gz`、兩份 scoring files、3,330 人的 `.afreq`、735 人 EAS keep list 與 Model A/B reference `.sscore`。第一次執行須明確 `freeze-model`，產生完整 SHA256 model lock。未通過 checksum 的模型會停止計分。

## 三個步驟（於伺服器執行）

```bash
PROJECT=/home/n102968/2026_INTERVENE_PRS
cd "$PROJECT/pipeline/PGS004879_integrated_v2"
cp -n config.example.json config.json
python3 -m unittest discover -s tests -v

# 只在初次部署執行一次：既有模型資源完整檢查與鎖定
python3 pipeline_v2.py freeze-model --config config.json

# 使用先前三案已完成的 raw audits 回歸驗證：應有 6/6 PASS
python3 pipeline_v2.py regression --config config.json \
  --raw-run pilot_fresh_02 --out-run method2_v2_regression_01

# 正式處理另一個單樣本（請換成真正的 sample ID 和 indexed gVCF）
python3 pipeline_v2.py run --config config.json \
  --mode fresh --sample SAMPLE_ID --gvcf /absolute/path/SAMPLE_ID.hard-filtered.gvcf.gz \
  --run-id method2_fresh_01
```

模型計分結果在 `$PROJECT/integrated_runs/<run-id>/<sample>/method2/`：`summary.json`、`scores.tsv`、`missing_loci.tsv.gz`、`report.md`、`status.json`、`inputs.lock.json`。

若先前 raw reconstruction 已完成，可使用 `--mode reuse --raw-source /absolute/path/to/completed/raw`，不用重新解析 gVCF。**指令和注意事項詳見 [`METHOD2_INTEGRATION_SOP.md`](METHOD2_INTEGRATION_SOP.md)。**

## 狀態與證據邊界

- 舊版三案 gVCF pipeline 的 server regression 及 Method 2 Python↔PLINK2 六組 scoring 已由使用者先前執行並驗證。
- 本次新整合入口已在本機執行 Python compile + synthetic tests；**尚未在使用者 server 實際跑完新入口的三案 regression**。因此部署後務必先跑 `regression` 再用於新的 WGS。
- 完成技術驗證不代表已具台灣人群 CHD outcome validation，EAS percentile 不是 CHD absolute risk。

**重要 ancestry 註記**：預設 raw-only 的 EAS percentile 僅表示相對於固定 EAS reference 的位置，並未替新樣本做 ancestry assignment；未知 ancestry 請加上 `--with-legacy-ancestry` 或另行完成 ancestry QC。
