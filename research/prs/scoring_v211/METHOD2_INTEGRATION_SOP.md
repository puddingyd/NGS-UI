# PGS004879 integrated v2：從 raw gVCF 到固定 Method 2 PRS 的單一入口

**科學方法基礎版本：v2.0（2026-10-31）；現行整合執行器：v2.1.1（文件更新 2026-11-02）。** Research use only；非 CHD absolute-risk／臨床驗證系統。

> **最新執行狀態：** 現已於伺服器通過三案 fresh gVCF（3/3）與固定 A/B（6/6）、兩案真正平行 fresh（2/2、4/4）、completed-run resume、post-roundtrip recovery。正式操作請用 `pipeline_v21.py`，見 [`V21_SOP.md`](V21_SOP.md)。以下保留 `pipeline_v2.py` 的 Model 2 科學邏輯與原始 v2 部署說明，供方法審查。

## 1. 要解決的問題與原則

原 `PGS004879_integrated_v1/pipeline.py`：

```text
raw gVCF → 已驗證 raw reconstruction / GQ, DP, FILTER, allele QC
         → strand-check / PGEN roundtrip / pgsc_calc raw score
         → target-dependent pgsc_calc ancestry/reference score → QC
```

Method 1 607,445-locus common3 只驗證固定交集能令 reference SUM 完全一致，**不是可向 50,000 人擴展的 canonical 模型**。Method 2 v0 則事先固定 PGS004879×HGDP+1kGP reference-compatible 的 **610,305 個**位點，避免因每個 target 的 QC/缺失不同而重算 reference distribution。

新版的預設流程：

```text
pipeline_v2.py run --mode fresh --gvcf SAMPLE.gvcf.gz
     │
     ├─ v1 vendor raw script（原算法未修改，保留 -resume）
     │   ├─ gVCF reference blocks + explicit ALT GT projection
     │   ├─ PASS, DP≥10, GQ≥20, dbSNP157 flip verification
     │   ├─ 22-chromosome PLINK2 round-trip / pgsc_calc raw SUM
     │   └─ raw/audit/all_score_locus_decisions.tsv
     │
     └─ method2_score.py（新）
         ├─ SHA256 model lock / 官方 scorefile identity
         ├─ 610,305 行 canonical source-row/REF/ALT/allele provenance crosswalk
         ├─ 合格 genotype：依 VCF GT/REF/ALT 重新核對 effect dosage
         ├─ 不合格 genotype：missing，代入 2×p_effect_3330
         ├─ Model A (610,305) 和 Model B (603,438)
         ├─ fixed reference EAS 735 ECDF percentile / sample SD Z
         └─ observed & imputed SUM / missing-locus audit / report
```

選用 `--with-legacy-ancestry` 時，會先執行**未修改的 v1 pipeline.py run**（包括 target-dependent ancestry + 原 QC），之後再執行 Method 2；預設則直接 raw → Method 2，**不需要再執行原先的 ancestry stage**。

## 2. 為何使用既有 canonical_v0 模型，而不是重新計算 reference？

現有 `/home/n102968/2026_INTERVENE_PRS/sensitivity/canonical_v0/` 已完成：

| 項目 | 結果 |
|---|---|
| 原始 PGS Catalog harmonized score | 610,677 rows；官方下載檔解壓 SHA256 與原始檔完全相同 |
| Reference-PVAR alleles compatible | Model A 610,305；其他 372 excluded |
| Fixed reference individuals | unrelated 3,330（原始 3,942、KING 排除 612）|
| EAS score distribution | unrelated EAS 735 |
| Model B | EAS 及 unrelated reference F_MISS ≤ 1%，603,438 loci |
| Scorefile SHA256 (A) | `358dc695afaee2387b14da0641f23a46b8d7cf33312f6ee14225eb9018d37f96` |
| Scorefile SHA256 (B) | `980541e5f1df2d5c5cd439f1f7020b68ce40a0af403bd6a6e44c02a62e675a42` |
| Allele provenance | 609,974 exact, 331 reverse-complement, 0 unexplained problems |
| Independent scoring | 三案 A/B 6/6 Python ↔ PLINK2 PASS；最大 SUM 誤差 4.58×10⁻⁷ |
| 官方來源解壓 SHA256 | `0ee9852635bc831bfbcbb144cb308157bed658e0e26fc4243737b3957aea4d3f` |
| 已有 `method2_v0_SHA256SUMS.txt` | 使用者驗證 17/17 檔案 OK |

因此新入口不重建 reference PGEN，不重新選位點，也不更動已確效的參數；它只**讀取固定的 master、frequency、scorefiles、reference `.sscore`**，計算新樣本。新 ZIP 內沒有這些大檔；需要原 server 的 `canonical_v0` 目錄。模型內容變更必須明確升版、重新建立新的模型鎖，不可直接覆蓋原 lock。

## 3. 部署與路徑

將 ZIP 解壓到：

```text
/home/n102968/2026_INTERVENE_PRS/pipeline/PGS004879_integrated_v2/
```

包內 `config.example.json` 已設定：

```json
"raw_script": "/home/n102968/2026_INTERVENE_PRS/pipeline/PGS004879_integrated_v2/vendor/PGS004879_FINAL_gVCF_resume.sh",
"method2_model_dir": "/home/n102968/2026_INTERVENE_PRS/sensitivity/canonical_v0",
"method2_lock": "/home/n102968/2026_INTERVENE_PRS/pipeline/PGS004879_integrated_v2/method2_model_lock.json"
```

```bash
PROJECT=/home/n102968/2026_INTERVENE_PRS
cd "$PROJECT/pipeline/PGS004879_integrated_v2"
cp -n config.example.json config.json
python3 -m json.tool config.json >/dev/null
python3 -m py_compile pipeline.py pipeline_v2.py method2_score.py qc_report.py
bash -n vendor/PGS004879_FINAL_gVCF_resume.sh
python3 -m unittest discover -s tests -v
```

其他 v1 dependencies：`nextflow`、`apptainer`、`bcftools`、`samtools`、indexed GRCh38 FASTA、dbSNP157、indexed gVCF、local PLINK2 image、Nextflow configs。若機器路徑與示例不同，編輯**本套件**的 `config.json`，不要改寫 vendor script。

## 4. 首次部署：固定 immutable model lock

```bash
python3 pipeline_v2.py freeze-model --config config.json
```

此指令先用 `sha256sum -c` 重查原 `canonical_v0/provenance/method2_v0_SHA256SUMS.txt` 17 個舊檔，再核對 source 原始內容、Model A/B scoring SHA、master/scorefile 一致性、reference 3330/EAS735 IID sets、無 provenance issues，產生新檔：

```text
pipeline/PGS004879_integrated_v2/method2_model_lock.json
```

`run` 每次會重查原始 source 的解壓 checksum 和**每個 lock 的檔案 SHA256**。若檔案已修改或被替換則 fail closed，不以新內容默默計算。首次產生 lock 時**不覆寫已存在 lock**。若原本固定 model 確有預期更動，應建立新目錄、新版本、新驗證而非直接更新 v0。

**注意**：freeze 依賴現有 17-file inventory 及模型資源。請不要清理 `canonical_v0`、`resources` 或正在使用的 Nextflow work files。

## 5. 先驗收三個已知樣本（無須重跑 raw gVCF）

```bash
python3 pipeline_v2.py regression --config config.json \
  --raw-run pilot_fresh_02 \
  --out-run method2_v2_regression_01
```

會讀取各樣本 `integrated_runs/pilot_fresh_02/<sample>/raw/audit/`，經目前新增的 canonical crosswalk / 2p imputation 重算 A/B；並比對**2026-10-31 已驗證的六組固定數值**：

| Sample | A SUM | A EAS percentile | B SUM | B EAS percentile |
|---|---:|---:|---:|---:|
| 26T00970 | 0.366424148908984 | 35.1020408163 | 0.363472617655232 | 35.3741496599 |
| 26T00969 | 0.474838562769280 | 60.5442176871 | 0.475804056110600 | 61.4965986395 |
| 26T00968 | 0.556804458044592 | 79.4557823129 | 0.552610418497840 | 79.4557823129 |

容許 `total_sum` 誤差 ≤ 2e-6（避免文字精度造成假錯誤），`percentile` 誤差 ≤ 1e-6 percentage point。最終應出現：

```text
METHOD2 THREE-SAMPLE REGRESSION: PASS (6/6)
```

結果獨立輸出至 `integrated_runs/method2_v2_regression_01/`，**不改動** `pilot_fresh_02`。**2026-11-02 狀態更新：** 此 v2 基準已完成實際 server regression 6/6；且 v2.1.1 三案 fresh、二案並行結果均重現相同 golden，詳 [`VALIDATION_REPORT_v211_2026-11-02.md`](VALIDATION_REPORT_v211_2026-11-02.md)。

## 6. 新樣本：一支程式從 raw gVCF 到 Method 2 PRS

```bash
python3 pipeline_v2.py run \
  --config config.json \
  --mode fresh \
  --sample SAMPLE_ID \
  --gvcf /absolute/path/SAMPLE_ID.hard-filtered.gvcf.gz \
  --run-id method2_fresh_01
```

重點：`SAMPLE_ID` 要與 gVCF 標頭內 sample ID 完全相同；gVCF 應有 `.tbi` 或 `.csi`。若使用已有 run ID，完整 raw output 必須已驗證而且 input SHA256 確認為同一個 gVCF；不會覆寫 incomplete raw directory。不同批次請使用**新的 run-id**。

每一個 sample 的輸出：

```text
$PROJECT/integrated_runs/method2_fresh_01/SAMPLE_ID/
├── raw/                         # v1 完整原始 genotype reconstruction / raw PGS
├── logs/raw_stage.log
└── method2/
    ├── inputs.lock.json         # raw audit checksum、model lock、source/config/code hashes
    ├── summary.json             # 完整結果及 genotype QC counters
    ├── scores.tsv               # A/B 兩筆 SUM、EAS Z、percentile、group
    ├── missing_loci.tsv.gz      # 每個 missing locus 的原因、beta、2p、imputed contribution
    ├── report.md                # 易讀摘要
    └── status.json              # PASS 或 FAILED
```

如也想輸出原 v1 的 ancestry percentiles/QC：

```bash
python3 pipeline_v2.py run \
  --config config.json --mode fresh \
  --sample SAMPLE_ID --gvcf /absolute/path/SAMPLE_ID.hard-filtered.gvcf.gz \
  --run-id method2_with_legacy_01 --with-legacy-ancestry
```

這會多執行 v1 `pgsc_calc --run_ancestry`，並把 `qc/qc.json` 的 ancestry/low-confidence 資訊寫入 Method 2 報告；**Method 2 的固定計分不依賴其 reference SUM**。

**重要**：預設 raw-only 模式不推論新 target ancestry；其 EAS percentile 是「假設以 EAS 735 人作比較」的條件性統計。若新樣本 ancestry 未知，建議使用 `--with-legacy-ancestry`，或採用經驗證的獨立 ancestry assignment，確認為 EAS 才作 ancestry-matched 解釋。

若已經有完整 raw stage，可以用另一個 run-id：

```bash
python3 pipeline_v2.py run --config config.json \
  --mode reuse --sample 26T00970 \
  --raw-source "$PROJECT/integrated_runs/pilot_fresh_02/26T00970/raw" \
  --run-id method2_reuse_01
```

## 7. Method 2 scoring 與嚴格檢查

逐行讀取 `all_score_locus_decisions.tsv`，以 chromosome+position 對應已經過 source-row / REF / ALT provenance 的 master（610,305 unique coordinates）。不只比對 raw contribution，**還重新核對 accepted `VCF_REF`、`VCF_ALT`、biallelic `VCF_GT`、ALT dosage 與 canonical effect allele dosage**。`CANDIDATE_DIRECT` 或 `FLIP_CONFIRMED` 才接受 observed。`LOW_GQ`、`LOW_DP`、`NON_PASS_FILTER`、`ALLELE_MISMATCH` 等皆當作 missing，而非 0/0。

計分：

```text
observed contribution = beta * genotype effect-allele dosage
missing contribution  = beta * 2*p_effect_unrelated3330
SUM = sum(observed) + sum(missing)
EAS percentile = 100 * count(EAS reference SUM <= target SUM) / 735
EAS Z = (SUM - EAS mean) / sample SD(EAS735)
```

Model A = 610,305 loci；Model B = reference missingness ≤1% 的 603,438 loci。全部 610,305 loci 必須在 raw audit 中找到；原始 scorefile 中**不在 Model A 的其他 372 loci**會單獨核對與 raw SUM 的貢獻關係，不當作 Method 2 observed/missing 的一部分。若 GVCF-derived Model A observed fraction 低於 99.5%，自動停止（目前 v0 的研究品質防線，未經臨床 cutoff validation）。

這份 package 的 research-output percentile 是之前經測試的 735-person EAS ECDF；不是 INTERVENE age/sex-based CHD absolute risk。

## 8. 常見問題與限制

- **新版跑不到 `method2_model_lock.json`**：先執行 `freeze-model`；它要有所有已驗證的 server 模型檔。
- **既有 model checksum mismatch**：停止；不要繞過 SHA 或直接更新 lock，先調查哪個模型檔被更動。
- **`raw audit ... missing`**：v1 raw 未完成或路徑錯誤；不能只用 variant-only VCF 計分。
- **原始 raw SUM 不等於 canonical observed SUM**：正常可能因固定 reference-compatible universe 排除了額外 loci；程式會獨立對帳兩邊 contribution。
- **VCF GT canonical allele mismatch**：不能自動 strand flip；保留錯誤並審查相關位點。對已驗證的 `rs496572`（chr11:54662726）00970 low quality，不可放寬 QC。
- **原 source metadata 的 C/T 與 canonical G/A**：331 個位點可合法 reverse complement；由 source-to-canonical provenance 和 v1 dbSNP157 共同支持；effect weight 不變。
- **原 v1 report 的 percentile 不等於 Method2 A/B percentile**：計分位點 universe 和 missing handling 不同，不應期待相等。
- **不提供新的 CHD absolute-risk**：尚未做台灣族群及 CHD phenotype outcome calibration。
- **Server status**：本次 ZIP 的 Python 與合成單元測試在作者 sandbox PASS，未替代使用者 server 的真實三案 regression。

## 9. 後續升級建議

目前已完成 6/6 server fresh regression 和 4/4 parallel regression，且 baseline/model lock 不變。接下來優先 **mid-workflow Nextflow partial-task interruption recovery**、5–10 樣本 throughput/resource test、shared image cache、封存 validation evidence；再進行獨立 ancestry QC、台灣 reference calibration、INTERVENE outcome-based absolute-risk。詳細驗收標準見 [`NEXT_STEPS_v211.md`](NEXT_STEPS_v211.md)。未經額外驗證，不應將 v0 視為可直接用於臨床 CHD 絕對風險報告的模型。

## 10. v2.1.1 科學核心與執行層分工（2026-11-02 補充）

- 本文 `pipeline_v2.py` 的固定 A/B 方法及 `canonical_v0` 模型均未修改。現行 production-candidate execution 由 `pipeline_v21.py batch/run` 負責樣本互斥、run fingerprint、安全 resume 與 report。
- `raw_work` 的實際伺服器設定已修正為 `/home/n102968/2026_INTERVENE_PRS`，不是任何一個特定樣本的資料夾；`raw/nextflow_work` 仍對每個 sample 隔離。
- V2.1.1 已有 27/27 unit tests；三案 fresh 6/6；兩案 parallel fresh 4/4。`verify-pilot` 是三案專用，兩案要使用四筆 `batch_summary_v21.tsv` 逐筆 golden compare。
- 現有唯一未實測的主要 recovery 情境為 **Nextflow task 完成一部分後故障，再真正使用 cache 恢復未完成部分**；不能將 first-start failure recovery 或 completed-run cache probe 當作此項驗證。
