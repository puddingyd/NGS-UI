# PGS004879（CHD）WGS → Raw PGS → HGDP + 1kGP Normalization → QC Report

**整合版操作手冊暨技術驗證、待解決事項清單**  
版本：`integrated-v1.0` ｜ 日期：2026-10-28 ｜ genome build：GRCh38 ｜研究用途（**尚非臨床驗證的 CHD 絕對風險工具**）

> **證據邊界**：這個封裝程式已在本機進行 Python／Shell 語法檢查、合成資料單元測試；**尚未在院內 server 實際執行新整合版的三個 gVCF**。本文列出的 26T00970／00969／00968 數值，是 2026-10 已完成的原始 pilot benchmark，**不是新整合器重新跑出的結果**。必須在 server 執行下列 regression 指令並通過比對，才算完成這次的新版本實測。不可將本機測試視為真實 WGS pipeline validation。

## 1. 目標與原則

本次不改寫已驗證的 gVCF 重建演算法，而是以 `vendor/PGS004879_FINAL_gVCF_resume.sh` 封裝：

1. **Raw stage**：gVCF → exact score loci／genotype QC／strand audit → corrected VCF／22 條染色體 PGEN round-trip → pgsc_calc v2.3.0 `--skip_ancestry` raw score；原腳本另外獨立對帳 `SUM`、`DENOM` 與 round-trip。
2. **Ancestry stage**：直接重用 Raw stage `roundtrip_22chr/pgen/chr1...chr22`，以相同 corrected scorefile 重新執行 `pgsc_calc v2.3.0 --run_ancestry`、HGDP + 1kGP v1、FRAPOSA OADP、`SuperPop`。
3. **QC/report stage**：自動擷取 pgsc_calc 的 `target SUM`、EAS percentile、Z-scores、reference counts、match summary、PCA LD-pruned markers、high-weight unmatched 變異、INTERVENE 七區間 boundary warning。
4. **Regression**：比較 26T00970、26T00969、26T00968 與 2026-10 已確認的 golden benchmarks；將錯誤、版本及參數保留在獨立結果目錄，不覆寫舊資料。

**不可混為一談**：原始 raw SUM 和 `--run_ancestry` 計算出的 SUM 各自使用不同的 scoring variant set，**不要求相等**，必須各自對照其 baseline。原先 reference score 00970 vs 00969 的 differential scoring 已以 1,861 SNP 解釋，不能把 reference SUM 不同誤認為 pgsc_calc 計算失敗。

### 檔案內容

| 檔案 | 功能 |
|---|---|
| `pipeline.py` | 主執行器：fresh / reuse、鎖定、輸入指紋、工作目錄、raw→ancestry→QC |
| `vendor/PGS004879_FINAL_gVCF_ORIGINAL.sh` | 過去已驗證的 raw script 原始副本；保留供審核／SHA256 comparison |
| `vendor/PGS004879_FINAL_gVCF_resume.sh` | 上述原始副本**僅新增** Nextflow `-resume`；作為本次實際 raw stage |
| `qc_report.py` | 建立每樣本 `qc.json`、`qc.tsv`、`report.md`、`missing_high_weight.tsv` |
| `compare_benchmarks.py`／`benchmarks.json` | 驗證 3 案例精確基準與容許誤差 |
| `discover_gvcfs.py` | 從原始 raw audit `source_sha256.txt` 找回三案 gVCF 真正檔案路徑，不猜測 Novaseq batch |
| `config.example.json` | 伺服器資源路徑與 QC、Nextflow 參數 |
| `manifest.example.tsv` | 三案例樣本表，包含既有 raw-dir 與（待填的）gVCF |
| `tests/test_logic.py` | 小型離線測試，不使用病人基因型 |

---

## 2. 伺服器安裝與檢查

建議存放在現有專案（下列程式假設整個資料夾複製到指定位置）：

```bash
PROJECT=/home/n102968/2026_INTERVENE_PRS
mkdir -p "$PROJECT/pipeline"
# 把下載的 zip 上傳到 server 後，在適當目錄 unzip，讓實際路徑為：
# $PROJECT/pipeline/PGS004879_integrated_v1/pipeline.py
cd "$PROJECT/pipeline/PGS004879_integrated_v1"

cp -n config.example.json config.json
python3 -m json.tool config.json >/dev/null
bash -n vendor/PGS004879_FINAL_gVCF_resume.sh
python3 -m py_compile pipeline.py qc_report.py compare_benchmarks.py discover_gvcfs.py
python3 -m unittest discover -s tests -v
```

**先確認 `config.json` 所有路徑是真實存在的檔案。** 尤其：

- 原始 **unmodified** scoring file：`resources/scorefiles/PGS004879_hmPOS_GRCh38.txt.gz`，610,677 loci；不得以 00970「已針對其 flips 修正過」的 scorefile 當成新樣本的 original scorefile。
- FASTA：`/home/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta` 和 `.fai`。
- dbSNP157：`resources/dbsnp157/GCF_000001405.40.gz` 及 `.tbi`/`.csi`。
- Raw Nextflow config：`resources/config/pgsc_calc_local.config`；Ancestry config：`resources/config/pgsc_calc_ancestry.config`。
- PGS reference：`resources/reference_panels/HGDP_1kGP/v1/pgsc_HGDP+1kGP_v1.tar.zst`。
- PLINK2 image：`resources/containers/ghcr.io-pgscatalog-plink2-2.00a5.10.img`。
- `nextflow`、`apptainer`、`bcftools`、`samtools` 均需在 PATH。

如原始 `resources/...` 與本機設定不同，**只改自己的 `config.json`**；不應直接更動 vendor 的 QC／allele reconstruction 程式。`raw_work` 指向舊腳本的歷史工具 root（`.../26T00970`），但 input FASTA、dbSNP、scorefile、image、config 都由本次 `config.json` 額外明確指定，不讀取 26T00970 的歷史 genotype 當其他樣本資料。

執行均在**前景**，不使用 tmux、nohup，也不需要保留瀏覽器／ChatGPT 連線。

---

## 3. 建議先做兩級實測

### 3A. Reuse raw：先驗證新封裝的 ancestry／QC（較快）

三個原始 raw directories 已經完成並通過 `audit/final_reconciliation.json`，先在**新** run-id 執行，避免影響已完成的 ancestry pilot：

```bash
cd /home/n102968/2026_INTERVENE_PRS/pipeline/PGS004879_integrated_v1
python3 pipeline.py regression \
  --config config.json \
  --manifest manifest.example.tsv \
  --expected benchmarks.json \
  --run-id pilot_reuse_01 \
  --mode reuse \
  --keep-going
```

這一步**沒有重新讀取 gVCF 或重跑 Raw PGS**；它使用已完成、已對帳的 `raw_source`，但重新執行 `pgsc_calc --run_ancestry`、QC 與數值 regression。`--keep-going` 讓其中一案失敗後仍測試另外兩案，最後照樣會標記整體失敗。

### 3B. Fresh full：三個病例從 gVCF 全新重跑（本次最終驗收）

先從三個既有 raw log 的 `audit/source_sha256.txt` 擷取**真實來源路徑**，不用猜測 00969、00968 所屬的 NovaSeq batch：

```bash
cd /home/n102968/2026_INTERVENE_PRS/pipeline/PGS004879_integrated_v1
python3 discover_gvcfs.py \
  --project /home/n102968/2026_INTERVENE_PRS \
  --out pilot_3samples_with_gvcf.tsv
cat pilot_3samples_with_gvcf.tsv
```

請確認各樣本 `gvcf` 欄為**非空且存在**的 indexed `.gvcf.gz`，並對照原始 audit 的 sample ID。若 `source_sha256.txt` 不存在，手動補上當初使用的**真實**來源；不要自行推測檔名。

```bash
python3 pipeline.py regression \
  --config config.json \
  --manifest pilot_3samples_with_gvcf.tsv \
  --expected benchmarks.json \
  --run-id pilot_fresh_01 \
  --mode fresh \
  --keep-going
```

`fresh` 會真正執行：gVCF → variant QC → strand handling → corrected scorefile → 22-chr PGEN round-trip → raw pgsc_calc → ancestry pgsc_calc → report → 三案例 benchmark。原本 3 案 Raw score / normalized score 皆會比對。**只有看到 `REGRESSION: PASS` 才能稱本次整合重跑成功。**

如果要對一個新 sample 獨立跑：

```bash
python3 pipeline.py run \
  --config config.json \
  --sample NEW_SAMPLE \
  --gvcf /absolute/path/NEW_SAMPLE.hard-filtered.gvcf.gz \
  --run-id batch_20261028_01 \
  --mode fresh
```

**若只想重檢原本 ancestry 的 Nextflow cache：**在同一 run-id 重跑 `run` 或加 `--recheck-ancestry`，Nextflow 一律帶 `-resume`。不會刪除或覆寫既有 raw results。

### 3C. Regression pass criteria

`benchmarks.json` 將以下數字當成 2026-10 pilot golden data（不是本次新跑結果）：

| Sample | Raw SUM | Raw matched | Raw DENOM | Ancestry SUM | EAS percentile | pgsc matched / excluded / unmatched | PCA markers |
|---|---:|---:|---:|---:|---:|---|---:|
| 26T00970 | 0.351907 | 609346 | 1218692 | 0.354661 | 35.6463 | 609032 / 314 / 1331 | 27772 |
| 26T00969 | 0.471545 | 609488 | 1218976 | 0.475192 | 60.6803 | 609162 / 326 / 1189 | 27803 |
| 26T00968 | 0.559173 | 609597 | 1219194 | 0.558069 | 79.5918 | 609272 / 325 / 1080 | 27764 |

此外三位均：`MostSimilarPop=EAS`、`RF_P_EAS=1.0`、`LowConfidence=False`、reference individuals 3330、EAS 735。相應 `Z_MostSimilarPop` 約 -0.34448、+0.3394、+0.8658。

比較規則：Raw / normalized `SUM` 絕對誤差不大於 `2×10^-6`、percentile 不大於 0.15 百分點、Z 不大於 0.0002；matched／excluded／unmatched、原始 matched、DENOM、reference counts、PCA markers 要**完全一致**。這是為**同版本、同設定、同輸入**設計的嚴格 regression。容忍度不應因測試失敗而任意放寬；應先排查 reference、軟體／設定或資料版本是否不同。

若 raw 源／新樣本 genotype 已改變，嚴格基準可能合理地失敗；必須評估影響後建立**新版本基準**，不能靜默更新 `benchmarks.json`。

---

## 4. 輸出目錄與重跑／失敗安全性

每個 run 獨立建立：

```text
$PROJECT/integrated_runs/<RUN_ID>/
├── regression_execution.tsv         # regression 逐案實際執行狀態
├── regression_comparison.tsv        # golden benchmark 詳細逐欄比對
├── regression_comparison.log
└── <SAMPLE>/
    ├── run_inputs.json              # 參數 + script/config SHA256 + gVCF/reference metadata
    ├── raw_inputs.json              # completed raw reconciliation/score/sheet SHA256 + PGEN stats
    ├── status.json                  # STARTED / RAW_VALID / ANCESTRY_FINISHED / COMPLETE / FAILED
    ├── raw/                         # fresh 模式：獨立新產生；reuse 模式：讀取舊目錄而不複製
    │   ├── audit/final_reconciliation.json
    │   ├── audit/all_score_locus_decisions.tsv
    │   ├── scorefiles/PGS004879_hmPOS_GRCh38.txt.gz
    │   ├── roundtrip_22chr/pgen/chr*.{pgen,pvar,psam}
    │   └── results/case<SAMPLE>/score/aggregated_scores.txt.gz
    ├── ancestry/
    │   ├── inputs/samplesheet.csv    # raw PGEN 的絕對路徑
    │   ├── work/                     # Nextflow -resume workspace
    │   └── results/case<SAMPLE>/
    │       ├── score/case<SAMPLE>_pgs.txt.gz
    │       ├── score/case<SAMPLE>_popsimilarity.txt.gz
    │       └── match/case<SAMPLE>_summary.csv
    ├── qc/
    │   ├── qc.json
    │   ├── qc.tsv
    │   ├── report.md
    │   └── missing_high_weight.tsv
    └── logs/{raw_stage,ancestry_stage,qc_stage}.log
```

**重跑原則**：

- 同一 `run-id + sample` 的參數／原始檔案 metadata，及已完成 raw 的 reconciliation／corrected scorefile／samplesheet SHA256、PGEN metadata 均會核對；變動就拒絕，要求另取新 run-id。
- 每樣本獨立 `flock`，避免同時兩個程序改動同一樣本工作區。
- Ancestry 永遠 `nextflow -resume`；成功輸出若已存在，可直接重新產出 QC，也可以 `--recheck-ancestry` 強制讓 Nextflow 檢查 cache。
- Raw vendor script 原始設計是 **new output only**。若 raw 工作目錄建立後失敗，因其內部非 checkpoint/resume 設計，封裝器**不會**強行續跑部分 raw 目錄（即使裡面某個 Nextflow command 有 `-resume`）；請換一個 run-id 重新跑 raw，保留失敗的 audit。
- `flock`、`run_inputs.json` 可避免一般誤覆寫；**不是**跨檔案系統備份或資料不可竄改保證。raw script `source_sha256.txt` 才有當時輸入 gVCF 的 full hash，之後要保存以建立更強 provenance。
- `reuse` 模式不能作為「從 gVCF 全新重跑」的驗證證據。
- 本套腳本**沒有**使用 tmux/nohup/背景工作。

---

## 5. QC 報告內容與判讀

`qc/report.md`／`qc.json` 記載：

- Raw `observed_sum`、`DENOM`、原始 matched；raw stage 已完成的 22-chr VCF→PGEN→VCF round-trip／獨立 SUM reconciliation。
- pgsc_calc `SUM`、`Z_MostSimilarPop`、`Z_norm1`、`Z_norm2`、`percentile_MostSimilarPop`。
- PCA `MostSimilarPop`、`RF_P_EAS`、`LowConfidence`、LD-pruned marker count。
- Match `matched`、`excluded`、`unmatched`；reference 3330 人與 EAS 735 人。
- `abs(effect_weight)>=0.001` 且沒有 matched 的 loci：列入 `missing_high_weight.tsv`；這是**人工檢視優先順序**而非臨床驗證門檻。00970 的 `rs10116277` 權重 0.01395，理應被標出。
- INTERVENE 七組 `<20, 20–40, 40–60, 60–80, 80–90, 90–95, >=95`；距離邊界 <2 百分點就 `boundary_review=True`。此 **2pp 是 pilot 的人工複核觸發，不是統計信賴區間**。
- `qc_status`: `FAIL`（硬性檢查失敗）／`PASS_WITH_REVIEW_FLAGS`（成功但有 QC／interpretation warning）／`PASS`。因為 high-weight missing 和 boundary 只是警示，**不能將 PASS_WITH_REVIEW_FLAGS 等同 FAILURE，也不可當成臨床報告核准**。

### 重要：為何三次 Reference SUM 不完全相同？

已完成的 dosage-level audit：

| 對照 | 原始 mean SUM 差異 | 說明 |
|---|---:|---|
| 00970 vs 00969 | +0.01292255 | 00970 的 `chr9:22081398 rs10116277` 為 G/T、GQ15、`LOW_GQ`，故沒有進入 PVAR；00969/00968 均 T/T、GQ48/50，正常 matched |
| 00970 vs 00968 | +0.01297388 | 相似情形，少數其他 matched variants 差異 |
| 00969 vs 00968 | +0.00005134 | Reference score Pearson r=0.999963，median percentile move 0.1361 pp |

734/735 位 EAS reference 有 rs10116277 的實際 G dosage；平均 G dosage = **0.63079019**，weight +0.01395，預測貢獻 +0.00879952，可解釋 00970 vs 00969 平均 reference shift 的約 **68.11%**；減去單位點後 difference SD 由 0.00945235 降為 0.00180309。

最後直接用兩份 corrected reference scorefiles 的 **1861 個 differential scoring variants** 重算，結果：Observed mean shift = `0.01292255`、Predicted = `0.01292637`、Residual mean = `−0.00000383`、Residual SD = `0.00000616`、Max abs residual = `0.0000904`。這提供了極強的**數值歸因證據**：Reference difference 主要是不同 target-specific matched scoring set，而不是 scoring weight、sample set 改變。小 residual 尚未做精確 missing-imputation/frequency 逐位點驗證。

**不要為了讓 reference SUM 接近就把 00970 的 GQ15 改成 PASS**。原本 GQ>=20 的 QC 有其意義；若需對低 GQ 作敏感度分析，須在獨立版本呈現並做 genotype-level validation。

---

## 6. 已知侷限、尚待解決的問題（建議保留於每一版 SOP）

### A. 近期正式上線前要處理

1. **固定分數可比性與 reference 定義**：目前 `pgsc_calc v2.3.0` 的 matched variants 依每位 target 改變，reference 被重算；相同 735 位 EAS 個體的 PGS 可能不同。先保留官方 target-specific reference normalization 作為 pilot；若要跨樣本、跨批次、跨疾病穩定比較，設計 versioned scoring variant universe、相同 genotype dosage/missingness 定義、可重複的 reference normalization。**不能單純固定 reference SUM，卻允許 target 以不同 SNP subset 計分。**
2. **High-impact missing variant handling**：`rs10116277` 等權重較大之 SNP 缺失會影響 reference distribution。需針對 high-weight missing count、絕對加權潛在影響、GQ/DP/AD/PL、IGV、strand、imputation sensitivity 設定一致 SOP；不建議無差別放寬 GQ。
3. **Boundary uncertainty**：00969 = 60.6803、00968 = 79.5918，皆距 INTERVENE 60/80th boundary <1pp。即使整體 reference rank 很穩定，也不保證個別人的 HR group 穩定；日後應評估 reference/scoring perturbation、bootstrap、imputation sensitivity。
4. **Cohort batch vs per-case scoring**：正式大量 WGS 需比較 cohort-input `pgsc_calc` 和現有 per-sample baseline，固定 missing-genotype policy、variable sites representation 與 output normalization。單純把多樣本放到同一次 run，不等於確保 score 可比性。
5. **Input／reference provenance**：固定 PGS Catalog source/版本/SHA256、GRCh38 FASTA、dbSNP157、PLINK2 image、pgsc_calc 2.3.0、reference archive、Nextflow config、GQ/DP/strand parameter。`run_inputs.json` 記錄小檔 SHA／大檔 size+mtime；若需完整審核級不可否認性，應為所有大 reference 資料建立 checksum manifest。後續請將本專案納入 Git version control；目前 zip 是不可變版本封存但非 Git commit。
6. **本整合器 fresh-mode roundtrip**：保留原 vendor full raw pipeline，raw 失敗不允許在半成品原目錄續跑；改 run-id 可完整重跑。未來可將 raw 改成有 stage-level atomic checkpoint 的 Nextflow workflow，但在重構之前不應宣稱 raw 具有完整 -resume。
7. **Reference normalization formula**：確認並版本化 `Z_MostSimilarPop`、`Z_norm1/2` 具體建模方式；不應將 empirical percentile 與 PCA-adjusted Z 混稱；執行多個 ancestry groups 前擴充。
8. **Reference allele frequency imputation**：PLINK2 `--score` 的 mean imputation 應核對實際 `--read-freq`，尤其 735 EAS 中 HGDP00782 在 rs10116277 missing；先前獨立差分計分有微量 residual，不能武斷視為 rounding only。

### B. 方法學與族群外推問題

9. **PCA markers 不是獨立 genome-wide panel**：三個案例約 27,764–27,803 LD-pruned SNP，來自 PGS-specific intersection；高信心 `MostSimilarPop=EAS` 支持 continental ancestry pilot，但未來需獨立 genome-wide PCA、relatedness/contamination check、reference projection stability。`RF_P_EAS=1.0` 是分類信心，不等於 100% EAS 血統比例。
10. **EAS vs Han-related vs 台灣 reference**：00970 既有 EAS percentile=35.65%，Han-related 252 人池 percentile=40.48%，可能跨 40th INTERVENE risk stratum。群體選擇對輸出有實質影響，HGDP+1kGP 不代表台灣代表性族群。
11. **Taiwan-specific population WGS**：院內約 1,500 WGS 可先探索性分析，但來源是否為罕病、兒科、家系或選擇性抽樣、relatedness 及 age/sex/ancestry 結構，均影響 reference 代表性；不能直接當作一般台灣人口分布。需要去親緣、ancestry、批次／測序平台 QC，並在獨立 cohorts 驗證。PRS percentile tail 需要更大有效樣本量。
12. **INTERVENE external HR transferability**：目前 20/40/60/80/90/95th 七組 HR 來源為原研究 cohort／model；不可假設直接搬用到台灣／East Asian族群仍校準。即使 PRS percentile 很準，也不等於能預測個人 CHD 的 absolute risk。
13. **Taiwan absolute CHD risk 待完成**：後續若用 GBD 台灣 `location_id=8`，須核對 IHD vs CHD 定義、性別年齡別 incidence/prevalence、competing other-cause mortality、5 年 interval hazard、baseline year、外推假設；需 local outcome cohort validation/recalibration，且不同資料來源的死亡／盛行／發生定義必須一致。**現階段不可把 INTERVENE 的 absolute risk 計算當成臨床確定風險。**
14. **其他疾病 PGS 推廣**：不能把 PGS004879 專用 scoring file、QC 基準、score range、reference interpretation 套到其他 PGS ID；將來需每疾病版本化 score manifest 和 phenotype-specific validation。

### C. Variant representation／統計細節

15. **Palindromic A/T、C/G**：原 raw 腳本為 historic compatibility 使用 `--palindromic-mode assume-direct`；此處仍存在 strand ambiguity 對照頻率與 dbSNP 資訊不足的可能。要臨床應用之前需嚴格審核或設計 `exclude` sensitivity benchmark；**不得自行改成 exclude 卻繼續用舊 golden SUM 宣稱 PASS**。
16. **Reference block 與低頻 genotype**：原腳本會從 gVCF reference block 推定某些 score loci 的 `0/0`，但 DIRECT block 的 identity audit 不逐位另行用 rsID/dbSNP 全確認（已驗證的 FLIP 則經 dbSNP）。這部分仍待進一步 evidence 強化；gVCF 的 `MIN_DP`、`GQ`、block completeness 與多等位 allele projection 均需持續抽查。
17. **Match count 的分母**：整套 scoring file 共 610,677 loci；pgsc_calc matched / excluded / unmatched 應按完整 scoring file 算，與 Raw stage 自己的 matched 609,346／609,488／609,597 **不是同一指標**。
18. **High-weight missing cutoff**：abs weight >=0.001 目前只是報表的 review threshold，不能代替正式 PGS uncertainty computation。應評估 single-locus potential contribution、missing allele frequency 和 group crossing risk。
19. **Reference rank stability**：00970 vs 00969 Pearson r=0.99820067、reference individuals changing group 50/735；00970 vs 00968 r=0.99818834、48/735；00969 vs 00968 r=0.99996319、10/735。這些是**同一批 HGDP reference 在不同 target-specific scoring sets 下的敏感度分析**，不是台灣人口風險分組的真正錯誤率。

---

## 7. 常見錯誤排除

| 症狀 | 判斷／處理 |
|---|---|
| `Missing configured file` | 確認 `config.json` 對照 server 位置，尤其 raw config 與 ancestry config 不同 |
| `Missing gVCF index` | 不要猜 index；確定 `.tbi` 或 `.csi`，必要時使用正確的 `bcftools index` |
| `Raw directory exists but is incomplete` | 原 raw script 拒絕覆寫；用新的 `--run-id`，保留舊 audit，勿自行刪原始工作目錄 |
| `Run inputs changed` | 同一 run-id 的設定／來源更改，封裝器拒絕混入；新版本用新 run-id |
| `Sample already has ... in progress` | 同 sample per-sample flock，檢查是否已有前景作業；不要用 `kill -9` 繞過正在跑的資料寫入 |
| `REFERENCE PGEN NOT FOUND` | 只有做 dosage-level differential audit 時才需要 reference PGEN；主整合流程只依官方 `--run_ancestry` 和 reference archive 自動處理 |
| `REF / dbSNP / gVCF index absent` | 需索引與 GRCh38 build 完整相容；不可只依 rsID 猜 build |
| normalized SUM 與 raw SUM 不相等 | **預期可能不相等**：兩者 matched sets 不同；應各自和 golden baseline 比對 |
| `LowConfidence=True` 或 non-EAS ancestry | 先人工 review；pilot golden EAS cutoffs 不適合直接套用不同 ancestry |
| 00970 rs10116277 high-weight unmatched | 已確認 `GQ=15`、GT=0/1、DP=74；維持原 QC，不要為通過 regression 直接強制計分 |
| `REGRESSION: FAIL` | 查看 `$PROJECT/integrated_runs/<run-id>/regression_comparison.tsv`，先檢查真正輸入版本／配置／流程是否一致，再決定是否有新 protocol version |

## 8. 驗收清單

- [ ] 在 server 執行 `bash -n`、`py_compile`、`unittest` 均 PASS。
- [ ] `config.json` 真實檔案、container、reference archive、input gVCF／indexes 都存在。
- [ ] `pilot_reuse_01` 的 `REGRESSION: PASS`，且三份 `qc/report.md` 可開啟。
- [ ] `pilot_fresh_01` 以**三個確實來源 gVCF**跑完，`REGRESSION: PASS`。
- [ ] 特別核對 00970 的 `rs10116277` 仍 `LOW_GQ`、不進入 raw/normalized matched genotype。
- [ ] 00969／00968 的 boundary review 為 True；00970 False。
- [ ] 保留 SHA256 provenance、`run_inputs.json`、`raw_inputs.json`、`status.json`、Nextflow work、QC JSON 和 regression TSV。
- [ ] 將 `integrated-v1.0` 與原始 source file hash 記錄到 Git 或院內版本控管。
- [ ] **完成上述各項之前，報告只稱為 technical pilot，不稱 clinically validated CHD risk pipeline。**

## 9. 相關方法與程式來源（供版本化核對）

- PGS Catalog / `pgsc_calc` v2.3.0：<https://github.com/PGScatalog/pgsc_calc/tree/v2.3.0>
- PLINK2 scoring：<https://www.cog-genomics.org/plink/2.0/score>
- PGS Catalog calculator ancestry documentation：<https://pgsc-calc.readthedocs.io/en/latest/how-to/ancestry.html>
- INTERVENE flagship public code：<https://github.com/intervene-EU-H2020/flagship>
- 本專案已驗證的 `PGS004879_FINAL_gVCF.sh`、原版 SOP、3 案例 ancestry pilot、rs10116277／1861-loci reference differential scoring audit（2026-10）；本手冊包含相關數值，但新整合腳本仍須執行上述 3-case regression。
