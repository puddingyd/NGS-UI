# GPN-Star V100 / M447：ClinVar missense variants 之 PP3/BP4 證據校正研究紀錄

**研究進度：Step 1–26；截至 2026-10-30**  
**文件版本：v0.3（Step 1–21 維持 v0.2 Archive audit；Step 22–26 依研究者提供的 server 終端輸出補充；非最終論文 Methods）**  
**下一階段：Step 27 — calibration-only 的連續 score-to-LR 方法、外部驗證與正式證據門檻規劃**  
**專案根目錄：** `/home/n102968/GPN-star_calibration_clinvar20261025`

> **2026-10-30 更新摘要（v0.3）**：主要分析維持 **expanded autosomal cohort，48,271 variants（P/LP 20,089；B/LB 28,182；5,992 genes）**、V100 signed 分數 `S=-llr_calibrated`、原 Step 19 gene-disjoint split。Step 21 探索性 validation 支持 **`2<S≤5` BP4 Moderate、`8<S≤9` PP3 Supporting、`9<S≤11` PP3 Moderate** 候選；兩端 Strong 因少數類別太少未確認。Step 22–26 的 AF、gene-universe 和 pathogenic-gene eligibility analyses 已完成，顯示 **any-consequence P/LP/P-only gene restriction 下上述主結論大致穩健；expanded ultra-rare（AF<10⁻⁴ + 高 AN inferred-unobserved）則降低 PP3 LR**。所有候選證據仍限研究用途，詳見第 17–23 節。

> **文件稽核範圍**：本版已直接解壓並檢查使用者提供的 `Archive.zip`，讀取其中 `scripts/`、`logs/`、`results/`，另以 Python 從封存輸出重新計算主要 cohort／gene／AUC／LR；47 項獨立一致性檢查全部通過，且 Step 20（36 列）與 Step 21（14 列）以同一程式、相同 seed 在隔離目錄重跑後逐欄一致。**未含於 ZIP 的原始 ClinVar VCF、gnomAD source VCF、all-sites AN BED、GPN-Star 原始 Parquet，仍無法做上游原始檔層級驗證**；Step 4 的產生指令、Step 6 的獨立執行步驟均未提供。所有「已確認」僅限於實際可追溯證據，不代表已完成外部臨床驗證。

---

## 1. 研究目的及基本概念

### 1.1 目的

使用同一批經品質控制的 ClinVar germline missense SNVs，比較 **GPN-Star V100 與 M447** 的 prediction performance，並嘗試在符合 ClinGen SVI／ACMG/AMP 解讀框架的條件下，估計 **PP3（偏向 pathogenic）及 BP4（偏向 benign）** 的候選證據強度與分數區間。

研究分兩層：

1. **Discrimination**：在 ClinVar P/LP 與 B/LB 之間比較 ROC AUC、average precision（AP）及 score distributions。
2. **Evidence calibration**：在明確定義的 reference cohorts 中估計 class-conditional score likelihood ratio（LR），利用 gene-disjoint calibration / validation 和 gene-cluster bootstrap 檢查可重現性。

本研究的 LR 定義為

\[
LR(I)=\frac{P(S\in I\mid\mathrm{P/LP})}{P(S\in I\mid\mathrm{B/LB})}
=\frac{n_{P,I}/N_P}{n_{B,I}/N_B}
\]

其中 \(I\) 是 **互不重疊的 score interval**，不是「大於某 cutoff 的所有變異」；\(N_P,N_B\) 是**同一分析 cohort、同一 split** 內 P/LP 與 B/LB 總數。此 LR 是在選定 ClinVar reference set 下估計的**證據比**，並非病患的 posterior pathogenic probability。

### 1.2 GPN-Star scores 與方向

本研究由本機 Parquet 讀取兩種欄位：

- `llr_calibrated`：signed model score；為方便臨床分析，定義 **\(S=-\mathrm{llr\_calibrated}\)**，使 \(S\) 越高代表越強的 predicted deleteriousness／constraint。
- `abs_llr_calibrated`：作者另行提供的 calibrated absolute/effect magnitude score；**不是** `abs(llr_calibrated)`。它是 secondary predictor，不可從 signed score 自行取絕對值代替。

V100 signed（\(S\)）經 Step 17 選為目前的 **primary predictor**；M447 signed、兩套 absolute scores 為 secondary/comparator analyses。這項模型選擇是在初步整體 AUC 結果已被查看後做出，屬研究過程中的探索性選擇；不可聲稱完全事前指定。

### 1.3 臨床證據門檻（研究框架）

| Evidence strength | Pathogenic-direction LR 下限 | Benign-direction LR 上限 |
| ----------------- | --------------------------:| ----------------------:|
| Supporting        | ≥ 2.406                    | ≤ 0.4156               |
| Moderate          | ≥ 5.790                    | ≤ 0.1727               |
| Strong            | ≥ 33.53                    | ≤ 0.0298               |
| Very Strong       | ≥ 1,124                    | ≤ 0.0009               |

上述為採用的 Bayesian/ClinGen SVI 參考 LR 門檻，**不是**已證實的 GPN-Star 臨床 score cutoffs。須以校正資料、保守不確定性、獨立驗證和適當 disease mechanism 來決定是否可賦予實際 PP3/BP4。後述 Step 21 僅標記 **candidate**。

---

## 2. 資料來源、版本及環境

| 資料／工具       | 版本／內容                                                 | 目前使用方式                                                                                                 |
| ----------- | ----------------------------------------------------- | ------------------------------------------------------------------------------------------------------ |
| ClinVar     | 2026-10-25 release、GRCh38                             | P/LP vs B/LB，review status ≥2 stars，missense SNV reference selection                                   |
| Ensembl VEP | v115                                                  | All-transcript annotation；ClinVar gene matching；MANE Select／MANE Plus Clinical／canonical fallback 相關旗標 |
| gnomAD      | v4.1 exomes + genomes                                 | 使用本機第三方 joint MANE-filtered VCF；`AC_e/AN_e`、`AC_g/AN_g` 計算 joint AF；另用 all-sites AN BED 檢視未觀察位點        |
| GPN-Star    | V100、M447，hg38，200M                                   | 預先計算的本機 per-chromosome Parquet，含 signed 和 absolute calibrated LLR                                      |
| 分析環境        | Linux server；Python virtualenv `.venv_gpnstar_scores` | Polars、NumPy、SciPy、scikit-learn；Python 3.12／Polars 2.1.0 曾確認可用                                         |

**GPN-Star Parquet 路徑：**

```text
/home/n102968/NGS_UI/biotools/gpn_star/
  gpn-star-hg38-v100-200m/llr/llr_chr{1..22,X,Y}.parquet
  gpn-star-hg38-m447-200m/llr/llr_chr{1..22,X,Y}.parquet
```

已確認 Parquet schema：`chrom: String`、`pos: Int64`、`ref: String`、`alt: String`、`llr_calibrated: Float32`、`abs_llr_calibrated: Float32`。研究使用 `CHROM:POS:REF:ALT` 組成的 `variant_key`，在 **GRCh38、1-based SNV 座標**下配對。

---

## 3. Step 1–14：原始方法、執行鏈及資料來源（Archive 已核對）

### 3.1 Step 1–21 完整 audit manifest

| Step   | 主要目的                                                                  | Archive 中原始 script／證據                                                                                   | 核對狀態                                                               |
| ------:| --------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ |
| **1**  | 解析 ClinVar GRCh38 missense SNV、P/LP vs B/LB、review stars              | `01_count_clinvar_stars.py`、log、`clinvar_20261025_{stars_summary.tsv,missense_variants.tsv.gz,qc.json}` | ✅ 腳本／輸出數字交叉驗證；原始 VCF 未入 ZIP                                        |
| **2**  | 初期 VEP AF 註解（gnomADe/gnomADg、dbNSFP joint）                            | `02_annotate_gnomad.sh`、test/full logs                                                                  | ✅ 命令與完成 log；完整註解 TSV 未入 ZIP                                        |
| **3**  | 比較 VEP exome/genome 與 dbNSFP joint AF 的覆蓋率、閾值一致性                      | `03_compare_gnomad_sources.py`、log、`03_af_coverage_full.tsv`                                            | ✅ 腳本／表／log 一致                                                      |
| **4**  | unresolved non-missense variants 中間抽樣／VCF                             | `04_unresolved_nonmissense.vcf`                                                                         | ⚠️ 有 5,499 筆 VCF 記錄，**沒有 Step 4 script/log，不能還原產生規則**              |
| **5**  | VEP v115 all-transcript annotation                                    | `05_transcript_vep.sh`、log、`05_clinvar_all_transcripts.tsv.gz`                                          | ✅ 腳本／輸出；後續 Step 7 讀入 4,456,611 transcript rows                     |
| **6**  | 可能的中間人工檢查／整併階段                                                        | ZIP 未見 `06_*` script、log、results                                                                        | ⚠️ **無直接紀錄，不能推定已執行獨立 Step 6**                                      |
| **7**  | 整合 ClinVar VCF、VEP transcript、初期 AF 成 variant master                  | `07_make_clinvar_master.py`、log、`07_*`                                                                  | ✅ 226,010 rows；primary ≥2 stars 59,642；89 筆 canonical fallback     |
| **8**  | 比較十種 reference cohort selection                                       | `08_summarize_clinvar_cohorts.py`、log、`08_*`                                                            | ✅ 10 組 counts 可追溯                                                  |
| **9**  | 在 primary missense ≥2 stars 變異內比較 AF 缺失與 discordance                  | `09_compare_strict_missense_af.py`、log、`09_*`                                                           | ✅ 主要 AF 缺失數字可追溯                                                    |
| **10** | 對 gnomAD v4.1 joint MANE-filtered VCF 精確等位基因配對                        | `10_annotate_gnomad_joint.sh`（初試失敗）、**`10b_annotate_gnomad_joint.py`（正式）**、logs                         | ✅ 修正方案可追溯；外部 VCF 原件未入 ZIP                                          |
| **11** | 將新 joint AF 合併回 master，移除舊 AF 欄位                                      | **`11b_merge_joint_master.py`**、log、`11_*`                                                              | ✅ primary 59,642、status/AF summary 一致                              |
| **12** | 以 exome + genome all-sites AN BED 檢視 site-level callable allele count | `12_annotate_all_sites_AN.py`、test/full logs、`12_*`                                                     | ✅ summary、226,010 rows、214,472 unique loci 的 log 可追溯；原始 BED 未入 ZIP |
| **13** | 查詢 joint subset VCF，分辨 same-locus other ALT vs no variant record      | `13_audit_unobserved_loci.py`、log、`13_*`                                                                | ✅ 10,878 audit rows、0 unexpected exact ALT                         |
| **14** | 產生三組 nested autosomal cohorts                                         | `14_build_calibration_cohorts.py`、log、`14_*`                                                            | ✅ 59,642 variant keys／三組總數與巢狀關係獨立核對                                |
| **15** | 本機 Parquet score annotation（V100/M447）                                | `15_annotate_gpnstar_local.py`、test/full logs、cache、`15_*`                                              | ✅ 59,641 paired-score candidates；原始 Parquet 未入 ZIP                 |
| **16** | V100/M447 variant-key paired merge                                    | `16_merge_gpnstar_models.py`、log、`16_*`                                                                 | ✅ 59,642 keys、59,641 paired present、1 MT                           |
| **17** | 3 cohorts × 2 models × signed/absolute scores 的 AUC/AP                | `17_compare_gpnstar_auc.py`、log、`17_*`                                                                  | ✅ **12 組 ROC AUC 從原始 paired rows 獨立重算一致**                          |
| **18** | Gene annotation QC、gene composition                                   | `18_gene_composition_qc.py`、log、`18_*`                                                                  | ✅ 5,992 genes、3 類組成獨立重算一致                                          |
| **19** | Gene-disjoint calibration/validation 80/20                            | `19_gene_disjoint_split.py`、log、`19_*`                                                                  | ✅ variant/gene counts、0 gene overlap、score consistency             |
| **20** | 1-point empirical LR profile；500 gene-cluster bootstrap               | `20_v100_local_lr_profile.py`、log、`20_*`                                                                | ✅ 36 列 **實際重跑逐欄一致**                                                |
| **21** | 七個 candidate bins、2,000 gene-cluster bootstrap                        | `21_v100_evidence_interval_qc.py`、log、`21_*`                                                            | ✅ 14 列 **實際重跑逐欄一致**；Strong 仍待獨立確認                                  |

> 請不要將編號 Step 4 或 Step 6 寫成已有完整可重現工作流程：ZIP 中沒有對應的程式／指令。本紀錄按**實際封存內容**逐項重建，而非由最終結果反推尚未提供的操作。

### 3.2 Step 1 — ClinVar 初始解析

**來源／輸入**：`clinvar_vcf/clinvar_20261025.vcf.gz`（GRCh38；原始 VCF 未打包）；`scripts/01_count_clinvar_stars.py`。

- 逐列讀 VCF，僅保留 `REF` 與 `ALT` 均為單一 A/C/G/T 的 biallelic SNV；初步 missense 依 ClinVar INFO/`MC` **包含 `missense_variant`** 過濾。這只是*初篩*，尚未確認 VEP 的 gene-matched transcript missense；該部分在 Step 5–7 解決。
- 嚴格 CLNSIG 分類：`pathogenic`、`likely_pathogenic`、`pathogenic/likely_pathogenic` → P/LP；`benign`、`likely_benign`、`benign/likely_benign` → B/LB；其餘組合（包括 VUS/conflicting）未納入。
- `CLNREVSTAT` 轉換為 0–4 stars；`criteria_provided,_single_submitter`=1；`criteria_provided,_multiple_submitters,_no_conflicts`=2；`reviewed_by_expert_panel`=3；`practice_guideline`=4。最終主要 cohort 只採 **≥2**。
- **重要歷史欄位**：Step 1 輸出 `variation_id` 曾暫填 ClinVar INFO/`ALLELEID`；**Step 7 從原始 VCF ID 欄位重新取回真正的 Variation ID，並分別保留 Allele ID**。不能把 Step 1 的 `variation_id` 直接當最終 Variation ID。

| 指標                                            | 原始 log／QC 結果 |
| --------------------------------------------- | ------------:|
| Total ClinVar VCF records                     | 4,599,696    |
| SNV records                                   | 4,269,569    |
| MC 包含 missense 的 SNVs                         | 2,581,130    |
| CLNSIG 屬 P/LP 或 B/LB 且成功納入                    | **226,010**  |
| ≥2 stars P/LP                                 | **25,971**   |
| ≥2 stars B/LB                                 | **43,022**   |
| ≥2 stars 兩類合計（尚未做 gene-matched transcript QC） | **68,993**   |

Stars 分布（全部 226,010 筆）如下：

| Stars | P/LP   | B/LB    | 合計      |
| -----:| ------:| -------:| -------:|
| 0     | 10,427 | 6,158   | 16,585  |
| 1     | 36,193 | 104,239 | 140,432 |
| 2     | 22,478 | 41,482  | 63,960  |
| 3     | 3,485  | 1,540   | 5,025   |
| 4     | 8      | 0       | 8       |

**核對**：封存 `clinvar_20261025_missense_variants.tsv.gz` 確實 226,010 筆、variant keys 唯一，且自行從 stars/label 重新計算 ≥2 stars 為 25,971／43,022；原始 VCF 的 4,599,696 計數僅能核對 log/JSON，未重新掃描 VCF。

### 3.3 Step 2 — 初期 gnomAD AF 註解（後續被 Step 10–11 取代）

**腳本**：`scripts/02_annotate_gnomad.sh`（test/full logs 都存在）。

- 將 Step 1 的 226,010 筆轉成 `gnomad_annotation/clinvar_missense.vcf`；前 100 筆轉成 `test100.vcf`。
- 使用 VEP v115/GRCh38 cached environment，輸出 `gnomADg_AF`、`gnomADg_EAS_AF`、`gnomADe_AF`、`gnomADe_EAS_AF`。若本機 `dbNSFP5.3a_with_pknn_grch38.gz` 可讀，加入 `gnomAD4.1_joint_AF`、`gnomAD4.1_joint_EAS_AF` 的 dbNSFP plugin annotation。
- logs 顯示 `Using dbNSFP 5.3a joint AF`；test 及 full 皆成功完成，full 輸出 `gnomad_annotation/clinvar_missense_gnomad.tsv`（該 TSV 未打包）。
- **重要區別**：此處的 `gnomAD4.1_joint_AF` 是初期 *dbNSFP annotation*；**不是** Step 10b 從 local VCF 的 `AC_e/AN_e/AC_g/AN_g` 重新推導的 final joint AF，兩者不可混用。

### 3.4 Step 3 — AF coverage／來源比較

**腳本／結果**：`03_compare_gnomad_sources.py`、`03_af_coverage_full.tsv`。完整執行分析 226,010 筆（VEP annotation rows 226,010）。

在 **≥2 stars** 但尚未 final gene-matched 過濾的子集：

| label | source       | n      | AF 有值  | AF coverage |
| ----- | ------------ | ------:| ------:| -----------:|
| P/LP  | dbNSFP joint | 25,971 | 11,810 | 45.47%      |
| P/LP  | VEP exomes   | 25,971 | 12,804 | 49.30%      |
| P/LP  | VEP genomes  | 25,971 | 6,853  | 26.39%      |
| B/LB  | dbNSFP joint | 43,022 | 38,617 | 89.76%      |
| B/LB  | VEP exomes   | 43,022 | 42,644 | 99.12%      |
| B/LB  | VEP genomes  | 43,022 | 41,782 | 97.12%      |

**研究意義**：P/LP 的 AF missingness 大於 B/LB；只採已取得 AF 的 variant 會系統性排除大量 P/LP，促成後來 all-sites AN 與 putatively-unobserved cohort 設計。AF 來源在 <0.001 處會產生 discordant classification，不應預設 dbNSFP joint 與 exome/genome AF 可互換。

### 3.5 Step 4 — unresolved non-missense VCF（追溯受限）

Archive 有 `results/04_unresolved_nonmissense.vcf`：2 行 header + **5,499 個 VCF variant records**（`wc -l` 為 5,501）。這是可確認的中間產物，但 `scripts/` 與 `logs/` **都沒有 `04_*`**。

- **可確認**：檔案存在、VCF 格式、記錄數量。
- **不可確認**：如何定義 unresolved、抽樣自何處、使用哪些過濾條件、是否曾用來更新主資料。
- 目前 Step 7 script 直接使用 Step 1／2／5 的產物，**沒有直接讀取 Step 4 VCF**；因此不能聲稱 Step 4 對 final reference cohort 有任何已證實的影響。

### 3.6 Step 5 — VEP all-transcript annotation

**腳本**：`05_transcript_vep.sh`。使用 `vep_115.sif`、GRCh38 offline cache，參數包括 `--tab --symbol --biotype --mane --canonical`，刻意**沒有** `--pick`，保留 all-transcript annotation；輸出欄位 `Uploaded_variation,Gene,SYMBOL,Feature,BIOTYPE,Consequence,MANE_SELECT,MANE_PLUS_CLINICAL,CANONICAL`。

**結果**：`results/05_clinvar_all_transcripts.tsv.gz`；`05_transcript_vep.log` 顯示輸出完成，**Step 7 載入 4,456,611 transcript rows**。該數量是 transcript annotations，不是不同的 genomic variant 數。

### 3.7 Step 6 — 無獨立可核對工作紀錄

封存的 `scripts/`、`logs/`、`results/` 均無 `06_*`，Step 7 的輸入又明確指向 Step 1／2／5。故不能確認是否真的執行一個編號為 6 的中間步驟。保留編號但標為**未提供執行依據**；未來若找到當時命令歷史或檔案才補入，不要為了湊齊 21 步而杜撰。

### 3.8 Step 7 — 建立 final ClinVar/VEP annotation master

**腳本／產物**：`07_make_clinvar_master.py`、`07_make_clinvar_master.log`、`07_clinvar_master_annotation.tsv.gz`、`07_master_qc.json`。以 `chrom:pos:ref:alt` matching，使用 SQLite 暫存 4,456,611 transcript rows，避免在記憶體載入所有 transcripts；保留 `all_transcripts_json`，一個 genomic variant 一列，總計 **226,010**。

真正的 `primary_missense` 定義（已直接確認 source code）：

1. ClinVar GENEINFO 必須**恰為單一 gene** 並能與 VEP `SYMBOL` matched；且 VEP transcript 的 `BIOTYPE=protein_coding`。
2. 至少一條對應 **MANE Select 或 MANE Plus Clinical** transcript 是 `missense_variant`；**或** 該 gene 完全沒有可用 MANE transcript，才允許使用 canonical transcript 的 missense 作 fallback。
3. 單純其他 transcript 是 missense，或 ClinVar 複數基因需人工判定者，**不能**自動成為 `primary_missense`。
4. **`strict_mane_missense` 與 `primary_missense` 不同**：Step 8–14 最終選樣使用後者，包含少量 canonical fallback。於 ≥2 stars 最終的 59,642 筆中，**89 筆為 canonical fallback**，另 **59,553 筆**以符合規則的 MANE missense 路徑納入。請避免 Methods 寫成「全部 MANE-only」。

QC：`star2_P/LP_primary_missense=22,589`、`star2_B/LB_primary_missense=37,053` → **59,642**。封存 master 已獨立掃描，226,010 rows、primary ≥2 兩類 counts 完全吻合，primary 變異 0 個 gene-match 違反、0 個 transcript rule 違反。

### 3.9 Step 8 — 十種 reference cohorts 的初始比較

**腳本／結果**：`08_summarize_clinvar_cohorts.py`、`08_cohort_summary.tsv`、`08_primary_rare_missense.tsv.gz`。最初偏重 **VEP exome AF <0.001** 的 rare selection，將 AF missingness 與其他篩選方式作比較。

| Cohort（原始 key）                           | P/LP       | B/LB       | Total      |
| ---------------------------------------- | ----------:| ----------:| ----------:|
| `01_primary_ge2_exome_0.001`             | 11,044     | 14,378     | 25,422     |
| `02_mane_only_ge2_exome_0.001`           | 11,029     | 14,355     | 25,384     |
| `03_broad_gene_missense_ge2_exome_0.001` | 11,087     | 14,890     | 25,977     |
| `04_ge1_exome_0.01`                      | 16,477     | 99,136     | 115,613    |
| `05_ge2_exome_0.01`                      | 11,115     | 27,059     | 38,174     |
| `06_ge2_exome_0.0001`                    | 10,441     | 6,218      | 16,659     |
| `07_ge3_exome_0.001`                     | 1,462      | 1,008      | 2,470      |
| `08_ge2_joint_0.001`                     | 10,317     | 17,884     | 28,201     |
| `09_ge2_all_AF`                          | **22,589** | **37,053** | **59,642** |
| `10_ge2_missing_exome`                   | 11,471     | 200        | 11,671     |

**注意**：名稱 `09_ge2_all_AF` 實際定義是 *≥2 stars + primary_missense、不套 AF cutoff*，不是所有變異都有 AF。`08_primary_rare_missense.tsv.gz` 是早期比較的 rare subset，**不是**後來 Step 14 的 expanded primary cohort。

### 3.10 Step 9 — primary missense cohort 的 AF missingness audit

**腳本**：`09_compare_strict_missense_af.py`；以 ≥2 stars 且 `primary_missense=1` 的 59,642 筆比較 VEP exome/genome 與 dbNSFP joint。

| label | source       | n      | AF 有值  | AF 缺失  | coverage |
| ----- | ------------ | ------:| ------:| ------:| --------:|
| P/LP  | VEP exomes   | 22,589 | 11,118 | 11,471 | 49.22%   |
| P/LP  | VEP genomes  | 22,589 | 6,013  | 16,576 | 26.62%   |
| P/LP  | dbNSFP joint | 22,589 | 10,377 | 12,212 | 45.94%   |
| B/LB  | VEP exomes   | 37,053 | 36,853 | 200    | 99.46%   |
| B/LB  | VEP genomes  | 37,053 | 36,159 | 894    | 97.59%   |
| B/LB  | dbNSFP joint | 37,053 | 34,905 | 2,148  | 94.20%   |

在 VEP exomes 與 dbNSFP joint **兩者均有 AF** 的 P/LP 9,748、B/LB 34,870 筆，AF<0.001 cutoff discordance 分別 2/9,748（0.02%）與 1,968/34,870（5.64%）。這是**不同 AF 來源的初期比較**，不是 Step 10 final joint AF 的一致性結果。

### 3.11 Step 10 / 10b — gnomAD v4.1 local joint AF 的修正流程

- **正式 source**：`gnomad_raw/joint_mane_v41/variantsGnomADv4.1Joint.vcf.bgz`（MANE-restricted 第三方衍生 source，非官方完整 variants VCF），有 `.tbi`；ZIP 不包含原始 VCF，不能從 Archive 再檢查下載 checksum 或被濾掉的 records。
- `10_annotate_gnomad_joint.sh` 初次測試期待 `AC_j/AN_j/AF_j`，log 明確失敗 `ERROR: Joint AC/AN/AF tags not found`；不能把這個初始腳本寫成最終 successful method。
- 正式跑成功的是 **`10b_annotate_gnomad_joint.py`**，讀取 `AC_e,AN_e,AC_g,AN_g`，對相符 ALT 以 joint `AC=(AC_e+AC_g)`、`AN=(AN_e+AN_g)`、`AF=AC/AN` 推導；EAS 同理。保留原 source FILTER／multi-allelic ALT index。若 AC/AN 所需 component 缺失，記為 incomplete，而非零。
- Test：matched 99/100、1 unresolved；Full：**226,010** queried，**166,705 exact matched**、**59,305 not_found_unresolved**；matched 中 joint AF calculated **166,661**、missing components **44**、PASS 有有效 AF **157,995**。
- `FILTER_PASS=158,032` 是 raw FILTER 分類數；`joint_AF_PASS=157,995` 是**PASS 且有有效 joint AF**，兩者分母／定義不同。計算時要保留此區別。

### 3.12 Step 11 / 11b — joint AF 合併回 master

**實際 script 名稱是** `11b_merge_joint_master.py`。將 Step 10b 的 `clinvar_full_joint.tsv` 回併 226,010-row master，**移除舊初期 VEP/dbNSFP AF 欄位**，保留新 `gnomad41_*`、status 等欄位，避免後續混讀不一致來源。

≥2 stars `primary_missense` = 59,642 筆（P/LP 22,589；B/LB 37,053）：

| Final joint status        | P/LP   | B/LB   | Total      |
| ------------------------- | ------:| ------:| ----------:|
| `observed_PASS`           | 10,892 | 36,450 | 47,342     |
| `observed_FILTERED`       | 983    | 435    | 1,418      |
| `matched_incomplete`      | 0      | 4      | 4          |
| `not_found_unresolved`    | 10,714 | 164    | 10,878     |
| `observed_PASS` 且 AF<0.01 | 10,890 | 30,161 | **41,051** |
| `observed_PASS` 且 AF≥0.01 | 2      | 6,289  | **6,291**  |

`AF<0.01` 只套在 **observed_PASS**，並非把所有缺失 AF 變異當成頻率 0。

### 3.13 Step 12 — gnomAD all-sites AN annotation

**腳本**：`12_annotate_all_sites_AN.py`。來源為 exomes/genomes gnomAD v4.1 all-sites AN BED（衍生 MANE-restricted tracks；原始 BED 不在 ZIP）。先將 226,010 SNV 組成 **214,472** 個 distinct genomic loci，再以 0-based BED `[POS−1,POS)` 與 tabix 查詢：exome 213,907 個 loci 有 AN；genome 213,945 個 loci 有 AN；每個 variant 保存 `AN_e_all_sites`、`AN_g_all_sites` 及其和（不是 ALT AC）。

≥2 stars primary cohort：

| AN QC                          | P/LP       | B/LB    | Total      |
| ------------------------------ | ----------:| -------:| ----------:|
| `not_found_unresolved`         | 10,714     | 164     | **10,878** |
| ALT absent、exome/genome AN 均有值 | 10,703     | 145     | 10,848     |
| ALT absent、部分或零 AN             | 10         | 10      | 20         |
| ALT absent、無 AN record         | 1          | 9       | 10         |
| 不限染色體，joint AN≥1,000,000       | **10,672** | **123** | **10,795** |

**重要技術限制**：gnomAD v4.1 官方 all-sites AN 是 genomic locus 的 callable allele number，可評估「此位置有多少 reference/alternate allele calls」；不能從 `AN≥1,000,000` 直接得到目標 ALT 的 `AC=0`。此外 all-sites AN 與 variant VCF 的 allele-specific AN 在某些 multiallelic loci 會不同，應保留來源差異。

### 3.14 Step 13 — unresolved ALT 的 locus-level audit

**腳本**：`13_audit_unobserved_loci.py`。對 Step 11/12 在 source VCF exact ALT 未匹配的 **10,878** variants 建立 0-based query BED，再以 `bcftools view -H -R` 回查**同一份第三方 MANE subset VCF**；區分沒有任何 record、同位置別的 ALT，以及不應出現的 exact ALT。

| Locus audit category            | P/LP  | B/LB | Total     |
| ------------------------------- | -----:| ----:| ---------:|
| `no_variant_record_at_position` | 7,872 | 110  | 7,982     |
| `other_variant_at_position`     | 2,842 | 54   | 2,896     |
| `UNEXPECTED_EXACT_ALT`          | 0     | 0    | **0**     |
| Autosomal（不論 AN）                | 9,544 | 150  | 9,694     |
| Autosomal AN≥1M                 | 9,504 | 119  | **9,623** |
| ├ 高 AN 且無任何 variant record      | 6,884 | 70   | 6,954     |
| └ 高 AN 且同位置有其他 ALT              | 2,620 | 49   | 2,669     |

這個「0 unexpected exact ALT」**只表示 source subset 的兩次查詢一致**，不是對官方全量 gnomAD 的 ALT absence proof。此點會影響最終 LR 的 generalizability。

### 3.15 Step 14 — 建立正式三組 autosomal reference cohorts

**腳本**：`14_build_calibration_cohorts.py`，輸出 `14_calibration_cohort_flags.tsv.gz`／summary。第一階段只留 ≥2 stars 且 `primary_missense=1` 的 **59,642** variants。`chromosome_class=autosome` 限 chr1–22。

巢狀定義：

1. **Observed-only**：原 source exact ALT 為 `observed_PASS`，且 derived joint AF **<0.01**。
2. **Conservative**：observed-only + `not_found_unresolved`、autosomal joint all-sites AN≥**1,000,000**、且 subset VCF **同 locus 無任何 variant record**。
3. **Expanded（主要校正集合）**：conservative + high AN、same locus 有**其他 ALT**，但找不到目標 ALT。

| Cohort                             | P/LP   | B/LB   | Total      | Nested selection                                  |
| ---------------------------------- | ------:| ------:| ----------:| ------------------------------------------------- |
| `autosomal_observed_only`          | 10,585 | 28,063 | **38,648** | Base                                              |
| `autosomal_no_other_ALT_candidate` | 17,469 | 28,133 | **45,602** | Base + 6,954 high-AN candidates                   |
| `autosomal_expanded_candidate`     | 20,089 | 28,182 | **48,271** | Conservative + 2,669 high-AN other-ALT candidates |

Expanded 內 `observed_PASS_rare`=38,648、`highAN_no_variant_record`=6,954、`highAN_other_ALT_at_site`=2,669；合計 **48,271**。本次獨立 audit 確認 59,642 個 variant_key 唯一、三組 P/LP/B/LB 計數一致、**零個違反 observed ⊆ conservative ⊆ expanded**。`observed_FILTERED`、AF≥0.01、sex/MT 或 AN 不足者均未直接進主要 autosomal cohort。

依後續 Step 15/16 scores，三類 `S=-V100 llr_calibrated` 中位數分別：

| Category                   | P/LP   | B/LB   | P/LP median S | B/LB median S |
| -------------------------- | ------:| ------:| -------------:| -------------:|
| `observed_PASS_rare`       | 10,585 | 28,063 | 8.3790        | 3.8220        |
| `highAN_no_variant_record` | 6,884  | 70     | 10.4580       | 2.3765        |
| `highAN_other_ALT_at_site` | 2,620  | 49     | 9.8360        | 2.3730        |

**研究決策／透明性**：選擇 expanded 作 primary，是為避免有充分 AN、但 source subset 沒有 exact ALT 的大量罕見 P/LP 被排除；選擇已在探索性 AUC 結果被查看後形成，**並非預先註冊的 reference-set choice**。Expanded vs observed LR 的變化可能包含 reference selection/AF bias，Step 22 必須保留敏感度分析。

---

## 4. Step 15–16：GPN-Star annotation 與 paired merge

### Step 15 — V100/M447 本機 Parquet annotation

- **腳本**：`scripts/15_annotate_gpnstar_local.py`；先以 chr21 100 筆／模型測試，兩套皆 **100/100** 非 null signed scores。
- Full mode 掃描 chr1–22、X、Y 對應 Parquet，按 `pos/ref/alt` 合併；per-chromosome caches 可支援續跑。
- Full QC：兩套模型皆 `score_present=59,641`、`score_missing=0`、`unsupported_chromosome=1`。
- 唯一不支援者為 **`MT:12706:T:C`（P/LP）**。MT 不在研究 primary autosomal cohort，不影響主要分析。
- 三個 autosomal cohorts 的 V100/M447 score coverage 均 **100%**。
- **輸出**：
  - `results/15_gpnstar_v100_full.tsv.gz`
  - `results/15_gpnstar_m447_full.tsv.gz`
  - 各模型 `15_gpnstar_{model}_qc_full.tsv`

### Step 16 — paired merge

- **腳本**：`scripts/16_merge_gpnstar_models.py`。
- 依 `variant_key` merge，不以列順序 merge；檢查來源 metadata 一致、無重複鍵、分數 finite。
- 加入 `gpnstar_v100_deleteriousness = -gpnstar_v100_llr_calibrated`、M447 對應欄位及 `paired_score_status`。
- QC：**59,642** keys 全部配對；**59,641** 雙模型有分數；**1** 筆 MT 無分數；三個 autosomal cohorts paired missing **0**。
- **輸出**：`results/16_gpnstar_v100_m447_paired.tsv.gz`、`results/16_gpnstar_paired_coverage.tsv`。

---

## 5. Step 17：V100 vs M447 predictive discrimination

**Archive 追加核對**：本次從 `16_gpnstar_v100_m447_paired.tsv.gz` 的 59,642 rows 重新以 Mann–Whitney U/rank statistic 計算全部 12 組 AUC（3 cohorts × 2 scores × 2 models），與 `17_gpnstar_auc_summary.tsv` **浮點數值完全一致**；observed-only 三組 AF strata 的 P/LP、B/LB counts 也逐一吻合。Bootstrap CI 為 original Step 17 的 500 次 *variant-level stratified paired bootstrap*，不能視為已控制 gene clustering。

- **腳本**：`scripts/17_compare_gpnstar_auc.py`；三個 cohorts × 兩種 score（signed、absolute）× 兩模型。
- AUC 使用 Mann–Whitney／rank-based 計算；AP 使用 average precision；對同一批 variants 使用 **500 次 stratified paired variant bootstrap** 估 AUC 95% CI 及 paired ΔAUC。
- 此時 CI 未以 gene clustering 處理同基因相依性；屬探索性 discrimination。

### 5.1 Signed score 的主要結果

| Cohort        | V100 AUC（95% variant-bootstrap CI） | M447 AUC | Paired ΔAUC V100−M447（95% CI） | V100 AP | M447 AP |
| ------------- | ---------------------------------- | --------:| ----------------------------- | -------:| -------:|
| Observed-only | **0.8898** (0.8862–0.8931)         | 0.8652   | **+0.0247** (0.0225–0.0270)   | 0.7286  | 0.6507  |
| Conservative  | **0.9184** (0.9160–0.9209)         | 0.8925   | **+0.0259** (0.0239–0.0278)   | 0.8719  | 0.8164  |
| Expanded      | **0.9218** (0.9194–0.9240)         | 0.8964   | **+0.0254** (0.0236–0.0273)   | 0.8903  | 0.8426  |

V100 absolute scores 的 AUC 依序為 **0.8796、0.9129、0.9176**；M447 absolute 為 **0.8642、0.8919、0.8959**。V100 signed 在三組均較其 absolute 分數略好，亦均優於 M447 signed；主要 predictive comparison 方向一致。

### 5.2 Observed-only 的 AF-stratified AUC（signed）

| Joint AF stratum       | P/LP   | B/LB   | V100 AUC   | M447 AUC |
| ---------------------- | ------:| ------:| ----------:| --------:|
| AF < 1×10^-4           | 10,061 | 6,426  | **0.8472** | 0.8054   |
| 1×10^-4 ≤ AF < 1×10^-3 | 464    | 10,496 | **0.8553** | 0.8366   |
| 1×10^-3 ≤ AF < 0.01    | 60     | 11,141 | **0.8443** | 0.8077   |

**解讀**：整體 observed-only AUC 0.8898 高於各 AF stratum 內的 AUC（約 0.84–0.86），提示 **AF-related confounding／case-mix effects**；分層結果支持 GPN-Star 仍有獨立於粗略 AF strata 的區辨訊號，但尚不能量化因果貢獻。高 AF stratum P/LP 只有 60 筆。各 cohort 的 P/LP proportion 不同，**不可直接用 AP 上升論證模型品質提高**。

- **輸出**：`17_gpnstar_auc_summary.tsv`、`17_gpnstar_paired_auc_delta.tsv`、`17_gpnstar_af_stratified_auc.tsv`。

---

## 6. Step 18：gene annotation 與 gene composition QC

- **腳本**：`scripts/18_gene_composition_qc.py`。
- 從 Step 12 master 的 `clinvar_genes`、`matched_vep_genes` 回接 Step 16 expanded variants；核對 gene-match 和單一 gene attribution。
- **48,271/48,271** expanded variants：`gene_match_exact`；多基因、部分重疊、缺失、gene 不一致皆為 **0**。
- **5,992** unique genes：
  - P/LP 與 B/LB 均有：**1,527**。
  - 僅 P/LP：**322**。
  - 僅 B/LB：**4,143**。
- Gene counts 分布不均，示例：`FBN1` P/LP 569、B/LB 20；`LDLR` 485/11；`PAH` 384/1；`KMT2D` 24/254。代表 label composition、disease mechanism 與基因別樣本量不均；不可單純把所有變異視為獨立同分布。
- **輸出**：`18_expanded_gene_counts.tsv`、`18_expanded_gene_annotation_audit.tsv.gz`、`18_expanded_gene_qc_summary.tsv`。

---

## 7. Step 19：固定 gene-disjoint calibration / validation

**Archive 追加核對**：`19_expanded_gene_split.tsv.gz` 中 48,271 個 variant_key 均唯一，與 Step 16 相應 score 完全一致；自行計數 5,992 genes，**0 個基因跨 calibration/validation**。`19_gene_split_manifest.tsv` 保存基因層級 assignment。此 split 從 2,000 次 candidate gene-stratified random partitions 中，以 *P/LP 和 B/LB variant counts* 與 20% 目標的偏差選擇，**沒有使用 GPN-Star score 作選分割依據**。但選取平衡 label composition 不代表已針對 gene-disease mechanism 或 AF 分布分層。

- **腳本**：`scripts/19_gene_disjoint_split.py`。
- 使用 `clinvar_genes` 為 gene grouping（Step 18 已確認每筆均對單一基因）。同一基因的所有 variants 分配到同一 split。
- 三層 gene composition：both／P/LP-only／B/LB-only；以 `seed=20261029`、`2000` candidate splits、`validation-fraction=0.20` 尋找 P/LP 與 B/LB counts 最接近 80/20 的 assignment，**不使用 GPN-Star score/AUC 來選擇 split**。
- 固定後對 expanded、conservative、observed-only 套用**相同 gene assignments**，gene leakage QC passed。

| Cohort        | Split       | P/LP       | B/LB       | Total      | Genes     | P/LP fraction |
| ------------- | ----------- | ----------:| ----------:| ----------:| ---------:| -------------:|
| Expanded      | Calibration | **16,079** | **22,548** | **38,627** | **4,794** | 41.63%        |
| Expanded      | Validation  | **4,010**  | **5,634**  | **9,644**  | **1,198** | 41.58%        |
| Conservative  | Calibration | 14,008     | 22,503     | 36,511     | 4,779     | 38.37%        |
| Conservative  | Validation  | 3,461      | 5,630      | 9,091      | 1,196     | 38.07%        |
| Observed-only | Calibration | 8,608      | 22,442     | 31,050     | 4,718     | 27.72%        |
| Observed-only | Validation  | 1,977      | 5,621      | 7,598      | 1,175     | 26.02%        |

- Expanded validation 約占所有 P/LP 的 **19.96%**、所有 B/LB 的 **19.99%**。
- **輸出**：`results/19_expanded_gene_split.tsv.gz`、`19_gene_split_manifest.tsv`、`19_gene_split_summary.tsv`。
- **重要**：此時之後不能依 AUC／LR 重新調整 split；此為 gene-disjoint **內部驗證**，不是完全外部、獨立的新 ClinVar/clinical reference set。

---

## 8. Step 20：V100 empirical local LR profile（探索性 1-point bins）

**Archive 重跑驗證**：用封存版 `20_v100_local_lr_profile.py`、相同固定 seed `20261029`、500 次 bootstrap，在隔離副本重新執行；36 列、全部欄位與封存 `20_v100_local_lr_profile.tsv` **完全一致**。這是「同一程式可重現」，不是獨立統計方法的外部交叉驗證。

- **腳本**：`scripts/20_v100_local_lr_profile.py`；只使用 expanded，分別分析 calibration 與 validation。
- 定義 \(S=-\text{V100 llr_calibrated}\)；以整數區間（如 `8<S<=9`）建立每 1 分寬的 bins，加上兩端 tail。
- LR 由各 split 內 P/LP／B/LB **class-conditional bin frequency 比率**取得；**500 次 gene-cluster bootstrap** 估 CI；原腳本要求每區間每類至少 20 筆才標 `estimable`，少數樣本為 `sparse_or_unstable`。

| V100 S bin | Calibration P/B | Calibration LR | Validation P/B | Validation LR |
| ---------- | ---------------:| --------------:| --------------:| -------------:|
| 4<S≤5      | 318 / 2,804     | 0.159          | 79 / 685       | 0.162         |
| 5<S≤6      | 768 / 2,368     | 0.455          | 153 / 630      | 0.341         |
| 6<S≤7      | 1,387 / 2,086   | 0.932          | 289 / 524      | 0.775         |
| 7<S≤8      | 2,129 / 1,553   | 1.922          | 461 / 382      | 1.696         |
| 8<S≤9      | 2,628 / 1,001   | 3.682          | 658 / 280      | 3.302         |
| 9<S≤10     | 2,751 / 558     | 6.914          | 750 / 146      | 7.217         |
| 10<S≤11    | 2,557 / 248     | 14.459         | 669 / 59       | 15.931        |
| 11<S≤12    | 1,960 / 69      | 39.834         | 541 / 13       | 58.469        |
| 12<S≤13    | 1,146 / 10      | 160.707（稀疏）    | 325 / 3        | 152.207（稀疏）   |

**觀察**：LR 隨 S 增加有明確、幾乎單調的變化；calibration/validation 趨勢近似。對非常高或低分的 tails，少數類別樣本量不足，點估計不宜直接轉成 evidence strength。此步驟只是探索性 profile，尚無正式 PP3/BP4 cutoff。

- **輸出**：`results/20_v100_local_lr_profile.tsv`。

---

## 9. Step 21：候選 evidence intervals 與 2,000 次 gene-cluster bootstrap

**Archive 重跑驗證**：用封存版 `21_v100_evidence_interval_qc.py`、相同 seed／2,000 次 gene-cluster bootstrap 重新執行；14 列所有欄位（包括 bootstrap CI）與封存表 **完全一致**。此外從 Step 19 variants **獨立按 score bins 計算** 14 個 P/LP、B/LB counts 與 LR，亦完全吻合。

**CI 解讀修正**：這裡是 **two-sided 95% percentile gene-cluster bootstrap CI**；Pejaver et al. 原始 threshold algorithm 採 one-sided 95% confidence bounds 與連續 score 附近的 local posterior/LR，兩者不能視為相同校正方法。Step 21 的 `candidate` labels 只是以當前分箱及 CI 與 LR 門檻比較，**不代表已按 Pejaver 的完整 procedure 完成 clinical calibration**。

### 9.1 固定比較的七個區間

從 Step 20 的已觀察分布選取七個**互不重疊**的 intervals：

1. `S<=2`
2. `2<S<=5`
3. `5<S<=6`
4. `6<S<=8`
5. `8<S<=9`
6. `9<S<=11`
7. `S>11`

**這組分箱是在先前查看過 calibration 與 validation profile 後決定，因此為 exploratory candidate mapping；不可稱作 untouched-validation、事前預註冊 cutoffs。**

### 9.2 方法

- **腳本**：`scripts/21_v100_evidence_interval_qc.py`。
- 使用 Step 19 已固定的 expanded split；各 split 分別計算 empirical bin LR。
- **2,000 次 gene-cluster bootstrap**：每次以 gene 為單位有放回抽樣，保留該 gene 的所有 variants，對 LR 計算 percentile 95% CI。
- Pathogenic direction 以 **LR CI 下界**對照 Supporting／Moderate／Strong 門檻；benign direction 以 **LR CI 上界**判斷。
- **原程式規則**：只要該 bin 的 P/LP 或 B/LB **<20**，即輸出 `sparse_confirmation`，不論 LR point estimate／CI 多好。Zero-count／無限 LR 也須小心處理。

### 9.3 全部 Step 21 結果

| S interval | Calibration P/LP / B/LB | Calibration LR（95% gene-bootstrap CI） | Validation P/LP / B/LB | Validation LR（95% CI）          | 程式回報（validation）                       |
| ---------- | -----------------------:| ------------------------------------- | ----------------------:| ------------------------------ | -------------------------------------- |
| **S≤2**    | 61 / 6,290              | **0.0136** (0.0092–0.0189)            | 7 / 1,526              | **0.0064** (0.0020–0.0126)     | **Sparse confirmation**（BP4 Strong 候選） |
| **2<S≤5**  | 538 / 8,365             | **0.0902** (0.0734–0.1100)            | 116 / 2,071            | **0.0787** (0.0436–0.1144)     | **BP4 Moderate candidate**             |
| **5<S≤6**  | 768 / 2,368             | 0.4548 (0.3853–0.5339)                | 153 / 630              | 0.3412 (0.2330–0.4501)         | Indeterminate                          |
| **6<S≤8**  | 3,516 / 3,639           | 1.3549 (1.1982–1.5251)                | 750 / 906              | 1.1631 (0.8773–1.4636)         | Indeterminate                          |
| **8<S≤9**  | 2,628 / 1,001           | **3.6816** (3.2448–4.1968)            | 658 / 280              | **3.3017** (2.5968–4.1086)     | **PP3 Supporting candidate**           |
| **9<S≤11** | 5,308 / 806             | **9.2352** (8.2661–10.3567)           | 1,419 / 205            | **9.7253** (8.0974–11.8814)    | **PP3 Moderate candidate**             |
| **S>11**   | 3,260 / 79              | **57.8681** (42.5595–81.2121)         | 907 / 16               | **79.6452** (44.1589–177.2135) | **Sparse confirmation**（PP3 Strong 候選） |

**最具可重現性的候選區間**：

- `2<S<=5`：cal/val CI 上界均低於 benign Moderate LR 0.1727，屬 **BP4 Moderate candidate**。
- `8<S<=9`：cal/val CI 下界分別 3.245、2.597，均高於 PP3 Supporting LR 2.406，屬 **PP3 Supporting candidate**，但 validation margin 較小。
- `9<S<=11`：cal/val CI 下界 8.266、8.097，均高於 PP3 Moderate LR 5.790，屬 **PP3 Moderate candidate**；計數及重現性目前最穩健。
- `S<=2`：calibration 達 BP4 Strong；validation P/LP **僅 7** 筆，標 sparse，**不可現在賦予 BP4 Strong**。
- `S>11`：calibration 達 PP3 Strong；validation B/LB **僅 16** 筆，標 sparse，**不可現在賦予 PP3 Strong**。

**輸出**：`results/21_v100_evidence_interval_qc.tsv`；log：`logs/21_v100_evidence_interval_qc.log`。

---

## 10. 截至 Step 21 的研究結論

1. 已建立 **59,642** 筆嚴格 ClinVar missense SNV reference data，並在 autosomes 形成 observed-only、conservative、expanded nested cohorts。
2. V100 與 M447 score annotation 覆蓋 **59,641/59,642**，唯一缺失為 MT SNV；expanded 的 **48,271** 筆全數有 paired scores。
3. V100 signed score 在三組 cohorts 的 ROC AUC 均高於 M447 signed；expanded AUC 為 **0.9218 vs 0.8964**，配對差 **+0.0254**。但 AUC 差受 variant-level bootstrap 與 reference composition 限制。
4. Expanded 額外 P/LP 的 V100 scores 較高，且整體 observed-only AUC 高於 within-AF AUC，顯示 **reference selection 與 AF confounding** 不可忽略。
5. Expanded 涵蓋 **5,992** genes，單一基因對應 QC 完整；Step 19 已固定 **4,794 gene calibration / 1,198 gene validation**，無 gene leakage。
6. V100 signed score 的 empirical LR 呈明顯遞增關係，且在 held-out genes 中大致重現。
7. 以 Step 21 七區間分箱法，**BP4 Moderate、PP3 Supporting、PP3 Moderate** 在探索性的 gene-disjoint validation 中具有候選支持；Strong tails 因少數類別筆數偏少而尚未確認。
8. **目前沒有正式可直接用於臨床報告的 PP3/BP4 cutoff**；上述候選 evidence mapping 仍須 sensitivity analysis、獨立外部 validation、機轉與 label circularity 評估。

---

## 11. 重要方法學限制與後續補強

| 風險                                       | 為何重要                                                                 | 建議／後續作法                                                                          |
| ---------------------------------------- | -------------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| gnomAD source 不完整                        | Joint MANE-filtered VCF 為第三方 subset；absence 不等於官方 gnomAD AC=0        | 對 expanded 中未匹配 ALT，以官方完整資料獨立交叉驗證，或清楚維持 `putatively_unobserved` 說法               |
| All-sites AN 與 variant-site AN 不等價       | 不同 filtering／variant representation／multi-allelic downcoding 可能影響可比性 | 保留來源欄位及差異 QC；不要直接等同 AF=0                                                         |
| AF／class imbalance                       | P/LP 在超罕見變異集中；AP 對 class prevalence 敏感                               | Step 22 observed-only 與 AF<1e-4 sensitivity；別把不同 cohort 的 AP 變化當作純模型提升           |
| Gene correlation／composition             | 遺傳性疾病 genes 的變異數及 P/LP/B/LB 比例差異甚大                                   | Gene-disjoint split、gene-cluster bootstrap；檢視 top-gene contribution、不同 mechanism |
| Validation 不是完全 untouched                | V100 與 expanded selection 及 Step 21 分箱，均在查看全體或 validation 描述性結果後決定   | 如實稱 exploratory gene-disjoint internal validation；最終門檻需獨立 external validation    |
| ClinVar circularity                      | ClinVar 原始 P/LP、B/LB 判定可能使用 PP3/BP4 或其他計算工具                          | 審查 evidence provenance／classification dependence；盡量使用符合獨立性條件的 reference data     |
| Mechanism-specific validity              | Missense score 對 LOF、GOF、dominant-negative 等機轉可能不同                   | 基因／疾病機轉分層，避免普遍套用全部 missense diseases                                             |
| Local bin LR 與 continuous calibration 不同 | 固定區間 LR 反映 bin frequency ratio；分箱選擇影響證據強度                            | 在 calibration-only 固定門檻後做 robust validation；必要時比較連續 smoothing/local LR           |
| Sparse tails／bootstrap                   | Strong tails 的少數類別非常少；百分位 CI 可能不穩定，數個 genes 即可造成影響                   | 暫不使用 Strong；檢查 per-gene source、zero counts、interval confidence stability，追加獨立資料  |
| 目前僅 autosomal                            | chrX / chrY / MT 的頻率及機制不同                                            | 核基因 X-linked 等另做獨立評估，不直接外推                                                       |

---

## 12. Step 22：當時的原定分析計畫（已於後續執行；結果見第 17 節）

**目的**：在**不修改 Step 19 gene assignment、Step 21 七分數區間**前提下，測試 evidence LR 的來源敏感性。

將使用同一份 `19_expanded_gene_split.tsv.gz`，分別計算 calibration／validation 內：

1. Expanded（primary）。
2. Conservative（排除同 locus 有其他 ALT 的 putatively unobserved subset）。
3. Observed-only（不使用 putatively unobserved ALT）。
4. Observed-only + joint AF < 1×10^-4（降低跨頻率差異）。

**歷史規劃的執行腳本**：`scripts/22_v100_cohort_sensitivity.py`。研究者已於 2026-10-30 提供實際輸出；見第 17 節。

**重點指標**：各 score interval 的 P/LP、B/LB counts、LR、gene-cluster bootstrap 95% CI、P/LP 和 B/LB 的 gene counts、單一 gene 最大貢獻比例與 `sparse` flag。優先檢查 `8<S<=9` 的 PP3 Supporting、`9<S<=11` 的 PP3 Moderate 是否仍通過門檻；再決定是否需要縮小或放寬 strong evidence 的候選範圍。

**原則**：Step 22 僅作 sensitivity analysis；不根據結果重切 genes、不重選 V100 模型、不悄悄調整七個 score bins。

---

## 13. 重要檔案索引（依 ZIP 實際存在檔案，另列未打包來源）

**Scripts（ZIP 實際存在）**：

```text
scripts/
  01_count_clinvar_stars.py
  02_annotate_gnomad.sh
  03_compare_gnomad_sources.py
  05_transcript_vep.sh
  07_make_clinvar_master.py
  08_summarize_clinvar_cohorts.py
  09_compare_strict_missense_af.py
  10_annotate_gnomad_joint.sh       # 初試失敗；非 final
  10b_annotate_gnomad_joint.py      # final gnomAD exact-allele extraction
  11b_merge_joint_master.py        # final merge
  12_annotate_all_sites_AN.py
  13_audit_unobserved_loci.py
  14_build_calibration_cohorts.py
  15_annotate_gpnstar_local.py
  16_merge_gpnstar_models.py
  17_compare_gpnstar_auc.py
  18_gene_composition_qc.py
  19_gene_disjoint_split.py
  20_v100_local_lr_profile.py
  21_v100_evidence_interval_qc.py
  README.txt
```

**Results（ZIP 實際存在）**：

```text
results/
  clinvar_20261025_qc.json
  clinvar_20261025_stars_summary.tsv
  clinvar_20261025_missense_variants.tsv.gz
  03_af_coverage_full.tsv
  04_unresolved_nonmissense.vcf
  05_clinvar_all_transcripts.tsv.gz
  07_clinvar_master_annotation.tsv.gz
  07_master_qc.json
  08_cohort_summary.tsv
  08_primary_rare_missense.tsv.gz
  09_af_coverage_strict_missense.tsv
  11_clinvar_master_gnomad41_joint_mane.tsv.gz
  11_joint_primary_summary.tsv
  12_clinvar_master_gnomad41_AN_full.tsv.gz
  12_AN_primary_qc_full.tsv
  13_unobserved_locus_audit.tsv
  13_unobserved_locus_summary.tsv
  14_calibration_cohort_flags.tsv.gz
  14_calibration_cohort_summary.tsv
  15_gpnstar_{v100,m447}_{test,full}.tsv.gz
  15_gpnstar_{v100,m447}_qc_{test,full}.tsv
  16_gpnstar_v100_m447_paired.tsv.gz
  16_gpnstar_paired_coverage.tsv
  17_gpnstar_auc_summary.tsv
  17_gpnstar_paired_auc_delta.tsv
  17_gpnstar_af_stratified_auc.tsv
  18_expanded_gene_counts.tsv
  18_expanded_gene_annotation_audit.tsv.gz
  18_expanded_gene_qc_summary.tsv
  19_expanded_gene_split.tsv.gz
  19_gene_split_manifest.tsv
  19_gene_split_summary.tsv
  20_v100_local_lr_profile.tsv
  21_v100_evidence_interval_qc.tsv
```

**Logs**：ZIP 含所有上面成功完成的主要 stages logs（含 Step 2 test/full、Step 10 初試失敗、10b test/full、Step 12 test/full、Step 15 test/full、Step 17–21）；**無 Step 4、Step 6 獨立 logs**。並有 Step 15 per-chromosome cache。

**主要原始來源但未打包進 ZIP（不可稱已直接讀取或 checksum 驗證）**：

- `clinvar_vcf/clinvar_20261025.vcf.gz`（ClinVar 2026-10-25 GRCh38 原始 VCF）。
- `gnomad_annotation/clinvar_missense_gnomad.tsv`、`gnomad_annotation/joint_mane_v41/clinvar_full_joint.tsv`（Step 2、Step 10b 的關鍵中間檔）。
- `gnomad_raw/joint_mane_v41/variantsGnomADv4.1Joint.vcf.bgz` 和 index（第三方 MANE subset）。
- `gnomad_raw/all_sites_AN_v41/gnomad.exomes.v4.1.AN.bed.gz`、`gnomad.genomes.v4.1.AN.bed.gz` 和 indexes（MANE-restricted AN tracks）。
- GPN-Star V100/M447 original genome-wide LLR Parquet：`/home/n102968/NGS_UI/biotools/gpn_star/gpn-star-hg38-{v100,m447}-200m/llr/llr_chr*.parquet`。

若要達到 full upstream reproducibility，應另保存 source URL、release version、MD5/SHA256、下載日期、完整 execution commands/container digests，並至少對 primary variant ALT absence 抽樣使用官方全量 gnomAD 独立交叉檢查。

---

## 14. 建議的 Methods／結果論文寫作重點（草稿）

**Methods 核心表述（暫定）**：

> We assembled a GRCh38 ClinVar reference set of biallelic missense single-nucleotide variants with at least two review stars and concordant single-gene ClinVar–VEP annotation, requiring missense consequence on a gene-matched MANE Select or MANE Plus Clinical transcript or, when the gene lacked MANE annotations, a canonical fallback transcript. Variants were annotated with GPN-Star V100 and M447 precomputed scores and evaluated in three nested autosomal cohorts defined by observed rare allele frequency or high allele-number loci with no matching alternate allele in the available gnomAD-derived subset. The expanded cohort (20,089 pathogenic/likely pathogenic and 28,182 benign/likely benign variants) was selected as the primary candidate reference set. A gene-disjoint split assigned 4,794 genes to calibration and 1,198 genes to validation. We estimated class-conditional score likelihood ratios for fixed, exploratory score intervals, with uncertainty estimated using two-sided percentile gene-cluster bootstrap. The selected score intervals were exploratory rather than strictly prespecified before seeing validation data. Because score intervals were refined after exploratory examination and the variant source was a filtered gnomAD subset, evidence mappings were considered provisional pending independent external validation.

**Results 核心表述（暫定）**：

> V100 signed scores discriminated ClinVar P/LP from B/LB missense variants better than M447 across the three study cohorts (expanded-cohort AUC 0.9218 versus 0.8964). Empirical interval likelihood ratios increased monotonically with V100 deleteriousness scores. In gene-disjoint validation, 2<S≤5 yielded an LR of 0.079 (95% gene-bootstrap CI 0.044–0.114), 8<S≤9 yielded an LR of 3.302 (2.597–4.109), and 9<S≤11 yielded an LR of 9.725 (8.097–11.881). Extreme-score intervals showed stronger point estimates but sparse minority-class counts; they were not considered validated Strong evidence categories.

*注意：上述英文段落是內部論文寫作草稿，正式投稿前必須再確認資料來源許可、樣本獨立性、方法說明及 reference citations。*

---

## 15. 尚須確認／尚未完成的上游驗證

- [x] 用 Archive 補齊 Step 1–3、5、7–10 的腳本、指令、輸出與 QC；明確記錄 Step 10 的失敗試跑與 10b 修正。
- [x] 確認 `primary_missense` / `strict_mane_missense` 並非同義，final ≥2 stars 集合有 **89 canonical fallback**。
- [x] 逐一核對 Step 11–21 原始表、counts、同基因資料洩漏、AUC、LR；重新執行 Step 20/21 完全一致。
- [ ] 補 Step 4 `04_unresolved_nonmissense.vcf` 的原始產生命令和用途；Step 6 是否曾獨立執行（若沒有，直接註明 omitted，不要補造）。
- [ ] 保存原始 ClinVar VCF／gnomAD VCF／AN BED／GPN-Star Parquet 的 **URL、MD5/SHA256、版本與下載日期**；本次 ZIP 未含這些大檔，無法做源檔層級再驗證。
- [ ] 使用**官方全量** gnomAD 針對 expanded 中的 putatively-unobserved ALT 獨立核對。不能因 `AN≥1M` 就宣稱 AC=0。
- [x] Step 22 sensitivity 已執行；Step 23b 另修正 ultra-rare 定義，納入通過高 AN/locus audit 的 inferred-unobserved ALT。詳細見第 17–19 節。
- [ ] 檢查 ClinVar classification evidence circularity、gene–disease mechanism、GOF/DN、基因間 label/AF 差異；設計真正獨立 external validation。
- [ ] 最終 threshold 需要 calibration-only rule、額外外部驗證及適當的 clinical lab SOP；目前**僅供研究，不是可直接使用的 PP3/BP4 Strong/Moderate cutoffs**。

## 16. 本版獨立稽核（Archive.zip，2026-10-29）

### 16.1 驗證範圍與結果

1. 解壓並**直接檢查** Archive 裡的 Step 1–21 scripts、logs、results；對照原先 `v0.1` 記錄補正腳本名稱、漏掉的步驟及 final source changes。
2. 以封存 variants/paired/split 表重新計算 **47 項** consistency checks，**47/47 PASS**，包括：ClinVar 226,010 keys 唯一與 star totals、Step 7 final `primary_missense`、Step 14 三組 nested cohorts、Step 15/16 pair coverage、Step 17 全部 12 組 AUC、Step 18/19 gene composition／0 leakage、Step 21 七分數區間的 class counts/LR。**Step 11/12 226,010-row full files 的數字以其原始 summary/log cross-check 為主，沒有從未打包的 source VCF/BED 再運算**。
3. 將 Step 20/21 script 指向 Archive 隔離副本、固定 seed 與 bootstrap 次數重跑：**Step 20 36/36 rows、Step 21 14/14 rows 每個輸出欄位相同（含 CI）**；這證實同版程式的可再現性，不證實 external generalizability。
4. 確認 `MT:12706:T:C` 是唯一 score unsupported variant；所有 **48,271** expanded autosomal variants 有 V100/M447 paired scores。

**本次受檢 Archive.zip SHA256**：`47a7792577d4b7b44e0058983a7e962c3567430fb22b9823a9215257724df080`（供之後確認是否為同一份封存來源）。

### 16.2 關鍵發現及需要更正的舊敘述

| 發現                                                                     | 文件處理                                                      |
| ---------------------------------------------------------------------- | --------------------------------------------------------- |
| Step 1 `variation_id` 最初實際由 `ALLELEID` 暫填                              | 明記 Step 7 改自 VCF ID 取真正 Variation ID，分開存 Allele ID        |
| Step 4 只有 VCF、Step 6 完全沒有獨立檔案                                          | 如實標為**不可完整還原**，不假裝已核實                                     |
| Step 10 初始 script 因缺 `AC_j/AN_j/AF_j` 而失敗                              | 改以 **10b** 從 `AC_e/AN_e/AC_g/AN_g` 建 final AF             |
| Final 使用 `11b_merge_joint_master.py` 而不是不存在的 11 腳本                     | 改列正確腳本                                                    |
| ≥2 stars `primary_missense` 不是全部 MANE-only                             | **89 canonical fallback**，需寫入 Methods                     |
| 59,642 不等於全數 autosomal cohort                                          | 48,271 才是 primary expanded autosomal                      |
| 同一份 gnomAD MANE-filtered VCF 的 exact ALT absence 不能證明 full gnomAD AC=0 | 始終稱 `putatively_unobserved`                               |
| Pejaver local one-sided CI 與本研究 two-sided percentile CI 不同             | Step 21 僅稱 empirical bin LR candidate，不宣稱完整複製原方法          |
| Step 21 Strong tail class counts 極少                                    | `S≤2` validation P/LP=7；`S>11` validation B/LB=16，仍暫緩臨床導入 |

### 16.3 方法學來源／可追溯的外部參考

- **Pejaver et al. (2022)**. *Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria.* *American Journal of Human Genetics* 109:2163–2177. https://doi.org/10.1016/j.ajhg.2022.10.013 ; open access https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/ 。Table 1 可確認 LR threshold（2.406、5.790、33.53、1124 及倒數）。其 threshold-estimation 使用 local neighborhood、one-sided confidence bound 與 bootstrapping，與本研究的固定 bins 有別。
- **gnomAD v4.1 release announcement**（2024-04-19）：https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/ 。官方說明 exome/genome all-sites AN 對 homozygous-reference / absent variation loci 的可用性；AN 是 site-level denominator，不等於特定 ALT 的 observed AC。
- **Song Lab GPN-Star canonical scores documentation**：https://huggingface.co/datasets/songlab/gpn-star-scores 。`llr_calibrated` 為 ALT-vs-REF calibrated LLR（更負較 constrained）；`abs_llr_calibrated` 為另行提供的 calibrated magnitude score，不能取 `abs(llr_calibrated)` 代替。
- **第三方 gnomAD v4.1 MANE-filtered AN tracks**：https://zenodo.org/records/21812932 。特別注意其衍生資料範圍，勿混同於官方 unrestricted all-sites download。

### 16.4 實務建議

本版建議存為 `docs/GPN-Star_Steps01-26_Research_Log_Audited_v0.3_20261030.md`。保留 v0.2 和 Step 19 的 gene manifest、Step 21 七個 interval boundaries 供版本追溯；Step 22–26 維持相同 split 和分箱。

**變更記錄**：

- 2026-10-30 `v0.3` — 保留 v0.2 的逐步 audit，補上 Step 22–26 使用者實際執行輸出、修正 inclusive ultra-rare cohort、gene eligibility sensitivity、C19ORF12 未解決的 4 筆 QC，新增跨 cohort 的候選 PP3/BP4 LR 門檻總表。Step 22–26 尚未由 Archive 原始檔獨立重跑。
- 2026-10-29 `v0.1` — 依先前 terminal 及對話建立初稿。
- 2026-10-29 `v0.2 audited` — 以 Archive 的實際 scripts/logs/results 補齊 Step 1–14、修正 10b/11b 與 MANE/canonical 定義，進行 47 項交叉檢查與 Step 20/21 實際重跑；留下原始未封存資料的追溯缺口和 Step 22 待續。

---

## 17. Step 22 — 原始 AF／cohort sensitivity（已完成）

**來源與稽核級別**：本節數值係研究者在 2026-10-30 提供之 `22_v100_cohort_sensitivity.tsv` 終端輸出；**未**重新取得伺服器原始 TSV 做獨立重跑。腳本 `scripts/22_v100_cohort_sensitivity.py`；使用 Step 19 固定 split、Step 21 七區間，1,000 次 gene-cluster bootstrap；欄位含每 bin 的 P/LP/B/LB、LR、95% CI、涉及基因數與 `sparse`。

四種最初使用的 reference sets：expanded、conservative、observed-only、**observed-only AF<10⁻⁴**。最後一項後來確認為**不完整的 ultra-rare 定義**：忽略高 AN inferred-unobserved ALT，因此僅保留為 historical sensitivity，不作研究主要 ultra-rare 結論（見 Step 23b）。

| Cohort / split                      | P/LP   | B/LB   | Total  |
| ----------------------------------- | ------:| ------:| ------:|
| Expanded / calibration              | 16,079 | 22,548 | 38,627 |
| Expanded / validation               | 4,010  | 5,634  | 9,644  |
| Conservative / calibration          | 14,008 | 22,503 | 36,511 |
| Conservative / validation           | 3,461  | 5,630  | 9,091  |
| Observed / calibration              | 8,608  | 22,442 | 31,050 |
| Observed / validation               | 1,977  | 5,621  | 7,598  |
| Observed-only AF<10⁻⁴ / calibration | 8,200  | 5,103  | 13,303 |
| Observed-only AF<10⁻⁴ / validation  | 1,861  | 1,323  | 3,184  |

**主要發現**：expanded 與 conservative 的 LR 很接近；observed-only 和僅 observed 的 ultra-rare 會改變 score-bin LR，尤其後者較高分 B/LB 的 class-conditional 比例增加。這不能僅憑相關比較分離 AF、gene composition、variant-selection 等原因；Step 23 和 23b 繼續處理。

## 18. Step 23 — Gene universe sensitivity（原版 ultra-rare 定義不完整）

**腳本**：`scripts/23_v100_gene_af_sensitivity.py`；輸出 `results/23_v100_gene_af_sensitivity.tsv`。比較同一 gene split 的三種 gene universes：`all_genes`、`expanded_both_label_genes`、`ultrarare_both_label_genes`；分別做 expanded、observed、observed-only AF<10⁻⁴。gene universe 在 calibration 和 validation **各自計算**，不改 Step 19 split。

- Expanded both-label genes：calibration **1,222**；validation **305**。
- 舊版 observed-ultra-rare both-label genes：calibration **429**；validation **106**。
- Observed-only ultra-rare validation `9<S≤11`：all genes LR **4.820**（95% CI 3.515–6.704）；若限制在這個不完整的 ultra-rare both-label gene universe，LR **6.092**（4.033–10.484）。
- 限制 both-label gene universe **不能等同於已充分控制 gene effects**，因為基因的納入資格受 ClinVar P/LP/B/LB labels 影響；score 分布和變異型別仍可能不同。

**重要修正**：只有 `observed AF<10⁻⁴` 的比較，不包含已通過高 AN/locus audit 的推定未觀察 ALT；不應用來否定 expanded ultra-rare 的證據等級。

## 19. Step 23b — Expanded-inclusive ultra-rare 定義與修正分析

**腳本**：`scripts/23b_v100_ultrarare_inclusive_sensitivity.py`，輸出 `results/23b_v100_ultrarare_inclusive_sensitivity.tsv`，保留 Step 23 的 cohort/gene-universe 架構並改變 ultra-rare membership。

定義（操作性；不把缺失 AF 欄位改寫成真實 AF=0）：

```text
Expanded ultra-rare (autosomal) =
    [expanded & observed_PASS & derived joint AF < 1e-4]
    UNION
    [expanded & high all-sites AN (>=1,000,000)
       & local joint subset does not match the exact ALT
       & locus audit shows no variant record or other ALT only]
```

包含同位置有其他 ALT、但沒有目標 ALT 的 locus audit 類型。第二組必須稱 **inferred-/putatively-unobserved**，因為 local gnomAD source 是第三方 MANE-filtered subset，高 AN 不足以獨立證明官方全量 VCF 的該 ALT AC=0。

| 來源                              | P/LP       | B/LB      | 合計         |
| ------------------------------- | ----------:| ---------:| ----------:|
| Observed joint AF<10⁻⁴          | 10,061     | 6,426     | 16,487     |
| High-AN inferred-unobserved ALT | 9,504      | 119       | 9,623      |
| **Expanded ultra-rare 合計**      | **19,565** | **6,545** | **26,110** |

固定 Step 19 split：calibration **15,671 P/LP + 5,209 B/LB = 20,880**；validation **3,894 P/LP + 1,336 B/LB = 5,230**。新的 both-label gene universe：calibration **575**、validation **154**（與舊版 429/106 不同）。

| Gene universe               | Split       | expanded ultra-rare 變異數 | `2<S≤5` LR (95% CI) | `8<S≤9` LR (95% CI) | `9<S≤11` LR (95% CI) | `S>11` LR (95% CI)      |
| --------------------------- | ----------- | -----------------------:| ------------------- | ------------------- | -------------------- | ----------------------- |
| All genes                   | Calibration | 20,880                  | 0.100 (0.082–0.123) | 2.419 (2.046–2.874) | 5.099 (4.287–6.048)  | 32.706 (21.976–52.173)  |
| All genes                   | Validation  | 5,230                   | 0.081 (0.046–0.120) | 2.252 (1.688–3.005) | 5.702 (4.258–7.872)  | 44.259 (21.227–173.529) |
| Expanded both-label genes   | Validation  | 4,348                   | 0.075 (0.039–0.116) | 2.271 (1.601–3.314) | 5.130 (3.774–7.392)  | 27.088 (12.272–111.481) |
| Ultra-rare both-label genes | Validation  | 3,292                   | 0.053 (0.015–0.101) | 1.987 (1.306–3.063) | 5.429 (3.846–8.000)  | 33.293 (15.452–107.544) |

**主要判讀**：BP4 `2<S≤5` 對 expanded ultra-rare 仍穩定支持 Moderate 候選。PP3 `8<S≤9` 的 validation CI 下界 1.688，**未達 Supporting 的 2.406**；`9<S≤11` validation CI 下界 4.258，**僅能在此 ultra-rare set 保守支持 Supporting，未支持 Moderate 的 5.790**。`S>11` validation B/LB=7，屬稀疏 tail，不賦予 Strong。Gene restrictions 不會完全消除差異。

將 `9<S≤11` 的兩類條件機率拆解：expanded validation P/LP **1,419/4,010=35.39%**，B/LB **205/5,634=3.64%**；expanded ultra-rare validation P/LP **1,396/3,894=35.85%**，B/LB **84/1,336=6.29%**。因此 LR 下降主要反映**該區間 ultra-rare B/LB 佔比上升**，而不是 P/LP 高分比例大幅下降；不能在此斷言造成現象的生物學因果機制。

## 20. Step 24 — High-scoring ultra-rare B/LB source QC（探索性）

**腳本**：`scripts/24_v100_ultrarare_benign_audit.py`；產出 `24_v100_class_conditional_rates.tsv` 與 `24_ultrarare_high_score_BLB.tsv`。定義 `S>8`、B/LB、expanded-inclusive ultra-rare，依原 Step 19 split 分組；不修改 label、不刪變異。

| Split       | S>8 B/LB 數量 | Observed AF<10⁻⁴ | Inferred-unobserved |
| ----------- | -----------:| ----------------:| -------------------:|
| Calibration | 725         | 716              | 9                   |
| Validation  | 187         | 184              | 3                   |
| **Total**   | **912**     | **900**          | **12**              |

Validation 187 筆分布在 `8<S≤9`: **96**、`9<S≤11`: **84**、`S>11`: **7**。其中 S>11 七筆為 `PAX3`（一筆 inferred-unobserved）、`HNF1A`（兩筆）、`FLNC`、`DCTN1`、`RHOBTB2`、`DNM2`。較多 validation high-scoring B/LB 的基因包括 RAI1 (13)、COL5A1 (9)、HNF1A (8)、EHMT1 (8)、TSC2 (8)。

**結論**：高分 ultra-rare B/LB **98.7%（900/912）屬 observed**，並非由 inferred-unobserved records 主導；不應根據 GPN-Star 高分就質疑或剔除其 ClinVar B/LB label。此步僅作 source composition QC，不再進一步逐筆展開 ClinVar mechanism evidence audit。

## 21. Step 25 — High-scoring B/LB 的 ClinVar master 回接（已完成，停止延伸個別審查）

**腳本**：`scripts/25_clinvar_high_score_BLB_review.py`；輸出 `25_high_score_BLB_clinvar_review.tsv`（912）與 `25_priority_BLB_review.tsv`（18）。優先清單是 validation 中 S>11 的 7 筆，聯集所有 high-scoring inferred-unobserved 12 筆，重疊 1 筆，總計 **18**。

原始 master 可回接欄位：`variation_id`、`clnsig`、`stars`、`review_status`、`geneinfo`、MANE flags、`strict_mane_missense`、`gnomad41_absence_review_status` 等。**Master 不包含 condition-specific SCV functional/segregation evidence**；Step 25 只是 metadata merge，不是獨立確認 B/LB classification。研究決策：**停止逐筆 evidence audit；焦點回到 primary expanded cohort 與 pathogenic-gene eligibility**。先前提出但未執行的 ClinVar VCV/SCV 擷取不得記成已完成。

## 22. Step 26 — Gene with ≥1 P/LP / Pathogenic variant（any consequence）eligibility

**目的**：與 primary expanded/all-genes 比較。不是要求本研究 missense subset 內同基因同時有 P/LP/B/LB，而是回到**同一版完整 ClinVar VCF**判斷該 gene 是否至少有一個 P/LP 或 Pathogenic，**不限變異後果**。其中 `Pathogenic-only, any consequence` 是更貼近 Pejaver et al. (2022) described gene inclusion rule 的敏感度定義；`P/LP any consequence` 是較寬版本。

**輸入**：`clinvar_vcf/clinvar_20261025.vcf.gz`（這次研究者已於 server 重新使用原始 VCF，唯未上傳到本對話）+ `results/19_expanded_gene_split.tsv.gz`；腳本 `scripts/26_gene_pathogenic_eligibility.py`；輸出 `26_gene_pathogenic_eligibility.tsv`、`26_gene_restricted_lr.tsv`。掃描完整 ClinVar VCF 4,599,696 records，其中符合腳本 P/LP parsing 且有 gene 的 records **365,790**（此處數字是 records，不是本研究 missense final cohort）。使用同一固定 Step 19 split、相同七個 bins、1,000 次 gene-cluster bootstrap。修復 shell 路徑設定時，`$PWD/clinvar_vcf/clinvar_20261025.vcf.gz` 可用，將 `~` 放在引號內不會自動展開。

| Gene mode（腳本欄位）              | Genes | Calibration P/LP | Calibration B/LB | Validation P/LP | Validation B/LB |
| ---------------------------- | -----:| ----------------:| ----------------:| ---------------:| ---------------:|
| `all_genes`                  | 5,992 | 16,079           | 22,548           | 4,010           | 5,634           |
| `primary_missense_PLP`       | 1,849 | 16,079           | 10,775           | 4,010           | 2,404           |
| `clinvar_any_type_PLP`       | 3,808 | 16,075           | 18,916           | 4,010           | 4,624           |
| `clinvar_any_type_P_only`    | 3,473 | 16,007           | 17,938           | 4,003           | 4,437           |
| `clinvar_any_type_PLP_2star` | 2,766 | 16,075           | 15,449           | 4,010           | 3,999           |

**解讀**：any-consequence P/LP restriction 保留 validation **8,634/9,644 variants**，病性（P/LP）4,010 筆全數保留，但排除一些 B/LB-only genes；P-only 與 ≥2-star gene subsets 也有類似模型 LR 分布。**不要把 P/LP gene universe 偷換成 P-only universe**：兩個定義不同，應在表格保留並列。

### 22.1 Step 26 特定 QC：C19ORF12 四筆 P/LP 未通過 any-type P/LP gene eligibility

研究者從 `results/26_gene_pathogenic_eligibility.tsv` 回查 Step 19，發現 **4 筆** P/LP missense 全部位於 `C19ORF12` 且屬 calibration：

```text
19:29702923:G:A  C19ORF12  calibration
19:29702966:C:T  C19ORF12  calibration
19:29702977:C:A  C19ORF12  calibration
19:29702977:C:T  C19ORF12  calibration
```

同一 ClinVar snapshot 本應至少含這四筆 P/LP，所以 `clinvar_any_type_PLP` 若精確追蹤基因識別，理論上至少應包含此基因。**目前將它列為未解決的 gene-symbol/GENEINFO 或 CLNSIG parser mismatch**（例如基因名稱大小寫或別名不一致），不能在未看原始 VCF 記錄前確定原因，更不能自行修改 label。結果影響 `clinvar_any_type_PLP` calibration P/LP 由 16,079 減為 16,075；`clinvar_any_type_PLP_2star` 也可能受同源 parser 影響。修正後須重跑 Step 26，並保留修前/修後結果；不必重跑 Step 1–25，也不必重切 genes。

### 22.2 Step 26 validation 的核心發現

`2<S≤5` 的 BP4 Moderate 在各 gene restriction 下均支持；`9<S≤11` PP3 Moderate 亦均支持。`8<S≤9` PP3 Supporting 對 confidence-bound choice 敏感：P/LP any-type **CI 下界=2.391**（略低於 Supporting LR=2.406），P-only **2.492** 和 ≥2-star P/LP any-type **2.451** 通過。S>11 雖有高 LR，但 15–16 個 B/LB reference，暫不賦予 Strong。

## 23. 主要研究結果與 PP3/BP4 候選 threshold comparison（截至 Step 26）

**所有資料均使用相同 V100 primary score `S = -llr_calibrated`。分類採 validation 95% gene-cluster bootstrap 的保守端點：BP4 用 CI 上界、PP3 用 CI 下界。** 各 row 的總分母只限於其 cohort/split/gene universe，區間 LR **不可當成簡單的 `S≥threshold` 累積 LR**。沒有顯著資料可確認 Strong 之前，不自動把 S>11 當成 PP3 Strong。

### 23.1 Evidence-strength 參考門檻

| Level      | PP3 LR 最低 | BP4 LR 最高 |
| ---------- | ---------:| ---------:|
| Supporting | 2.406     | 0.4156    |
| Moderate   | 5.790     | 0.1727    |
| Strong     | 33.53     | 0.0298    |

### 23.2 **Primary expanded / all genes**：全部七區間

| V100 score interval | Calibration LR (95% CI)   | Validation LR (95% CI)     | Interpretation                                   |
| ------------------- | ------------------------- | -------------------------- | ------------------------------------------------ |
| S≤2                 | 0.0136 (0.0092–0.0189)    | 0.0064 (0.0020–0.0126)     | BP4 Strong 數值候選，validation P/LP=7，sparse，**不採用** |
| 2<S≤5               | 0.0902 (0.0734–0.1100)    | 0.0787 (0.0436–0.1144)     | **BP4 Moderate candidate**                       |
| 5<S≤6               | 0.4548 (0.3853–0.5339)    | 0.3412 (0.2330–0.4501)     | No assignment                                    |
| 6<S≤8               | 1.3549 (1.1982–1.5251)    | 1.1631 (0.8773–1.4636)     | No assignment                                    |
| 8<S≤9               | 3.6816 (3.2448–4.1968)    | 3.3017 (2.5968–4.1086)     | **PP3 Supporting candidate**                     |
| 9<S≤11              | 9.2352 (8.2661–10.3567)   | 9.7253 (8.0974–11.8814)    | **PP3 Moderate candidate**                       |
| S>11                | 57.8681 (42.5595–81.2121) | 79.6452 (44.1589–177.2135) | PP3 Strong 數值候選，validation B/LB=16，sparse        |

*註：上表來自 Step 21 **2,000 次** gene-bootstrap；下兩表 Step 22/23b/26 使用 **1,000 次**，因此 primary 的 CI 端點在不同表格稍有變動，非統計錯誤。*

### 23.3 Gene eligibility sensitivity：固定 expanded cohort 的 validation LR

| Gene universe                      | Genes | Validation P/LP / B/LB | BP4: 2<S≤5 (LR; CI; 判定)              | PP3: 8<S≤9 (LR; CI; 判定)                  | PP3: 9<S≤11 (LR; CI; 判定)              | S>11 LR（sparse）    |
| ---------------------------------- | -----:| ----------------------:| ------------------------------------ | ---------------------------------------- | ------------------------------------- | ------------------ |
| All genes (primary)                | 5,992 | 4,010/5,634            | 0.079 (0.043–0.117); BP4 Moderate 候選 | 3.302 (2.585–4.060); PP3 Supporting 候選   | 9.725 (8.138–11.819); PP3 Moderate 候選 | 79.645（B/LB=16，稀疏） |
| ≥1 P/LP missense in primary        | 1,849 | 4,010/2,404            | 0.078 (0.044–0.114); BP4 Moderate 候選 | 3.131 (2.329–4.268); 不達 Supporting 的保守界限 | 7.595 (6.068–10.051); PP3 Moderate 候選 | 36.250（B/LB=15，稀疏） |
| ≥1 P/LP, any consequence           | 3,808 | 4,010/4,624            | 0.078 (0.043–0.114); BP4 Moderate 候選 | 3.072 (2.391–3.797); 不達 Supporting 的保守界限 | 8.658 (7.036–10.814); PP3 Moderate 候選 | 69.725（B/LB=15，稀疏） |
| ≥1 Pathogenic, any consequence     | 3,473 | 4,003/4,437            | 0.078 (0.042–0.112); BP4 Moderate 候選 | 3.153 (2.492–4.015); PP3 Supporting 候選   | 8.701 (7.144–10.733); PP3 Moderate 候選 | 67.022（B/LB=15，稀疏） |
| ≥1 P/LP, ≥2 stars, any consequence | 2,766 | 4,010/3,999            | 0.078 (0.043–0.114); BP4 Moderate 候選 | 3.155 (2.451–3.973); PP3 Supporting 候選   | 8.629 (6.942–10.732); PP3 Moderate 候選 | 60.301（B/LB=15，稀疏） |

### 23.4 AF/cohort sensitivity：固定 all-genes gene-disjoint validation LR

| Cohort                                                              | BP4: 2<S≤5 (LR; CI; 判定)                | PP3: 8<S≤9 (LR; CI; 判定)                  | PP3: 9<S≤11 (LR; CI; 判定)               | S>11 LR（sparse）                         |
| ------------------------------------------------------------------- | -------------------------------------- | ---------------------------------------- | -------------------------------------- | --------------------------------------- |
| Expanded primary / all genes                                        | 0.079 (0.044–0.114); BP4 Moderate 候選   | 3.302 (2.597–4.109); PP3 Supporting 候選   | 9.725 (8.097–11.881); PP3 Moderate 候選  | 79.645 (44.159–177.213); **B/LB=16，稀疏** |
| Conservative                                                        | 0.084 (0.047–0.127); BP4 Moderate 候選   | 3.381 (2.706–4.215); PP3 Supporting 候選   | 9.570 (7.922–11.818); PP3 Moderate 候選  | 75.540 (40.895–163.804); **B/LB=16，稀疏** |
| Observed-only (AF<0.01)                                             | 0.146 (0.094–0.193); BP4 Supporting 候選 | 4.265 (3.578–5.144); PP3 Supporting 候選   | 8.109 (6.731–10.096); PP3 Moderate 候選  | 37.530 (20.886–79.047); **B/LB=15，稀疏**  |
| Expanded ultra-rare: observed AF<10⁻⁴ + high-AN inferred-unobserved | 0.081 (0.046–0.120); BP4 Moderate 候選   | 2.252 (1.688–3.005); 不達 Supporting 的保守界限 | 5.702 (4.258–7.872); PP3 Supporting 候選 | 44.259 (21.227–173.529); **B/LB=7，稀疏**  |
| Observed ultra-rare ONLY: AF<10⁻⁴                                   | 0.153 (0.100–0.202); BP4 Supporting 候選 | 2.903 (2.212–3.734); 不達 Supporting 的保守界限 | 4.820 (3.515–6.704); PP3 Supporting 候選 | 22.986 (9.972–88.537); **B/LB=6，稀疏**    |

**比較重點**：

1. **Primary expanded**：`2<S≤5` BP4 Moderate；`8<S≤9` PP3 Supporting；`9<S≤11` PP3 Moderate。Strong tails 暫緩。
2. **Gene with ≥1 P/LP, any consequence**：BP4 Moderate 與 PP3 Moderate 保持；PP3 Supporting 的 CI 下界 **2.391** 只比門檻 **2.406** 低 **0.015**，因此是 borderline，不宜四捨五入後宣稱完全通過。**Pathogenic-only any type** 的 PP3 Supporting lower CI **2.492** 則通過。C19ORF12 4 筆 gene parser QC 尚待修正。
3. **Conservative**：三個核心 evidence tiers 與 primary 一致，說明 same-locus other-ALT candidate 的增減沒有主導結論。
4. **Observed-only AF<0.01**：BP4 2–5 的 CI 上界約 **0.193**，只支持 Supporting；兩個 PP3 區間仍分別支持 Supporting、Moderate。
5. **Expanded inclusive ultra-rare**：BP4 Moderate 維持；`8<S≤9` 未穩健達 PP3 Supporting；`9<S≤11` 僅穩健支持 Supporting 而**未**達 Moderate。這是與 primary calibration 最重要的差別，不可忽略。
6. **舊版 observed-only ultra-rare**：不包括 high-AN inferred-unobserved variants，是比較資料來源的 historical sensitivity，**不是本研究採納的最終 ultra-rare 定義**。
7. **Gene / AF 表不是完全因子設計**：Step 26 的 pathogenic-gene restriction 是 expanded cohort 下的比較，Step 23b 的 ultra-rare gene universes 是另一套 comparisons；不能任意把兩表的 LR 組合成「P/LP any-type gene ∩ expanded ultra-rare」的結果。若未來需要真正交叉分析，必須執行相應的 subset LR 計算。

**不宜直接宣稱數值 cutoff `S≥9` 等於 PP3 Moderate**：目前是 `9<S≤11` *bin-specific LR*；對 `S>11` 的 clinical Strong 尚未確認，是否至少可用 Moderate 也需要以合併高分區間或連續 local LR 另行估計，不能偷換分箱條件。

## 24. 方法學狀態、未解決項目與接續方向

### 24.1 已確定的研究決策

- **Primary cohort**：expanded, 48,271 autosomal variants，ClinVar ≥2 stars、gene-matched MANE/canonical missense + observed AF<0.01 或 high-AN inferred-unobserved ALT。
- **Primary predictor**：V100 signed `S=-llr_calibrated`；M447 和 absolute scores secondary。
- **Primary split**：Step 19 原 gene-disjoint 80/20，不重新切分。
- **Sensitivity**：conservative、observed-only、expanded inclusive ultra-rare，和 any-consequence P/LP/P-only pathogenic-gene eligibility；both-label-gene analysis 只作探索性。
- **Step 24/25**：只保留高分 B/LB composition QC，停止個別 ClinVar case review；不因分數高剔除 B/LB。

### 24.2 尚待完成／不可誤稱已驗證

- [ ] **C19ORF12 四筆 P/LP** 的 VCF `CLNSIG`/`GENEINFO`、gene symbol alias/case normalization mismatch：確定原因並重跑 Step 26 gene eligibility；目前只確認輸出中的不一致，**沒有原始 VCF row-level 查核結果**。
- [ ] 官方完整 gnomAD variant callset 上對 inferred-unobserved ALT absence 的 cross-check；all-sites AN 本身不證明 AC=0。
- [ ] ClinVar reference classification 與 PP3/BP4 其他計算證據的可能 circularity；不同 gene–disease mechanism 的適用性。
- [ ] Step 21 fixed-bin mapping 來自查看 calibration 和 validation 描述性結果後的探索性分箱，**現有 validation 並非完全 untouched**；最終臨床可用門檻需真正 external validation。
- [ ] 若採 Pejaver-style continuous local LR：須在 calibration-only 設定平滑、鄰域、one-sided conservative-bound 等規則；保留 held-out genes 不用於調參，另行外部驗證。不可把目前的 two-sided percentile bin CI 寫成完整 Pejaver procedure。
- [ ] BP4/PP3 Strong tails 的少數反向類別樣本太少，即使 point estimate/CI 看似符合 Strong，也暫不當作臨床 cutoff。

### 24.3 主要檔案與分析版本追蹤

```text
scripts/22_v100_cohort_sensitivity.py
results/22_v100_cohort_sensitivity.tsv
scripts/23_v100_gene_af_sensitivity.py
results/23_v100_gene_af_sensitivity.tsv
scripts/23b_v100_ultrarare_inclusive_sensitivity.py
results/23b_v100_ultrarare_inclusive_sensitivity.tsv
scripts/24_v100_ultrarare_benign_audit.py
results/24_v100_class_conditional_rates.tsv
results/24_ultrarare_high_score_BLB.tsv
scripts/25_clinvar_high_score_BLB_review.py
results/25_high_score_BLB_clinvar_review.tsv
results/25_priority_BLB_review.tsv
scripts/26_gene_pathogenic_eligibility.py
results/26_gene_pathogenic_eligibility.tsv
results/26_gene_restricted_lr.tsv
```

**本 v0.3 稽核聲明**：Step 1–21 的 47/47、36/36、14/14 重跑與 Archive audit 是 **v0.2 的歷史結果**；Step 22–26 的數值來自研究者直接貼出的終端執行結果，這次根據跨表 count/split/threshold 重新做了有限的**算術一致性核對**，**沒有**直接存取使用者伺服器、原始 VCF、gnomAD 或 Step 22–26 原始 TSV，更沒有獨立重跑 Step 22–26。對未呈現的 score bins、gene∩AF 交叉子集、尚未解決的 C19ORF12 成因，一律不補造數字。

---

**更新日期**：2026-10-30。**版本**：v0.3；完整保存 v0.2 稽核版正文並將 Step 22–26 補附於第 17–24 節。
