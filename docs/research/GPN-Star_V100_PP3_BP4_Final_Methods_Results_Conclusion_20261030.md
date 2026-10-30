# GPN-Star V100 於 ClinVar missense SNVs 的 PP3/BP4 候選證據校正：最終採用流程

**文件類型：Methods → Results → Conclusion（研究成果整理稿）**  
**資料截至：2026-10-30｜分析進度：Step 27A–E｜版本：v1.0（research-only）**  
**研究專案：** `/home/n102968/GPN-star_calibration_clinvar20261025`  
**主要研究目標：** 以 ClinVar P/LP、B/LB 參考變異評估 GPN-Star V100 的辨識能力，並對連續 signed score 執行 local likelihood-ratio calibration，估計 PP3/BP4 候選證據等級的 score intervals。

> **本文只呈現最後採納的分析鏈**：ClinVar 基因對應／consequence QC、gnomAD v4.1 joint-derived AF + all-sites AN/locus audit、expanded primary cohort、V100/M447 paired discrimination、固定 gene-disjoint split、**adaptive continuous local LR**、fixed-interval validation 及 AF/gene-eligibility/bandwidth sensitivity。其他早期探索流程與研究歷史另存於完整 v0.4 Research Log；本文件不將那些探索性結果列入 final Methods。

---

## Methods

### 1. Study design and reference classification

本研究為 ClinVar-based retrospective computational evidence calibration。主要參考來源為 **2026-10-25 ClinVar GRCh38 snapshot**（原始 VCF `clinvar_vcf/clinvar_20261025.vcf.gz`）；初篩基因變異分類為 **Pathogenic/Likely pathogenic（P/LP）** 與 **Benign/Likely benign（B/LB）**，其他 conflicting/VUS 類別不納入最終 reference set，並要求 **ClinVar review status ≥2 stars**。研究僅納入 **autosomal single-nucleotide missense variants**。不將 ClinVar 星等視為對 gene–disease mechanism 或分類獨立性的完整保證。

### 2. Missense consequence and gene annotation

使用 **Ensembl VEP v115**（GRCh38）進行 **all-transcript consequence annotation**，以 ClinVar 單一 gene designation 與 VEP `SYMBOL` 對應。入選 missense variant 須位於該 gene 的 protein-coding **MANE Select** 或 **MANE Plus Clinical** transcript；若該 gene **無可用 MANE transcript**，才允許 canonical transcript missense 作 fallback。最終在 review-status ≥2-stars 的 **59,642** 筆 primary-missense SNVs 中，**59,553** 由 MANE 對應判定，**89** 為 canonical fallback；P/LP **22,589**、B/LB **37,053**。

**VEP 在最後採納流程中的角色是 gene/transcript/consequence 判定，而不是作為最終 gnomAD AF annotation source。**

### 3. Allele-frequency annotation and source audit

**最終 AF source** 使用本地 **gnomAD v4.1 joint-derived MANE-filtered VCF**，依同一 GRCh38 `chrom:pos:ref:alt` 精確比對目標 alternate allele；在對應 INFO 組件有效時，計算：

\[
\mathrm{AF}_{joint}=\frac{AC_e+AC_g}{AN_e+AN_g}.
\]

原始 exome/genome 對應 ALT AC/AN 組件不完整者不得自行補零；只將符合來源 `FILTER=PASS` 且有有效 joint AF 的 variants 當作 observed high-confidence calls。後續使用 **gnomAD v4.1 exome/genome all-sites AN** 的 site-level tracks，判斷在 joint-derived subset 中找不到目標 ALT 的 loci 是否具有充分 callability，並再次於同一 variant subset VCF 執行 locus-level 查詢，區分：

- **High-AN no-variant-record locus**：在來源 subset 內同位置無 variant record；
- **High-AN other-ALT locus**：同位置有其他 ALT，但未找到本研究的目標 ALT。

設定 summed joint all-sites AN **≥1,000,000**。此程序可將 unmatched ALT 判定為 **putatively unobserved in the available source**，**不能證明在官方完整 gnomAD 全量 VCF 中目標 ALT 的 AC=0**。All-sites AN 是位點層級 denominator，不能直接取代 ALT-specific AC。

### 4. Primary autosomal reference cohort

建立三個巢狀 cohorts：

| Cohort | 定義 | P/LP | B/LB | N |
|---|---|---:|---:|---:|
| Observed-only | Joint-derived VCF 中 PASS 且 **AF<0.01** | 10,585 | 28,063 | **38,648** |
| Conservative | Observed-only + high-AN、同位置無 variant record、目標 ALT inferred-unobserved | 17,469 | 28,133 | **45,602** |
| **Expanded（primary）** | Conservative + high-AN、同位置有別的 ALT 但無目標 ALT | **20,089** | **28,182** | **48,271** |

選擇 **expanded** 作 primary，目的為避免單純因 gnomAD-derived VCF 未收錄 extremely rare alternate alleles，系統性遺失大量 ClinVar P/LP variants。此選擇係在研究探索過程中形成，並非 pre-registered，故設計 nested cohort 與 AF sensitivity checks。所有 expanded autosomal variants 都具有 GPN-Star 兩版本的 score。

### 5. GPN-Star score annotation and predictive discrimination

從本地 GRCh38 GPN-Star **V100** 與 **M447** precomputed per-chromosome Parquet 以 genomic allele key 配對。主要 predictor 選用 V100 **signed score**：

\[
S=-\mathrm{llr\_calibrated},\quad S\uparrow\text{ indicates increasing predicted deleteriousness}.
\]

M447 signed 作為主要 discrimination comparator；兩模型的其他 score 欄位為 secondary descriptors，不用於本次最終 PP3/BP4 threshold selection。比較 ROC AUC、average precision（AP）及 paired AUC difference。Step 17 的 AUC 95% CI 使用 **500 次 variant-level stratified paired bootstrap**，須與後續 LR 的 **gene-cluster bootstrap** 區分。

### 6. Gene-disjoint calibration and internal validation

依 `calibration_gene` 固定 **gene-disjoint calibration/validation split**；使用原 Step 19 gene manifest、不在後續依 LR 重新切分；**0 genes overlap**。

| Split | Genes | P/LP | B/LB | Variants |
|---|---:|---:|---:|---:|
| Calibration | **4,794** | **16,079** | **22,548** | **38,627** |
| Validation | **1,198** | **4,010** | **5,634** | **9,644** |
| Total | **5,992** | **20,089** | **28,182** | **48,271** |

因研究者先前曾檢視 validation score profile 和早期區間結果，這裡嚴格稱為 **exploratory gene-disjoint internal validation**，不稱完全 untouched held-out、獨立外部驗證。

### 7. Continuous local LR calibration（primary threshold-development method）

**Step 27A calibration only**，對連續 V100 `S` 的 **0.1-score grid**（實際 −7.2 至 13.4，共 207 points）估計局部 class-conditional likelihood ratio：

\[
\widehat{LR}_{local}(s)=
\frac{n_{P}(s-\epsilon(s)\leq S\leq s+\epsilon(s))/N_{P}}
{n_{B}(s-\epsilon(s)\leq S\leq s+\epsilon(s))/N_{B}}.
\]

Adaptive symmetric windows 的半寬 `ε(s)` 為使局部至少包含 **300 variants**、**25 P/LP** 與 **25 B/LB** 所需的最小距離；事後 reliability criteria 包括 **每類 ≥10 genes** 和 **maximum halfwidth 1.5**。使用 **1,000 gene-cluster bootstrap replicates**（seed 20261030）取得**逐點單側 95%** lower/upper bounds，僅在主要最長 supported continuous domain（**−0.7–12.9**）內依保守端點搜尋 evidence thresholds。BP4 要求 local **upper bound ≤** 相應 LR 目標，PP3 要求 local **lower bound ≥** 相應 LR 目標，並要求在目前 **supported domain** 中向相應尾端的連續 grid points 維持條件成立。

採用 evidence LR 參考目標：

| Strength | PP3 LR lower threshold | BP4 LR upper threshold |
|---|---:|---:|
| Supporting | **2.406** | **0.4156** |
| Moderate | **5.790** | **0.1727** |
| Strong | **33.53** | **0.0298** |

這是一個 **Pejaver-inspired adaptive local LR** 實作，**不是** Pejaver et al. 完整的 reference-background/prior calibration。此 bootstrap **固定**由原 calibration variants 選出的 windows，未量化 window-selection uncertainty；逐點界限不是 entire-curve simultaneous confidence band。Grid step=0.1 **不是** cutoff uncertainty=±0.1。

### 8. Frozen interval LR and evidence validation

使用 Step 27A 的 original bandwidth 設定，在 calibration 中固定以下六個 boundaries：**1.7、3.8、5.1、8.4、9.8、11.6**，並劃分七個互斥 score intervals。依 Step 27B 在 validation set 直接計算每個互斥 interval 的：

\[
\widehat{LR}(I)=\frac{n_{P,I}/N_P}{n_{B,I}/N_B}.
\]

報告 **2,000 次 gene-cluster bootstrap** 的 two-sided 95% CI 與 one-sided 95% bounds。BP4 以 **one-sided upper** 比較 upper LR 目標；PP3 以 **one-sided lower** 比較 lower LR 目標。另設研究可靠性規則：每區間 P/LP、B/LB 各 **≥20 variants** 且各 **≥10 genes**，不足標記 **sparse，不賦予 Strong**。Cumulative `S≤t` 或 `S≥t` tail LR 僅作 secondary QC，**不等同**互斥 evidence interval LR。

### 9. Sensitivity analyses（均固定原本 score cutoffs）

**AF/cohort sensitivity（Step 27C）**：expanded（primary）、conservative、observed-only（AF<0.01）、expanded ultra-rare（**observed AF<10⁻⁴ + high-AN inferred-unobserved ALT**）。

**Gene-eligibility sensitivity（Step 26/27C）**：在相同的 **2026-10-25 完整 ClinVar VCF** 中判斷 gene 是否有 **≥1 P/LP（any consequence）**、**≥1 Pathogenic（any consequence）**、或 **≥1 P/LP with ≥2 review stars（any consequence）**。此選取與「本研究中必須同時有 P/LP、B/LB missense」不同。原 Step 26 對 `C19ORF12` 有四筆 P/LP calibration variants 未通過 any-type P/LP gene membership，列為尚待原 VCF `GENEINFO`/symbol/CLNSIG parser QC 的限制。

**Bandwidth sensitivity（Step 27D–E）**：只在 calibration data 比較 smaller（minimum pooled **200**、each class **20**）、original（**300/25**）和 larger（**500/40**）adaptive windows，其他配置相同；將各自得出的六個 cutoffs **不再微調**地套回相同 gene-disjoint validation set。Step 27C/E 的 LR 使用 **2,000 gene-cluster bootstrap replicates**；以上 cohort/gene/bandwidth comparisons 不是彼此獨立的 external validations。

---

## Results

### 1. Reference set and model discrimination

原始 ClinVar VCF **4,599,696** records；經 final consequence、review status 選樣所得 ≥2-stars primary missense SNVs **59,642**，其中 expanded autosomal primary **48,271**（P/LP 20,089，B/LB 28,182），涵蓋 **5,992 genes**。本地 GPN-Star V100 與 M447 對 expanded primary variants 均有分數（coverage **100%**）。

| Cohort | V100 signed AUC（95% variant bootstrap CI） | M447 signed AUC | Paired ΔAUC V100–M447（95% CI） |
|---|---|---:|---|
| Observed-only | **0.8898** (0.8862–0.8931) | 0.8652 | +0.0247 (0.0225–0.0270) |
| Conservative | **0.9184** (0.9160–0.9209) | 0.8925 | +0.0259 (0.0239–0.0278) |
| **Expanded primary** | **0.9218** (0.9194–0.9240) | **0.8964** | **+0.0254** (0.0236–0.0273) |

V100 signed 在三個 cohorts 均優於 M447 signed。但不同 cohort 的 class composition 與 frequency distribution 改變，AUC/AP 不可解讀成純粹模型本身的因果改善。

### 2. Calibration-derived continuous thresholds

Step 27A original configuration 在 calibration set 找到六個候選 score cutoffs。這裡顯示每個 **cutoff local score neighborhood** 的 LR，非後續 evidence interval LR。

| Evidence cutoff | Score criterion | Local LR | One-sided conservative bound | Reference target |
|---|---|---:|---:|---:|
| BP4 Strong | `S≤1.7` | 0.0201 | **upper 0.0290** | ≤0.0298 |
| BP4 Moderate | `S≤3.8` | 0.0831 | **upper 0.1178** | ≤0.1727 |
| BP4 Supporting | `S≤5.1` | 0.3029 | **upper 0.3915** | ≤0.4156 |
| PP3 Supporting | `S≥8.4` | 4.8875 | **lower 3.9560** | ≥2.406 |
| PP3 Moderate | `S≥9.8` | 7.7433 | **lower 6.0130** | ≥5.790 |
| PP3 Strong | `S≥11.6` | 61.5341 | **lower 43.7351** | ≥33.53 |

Calibration score grid 207 points（−7.2–13.4）；最長 supported local domain **−0.7 至 12.9**。在 score tails，adaptive halfwidth 可能顯著擴張（如 S=12.0 時 ε=0.596），所以 0.1 grid 不是生物學／統計 cutoff 的實際精確度。

### 3. Fixed-interval gene-disjoint validation（primary evidence-mapping result）

六個 cutoffs 固定後，於 gene-disjoint validation（P/LP 4,010；B/LB 5,634）得到：

| Mutually exclusive V100 score interval | Candidate interpretation | Validation P/LP / B/LB | Empirical interval LR（two-sided 95% CI） | One-sided bound / target | Decision |
|---|---|---:|---|---|---|
| `S ≤1.7` | BP4 Strong | **6 / 1,358** | 0.006 (0.002–0.012) | upper 0.011 / ≤0.0298 | **Sparse；not confirmed** |
| **`1.7 < S ≤3.8`** | **BP4 Moderate** | **37 / 1,400** | **0.037 (0.021–0.053)** | **upper 0.051 ≤0.1727** | **Internally supported** |
| **`3.8 < S ≤5.1`** | **BP4 Supporting** | **88 / 907** | **0.136 (0.069–0.211)** | **upper 0.199 ≤0.4156** | **Internally supported** |
| `5.1 < S <8.4` | No evidence | 1,126 / 1,605 | 0.986 (0.752–1.204) | N/A | No assignment |
| **`8.4 ≤ S <9.8`** | **PP3 Supporting** | **1,020 / 270** | **5.308 (4.416–6.338)** | **lower 4.561 ≥2.406** | **Internally supported** |
| **`9.8 ≤ S <11.6`** | **PP3 Moderate** | **1,173 / 87** | **18.943 (14.473–25.657)** | **lower 15.164 ≥5.790** | **Internally supported** |
| `S ≥11.6` | PP3 Strong | **560 / 7** | 112.399 (48.378–378.237) | lower 56.614 / ≥33.53 | **Sparse；not confirmed** |

**QC**：Independent recount from Step 19 與 Step 27B TSV 顯示 `S≤1.7` B/LB **1,358 across 545 genes**，`8.4≤S<9.8` B/LB **270 across 167 genes**；Step 27B label/gene counts 可追溯。PP3 Supporting interval LR point estimate 5.308 **不表示**達 Moderate：其 one-sided lower bound **4.561<5.790**，故僅支持 Supporting。PP3 Moderate 也未達 Strong 的 33.53 lower-bound requirement。

### 4. AF and gene-eligibility sensitivity（fixed original thresholds）

Step 27C 使用相同六個 score boundaries，對七個**重疊** reference-set definitions 分別重算 interval LR；下表每格為 **LR [one-sided conservative bound]**（BP4 upper、PP3 lower）。四個候選 tiers 在七組均通過相對應 evidence LR target 和每類樣本／基因數 QC。

| Cohort / gene restriction | BP4 Moderate | BP4 Supporting | PP3 Supporting | PP3 Moderate |
|---|---|---|---|---|
| **Expanded primary** | **0.037 [0.051]** | **0.136 [0.199]** | **5.308 [4.561]** | **18.943 [15.164]** |
| Conservative | 0.041 [0.055] | 0.140 [0.203] | 5.374 [4.627] | 18.286 [14.664] |
| Observed-only AF<0.01 | 0.071 [0.092] | 0.238 [0.326] | 5.920 [5.155] | 12.843 [10.166] |
| **Expanded ultra-rare** | 0.037 [0.052] | 0.140 [0.202] | **3.509 [2.837]** | **11.036 [7.859]** |
| ≥1 ClinVar P/LP, any consequence | 0.038 [0.051] | 0.133 [0.195] | 4.901 [4.175] | 17.566 [13.840] |
| ≥1 ClinVar Pathogenic, any consequence | 0.038 [0.052] | 0.133 [0.193] | 4.930 [4.181] | 17.321 [13.498] |
| ≥1 ClinVar P/LP ≥2 stars, any consequence | 0.038 [0.052] | 0.131 [0.191] | 4.962 [4.149] | 16.711 [13.184] |

特別是 clinically relevant 的 **expanded ultra-rare**（observed AF<10⁻⁴ + high-AN inferred-unobserved）保留了 **26,110** variants（calibration 20,880；validation **5,230**，其中 P/LP **3,894**、B/LB **1,336**）：PP3 Supporting lower **2.837≥2.406**、PP3 Moderate lower **7.859≥5.790**。**原始高 AN candidate source 的 ALT absence 尚未以官方全量 gnomAD 獨立交叉驗證**。

### 5. Adaptive-window robustness（original thresholds preserved）

Step 27D 在 calibration-only 比較三組 adaptive windows；Step 27E 不經重新微調地將相應 cutoffs 固定用於既有 validation genes。

| Evidence cutoff | Smaller window | **Original window** | Larger window |
|---|---:|---:|---:|
| BP4 Strong | 1.7 | **1.7** | 1.6 |
| BP4 Moderate | 3.8 | **3.8** | 4.2 |
| BP4 Supporting | 5.1 | **5.1** | 5.1 |
| PP3 Supporting | 8.4 | **8.4** | 8.3 |
| PP3 Moderate | 9.9 | **9.8** | 9.5 |
| PP3 Strong | 11.6 | **11.6** | 11.5 |

四個非 Strong intervals 在三組 settings 均維持對應 evidence tier：

| Evidence | Smaller LR [bound] | **Original LR [bound]** | Larger LR [bound] | Evidence target |
|---|---|---|---|---|
| BP4 Moderate | 0.037 [upper 0.051] | **0.037 [0.051]** | 0.042 [0.055] | ≤0.1727 |
| BP4 Supporting | 0.136 [upper 0.199] | **0.136 [0.199]** | 0.168 [0.250] | ≤0.4156 |
| PP3 Supporting | 5.449 [lower 4.701] | **5.308 [4.561]** | 4.553 [3.852] | ≥2.406 |
| PP3 Moderate | 20.199 [lower 15.875] | **18.943 [15.164]** | 14.983 [12.421] | ≥5.790 |

Bandwidth selection 使 BP4 Moderate cutoff 變動約 0.4 分、PP3 Moderate cutoff 變動約 0.4 分，顯示 cutoff boundary location 有估計不確定性，**不是**所有參數都得到完全相同的 score threshold。原 original thresholds 已保留，不依 validation 選擇較有利設定。兩端 Strong 的少數 label 始終 sparse（BP4 Strong validation P/LP 6；PP3 Strong validation B/LB 7）。

---

## Conclusion

### Main findings

在 ClinVar 2026-10-25、GRCh38、≥2-star、gene-matched missense SNV reference set 中，V100 signed score 於 expanded primary cohort 的 discriminatory performance **優於 M447**（ROC AUC **0.9218** vs **0.8964**）。以 expanded gene-disjoint calibration set 建立的 **adaptive continuous local LR**，在 **0.1-score grid** 下得到 **1.7/3.8/5.1/8.4/9.8/11.6** 六個候選界線；將這些界線固定於原 gene-disjoint validation 後，四個**互斥** evidence intervals 具有本研究預先設定的樣本量／基因數支持，以及達相應門檻的單側 bootstrap confidence bounds。

**主要研究候選 mapping（V100 `S=-llr_calibrated`）：**

| V100 score | Candidate evidence level | Status |
|---|---|---|
| **`1.7 < S ≤3.8`** | **BP4 Moderate** | Internally supported |
| **`3.8 < S ≤5.1`** | **BP4 Supporting** | Internally supported |
| **`8.4 ≤ S <9.8`** | **PP3 Supporting** | Internally supported |
| **`9.8 ≤ S <11.6`** | **PP3 Moderate** | Internally supported |
| `S≤1.7` / `S≥11.6` | BP4 Strong / PP3 Strong | **Not confirmed: sparse opposite-label reference variants** |
| `5.1<S<8.4` | No evidence | No assignment |

**Robustness**：四個 candidate tiers 在七種互相重疊的 AF/gene-eligibility reference definitions、三種 adaptive bandwidth settings 下均維持同一證據等級；其中 ultra-rare cohort 的 PP3 LR 下降，但單側保守界限仍達相應 Supporting/Moderate 門檻。

### Limitations and clinical implications

此 mapping **屬探索性、研究用途的 candidate calibration，而非可直接作為 ACMG/AMP clinical variant classification 的核准 cutoff**。研究採 ClinVar labels，存在 computational-evidence circularity；variants 可能不獨立，已以 gene-cluster bootstrap 部分處理，但 disease mechanism / gene-specific validity 未確效。過往研究已檢視 validation score 相關資訊，故目前 validation 不是真正 untouched external dataset；score threshold selection 也只採 **pointwise one-sided bootstrap**，不是 simultaneous band，且 adaptive-window-selection uncertainty 未完整納入。Expanded cohort 的高 AN inferred-unobserved ALT 只對來源 subset 有效，尚未於官方全量 gnomAD 驗證；`C19ORF12` gene-eligibility parser QC 未解決。Strong evidence tails 的 minority-class reference 太少，**目前不應升級 Strong**。

### Recommended next stage

凍結 original Step 19 gene manifest、Step 27A original threshold table、Step 27B–E outputs、對應 scripts/logs 與資料來源 checksum；下一步優先進行**真正獨立的 clinical / external reference validation**，並以 gene–disease mechanism、evidence provenance 和 full-source gnomAD ALT cross-check 評估臨床適用性。在外部結果取得前，不因重複使用同一批 validation records 而繼續微調 0.1-score cutoffs。

---

## Data provenance / reproducibility notes

- **ClinVar**：2026-10-25 GRCh38 source；`clinvar_vcf/clinvar_20261025.vcf.gz`；gene/transcript consequence 由 VEP v115 完成。來源是既有分析的本地 snapshot，不代表本次文件製作重新下載或驗證 source checksum。
- **Retained scripts/results**：`07_make_clinvar_master.py`、`10b_annotate_gnomad_joint.py`、`11b_merge_joint_master.py`、`12_annotate_all_sites_AN.py`、`13_audit_unobserved_loci.py`、`14_build_calibration_cohorts.py`、`15_annotate_gpnstar_local.py`、`16_merge_gpnstar_models.py`、`17_compare_gpnstar_auc.py`、`18_gene_composition_qc.py`、`19_gene_disjoint_split.py`、`26_gene_pathogenic_eligibility.py`、`27a_v100_adaptive_local_lr.py`、`27b_v100_fixed_interval_validation.py`、`27c_v100_fixed_cutoff_sensitivity.py`、`27d_bandwidth_sensitivity.py`、`27e_v100_{smaller,larger}_validation.py`，及其 `results/`、`logs/` 同名輸出。
- **稽核區分**：Steps 1–21 曾對原 Archive 執行 47 項一致性檢查（詳細見完整 v0.4 history）。Steps 22–27E 的新增數字以研究者直接提供的 server outputs 為基礎，本文製作時做了 cross-table count/LR/boundary consistency checks，**沒有存取研究者 live server 重新計算 Step 27A–E**。
- **主要方法學參考**：Pejaver V et al. *Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria*. *Am J Hum Genet*. 2022;109:2163–2177. DOI: https://doi.org/10.1016/j.ajhg.2022.10.013. 本研究為 *Pejaver-inspired* local LR adaptation，不能宣稱完整複製該文 calibration pipeline。

**本文件與完整研究歷程的區別**：本文件僅保留最終採用 workflow（含必要的 sensitivity/QC），適合後續延伸成論文 Methods/Results 或工作報告；完整舊方法和歷史偏差檢查仍保留於 `GPN-Star_ClinVar_PP3_BP4_Steps01-27E_Research_Log_Audited_v0.4_20261030.md`。
