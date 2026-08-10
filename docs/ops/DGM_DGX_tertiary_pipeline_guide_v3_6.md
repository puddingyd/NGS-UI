# 臨床三級分析 Pipeline 使用說明

**版本：v3.6**
**更新日期：2026-08-03**
**負責人：林伯昱（p88124019@gs.ncku.edu.tw）**

> v3.6 更新：新增 DGX-2 執行方式（檔名移除 DGM，本說明同時適用 DGM 與 DGX-2）、ClinGen 專家判讀對照（ERepo）、`--academic_dbnsfp`（dbNSFP 5.3a + 新 in-silico
> 工具）、PVS1 改為 ClinGen SVI 決策樹分級、新增 DGX-2 profile（含 GPU lock）、
> Pangolin 可切 CPU/GPU。SNV/indel 表由 65 欄擴充為 **81 欄**（新欄位一律附加在最後，
> 既有欄位位置不變）。
>
> v3.5 更新：新增 `STRAND_BIAS` 欄位（SNV/indel 表）、DRAGEN 原生 ploidy QC
> （`00_prepare/{SAMPLE_ID}.ploidy_qc.txt`），並修正 DRAGEN mito 已納入 annotation 的說明。

---

## 目前開發進度

### ✅ Phase 1 已完成（可使用）

| 步驟 | 功能 | 狀態 |
|------|------|------|
| PREPARE_VCF | NCKUH ensemble VCF 前處理（CALLERS tag、PASS 過濾）| ✅ 測試通過 |
| PREPARE_VCF_DRAGEN | DRAGEN VCF 前處理（CALLERS=DRAGEN、chrM 分流、自動建 tabix index）| ✅ 測試通過 |
| VEP_ANNOTATE | VEP 115 annotation（dbNSFP、LOFTEE、ClinVar、gnomAD、1000G）| ✅ 測試通過 |
| PANGOLIN_SCORE | Splice variant GPU inference | ✅ 測試通過 |
| PARSE_CSQ | VEP CSQ 解析 + ClinVar lookup + 輸出 TSV（61 欄，含 STRAND_BIAS、HGNC_ID）| ✅ 測試通過 |
| ACMG_CLASSIFY | ACMG/AMP Phase 1 自動分類（ClinGen SVI 2022）| ✅ 測試通過 |
| MITO_ANNOTATE | mtDNA annotation（VEP 輕量 + gnomAD mito v3.1（CC0）+ ClinVar，NCKUH/DRAGEN 雙支援）| ✅ 測試通過 |
| STR_ANNOTATE | STR threshold 分類（STRchive，NCKUH GangSTR/DRAGEN ExpansionHunter）| ✅ 測試通過 |

### ✅ Phase 2–3 已完成

| 步驟 | 功能 | 狀態 |
|------|------|------|
| ANNOTSV CNV/SV | AnnotSV 3.5.10 annotation（NCKUH WES/WGS + DRAGEN CNV/SV）| ✅ 測試通過 |
| PGx | PharmCAT 3.2.0 + StellarPGx（CYP2D6）+ OptiType 1.3.5（HLA-A/B）| ✅ 測試通過 |

### 🔲 後續 Phase（開發中）

| Phase | 功能 |
|-------|------|
| Phase 4 | WhatsHap phasing、Evo2 non-coding score、報告產生 |
| Phase 5 | Phenotype module（Exomiser + LIRICAL）、ROH |

---

## 執行方式

### 執行環境

目前有兩台機器可以跑，差別只在 `-profile` 與路徑：

| 機器 | profile | 程式碼位置 | 輸出位置 | 備註 |
|------|---------|-----------|---------|------|
| **DGM Server** | `dgm` | `/home/pipeline/tertiary_code` | `/home/pipeline/tertiary_output` | 32 核、獨佔 GPU |
| **DGX-2** | `dgx` | `/datalake_Intermediate/pipeline/tertiary_code` | `/datalake_Intermediate/pipeline/nextflow_output` | 48 核、V100×6，**與二級共用 GPU → 自動搶卡** |

下面的指令以 **DGM** 為例；要在 DGX-2 跑，把 `-profile dgm` 換成 `-profile dgx`、
路徑換成上表的 DGX-2 欄位即可（見文末「在 DGX-2 執行」）。

### 環境準備

```bash
# DGM
source /home/pipeline/pipeline_code/DGM_NGS2ndAnalysis.sh

# DGX-2：登入後直接用 nextflow（不需 source）
ssh n101569@10.11.33.75
```

### Sample Sheet 格式

從 v3.1 開始改為 **sample sheet 批次輸入**，支援一次跑多個樣本。

建立 CSV 檔案（例如 `/home/pipeline/samplesheet.csv`）：

```csv
sample_id,pipeline_type,input_dir,seq_type,hpo
26WE0001,nckuh,/home/pipeline/nextflow_output/26WE0001/26WE0001,WES,
26WE0002,nckuh,/home/pipeline/nextflow_output/26WE0002/26WE0002,WES,HP:0001250|HP:0001263
26WG0001,dragen,/home/datalake_Raw/Novaseq/20260428_LH00873_0015_B23NG3WLT4,WGS,
```

**欄位說明：**

| 欄位 | 必填 | 說明 |
|------|------|------|
| `sample_id` | ✅ | 樣本 ID，不可重複 |
| `pipeline_type` | ✅ | `nckuh`（NCKUH 二級分析）或 `dragen`（Illumina DRAGEN）|
| `input_dir` | ✅ | 輸入資料夾（見下方路徑規則）|
| `seq_type` | ✅ | `WES` 或 `WGS` |
| `hpo` | ❌ | HPO term，多個用 `\|` 分隔，無則留空 |

**路徑規則（pipeline 自動從 input_dir + sample_id 組合路徑）：**

```
nckuh：{input_dir}/04_snv_indel/{sample_id}.ensemble.fixed.vcf.gz
dragen：{input_dir}/vcf.gz/{sample_id}.hard-filtered.vcf.gz
```

> **注意：** 同一個 sample sheet 的 `pipeline_type` 必須一致。若要在同一個 sample sheet 混放不同類型，執行時加 `--pipeline_type` 過濾（見下方進階用法）。

---

### 執行指令

#### NCKUH 批次執行 (如果裡面全部都來自自己的pipeline)

```bash
nextflow -c /home/pipeline/tertiary_code/nextflow_tertiary.config \
    run /home/pipeline/tertiary_code/main_tertiary.nf \
    -profile dgm \
    --samplesheet /home/pipeline/samplesheet_nckuh.csv \
    --out_dir /home/pipeline/tertiary_output \
    -resume
```

#### DRAGEN 批次執行 (如果裡面全部都來自dragen)

```bash
nextflow -c /home/pipeline/tertiary_code/nextflow_tertiary.config \
    run /home/pipeline/tertiary_code/main_tertiary.nf \
    -profile dgm \
    --samplesheet /home/pipeline/samplesheet_dragen.csv \
    --out_dir /home/pipeline/tertiary_output \
    -resume
```

> **DRAGEN 注意事項：**
> - `input_dir` 指向 DRAGEN 輸出的**上層資料夾**（含多個樣本的 VCF），pipeline 會自動用 `{sample_id}.hard-filtered.vcf.gz` 搜尋
> - `.tbi` index 如果不存在，pipeline 會自動建立，不需要手動準備
> - chrM variant 會自動分流到獨立的 mito VCF（`00_prepare/{sample_id}.mito_for_annotation.vcf.gz`），
>   並**進入** mito annotation pipeline → `04_mito/{sample_id}.mito.tsv`
> - 若 DRAGEN 輸出含原生 `{sample_id}.ploidy.vcf.gz`，pipeline 會自動產生性別/倍體 QC
>   `00_prepare/{sample_id}.ploidy_qc.txt`（找不到則 warn 後跳過，不影響其他分析）

#### PGx（Pharmacogenomics）選項

PGx module 預設開啟（`run_pgx = true`，`run_pgx_cyp2d6 = true`，`run_pgx_hla = true`）。

**速度說明：**

| Module | 大約時間（WGS）| 說明 |
|--------|--------------|------|
| GATK gVCF（PGX_GVCF）| 3-8 分鐘 | 只跑 1207 個 PGx 位點，比全基因組快很多 |
| StellarPGx（CYP2D6）| 15-30 分鐘 | Graphtyper2 graph re-genotyping，不可加速 |
| OptiType（HLA-A/B）| 5-10 分鐘 | razers3 比對 + ILP solver |
| PharmCAT | 2-5 分鐘 | VCF preprocessing + Java solver |
| MT-RNR1 mpileup | < 1 分鐘 | bcftools mpileup 只跑 chrM 3 個位點 |

> ⚠️ 若追求速度（如測試或緊急出報告），可暫時關閉 PGx，等其他分析完成後再補跑：

```bash
# 關閉全部 PGx（最快，適合第一輪測試）
nextflow -c /home/pipeline/tertiary_code/nextflow_tertiary.config \
    run /home/pipeline/tertiary_code/main_tertiary.nf \
    -profile dgm \
    --samplesheet /home/pipeline/samplesheet_nckuh.csv \
    --out_dir /home/pipeline/tertiary_output \
    --run_pgx false \
    -resume

# 只關閉 StellarPGx（保留 PharmCAT + HLA，速度較快）
nextflow -c /home/pipeline/tertiary_code/nextflow_tertiary.config \
    run /home/pipeline/tertiary_code/main_tertiary.nf \
    -profile dgm \
    --samplesheet /home/pipeline/samplesheet_nckuh.csv \
    --out_dir /home/pipeline/tertiary_output \
    --run_pgx_cyp2d6 false \
    -resume

# 補跑 PGx（--resume 讓其他 module 從 cache 讀取）
nextflow -c /home/pipeline/tertiary_code/nextflow_tertiary.config \
    run /home/pipeline/tertiary_code/main_tertiary.nf \
    -profile dgm \
    --samplesheet /home/pipeline/samplesheet_nckuh.csv \
    --out_dir /home/pipeline/tertiary_output \
    -resume
```

**WES 注意：** WES 樣本無 BAM，StellarPGx 和 OptiType 自動跳過，只跑 PharmCAT（VCF-only 模式），速度快很多。

#### `--academic_dbnsfp`（換用 dbNSFP 5.3a，預設關閉）

預設用 **dbNSFP 4.9c**，裡面的工具全部可商用。加上 `--academic_dbnsfp true` 會改用
**dbNSFP 5.3a**，並多抓四個預測工具：**REVEL、MutPred2、VEST4、CADD_phred**。

> ⚠️ 這四個工具多為「**學術免費、商業使用需另行授權**」（CADD 尤其明確），所以不放在預設路徑。
> 收費臨床服務請維持關閉。

```bash
nextflow -c /home/pipeline/tertiary_code/nextflow_tertiary.config \
    run /home/pipeline/tertiary_code/main_tertiary.nf \
    -profile dgm \
    --academic_dbnsfp true \
    --samplesheet /home/pipeline/samplesheet_nckuh.csv \
    --out_dir /home/pipeline/tertiary_output \
    -resume
```

**兩種模式的差別：**

| | 預設（4.9c）| `--academic_dbnsfp true`（5.3a）|
|---|---|---|
| ACMG 計分用的族群頻率 | gnomAD 2.1.1 exomes（整體）| gnomAD 2.1.1 exomes（non_cancer 子集）|
| P-KNN（GUI 排序主訊號）| 有 | 有，且覆蓋更完整（P-KNN 原生就是 5.3 產生）|
| REVEL / MutPred2 / VEST4 / CADD | 空值 `.` | 有值 |
| gnomAD 4.1 參考欄位 | 空值 `.` | 有值 |
| `DBNSFP_VERSION` 欄 | `4.9c` | `5.3a` |

**重點：ACMG 判讀基準不會因為換版而改變。** 兩種模式的 ACMG 都用 gnomAD 2.1.1，新加的工具與
gnomAD 4.1 **只是參考欄位，不參與計分**。執行時 banner 會印出實際使用的 dbNSFP 檔案路徑，
可用來核對。

#### Pangolin 的 GPU / CPU

Pangolin 是唯一需要 GPU 的步驟。沒有 GPU 的機器加上 `--use_gpu_pangolin false` 即可改走 CPU
（較慢，結果相同），其他步驟完全不受影響。

#### 在 DGX-2 執行

DGX-2 與二級分析共用 reference、容器與 GPU，所以：

- **會自動搶卡**：`-profile dgx` 已啟用 GPU lock，每個 Pangolin 自己搶一張空閒 V100，跑完
  歸還（即使中途失敗也會還）。**不需要也不要手動指定卡號**。
- **多樣本會並行**：最多 6 個 Pangolin 同時跑（對應 6 張 V100）。單一樣本只會用到一張。
- **DRAGEN 原始資料**在 `/datalake_Raw`，profile 已掛載，`input_dir` 直接寫該路徑即可。

```bash
ssh n101569@10.11.33.75

nextflow -c /datalake_Intermediate/pipeline/tertiary_code/nextflow_tertiary.config \
    run /datalake_Intermediate/pipeline/tertiary_code/main_tertiary.nf \
    -profile dgx \
    --samplesheet /datalake_Intermediate/pipeline/samplesheet_nckuh.csv \
    --out_dir /datalake_Intermediate/pipeline/nextflow_output \
    -resume
```

其餘所有選項（`--pipeline_type`、`--run_pgx`、`--academic_dbnsfp` …）用法與 DGM 完全相同。

#### 進階：用 `--pipeline_type` 過濾混合 sample sheet

如果 sample sheet 混有 `nckuh` 和 `dragen` 兩種，可以用 `--pipeline_type` 指定這次只跑哪種，不符合的 row 會 warn 後跳過：

```bash
# 只跑 dragen 的 row
nextflow -c /home/pipeline/tertiary_code/nextflow_tertiary.config \
    run /home/pipeline/tertiary_code/main_tertiary.nf \
    -profile dgm \
    --pipeline_type dragen \
    --samplesheet /home/pipeline/samplesheet_all.csv \
    --out_dir /home/pipeline/tertiary_output \
    -resume
```

---

## 輸出檔案說明

執行完成後，輸出位於：

```
/home/pipeline/tertiary_output/{SAMPLE_ID}/
```

### v3.1 輸出結構

```
{SAMPLE_ID}/
├── 00_prepare/
│   ├── {SAMPLE_ID}.snv_for_annotation.vcf.gz      ← 前處理後的 SNV（中間檔）
│   └── {SAMPLE_ID}.snv_for_annotation.vcf.gz.tbi
│   （DRAGEN 額外輸出）
│   ├── {SAMPLE_ID}.mito_for_annotation.vcf.gz     ← chrM variants（→ 04_mito）
│   ├── {SAMPLE_ID}.mito_for_annotation.vcf.gz.tbi
│   └── {SAMPLE_ID}.ploidy_qc.txt                  ← 性別/倍體 QC（DRAGEN 原生 ploidy.vcf）
├── 01_vep/
│   ├── {SAMPLE_ID}.vep.vcf.gz                     ← VEP annotation 結果（中間檔）
│   └── {SAMPLE_ID}.vep.vcf.gz.tbi
├── 02_pangolin/
│   ├── {SAMPLE_ID}.pangolin.vcf.gz                ← Splice variant 分數（中間檔）
│   └── {SAMPLE_ID}.pangolin.vcf.gz.tbi
├── 03_acmg/
│   └── {SAMPLE_ID}.snv_indel.acmg.tsv             ← ★ SNV/Indel 最終輸出（81 欄）
├── 04_mito/                                        ← ★ v3.2 新增
│   ├── {SAMPLE_ID}.mito.tsv                       ← mtDNA 輸出（21 欄）
│   └── {SAMPLE_ID}.mito.vep.vcf.gz               ← VEP 中間檔
├── 05_str/                                         ← ★ v3.2 新增
│   └── {SAMPLE_ID}.str.tsv                        ← STR 輸出（22 欄）
├── 06_cnv_sv/                                      ← ★ v3.2 新增
│   ├── {SAMPLE_ID}.cnv.annotated.tsv              ← CNV AnnotSV 輸出
│   ├── {SAMPLE_ID}.cnv.unannotated.tsv
│   ├── {SAMPLE_ID}.sv.annotated.tsv               ← SV AnnotSV 輸出
│   └── {SAMPLE_ID}.sv.unannotated.tsv
└── 07_pgx/                                         ← 藥物基因組學（PGx）
    ├── {SAMPLE_ID}.pgx.tsv                        ← ★ PGx 報告（16 欄，CPIC Level A）
    ├── {SAMPLE_ID}.pharmcat.report.json           ← PharmCAT 完整報告
    ├── {SAMPLE_ID}.outside_calls.tsv              ← PharmCAT outside calls
    ├── {SAMPLE_ID}.stellarpgx.tsv                 ← CYP2D6 diplotype（WGS only）
    └── {SAMPLE_ID}.optitype.tsv                   ← HLA-A/B/C（WGS only）
```

> **v3.2 新增：** `04_mito/`、`05_str/`、`06_cnv_sv/` 目錄。STR TSV 只包含有對應 STRchive 記錄的 locus（已知致病 STR），DRAGEN 約 53 筆，NCKUH WES 約 5 筆。

---

### 主要輸出欄位（snv_indel.acmg.tsv，81 欄）

> v3.6 起由 65 欄擴充為 81 欄。**新欄位一律附加在最後**，既有欄位的位置沒有變動，
> 用欄位編號取值的舊腳本不受影響（但新欄位的編號請以本節為準）。

#### 位置資訊（欄 1–5）
| 欄位 | 說明 |
|------|------|
| CHROM, POS, REF, ALT | 變異座標 |
| RS_ID | dbSNP rsID（如 rs72631890）|

#### Transcript 資訊（欄 6–14）
| 欄位 | 說明 |
|------|------|
| GENE | HGNC gene symbol |
| TRANSCRIPT | Ensembl transcript ID |
| TRANSCRIPT_TYPE | MANE_SELECT / MANE_PLUS_CLINICAL / CANONICAL / APPRIS_P1 / BEST_CONSEQUENCE |
| HGVS_C, HGVS_P | HGVS 命名 |
| CONSEQUENCE | VEP consequence（如 missense_variant）|
| IMPACT | HIGH / MODERATE / LOW / MODIFIER |
| EXON, INTRON | Exon/Intron 編號（格式：2/19）|

#### Caller 資訊（欄 15–23）
| 欄位 | 說明 |
|------|------|
| CALLERS | `DV+HC` / `DV` / `HC`（NCKUH）或 `DRAGEN` |
| DP_DV, AD_DV, VAF_DV | DeepVariant（DRAGEN 樣本為 DRAGEN）read depth、allelic depth、VAF |
| DP_HC, AD_HC | HaplotypeCaller read depth、allelic depth（DRAGEN 無 HC，此二欄為空）|
| ZYGOSITY | het / hom / hemizygous / unknown |
| GT_DV, GT_HC | Genotype |

#### Strand bias（欄 24）
| 欄位 | 說明 |
|------|------|
| STRAND_BIAS | 股偏警示：`PASS` / `WARN(FS=..,SOR=..)`（GATK 門檻：SNV FS>60/SOR>3.0，indel FS>200/SOR>10.0）；`.` = 無 FS/SOR（DeepVariant-only 位點）→ 需人工複核。germline 只標記、不硬刪 |

#### 族群頻率（欄 25–31）
| 欄位 | 說明 |
|------|------|
| GNOMAD_G_AF, GNOMAD_G_EAS_AF | gnomAD genome AF（全體 + 東亞）|
| GNOMAD_E_AF, GNOMAD_E_EAS_AF | gnomAD exome AF（全體 + 東亞）|
| GNOMAD_E_AF_DBNSFP, GNOMAD_E_EAS_AF_DBNSFP | dbNSFP 版本 gnomAD exome AF |
| TG_EAS_AF | 1000 Genomes 東亞次族群 AF |

#### ClinVar（欄 32–36）
| 欄位 | 說明 |
|------|------|
| CLINVAR_SIG | ClinVar 判定（Pathogenic / Likely_pathogenic / ...）|
| CLINVAR_STARS | ClinVar 星級（0–4）|
| CLINVAR_DN | ClinVar 相關疾病名稱 |
| CLINVAR_SIGCONF | ClinVar conflicting interpretations 明細 |
| CLINVAR_VARIATION_ID | ClinVar Variation ID（GUI 組 URL：`https://www.ncbi.nlm.nih.gov/clinvar/variation/{ID}/`）|

#### OMIM（欄 37）
| 欄位 | 說明 |
|------|------|
| OMIM_IDS | OMIM 疾病 ID，逗號分隔（GUI 組 URL：`https://www.omim.org/entry/{ID}`）|

#### LOFTEE（欄 38–41）
| 欄位 | 說明 |
|------|------|
| LOFTEE | HC / LC / .（LoF 信心度）|
| LOFTEE_FILTER | LOFTEE filter 原因（. = 通過）|
| LOFTEE_FLAGS | LOFTEE 額外旗標 |
| LOFTOOL | LoF tolerance score |

#### In Silico 預測（欄 42–54）
| 欄位 | 說明 |
|------|------|
| BAYESDEL_NOAF, BAYESDEL_NOAF_PRED | BayesDel（ClinGen SVI 推薦）|
| ALPHAMISSENSE, ALPHAMISSENSE_PRED | AlphaMissense |
| ESM1B, ESM1B_PRED | ESM1b 蛋白質語言模型（★ 越負越致病）|
| VARITY_R | VARITY_R（ClinGen SVI 推薦）|
| SIFT, SIFT_PRED | SIFT |
| DANN | DANN |
| PHACTBOOST | PHACTboost |
| PHYLOP100 | phyloP 100-way vertebrate conservation |
| GERP | GERP++ conservation |

#### P-KNN（欄 55–56）
| 欄位 | 說明 |
|------|------|
| PKNN_LLR | P-KNN joint calibration LLR（missense only）|
| PKNN_EVIDENCE | PP3/BP4 evidence 強度 |

#### Splice（欄 57–58）
| 欄位 | 說明 |
|------|------|
| PANGOLIN_SCORE | Pangolin splice score（最大值）|
| PANGOLIN_DETAIL | Pangolin 完整輸出字串 |

#### 蛋白質（欄 59–60）
| 欄位 | 說明 |
|------|------|
| DOMAINS | 蛋白質 domain 資訊 |
| SWISSPROT | UniProt/SwissProt ID（GUI 組 URL：`https://www.uniprot.org/uniprot/{ID}`）|

#### Gene identifier（欄 61）
| 欄位 | 說明 |
|------|------|
| HGNC_ID | HGNC 基因識別碼（GUI 組 URL：`https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{ID}`）|

#### ClinGen 專家判讀對照（欄 62–64，v3.6 新增）
ClinGen Evidence Repository（ERepo）收錄各 **VCEP 專家小組**的變異判讀，含他們實際套用了哪些
ACMG criteria。以 ClinVar Variation ID 對照到我們的表。

| 欄位 | 說明 |
|------|------|
| CLINGEN_VCEP_CLASS | 專家小組的判讀結論 |
| CLINGEN_VCEP_CRITERIA | 專家實際套用的 criteria（如 `PVS1,PM2_Supporting,PP3`）|
| CLINGEN_VCEP_PANEL | 哪個 VCEP 判的 |

> **只作對照、不進計分。** ClinGen SVI 2018 建議不要用 PP5/BP6（拿他人判讀當證據會變成循環
> 論證），本 pipeline 也未實作 PP5/BP6。實測 NA12878 WES 有 77 個變異被專家判讀過。

#### dbNSFP 5.3a 專屬工具與 gnomAD 4.1（欄 65–72，v3.6 新增）
| 欄位 | 說明 |
|------|------|
| REVEL / MUTPRED2 / MUTPRED2_PRED / VEST4 / CADD_PHRED | 只有 `--academic_dbnsfp true` 才有值 |
| GNOMAD41_JOINT_AF / GNOMAD41_JOINT_EAS_AF | gnomAD 4.1 joint（exomes+genomes，~80 萬人）**參考用，不進計分** |
| DBNSFP_VERSION | 這批分數來自 `4.9c` 或 `5.3a` |

> 覆蓋率說明：這些工具只涵蓋 **missense SNV**（dbNSFP 本身就只收 nsSNV），synonymous /
> intron / UTR / indel 一律為 `.`，屬正常。實測在 missense SNV 母體內：CADD 99.5%、
> P-KNN 99.5%、MutPred2 98.9%、AlphaMissense 98.8%、REVEL 95.9%。

#### PVS1 決策樹輸入（欄 73–74，v3.6 新增）
| 欄位 | 說明 |
|------|------|
| NMD | VEP NMD plugin：預測**逃過** NMD 時才有值；空值代表會被 NMD 降解 |
| PROTEIN_POSITION | `123/456` 格式，用來算截斷掉多少比例的蛋白質 |

#### ACMG 分類（欄 75–78）
| 欄位 | 說明 |
|------|------|
| ACMG_CRITERIA | 觸發的所有 criteria，逗號分隔（如 `PVS1,PM2_Supporting`）|
| ACMG_SCORE | 數值化分數（越高越可能致病，用於排序）|
| ACMG_CLASS | `Pathogenic` / `Likely_Pathogenic` / `VUS` / `Likely_Benign` / `Benign` |
| ACMG_NOTES | 觸發原因說明，含數值依據（如 `PVS1:LOFTEE=HC,gene=MECP2,HI=3`）|

#### 與專家判讀的一致性 + PVS1 分級（欄 79–81，v3.6 新增）
| 欄位 | 說明 |
|------|------|
| CLINGEN_AGREEMENT | `AGREE`（同一級）/ `DIFFER_TIER`（方向相同、強度不同）/ `DIFFER`（方向不同 → **優先人工複核**）/ `.`（無專家判讀）|
| PVS1_STRENGTH | `PVS1` / `PVS1_Strong` / `PVS1_Moderate` / `PVS1_Supporting` / `.` |
| PVS1_REASON | 決策樹走了哪條分支，供人工複核 |

**PVS1 為什麼會分級？**（ClinGen SVI, Abou Tayoun 2018）

舊版只要「LOFTEE 高信心 + 該基因 LoF 致病」就給滿分 8 分。但 SVI 指出這樣過度樂觀：即使是
predicted LoF，若**逃過 NMD**、只截掉蛋白尾端、或落在非關鍵區域，實際影響有限，應該降級。

| 情境 | 判定 | 分數 |
|------|------|------|
| nonsense / frameshift / splice±1,2 且**會被 NMD 降解** | `PVS1` | 8 |
| 逃過 NMD，但落在**已知功能域** | `PVS1_Strong` | 4 |
| 逃過 NMD，且**移除 >10% 蛋白質** | `PVS1_Strong` | 4 |
| 逃過 NMD，只截尾一點點 | `PVS1_Moderate` | 2 |
| 起始密碼子 `start_lost`（SVI 上限） | `PVS1_Moderate` | 2 |
| LOFTEE 低信心、或該基因 LoF 非致病機轉（ClinGen HI≠3）| 不觸發 | 0 |

> **這會改變判讀結果**：逃過 NMD 又只截尾的變異，可能從 `Likely_Pathogenic` 降到 `VUS`。
> 這正是決策樹的目的，但看到判讀與舊版不同時請以 `PVS1_REASON` 核對原因。
>
> 另外，`PVS1_STRENGTH` 有值但 `ACMG_CLASS` 是 `Benign` **不是矛盾** —— BA1（族群 AF > 5%）
> 是 stand-alone benign，會蓋過一切；`PVS1_STRENGTH` 仍獨立記錄決策樹的評估（「是 LoF 但族群
> 常見」本身就是有用的資訊）。

---

### mito.tsv 欄位（21 欄）

| 欄位 | 說明 |
|------|------|
| CHROM, POS, REF, ALT | 位置 |
| GENE | MT gene（MT-ND1, MT-RNR2 等）|
| HGVS_C, HGVS_P | HGVS 命名（upstream variant 無值）|
| CONSEQUENCE, IMPACT, BIOTYPE | VEP 後果 |
| GENOTYPE, DP | 樣本 GT 和 read depth |
| AF_SAMPLE | Heteroplasmy level（FORMAT/AF，0-1）|
| CLINVAR_SIG, CLINVAR_DN | ClinVar 致病性和疾病名稱 |
| GNOMAD_MITO_AF_HOM | gnomAD mito v3.1 同質性 AF（heteroplasmy ≥ 0.95，CC0）|
| GNOMAD_MITO_AF_HET | gnomAD mito v3.1 異質性 AF（heteroplasmy 0.10-0.95）|
| GNOMAD_MITO_AN | gnomAD mito v3.1 總樣本數（56,434）|
| CLINVAR_VARIATION_ID | ClinVar Variation ID |
| OMIM_IDS | OMIM 疾病 ID |
| PIPELINE | nckuh / dragen |

> **注意：** NCKUH 只輸出 PASS variant；DRAGEN 輸出所有 chrM variant（含 non-PASS），請依 FILTER 欄位篩選。
> **License：** 舊版使用 MITOMAP（CC BY-NC，收費臨床屬商業用途）。現改用 gnomAD mito v3.1（CC0）+ ClinVar（開放）。

### str.tsv 欄位（22 欄）

| 欄位 | 說明 |
|------|------|
| CHROM, POS, END | 位置 |
| STR_ID | STRchive locus ID（如 HD_HTT）|
| GENE | 基因名稱 |
| MOTIF | repeat 單元（如 CAG）|
| LOCUS_STRUCTURE | repeat 結構（如 (CAG)*CAACAG(CCG)*）|
| TYPE | locus 位置（Coding / 5' UTR 等）|
| REPCN_A1, REPCN_A2 | 兩個 allele 的 repeat count |
| DP | Read depth（NCKUH: FORMAT/DP；DRAGEN: FORMAT/LC locus coverage）|
| REPCI | Confidence interval |
| BENIGN_MIN, BENIGN_MAX | 正常範圍（min 用於缺失型 locus 如 VWA1）|
| PATHOGENIC_MIN, PATHOGENIC_MAX | 致病範圍 |
| INTERMEDIATE_MIN, INTERMEDIATE_MAX | 中間範圍（部分 locus 有）|
| CLASSIFICATION | normal / intermediate / borderline / pathogenic / no_threshold |
| DISEASE | 疾病名稱 |
| INHERITANCE | AD / AR / XL |
| PIPELINE | nckuh / dragen |

---

### cnv.annotated.tsv / sv.annotated.tsv 欄位（AnnotSV both 模式）

AnnotSV 輸出兩種行（`Annotation_mode` 欄位）：
- **full**：整個 SV 的 annotation（每個 SV 一行）
- **split**：按 gene 拆分（一個 SV 跨多個 gene 時，每個 gene 一行）

臨床判讀建議：先看 `full` 行的 `ACMG_class`，再用 `split` 行展開各基因詳情。

#### 關鍵欄位說明

| 欄位 | 說明 |
|------|------|
| AnnotSV_ID | 格式：chrom_start_end_SVTYPE_n |
| SV_chrom, SV_start, SV_end | SV 座標 |
| SV_length | SV 長度（bp），deletion 為負值 |
| SV_type | DEL / DUP / INS / INV / BND |
| Samples_ID | 樣本 ID（CNVkit BED 輸入時若無 sample column 會是 NA）|
| Annotation_mode | full / split |
| Gene_name | 覆蓋的基因（sorted by genomic coordinate）|
| ACMG_class | 1=Benign / 2=Likely benign / 3=VUS / 4=Likely pathogenic / 5=Pathogenic / NA |
| AnnotSV_ranking_score | ACMG/ClinGen 評分（≥0.99=Pathogenic，≤-0.99=Benign）|
| AnnotSV_ranking_criteria | 評分依據說明 |
| HI | ClinGen Haploinsufficiency score（3=sufficient evidence）|
| TS | ClinGen Triplosensitivity score（3=sufficient evidence）|
| OMIM_morbid | yes = 覆蓋 OMIM morbid gene |
| OMIM_phenotype | OMIM 表型（split 行）|
| OMIM_inheritance | 遺傳模式（split 行）|
| P_gain_phen, P_gain_source | 已知致病 gain（duplication）資料庫命中 |
| P_loss_phen, P_loss_source | 已知致病 loss（deletion）資料庫命中 |
| B_gain_source, B_gain_AFmax | 已知良性 gain，最大 AF |
| B_loss_source, B_loss_AFmax | 已知良性 loss，最大 AF |
| GnomAD_pLI | gnomAD LoF intolerance score（≥0.9=intolerant）|
| ExAC_delZ, ExAC_dupZ | ExAC CNV intolerance Z score |
| CytoBand | 細胞遺傳帶（如 22q11.2）|
| RE_gene | 覆蓋的 regulatory element 所調控的基因（full 行）|
| GenCC_classification | GenCC 基因-疾病證據等級 |

> **臨床過濾建議：**
> 1. 先篩 `Annotation_mode == "full"` 的行，看 `ACMG_class` 4 或 5
> 2. NCKUH WES CNV 約 112 行，DRAGEN CNV 約 348 行（PASS + non copy-neutral）
> 3. DRAGEN SV 約 30,068 行（含大量 INS/DEL，臨床先關注 class 4-5 的 DEL/DUP）



---

### pgx.tsv 欄位（16 欄）

`07_pgx/{SAMPLE_ID}.pgx.tsv` — 臨床用藥基因組學報告，只包含 CPIC Level A 基因。

#### 欄位說明

| 欄位 | 說明 |
|------|------|
| SAMPLE_ID | 樣本 ID |
| PIPELINE | `nckuh` 或 `dragen` |
| GENE | 基因名稱（CYP2D6、CYP2C19、HLA-B 等 CPIC Level A 基因）|
| DIPLOTYPE | Star allele diplotype（如 `*1/*5`）或 HLA allele（如 `*03:01/*68:01`）|
| ACTIVITY_SCORE | CYP2D6、CYP2C9 等有 activity score 的基因填入數值；其他基因為 `.` |
| PHENOTYPE | 代謝表型（Normal Metabolizer / Intermediate Metabolizer / Poor Metabolizer / Ultrarapid Metabolizer）|
| DRUG | 藥物名稱（英文小寫，如 `clopidogrel`）|
| GUIDELINE_SOURCE | 指引來源（`CPIC` / `DPWG` / `FDA`）|
| RECOMMENDATION | 劑量調整或替代藥物建議（來自 PharmCAT report.json，HTML entities 已解碼）|
| IMPLICATION | 基因型對藥物代謝影響說明 |
| CPIC_LEVEL | CPIC 證據等級（`Strong` / `Moderate` / `.`）|
| DPWG_LEVEL | DPWG 證據等級（`Strong` / `Moderate` / `.`）|
| OUTSIDE_CALLER | outside call 來源（`PharmCAT-outside`：StellarPGx/OptiType 提供）|
| MTRN1_RISK | MT-RNR1 aminoglycoside 風險（`HIGH` / `LOW` / Unknown）。`LOW` 需要 mpileup 證實 chrM:827/1494/1555 至少一個位點 `DP ≥ 10`；三個位點都沒深度時**不輸出 MT-RNR1 列**（＝Unknown，代表「沒測到」而非「陰性」）|
| NOTES | 補充說明（heteroplasmy AF、ClinVar sig；MT-RNR1 列會標明哪些位點已確認覆蓋、哪些未評估）|
| EVIDENCE_STRENGTH | 整體證據強度（`Strong` / `Moderate` 等）|

> **DRAGEN 交叉註記（僅 DRAGEN 樣本，v3.5）：** `NOTES` 欄會附上 DRAGEN 原生 PGx 判讀
> （`other/{sample}/germline_seq/{sample}.targeted.json`）與我們的比對，格式
> `DRAGEN <判定>: <DRAGEN 原始 genotype>`：
> - `一致`：正規化後與我們的 diplotype 相同（含落在 DRAGEN 模糊多重解 `;` 候選集內）；
> - `不一致`：同命名系統但不同（如 CYP2D6、DPYD）→ 建議人工複核；
> - `未比對`：命名系統不同（star vs HGVS/rs），不下判定、僅附 DRAGEN 原文。
>
> reference 跨寫法（`*1`／`Reference`／`B(wildtype)`）視為相同；DRAGEN 沒有的基因（如 HLA）不註記。
> **欄位不變**，只是把 DRAGEN 結果寫進既有的 `NOTES`。找不到 `targeted.json` → 該樣本保留原
> `NOTES`、不報錯。

#### CPIC Level A 基因清單

| 基因 | Outside caller | 主要臨床意義 |
|------|--------------|-------------|
| CYP2D6 | StellarPGx（WGS）| 抗憂鬱劑、止痛藥、抗精神病藥 |
| CYP2C19 | PharmCAT（VCF）| clopidogrel、PPI、抗憂鬱劑 |
| CYP2C9 | PharmCAT（VCF）| warfarin、NSAIDs、phenytoin |
| DPYD | PharmCAT（VCF）| fluoropyrimidine（5-FU）毒性 |
| TPMT | PharmCAT（VCF）| thiopurine（azathioprine）毒性 |
| NUDT15 | PharmCAT（VCF）| thiopurine（azathioprine）毒性（亞洲族群重要）|
| SLCO1B1 | PharmCAT（VCF）| simvastatin 肌肉毒性 |
| HLA-A | OptiType（WGS）| abacavir 過敏（*31:01 positive）|
| HLA-B | OptiType（WGS）| abacavir（*57:01）、carbamazepine（*15:02）、allopurinol（*58:01）過敏 |
| UGT1A1 | PharmCAT（VCF）| irinotecan 毒性 |
| G6PD | PharmCAT（VCF）| 多種藥物溶血風險 |
| MT-RNR1 | mito pipeline + BAM mpileup | aminoglycoside 致聾風險（CPIC Level A）|

> **WES 注意：** WES 有 BAM，會跑 GATK gVCF（大幅改善 PharmCAT 結果），但 StellarPGx（CYP2D6）和 OptiType（HLA-A/B）自動跳過（WGS only）。WES Called 基因數約 17/23。

> **HLA-B 準確度限制：** OptiType 在 reads 不夠多時容易將 heterozygous 誤判為 homozygous。NA12878 ground truth 為 `B*07:02/B*40:02`，pipeline 目前 call 為 `B*08:01/B*08:01`。臨床最重要的 `*57:01`（abacavir）和 `*58:01`（allopurinol）negative 結果正確，臨床安全性無虞。如需更高 HLA 準確度，建議用專用 HLA typing 服務。

#### 已知不準確基因（WGS/WES 均適用）

| 基因 | 原因 |
|------|------|
| CYP2D6（WES）| WES capture 不完整，VCF-based calling 準確度低 |
| HLA-A/B（WES）| 需要 BAM，WES 無 BAM |
| G6PD | X-linked，hemizygous calling 需要額外處理 |

---

### PGx 輸出驗證指令

```bash
SAMPLE_ID=26WG0001
PGX=/home/pipeline/tertiary_output/${SAMPLE_ID}/07_pgx/${SAMPLE_ID}.pgx.tsv

# 確認輸出目錄
ls -lh /home/pipeline/tertiary_output/${SAMPLE_ID}/07_pgx/

# CYP2D6 diplotype + activity score
grep "CYP2D6\|GENE" $PGX | cut -f3-6 | head -5
# WGS 預期：CYP2D6  *X/*X  1.0（或其他數值）  Intermediate/Normal/Poor Metabolizer

# HLA-A/B（WGS only）
grep "HLA" $PGX | cut -f3,4,6 | head -5
# WGS 預期：HLA-B  *07:02/*73:01  *57:01 negative; *15:02 negative

# 總筆數（WGS 約 140 行，WES 約 100 行）
wc -l $PGX

# 確認 ACTIVITY_SCORE 有值（CYP2D6/CYP2C9 應有數值，其他為 .）
awk -F'\t' 'NR>1 && $5!="."' $PGX | cut -f3,5 | sort -u

# 確認 outside caller 欄位
awk -F'\t' 'NR>1 && $13!=""' $PGX | cut -f3,13 | sort -u
```

---

## 驗證指令

收到新版輸出後，請執行以下指令確認正確：

```bash
# 設定路徑（替換 SAMPLE_ID）
SAMPLE_ID=26WE0001
TSV=/home/pipeline/tertiary_output/${SAMPLE_ID}/03_acmg/${SAMPLE_ID}.snv_indel.acmg.tsv

echo "========================================="
echo "Step 1：欄位數（應為 81）"
echo "========================================="
head -1 $TSV | tr '\t' '\n' | wc -l

echo ""
echo "========================================="
echo "Step 2：最後 10 個欄位名稱"
echo "（確認 HGNC_ID 和 4 個 ACMG 欄位存在）"
echo "========================================="
head -1 $TSV | tr '\t' '\n' | tail -10

echo ""
echo "========================================="
echo "Step 3：總 variant 數"
echo "========================================="
awk 'NR>1' $TSV | wc -l

echo ""
echo "========================================="
echo "Step 4：CALLERS 欄位確認（第 15 欄）"
echo "（NCKUH 應為 DV+HC/DV/HC；DRAGEN 應為 DRAGEN）"
echo "========================================="
awk -F'\t' 'NR>1 {print $15}' $TSV | sort | uniq -c | sort -rn | head -5

echo ""
echo "========================================="
echo "Step 5：ClinVar 分布（應有非 . 的值）"
echo "========================================="
awk -F'\t' 'NR>1 {print $32}' $TSV | sort | uniq -c | sort -rn | head -10

echo ""
echo "========================================="
echo "Step 6：ACMG 分類分布"
echo "========================================="
awk -F'\t' 'NR>1 {print $77}' $TSV | sort | uniq -c | sort -rn

echo ""
echo "========================================="
echo "Step 7：P/LP variant 詳細資訊"
echo "========================================="
awk -F'\t' 'NR>1 && ($77=="Pathogenic" || $77=="Likely_Pathogenic")' $TSV \
    | cut -f1,2,6,11,75,76,77,78,80 | head -10   # GENE/HGVS/ACMG 4 欄 + PVS1_STRENGTH

echo ""
echo "========================================="
echo "驗證完成"
echo "========================================="
```

### Mito 輸出確認

```bash
SAMPLE_ID=26WE0001
MITO=/home/pipeline/tertiary_output/${SAMPLE_ID}/04_mito/${SAMPLE_ID}.mito.tsv

# 確認輸出存在
ls -lh /home/pipeline/tertiary_output/${SAMPLE_ID}/04_mito/

# 總 variant 數（NCKUH 約 30-60 筆，DRAGEN 約 40-60 筆）
wc -l $MITO

# gnomAD mito 命中筆數（$14 = GNOMAD_MITO_AF_HOM）
awk -F'\t' 'NR>1 && $14!="."' $MITO | wc -l

# CONSEQUENCE 分布
awk -F'\t' 'NR>1 {print $8}' $MITO | sort | uniq -c | sort -rn
```

### STR 輸出確認

```bash
STR=/home/pipeline/tertiary_output/${SAMPLE_ID}/05_str/${SAMPLE_ID}.str.tsv

# 確認輸出存在
ls -lh /home/pipeline/tertiary_output/${SAMPLE_ID}/05_str/

# STRchive 命中筆數（DRAGEN: 約 53；NCKUH WES: 約 5）
wc -l $STR

# 分類分布
awk -F'\t' 'NR>1 {print $19}' $STR | sort | uniq -c | sort -rn

# Pathogenic / Intermediate 詳細資訊
awk -F'\t' 'NR>1 && ($19=="pathogenic" || $19=="intermediate")' $STR \
    | cut -f4,5,9,10,15,16,19,20
```

### CNV/SV 輸出確認

```bash
SAMPLE_ID=26WG0001
CNV_DIR=/home/pipeline/tertiary_output/${SAMPLE_ID}/06_cnv_sv

# 確認輸出存在
ls -lh $CNV_DIR/

# CNV 行數（NCKUH WES 約 100 行；NCKUH WGS 約數百行；DRAGEN 約數百行）
wc -l $CNV_DIR/${SAMPLE_ID}.cnv.annotated.tsv

# SV 行數（NCKUH/DRAGEN WGS 可能超過 10 萬行，正常）
wc -l $CNV_DIR/${SAMPLE_ID}.sv.annotated.tsv

# CNV：確認有 AnnotSV 分類欄位
head -1 $CNV_DIR/${SAMPLE_ID}.cnv.annotated.tsv | tr '	' '
' | grep -n "AnnotSV_ranking\|ACMG\|Pathogenic"

# SV：確認 unannotated 筆數（沒有對應基因的 SV，通常為少數）
wc -l $CNV_DIR/${SAMPLE_ID}.sv.unannotated.tsv
```

### DRAGEN ploidy QC 確認（僅 DRAGEN 樣本）

```bash
SAMPLE_ID=VAL-10
cat /home/pipeline/tertiary_output/${SAMPLE_ID}/00_prepare/${SAMPLE_ID}.ploidy_qc.txt
# 預期：estimated_sex_karyotype 與 samplesheet 宣告一致 → sex_check: OK
# 每條 contig 的 NDC 應 ~1.0；偏離過多的 contig 會列在 WARNINGS（疑似非整倍體，需人工確認）
```

### PGx 輸出確認

```bash
SAMPLE_ID=26WG0001
PGX_DIR=/home/pipeline/tertiary_output/${SAMPLE_ID}/07_pgx

# 確認輸出檔案
# WGS：6 個（pgx.tsv, pharmcat.report.json, outside_calls.tsv,
#            stellarpgx.tsv, optitype.tsv, mtrn1.vcf.gz）
# WES：4 個（pgx.tsv, pharmcat.report.json, outside_calls.tsv, mtrn1.vcf.gz）
ls -lh $PGX_DIR/

# pgx.tsv 總行數（WGS 約 180-200 行，WES 約 130-150 行）
wc -l $PGX_DIR/${SAMPLE_ID}.pgx.tsv

# Called 基因數（PharmCAT）
# WGS 預期：22/23（只有 MT-RNR1 不在 PharmCAT 內）
# WES 預期：17/23（CYP2D6/HLA-A/HLA-B/IFNL3/MT-RNR1/VKORC1 可能 Unknown）
python3 -c "
import json
with open('$PGX_DIR/${SAMPLE_ID}.pharmcat.report.json') as f:
    d = json.load(f)
genes = d.get('genes', {})
called = [g for g, v in genes.items()
          if v.get('sourceDiplotypes', [{}])[0].get('label', '') not in
          ('Unknown/Unknown', 'Unknown', 'No Result', '')]
print(f'PharmCAT called: {len(called)}/{len(genes)} genes')
print('Called:', sorted(called))
"

# CYP2D6（WGS：StellarPGx outside call；WES：Unknown）
echo "--- CYP2D6 ---"
grep "CYP2D6" $PGX_DIR/${SAMPLE_ID}.pgx.tsv | cut -f3-6 | sort -u | head -3

# HLA-A/HLA-B（WGS only）
echo "--- HLA ---"
grep "HLA" $PGX_DIR/${SAMPLE_ID}.pgx.tsv | cut -f3,4,6 | sort -u | head -5

# MT-RNR1（所有樣本，mpileup 確認）
echo "--- MT-RNR1 ---"
grep "MT-RNR1" $PGX_DIR/${SAMPLE_ID}.pgx.tsv | cut -f3,4,13,14,15
# 預期：Reference  LOW  Coverage confirmed at chrM positions: 827,1494,1555
# 若有風險：hgvs（m.1555A>G 等）  HIGH

# 臨床警示
echo "--- B*57:01 / B*58:01 ---"
grep "HLA-B" $PGX_DIR/${SAMPLE_ID}.pgx.tsv | cut -f3,4,6 | head -3
```

> **PGx 結果解讀：** HLA-B `*57:01 negative` 表示無 abacavir 過敏風險；`*58:01 negative` 表示無 allopurinol 嚴重過敏風險；`*15:02 negative` 表示無 carbamazepine 嚴重皮膚反應風險。這三個是台灣臨床最常用的 CPIC Level A HLA 警示。

---

## 常見問題

**Q：`-resume` 沒有效果，從頭重跑？**

代表 work 目錄不存在或被清除。加上 `-resume` 即可讓已完成的 process 不重跑。

**Q：Pangolin 輸出是空的？**

該樣本可能沒有 splice region variant，這是正常現象。`PANGOLIN_SCORE` 欄位全部為 `.`。

**Q：CLINVAR_SIG 全部都是 `.`？**

請確認 `/home/pipeline/reference/hg38/tertiary/clinvar/` 下的 ClinVar VCF contig 格式是否為 `chr1`（而非 `1`）。NCBI 官方下載的 ClinVar VCF 預設不含 `chr` 前綴，需要用 `bcftools annotate --rename-chrs` 轉換後才能與 hg38 pipeline 正確 match。

**Q：DRAGEN 樣本找不到 VCF？**

請確認 `input_dir` 是 DRAGEN run 的**上層資料夾**，pipeline 會自動去 `{input_dir}/vcf.gz/` 下找 `{sample_id}.hard-filtered.vcf.gz`。例如：

```
input_dir = /home/datalake_Raw/Novaseq/20260428_LH00873_0015_B23NG3WLT4
pipeline 會去找：/home/datalake_Raw/Novaseq/20260428_LH00873_0015_B23NG3WLT4/vcf.gz/VAL-10.hard-filtered.vcf.gz
```

**Q：DRAGEN 缺少 `.tbi` index？**

不需要手動建，pipeline 的 `ADD_DRAGEN_TAG` process 會自動偵測並建立。

**Q：mito.tsv 的 GNOMAD_MITO_AF_HOM 和 AF_HET 是什麼？**

來自 gnomAD v3.1 mito（CC0，56,434 個 WGS）的本機查表。`AF_HOM` 是同質性頻率（heteroplasmy ≥ 95%），`AF_HET` 是異質性頻率（10-95%）。AF 高的 variant 通常是族群多型性。GNOMAD_MITO_AN = 56,434 表示全部樣本均有資料。

**Q：str.tsv 的筆數比預期少？**

STR TSV 只輸出有對應 STRchive 記錄的 locus（已知致病 STR，目前約 73 個 locus）。DRAGEN WGS 約命中 53 筆，NCKUH WES 因 capture 範圍限制約 5 筆，其餘 STR locus 不在輸出中。

**Q：str.tsv 出現 `no_threshold` 分類？**

部分 STRchive locus 的 pathogenic_min 設定非常高（如 FAME 相關疾病 > 100 copies），正常人的 repeat count 也可能在此範圍，屬於正常變異。`no_threshold` 表示 STRchive 對此 locus 的分類標準特殊，需人工判讀。

**Q：同一個 sample sheet 可以混放 nckuh 和 dragen 嗎？**

可以放，但執行時需要加 `--pipeline_type` 指定要跑哪種，另一種 row 會被跳過：

```bash
--pipeline_type nckuh   # 只跑 nckuh 的 row
--pipeline_type dragen  # 只跑 dragen 的 row
```

若不加 `--pipeline_type`，而 sample sheet 內有兩種 `pipeline_type`，pipeline 會報錯提示你拆開或加上過濾參數。

---
