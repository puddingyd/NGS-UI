# 臨床三級分析 Pipeline 使用說明

**版本：v3.2**
**更新日期：2026-06-06**
**負責人：林伯昱（p88124019@gs.ncku.edu.tw）**

---

## 目前開發進度

### ✅ Phase 1 已完成（可使用）

| 步驟 | 功能 | 狀態 |
|------|------|------|
| PREPARE_VCF | NCKUH ensemble VCF 前處理（CALLERS tag、PASS 過濾）| ✅ 測試通過 |
| PREPARE_VCF_DRAGEN | DRAGEN VCF 前處理（CALLERS=DRAGEN、chrM 分流、自動建 tabix index）| ✅ 測試通過 |
| VEP_ANNOTATE | VEP 115 annotation（dbNSFP、LOFTEE、ClinVar、gnomAD、1000G）| ✅ 測試通過 |
| PANGOLIN_SCORE | Splice variant GPU inference | ✅ 測試通過 |
| PARSE_CSQ | VEP CSQ 解析 + ClinVar lookup + 輸出 TSV（61 欄，含 HGNC_ID）| ✅ 測試通過 |
| ACMG_CLASSIFY | ACMG/AMP Phase 1 自動分類（ClinGen SVI 2022）| ✅ 測試通過 |
| MITO_ANNOTATE | mtDNA annotation（VEP 輕量 + MITOMAP 查表，NCKUH/DRAGEN 雙支援）| ✅ 測試通過 |
| STR_ANNOTATE | STR threshold 分類（STRchive，NCKUH GangSTR/DRAGEN ExpansionHunter）| ✅ 測試通過 |

### 🔲 後續 Phase（開發中）

| Phase | 功能 |
|-------|------|
| Phase 2 | AnnotSV CNV/SV annotation |
| Phase 3 | Exomiser/LIRICAL 表型驅動排序、PGx（Aldy）|
| Phase 4 | WhatsHap phasing、Evo2 non-coding score、報告產生 |

---

## 執行方式

### 環境準備

```bash
source /home/pipeline/pipeline_code/DGM_NGS2ndAnalysis.sh
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
> - chrM variant 會自動分流到獨立的 mito VCF（備用，目前不進入 annotation pipeline）

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
│   ├── {SAMPLE_ID}.mito_for_annotation.vcf.gz     ← chrM variants（中間檔）
│   └── {SAMPLE_ID}.mito_for_annotation.vcf.gz.tbi
├── 01_vep/
│   ├── {SAMPLE_ID}.vep.vcf.gz                     ← VEP annotation 結果（中間檔）
│   └── {SAMPLE_ID}.vep.vcf.gz.tbi
├── 02_pangolin/
│   ├── {SAMPLE_ID}.pangolin.vcf.gz                ← Splice variant 分數（中間檔）
│   └── {SAMPLE_ID}.pangolin.vcf.gz.tbi
├── 03_acmg/
│   └── {SAMPLE_ID}.snv_indel.acmg.tsv             ← ★ SNV/Indel 最終輸出（65 欄）
├── 04_mito/                                        ← ★ v3.2 新增
│   ├── {SAMPLE_ID}.mito.tsv                       ← mtDNA 輸出（24 欄）
│   └── {SAMPLE_ID}.mito.vep.vcf.gz               ← VEP 中間檔
└── 05_str/                                         ← ★ v3.2 新增
    └── {SAMPLE_ID}.str.tsv                        ← STR 輸出（22 欄）
```

> **v3.2 新增：** `04_mito/` 和 `05_str/` 目錄。STR TSV 只包含有對應 STRchive 記錄的 locus（已知致病 STR），DRAGEN 約 53 筆，NCKUH WES 約 5 筆。

---

### 主要輸出欄位（snv_indel.acmg.tsv，65 欄）

#### 位置資訊（欄 1–5）
| 欄位 | 說明 |
|------|------|
| CHROM, POS, REF, ALT | 變異座標 |
| RS_ID | dbSNP rsID（如 rs72631890）|

#### Transcript 資訊（欄 6–15）
| 欄位 | 說明 |
|------|------|
| GENE | HGNC gene symbol |
| TRANSCRIPT | Ensembl transcript ID |
| TRANSCRIPT_TYPE | MANE_SELECT / MANE_PLUS_CLINICAL / CANONICAL |
| HGVS_C, HGVS_P | HGVS 命名 |
| CONSEQUENCE | VEP consequence（如 missense_variant）|
| IMPACT | HIGH / MODERATE / LOW / MODIFIER |
| EXON, INTRON | Exon/Intron 編號（格式：2/19）|
| MANE_ALL | JSON：所有 MANE transcript 的後果（供 GUI 展開）|

#### Caller 資訊（欄 16–24）
| 欄位 | 說明 |
|------|------|
| CALLERS | `DV+HC` / `DV` / `HC`（NCKUH）或 `DRAGEN` |
| DP_DV, AD_DV, VAF_DV | DeepVariant read depth、allelic depth、VAF |
| DP_HC, AD_HC | HaplotypeCaller read depth、allelic depth |
| ZYGOSITY | het / hom / hemizygous / unknown |
| GT_DV, GT_HC | Genotype |

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

#### ACMG 分類（欄 62–65）
| 欄位 | 說明 |
|------|------|
| ACMG_CRITERIA | 觸發的所有 criteria，逗號分隔（如 `PVS1,PM2_Supporting`）|
| ACMG_SCORE | 數值化分數（越高越可能致病，用於排序）|
| ACMG_CLASS | `Pathogenic` / `Likely_Pathogenic` / `VUS` / `Likely_Benign` / `Benign` |
| ACMG_NOTES | 觸發原因說明，含數值依據（如 `PVS1:LOFTEE=HC,gene=MECP2,HI=3`）|

---

### mito.tsv 欄位（24 欄）

| 欄位 | 說明 |
|------|------|
| CHROM, POS, REF, ALT | 位置 |
| GENE | MT gene（MT-ND1, MT-RNR2 等）|
| HGVS_C, HGVS_P | HGVS 命名（upstream variant 無值）|
| CONSEQUENCE, IMPACT, BIOTYPE | VEP 後果 |
| GENOTYPE, DP | 樣本 GT 和 read depth |
| AF_SAMPLE | Heteroplasmy level（FORMAT/AF，0-1）|
| GNOMAD_MITO_AF | gnomAD mito AF（目前全為 `.`，VEP cache 限制）|
| CLINVAR_SIG, CLINVAR_DN | ClinVar 致病性和疾病名稱 |
| MITOMAP_LOCUS | MITOMAP locus（如 MT-CYB）|
| MITOMAP_DISEASE | MITOMAP 疾病名稱 |
| MITOMAP_STATUS | Cfrm（確認）/ Reported / Conflicting reports |
| MITOMAP_HOMO, MITOMAP_HETERO | Homoplasmy/Heteroplasmy 報告（+/-/nr）|
| MITOMAP_MITOTIP | MitoTIP score（tRNA/rRNA only）|

> **注意：** NCKUH 只輸出 PASS variant；DRAGEN 輸出所有 chrM variant（含 non-PASS），請依 FILTER 欄位篩選。

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

## 驗證指令

收到新版輸出後，請執行以下指令確認正確：

```bash
# 設定路徑（替換 SAMPLE_ID）
SAMPLE_ID=26WE0001
TSV=/home/pipeline/tertiary_output/${SAMPLE_ID}/03_acmg/${SAMPLE_ID}.snv_indel.acmg.tsv

echo "========================================="
echo "Step 1：欄位數（應為 65）"
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
echo "Step 4：CALLERS 欄位確認"
echo "（NCKUH 應為 DV+HC/DV/HC；DRAGEN 應為 DRAGEN）"
echo "========================================="
awk -F'\t' 'NR>1 {print $16}' $TSV | sort | uniq -c | sort -rn | head -5

echo ""
echo "========================================="
echo "Step 5：ClinVar 分布（應有非 . 的值）"
echo "========================================="
awk -F'\t' 'NR>1 {print $32}' $TSV | sort | uniq -c | sort -rn | head -10

echo ""
echo "========================================="
echo "Step 6：ACMG 分類分布"
echo "========================================="
awk -F'\t' 'NR>1 {print $64}' $TSV | sort | uniq -c | sort -rn

echo ""
echo "========================================="
echo "Step 7：P/LP variant 詳細資訊"
echo "========================================="
awk -F'\t' 'NR>1 && ($64=="Pathogenic" || $64=="Likely_Pathogenic")' $TSV \
    | cut -f1,2,6,11,62,63,64,65 | head -10

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

# MITOMAP 命中筆數
awk -F'\t' 'NR>1 && $21!="."'  $MITO | wc -l

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

**Q：mito.tsv 的 GNOMAD_MITO_AF 全部是 `.`？**

這是已知限制。VEP 115 cache 的 `--af_gnomadg` 不包含 chrM 的 gnomAD mito 資料，所以這個欄位目前無法填入。臨床判讀請參考 MITOMAP_STATUS 和 AF_SAMPLE（heteroplasmy level）。

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
