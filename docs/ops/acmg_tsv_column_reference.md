# 三級分析 ACMG TSV 欄位說明文件

**檔案名稱：** `{SAMPLE_ID}.snv_indel.acmg.tsv`  
**產生程式：** `acmg_classifier.py`（輸入為 `parse_vep_csq.py` 產生的 annotated TSV）  
**總欄位數：** 65 欄（61 個 annotation 欄位 + 4 個 ACMG 欄位）  
**缺失值表示：** `.`（一律用英文句點）  
**分隔符號：** Tab（`\t`）  
**版本：** v3.1（2026-05）

---

## 一、位置資訊（欄 1–5）

| # | 欄位名稱 | 資料型別 | 範例值 | 說明 |
|---|---------|---------|-------|------|
| 1 | `CHROM` | String | `chr1` | 染色體，UCSC 格式（chr 前綴） |
| 2 | `POS` | Integer | `94730790` | 位置，1-based，VCF 標準 |
| 3 | `REF` | String | `C` | Reference allele |
| 4 | `ALT` | String | `T` | Alternate allele（multi-allelic 用逗號分隔，少見） |
| 5 | `RS_ID` | String | `rs1801133` | dbSNP rsID；優先用 ClinVar lookup，備用 VEP Existing_variation；無則 `.` |

---

## 二、Transcript 資訊（欄 6–15）

| # | 欄位名稱 | 資料型別 | 範例值 | 說明 |
|---|---------|---------|-------|------|
| 6 | `GENE` | String | `MTHFR` | HGNC 官方基因符號（VEP `--symbol`） |
| 7 | `TRANSCRIPT` | String | `ENST00000376592.5` | 代表 transcript 的 Ensembl ID（含版本號） |
| 8 | `TRANSCRIPT_TYPE` | String | `MANE_SELECT` | 代表 transcript 類型：`MANE_SELECT` / `MANE_PLUS_CLINICAL` / `CANONICAL` / `OTHER` |
| 9 | `HGVS_C` | String | `ENST00000376592.5:c.665C>T` | cDNA 命名（HGVS nomenclature） |
| 10 | `HGVS_P` | String | `ENSP00000365783.2:p.Ala222Val` | 蛋白質命名；同義突變或非編碼為 `.` |
| 11 | `CONSEQUENCE` | String | `missense_variant` | VEP SO consequence；多個 consequence 用 `&` 連接，例如 `missense_variant&splice_region_variant` |
| 12 | `IMPACT` | String | `MODERATE` | VEP impact 等級：`HIGH` / `MODERATE` / `LOW` / `MODIFIER` |
| 13 | `EXON` | String | `5/11` | 所在 exon（`exon編號/總exon數`）；intronic variant 為 `.` |
| 14 | `INTRON` | String | `.` | 所在 intron（`intron編號/總intron數`）；exonic variant 為 `.` |
| 15 | `MANE_ALL` | JSON String | `[{"tx":"NM_005957.5","enst":"ENST00000376592","type":"MANE_SELECT","consequence":"missense_variant","hgvsc":"...","hgvsp":"...","impact":"MODERATE"}]` | 所有 MANE Select / MANE Plus Clinical transcript 的 JSON 陣列；GUI 可展開顯示所有 MANE transcript；無 MANE 時為 `[]` |

> **MANE_ALL JSON 欄位說明：**
> 每個元素包含：`tx`（RefSeq NM ID）、`enst`（Ensembl ID）、`type`（MANE_SELECT 或 MANE_PLUS_CLINICAL）、`consequence`、`hgvsc`、`hgvsp`、`impact`

---

## 三、Caller 資訊（欄 16–24）

| # | 欄位名稱 | 資料型別 | 範例值 | 說明 |
|---|---------|---------|-------|------|
| 16 | `CALLERS` | String | `DV+HC` | 偵測到此 variant 的 caller：`DV+HC`（兩者都有，最高信心）/ `DV`（只有 DeepVariant）/ `HC`（只有 HaplotypeCaller） |
| 17 | `DP_DV` | Integer | `45` | DeepVariant read depth（FORMAT/DP） |
| 18 | `AD_DV` | String | `28,17` | DeepVariant allelic depth，格式 `REF,ALT`（FORMAT/AD） |
| 19 | `VAF_DV` | Float | `0.3778` | DeepVariant variant allele fraction（FORMAT/VAF），四捨五入至小數點後 4 位 |
| 20 | `DP_HC` | Integer | `43` | HaplotypeCaller read depth |
| 21 | `AD_HC` | String | `27,16` | HaplotypeCaller allelic depth，格式 `REF,ALT` |
| 22 | `ZYGOSITY` | String | `het` | 接合性：`het`（雜合）/ `hom`（純合）/ `hemizygous`（半合，chrX/Y）/ `ref`（參考序列）/ `unknown` |
| 23 | `GT_DV` | String | `0/1` | DeepVariant genotype |
| 24 | `GT_HC` | String | `0/1` | HaplotypeCaller genotype |

---

## 四、族群頻率（欄 25–31）

所有頻率欄位為 Float（0–1），資料庫中未收錄的 variant 為 `.`。

| # | 欄位名稱 | 來源 | 範例值 | 說明 |
|---|---------|------|-------|------|
| 25 | `GNOMAD_G_AF` | gnomAD genome v4 | `0.000123` | 全人群 allele frequency |
| 26 | `GNOMAD_G_EAS_AF` | gnomAD genome v4 | `0.000089` | 東亞族群 AF |
| 27 | `GNOMAD_E_AF` | gnomAD exome v4 | `0.000145` | 全人群 AF |
| 28 | `GNOMAD_E_EAS_AF` | gnomAD exome v4 | `0.000067` | 東亞族群 AF |
| 29 | `GNOMAD_E_AF_DBNSFP` | dbNSFP 4.9c（gnomAD exome） | `0.000140` | 全人群 AF（dbNSFP 版本，補充用） |
| 30 | `GNOMAD_E_EAS_AF_DBNSFP` | dbNSFP 4.9c（gnomAD exome） | `0.000060` | 東亞族群 AF（dbNSFP 版本） |
| 31 | `TG_EAS_AF` | 1000 Genomes EAS | `0.0` | 東亞族群 AF（1000 Genomes Phase 3） |

> **GUI 建議：** 顯示時優先用 `GNOMAD_G_AF`，其次 `GNOMAD_E_AF`；若要顯示 EAS 頻率則用 `GNOMAD_G_EAS_AF`。

---

## 五、ClinVar（欄 32–36）

| # | 欄位名稱 | 資料型別 | 範例值 | 說明 |
|---|---------|---------|-------|------|
| 32 | `CLINVAR_SIG` | String | `Pathogenic` | ClinVar 臨床意義；原始值，例如 `Pathogenic` / `Likely_pathogenic` / `Benign` / `Uncertain_significance` / `Conflicting_classifications_of_pathogenicity` |
| 33 | `CLINVAR_STARS` | Integer | `2` | ClinVar review status 星星數（0–4）；`4`=practice guideline，`3`=expert panel，`2`=多機構一致，`1`=單機構或衝突，`0`=無標準 |
| 34 | `CLINVAR_DN` | String | `Homocystinuria,methylcobalamin...` | ClinVar 疾病名稱（`&_` 分隔多個疾病） |
| 35 | `CLINVAR_SIGCONF` | String | `Pathogenic(2)&Uncertain_significance(1)` | 有衝突分類時的詳細記錄；無衝突則為 `.` |
| 36 | `CLINVAR_VARIATION_ID` | Integer | `1265` | ClinVar Variation ID；GUI 組連結用：`https://www.ncbi.nlm.nih.gov/clinvar/variation/{ID}/` |

---

## 六、OMIM（欄 37）

| # | 欄位名稱 | 資料型別 | 範例值 | 說明 |
|---|---------|---------|-------|------|
| 37 | `OMIM_IDS` | String | `236250,609325` | OMIM phenotype ID，逗號分隔（一個 variant 可對應多個疾病）；GUI 組連結用：`https://www.omim.org/entry/{ID}` |

---

## 七、LOFTEE（欄 38–41）

| # | 欄位名稱 | 資料型別 | 範例值 | 說明 |
|---|---------|---------|-------|------|
| 38 | `LOFTEE` | String | `HC` | LOFTEE LoF 分級：`HC`（High Confidence）/ `LC`（Low Confidence）/ `.`（非 LoF variant） |
| 39 | `LOFTEE_FILTER` | String | `.` | LOFTEE filter 原因；`.` 代表通過所有 filter（HC 且無問題）；有值代表 LC 的原因，例如 `END_TRUNC`、`INCOMPLETE_CDS` |
| 40 | `LOFTEE_FLAGS` | String | `.` | LOFTEE 額外旗標，例如 `PHYLOCSF_WEAK`；`.` 代表無 |
| 41 | `LOFTOOL` | Float | `0.0841` | LoFtool 基因 LoF 不耐受分數（0–1，越小越不耐受）；**基因層級分數**，同一基因所有 variant 相同 |

---

## 八、In Silico 預測分數（欄 42–54）

所有分數欄位為 Float，預測欄位為 String（`D`=Deleterious/Damaging、`T`=Tolerated/Benign），缺失為 `.`。

| # | 欄位名稱 | 範圍 | 範例值 | 說明 |
|---|---------|------|-------|------|
| 42 | `BAYESDEL_NOAF` | −1 ~ +1 | `0.285` | BayesDel（不含 AF 特徵）；> 0.5 = pathogenic（ClinGen 校準閾值另見 ACMG_NOTES） |
| 43 | `BAYESDEL_NOAF_PRED` | D/T | `D` | BayesDel 預測：`D`（Deleterious）/ `T`（Tolerated） |
| 44 | `ALPHAMISSENSE` | 0 ~ 1 | `0.923` | AlphaMissense；> 0.564 = likely pathogenic（開發者閾值）；ClinGen 校準閾值見 ACMG_NOTES |
| 45 | `ALPHAMISSENSE_PRED` | String | `likely_pathogenic` | AlphaMissense 預測：`likely_pathogenic` / `ambiguous` / `likely_benign` |
| 46 | `ESM1B` | 約 −30 ~ +10 | `−11.3` | ESM1b 蛋白質語言模型分數；**越負越致病**（方向與其他工具相反） |
| 47 | `ESM1B_PRED` | D/T | `D` | ESM1b 預測 |
| 48 | `VARITY_R` | 0 ~ 1 | `0.721` | VARITY_R（rare variant 訓練版）；> 0.5 = pathogenic（開發者閾值） |
| 49 | `SIFT` | 0 ~ 1 | `0.02` | SIFT 分數；**越小越致病**（< 0.05 = deleterious） |
| 50 | `SIFT_PRED` | D/T | `D` | SIFT 預測：`deleterious` / `tolerated` |
| 51 | `DANN` | 0 ~ 1 | `0.9821` | DANN deep learning 分數；越大越致病 |
| 52 | `PHACTBOOST` | 0 ~ 1 | `0.634` | PHACTboost；越大越致病 |
| 53 | `PHYLOP100` | 約 −20 ~ +10 | `5.12` | phyloP100way vertebrate 保守性分數；**> 3.58 = 高度保守**（BP7/BP3 判斷閾值） |
| 54 | `GERP` | 約 −12 ~ +6 | `4.87` | GERP++ RS 分數；> 2 = 保守性位置 |

---

## 九、P-KNN（欄 55–56）

| # | 欄位名稱 | 資料型別 | 範例值 | 說明 |
|---|---------|---------|-------|------|
| 55 | `PKNN_LLR` | Float | `2.847` | P-KNN joint calibration log likelihood ratio；正值支持致病，負值支持良性；只有 missense variant 有值（~30% 命中率） |
| 56 | `PKNN_EVIDENCE` | String | `PP3_Moderate` | 由 `PKNN_LLR` 換算的 ACMG evidence 強度：`PP3_Supporting/Moderate/Strong` / `BP4_Supporting/Moderate/Strong` / `.` |

---

## 十、Splice（欄 57–58）

| # | 欄位名稱 | 資料型別 | 範例值 | 說明 |
|---|---------|---------|-------|------|
| 57 | `PANGOLIN_SCORE` | Float | `0.832` | Pangolin splice effect 最大分數（絕對值）；≥ 0.1 表示可能影響 splicing；只有 splice candidate 有值 |
| 58 | `PANGOLIN_DETAIL` | String | `BRCA1\|ACCEPTOR_GAIN:0.832:123456789\|DONOR_LOSS:...` | Pangolin 詳細輸出，包含 gain/loss 類型和位置 |

---

## 十一、蛋白質資訊（欄 59–60）

| # | 欄位名稱 | 資料型別 | 範例值 | 說明 |
|---|---------|---------|-------|------|
| 59 | `DOMAINS` | String | `Pfam_domain:PF00001&SMART_domains:SM00320` | VEP DOMAINS，`&` 分隔多個 domain；可用於 PM1（functional domain）和 BP3（repeat region）判斷 |
| 60 | `SWISSPROT` | String | `P35557` | UniProt/SwissProt accession；GUI 組連結用：`https://www.uniprot.org/uniprot/{ID}` |

---

## 十二、基因識別碼（欄 61）

| # | 欄位名稱 | 資料型別 | 範例值 | 說明 |
|---|---------|---------|-------|------|
| 61 | `HGNC_ID` | String | `HGNC:7436` | HGNC 基因識別碼（含 `HGNC:` 前綴）；GUI 組連結用：`https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{ID}` |

---

## 十三、ACMG 分類結果（欄 62–65）

這四欄由 `acmg_classifier.py` 計算產生，是 GUI 最主要顯示的欄位。

| # | 欄位名稱 | 資料型別 | 範例值 | 說明 |
|---|---------|---------|-------|------|
| 62 | `ACMG_CRITERIA` | String | `PVS1,PM2_Supporting` | 觸發的所有 criteria，逗號分隔；`.` 表示無 criteria 觸發 |
| 63 | `ACMG_SCORE` | Integer | `9` | 數值化總分（Tavtigian 2020 點數系統）；可用於排序：分數越高越可能致病；分數越低越傾向良性 |
| 64 | `ACMG_CLASS` | String | `Likely_Pathogenic` | 最終 ACMG 分類，共 5 類（見下表） |
| 65 | `ACMG_NOTES` | String | `PVS1:LOFTEE=HC,gene=SMC1A,HI=3\|PM2_Supporting:gene=SMC1A,MOI=Unknown,gnomAD_AF=absent` | 觸發原因說明，`\|` 分隔各 criteria；包含觸發的數值依據，供臨床人員審閱 |

### ACMG_CLASS 分類對應

| `ACMG_CLASS` 值 | ACMG_SCORE 範圍 | 說明 |
|----------------|----------------|------|
| `Pathogenic` | ≥ 10 | 致病 |
| `Likely_Pathogenic` | 6 ~ 9 | 可能致病 |
| `VUS` | 0 ~ 5 | Variant of Uncertain Significance |
| `Likely_Benign` | −1 ~ −6 | 可能良性 |
| `Benign` | ≤ −7 或觸發 BA1 | 良性（觸發 BA1 直接判為 Benign，不看分數） |

### ACMG_CRITERIA 常見值說明

| Criteria | 點數 | 觸發條件 |
|---------|------|---------|
| `BA1` | −8 | gnomAD AF > 5%（直接 Benign） |
| `BS1` | −4 | gnomAD AF > 1% |
| `PVS1` | +8 | LOFTEE HC + ClinGen HI score = 3 |
| `PM4` | +2 | inframe indel 或 stop_lost |
| `PM2_Supporting` | +1 | EAS AF 極罕見（< 0.001%）或 AD 基因全人群 AF < 0.01% |
| `PP3_Supporting` | +1 | 計算工具預測致病（Supporting 強度） |
| `PP3_Moderate` | +2 | 計算工具預測致病（Moderate 強度） |
| `PP3_P3` | +3 | 計算工具預測致病（+3 強度） |
| `PP3_Strong` | +4 | 計算工具預測致病（Strong 強度） |
| `BP3` | −1 | inframe indel 在 repeat region + 低保守性 |
| `BP4_Supporting` | −1 | 計算工具預測良性（Supporting 強度） |
| `BP4_Moderate` | −2 | 計算工具預測良性（Moderate 強度） |
| `BP4_M3` | −3 | 計算工具預測良性（−3 強度） |
| `BP4_Strong` | −4 | 計算工具預測良性（Strong 強度） |
| `BP7` | −1 | synonymous + 低 splice 影響 + 低保守性 |

> **PP3/BP4 計算工具優先順序（cascade fallback）：**
> P-KNN → AlphaMissense → ESM1b → VARITY_R → BayesDel_noAF
> 只用第一個有值且不在 indeterminate 範圍的工具；觸發的工具名稱記錄在 `ACMG_NOTES`。

> **注意：** `ACMG_NOTES` 有時會包含 `PVS1_NOT_TRIGGERED:...` 字樣，這代表 variant 是 LOFTEE HC 但因為 ClinGen HI score 不夠（例如 HI=2 或不在 ClinGen 資料庫），PVS1 未觸發，需要人工審閱。

---

## 十四、Annotation database version sidecar

ClinVar release date 不應重複寫在 TSV 每一列，也不能由 TSV mtime 推測。三級 pipeline 在完成 `03_acmg/{source}.snv_indel.acmg.tsv` 時，應一起原子寫入：

`03_acmg/{source}.annotation_versions.json`

```json
{
  "schema_version": 1,
  "pipeline_version": "3.5.0",
  "databases": {
    "clinvar": {
      "release_date": "2026-05-10",
      "source": "clinvar_20260510.vcf.gz",
      "sha256": "<optional source checksum>"
    }
  }
}
```

`release_date` 必須是 annotation run 實際使用的 ClinVar snapshot 日期（ISO `YYYY-MM-DD`），不是執行日期或下載日期。NGS-UI 也相容 `{tsv filename}.meta.json`、`03_acmg/annotation_versions.json`、`08_postprocessing/annotation_versions.json` 與 `pipeline_source.json.annotation_versions`，但新 pipeline 應優先使用上述 per-source sidecar。缺少 metadata 時 UI 顯示「三級輸出未提供」，不會以 hard-coded 日期代替。

---

## 附錄：GUI 建議的顯示優先順序

1. **主列表（需顯示）：** `GENE`, `HGVS_C`, `HGVS_P`, `CONSEQUENCE`, `ZYGOSITY`, `GNOMAD_G_EAS_AF`, `CLINVAR_SIG`, `CLINVAR_STARS`, `ACMG_CRITERIA`, `ACMG_CLASS`, `ACMG_SCORE`
2. **展開詳情（點擊後顯示）：** 族群頻率全部欄位、所有 in silico 分數、`MANE_ALL` JSON、`ACMG_NOTES`
3. **連結組成：**
   - ClinVar：`https://www.ncbi.nlm.nih.gov/clinvar/variation/{CLINVAR_VARIATION_ID}/`
   - OMIM：`https://www.omim.org/entry/{每個OMIM_ID}`
   - UniProt：`https://www.uniprot.org/uniprot/{SWISSPROT}`
   - HGNC：`https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{HGNC_ID}`
4. **預設排序：** `ACMG_SCORE` 降冪（最可能致病的排在最上面）
5. **預設過濾：** 可考慮預設隱藏 `ACMG_CLASS = Benign`，讓臨床人員專注在 P/LP/VUS/LB
