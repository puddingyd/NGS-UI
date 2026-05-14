# 三級分析 Pipeline 反饋 — `puddingyd/tertiary_code`

> 對接 GUI 端 (`puddingyd/NGS-UI`) 過程中對新版三級 pipeline 的整理觀察。  
> 驗證對象：**VAL-58**（WGS）  
> Secondary 輸入：`/home/datalake_Intermediate/n102968/nextflow_output/WGS_VAL/VAL-58/`  
> Tertiary 輸出：`tertiary_output/VAL-58/VAL-58.snv_indel.annotated.tsv`  
> 截至：2026-05-14

排序：**修起來最會痛的先講**。每條都附 (a) 觀察 (b) 影響 (c) 建議。中段大半是「我這邊已先寫 stop-gap script 補上、等你補完可拔除」的部分。

---

## 🔴 P0 — 資料正確性，必修

### 1. `VAL-58.ensemble.fixed.vcf.gz` 的 HC sample column 名稱錯誤

**觀察**

`#CHROM` 行只有兩個 sample column：

```
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  VAL-58_DV  VAL-26
```

第二個應該是 `VAL-58_HC`（HaplotypeCaller 那條沒被 reheader）。`scripts/add_callers_tag.py` 的 docstring 明寫期待 `{sample_id}_DV` / `{sample_id}_HC`。

**影響**

tertiary pipeline 第一步 `PREPARE_VCF:ADD_CALLERS_TAG` 直接 fail：

```
[ERROR] 找不到 HC sample column：VAL-58_HC
```

只能手動 `bcftools reheader -s rename.tsv` 後重 stage 一份 input dir 才能跑。

**建議**

Secondary pipeline ensemble 合併那一步，HC VCF 在 merge 前先 `bcftools reheader -s` 成 `${sample}_HC`（DV 那邊已經做對了）。順便驗 sample-sheet → BAM SM 對應，`VAL-26` 看起來像某條 sample 寫錯混進來。

---

### 2. Apptainer / Nextflow temp 目錄非建置者寫不進去

**觀察**

```
/home/pipeline/nextflow_temp/        drwxr-sr-x  n101569    其他人寫不進
/home/pipeline/apptainer_temp/       drwxr-sr-x  n101569    同上
/home/pipeline/nextflow_containers/  drwxr-sr-x  n101569    APPTAINER_CACHEDIR 也指這
```

`NGS2ndAnalysis_env.sh` 把 `NXF_TEMP` / `APPTAINER_TMPDIR` / `APPTAINER_CACHEDIR` 三個 env var 都 export 到上述位置。

**影響**

任何不是原 owner 的人 source env 後跑 nextflow / apptainer，第一步必炸：

```
java.nio.file.AccessDeniedException: /home/pipeline/nextflow_temp/nxf-...
FATAL: mkdir /home/pipeline/apptainer_temp/build-temp-...: permission denied
```

每個新使用者第一次都要自己摸清楚要 override 哪三個 env var。

**建議**

擇一：

- env script 改成 `NXF_TEMP=$HOME/.nextflow_tmp`、`APPTAINER_TMPDIR=$HOME/.apptainer/tmp` 之類 user-local default
- 那三個目錄 chmod `drwxrwxrwt`（group write + sticky）並 chgrp 到 shared group
- 至少在 env script 開頭 print「請先 export 以下變數到自己家目錄」

---

## 🟠 P1 — GUI 整合卡點 / 我們被迫做 stop-gap

### 3. Phase-1 TSV 完全未過濾，5.35M 列 / 1.6 GB

**觀察**

`{sample}.snv_indel.annotated.tsv` 把 VEP-PICK 完所有變異 dump 出來，包含 `IMPACT=MODIFIER` 的 `upstream_gene_variant`、`intron_variant` 等臨床完全不報的位點。

```
$ wc -l VAL-58.snv_indel.annotated.tsv
5350698  → 1.6 GB

$ awk -F'\t' '...' IMPACT 分佈
5319084 MODIFIER
  19650 LOW
  11355 MODERATE
    608 HIGH
```

**影響**

- GUI adapter 把整份 TSV 讀進 dict → OOM / service 卡死（我這邊撞到過一次）
- 已寫 stop-gap `scripts/filter_snv_tsv.py`：rare (`GNOMAD_G_AF ≤ 0.01`) + `IMPACT ∈ {HIGH, MODERATE}` + 排除 alt-contig + 排除 `*`-allele + 保留 ClinVar P/LP rescue → 5.35M 列 → ~1,363 列

**建議**

在 `PARSE_VEP_CSQ` 之後加 `FILTER` step，輸出兩份：

| 檔名 | 目的 |
|---|---|
| `{sample}.snv_indel.full.annotated.tsv` | 保留現況、archive / re-analysis 用 |
| `{sample}.snv_indel.annotated.tsv` | 過濾完，給 GUI 用 |

過濾規則參考 `scripts/filter_snv_tsv.py` 的 default（可以再放寬到也保留 `LOW`），但**至少要砍掉 MODIFIER + alt contig + `*`-allele**。

---

### 4. CLINVAR_* 四欄全是 placeholder

**觀察**

VEP 設定有 reference 到 ClinVar VCF (`clinvar_20260418.vcf.gz`)，但 TSV 裡四欄都沒值：

```
CLINVAR_SIG       = .   (所有列)
CLINVAR_STARS     = 0   (所有列)
CLINVAR_DN        = .
CLINVAR_SIGCONF   = .
```

VAL-58 整份 TSV `CLINVAR_SIG ~ /Pathogenic|Likely_pathogenic/` 命中數 = **0**。

**影響**

- 所有 tier 分類錯位：`1A` (ClinVar P/LP ≥ 1★) 永遠空、`1B` (Frameshift LOFTEE HC) 雖有但少
- Reviewer 看不到 ClinVar 已知致病資訊
- 已寫 stop-gap `scripts/annotate_clinvar.py`：直接讀本地 `clinvar_20260418.vcf.gz`，by `(CHROM,POS,REF,ALT)` join，backfill 四欄。VAL-58 跑下來 backfill 48,514 列

**建議**

VEP 加：

```bash
--custom file=clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNSIGCONF
```

或在 PARSE_CSQ 後用 bcftools annotate 補。確保產出：

| 欄位 | 來源 / 對應 |
|---|---|
| `CLINVAR_SIG` | `INFO/CLNSIG` (Pathogenic / Likely_pathogenic / ...) |
| `CLINVAR_STARS` | `INFO/CLNREVSTAT` → NCBI 0-4 階層（practice_guideline=4, reviewed_by_expert_panel=3, multiple_submitters_no_conflicts=2, single_submitter / conflicting=1, no_assertion=0） |
| `CLINVAR_DN` | `INFO/CLNDN`，`%2C` URL-decode |
| `CLINVAR_SIGCONF` | `INFO/CLNSIGCONF`（Conflicting 才有） |

**Sentinel 一致性 bug**：目前 `CLINVAR_STARS` 預設 `0` 但 `CLINVAR_SIG` 預設 `.` → 混合 sentinel。GUI 端「該欄是否有資料」邏輯（`!_has_value(row, ...)`）會被 `"0"` 誤判為「已有資料」，導致 backfill 漏寫 STARS。建議「沒做 ClinVar 註解時四欄統一 `.` / 統一空白」。

---

### 5. ACMG_CLASSIFY 整步未實作

**觀察**

```
ACMG_EVIDENCE = (empty)
ACMG_POINTS   = (empty)
ACMG_CLASS    = (empty)
```

`main_tertiary.nf` 進度註解：`🔲 ACMG_CLASSIFY - ACMG evidence 收集與分類（Phase 1）`。

**影響**

Tier `1C`（ACMG ≥ 4 points）規則完全不生效；reviewer 看不到 baseline ACMG class。

我用 GeneBe REST API 做 stop-gap (`scripts/annotate_acmg_genebe.py`)：rare variants → sites VCF → GeneBe annotate → 取 `INFO/acmg_score` 跟 `INFO/acmg_criteria` 回填三欄。VAL-58 跑下來填了 1,357 / 1,363 列。

**建議**

依你原計畫 Phase 1 上 ACMG_CLASSIFY 即可。GUI adapter 是「pipeline 自己有 ACMG 值 → 直接用；沒值 → 不接 GeneBe」的設計（其實 GeneBe 是回填到 TSV，pipeline 補上後 stop-gap script 自然停用）。

---

### 6. CNV/SV 同事件被多 caller 重複報告

**觀察**

VAL-58 的 SV tier 上看到同一個 DEL/DUP 被 Manta / DELLY / cnvkit 等 caller 各報一次，breakpoint 差 100-500 bp。AnnotSV 對每筆獨立算 `ranking_score`，但它們其實是同一事件。

**影響**

- UI 上同位點看到 4-5 張幾乎一樣的 SV 卡，reviewer 工作量大
- AnnotSV ranking 被汙染（同事件分數加成）
- 前端已做視覺 collapse（同 `CHROM` + 同 `sv_type` + reciprocal overlap ≥ 0.8 合一張，alias 收進 expander），但只是視覺層

**建議**

`05_cnv_sv` 多 caller union 之後、進 AnnotSV **之前**跑 SV merge。兩個業界標準工具：

**Truvari**（推薦）：

```bash
bcftools sort sv.union.vcf.gz | bgzip > sv.sorted.vcf.gz
tabix -p vcf sv.sorted.vcf.gz
truvari collapse \
  -i sv.sorted.vcf.gz \
  -o sv.merged.vcf.gz -c sv.removed.vcf.gz \
  --refdist 500 \
  --pctseq 0 \
  --pctsize 0.7 \
  --pctovl 0.8 \
  --keep maxqual
```

**SURVIVOR**：

```bash
SURVIVOR merge vcf_list.txt 500 1 1 1 0 30 sv.merged.vcf
# 500bp breakpoint tolerance, ≥1 caller, type match, strand match,
# ignore length, min 30bp
```

合完 SUPP_VEC 還能保留「哪些 caller 支持」資訊，下游 ranking 才不會被同事件加成。pipeline 修完後我這邊前端 cluster 邏輯會自然 no-op。

---

### 7. 欄名漂移（舊 → 新 pipeline 未對齊）

新版 TSV 跟舊 pipeline 的 contract 有 7 個 rename + 多個拆分 + 多個遺失。GUI adapter 我用 `_coalesce` / `_max_multi` 都接好了，但 long-term 應該對齊：

| 邏輯欄位 | 舊 | 新 | 我的 stop-gap |
|---|---|---|---|
| AD（per caller） | `AD` | `AD_DV` / `AD_HC` | `_coalesce(AD, AD_DV, AD_HC)` |
| VAF | `VAF` | `VAF_DV` / `VAF_HC` | `_coalesce(VAF, VAF_DV, VAF_HC)` |
| ClinVar conflicting | `CLINVAR_CONF` | `CLINVAR_SIGCONF` | `_coalesce` |
| LOFTEE HC 判定 | `LOFTEE_HC`（HC / blank） | `LOFTEE`（HC / LC / .） | `_coalesce` |
| BayesDel | `BAYESDEL` | `BAYESDEL_NOAF` | `_coalesce`（注意 with-AF vs no-AF 是不同變體，數值不完全等價） |
| Splice impact | `SPLICEAI_MAX` | `PANGOLIN_SCORE` | 兩個 payload field 並存（Pangolin 還有正負號，代表 gain/loss） |
| LM path-pred | `ESM2_SCORE` | `ESM1B` | `_coalesce`（**ESM1b 比 ESM2 早一代，請確認是不是 typo / regression**） |

**建議**：對齊到舊欄名（讓既有 consumer 不用改），或正式 document 新 schema + bump 一個 version field 讓 consumer 顯式對齊。

---

## 🟡 P2 — 品質 / 顯示

### 8. VEP `--pick` 策略可能選錯 transcript

**觀察**

VAL-58 raw TSV 的 LoF consequence 統計（前幾名）：

```
391 splice_donor_region_variant&intron_variant
258 splice_donor_region_variant&intron_variant&non_coding_transcript_variant
146 splice_donor_variant&non_coding_transcript_variant
134 frameshift_variant
104 splice_donor_5th_base_variant&intron_variant
103 splice_acceptor_variant&non_coding_transcript_variant
 87 splice_donor_5th_base_variant&intron_variant&non_coding_transcript_variant
 84 stop_gained
```

`&non_coding_transcript_variant` suffix 出現多次 — PICK 選到 lncRNA / pseudogene transcript，可能蓋掉同位置 protein-coding 的 HIGH consequence。

**建議**

把 `--pick_order` 改成 protein-coding biotype 優先：

```
--pick_order biotype_rank,canonical,mane_select,tsl,appris,length
```

或更狠：完全不 PICK，輸出 per-variant × per-transcript 全展開，下游 PARSE_CSQ 自己挑 — 檔案大 5-10 倍但保留完整訊息。

---

### 9. LOFTEE 覆蓋率偏低

**觀察**

VAL-58 全 5.35M 列裡 LOFTEE 欄位分佈：

```
1299  .    （unannotatable 或 non-LoF）
  56  HC
   8  LC
```

對照 IMPACT=HIGH 的 ~80 列（frameshift / stop_gained / splice_donor / splice_acceptor），LOFTEE 只 cover 56。

**疑問**

- 是不是只有 PICK 選到的 transcript 才跑 LOFTEE？多 transcript 變異漏判？
- `--plugin LoF,...` 的 GERP / phyloP 路徑、conservation cutoffs 設定對嗎？

**建議**

確認 LOFTEE plugin 對所有 protein-coding transcript 都跑得到；如果 PICK 之後才跑，多 transcript 的 LoF 會漏。

---

### 10. 缺失欄位（舊 pipeline 有、新版沒）

按臨床用處排：

| 欄位 | 用處 | 重要度 |
|---|---|---|
| `TWB_AF` (Taiwan BioBank AF) | 台灣個案 BS1 / BA1 不能光用 gnomAD，TWB 必要 | **高** |
| `IN_ROH` | ROH 區內變異 → AR 候選關鍵旗標 | **中** |
| `IN_BLACKLIST` | QC blacklist 區域 (false-positive hotspot) | **中** |
| `REVEL` / `METARNN` / `CADD_PHRED` | 跟新版 AlphaMissense / VARITY_R 部分重疊但各 cohort calibration 不同 | 中 |
| `EVO2_SCORE` | LM-based path 預測；ESM1B 補上後可選 | 低 |
| `PKNN_LLR` | 之前計畫書寫的「primary missense PP3/BP4」metric | 看判讀策略 |
| `PHASE_GROUP` / `PHASE_RESULT` | compound-het 在 cis/trans 判定 | 高（但屬 Phase 4） |
| `OMIM_LINK` / `GNOMAD_LINK` / `CLINVAR_LINK` | GUI 有 fallback 拼網址 | **不必補** |
| `REPORT_CLASS` | legacy | **不必補** |

---

### 11. 新增欄位 GUI 還沒完全消化

新 pipeline 多出來的欄位（GUI 端我大多已接到 SNV 卡片 in-silico 列，但有些更積極呈現可期待）：

`IMPACT`（filter script 用、payload 沒帶到 UI）、`DP_DV` / `DP_HC`、`GNOMAD_E_AF_DBNSFP` / `_EAS_AF_DBNSFP`、`LOFTOOL`（gene-level intolerance）、`BAYESDEL_NOAF_PRED` / `ALPHAMISSENSE_PRED` / `ESM1B_PRED` / `SIFT_PRED`（categorical preds）、`VARITY_R` / `SIFT` / `DANN` / `PHACTBOOST` / `PHYLOP100` / `GERP`、`PANGOLIN_DETAIL`、`DOMAINS` / `SWISSPROT`。

這部分純粹 GUI 端的呈現，pipeline 不用改。

---

## 🟢 P3 — 已知 Phase 2-4 待補（你計畫書已列）

`main_tertiary.nf` 的 TODO 跟我的觀察一致：

- 🔲 **CNV/SV (AnnotSV)** — Phase 2
- 🔲 **Mitochondria** — 我這邊用 `scripts/annotate_mito_vcf.sh` (MITOMAP-only) 補上當 stop-gap，Phase 3 上線可直接 replace
- 🔲 **STR (STRchive)** — Phase 3
- 🔲 **ROH** — 跟欄位 `IN_ROH` 連動；目前 GUI 的 ROH summary card 顯示「(無資料)」
- 🔲 **Phenotype (Exomiser + LIRICAL)** — Phase 2；UI 上 reviewer 按 ▶ 開始分析 會 trigger，但需要 pipeline 把 yml input 準備好
- 🔲 **WhatsHap phasing** — 跟 `PHASE_*` 欄位連動，Phase 4
- 🔲 **PGx (Aldy)** — Phase 3；UI 上 PharmCAT 卡空著
- 🔲 **ACMG SF v3.2** — Phase 3；secondary findings 區段空著

---

## 📎 附錄 — 我這邊的三隻 stop-gap script

GUI repo `puddingyd/NGS-UI` 的 `scripts/`：

| Script | 對應 pipeline 缺失 | 何時可以拔 |
|---|---|---|
| `annotate_clinvar.py` | §4 (ClinVar 註解) | pipeline 補上 VEP `--custom` 或 bcftools annotate ClinVar 後 |
| `filter_snv_tsv.py` | §3 (TSV 未過濾) | pipeline FILTER step 上線後 |
| `annotate_acmg_genebe.py` | §5 (ACMG_CLASSIFY) | pipeline ACMG_CLASSIFY step 上線後 |

三隻邏輯都是「**只填空白格**」— pipeline 真的補上、TSV 該欄非空時，script 變成 no-op。可以直接吸進 pipeline（用 nextflow `process` 包），參數 / 邏輯 / 預設值都有清楚的 docstring。

---

## 💬 建議的優先序

1. **馬上**：P0 那兩條（VAL-58 sample naming + 共用 temp dir 權限）— 不修任何人都跑不動
2. **下一個 sprint**：P1 §3 / §4 / §5（unfiltered TSV / ClinVar / ACMG）— GUI 立刻能呈現完整資訊
3. **接著**：P1 §6 (SV multi-caller merge)、§7 (欄名對齊)
4. **有空時**：P2、P3 按計畫推進

任何問題我都可以一起 debug。前端 stop-gap 在 GUI repo 分支 `claude/plan-ngs-ui-RQW8J`，commit history 是 1:1 對應修哪個問題。

— 整理人：Claude（GUI 端維護者）
