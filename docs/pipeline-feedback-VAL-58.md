# 三級分析 Pipeline 反饋 — 

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

### 3. Phase-1 TSV 未過濾，5.35M 列 / 1.6 GB

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

- GUI adapter 把整份 TSV 讀進 dict → OOM / GUI 要讀很久
- 先用了 filter: rare (`GNOMAD_G_AF ≤ 0.01`) + `IMPACT ∈ {HIGH, MODERATE}` + 排除 alt-contig + 排除 `*`-allele + 保留 ClinVar P/LP rescue → 5.35M 列 → ~1,363 列
- Impact = LOW 或許可以留著？？

**建議**

在 `PARSE_VEP_CSQ` 之後加 `FILTER` step，輸出兩份：

| 檔名 | 目的 |
|---|---|
| `{sample}.snv_indel.full.annotated.tsv` | 保留現況、archive / re-analysis 用 |
| `{sample}.snv_indel.annotated.tsv` | 過濾完，給 GUI 用 |

可能考慮濾掉 AF > 0.01 + Impact = MODIFIER + alt contig

---

### 4. VEP `--pick` 策略可能選錯 transcript

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
