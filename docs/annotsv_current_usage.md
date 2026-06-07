# NGS-UI 目前的 AnnotSV 使用方式

這份文件只描述 NGS-UI 目前怎麼跑 CNV/SV 的 AnnotSV，以及輸出如何接回 GUI。

## 安裝位置

AnnotSV 目前預設放在 DGM/NGS-UI 本機資料樹底下：

- 程式：`/home/n102968/NGS_UI/biotools/AnnotSV/bin/AnnotSV`
- AnnotSV annotation database：`/home/n102968/NGS_UI/biotools/AnnotSV/share/AnnotSV`


NGS-UI 自己的呼叫程式在：

- `scripts/run_annotsv_cnv_sv.sh`：統一入口，依來源分派 DRAGEN 或 in-house。
- `scripts/annotate_dragen_cnv_annotsv.sh`：DRAGEN CNV/SV VCF。
- `scripts/annotate_inhouse_cnv_sv_annotsv.sh`：in-house gCNV + Delly VCF。


## 輸入

DRAGEN 模式給 hard-filtered SNV VCF，script 會從同一個 `vcf.gz/` 目錄找 sibling 檔案：

- `{sample}.cnv.vcf.gz` → CNV
- `{sample}.sv.vcf.gz` → SV

in-house 模式直接給兩個 VCF：

- gCNV：`{sample}.gcnv.vcf.gz`
- Delly：`{sample}.delly.vcf.gz`

## 執行範例

DRAGEN：

```bash
scripts/run_annotsv_cnv_sv.sh \
  --sample 26WG0001 \
  --out-dir $HOME/NGS_UI/tertiary_output/26WG0001 \
  --dragen-cnv-source /path/to/vcf.gz/26WG0001.hard-filtered.vcf.gz
```

in-house：

```bash
scripts/run_annotsv_cnv_sv.sh \
  --sample 26WE0001 \
  --out-dir $HOME/NGS_UI/tertiary_output/26WE0001 \
  --inhouse-cnv-vcf /path/to/26WE0001.gcnv.vcf.gz \
  --inhouse-sv-vcf  /path/to/26WE0001.delly.vcf.gz
```

## AnnotSV 參數

兩種來源最後都會呼叫 AnnotSV，主要參數一致：

```bash
AnnotSV \
  -SVinputFile <input.vcf.gz> \
  -outputDir <sample_out>/_annotsv_<cnv|sv> \
  -outputFile <sample>.<cnv|sv>.annotated.tsv \
  -genomeBuild GRCh38 \
  -annotationsDir "$ANNOTSV_ANNOTATIONS" \
  -SVinputInfo 1
```

## 輸出

NGS-UI GUI 固定讀 sample 目錄下這兩個檔案：

```text
$HOME/NGS_UI/tertiary_output/{sample}/cnv.annotated.tsv
$HOME/NGS_UI/tertiary_output/{sample}/sv.annotated.tsv
```

TSV 是 AnnotSV 標準輸出。GUI adapter 會讀 `Annotation_mode=full` 的事件列，並聚合 `split` gene rows，主要使用座標、SV type、cytoband、gene、pathogenic/benign overlap、OMIM/GenCC、`AnnotSV_ranking_score` 與 `ACMG_class` 等欄位。
