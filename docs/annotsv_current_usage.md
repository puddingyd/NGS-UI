# NGS-UI 目前的 AnnotSV 使用方式

這份文件只描述 NGS-UI 目前 CNV/SV 的 AnnotSV 輸出從哪裡來，以及如何接回 GUI。

## 主要來源：三級 pipeline（Nextflow）

CNV/SV 的 AnnotSV 註解由正式三級 pipeline 負責，輸出在：

```text
/home/pipeline/tertiary_output/{sample}/06_cnv_sv/{source_sample_id}.cnv.annotated.tsv
/home/pipeline/tertiary_output/{sample}/06_cnv_sv/{source_sample_id}.sv.annotated.tsv
```

NGS-UI worker 在 Nextflow 結束後驗證每個 sample 的 staging 輸出；**缺 CNV 或 SV TSV 會直接讓 job 失敗**，
不再有 NGS-UI 本機 AnnotSV fallback（舊的 `run_annotsv_cnv_sv.sh` /
`annotate_dragen_cnv_annotsv.sh` / `annotate_inhouse_cnv_sv_annotsv.sh` 已移除）。

## DRAGEN integrated CNV rescue

DRAGEN 模式另外有 `scripts/rescue_dragen_cnv.py`（post-processing 的 `cnv-rescue` 步驟），
會從 `{sample}.cnv_sv.vcf.gz` 取 integrated PASS + Rule B CNV 重新跑 AnnotSV，
輸出 `08_postprocessing/{sample}.cnv.review.tsv`。它使用的 AnnotSV runtime 由
`backend/app/services/dragen_cnv_rescue.py::annotsv_command()` 決定：

- 程式：`ANNOTSV_BIN`，或 `$NGS_UI_HOME/biotools/AnnotSV/bin/AnnotSV`，否則使用 `annotsv_3.5.10.sif`（`NGS_UI_ANNOTSV_SIF` 可覆寫）
- annotation database：`ANNOTSV_ANNOTATIONS`，或 `/home/pipeline/reference/hg38/tertiary/annotsv_annotations/share/AnnotSV`、`$NGS_UI_HOME/biotools/AnnotSV/share/AnnotSV`

In-house（NCKUH）模式不跑 rescue。

## 輸出

TSV 是 AnnotSV 標準輸出。GUI adapter 會讀 `Annotation_mode=full` 的事件列，並聚合 `split` gene rows，
主要使用座標、SV type、cytoband、gene、pathogenic/benign overlap、OMIM/GenCC、`AnnotSV_ranking_score`
與 `ACMG_class` 等欄位。
