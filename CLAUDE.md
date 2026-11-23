# AGENTS.md

如果對話窗輸入 go，請做以下：

請只讀 **「待處理事項.md」**，先不要讀其他文件，先判斷哪些是必要讀的文件，然後再去讀和處理「待處理事項.md」中的事項。完成後請在 **「處理紀錄.md」** 摘要你做了什麼，並從「待處理事項.txt」 清掉已經完成的項目（若尚未完成，請留著那個項目，並加註尚未完成的部分）。

**有任何架構變動或功能增減，請順手更新這份 AGENTS.md、README.md 和 frontend/VERSION.md。**

---

## 0. 一句話總覽

成大醫院基因醫學部 **NGS 分析平台**（repo `puddingyd/NGS-UI`）。這份文件記錄整個系統的流程、結構、慣例與踩雷點，讓新的 session / auto-compact 後也能完全接續工作。

醫院內部的 NGS 三級分析判讀工具：FastAPI 後端 + 原生 JS 前端（無 build step），部署在內網 `192.168.84.91:8765`。Nextflow 與 UI state 統一放在 `/home/datalake_Intermediate/pipeline/tertiary_output/{LIS_ID}/`：pipeline 使用 00-07，UI 衍生/判讀狀態使用 `08_postprocessing/`。reviewer 在 UI 裡看 SNV/Indel + CNV/SV + Mitochondria + STR + PGx、標記與匯出報告；另有 `/phenotype/` 與個案清單工具。

---

## 1. Git 工作流（重要）

- 直接在遠端 default branch **`claude/plan-ngs-ui-RQW8J`** 開發並 push（舊的 `Codex/plan-ngs-ui-RQW8J` 已不再使用）；最新、完整的專案說明以 `AGENTS.md` 為準。
- dev 機（`n102968@server`，repo 在 `~/NGS_UI/NGS-UI`，remote 設成 SSH `git@github.com:puddingyd/NGS-UI.git`）從那個分支 `git pull`；如需進 `main` 由使用者在 GitHub 開 PR 合併。
- push 失敗（網路）retry 最多 4 次 exponential backoff（2/4/8/16s）。不要建 PR 除非使用者明說。

---

## 2. 目錄佈局 / env var

`config.py` 裡所有路徑都從 `NGS_UI_HOME` 推導（env `NGS_UI_HOME` 優先；未設定時，只有典型 production checkout `NGS_UI/NGS-UI` 會用 `REPO_ROOT.parent`，一般 standalone checkout 會用 `REPO_ROOT`）。每個都有自己的 env override。在 sandbox/桌面 checkout 裡 `NGS_UI_HOME` = repo root，在 dev 機 `NGS_UI_HOME` = `~/NGS_UI`。

| 內容 | 預設路徑 | env override |
|---|---|---|
| 程式碼（這個 git checkout）= `REPO_ROOT` | `NGS_UI_HOME/NGS-UI/` | — |
| 新三級統一 root（Nextflow 00-07 + UI 08） | `/home/datalake_Intermediate/pipeline/tertiary_output/{LIS_ID}/` | `NGS_UI_TERTIARY_ROOT`（舊 `NGS_UI_PIPELINE_OUT_ROOT` alias） |
| 舊 UI sample root（遷移 fallback；**目錄已移除**，程式只在目錄存在時讀取） | `NGS_UI_HOME/tertiary_output/{LIS_ID}/` | `NGS_UI_LEGACY_TERTIARY_OUTPUT_ROOT`（舊 `TERTIARY_OUTPUT_ROOT` alias） |
| 舊 Nextflow root（遷移 fallback；**目錄已移除**，程式只在目錄存在時讀取） | `/home/pipeline/tertiary_output/{LIS_ID}/` | `NGS_UI_LEGACY_PIPELINE_OUT_ROOT` |
| 未登錄個案共用索引（不含 phenotype 或 roster 個資） | `NGS_UI_HOME/data/unregistered_samples.json` | 隨 `NGS_UI_DATA_ROOT`；啟動預建、5 分鐘過期背景更新，註冊／取消註冊／pipeline 完成或刪除更新單筆 |
| 樣本清單快取（每次 list 都被重寫） | `/home/datalake_Intermediate/pipeline/tertiary_output/_index.json` | `NGS_UI_INDEX_PATH` |
| 伺服器執行狀態（`users.db`、`jobs/`） | `NGS_UI_HOME/data/` | `NGS_UI_DATA_ROOT` |
| 共用個案分類清單／批次復原紀錄 | `NGS_UI_HOME/data/categories.json`、`categories.pending.json`、`categories.lock` | 隨 `NGS_UI_DATA_ROOT`；管理分類儲存即時生效，改名／移除同步套用既有個案 |
| Manual ACMG + active SNV observations SQLite | `NGS_UI_HOME/data/manual_acmg.sqlite` | `NGS_UI_MANUAL_ACMG_DB` |
| SNV review payload in-memory LRU 筆數 | `2` | `NGS_UI_SNV_CACHE_MAX`（`0`=關閉） |
| 完整 SNV annotation fallback 可進 cache 的大小上限 | `100 MB` | `NGS_UI_SNV_CACHE_MAX_RAW_MB`（大型 WGS raw TSV 預設不 cache；gene search 缺 index 時串流掃描） |
| 病患 phenotype.txt + Clinical presentation sidecar | `NGS_UI_HOME/patient_phenotype/`；主檔為 `{MRN}_phenotype.txt` / `{MRN}_clinical_presentation.txt`，舊 `{LIS}_{MRN}` / `{LIS}` phenotype 只作 fallback | `NGS_UI_PHENOTYPE_DIR` |
| MRN-level 病歷附件（PDF/JPG/PNG/TIFF）+ SQLite catalog/audit + trash | `NGS_UI_HOME/patient_documents/` | `NGS_UI_PATIENT_DOCUMENTS_DIR`；串流上傳不設 per-file cap，但保留 `NGS_UI_PATIENT_DOCUMENTS_MIN_FREE_GB`（預設 10 GB，`0`=停用 reserve） |
| 上傳的個案清單 xlsx + `roster.json` | `NGS_UI_HOME/patient_list/` | `NGS_UI_PATIENT_LIST_DIR` |
| HPO reference（`hp.obo`、`phenotype_to_genes.txt` 等大型檔） | `NGS_UI_HOME/phenotype_data/` | `NGS_UI_PHENO_DATA_DIR` |
| GenCC / ClinGen / MONDO gene-disease SQLite index | `NGS_UI_HOME/phenotype_data/gene_disease/gene_disease.sqlite` | `NGS_UI_GENE_DISEASE_DB` |
| GenCC / ClinGen / MONDO raw downloads | `NGS_UI_HOME/phenotype_data/gene_disease/raw/` | `NGS_UI_GENE_DISEASE_RAW_DIR` |
| 固定 WES-I / WES-II / WGS gene panels | `REPO_ROOT/phenotype_data/gene_panels/` | `NGS_UI_GENE_PANELS_DIR` |
| 固定 panel UI index | `REPO_ROOT/phenotype_data/fixed_panels/` | `NGS_UI_FIXED_PANELS_DIR` |
| custom gene panels + 名稱／來源／輸出名稱對照表 | `REPO_ROOT/phenotype_data/custom_panels/`（`panel_metadata.tsv`） | `NGS_UI_CUSTOM_GENE_PANELS_DIR` |
| expanded reportable gene list / HGNC alias / dead-zone tables | `REPO_ROOT/ngs_panel_deadzone/` | `NGS_UI_PANEL_DEADZONE_DIR` |
| OMIM.xlsx | `NGS_UI_HOME/OMIM/OMIM.xlsx` | `NGS_UI_OMIM_XLSX` |
| OMIM 擴充功能擷取原文（append-only gzip JSON + `index.jsonl`） | `NGS_UI_HOME/OMIM/captures/` | `NGS_UI_OMIM_CAPTURE_DIR` |
| 平台 OMIM 資料庫（存在時取代 OMIM.xlsx 作為來源） | `NGS_UI_HOME/OMIM/omim.sqlite` | `NGS_UI_OMIM_DB` |
| GIAB stratification BED + `strata_manifest.json` | `NGS_UI_HOME/biotools/giab_stratification/` | `NGS_UI_GIAB_STRAT_DIR` |
| GeneBe ACMG 主 DB（bgzip TSV；lazy 產生同目錄 SQLite cache；三級大批次優先，Somatic 小批次作 live API fallback） | `NGS_UI_HOME/biotools/genebe/genebe_hg38.tsv.gz` | `NGS_UI_GENEBE_DB` |
| GeneBe live API 成功／negative cache（三級只查主 DB miss 且符合 review filter 的點；Somatic 預設 ≤10,000 個候選時查全部具體 allele） | `NGS_UI_HOME/biotools/genebe/genebe_api_cache.sqlite` | `NGS_UI_GENEBE_API_CACHE` |
| GeneBe API 待匯入 7 欄 TSV + JSON sidecar | `NGS_UI_HOME/biotools/genebe/api_pending/` | `NGS_UI_GENEBE_API_PENDING_DIR` |
| LitVar2 官方 bulk JSON + slim SQLite（只供 reviewer 文獻脈絡，不影響 tier/ACMG/report） | `NGS_UI_HOME/biotools/litvar2/{litvar2_variants.json.gz,litvar2.sqlite}` | `NGS_UI_LITVAR2_DIR`、`NGS_UI_LITVAR2_BULK_PATH`、`NGS_UI_LITVAR2_DB`、`NGS_UI_LITVAR2_BULK_URL` |
| 每週最新版 ClinVar UI 比對庫 | `NGS_UI_HOME/biotools/clinvar_latest/{clinvar.vcf.gz,clinvar_latest.sqlite}` | `NGS_UI_CLINVAR_LATEST_DIR`、`NGS_UI_CLINVAR_LATEST_DB`、`NGS_UI_CLINVAR_LATEST_VCF`、`NGS_UI_CLINVAR_LATEST_URL`；每週三 04:00 timer 原子更新，只影響 UI overlay |
| Research-only SpliceAI reference | `NGS_UI_HOME/biotools/spliceai/spliceai_scores.raw.{snv,indel}.hg38.vcf.gz` | 兩檔與各自 `.tbi` 都必須存在；可由 `run_stopgaps.sh --spliceai-snv/--spliceai-indel` 覆寫 |
| Research-only AutoMap image | `/home/datalake_Intermediate/pipeline/nextflow_containers/automap_1.3.sif`（舊 `/home/pipeline/nextflow_containers/automap_1.3.sif` fallback） | `NGS_UI_AUTOMAP_SIF`；只在 in-house 缺既有 AutoMap 且三級勾 Research-only 時使用 |
| GPN-Star V/M GRCh38 canonical LLR（僅 Research-only） | `NGS_UI_HOME/biotools/gpn_star/`：M `manifest.json`、V `manifest.v.json` + 各自 `gpn-star-hg38-{m447,v100}-200m/llr/llr_chr{1-22,X,Y}.parquet` | `NGS_UI_GPN_STAR_DIR`；僅明確 Research-only 才由 Polars join compact review SNV，不進 raw/overlay/gene search；固定 revision，無 timer；缺資料／Polars／查詢失敗只 warning，不阻擋三級 |
| GPN-MSA（所有三級 case 與 Somatic） | `NGS_UI_HOME/biotools/gpn_msa/scores.tsv.bgz` + `.tbi` | `NGS_UI_GPN_MSA_DB`；三級與 Somatic 都維持 tabix best-effort annotation；MSA、Star 各自欄位，並存顯示 |
| 本院 AF sites VCF（院內 WGS cohort；核基因組用 `INHOUSE_AC/AN/AF`，chrM 另以 carrier frequency 定義並帶 `INHOUSE_NHOM/HET_MT`；缺檔＝本院 AF 靜默關閉） | `NGS_UI_HOME/biotools/inhouse_af/inhouse_af.hg38.vcf.gz` | `NGS_UI_INHOUSE_AF_DB` |
| Exomiser/LIRICAL CLI | `NGS_UI_HOME/biotools/` | `NGS_UI_BIOTOOLS_DIR` |
| VCF（per-sample） | `NGS_UI_HOME/vcf/` | `NGS_UI_VCF_DIR` |
| 三級分析 DRAGEN VCF 搜尋 roots | `/home/datalake_Raw/Novaseq:/home/datalake_Intermediate/n102968` | `NGS_UI_DRAGEN_VCF_ROOTS`（可用 `:` 分隔多個 root/glob；尋找 `*/vcf.gz/*hard-filtered.vcf.gz`） |
| IGV in-house BAM 搜尋根目錄 | `/home/datalake_Intermediate/pipeline/nextflow_output` | `NGS_UI_INHOUSE_BAM_ROOT`（可用 `:` 分隔多個 root；舊 `NGS_UI_BAM_ROOT` 仍作 fallback） |
| IGV DRAGEN BAM 搜尋根目錄 | `/home/datalake_Raw/Novaseq:/home/datalake_Intermediate/n102968` | `NGS_UI_DRAGEN_BAM_ROOT`（可用 `:` 分隔多個 root；尋找 `<run>/bam/{sample}.bam`，排除 `{sample}.repeats.bam`；systemd 若有設定此 env 需同步併入測試 root） |
| IGV hg38 reference（FASTA + `.fai`） | `/home/pipeline/reference/hg38` | `NGS_UI_IGV_REF_DIR` |
| 二級分析 WES FASTQ 搜尋 roots | `/home/datalake_Raw/NextSeq2000:/home/datalake_Raw/Other/Reanalysis:/datalake_Raw/Other/Reanalysis` | `NGS_UI_SECONDARY_WES_FASTQ_ROOTS` |
| 二級分析 WGS FASTQ 搜尋 roots | `/home/datalake_Raw/Novaseq`（WGS 下拉選單依主 sample 聚合並顯示 lane/FASTQ 數；samplesheet 展開為每 lane 一列，不用 merged） | `NGS_UI_SECONDARY_WGS_FASTQ_ROOTS` |
| 二級分析 samplesheet server output | `/home/datalake_Intermediate/pipeline/nextflow_output` | `NGS_UI_SECONDARY_OUTPUT_ROOT` |
| 二級分析 samplesheet staging | `/home/datalake_Intermediate/pipeline/nextflow_samplesheet_staging` | `NGS_UI_SECONDARY_SAMPLESHEET_STAGING_ROOT`（DGX2 對應 `NGS_UI_SECONDARY_DGX_SAMPLESHEET_STAGING_ROOT`） |
| 二級分析 DGX2 raw root（對應 DGM `/home/datalake_Raw`） | `/datalake_Raw/datalake_Raw` | `NGS_UI_SECONDARY_DGX_RAW_ROOT` |
| 二級分析 DGX2 output / launch / work roots | `/datalake_Intermediate/pipeline/nextflow_output`、`/datalake_Intermediate/pipeline/nextflow_launch`、`/raid/DGM/work` | `NGS_UI_SECONDARY_DGX_OUTPUT_ROOT`、`NGS_UI_SECONDARY_DGX_LAUNCH_ROOT`、`NGS_UI_SECONDARY_DGX_WORK_ROOT` |
| 三級分析 Nextflow cache/work root | `NGS_UI_HOME/nf_work/` | `NGS_UI_TERTIARY_NF_WORK_ROOT`；DRAGEN 與 NCKUH 分別共用 `contexts/shared-dragen/`、`contexts/shared-nckuh/` 的 launch/session lineage；新 lineage 預設用其 `work/`，承接舊 session 時以 `work-dir.txt` 固定沿用原 work root |
| PRS 個案重建／Nextflow／工具快取 | `NGS_UI_HOME/nf_work/prs/{個案ID}/`（預設沿用三級 Nextflow work root） | `NGS_UI_PRS_WORK_ROOT`；未設定時跟隨 `NGS_UI_TERTIARY_NF_WORK_ROOT/prs/`，每案獨立 project；正式計分／QC／風險快照存08；三級暫存清理明確保留PRS，完整刪個案才清該案work |
| Exomiser/LIRICAL 輸入模板 | `REPO_ROOT/phenotype_reference/`（**在 repo 裡，不在 NGS_UI_HOME**） | — |
| 前端靜態檔 | `REPO_ROOT/frontend/` | `FRONTEND_DIR` |
| 首頁歡迎文字與版本紀錄 | `REPO_ROOT/frontend/VERSION.md`（前端啟動時讀取並顯示在未載入個案的首頁） | — |
| 首頁 NGS 分析流程圖文字 | `REPO_ROOT/frontend/ANALYSIS_FLOW.md`（固定 section/key schema；可編輯 `|` 後文字） | — |
| EMR client id（NCKU intranet） | — | `NGS_UI_EMR_CLIENT_ID`（空 = 整套 EMR 功能關閉） |
| Redis（job queue） | `redis://127.0.0.1:6379/0` | `REDIS_URL` |
| Java / Exomiser 路徑等 | 見 `config.py` | `EXOMISER_HOME`、`LIRICAL_HOME`、`JAVA_BIN`… |

`.gitignore`：`tertiary_output/`、`data/`、`patient_list/`、`patient_documents/`、`biotools/`、`phenotype_data/*`（但例外追蹤 `phenotype_data/gene_panels/*.txt`、`phenotype_data/custom_panels/*.txt`、`phenotype_data/fixed_panels/**/*.txt` 與 `phenotype_data/fixed_panels/index.json`）、`_index.json`、`__pycache__/`、`*.pyc`、`.venv/`、`node_modules/`。所以 fixed/custom panel 檔在 git 裡，GeneBe DB/API cache/pending TSV、病歷附件與 audit DB、HPO 大型 reference、樣本資料、roster 都不在 git 裡，部署時要自己放到 `NGS_UI_HOME/` 底下。`OMIM.xlsx` 原本在 repo 根，已 `git rm` 移出（dev 機放在 `NGS_UI_HOME/OMIM/`）。

**每樣本目錄 `/home/datalake_Intermediate/pipeline/tertiary_output/{LIS_ID}/` 內容：**
```
03_acmg/{source}.snv_indel.acmg.tsv    唯一完整 SNV source of truth（immutable）
03_acmg/{source}.annotation_versions.json  annotation DB version sidecar（至少 ClinVar release_date）
04_mito/{source}.mito.tsv              Mito source，UI 直接讀
05_str/{source}.str.tsv                STR source，UI 直接讀
06_cnv_sv/{source}.{cnv,sv}.annotated.tsv  CNV/SV source，UI 直接讀
07_pgx/{source}.pgx.tsv / PharmCAT JSON    PGx source，UI 直接讀
08_postprocessing/
  {LIS_ID}.cnv.review.tsv      DRAGEN 整合檔 PASS＋Rule B 全部重新 AnnotSV；完成 manifest 存在時優先讀取
  {LIS_ID}.cnv.rescued.annotated.tsv  僅舊 v1 存在；v2 不再產生，重跑由 managed promotion 清除
  {LIS_ID}.cnv_rescue.json     v2 規則版本、整合來源 SHA256、整合後座標與 FILTER/QUAL、支持 SV、計數
  {LIS_ID}.layout.json        layout v3 activation marker；最後原子寫入
  {LIS_ID}.pipeline_source.json  source_path/source_sample_id/source_vcf_path/source_vcf_size/annotated_at
  {LIS_ID}.specialized.manifest.json  四 locus 資料包版本／checksum／generation／來源追溯
  {LIS_ID}.specialized.source.json  完整小型 targeted JSON 快照
  {LIS_ID}.specialized.evidence.jsonl.gz  四 locus 原始 VCF headers／records／同源區證據
  {LIS_ID}.specialized.calls.json  GBA1／CYP21A2／SMN／HBA 正規化 findings 與 locus 狀態
  {LIS_ID}.snv_annotations.sqlite  ClinVar latest/GeneBe/GIAB/院內 AF/MANE/SpliceAI/LitVar2 sparse overlay
  {LIS_ID}.clinvar_comparison.json baseline/latest release、比對統計與升降規則 marker
  {LIS_ID}.snv_indel.review.tsv    主畫面衍生檔；WES DP≥20 + BED/AF + ClinVar rescue
  {LIS_ID}.snv_indel.review.tsv.source.json  raw/overlay/filter signature + GPN-MSA／Star 狀態 + gpn_star_research_only
  {LIS_ID}.snv_gene_index.sqlite   gene / variant id → 03_acmg raw byte offset
  {LIS_ID}.litvar2_annotation.json  filtered post-processing 完成 marker + DB 版本
  {LIS_ID}.litvar2_on_demand.sqlite gene search / 手動單點查詢的版本化正負結果快取
  {LIS_ID}.sample_metadata.json / {LIS_ID}.case_summary.json
  {LIS_ID}.mito.annotated.tsv      只有確實補到 MITOMAP 時才保留的小型衍生檔
  {LIS_ID}.ploidy.vcf.gz / {LIS_ID}.ploidy_qc.txt
  {LIS_ID}.qc_summary.json / {LIS_ID}.vcf_from_tsv.vcf.gz
  {LIS_ID}.roh_regions.tsv / {LIS_ID}.roh_summary.json
  {LIS_ID}.roh.source.{automap.tsv,automap.pdf,bcftools.txt,dragen.bed,dragen_metrics.csv}
  analyses/{ver}/{LIS_ID}.analysis.json, {LIS_ID}.pheno_score.tsv, prefixed results/analysis_files
```
layout v3 的 sample-owned 檔案一律用 `{LIS_ID}.<name>`；resolver 固定 prefixed exact name 優先，再 fallback 舊 unprefixed exact name，不可用 broad wildcard。舊 layout v2 `layout.json` 及未加前綴檔仍可讀寫；完整重跑會非破壞性複製 v2 state，最後才寫 prefixed v3 marker，不刪舊檔。新樣本只有 marker 存在且沒有 metadata 時才列為未登錄，避免 post-processing 未完成就被載入。

SNV post-processing 允許執行期間建立隱藏 working TSV，既有 in-place annotator 跑完後比較 03 raw 產生 sparse overlay，worker 的 `finally` 一定刪除 working TSV；不得在 08 永久留下第二份 `snv_indel.annotated.tsv`。`REF/ALT=*` 與非 primary contig 在 review/index/search/VCF consumer 端排除，不改寫 raw。

新版三級 v3.8 的 SNV/Indel source schema 為 82 欄，新增 `HAPLOID_HET` 並使 `ACMG_CLASS` 移至第 78 欄；所有 GUI consumer 必須依欄名解析，不可依位置。嚴格 staging 驗證要求 `HAPLOID_HET`；`CALLERS=NONE` 在 review、gene index、VCF、gene search 與 adapter 共用的 raw-row eligibility 中防禦性排除。非空 `HAPLOID_HET`（男性 chrX 非 PAR 原始 het）保留為候選並在卡片顯示紅色人工複核 badge，含 DV／HC 來源；不因此修改 tier、ACMG 或報告選取。AD 的 `.` 保持缺值，ZYGOSITY 直接採 pipeline 依實際 ALT caller 與 GT ploidy 算出的值；GUI 不自行推測性染色體 zygosity。因 `add_callers_tag.py`／`parse_vep_csq.py` 的內容變更不一定影響 Nextflow task hash，v3.8 首次部署後若 strict validation 收到舊 81 欄結果，錯誤訊息須指出 stale `-resume` cache；在沒有 active job 時清理一次三級 Nextflow 暫存再重跑。二級 ensemble 已重跑、輸入檔已改變時通常會自行失效。 最新版 DRAGEN compound 規則不新增欄位：只有 phase 已知且同 phase set 的 het 合成（有 `|` 但無 PS 依 VCF 規格視為同一 set），重疊只允許落在 indel 錨定鹼基；其他重疊與 symbolic ALT 由 pipeline 原樣保留。同一鹼基由 pipeline 先套 SNV 再套 indel；GUI 不做 overlap／phase-based 合成或去重，必須保留不同 variant ID。DRAGEN 多等位拆分後若殘留 `ZYGOSITY=ref`／`GT=0/0` 假列，raw-row eligibility 與 `CALLERS=NONE` 一樣排除；不要單看某一 caller 的 `GT_DV=0/0`，NCKUH 可能由另一 caller 真正 call 到 ALT。

---

## 3. 整體資料流

```
三級 pipeline (Nextflow)
  → 每次都啟動 `-resume`；同 pipeline 類型共用 session，跨 batch 沿用相同 sample/input cache
  → `.staging/{job_id}/{source_sample_id}/00-07`
  → 驗證 00-06 + SNV/Mito/STR/CNV/SV（勾選時再驗證 07/PGx）
  → staging 內完成 08 post-processing
  → 可 rollback promotion 到 `/home/datalake_Intermediate/pipeline/tertiary_output/{LIS_ID}/`
  → NGS-UI worker 直接讀 03-07；不複製 SNV/CNV/SV/STR/PGx
  → SNV working TSV（暫存；SpliceAI、GeneBe API、LitVar2 共用 final-review predicate）→ sparse overlay + LitVar2 marker → final review TSV 以 tabix 補 GPN-MSA，Research-only 再以 Polars 補 GPN-Star V/M（各自 best-effort）+ prefixed index → 刪 working TSV
  （2026-11-19 起整批「逐步驟」執行：`scripts/postprocess_batch.py` 對一個 chunk 的所有 sample 每個資源只查一次，再平行寫回各 sample，詳見 §11）
  → 08_postprocessing/{LIS_ID}.layout.json 最後寫入，啟用 UI 讀取
  → Mito 只有 MITOMAP 確實補值時才留 08_postprocessing/{LIS_ID}.mito.annotated.tsv
  → ROH 只寫 08：in-house 優先二級 AutoMap；缺 AutoMap 時 Research-only 補跑，否則/失敗用 BCFtools；DRAGEN 永遠用原生 ROH
  → (HPO/panel: 載入新個案、/phenotype/ 與主畫面 default 組合共用 patient_phenotype/{MRN}_phenotype.txt；主畫面其他 analysis version 不覆蓋主檔；舊 LIS 檔只作 fallback)
  → (Clinical presentation: /phenotype/ 工具依 MRN 產生 patient_phenotype/{MRN}_clinical_presentation.txt；主畫面 reviewer 修改後會同步寫回)
  → (Documents: 主畫面與 /phenotype/ 登入後依 MRN 共用 patient_documents/files/{MRN}/{UUID}；SQLite 管理顯示檔名、hash 與 audit，刪除移入 trash)
  → (個案清單: 上傳 xlsx → patient_list/roster.json)

UI 流程：
  載入新個案 (modal；未登錄個案欄位支援 LIS ID/source sample/姓名/MRN typeahead 搜尋，清單先顯示快取並背景更新，選取後才讀個案詳細資料，可按「更新清單」強制掃描) → POST /api/samples → 寫 sample_metadata.json + analyses/default/analysis.json
                                          + (write_version side-effect 寫 pheno_score.tsv；SNV in_panel 動態由 pheno_score.tsv 補值)
                                          + 有至少一個 HPO term 才排 Exomiser/LIRICAL；缺 VCF 時由 worker 背景建立
  個案清單 (modal) → DELETE /api/samples/{id}?delete_pipeline_output=false
                  → 只刪 sample_metadata.json、case_summary.json 與 analyses/，保留 annotation/PGx/ploidy 等 local tertiary outputs
                  → sample 從已載入個案清單消失，並回到「載入新個案」的未登錄清單；真正刪 pipeline output 走三級分析清單 / DELETE /api/dragen/outputs/{id}
  選樣本 (combobox) → GET /api/samples/{id}  (核心 payload，aux_pending=true)
                    → 背景分別 GET /api/samples/{id}/secondary-snv、/cnv、/sv 和 /mito（SNV 先顯示；Secondary findings、CNV 可各自載完先顯示）
                    → GET /api/samples/{id}/report  (reviewer 編輯狀態)
  reviewer 在卡片標 1/2/C/0、secondary findings 勾 ✓、寫 comment、勾 disease/gene checkbox、用 structured modal 改 ACMG
  自動儲存 (1.5s debounce) → PUT /api/samples/{id}/report
  ▶ 開始分析 → POST /api/samples/{id}/phenotype (存 HPO/panels + 算 pheno_score)
            → POST /api/samples/{id}/jobs/exomiser_lirical (enqueue rerun worker)
  匯出診斷報告 → GET /api/samples/{id}/report.docx
  輸出 PDF → 瀏覽器列印視窗（報告區 causative / other / candidate 卡片摘要 + 本次檢測基因清單）
```

舊個案遷移：先 `python scripts/migrate_tertiary_output_layout.py --sample ID` dry-run，再對 1–2 個 case 加 `--apply`。驗證 UI、gene search、報告與 metadata checksum 後才 `--all --apply`。工具不刪 legacy tree，且只在所有衍生檔驗證通過後寫 marker；`--sample ID --rollback --apply` 只停用 marker，立即回 legacy。遷移已完成：兩個舊目錄（`~/NGS_UI/tertiary_output`、`/home/pipeline/tertiary_output`）已移除，`--rollback` 回 legacy 已不可行；所有個案都在 `/home/datalake_Intermediate/pipeline/tertiary_output/{LIS_ID}/`（例如 `26T00970-dragen/03_acmg/26T00970.snv_indel.acmg.tsv`）。

---

## 4. Adapter / tier 結構

每種 variant type 一個 adapter（`backend/app/adapters/`），各回傳 `(variants_dict, categories_dict)`。`sample_loader.load_sample()` 把它們全部塞進一個 payload —— id namespace 不衝突（prefix 不同）。

| Type | adapter | id 格式 | tiers（payload key） | tier 規則 | 排序 |
|---|---|---|---|---|---|
| **SNV/Indel** | `snv_tsv.py` | `chr{N}-{pos}-{ref}-{alt}` | `1A 1B 1C 2 3`（互斥） | `classify_tier`：1A = ClinVar P/LP ≥1★；1B = ClinVar P/LP（任何）或 LOFTEE HC；1C = ACMG points ≥4；2 = ACMG points 1-3；3 = 其餘。`_normalize_acmg_class` 把 `VUS`→`Uncertain significance` 等 | 各 tier 內 `total_score`（= geno_score + pheno_score）desc，tie-break by id |
| **CNV** | `annotsv_tsv.py`（`source="cnv"`） | AnnotSV_ID | `CNV-1A`(Clinical) `CNV-1B`(Pathogenic)（**獨立分區，可重複**） | 1A = SV 涵蓋的任一 gene 在 pheno set（score>0）；1B = AnnotSV `ACMG_class` ∈ {4,5} | 各 tier 依 `max_pheno_score + scaled AnnotSV_ranking_score` desc → id；報告區與 DOCX 的 CNV/SV 也依此 combined score 排。基因表 `genes` 切到前 10（後端，`genes_overflow` = in-panel 溢出、`genes_compact` = 其餘只含 `{gene,omim_id,in_panel}`） |
| **SV** | `annotsv_tsv.py`（`source="sv"`） | AnnotSV_ID | `SV-2A` `SV-2B` | 同 CNV，用 sv 檔 | 同 CNV |
| **Mito** | `mito_tsv.py` | `chrM-{pos}-{ref}-{alt}` | `MITO-1`(Pathogenic) `MITO-2`(Rare / reported mtDNA variant) `MITO-3`(Other variant)（互斥） | 支援舊 `mito.annotated.tsv` 與 v3.2 `04_mito/{sample}.mito.tsv` 欄位；worker 複製 v3.2 mito TSV 後會靜默回補本地 MITOMAP 欄位（只供分類，不顯示在卡片、不額外寫 log）。**若有 `FILTER` 欄則先 `FILTER=PASS` 過濾**，無 `FILTER` 欄則不做 UI filter。1 = ClinVar P/LP（排除 conflicting）或 MITOMAP `Cfrm`/`Confirmed`/`[P]`/`[LP]`/明確 pathogenic；2 = 非 ClinVar B/LB 且 gnomAD-mito 缺值或 rare（`max(GNOMAD_MITO_AF_HOM, GNOMAD_MITO_AF_HET, GNOMAD_MITO_AF) < 0.01`），或有 ClinVar/MITOMAP disease 記錄；3 = 其他 PASS mtDNA variants。 | `in_panel`、heteroplasmy、depth、gnomAD-mito AF、position 排序；heteroplasmy 只作 tie-break，不作致病性門檻。 |

**報告區**（`REPORT_SECTION_DEFS`，全在 `frontend/index.html` 的 `#report-sections` + Secondary-findings 折疊群組）：
- Causative（status `1`）、Other（`2`）、**Candidate（`C`）** —— 三段 default open，有 disease checkbox、可「＋ 新增 variant」。三段先按 `total_score` desc 排，再把**同基因的 cluster 在一起**（最高分基因的整組排最前；手動新增的無 gene_symbol 留原位）。
- ACMG SF / Proactive / Carrier / PharmCat —— 收在「Secondary findings」折疊群組（純文字標題 + 三角形鈕，無卡片框）。
- status 圓形按鈕：**`1/2/C/0`**（C → Candidate 區）。`C` 與 `0` 可並存（metadata 存成 `C,0`），`1` 或 `2` 會反選其他狀態；再次點擊已選項目可清空；同一 variant 在分析區、報告區與搜尋 modal 的按鈕會同步上色。
- 「輸出 PDF」在三組匯出按鈕旁，使用瀏覽器列印視窗輸出報告區卡片摘要：標題 `{LIS_ID} Report`，右上固定 `{LIS_ID}_report_{YYYY/M/D HH:mm}`、右下頁碼；含 causative/other/candidate，不含 comment、More 展開內容、Secondary findings、操作按鈕與 CNV/SV 已知致病/良性區域重疊。OMIM disease summary 只留到遺傳模式括號。causative/other 空區仍保留，candidate 空區省略。
- DOCX 診斷報告依序輸出第一類、第二類、固定家族檢測建議，再集中輸出各 variant 的「參考資料」。CNV/SV 的 Disease override 在多基因或無 OMIM gene 片段中會直接併入第 1 點位置描述（`...，與 {Disease} 相關。`），不另列一點，因此後續編號前移；單一 gene 仍使用含疾病與遺傳模式的詳細句型。SNV/Indel 表格會帶 TSV `RS_ID`，RS ID 欄預留尾端空格避免黏到結構欄；SNV/Indel 與 Mito 的核苷酸欄內容每行最多 13 字元，蛋白質 `(p....)` 固定另起一行。遺傳模式 `XL / XLD / XLR` 分別顯示為「性聯遺傳 / 性染色體顯性遺傳 / 性染色體隱性遺傳」，`Likely pathogenic` 中文為「疑似致病性」。§五.4 grouped gene list 不縮排、HPO 標題不再重複顯示 `(HP:...)`，且各 HPO/panel 區塊之間留空行；§五.4 只列 `ngs_panel_deadzone/panel/panel_loose.hgnc_canonical.txt` 內的 disease-associated genes。WGS 報告會把 dead-zone genes 以 `GENE（exon 2, 4-6 <15X）` 短括號標註，基因清單後空一行並以無縮排註解「註：括號中標示之 exon 為 cohort dead-zone，代表該 exon coverage 低於本檢測判讀門檻（<15X）。」說明；WES 報告不輸出 dead-zone 標註，僅在主畫面顯示。

**「in_panel」概念**：`pheno_score.tsv` 裡 score>0 的基因 = 病人 HPO/panel 相關的基因（`phenotype_scorer.compute_pheno_match` 回 `{gene: matched_weight}` + total_weight；`compute_pheno_score` = 之後做 `100*matched/total` 正規化）。CNV/SV 的 Clinical tier、Mito 的排序 tie-breaker、CNV/SV 卡片基因表的 ⭐ 標記都用這個。CNV/SV 卡片「Pheno」欄顯示成 `matched/total`（乘 100 前的原始狀態）。`has_phenotype` = bool(hpo or panels)，前端用來在 CNV/SV Clinical 區空白時顯示「請先設定 phenotype」提示。`Disease-associated` 是另一層報告基因範圍，來自 `panel_loose(_plus_clinical).hgnc_canonical.txt`，v10.25 起再聯集平台 OMIM 資料中有疾病 MIM 的基因（`{ }`、`?` 算，`[ ]` 不算）；SNV 主畫面預設用它過濾，DOCX §五.4 也只列這個範圍內的 HPO/panel genes。Mito disease 來自 `CLINVAR_DN`，前端用 `&` 拆成 checkbox，勾選者才進 DOCX 疾病句。`panel_deadzone.py` 會用 `HGNC_ID` 優先把 VEP symbol 轉成 HGNC-current 名稱，並提供 WES 20X / WGS DRAGEN 15X dead-zone exon 查詢。

---

## 5. 前端版面 / 卡片（`frontend/index.html` + `app.js` + `style.css`）

- **Topbar**（深色 `#24292f`，z-index 110 蓋過 login modal）：左 hamburger（toggle `#sidebar`），中標題「成大醫院基因醫學部 NGS 分析平台」，右 `登入/登出`（同一顆鈕 toggle：`data-loggedIn` 切 handler）·`輸入臨床表徵 (HPO/panel)`（`<a href="/phenotype/" target="_blank">`，**未登入也顯示**）·`上傳個案清單`（xlsx upload）。`.btn` 用 `text-decoration:none` + inline-flex，所以 `<a class="btn">` 跟 `<button class="btn">` 一樣。
- **首頁歡迎 / 版本紀錄**：未載入個案時，搜尋卡片下方顯示 `#welcome-card`，內容由 `frontend/VERSION.md` 讀取並用前端的簡易 Markdown renderer 呈現；載入任一 sample 後自動隱藏。之後若有影響判讀流程、報告輸出、資料載入或主要工具入口的更新，要評估是否同步更新 `frontend/VERSION.md`。
- **登入 modal**：「申請帳號：PYTHONPATH=backend python -m app create-user」那段提示字用 `.login-hint`/`.login-hint-code`（白底白字，反白才看得到）。
- **三個分析卡片**（`#card-snv`、`#card-cnv-sv`、`#card-mito`）各有 tier-tab bar + tier panels。`renderTierTabBar`/`renderCnvSvTabBar`/`renderMitoTabBar`；tab-click dispatch 統一處理三組（用 `data-tier` 判斷）。tier-panel 顏色：SNV 紅/黃系、CNV 藍 `#bfdbfe`、SV 紫 `#ddd6fe`、Mito teal `#99f6e4`。卡片要包在 `.block-body` 裡才有 inset 效果（`.tier-panel > .block-body { padding-top: 8px }`）。Mito/STR/ROH 卡片：STR 跟 ROH 還是「（無資料）」placeholder。
- **Staged loading**：`GET /samples/{id}` 只回核心（meta + reports + review TSV 的 SNV + analyses + has_phenotype，`aux_pending: true`，CNV/SV/Mito 是空 dict）。前端 `loadSample` render 完後背景分別打 `GET /samples/{id}/cnv`、`/sv` 和 `/mito`，CNV 載完會先顯示，SV 可繼續在背景載入；回來後 merge 進 `state.data` + re-render 那張卡。`state._auxLoadToken` 防 race（切樣本後晚到的回應丟掉）。等待時對應 CNV/SV、Mito panel 顯示「載入中…」、tab count「…」。SNV tier 只 render active tab，切 tab 時才建立該 tier 卡片 DOM。
- 載入或重新載入個案時會重設 sample-scoped UI state：SNV/CNV/Mito/STR active tab 回預設、主畫面 gene search input 與 gene search modal 舊結果清空。
- **SNV 顯示 filter**：分析區標題列依序是 `Disease-associated`（預設勾；限 `panel_loose`，只影響 SNV 主畫面，不限制 CNV/SV 或 gene search modal）、`In panel only`（預設勾）、`gnomAD_G_AF < 0.01`（預設勾；缺值保留）、`VAF < 0.2`（預設不勾；主畫面預設隱藏低 VAF 或 ref genotype 點位，勾後才顯示）、`impact=MODIFIER`（預設不勾；勾後才顯示 MODIFIER）、OMIM 顯示 toggle。這些都是 UI-only，不刪原始 TSV；主畫面 review TSV 預載層只保留 `NGS_UI_CDS_CANDIDATE_BED` 內且 AF<0.01/AF缺值的點，ClinVar P/LP 與 reviewer 已標記點不受 BED/AF 限制，其他超過者仍可用 gene search 查回。`IMPACT=LOW` 一律保留。tier tab count 跟著 filter 重算。若 reviewer 標記不在 `panel_loose` 的 SNV，分析區會提示 DOCX §五.4 不會列該 gene，但不阻擋標記。
- **ESM1b 配色**：依 Bergquist et al. 2025 ClinGen SVI 校準區間壓成五色；ESM1b 與一般 predictor 方向相反，分數越低越偏致病（`≤-24` P、`≤-10.7` LP、`≤-6.2` VUS、`≤8.7` LB、其餘 B）。多 transcript TSV 用最低分作 worst case。
- **GIAB stratification badge**：SNV/Indel 卡片那排 `MANE_SELECT / DRAGEN / In panel / LOFTEE HC` 標籤後面，會依變異落在哪些 GIAB 困難區（homopolymer / tandem_repeat / segdup / low_mappability / gc_extreme / other_difficult）多畫琥珀色 `.badge-giab`。資料來自 TSV `GIAB_STRATA` 欄（逗號分隔 label）→ adapter `giab_strata: [...]` → 前端 `GIAB_STRATA_DISPLAY` map（不認得的 label 直接顯示原字串，加 strata 不用改前端）。純 UI QC 提示，不影響 tier 排序、不進 DOCX。`GIAB_STRATA` 由 `scripts/annotate_giab_strata.py` 在三級分析尾段寫入（見 §7），review TSV 與 gene search 都會帶這欄。BED 不在 git，部署時用 `scripts/download_giab_strata.sh` 放到 `NGS_UI_GIAB_STRAT_DIR`。
- **SNV OMIM 外部連結**：SNV/Indel 卡片的 OMIM 按鈕優先用 TSV `OMIM_LINK` / `OMIM_id` / `omim_ids`；若沒有 OMIM disease/ID 但有 `gene_symbol`，仍顯示 OMIM 並連到 OMIM geneMap 搜尋該 gene，方便確認資料庫是否更新。
- **個案搜尋 filter**：主畫面搜尋框上方是可複選、可取消的 `WES` / `WGS` 圓形 chips；`matchSamples()` 與 Enter 精準查找都只在勾選的 test type 內搜尋。提示文字放在 input placeholder。
- **個案清單管理**：modal 內有獨立的 `WES` / `WGS` 圓形 chips 與全文搜尋框，會搜尋 ID、姓名、MRN、variant、疾病、comment 與日期等所有表格內容。首頁 `/api/samples` 只回輕量搜尋索引；打開個案清單才呼叫 `/api/samples/case-summary` 載入摘要。表格摘要 status `1` 的 causative SNV/CNV/SV、實際勾選的 SNV OMIM diseases 與 CNV/SV reviewer `Disease`、status `2` 的 other SNV/CNV/SV、主畫面 comment、簽收時間與載入時間；多個 variant / disease 逐行顯示，長 HGVS 可在任意字元換行。摘要依 metadata / TSV/index / OMIM signature 使用 process LRU + per-sample `case_summary.json` 持久化快取；SNV 摘要優先用 `snv_gene_index.sqlite` 依 variant id seek，CNV/SV 摘要只讀目標 id，避免 DRAGEN/WGS 大型 TSV 在首頁或重啟後被完整重掃。
- **Gene 搜尋**：SNV/Indel 與 CNV/SV 標題列右側各有 gene 搜尋框；SNV 打 `GET /samples/{id}/snv-search?genes=...`，以 prefixed `snv_gene_index.sqlite` seek `03_acmg` raw rows，再合併 prefixed `snv_annotations.sqlite`，不受 review TSV 限制。搜尋 modal 自己保留 `gnomAD_G_AF < 0.01` toggle。若個案有有效的 `{LIS_ID}.litvar2_annotation.json` marker，gene search 會把尚無 LitVar2 狀態的結果一次批次送到本地 SQLite 查詢，不逐卡發 request；正／負結果依目前 LitVar2 DB fingerprint 快取到 `{LIS_ID}.litvar2_on_demand.sqlite`，DB 更新後舊 cache 自動失效。每張 SNV 卡片 LitVar2 標題旁另有 refresh icon，可強制以最新本地 DB 重查單點；舊個案也允許手動查詢，但未手動更新的其他點仍顯示 `NA (請重跑三級)`。舊 unprefixed sidecar 與缺 index 的 bounded-memory raw scan 均保留 fallback；legacy 未遷移個案仍可讀舊 `snv_indel.annotated.tsv`。
- **IGV.js modal**：SNV/Indel、Mitochondria 與 CNV/SV 卡片的 ext-links 最左側有 `IGV`；基本資料性別下方另有「☑ 已確認」核取方塊與「確認 SRY」，以 hg19/hg38 對應 gene region 開啟相同 modal。SRY 模式初始只帶當前 sample，仍可手動加入 sibling，alignment coverage 固定 `autoscale:false,max:100`；勾選狀態存 `reports.sry_confirmed`。一般 modal 點擊後標題顯示 sample、padded locus、variant HGVS/類型與原始 `[build]` 座標；alignment track 預設 `displayMode: "SQUISHED"`、`visibilityWindow: 5000000`，SNV/Indel 與 Mito track height 固定 300，CNV/SV track height 固定 50。SNV/Indel 與 Mito 顯示前後 100bp，CNV/SV 原則上顯示原區間前後各 20% flanking area；若事件本身 ≤5 Mb，padding 會縮小以保持初始視窗 ≤5 Mb，直接載入 BAM coverage 而不先顯示 zoom-in 提示。browser-level `roi` 仍標出未加 padding 的實際 CNV/SV 區間；CNV/SV modal 內所有 BAM alignment coverage tracks 使用同一個 `autoscaleGroup`，讓 y-axis data range 一致以便比較 deletion/duplication；染色體先經 `_normalizeChrom()` 轉成 `chrN`。modal 先列 primary BAM + 同 batch 最多 2 個 sibling，可增刪 track，確認後才按「載入 IGV」。三級分析 job state 與 `pipeline_source.json` 保存原始 `source_sample_id` 與 `pipeline_type`；`GET /api/igv/bams` 會依來源分流，sample ID `-dragen` suffix 強制 DRAGEN roots、`-nckuh` / legacy `-inhouse` 強制 in-house roots，無 suffix hint 才看 sidecar。DRAGEN 找 `NGS_UI_DRAGEN_BAM_ROOT/<run>/bam/{sample}.bam`，預設 roots 包含 `/home/datalake_Raw/Novaseq` 與 `/home/datalake_Intermediate/n102968`，並排除 `{sample}.repeats.bam`；in-house 找 `NGS_UI_INHOUSE_BAM_ROOT/<batch>/{sample}/02_alignment/`。若 reviewer 自訂輸出 sample ID，會用 sidecar 回查原始 BAM ID；舊資料則移除已知 suffix 後 fallback 到原始 sample ID，並可容忍 in-house alignment 目錄中唯一的 `*.aligned.sorted.bam` 檔名不完全等於 sample ID。IGV modal 另有「其他路徑」下拉，呼叫 `/api/igv/bam-folders` 列出 DRAGEN `<run>/bam/` 與 in-house `<batch>`，選資料夾後再用 `/api/igv/bam-folder?dir=...` 列出 BAM（DRAGEN 取直下一層、in-house batch 取各 sample 的 `02_alignment`，排除 `.repeats.bam`）；按「+ 加入」會把選定 BAM 附加在現有 tracks 後方供比較，不取代原個案或同 batch 清單；相同路徑不可重複加入。`routers/igv.py` 的 `/bams`、`/batch-samples`、`/bam-folders`、`/bam-folder`、`/genome`、`/file` 全部需登入；`/file` 支援 HTTP range，且只允許 BAM roots 下的 BAM/BAI/CRAM/CRAI 與 `NGS_UI_IGV_REF_DIR` 下的 reference 檔。本機 hg38 FASTA + `.fai` 避免 igv.js 連 hospital intranet 擋住的 AWS S3。
- **載入效能**：`snv_review.ensure_review_tsv()` 自動維護 `snv_indel.review.tsv`，SNV/Indel 預載與搜尋都排除 chrM/MT，WES 先直接濾 DP<20，再讓主畫面只預載 `NGS_UI_CDS_CANDIDATE_BED` 內且 `GNOMAD_G_AF < 0.01` 或 AF 缺值的點；baseline ClinVar P/LP 及 weekly ClinVar `UP_TO_PLP` / `DOWN_FROM_PLP` rescue 不受 BED/AF 限制但仍受 WES DP≥20 限制。這個權威判斷是 `snv_review.is_review_retained()`，review materialization、SpliceAI、GeneBe live API fallback、LitVar2 都必須共用。review TSV 產生末段對 GRCh38 primary-contig SNV 以固定 BGZF/TBI tabix 補 `GPN_MSA_SCORE`（所有 case），只有明確 Research-only 才以本機 Parquet/Polars 另補 `GPN_STAR_SCORE`（M，相容舊欄位）及 `GPN_STAR_V_SCORE`（V）。兩者依 1-based `(chrom,pos,ref,alt)` 精確 join，不寫 raw／overlay，raw gene search 不補 GPN。Star 不對 indel/MNV/ref=alt 推分，缺值不設為 0。`dragen_run.py` 將 checkbox 傳成 `run_stopgaps.sh --research-only` → `build_snv_review_tsv.py --research-only`；CLI 未帶旗標記 false。manifest 的 `gpn_star_research_only` 是持久 opt-in，UI只依明確true顯示／刷新Star，v9.76恢復資料與套件signature變動時更新；舊 manifest 沒旗標不自動啟用 Star。關閉時不查 Star 資源並記 `skipped_not_research_only`，v9.75起UI亦隱藏殘留Star分數；完整重建只產生當次模式允許的註解。MSA／Star-M／Star-V 各有獨立 signature、`*_annotation` 狀態與 `*_values` 成功來源，任一缺資料／查詢失敗不中斷另一個或三級。Star 缺 manifest／shard／Polars 分別記 `skipped_missing_db`／`skipped_missing_shard`／`skipped_missing_polars`，其他失敗記 `failed`；失敗 refresh 保留成功分數來源。明確執行補建且只有單一GPN資源變更時，只補該分數，不重掃raw或重查另一套；review mtime 使 SNV LRU 自動失效。reviewer 已標記但不在 review TSV 的 SNV 由 `snv_gene_index.sqlite` 依 variant id 補入，同樣套 WES DP≥20，不會讓 review TSV 因 status 變動而重建；舊樣本缺少或尚未更新 index 時，會對少量已標記 id 做 raw TSV fallback scan，讓 gene search 標入 1/2/C 的 variant 在報告區仍能顯示。三級分析 post-processing 階段先執行 `annotate_mane_refseq.py`、`build_snv_review_tsv.py` 與 `build_snv_gene_index.py`；載入時 lazy rebuild 僅作為舊樣本 fallback。Exomiser/LIRICAL 的 `vcf_from_tsv.vcf.gz` 由完整 TSV 產生，但 v3.5 多 transcript rows 會先依 `(CHROM,POS,REF,ALT)` 去重；`vcf_from_tsv.vcf.gz.source.json` 記錄 writer version + TSV signature，讓舊未去重 VCF 自動重建。`sample_loader` 的 process-local bounded LRU 依 TSV / analysis / pheno / Exomiser / LIRICAL / OMIM signature 自動失效；完整 raw SNV TSV gene search 與個案清單單列摘要優先走 `snv_gene_index.sqlite`，不在載入個案後自動預熱。`/api/samples` 不再計算 case-list summary；`/api/samples/case-summary` 優先讀 `tertiary_output/_case_table.json`，只有總表缺失或 sample 集合不一致才重建，單列摘要仍以 `case_summary.json` 跨重啟保留；PDF gene list 延後到 `/api/samples/{id}/report-gene-list`，避免樣本初載同步解析 phenotype gene map。刪除個案只 invalidate 該 sample 的 SNV 與個案摘要 cache，不同步重掃全部清單。大型 JSON response 用 gzip。OMIM enrichment 每批只檢查一次 xlsx mtime。FastAPI startup 用 daemon thread 背景預熱 HPO、phenotype gene map、OMIM 與 mito ClinVar，避免大型 `phenotype_to_genes.txt` 擋住 HTTP port。
- **三級分析 modal / job log**：in-house 與 DRAGEN VCF 各自使用 typeahead，支援同來源批次；「加入同目錄全部檢體」位於「加入批次」左側，以索引 `batch_dir` 將 DRAGEN 同 run/input 目錄或 NCKUH 同一上層 batch 目錄的檢體一次加入，預設輸出 ID 為 `{source}-dragen` / `{source}-nckuh`，並沿用 mode 與 sample/source/path 去重；舊索引缺 `batch_dir` 時由前端既有路徑推導。VCF discovery index 由 `deploy/systemd/ngs-ui-tertiary-index-update.timer` 每天台北時間 02:00 自動刷新，人工更新與 24 小時 stale fallback 保留，`refresh_index()` 以跨程序 file lock 防止 timer 與 UI 同時寫入。worker 保存 `source_sample_id`，每次都用 v3.x samplesheet執行 Nextflow `-resume`，正式 output 是否已有檔案只寫 log，不再作 skip gate。DRAGEN 與 NCKUH 分別固定使用 `NGS_UI_TERTIARY_NF_WORK_ROOT/contexts/shared-{dragen,nckuh}/launch` 的共享 lineage；新 lineage 預設用 sibling `work/`，承接舊 session 時則以 `work-dir.txt` 固定沿用該 session 的原 work root，確保 cache metadata 與 task outputs 同時存在。因此 batch 組成改變時，相同 sample/input 的 task 仍可沿用 cache，input path/mtime/size、script、container 等 hash 條件改變時由 Nextflow 自行重跑。同一模式用 `TERTIARY_JOBS_DIR/.nextflow_locks/{dragen,nckuh}.lock` blocking `flock` 保護共享 LevelDB；等待時 step 為 `waiting-nextflow-cache`（2%），DRAGEN/NCKUH 彼此可並行。每次 Nextflow 成功、失敗或終止後，worker 都會在釋放 shared cache lock 前，把該次 shared launch 的 `.nextflow.log` 原子快照到 `TERTIARY_JOBS_DIR/{job_id}/.nextflow.log`，避免下一批輪替後下載到錯誤 run；三級分析清單的「查看 Log」讀完整 `log.txt`（只有執行中 modal polling 仍取 tail），旁邊的「下載 .nextflow.log」回傳該 job 快照，正在 Nextflow 的 job 可讀目前 live log，舊 job 沒有快照時不 fallback 到共享目錄。共享 lineage 第一次建立時，會掃 `REPO_ROOT/.nextflow` 與舊 `contexts/*/launch/.nextflow`，選本次 sample overlap 最大、再依 OK/ERR 與時間排序的一個 session 承接；不同 session ID cache 不合併。Nextflow `--out_dir` 固定指到正式 root 下 `.staging/{job_id}`，完整驗證 00-06 + SNV/Mito/STR/CNV/SV（勾選 PGx 時另驗證 07/PGx）後才在 staging 執行 post-processing，建立 prefixed sparse overlay/review/index、MITOMAP/ploidy 衍生檔與 marker。最後 `promote-output` 以同 filesystem rename 切換正式 00-07 與 worker-owned 08 檔，保留 live reviewer metadata/case summary/analyses，marker 最後切換；失敗會 rollback，原 live generation 不變。取消 PGx 的成功重跑會移除舊 07，避免 stale PGx。新樣本 marker 完成前不列入「載入新個案」。進度為 Nextflow 3-82%、prepare/post-processing/overlay/review/index 82-99%、promotion 99%；清單忽略 `.staging`/`.rollback` 並合併新 root、legacy pipeline root 與 job state。完整刪除才會刪 pipeline/UI tree 與該 job 的兩種 log；個案清單刪除只 unregister。
- **三級分析 output 來源分流 / promotion**：正式統一 root 以 UI sample ID（如 `{source}-dragen` / `{source}-nckuh`）管理；staging 以 samplesheet 的原始 source ID 接收 Nextflow 產物，驗證後 promotion 到 UI ID。`-WES` / `-WGS` 是 test-type ID，不可當 caller alias 剝除。legacy `/home/pipeline/tertiary_output` 已移除；程式保留的 fallback 只在目錄存在時讀取，不參與 cache／promotion。
- **三級分析 timing log**：worker-owned log timestamp 改成後綴 `[YYYY-MM-DD HH:MM:SS]`，不含 `T` / `+08:00`。Nextflow 以 `-ansi-log false` 執行，job log 是每個 task 一行的純文字（不再是數萬行 ANSI 重繪畫面）；worker 另寫 `[nextflow-step] {PROCESS} start` / `done n/total elapsed=…m`，state 內 `nextflow_step_history` 記錄 start/done、`done/total` 與 elapsed。post-processing summary log 只把單獨的 `skipped` 狀態字當成 skipped，GeneBe 的 `'*'-skipped 0` 這類計數不再讓 genebe 被誤標成 skipped。
- **三級 post-processing**：整批逐步驟執行（`scripts/postprocess_batch.py`，`run_stopgaps.sh` 為單一 sample wrapper），每個資源對整批只查一次、各 sample 平行寫回；詳見 `AGENTS.md` §7／§11。
- **三級分析進度 / 取消**：時間加權進度條（`backend/app/services/tertiary_progress.py`）：Nextflow 用 `-ansi-log false -with-trace`，依 task 完成與每個 process 的預估秒數推進；Nextflow / post-processing / promote 依「固定成本 + 每例成本 × n」模型分配比例，profile 存 `data/jobs/tertiary/_timing_profile.json`；前端顯示 `progress_pct`、預估剩餘時間，Nextflow 階段另顯示 console 風格即時狀態區塊（`nextflow_processes`）。詳見 `AGENTS.md`。進度面板「終止」按鈕會 `POST /api/dragen/jobs/{job_id}/cancel`，後端向 worker process group 送 SIGTERM 並標記 `cancelled`。
- **DOCX CNV/SV 表格**：「變異位置」欄使用 buffered wrap，內容寬度比欄寬少一格，避免座標字串貼到後面的拷貝數欄。
- Mito TSV 複製完成後會呼叫 `mitomap_mito.annotate_mito_tsv()` 靜默回補 `MITOMAP_*` 欄位；這不是進度 step，不會在 job log 額外顯示。回補欄位只給 `mito_tsv.py` 做 tier 分類，前端卡片不呈現 MITOMAP。
- **共用 class**：`.variant-head`（`#index` span + `.status-select` + ...）、`.btn-copy`（`COPY_ICON_SVG`）、`.ext-links`（Varsome 那種按鈕樣式）、`.cnv-sv-detail-box`（灰底兩行 + 折疊區）、`.cnv-sv-reasoning`、`.cnv-sv-comment-text`、`.acmg-class` + `.sig-p/.sig-lp/.sig-vus/.sig-lb/.sig-b` 五級色。**`.modal-card input[type=text]` 那條 catch-all（width:100%）被 `.variant-card`/`.cnv-sv-card` 裡的 input 排除**（不然 gene-search modal 裡的 ACMG points 那 3em 格會被撐爆）。
- **CNV/SV 卡片**（`renderCnvSvCard`）：每張左側無色條、inset 在 tier-panel 上；header = `#index` + status select（1/2/C/0）+ `CNV/SV` tag + SV-type pill（DEL紅/DUP藍/INV橙/INS紫/TRA灰）+ 座標 + 複製鈕 + `{chromN}{cytoband}`（如 `12p11.21-q24.33`）+ ext-links（UCSC/DECIPHER/dbVar/GeneCards）靠右；detail box 第一行 = ACMG 五級下拉（`.cnv-sv-acmg-select` + `sig-*` 色；reviewer override 存 `state.reports.edits[id].ACMG_class_sv`，跟 SNV ACMG 分開）+「涵蓋基因數: 1518（疾病相關：28）」+ 基因型 + Filter + Qual，第二行 = `AnnotSV 評分依據`（折疊）+ Score；基因表預設只顯示 phenotype 相關（`in_panel=true`）基因（首欄 checkbox → `state.reports.edits[id].report_genes`；Phenotype 與 Inheritance cell 點擊展開 `.gene-clip-cell`），按小三角形才 lazy-render phenotype score 為 0 的其餘完整列 / compact chips；「已知致病區域重疊」「已知良性區域重疊」兩段（DEL→只 P_loss/B_loss、DUP→只 gain、其他→全部；內容 CSS `-webkit-line-clamp:2` + 「展開全部」鈕；無資料顯示提示而非整段消失）；最後 Disease textarea（→ `state.reports.edits[id].disease`，DOCX 與個案清單使用）與 Comment textarea（→ `state.reports.edits[id].comment`）。同來源、同染色體、同 DEL/DUP 且相鄰 gap ≤250 kb 的片段預設自動整合；copy number 差異不再阻擋視覺整合，原始 segments 仍可展開查看各自 CN。UI / DOCX / 個案清單改用 parent CNV/SV。parent 取合併片段中最高 combined score 的代表資料，分析區與報告區皆依 `max_pheno_score + scaled AnnotSV score` 排序；前端依 sample payload 快取整合 view，render tiers 與 cards 時不會反覆排序及重建 parent。
- **Mito 卡片**（`renderMitoCard`）：header = `#index` + status select + locus pill（protein紅/tRNA橙/rRNA藍/control灰）+ mitochondrial `m.HGVS`（v3.2 由 `POS/REF/ALT` 產生，不使用 transcript `HGVS_C` 的 `c.`）+ 複製鈕 + gene + heteroplasmy 徽章（teal）+ gnomAD-MT / ClinVar 連結；detail box 第一行 = ACMG 五級下拉（`ACMG_classification_mito` reviewer override，報告區同步重 render）+ `ref→alt` / 類型 / `Heteroplasmy (AD·DP)` / genotype，只有來源 TSV 真的有 `FILTER` 時才顯示 Filter；第二行 = Consequence + Impact + Biotype + Protein change（v3.2 `HGVS_P` 會 URL decode 並移除 protein transcript，只留報告用 `p.xxx`）+ ClinVar + gnomAD-mito AF + TLOD（缺值顯示 `—`）。MITOMAP 欄位只供後端 tier 分類，MITOMAP 區塊與 MITOMAP 連結不顯示；`CLINVAR_DN` 以 `&` 拆成 ClinVar disease checkbox，勾選者才進 DOCX。
- **OMIM disease list**（SNV 卡片，`renderDiseaseList`）：`<details>` 列出完整 `Disease1..16`（跳過 NA），summary 有報告勾選 checkbox（→ `state.reports.edits[id].report_disease_picks`：`{phenotypeMIM}|{normalized name}` key + label 快照，OMIM 改名／調序不會換錯疾病；舊 `report_diseases` slot 只作 fallback，`scripts/migrate_report_disease_picks.py` 一次性轉換），缺說明的 `*` 連到 OMIM entry 頁，展開的黃底框（`.disease-detail`）底部有「▴ 收合」鈕。只有「個案清單」摘要會呼叫 `omim_store.compact_disease_label()`，把 curator 長篇說明截到第一個可辨識的遺傳模式括號（如 `(AD)`、`(AR)`）結束；variant 卡片、原始 `OMIM.xlsx` 與 `OMIM_disease` 都保留完整內容。

---

## 6. pheno_score 自動寫入

`analyses_store.write_version()` 寫完 `analysis.json` 後 side-effect 算 `compute_pheno_score()` 並寫 `pheno_score.tsv`（HPO/panels 為空就刪掉舊檔）。所以 `register` 新個案、編輯 phenotype（`routers/phenotype.py`）、複製/重命名 version（`routers/analyses.py`）都會自動產生 pheno_score.tsv，不用等「▶ 開始分析」。`sample_loader` 還有 lazy backfill：載入時 pheno_score.tsv 缺失或比 analysis.json 舊就即時重算。SNV 的 `in_panel` 不再寫回 raw TSV 的 `IN_PANEL` 欄；主畫面載入與 gene search 都用當前 active analysis 的 `pheno_score.tsv` 動態補 `in_panel` / `pheno_score` / `total_score`。`write_pheno_table(sample_id, scores, target_dir=...)` 支援指定 version 目錄（不一定是 active 的）。

---

## 7. 轉換 / annotation scripts（`scripts/`）

| script | 用途 |
|---|---|
| `convert_anno_combined_to_tertiary_tsv.py` | 舊 R pipeline 的 `anno_combined.txt.gz` → `snv_indel.annotated.tsv`（去重 by 變異留最佳 transcript：MANE_SELECT > MANE_PLUS_CLINICAL > CANONICAL > any；缺的欄留空）。用法：`--in <file> --out <輸出目錄>/snv_indel.annotated.tsv`（舊 layout 轉檔用；新 layout 的完整 SNV 來源是 `/home/datalake_Intermediate/pipeline/tertiary_output/{LIS_ID}/03_acmg/`） |
| `convert_old_json_to_tertiary_tsv.py` | 舊 webdata JSON → `snv_indel.annotated.tsv` |
| `filter_snv_tsv.py` | GUI TSV 的結構性清理：只移除 `REF/ALT=*` 與非 primary contig；AF、VAF、IMPACT 留在 TSV，由 UI 顯示 filter 控制。v3.x 正式 pipeline 自行處理 NCKUH/DRAGEN 前處理（舊 staging 流程已移除） |
| `build_snv_review_tsv.py` / `services/gpn_msa.py` / `services/gpn_star.py` | 從 03 raw + optional `--overlay` 預建 08 review TSV；統一套 `is_review_retained()` 後，以 `NGS_UI_GPN_MSA_DB` + tabix 補 MSA，`--research-only` 時再用 `NGS_UI_GPN_STAR_DIR` + Polars 補 Star。正式 chain 一律 best-effort；缺資源或 query failure 記 manifest status + warning 而不中止，`--require-gpn-msa`／`--require-gpn-star` 只給人工 strict 診斷，Star 仍須 `--research-only`。manifest 的 `test_type` 用 `test_types.normalize_test_type` 記錄成後端 loader 會要求的值（年份+T 如 `26T…` 一律 `TITAN-WGS`，即使 worker 傳 `--test-type WGS`；篩選行為與 WGS 相同），否則 UI 第一次開個案會因 manifest 不符整份重建；`run_stopgaps.sh` 傳 `--sample "$SID"`，因為 staging 時 `08_postprocessing` 的上層不是 sample ID |
| `build_snv_gene_index.py` | 從 03 raw 預建 `--out 08_postprocessing/{LIS_ID}.snv_gene_index.sqlite`，保存 gene/variant id 與 raw byte offset |
| `annotate_giab_strata.py` | 在暫存 working TSV 補 `GIAB_STRATA`，最後存入 sparse overlay；不寫回 03 raw。strata 由 `NGS_UI_GIAB_STRAT_DIR/strata_manifest.json` 定義 |
| `download_giab_strata.sh` | 在 dev 機跑一次：從 GIAB `@all` GRCh38 stratification tarball（預設 v3.6 `genome-stratifications-GRCh38@all.tar.gz`）下載/解壓，再用**穩定 substring glob**（`*AllHomopolymers*`、`*AllTandemRepeats*`、`*segdups*`、`*lowmappabilityall*`、`*gclt*orgt*_slop50*`、`*allOtherDifficult*`；多筆取最短檔名＝union）挑出 homopolymer / tandem repeat / segdup / low mappability / GC extreme / other difficult 六個子集到 `NGS_UI_GIAB_STRAT_DIR`，並依實際解出的檔名寫 `strata_manifest.json`。pattern-based 所以跨 GIAB 版本（v3.1/v3.6…檔名 slop/threshold 後綴會變）都適用。BED 大、不進 git。用法 `scripts/download_giab_strata.sh [--dir <out>] [--tarball <local.tar.gz>] [--url <tarball-url>] [--keep] [--force]`，可用 `GIAB_STRAT_TARBALL_URL` 覆寫來源。 |
| `backfill_giab_strata.sh` | 對「GIAB 步驟之前就分析完」的舊樣本補 `GIAB_STRATA`：逐 sample 跑 annotate → 重建 review TSV → **重建 gene index**（annotate 原子改寫整個 TSV 會位移所有 byte offset，而 gene index 只在缺檔時才自動重建，故必須一起重建，否則 gene 搜尋會 seek 到錯位 bytes）。idempotent；不帶參數＝處理 `TERTIARY_OUTPUT_ROOT` 下所有樣本，帶 sample ID 則只處理那幾個。新分析不需要（`run_stopgaps.sh` 已內含 annotate 步驟）。 |
| `annotate_acmg_genebe.py` | GeneBe ACMG 第二意見：主 DB 查整張 TSV，主 DB/API cache 都 miss 且符合 review filter 的點才用 pygenebe Apptainer live API。API 成功結果寫中央 cache 與精確 7 欄、去重的 `api_pending/*.tsv`；no-result 預設 cache 30 天。缺 credentials/SIF、timeout 或 batch 失敗只警告，不中止三級分析；主 DB 與 pipeline ACMG 仍可用。 |
| `deploy_genebe_db.sh` | 一鍵部署 GeneBe DB：驗 bgzip 完整性 → 濾掉 pos 非整數的壞行（`--sort` 可一併排序）→ bgzip → `tabix -s1 -b2 -e2 -c '#'` 重建索引 → 探測 → **原子換檔**安裝到 `NGS_UI_GENEBE_DB`。給 daily rebuild 後一行部署、避免手動 clean+reindex / stale index。NGS-UI 本身串流讀 `.gz`、不靠 `.tbi`，但索引仍給 tabix 類工具用。用法 `scripts/deploy_genebe_db.sh [SRC.tsv.gz] [--dest DEST] [--sort]`。DB 建置端要求見 `docs/genebe_db_requirements.md`。 |
| `compare_genebe_spliceai_coverage.py` | 診斷 GeneBe 本機資料庫的 SpliceAI 完整性：從 `genebe_hg38.tsv.gz` 抽樣 variant，使用目前 extra-VEP 的 VEP SpliceAI plugin 路徑重算 `SPLICEAI_MAX`，輸出 coverage summary、mismatch 明細與 run metadata，供評估是否能用本機 GeneBe DB 取代 extra-VEP / GeneBe API；正式估計建議 `--max-sites 100000 --sample-mode chrom-balanced`。 |
| `import_fixed_panels.py` | `fixed_panel/WES-I.xlsx`、`WES-II.xlsx` 與 `other_panel/` → `phenotype_data/fixed_panels/index.json` + `gene_panels/*.txt`。WES Excel 只讀 `gene panel list` 標記列起始的基因區塊，不可把疾病名或資料來源列當成基因；顯示名稱只移除科別/級別前綴，需保留 `Non-syndromic` 這類 panel 名內的連字號。持久化 key 為相容既有 analysis metadata 仍沿用舊規則，和顯示名稱分開。更新來源 Excel 後要重新執行匯入。 |
| `annotate_dragen_mito_vcf.sh` | DRAGEN hard-filtered VCF 內含 chrM calls 時的 legacy mito wrapper。直接把 `{sample}.hard-filtered.vcf.gz` 交給 `parse_mito_vcf.py`，只取 `chrM/MT/chrMT` rows，`FORMAT/AF` 當 heteroplasmy、`FORMAT/SQ` 在缺 `INFO/TLOD` 時填入 TLOD 欄，輸出 `<outdir>/mito.annotated.tsv`。用法：`scripts/annotate_dragen_mito_vcf.sh --in /path/{sample}.hard-filtered.vcf.gz --sample {sample} --outdir ~/legacy_mito/{sample}`（獨立輸出，不會被 UI 讀取）。 |
| `annotate_mito_vcf.sh` + `parse_mito_vcf.py` | Legacy/手動補跑 MITOMAP-only mito 轉檔。正式流程直接讀 `04_mito/{sample}.mito.tsv`；只有 MITOMAP 確實補值時 worker 才保留 `08_postprocessing/mito.annotated.tsv`。 |
| `migrate_to_versioned_layout.py` / `migrate_vcf_path.py` / `rewrite_vcf_paths.py` | 一次性的舊→新佈局遷移 |
| `migrate_report_disease_picks.py` | 一次性把舊 SNV 疾病勾選（slot 編號）轉成 `report_disease_picks`；預設 dry-run，`--apply` 寫入 |
| `omim_db.py` | 平台 OMIM 資料庫：`init`（固定舊勾選後從 OMIM.xlsx 建庫）、`status`、`reviews`、`export --out` |
| `probe_emr_api.py` | NCKU EMR API 診斷（urllib，dump 到 /tmp/emr_probe/；用內建 14 筆 MRN 清單，不讀 argv） |
| `docs/ops/pangolin_cache/` | **正式 pipeline 用**（不在 NGS-UI 內執行，部署到 `/home/pipeline/tertiary_code/`）：`pangolin_cache.py` 跨 sample 的 Pangolin 分數 SQLite 快取（split → 只送沒算過的變異給 GPU → merge 依原順序組回；依容器/DB/參考序列/參數分區；`seed` 匯入既有 `02_pangolin` 輸出並跳過衝突值、`compare` 驗證、`stats`）+ 替換用 `PANGOLIN_SCORE.nf`（以正式版為底的最小改動；`params.pangolin_cache_db` 空字串＝舊行為；同時把 `params.containsKey` 改成直接讀屬性）+ `apply_to_pipeline.py`（預設只預覽 diff，`--install` 才寫入並留 `*.bak-時間` 備份）+ README 部署/驗證步驟。cohort 重複率 98.4%，WES 每例 Pangolin 由 ~8 分鐘降到 <1 分鐘。**2026-11-21 已部署到 DGM 正式 pipeline**：快取在 `/home/pipeline/reference/hg38/tertiary/pangolin/cache/pangolin_scores.sqlite`（`seed --since 2026-08-03` 匯入 147,916 筆；08-03 前的輸出有 GPU 浮點差異造成的衝突，已排除），安裝備份為 `modules/snv_annotation.nf.bak-*` 與 `nextflow_tertiary.config.bak-*`；24WE 批次實測每例命中 99%+、只算 13–31 個新變異，與舊輸出 `compare` 完全一致。測試：`tests/test_pangolin_cache.py` |

`mito.annotated.tsv` adapter 同時支援舊欄位（22）：`CHROM POS REF ALT HGVS_M GENE LOCUS_TYPE CONSEQUENCE AA_CHANGE HETEROPLASMY AD DEPTH FILTER TLOD MITOMAP_DISEASE MITOMAP_STATUS MITOMAP_PLASMY MITOMAP_GB_FREQ MITOMAP_GB_SEQS MITOMAP_REFS MITOTIP_SCORE MITOMAP_ALLELE`，以及 v3.2 pipeline mito 欄位：`CHROM POS REF ALT GENE HGVS_C HGVS_P CONSEQUENCE IMPACT BIOTYPE GENOTYPE DP AF_SAMPLE GNOMAD_MITO_AF CLINVAR_SIG CLINVAR_DN GNOMAD_MITO_AF_HOM GNOMAD_MITO_AF_HET GNOMAD_MITO_AN CLINVAR_VARIATION_ID OMIM_IDS PIPELINE`。`AF_SAMPLE→heteroplasmy`、`DP→depth`、`HGVS_P→aa_change`（清成 `p.xxx`），`HGVS_C` 是 transcript `c.`，不拿來當卡片標題；adapter 會用 `POS/REF/ALT` 產生 `m.{pos}{ref}>{alt}`。若 TSV 含 `FILTER`，adapter 會照舊先保留 PASS；若沒有 `FILTER`，不做後端 FILTER 篩選，UI 也不顯示 Filter 欄。v3.2 TSV 複製後可由 `mitomap_mito.annotate_mito_tsv()` 精確 `(POS,REF,ALT)` 回補 MITOMAP 欄位；UI/報告卡片不呈現 MITOMAP，但 tier 1 會使用 MITOMAP confirmed/pathogenic 狀態。

`mito.vcf.gz` 是 GATK Mutect2 `--mitochondria-mode` 輸出（FilterMutectCalls + 黑名單 mask；chrM = GRCh38 chrM = rCRS/NC_012920.1；`FORMAT/AF` = heteroplasmy fraction，`FORMAT/DP` = depth，`FORMAT/AD`）。

完整 SNV source 是 `03_acmg/{SAMPLE_ID}.snv_indel.acmg.tsv`。v3.4 以前含 `MANE_ALL`（65 欄）；v3.5 移除 `MANE_ALL`、新增 `STRAND_BIAS`，仍為 65 欄且同一 genomic variant 可有多 transcript rows。NGS-UI post-processing 補的 `REFSEQ_*`、`GENEBE_*`、GIAB、院內 AF 等只存 sparse overlay，不寫回 03；adapter 合併 raw+overlay 後再依 `CHROM-POS-REF-ALT` 合卡。
CNV/SV 是 AnnotSV 標準輸出（128 欄；`Annotation_mode` full=一個 SV 一列、split=每 gene 一列；adapter `annotsv_tsv.py` 用 index-based 解析只取 ~30 欄、聚合 full+split）。

`tertiary_code/` 裡有次級 pipeline 的 Nextflow（`main_tertiary.nf`、`modules/`、`scripts/parse_vep_csq.py` 等）+ config（`dgm` profile = NCKU 正式環境 `/home/pipeline/reference/hg38`、`/home/pipeline/nextflow_containers`、`--bind /home`；`local` profile = `/scratch/pylin1991/...`、`/data/pylin1991/nf-containers`）—— 純參考用。UI worker 使用正式環境 `/home/pipeline/tertiary_code/{main_tertiary.nf,nextflow_tertiary.config,scripts}`（config 可用 `NGS_UI_TERTIARY_CONFIG` 覆寫）；PGx checkbox 預設勾選，取消時 Nextflow 加 `--run_pgx false`；repo 內的 workflow/scripts 版本可能落後，不可整包覆寫正式 `scripts_dir`。

---

## 8. 輸入臨床表徵工具（`/phenotype/`）

`frontend/phenotype/`（從舊的 GitHub-backed hpo-docs 改寫，砍掉 GitHub OAuth/terminal/run-analysis）。**不需登入**，由 NGS-UI 伺服器靜態服務在 `/phenotype/`（`main.py` 加 `GET /phenotype` → 307 redirect `/phenotype/`，註冊在 StaticFiles mount 之前）。功能：HPO term 搜尋（Fuse.js + 本地 `hpo_data.json` 3.5MB，**在 repo 裡** `frontend/phenotype/`）；Gene Panels 搜尋（打 `GET /api/phenotype-tool/panels`）；**Custom panel**（名稱 + 基因清單 textarea + weight；按「產生 phenotype.txt」時 POST `/api/phenotype-tool/custom-panel` 建檔到 `gene_panels/`、即時更新 `phenotype_scorer` 記憶體、名稱自動清理成 `[A-Za-z0-9_-]{1,64}`、衝突 409、基因**不大寫**（`C7orf50` 保留小寫）、case-sensitive 去重）。「產生 phenotype.txt」一鍵：建 custom panel → 組 TSV → POST `/api/phenotype-tool/save` 寫到 `patient_phenotype/`。MRN 或 LIS_ID 至少填一個；檔名：兩個都填 `{code}_{mrn}_phenotype.txt`、只 LIS_ID `{code}_phenotype.txt`、只 MRN `{mrn}_phenotype.txt`。「載入既有資料」用 `GET /api/phenotype-tool/load?code=&mrn=`。phenotype.txt 格式：`phenotype\thpo_name\tweight` 表頭 + `HP:xxxxxxx\t<name>\t<weight>` / `<panel_name>\t\t<weight>` 列（`phenotype_io.parse` 讀這個）。

主畫面的 Patient phenotype card 與「載入新個案」modal 都有 `WES-I / WES-II / WGS / Other panel` tabs，預設展開 `Other panel`；固定 panel chip 由 `GET /api/phenotype-tool/fixed-panels` 載入，勾選後以 weight=1 寫入各自的 panels working copy，Other panel typeahead 則排除固定 panel key，避免同一個 panel 出現兩個入口。Patient phenotype、新個案 modal 與 `/phenotype/` 工具的 HPO / panel 下拉都支援鍵盤上下鍵與 Enter 選取；dropdown 開啟時 Enter 只會選項目，不會送出載入個案或開始分析。Patient phenotype 和 Comment 中間的 Dead zone card 依目前 HPO + panel gene set 顯示 cohort dead exons；WES 用 20X，WGS（in-house/DRAGEN 都一樣）用 DRAGEN 15X；主畫面預設只顯示前 10 列，可用小三角形展開全部。平台啟動讀取索引、載入個案核心資料或登錄新個案期間會顯示不可用背景點擊或 ESC 誤關閉的「資料載入中」遮罩；Exomiser/LIRICAL job 完成後的 sample refresh 使用 silent reload，不再用全頁遮罩阻擋判讀。

`routers/phenotype_tool.py`：`GET /api/phenotype-tool/panels`、`POST /api/phenotype-tool/save`、`GET /api/phenotype-tool/load`、`POST /api/phenotype-tool/custom-panel` —— **全公開無 auth**（intranet 信任 + 嚴格驗證：token 限 `[A-Za-z0-9_-]{1,32}`、檔名從驗證過的 token 拼、內容 ≤64KB、panel 基因 ≤5000）。

---

## 9. 上傳個案清單（roster）

`patient_list_store.py`：上傳 NCKU「未完成報告清單」xlsx（`POST /api/patient_list`，**需登入**），原始檔存到 `patient_list/{ts}_{name}.xlsx`，merge 進 `patient_list/roster.json`（**additive**，不刪舊的）。前端 file input 支援一次多選，必須逐檔依序 POST（不可平行改寫 `roster.json`）；modal 只顯示當次批次每個檔案的成功／失敗與統計，歷次紀錄用預設收合的 `<details>`，展開時才讀 `/api/patient_list/uploads`。解析器逐一掃描所有工作表，找任一列／欄內含 `檢體編號` 的標題列；`檢驗名稱` 與新版表單的 `檢驗項目` 視為同欄，跨工作表 by LIS_ID 去重且第一筆優先。砍 `8BB1` 前綴得 LIS_ID（`8BB126WE0092`→`26WE0092`），檢驗名稱／項目判斷 WES/WGS，其他可讀欄位為 `病歷號|姓名|科別|開單醫師|簽收時間`，缺少的 optional 欄位留空。`unregistered_samples.listing()` / `detail()` 用 roster 自動填「載入新個案」modal 的 MRN/姓名/Test type（科別只當提示文字）；若三級分析 UI sample ID 帶 `{LIS_ID}-dragen` / `{LIS_ID}-nckuh` / `{LIS_ID}-inhouse` 後綴，roster 與 phenotype lookup 會先試完整 sample ID，再試去掉 caller 後綴的 LIS_ID，並可用 `pipeline_source.json.source_sample_id` 作為 fallback；若 `pipeline_source.json` 顯示 DRAGEN 來源，或 in-house source VCF 大於 100 MB，前端會把 Test type 預設為 WGS 並覆蓋 roster 預設。phenotype 檔查找順序固定為 `{roster_mrn}_phenotype.txt` → exact `{candidate_lis_id}_{roster_mrn}_phenotype.txt` → `{candidate_lis_id}_phenotype.txt` → `{candidate_lis_id}_*_phenotype.txt`（glob）；MRN-only 即使是 header-only 也為權威，不再落到舊資料。Clinical presentation sidecar 同樣先用 MRN-only `{roster_mrn}_clinical_presentation.txt`，再 fallback 舊 LIS/MRN/candidate suffix。`GET /api/patient_list` = 看目前 roster（debug）。

> **`_index.json` 不要拿來放 roster** —— 它是「已登錄樣本清單快取」，`list_index()` 每次都重寫。roster 用獨立檔。

---

## 10. 其他

- **認證**：SQLite `data/users.db` + bcrypt（`users.py`）。建帳號：`PYTHONPATH=backend python -m app create-user [username]`（從 `backend/` 的 parent dir 跑；不用重啟服務）。`PYTHONPATH=backend python -m app list-users`。8h session cookie（SameSite=Lax）。沒有改密碼/刪帳號的指令。
- **OMIM annotation**：`omim_store.py` 啟動時讀 `OMIM.xlsx`（`_warm_caches` 預載；mtime 變了自動 reload；找不到檔就靜默 disable）。`sample_loader` 每個 SNV 變異 join `Disease1..16`/`OMIM_id`/`OMIM_disease`/`Inheritance`（OMIM_LINK 解析出的 OMIM_id 優先、gene_symbol fallback）。OMIM.xlsx 欄位：`OMIM_id | gene_symbol | OMIM_disease | Inheritance | Disease1..16 | Done`（17824 列，`OMIM_disease` 多行文字、每行 `<病名> (繼承碼)`）。`Disease1..16` 全數缺失時，`omim_store` 會依序從 `OMIM_disease` 合成最多 16 個單行 fallback；已有任一 curator slot 時不補其餘空欄。
- **OMIM 擷取（PoC）**：`reference/omim_capture/`（Chrome/Edge MV3）在 `omim.org/entry/*` 顯示面板，按 Capture 才讀頁面並自動展開 Clinical Synopsis，背景 worker 以個人 capture token（`Authorization: Bearer ngsc_…`，平台 `/omim-capture/` 產生，只存 hash）送到 `POST /api/omim/captures`。`services/omim_capture.py` 解析頁面文字、存 `OMIM_CAPTURE_DIR/{MIM}/*.json.gz`、與 `omim_store` 工作簿比對（缺說明／OMIM Edit History 比 Done 新／已是最新）。平台 OMIM 資料庫（`services/omim_db.py`，`scripts/omim_db.py init` 建立）存在時 `omim_store` 改讀它：骨架＝xlsx `OMIM_disease` 各行、保留舊 Disease 欄原文、遷移配不上者進待確認；面板／`/omim-capture/` 可「採用到平台」（`POST /api/omim/captures/{mim}/{name}/adopt`），基因頁採用只同步骨架、不自動刪除；OMIM.xlsx 改由 `/api/omim/export.xlsx` 匯出。`services/omim_queue.py` 以基因為單位產生待補清單（先基因頁同步、再補缺說明／逾一年未更新的疾病；進行中→疾病相關基因的基因頁→其他；待確認併入各基因），`GET /api/omim/queue` 與擴充功能面板「前往」（`/api/omim/queue/next`，同一個基因做完才換下一個，排除已擷取未採用、7 天內採用與略過的頁面或基因）。`services/omim_lookup.py` 提供基因查詢（`GET /api/omim/genes/search`、`GET /api/omim/gene?q=` symbol／別名／gene MIM／phenotype MIM），顯示卡片實際使用的內容、來源、採用紀錄、待補與待確認。
- **自動儲存**：reviewer 編輯後 1.5s debounce → PUT `/samples/{id}/report`，包含 CNV/SV Disease 與 Comment textarea；開啟個案清單或匯出 DOCX 前會先 `flushPendingSave()`，避免 reviewer 剛輸入完成但 debounce 尚未觸發時漏掉最新內容。三個位置的「儲存」按鈕（top/mid/bottom，class `.js-btn-save`/`.js-save-hint`）；存成功後 hint 顯示 `已儲存（HH:MM:SS）`；`beforeunload` 在 dirty/inflight 時警告；`_lastSavedAt` 切樣本時重設 null。
- **EMR 整合**（`emr_client.py` + `routers/emr.py`）：NCKU intranet 兩支 API —— GetPhenotypeList（broken JSON，需修復）+ APIM easyform/getdata（X-IBM-Client-Id header）。reviewer txt > EMR 的 HPO 優先序；EMR 的 sex 覆寫 reviewer 打的；新 `genetic_counseling` 欄。`NGS_UI_EMR_CLIENT_ID` 空 = 整套關閉。sample-card 上有「🔗 EMR」連結 + 「EMR 同步」按鈕；載入新個案 modal 也有 EMR 同步（但 EMR 不回姓名，姓名只能手打）。
- **Exomiser/LIRICAL rerun worker**（`workers/exomiser_lirical.py`）：渲染 `phenotype_reference/exomiser_input.yml`/`lirical_input.yaml` 模板 → java -jar 跑 → 結果寫 `analyses/{ver}/exomiser_results.tsv`/`lirical_results.tsv`。**不算 pheno_score**（那是 `write_version` 的事）。RQ/Redis job queue。
- **`docs/`**：`ACMG_SF_v3.3.txt`、`carrier_mackenzie_1300+.txt`、`proactive.txt`（panel 文字檔）、舊版分析網頁的 `app.js/index.html/style.css`（參考用）。

---

## 11. 已知踩雷 / 慣例

- `REPO_ROOT / "phenotype_data"` 已改成 `config.PHENO_DATA_DIR`（= `NGS_UI_HOME/phenotype_data`，**無 fallback** —— dev 機部署時一定要 `mv ~/NGS_UI/NGS-UI/phenotype_data ~/NGS_UI/phenotype_data`，不然 HPO 搜尋空、pheno_score 全 0）。`hpo_ontology.py` 跟 `phenotype_scorer.py` 都 import 它。
- MITOMAP 兩個 TSV 是 **Latin-1**，不是 UTF-8（`0xa0` nbsp），loader 用 `encoding="latin-1"`。
- `parse_mito_vcf.py` 不做 POS-only MITOMAP fallback（不然 `m.114C>A` 會被配到 `m.114C>T` 的 disease，等等）。
- `.modal-card input[type=text]` catch-all 已排除 `.variant-card`/`.cnv-sv-card` 裡的 input。
- ACMG_CLASS 在 `snv_tsv._normalize_acmg_class` 正規化（`VUS`/`uncertain_significance`/各種大小寫 → `Uncertain significance`；認不出的留原樣 → UI 顯示 `—`）。
- SNV/Indel ACMG 顯示與輸出優先序：reviewer override → GeneBe `GENEBE_ACMG_CLASS` → pipeline `ACMG_CLASS`；criteria / score 亦同。
- GeneBe `GENEBE_ACMG_*` 採主 DB → API cache → review-filtered live API；主 DB 永遠優先。API 成功結果保存為中央 SQLite cache 與正式 7 欄 pending TSV，API 失敗不中斷分析。主 DB 查詢預設用 signature-aware lazy SQLite，失敗才 streaming；完整 working TSV 仍只暫存並由 worker `finally` 刪除。
- CNV/SV `genes` 陣列在 adapter 裡切到前 10（不然跨染色體的大 DEL 會塞 1500+ gene record 到 payload）；前端 overflow chip lazy-render（`<details>` 的 toggle 事件不 bubble，用 capture-phase listener）。
- `tier-tab` 的 click dispatch 用 `data-tier` 值判斷 SNV / CNV-SV / Mito 哪一組（`CNV_SV_TIER_ORDER`、`MITO_TIER_ORDER` 的 includes 檢查）。
- ACMG 五級色 class `.sig-p/.sig-lp/.sig-vus/.sig-lb/.sig-b` 是全域的；CNV/SV 的 ACMG 下拉用它們但**不要**加 `.acmg-class` class（會觸發 SNV 的 change handler 寫錯 state）。
- Mito tier 1/2/3 是「有 FILTER 則先 PASS 過濾，無 FILTER 則全收」；MITOMAP 只回補欄位與輔助 tier 分類，不做 disease-relevant 預篩，也不顯示在卡片。其他 PASS mtDNA variants 會進 `MITO-3`。
- Nextflow process 的 script 區塊裡**不要對 `params` 呼叫方法**（`params.containsKey(...)`、`params.get(...)`、`params[...]`）：Nextflow 會把整個 params 算進 task hash，而 NGS-UI 每個 job 的 `out_dir` / `samplesheet` 都不同 → 該 process 永遠無法 `-resume`（正式 PANGOLIN_SCORE 曾因此每次重跑）。直接讀 `params.xxx ?: default`。
- 別 commit 沒被 `.gitignore` 排除的患者資料 / 大檔。

---

## 12. TODO / 還沒做

- **候選診斷整體排序（表型 + 變異 + 遺傳模式）**：取代現行 `total_score = geno_score + pheno_score` 與 `phenotype_scorer` 的有或無查表。研究與實作計畫見 `docs/research/NGS_Integrated_Prioritization_Plan_v3.md`（第 4 章為 HPO → 疾病 → 基因表型引擎；第 5.2 節為在 NGS-UI 後處理自建的 ACMG 變異規則引擎，取代以 GeneBe／pipeline ACMG 作為排序點數來源），尚未實作。計畫 §10.7 指出現行 `phenotype_to_genes.txt` 含 `NOT`（明確排除）與非表型 term，會被計入 pheno_score，可獨立先修。
- **STR / ROH 卡片**（目前還是「（無資料）」placeholder；STRchive / ROH summary 還沒接）。
- mito **haplogroup**（Haplogrep2 sidecar）—— 沒做。
- PharmCAT / PGx 卡片（payload 裡 `pharmcat: {}` 是空的）。
