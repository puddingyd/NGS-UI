# NGS 分析平台 — 實作筆記 / Implementation Notes

成大醫院基因醫學部 NGS 分析平台（`puddingyd/NGS-UI`，分支 `claude/plan-ngs-ui-RQW8J`）的關鍵設計決策。這份是「給未來的自己/接手者」的速查表，不是教學文件。

---

## 1. 目錄佈局 / env var

`config.py` 裡所有路徑都從 `NGS_UI_HOME` 推導（env `NGS_UI_HOME`，否則 `REPO_ROOT.parent`，再否則 `REPO_ROOT`）。每個都有自己的 env override：

| 內容 | 路徑（預設） | env |
|---|---|---|
| 程式碼（這個 git checkout） | `NGS_UI_HOME/NGS-UI/`（= `REPO_ROOT`） | — |
| 每樣本的 TSV + sidecar | `NGS_UI_HOME/tertiary_output/{LIS_ID}/` | `TERTIARY_OUTPUT_ROOT` |
| 樣本清單快取（每次 list 都重寫） | `NGS_UI_HOME/tertiary_output/_index.json` | `NGS_UI_INDEX_PATH` |
| 伺服器執行狀態（users.db、jobs/） | `NGS_UI_HOME/data/` | `NGS_UI_DATA_ROOT` |
| 病患 phenotype.txt | `NGS_UI_HOME/patient_phenotype/` | `NGS_UI_PHENOTYPE_DIR` |
| 上傳的個案清單 xlsx + roster.json | `NGS_UI_HOME/patient_list/` | `NGS_UI_PATIENT_LIST_DIR` |
| HPO/panel 參考資料（hp.obo、phenotype_to_genes.txt、gene_panels/） | `NGS_UI_HOME/phenotype_data/` | `NGS_UI_PHENO_DATA_DIR` |
| gene panels（含自訂） | `NGS_UI_HOME/phenotype_data/gene_panels/` | `NGS_UI_GENE_PANELS_DIR` |
| OMIM.xlsx | `NGS_UI_HOME/OMIM/OMIM.xlsx` | `NGS_UI_OMIM_XLSX` |
| Exomiser/LIRICAL CLI | `NGS_UI_HOME/biotools/` | `NGS_UI_BIOTOOLS_DIR` |
| Exomiser/LIRICAL 輸入模板 | `REPO_ROOT/phenotype_reference/`（**在 repo 裡，不在 NGS_UI_HOME**） | `FRONTEND_DIR`等另計 |
| EMR client id（NCKU intranet） | — | `NGS_UI_EMR_CLIENT_ID`（空 = 關閉 EMR 功能） |

`.gitignore`：`tertiary_output/`、`data/`、`patient_list/`、`phenotype_data/`、`_index.json`。所以 panel 檔、phenotype 資料、樣本資料都**不在 git 裡**，部署時要自己放到 `NGS_UI_HOME/` 底下。

每樣本目錄 `tertiary_output/{LIS_ID}/` 裡面：
```
sample_metadata.json          patient-level：基本資料 + reviewer 編輯 + active_analysis
snv_indel.annotated.tsv       SNV/Indel（各 analysis version 共用同一份 TSV）
cnv.annotated.tsv             CNV（AnnotSV 輸出；pipeline 丟）
sv.annotated.tsv              SV（AnnotSV 輸出）
mito.annotated.tsv            粒線體（由 scripts/annotate_mito_vcf.sh 產生）
qc_summary.json roh_summary.json
analyses/{ver}/
  analysis.json               hpo + selected_panels + note（version-level）
  pheno_score.tsv             gene → 0-100 分數（write_version 的 side effect）
  exomiser_results.tsv lirical_results.tsv  （rerun worker 寫）
```

---

## 2. Adapter / tier 結構

每種 variant type 一個 adapter（`backend/app/adapters/`），各自回傳 `(variants_dict, categories_dict)`，sample_loader 把它們全部塞進同一個 payload（id namespace 不衝突，因為 prefix 不同）。

| Type | adapter | id 格式 | tiers（payload key） | tier 規則 | 預設排序 |
|---|---|---|---|---|---|
| SNV/Indel | `snv_tsv.py` | `chr{N}-{pos}-{ref}-{alt}` | `1A 1B 1C 2 3`（互斥） | ClinVar / LOFTEE / ACMG points 分層 | 各 tier 內 total_score desc |
| CNV | `annotsv_tsv.py`（`source="cnv"`） | AnnotSV_ID | `CNV-1A`(Clinical) `CNV-1B`(Pathogenic)（**獨立分區，可重複**） | 1A = 任一 gene 在 pheno set；1B = AnnotSV ACMG_class ∈ {4,5} | 1A: max pheno_score desc；1B: ranking_score desc |
| SV | `annotsv_tsv.py`（`source="sv"`） | AnnotSV_ID | `SV-2A` `SV-2B`（同 CNV） | 同 CNV，用 sv 檔 | 同 CNV |
| Mito | `mito_tsv.py` | `chrM-{pos}-{ref}-{alt}` | `MITO-1`(Pathogenic) `MITO-2`(Disease-associated)（互斥） | **先過濾 FILTER=PASS**；1 = MITOMAP Cfrm/Confirmed/[P]/[LP] 或 MitoTIP (likely-)pathogenic；2 = 有 MITOMAP_DISEASE。**沒 MITOMAP 紀錄的 polymorphism/haplogroup 完全不列** | disease-relevance key：`(status_rank, mitotip_rank, in_panel, -refs, -heteroplasmy, pos)` |

報告區（`REPORT_SECTION_DEFS`）：Causative(1) / Other(2) / **Candidate(C)** / 然後 ACMG SF / Proactive / Carrier / PharmCat 收在「Secondary findings」折疊群組裡。status 下拉值：`—/1/2/C/0/X`（C → Candidate 區）。報告區的 1/2/C 三段會先按 total_score 排，再把**同基因的 cluster 在一起**（最高分基因的整組排最前）。

「in_panel」概念：`pheno_score.tsv` 裡 score>0 的基因 = 病人 HPO/panel 相關的基因。CNV/SV 的 Clinical tier、Mito 的排序 tie-breaker 都用這個。CNV/SV 卡片的「Pheno」欄顯示成 `matched/total`（用 `phenotype_scorer.compute_pheno_match`，乘 100 前的原始狀態）。

---

## 3. 前端版面 / 卡片

- Topbar（深色 `#24292f`）：左 hamburger（toggle sidebar），中標題，右 `登入/登出`（toggle 同一顆鈕）·`輸入臨床表徵 (HPO/panel)`（連到 `/phenotype/`，**不需登入也顯示**）·`上傳個案清單`（xlsx upload）。`.topbar` z-index 110（蓋過 login modal，所以登入前 phenotype 連結仍可點）。
- 三個分析卡片各有 tier-tab bar + tier panel（`renderTierTabBar`/`renderCnvSvTabBar`/`renderMitoTabBar`）。tab-click dispatch 統一處理三組（用 `data-tier` 判斷哪一組）。tier-panel 顏色：SNV 紅/黃系、CNV 藍 `#bfdbfe`、SV 紫 `#ddd6fe`、Mito teal `#99f6e4`。卡片包在 `.block-body` 裡才有 inset 效果。
- SNV/Indel 與 CNV/SV 標題列右側各有一個 gene 搜尋框 + SNV 那邊還有 LIRICAL / Exomiser 兩顆按鈕 → 跳 `#gene-search-modal`（重用 `renderVariantCard`/`renderCnvSvCard`，所以卡片完整可互動）。modal `max-width: 1100px`。
- 卡片共用 class：`.variant-head`（#index + status select + ...）、`.btn-copy`（COPY_ICON_SVG）、`.ext-links`（Varsome 那種按鈕樣式）、`.cnv-sv-detail-box`（灰底兩行 + 折疊區）、`.cnv-sv-comment-text`、`.acmg-class` + `.sig-p/.sig-lp/.sig-vus/.sig-lb/.sig-b` 五級色。
- CNV/SV 卡片：每張左側無色條，inset 在 tier-panel 色塊上；header 有 SV-type pill（DEL紅/DUP藍/INV橙/INS紫/TRA灰）、座標+複製鈕、`{chromN}{cytoband}`（如 `12p11.21-q24.33`）、ext links 靠右；detail box 第一行 ACMG 下拉 + 涵蓋基因數（`1518（疾病相關：28）`）+ 基因型 + Filter + Qual，第二行 AnnotSV 評分依據（折疊）+ Score；基因表預設 10 列（後端 `genes_overflow`/`genes_compact` 切，前端 lazy-render chip overflow），首欄 checkbox（→ `report_genes`），Phenotype 與 Inheritance cell 點擊展開（`.gene-clip-cell`）；「已知致病區域重疊」「已知良性區域重疊」兩段（DEL→只 P_loss/B_loss、DUP→只 gain、其他→全部；內容 CSS line-clamp 2 行 + 展開鈕；無資料顯示提示而非消失）。
- Mito 卡片：header 有 locus pill（protein紅/tRNA橙/rRNA藍/control灰）、`m.HGVS`+複製鈕、gene、heteroplasmy 徽章（teal）、gnomAD-MT/MITOMAP 連結；detail box 第一行 `ref→alt`/類型/`Heteroplasmy (AD·DP)`/Filter（有 `ⓘ` tooltip 解釋 Mutect2 旗標），第二行 Consequence/HGVSc/HGVSp/AA/TLOD（`ⓘ` tooltip）；折疊 MITOMAP 區（disease/status/plasmy/GenBank freq/MitoTIP/refs）；comment textarea。

---

## 4. pheno_score 自動寫入

`analyses_store.write_version()` 寫完 analysis.json 後會 side-effect 算 `compute_pheno_score()` 並寫 `pheno_score.tsv`（HPO/panels 為空就刪掉舊檔）。所以 register 新個案、編輯 phenotype、複製 version 都會自動產生 pheno_score.tsv，不用等「▶ 開始分析」。`sample_loader` 還有 lazy backfill：pheno_score.tsv 缺失或比 analysis.json 舊就即時重算。`patient_store.register` 額外重寫 SNV TSV 的 `IN_PANEL` 欄。`compute_pheno_score` = `compute_pheno_match`（回 `{gene: matched_weight}` + total_weight）後做 `100*matched/total` 正規化。

---

## 5. 轉換 / annotation script（scripts/）

| script | 用途 |
|---|---|
| `convert_anno_combined_to_tertiary_tsv.py` | 舊 R pipeline 的 `anno_combined.txt.gz` → `snv_indel.annotated.tsv`（去重 by 變異，留最佳 transcript；缺的欄留空） |
| `convert_old_json_to_tertiary_tsv.py` | 舊 webdata JSON → `snv_indel.annotated.tsv` |
| `annotate_mito_vcf.sh` | GATK Mutect2-mito VCF → `mito.annotated.tsv`：`bcftools norm` → VEP（chrM 自動用粒線體密碼子表，不掛 plugin）→ `parse_mito_vep.py`。輸出最終 TSV 不帶 sample 前綴（直接放 `tertiary_output/{LIS_ID}/mito.annotated.tsv`），中間 VEP VCF 帶前綴。路徑 env override，預設對齊 tertiary pipeline 的 `dgm` profile（`/home/pipeline/reference/hg38/...`、`/home/pipeline/nextflow_containers/vep_115.sif`、`bcftools_1.23.1.sif`、`--bind /home`） |
| `parse_mito_vep.py` | VEP VCF + 本地 MITOMAP 兩個 TSV（`mitomap_mutations_coding_control.tsv`、`mitomap_mutations_rna.tsv`，env `MITOMAP_DIR`）→ `mito.annotated.tsv`。CSQ 值會 URL-decode（`%3D`→`=`）；locus type 用 rCRS 基因座標表判（D-loop→`MT-CR`/control，gap→intergenic，OriL→`MT-OLR`）；MITOMAP 只做精確 `(pos,ref,alt)` 比對（不做 POS-only fallback，避免配錯 allele）；dedupe by `(pos,ref,alt)` 留 TLOD 最高的。Heteroplasmy = `FORMAT/AF`，depth = `FORMAT/DP` |
| `migrate_to_versioned_layout.py` / `migrate_vcf_path.py` / `rewrite_vcf_paths.py` | 一次性的舊→新佈局遷移 |
| `probe_emr_api.py` | NCKU EMR API 診斷（urllib，dump 到 /tmp/emr_probe/） |

`mito.annotated.tsv` 欄位：`CHROM POS REF ALT HGVS_M GENE LOCUS_TYPE CONSEQUENCE HGVS_C HGVS_P AA_CHANGE HETEROPLASMY AD DEPTH FILTER TLOD MITOMAP_DISEASE MITOMAP_STATUS MITOMAP_PLASMY MITOMAP_GB_FREQ MITOMAP_GB_SEQS MITOMAP_REFS MITOTIP_SCORE MITOMAP_ALLELE`。

> mito 判讀標準：**FILTER=PASS only**（GATK Mutect2-mito best-practices）。TLOD 已內含在 FILTER 判定裡，不另設門檻。

---

## 6. 輸入臨床表徵工具（/phenotype/）

`frontend/phenotype/`（從原本 GitHub-backed 的 hpo-docs 改寫，砍掉 GitHub/terminal/run-analysis）。**不需登入**，由 NGS-UI 伺服器靜態服務在 `/phenotype/`（`main.py` 加 `/phenotype` → `/phenotype/` redirect）。功能：HPO term 搜尋（Fuse.js + 本地 `hpo_data.json` 3.5MB，在 repo 裡）、Gene Panels 搜尋（打 `GET /api/phenotype-tool/panels`）、**Custom panel**（名稱 + 基因清單 textarea + weight；產生時 POST `/api/phenotype-tool/custom-panel` 建檔到 `gene_panels/`、即時更新 `phenotype_scorer` 記憶體、名稱自動清理成 `[A-Za-z0-9_-]`、衝突 409、基因不大寫不驗證、case-sensitive 去重）。「產生 phenotype.txt」一鍵：建 custom panel → 組 TSV → POST `/api/phenotype-tool/save` 寫到 `patient_phenotype/`。MRN 或 LIS_ID 至少填一個；檔名：兩個都填 `{code}_{mrn}_phenotype.txt`、只 LIS_ID `{code}_phenotype.txt`、只 MRN `{mrn}_phenotype.txt`。

`POST/GET /api/phenotype-tool/*` 在 `routers/phenotype_tool.py`，**公開無 auth**（intranet 信任 + 嚴格驗證）。

---

## 7. 上傳個案清單（roster）

`patient_list_store.py`：上傳 NCKU 的「未完成報告清單」xlsx（`POST /api/patient_list`，需登入），原始檔存到 `patient_list/{ts}_{name}.xlsx`，merge 進 `patient_list/roster.json`（additive，不刪舊的）。檔案格式：找到 col 0 == `檢體編號` 的標題列，砍 `8BB1` 前綴得 LIS_ID（`8BB126WE0092`→`26WE0092`），`檢驗名稱`→WES/WGS，by LIS_ID 去重。`list_unregistered()` 用 roster 自動填「載入新個案」modal 的 MRN/姓名/Test type（科別只當提示）；phenotype 檔查找順序：`{lis_id}_{roster_mrn}_phenotype.txt` → `{lis_id}_phenotype.txt` → `{lis_id}_*_phenotype.txt` → `{roster_mrn}_phenotype.txt`。

> **`_index.json` 不要拿來放 roster** — 它是「已登錄樣本清單快取」，`list_index()` 每次都重寫。roster 用獨立檔。

---

## 8. 其他

- 認證：SQLite `data/users.db` + bcrypt。建帳號：`PYTHONPATH=backend python -m app create-user [username]`（不用重啟）。`python -m app list-users`。
- OMIM annotation：`omim_store.py` 啟動時讀 OMIM.xlsx（`_warm_caches` 預載，mtime 變了自動 reload），sample_loader 每個 SNV 變異 join `Disease1..5`/`OMIM_id`/`OMIM_disease`/`Inheritance`（OMIM_LINK 的 id 優先、gene_symbol fallback）。前端 `renderDiseaseList` 渲染（含報告勾選 checkbox + 黃底框 + 「▴ 收合」鈕）。
- ACMG_CLASS 正規化：snv adapter 的 `_normalize_acmg_class`（`VUS`→`Uncertain significance` 等）。
- 自動儲存：reviewer 編輯後 1.5s debounce 自動 PUT `/samples/{id}/report`；三個位置的「儲存」按鈕（top/mid/bottom，class `.js-btn-save`/`.js-save-hint`）；存成功後 hint 顯示 `已儲存（HH:MM:SS）`。
- git 工作流：開發在 sandbox，推到 `claude/plan-ngs-ui-RQW8J` 分支（推 `main` 會被 proxy 擋 403）；dev 機從那個分支 pull。

---

## TODO / 還沒做

- CNV/SV / Mito 的 docx 報告匯出（目前 docx 只支援 SNV；CNV/SV 的 `report_genes`、`ACMG_class_sv` 跟 Mito 的 reviewer 編輯都還沒接到 docx）。
- STR / ROH 卡片（目前還是「（無資料）」placeholder）。
- mito 的 haplogroup（Haplogrep2 sidecar）— 目前沒做。
- gnomAD-mito / HelixMTdb 族群頻率（目前 mito 只有 MITOMAP GenBank freq）。
