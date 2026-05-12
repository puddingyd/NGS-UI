# CLAUDE.md — 給 Claude（與接手者）的權威參考

成大醫院基因醫學部 **NGS 分析平台**（repo `puddingyd/NGS-UI`）。這份文件記錄整個系統的流程、結構、慣例與踩雷點，讓新的 session / auto-compact 後也能完全接續工作。**有任何架構變動，順手更新這份。**

---

## 0. 一句話總覽

醫院內部的 NGS 三級分析判讀工具：FastAPI 後端 + 原生 JS 前端（無 build step），部署在內網 `192.168.84.91:8765`。次級 pipeline（Nextflow，在另一台 compute cluster）產出 per-sample 的 TSV 丟進 `tertiary_output/{LIS_ID}/`；reviewer 在這個 UI 裡載入個案、看 SNV/Indel + CNV/SV + Mitochondria 的變異卡片、標記 causative/other/candidate、寫 comment、匯出診斷報告 (docx)。另有一個獨立的「輸入臨床表徵 (HPO/panel)」工具掛在 `/phenotype/`，和「上傳個案清單」功能建立 LIS_ID↔MRN↔姓名 對應。

---

## 1. Git 工作流（重要）

- 開發在這個 sandbox；推到分支 **`claude/plan-ngs-ui-RQW8J`**。
- **推 `main` 會被 proxy 擋 HTTP 403** —— 一律 `git push origin claude/plan-ngs-ui-RQW8J`。
- dev 機（`n102968@server`，repo 在 `~/NGS_UI/NGS-UI`，remote 設成 SSH `git@github.com:puddingyd/NGS-UI.git`）從那個分支 `git pull`；如需進 `main` 由使用者在 GitHub 開 PR 合併。
- push 失敗（網路）retry 最多 4 次 exponential backoff（2/4/8/16s）。不要建 PR 除非使用者明說。

---

## 2. 目錄佈局 / env var

`config.py` 裡所有路徑都從 `NGS_UI_HOME` 推導（env `NGS_UI_HOME`，否則 `REPO_ROOT.parent`（若 `../NGS-UI` 存在），再否則 `REPO_ROOT`）。每個都有自己的 env override。在 sandbox 裡 `NGS_UI_HOME` = repo root（沒有 parent NGS-UI），在 dev 機 `NGS_UI_HOME` = `~/NGS_UI`。

| 內容 | 預設路徑 | env override |
|---|---|---|
| 程式碼（這個 git checkout）= `REPO_ROOT` | `NGS_UI_HOME/NGS-UI/` | — |
| 每樣本 TSV + sidecar | `NGS_UI_HOME/tertiary_output/{LIS_ID}/` | `TERTIARY_OUTPUT_ROOT` |
| 樣本清單快取（每次 list 都被重寫） | `NGS_UI_HOME/tertiary_output/_index.json` | `NGS_UI_INDEX_PATH` |
| 伺服器執行狀態（`users.db`、`jobs/`） | `NGS_UI_HOME/data/` | `NGS_UI_DATA_ROOT` |
| 病患 phenotype.txt | `NGS_UI_HOME/patient_phenotype/` | `NGS_UI_PHENOTYPE_DIR` |
| 上傳的個案清單 xlsx + `roster.json` | `NGS_UI_HOME/patient_list/` | `NGS_UI_PATIENT_LIST_DIR` |
| HPO/panel 參考資料（`hp.obo`、`phenotype_to_genes.txt`、`gene_panels/`） | `NGS_UI_HOME/phenotype_data/` | `NGS_UI_PHENO_DATA_DIR` |
| gene panels（含使用者自訂的） | `NGS_UI_HOME/phenotype_data/gene_panels/` | `NGS_UI_GENE_PANELS_DIR` |
| OMIM.xlsx | `NGS_UI_HOME/OMIM/OMIM.xlsx` | `NGS_UI_OMIM_XLSX` |
| Exomiser/LIRICAL CLI | `NGS_UI_HOME/biotools/` | `NGS_UI_BIOTOOLS_DIR` |
| VCF（per-sample） | `NGS_UI_HOME/vcf/` | `NGS_UI_VCF_DIR` |
| Exomiser/LIRICAL 輸入模板 | `REPO_ROOT/phenotype_reference/`（**在 repo 裡，不在 NGS_UI_HOME**） | — |
| 前端靜態檔 | `REPO_ROOT/frontend/` | `FRONTEND_DIR` |
| EMR client id（NCKU intranet） | — | `NGS_UI_EMR_CLIENT_ID`（空 = 整套 EMR 功能關閉） |
| Redis（job queue） | `redis://127.0.0.1:6379/0` | `REDIS_URL` |
| Java / Exomiser 路徑等 | 見 `config.py` | `EXOMISER_HOME`、`LIRICAL_HOME`、`JAVA_BIN`… |

`.gitignore`：`tertiary_output/`、`data/`、`patient_list/`、`phenotype_data/`、`_index.json`、`__pycache__/`、`*.pyc`、`.venv/`、`node_modules/`。所以 panel 檔、phenotype 參考資料、樣本資料、roster 都**不在 git 裡**，部署時要自己放到 `NGS_UI_HOME/` 底下。`OMIM.xlsx` 原本在 repo 根，已 `git rm` 移出（dev 機放在 `NGS_UI_HOME/OMIM/`）。

**每樣本目錄 `tertiary_output/{LIS_ID}/` 內容：**
```
sample_metadata.json          patient-level：基本資料 + reviewer 編輯狀態 + active_analysis 指標 + tags/comment/status/edits/panels/manual_variants/clinical_description/genetic_counseling
snv_indel.annotated.tsv       SNV/Indel（pipeline 丟；各 analysis version 共用同一份）
cnv.annotated.tsv             CNV（AnnotSV 輸出；pipeline 丟）
sv.annotated.tsv              SV（AnnotSV 輸出；pipeline 丟）
mito.annotated.tsv            粒線體（由 scripts/annotate_mito_vcf.sh 產生）
qc_summary.json roh_summary.json   （前端的 QC 警告卡 / ROH）
analyses/{ver}/
  analysis.json               hpo + selected_panels + note（version-level）
  pheno_score.tsv             gene → 0-100 分數（write_version 的 side effect）
  exomiser_results.tsv lirical_results.tsv  （rerun worker 寫；可能不存在）
  analysis_files/             Exomiser/LIRICAL 的 run 目錄
  {LIS_ID}_{MRN}_phenotype.txt   audit copy（register 時若有 phenotype）
```
未登錄的樣本 = `tertiary_output/{X}/` 有 `snv_indel.annotated.tsv` 但沒有 `sample_metadata.json`；「載入新個案」modal 列這些。

---

## 3. 整體資料流

```
次級 pipeline (Nextflow, 別台 cluster)
  → tertiary_output/{LIS_ID}/{snv_indel,cnv,sv}.annotated.tsv
  → (mito: 在 dev 機跑 scripts/annotate_mito_vcf.sh --in <mito.vcf.gz> --outdir tertiary_output/{LIS_ID})
  → (HPO/panel: 用 /phenotype/ 工具產生 patient_phenotype/{...}_phenotype.txt)
  → (個案清單: 上傳 xlsx → patient_list/roster.json)

UI 流程：
  載入新個案 (modal) → POST /api/samples → 寫 sample_metadata.json + analyses/default/analysis.json
                                          + (write_version side-effect 寫 pheno_score.tsv) + 重寫 SNV TSV 的 IN_PANEL 欄
  選樣本 (combobox) → GET /api/samples/{id}  (核心 payload，aux_pending=true)
                    → 背景 GET /api/samples/{id}/cnv-sv 和 /mito
                    → GET /api/samples/{id}/report  (reviewer 編輯狀態)
  reviewer 在卡片標 1/2/C/0/X、寫 comment、勾 disease/gene checkbox、改 ACMG
  自動儲存 (1.5s debounce) → PUT /api/samples/{id}/report
  ▶ 開始分析 → POST /api/samples/{id}/phenotype (存 HPO/panels + 算 pheno_score)
            → POST /api/samples/{id}/jobs/exomiser_lirical (enqueue rerun worker)
  匯出診斷報告 → GET /api/samples/{id}/report.docx
```

---

## 4. Adapter / tier 結構

每種 variant type 一個 adapter（`backend/app/adapters/`），各回傳 `(variants_dict, categories_dict)`。`sample_loader.load_sample()` 把它們全部塞進一個 payload —— id namespace 不衝突（prefix 不同）。

| Type | adapter | id 格式 | tiers（payload key） | tier 規則 | 排序 |
|---|---|---|---|---|---|
| **SNV/Indel** | `snv_tsv.py` | `chr{N}-{pos}-{ref}-{alt}` | `1A 1B 1C 2 3`（互斥） | `classify_tier`：1A = ClinVar P/LP ≥1★；1B = ClinVar P/LP（任何）或 LOFTEE HC；1C = ACMG points ≥4；2 = ACMG points 1-3；3 = 其餘。`_normalize_acmg_class` 把 `VUS`→`Uncertain significance` 等 | 各 tier 內 `total_score`（= geno_score + pheno_score）desc，tie-break by id |
| **CNV** | `annotsv_tsv.py`（`source="cnv"`） | AnnotSV_ID | `CNV-1A`(Clinical) `CNV-1B`(Pathogenic)（**獨立分區，可重複**） | 1A = SV 涵蓋的任一 gene 在 pheno set（score>0）；1B = AnnotSV `ACMG_class` ∈ {4,5} | 1A：trigger gene 的 max pheno_score desc → ranking_score desc → id；1B：ranking_score desc → id。基因表 `genes` 切到前 10（後端，`genes_overflow` = in-panel 溢出、`genes_compact` = 其餘只含 `{gene,omim_id,in_panel}`） |
| **SV** | `annotsv_tsv.py`（`source="sv"`） | AnnotSV_ID | `SV-2A` `SV-2B` | 同 CNV，用 sv 檔 | 同 CNV |
| **Mito** | `mito_tsv.py` | `chrM-{pos}-{ref}-{alt}` | `MITO-1`(Pathogenic) `MITO-2`(Disease-associated)（互斥） | **先 `FILTER=PASS` 過濾**（GATK Mutect2-mito best-practices；非 PASS 全是 artifact/雜訊；TLOD 已內含在 FILTER 判定裡，不另設門檻）；**再只留 disease-relevant**（pathogenic 或有 `MITOMAP_DISEASE`）—— 沒 MITOMAP 紀錄的 polymorphism/haplogroup 變異完全不列。1 = MITOMAP status `Cfrm`/`Confirmed`/`[P]`/`[LP]` 或 MitoTIP (likely-)pathogenic；2 = 其餘有 disease | disease-relevance key：`(status_rank, mitotip_rank, in_panel(0/1), -refs, -heteroplasmy, pos)`（不再以 heteroplasmy 為主排序 —— heteroplasmy 高 ≠ 致病；只當 tie-breaker） |

**報告區**（`REPORT_SECTION_DEFS`，全在 `frontend/index.html` 的 `#report-sections` + Secondary-findings 折疊群組）：
- Causative（status `1`）、Other（`2`）、**Candidate（`C`）** —— 三段 default open，有 disease checkbox、可「＋ 新增 variant」。三段先按 `total_score` desc 排，再把**同基因的 cluster 在一起**（最高分基因的整組排最前；手動新增的無 gene_symbol 留原位）。
- ACMG SF / Proactive / Carrier / PharmCat —— 收在「Secondary findings」折疊群組（純文字標題 + 三角形鈕，無卡片框）。
- status 下拉值：**`—/1/2/C/0/X`**（C → Candidate 區）。SNV 卡片用 `statusOptions("candidate")`、CNV/SV/Mito 卡片 hardcode 同一組。

**「in_panel」概念**：`pheno_score.tsv` 裡 score>0 的基因 = 病人 HPO/panel 相關的基因（`phenotype_scorer.compute_pheno_match` 回 `{gene: matched_weight}` + total_weight；`compute_pheno_score` = 之後做 `100*matched/total` 正規化）。CNV/SV 的 Clinical tier、Mito 的排序 tie-breaker、CNV/SV 卡片基因表的 ⭐ 標記都用這個。CNV/SV 卡片「Pheno」欄顯示成 `matched/total`（乘 100 前的原始狀態）。`has_phenotype` = bool(hpo or panels)，前端用來在 CNV/SV Clinical 區空白時顯示「請先設定 phenotype」提示。

---

## 5. 前端版面 / 卡片（`frontend/index.html` + `app.js` + `style.css`）

- **Topbar**（深色 `#24292f`，z-index 110 蓋過 login modal）：左 hamburger（toggle `#sidebar`），中標題「成大醫院基因醫學部 NGS 分析平台」，右 `登入/登出`（同一顆鈕 toggle：`data-loggedIn` 切 handler）·`輸入臨床表徵 (HPO/panel)`（`<a href="/phenotype/" target="_blank">`，**未登入也顯示**）·`上傳個案清單`（xlsx upload）。`.btn` 用 `text-decoration:none` + inline-flex，所以 `<a class="btn">` 跟 `<button class="btn">` 一樣。
- **登入 modal**：「申請帳號：PYTHONPATH=backend python -m app create-user」那段提示字用 `.login-hint`/`.login-hint-code`（白底白字，反白才看得到）。
- **三個分析卡片**（`#card-snv`、`#card-cnv-sv`、`#card-mito`）各有 tier-tab bar + tier panels。`renderTierTabBar`/`renderCnvSvTabBar`/`renderMitoTabBar`；tab-click dispatch 統一處理三組（用 `data-tier` 判斷）。tier-panel 顏色：SNV 紅/黃系、CNV 藍 `#bfdbfe`、SV 紫 `#ddd6fe`、Mito teal `#99f6e4`。卡片要包在 `.block-body` 裡才有 inset 效果（`.tier-panel > .block-body { padding-top: 8px }`）。Mito/STR/ROH 卡片：STR 跟 ROH 還是「（無資料）」placeholder。
- **Staged loading**：`GET /samples/{id}` 只回核心（meta + reports + SNV + analyses + has_phenotype，`aux_pending: true`，CNV/SV/Mito 是空 dict）。前端 `loadSample` render 完後背景 `GET /samples/{id}/cnv-sv` 和 `/mito`，回來後 merge 進 `state.data` + re-render 那張卡。`state._auxLoadToken` 防 race（切樣本後晚到的回應丟掉）。等待時 CNV/SV、Mito panel 顯示「載入中…」、tab count「…」。
- **Gene 搜尋**：SNV/Indel 與 CNV/SV 標題列右側各有一個 gene 搜尋框，SNV 那邊還多 LIRICAL / Exomiser 兩顆按鈕 → 跳 `#gene-search-modal`（`max-width:1100px`，重用 `renderVariantCard`/`renderCnvSvCard`，所以卡片完整可互動）。
- **共用 class**：`.variant-head`（`#index` span + `.status-select` + ...）、`.btn-copy`（`COPY_ICON_SVG`）、`.ext-links`（Varsome 那種按鈕樣式）、`.cnv-sv-detail-box`（灰底兩行 + 折疊區）、`.cnv-sv-reasoning`、`.cnv-sv-comment-text`、`.acmg-class` + `.sig-p/.sig-lp/.sig-vus/.sig-lb/.sig-b` 五級色。**`.modal-card input[type=text]` 那條 catch-all（width:100%）被 `.variant-card`/`.cnv-sv-card` 裡的 input 排除**（不然 gene-search modal 裡的 ACMG points 那 3em 格會被撐爆）。
- **CNV/SV 卡片**（`renderCnvSvCard`）：每張左側無色條、inset 在 tier-panel 上；header = `#index` + status select（1/2/C/0/X）+ `CNV/SV` tag + SV-type pill（DEL紅/DUP藍/INV橙/INS紫/TRA灰）+ 座標 + 複製鈕 + `{chromN}{cytoband}`（如 `12p11.21-q24.33`）+ ext-links（UCSC/DECIPHER/dbVar/GeneCards）靠右；detail box 第一行 = ACMG 五級下拉（`.cnv-sv-acmg-select` + `sig-*` 色；reviewer override 存 `state.reports.edits[id].ACMG_class_sv`，跟 SNV ACMG 分開）+「涵蓋基因數: 1518（疾病相關：28）」+ 基因型 + Filter + Qual，第二行 = `AnnotSV 評分依據`（折疊）+ Score；基因表預設 10 列（首欄 checkbox → `state.reports.edits[id].report_genes`；Phenotype 與 Inheritance cell 點擊展開 `.gene-clip-cell`；overflow `<details>` lazy-render chips + in-panel 溢出用完整表格）；「已知致病區域重疊」「已知良性區域重疊」兩段（DEL→只 P_loss/B_loss、DUP→只 gain、其他→全部；內容 CSS `-webkit-line-clamp:2` + 「展開全部」鈕；無資料顯示提示而非整段消失）；最後 Comment textarea（→ `state.reports.edits[id].comment`）。
- **Mito 卡片**（`renderMitoCard`）：header = `#index` + status select + locus pill（protein紅/tRNA橙/rRNA藍/control灰）+ `m.HGVS` + 複製鈕 + gene + heteroplasmy 徽章（teal）+ gnomAD-MT/MITOMAP 連結；detail box 第一行 = `ref→alt` / 類型 / `Heteroplasmy (AD·DP)` / Filter（有 `ⓘ` tooltip：`MITO_FILTER_GLOSS` 解釋 Mutect2 旗標）；第二行 = Consequence（MITOMAP-only 的可讀標籤：missense / synonymous / stop_gained / non-coding (tRNA) / …）+ Protein change（MITOMAP `Amino Acid Change` 欄，如 `A52T`）+ TLOD（`ⓘ` tooltip）；折疊 MITOMAP 區（disease / status / plasmy reports / GenBank freq / MitoTIP / refs / allele）；Comment textarea。
- **OMIM disease list**（SNV 卡片，`renderDiseaseList`）：`<details>` 列出 `Disease1..5`（跳過 NA），summary 有報告勾選 checkbox（→ `state.reports.edits[id].report_diseases`），展開的黃底框（`.disease-detail`）底部有「▴ 收合」鈕。

---

## 6. pheno_score 自動寫入

`analyses_store.write_version()` 寫完 `analysis.json` 後 side-effect 算 `compute_pheno_score()` 並寫 `pheno_score.tsv`（HPO/panels 為空就刪掉舊檔）。所以 `register` 新個案、編輯 phenotype（`routers/phenotype.py`）、複製/重命名 version（`routers/analyses.py`）都會自動產生 pheno_score.tsv，不用等「▶ 開始分析」。`sample_loader` 還有 lazy backfill：載入時 pheno_score.tsv 缺失或比 analysis.json 舊就即時重算。`patient_store.register` 額外重寫 SNV TSV 的 `IN_PANEL` 欄。`write_pheno_table(sample_id, scores, target_dir=...)` 支援指定 version 目錄（不一定是 active 的）。

---

## 7. 轉換 / annotation scripts（`scripts/`）

| script | 用途 |
|---|---|
| `convert_anno_combined_to_tertiary_tsv.py` | 舊 R pipeline 的 `anno_combined.txt.gz` → `snv_indel.annotated.tsv`（去重 by 變異留最佳 transcript：MANE_SELECT > MANE_PLUS_CLINICAL > CANONICAL > any；缺的欄留空）。用法：`--in <file> --out tertiary_output/{LIS}/snv_indel.annotated.tsv` |
| `convert_old_json_to_tertiary_tsv.py` | 舊 webdata JSON → `snv_indel.annotated.tsv` |
| `annotate_mito_vcf.sh` + `parse_mito_vcf.py` | **MITOMAP-only（無 VEP）**。GATK Mutect2-mito VCF → `mito.annotated.tsv`：純 Python 讀 .vcf.gz、Python 端拆 multiallelic、HGVS_M 本地算（SNV `m.{pos}{ref}>{alt}`；indel 簡化 del/ins/dup）、gene/locus 用 rCRS 座標表（`_MT_GENES`/`_gene_at`，D-loop→`MT-CR`/control、OriL gap→`MT-OLR`、其他 gap→intergenic）、consequence + AA change 從 MITOMAP 的 `Amino Acid Change` 欄推、MITOMAP 只做精確 `(pos,ref,alt)` 比對（cc 用 `Nucleotide Change`、rna 用 `<ref><pos><alt>` 的 `Allele`；不做 POS-only fallback）、dedupe by `(pos,ref,alt)` 留 TLOD 最高。`MITOMAP_DIR` env（預設 `${REF_DIR:-/home/pipeline/reference/hg38}/tertiary/mitomap`）要有 `mitomap_mutations_coding_control.tsv`、`mitomap_mutations_rna.tsv`（**Latin-1 編碼**，loader 用 latin-1 讀）。輸出 `mito.annotated.tsv`（不帶 sample 前綴，直接放 `tertiary_output/{LIS}/`）。用法：`scripts/annotate_mito_vcf.sh --in <mito.vcf.gz> --sample {LIS} --outdir tertiary_output/{LIS}`。批次跑：`for S in ...; do IN=$(ls /home/datalake_Intermediate/.../*/$S/07_mitochondria/$S.mito.vcf.gz \| head -1); scripts/annotate_mito_vcf.sh --in "$IN" --sample "$S" --outdir ~/NGS_UI/tertiary_output/$S; done` |
| `migrate_to_versioned_layout.py` / `migrate_vcf_path.py` / `rewrite_vcf_paths.py` | 一次性的舊→新佈局遷移 |
| `probe_emr_api.py` | NCKU EMR API 診斷（urllib，dump 到 /tmp/emr_probe/；用內建 14 筆 MRN 清單，不讀 argv） |

`mito.annotated.tsv` 欄位（22）：`CHROM POS REF ALT HGVS_M GENE LOCUS_TYPE CONSEQUENCE AA_CHANGE HETEROPLASMY AD DEPTH FILTER TLOD MITOMAP_DISEASE MITOMAP_STATUS MITOMAP_PLASMY MITOMAP_GB_FREQ MITOMAP_GB_SEQS MITOMAP_REFS MITOTIP_SCORE MITOMAP_ALLELE`。

`mito.vcf.gz` 是 GATK Mutect2 `--mitochondria-mode` 輸出（FilterMutectCalls + 黑名單 mask；chrM = GRCh38 chrM = rCRS/NC_012920.1；`FORMAT/AF` = heteroplasmy fraction，`FORMAT/DP` = depth，`FORMAT/AD`）。

`snv_indel.annotated.tsv` 主要欄位：`CHROM POS REF ALT GENE TRANSCRIPT TRANSCRIPT_TYPE HGVS_C HGVS_P CONSEQUENCE MANE_ALL CALLERS ZYGOSITY GT_DV GT_HC AD VAF CLINVAR_SIG CLINVAR_STARS CLINVAR_DN CLINVAR_CONF GNOMAD_G_AF GNOMAD_G_EAS_AF GNOMAD_E_AF GNOMAD_E_EAS_AF TWB_AF PKNN_LLR REVEL BAYESDEL ALPHAMISSENSE METARNN ESM2_SCORE EVO2_SCORE SPLICEAI_MAX CADD_PHRED LOFTEE_HC LOFTEE_FILTER LOFTEE_FLAGS ACMG_EVIDENCE ACMG_POINTS ACMG_CLASS PHASE_GROUP PHASE_RESULT IN_ROH IN_PANEL IN_BLACKLIST OMIM_LINK GNOMAD_LINK CLINVAR_LINK REPORT_CLASS`。
CNV/SV 是 AnnotSV 標準輸出（128 欄；`Annotation_mode` full=一個 SV 一列、split=每 gene 一列；adapter `annotsv_tsv.py` 用 index-based 解析只取 ~30 欄、聚合 full+split）。

`tertiary_code/` 裡有次級 pipeline 的 Nextflow（`main_tertiary.nf`、`modules/`、`scripts/parse_vep_csq.py` 等）+ config（`dgm` profile = NCKU 正式環境 `/home/pipeline/reference/hg38`、`/home/pipeline/nextflow_containers`、`--bind /home`；`local` profile = `/scratch/pylin1991/...`、`/data/pylin1991/nf-containers`）—— 純參考用，這個 repo 不跑它。

---

## 8. 輸入臨床表徵工具（`/phenotype/`）

`frontend/phenotype/`（從舊的 GitHub-backed hpo-docs 改寫，砍掉 GitHub OAuth/terminal/run-analysis）。**不需登入**，由 NGS-UI 伺服器靜態服務在 `/phenotype/`（`main.py` 加 `GET /phenotype` → 307 redirect `/phenotype/`，註冊在 StaticFiles mount 之前）。功能：HPO term 搜尋（Fuse.js + 本地 `hpo_data.json` 3.5MB，**在 repo 裡** `frontend/phenotype/`）；Gene Panels 搜尋（打 `GET /api/phenotype-tool/panels`）；**Custom panel**（名稱 + 基因清單 textarea + weight；按「產生 phenotype.txt」時 POST `/api/phenotype-tool/custom-panel` 建檔到 `gene_panels/`、即時更新 `phenotype_scorer` 記憶體、名稱自動清理成 `[A-Za-z0-9_-]{1,64}`、衝突 409、基因**不大寫**（`C7orf50` 保留小寫）、case-sensitive 去重）。「產生 phenotype.txt」一鍵：建 custom panel → 組 TSV → POST `/api/phenotype-tool/save` 寫到 `patient_phenotype/`。MRN 或 LIS_ID 至少填一個；檔名：兩個都填 `{code}_{mrn}_phenotype.txt`、只 LIS_ID `{code}_phenotype.txt`、只 MRN `{mrn}_phenotype.txt`。「載入既有資料」用 `GET /api/phenotype-tool/load?code=&mrn=`。phenotype.txt 格式：`phenotype\thpo_name\tweight` 表頭 + `HP:xxxxxxx\t<name>\t<weight>` / `<panel_name>\t\t<weight>` 列（`phenotype_io.parse` 讀這個）。

`routers/phenotype_tool.py`：`GET /api/phenotype-tool/panels`、`POST /api/phenotype-tool/save`、`GET /api/phenotype-tool/load`、`POST /api/phenotype-tool/custom-panel` —— **全公開無 auth**（intranet 信任 + 嚴格驗證：token 限 `[A-Za-z0-9_-]{1,32}`、檔名從驗證過的 token 拼、內容 ≤64KB、panel 基因 ≤5000）。

---

## 9. 上傳個案清單（roster）

`patient_list_store.py`：上傳 NCKU「未完成報告清單」xlsx（`POST /api/patient_list`，**需登入**），原始檔存到 `patient_list/{ts}_{name}.xlsx`，merge 進 `patient_list/roster.json`（**additive**，不刪舊的）。xlsx 格式：找 col 0 == `檢體編號` 的標題列，砍 `8BB1` 前綴得 LIS_ID（`8BB126WE0092`→`26WE0092`），`檢驗名稱`→WES/WGS，by LIS_ID 去重，欄位 `檢體編號|病歷號|姓名|檢驗名稱|...|科別|...`。`sample_loader.list_unregistered()` 用 roster 自動填「載入新個案」modal 的 MRN/姓名/Test type（科別只當提示文字）；phenotype 檔查找順序：`{lis_id}_{roster_mrn}_phenotype.txt` → `{lis_id}_phenotype.txt` → `{lis_id}_*_phenotype.txt`（glob）→ `{roster_mrn}_phenotype.txt`。`GET /api/patient_list` = 看目前 roster（debug）。

> **`_index.json` 不要拿來放 roster** —— 它是「已登錄樣本清單快取」，`list_index()` 每次都重寫。roster 用獨立檔。

---

## 10. 其他

- **認證**：SQLite `data/users.db` + bcrypt（`users.py`）。建帳號：`PYTHONPATH=backend python -m app create-user [username]`（從 `backend/` 的 parent dir 跑；不用重啟服務）。`PYTHONPATH=backend python -m app list-users`。8h session cookie（SameSite=Lax）。沒有改密碼/刪帳號的指令。
- **OMIM annotation**：`omim_store.py` 啟動時讀 `OMIM.xlsx`（`_warm_caches` 預載；mtime 變了自動 reload；找不到檔就靜默 disable）。`sample_loader` 每個 SNV 變異 join `Disease1..5`/`OMIM_id`/`OMIM_disease`/`Inheritance`（OMIM_LINK 解析出的 OMIM_id 優先、gene_symbol fallback）。OMIM.xlsx 欄位：`OMIM_id | gene_symbol | OMIM_disease | Inheritance | Disease1..5 | Done`（17822 列，`OMIM_disease` 多行文字、每行 `<病名> (繼承碼)`）。`Disease1..5` 缺失時 `omim_store` 會從 `OMIM_disease` 的每一行合成。
- **自動儲存**：reviewer 編輯後 1.5s debounce → PUT `/samples/{id}/report`；三個位置的「儲存」按鈕（top/mid/bottom，class `.js-btn-save`/`.js-save-hint`）；存成功後 hint 顯示 `已儲存（HH:MM:SS）`；`beforeunload` 在 dirty/inflight 時警告；`_lastSavedAt` 切樣本時重設 null。
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
- CNV/SV `genes` 陣列在 adapter 裡切到前 10（不然跨染色體的大 DEL 會塞 1500+ gene record 到 payload）；前端 overflow chip lazy-render（`<details>` 的 toggle 事件不 bubble，用 capture-phase listener）。
- `tier-tab` 的 click dispatch 用 `data-tier` 值判斷 SNV / CNV-SV / Mito 哪一組（`CNV_SV_TIER_ORDER`、`MITO_TIER_ORDER` 的 includes 檢查）。
- ACMG 五級色 class `.sig-p/.sig-lp/.sig-vus/.sig-lb/.sig-b` 是全域的；CNV/SV 的 ACMG 下拉用它們但**不要**加 `.acmg-class` class（會觸發 SNV 的 change handler 寫錯 state）。
- Mito tier 1/2 是「先 PASS 過濾，再只留 disease-relevant」—— 沒 MITOMAP 紀錄的 polymorphism 完全不進 payload。
- 別 commit 沒被 `.gitignore` 排除的患者資料 / 大檔。

---

## 12. TODO / 還沒做

- CNV/SV / Mito 的 **docx 報告匯出**（目前 docx 只支援 SNV；CNV/SV 的 `report_genes`/`ACMG_class_sv`、Mito 的 reviewer 編輯都還沒接到 docx）。
- **STR / ROH 卡片**（目前還是「（無資料）」placeholder；STRchive / ROH summary 還沒接）。
- mito **haplogroup**（Haplogrep2 sidecar）—— 沒做。
- mito **gnomAD-mito / HelixMTdb 族群頻率**（目前只有 MITOMAP GenBank freq）。
- PharmCAT / PGx 卡片（payload 裡 `pharmcat: {}` 是空的）。
