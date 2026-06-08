# AGENTS.md

如果對話窗輸入 go，請做以下：

請只讀 **「待處理事項.md」**，先不要讀其他文件，先判斷哪些是必要讀的文件，然後再去讀和處理「待處理事項.md」中的事項。完成後請在 **「處理紀錄.md」** 摘要你做了什麼，並從「待處理事項.txt」 清掉已經完成的項目（若尚未完成，請留著那個項目，並加註尚未完成的部分）。

**有任何架構變動或功能增減，請順手更新這份 AGENTS.md、README.md 和 frontend/VERSION.md。**

---

## 0. 一句話總覽

成大醫院基因醫學部 **NGS 分析平台**（repo `puddingyd/NGS-UI`）。這份文件記錄整個系統的流程、結構、慣例與踩雷點，讓新的 session / auto-compact 後也能完全接續工作。

醫院內部的 NGS 三級分析判讀工具：FastAPI 後端 + 原生 JS 前端（無 build step），部署在內網 `192.168.84.91:8765`。次級 pipeline（Nextflow，在另一台 compute cluster）產出 per-sample 的 TSV 丟進 `tertiary_output/{LIS_ID}/`；reviewer 在這個 UI 裡載入個案、看 SNV/Indel + CNV/SV + Mitochondria 的變異卡片、標記 causative/other/candidate、寫 comment、匯出診斷報告 (docx)。另有一個獨立的「輸入臨床表徵 (HPO/panel)」工具掛在 `/phenotype/`，和「上傳個案清單」功能建立 LIS_ID↔MRN↔姓名 對應。

---

## 1. Git 工作流（重要）

- 開發在這個 sandbox；推到分支 **`Codex/plan-ngs-ui-RQW8J`**。
- **推 `main` 會被 proxy 擋 HTTP 403** —— 一律 `git push origin Codex/plan-ngs-ui-RQW8J`。
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
| panel_loose / HGNC alias / dead-zone tables | `REPO_ROOT/ngs_panel_deadzone/` | `NGS_UI_PANEL_DEADZONE_DIR` |
| OMIM.xlsx | `NGS_UI_HOME/OMIM/OMIM.xlsx` | `NGS_UI_OMIM_XLSX` |
| Exomiser/LIRICAL CLI | `NGS_UI_HOME/biotools/` | `NGS_UI_BIOTOOLS_DIR` |
| VCF（per-sample） | `NGS_UI_HOME/vcf/` | `NGS_UI_VCF_DIR` |
| IGV in-house BAM 搜尋根目錄 | `/home/datalake_Intermediate/pipeline/nextflow_output` | `NGS_UI_INHOUSE_BAM_ROOT`（可用 `:` 分隔多個 root；舊 `NGS_UI_BAM_ROOT` 仍作 fallback） |
| IGV DRAGEN BAM 搜尋根目錄 | `/home/datalake_Raw/Novaseq` | `NGS_UI_DRAGEN_BAM_ROOT`（可用 `:` 分隔多個 root；尋找 `<run>/bam/{sample}.bam`，排除 `{sample}.repeats.bam`） |
| IGV hg38 reference（FASTA + `.fai`） | `/home/pipeline/reference/hg38` | `NGS_UI_IGV_REF_DIR` |
| Exomiser/LIRICAL 輸入模板 | `REPO_ROOT/phenotype_reference/`（**在 repo 裡，不在 NGS_UI_HOME**） | — |
| 前端靜態檔 | `REPO_ROOT/frontend/` | `FRONTEND_DIR` |
| 首頁歡迎文字與版本紀錄 | `REPO_ROOT/frontend/VERSION.md`（前端啟動時讀取並顯示在未載入個案的首頁） | — |
| EMR client id（NCKU intranet） | — | `NGS_UI_EMR_CLIENT_ID`（空 = 整套 EMR 功能關閉） |
| Redis（job queue） | `redis://127.0.0.1:6379/0` | `REDIS_URL` |
| Java / Exomiser 路徑等 | 見 `config.py` | `EXOMISER_HOME`、`LIRICAL_HOME`、`JAVA_BIN`… |

`.gitignore`：`tertiary_output/`、`data/`、`patient_list/`、`phenotype_data/`、`_index.json`、`__pycache__/`、`*.pyc`、`.venv/`、`node_modules/`。所以 panel 檔、phenotype 參考資料、樣本資料、roster 都**不在 git 裡**，部署時要自己放到 `NGS_UI_HOME/` 底下。`OMIM.xlsx` 原本在 repo 根，已 `git rm` 移出（dev 機放在 `NGS_UI_HOME/OMIM/`）。

**每樣本目錄 `tertiary_output/{LIS_ID}/` 內容：**
```
sample_metadata.json          patient-level：基本資料 + reviewer 編輯狀態 + active_analysis 指標 + tags/comment/status/edits/panels/manual_variants/clinical_description/genetic_counseling
case_summary.json             個案清單摘要持久化快取：causative / disease / other variant / comment / sign_received_at，依 metadata、TSV/index、OMIM signature 自動失效
pipeline_source.json          三級分析來源 audit：source_path/source_sample_id/source_vcf_path/copied_at
snv_indel.annotated.tsv       SNV/Indel 完整來源（pipeline 丟；各 analysis version 共用）
snv_indel.review.tsv          主畫面用衍生檔（自動重建；candidate BED 內 AF<0.01 / AF缺值 / ClinVar P/LP；reviewer 已標記 SNV 另由 index 補入）
snv_indel.review.tsv.source.json  衍生檔 manifest（raw mtime/size + AF 門檻 + candidate BED signature；不含 reviewer keep_ids）
snv_gene_index.sqlite         完整 TSV gene search index（gene → raw TSV byte offset；tertiary job/stopgaps 預建）
cnv.annotated.tsv             CNV（AnnotSV 輸出；pipeline 丟）
sv.annotated.tsv              SV（AnnotSV 輸出；pipeline 丟）
mito.annotated.tsv            粒線體（三級 pipeline v3.2 的 `04_mito/{source_sample_id}.mito.tsv` 複製而來；舊樣本可能仍是 `scripts/annotate_mito_vcf.sh` 輸出格式）
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
  → (mito: v3.2 pipeline 輸出 04_mito/{SAMPLE_ID}.mito.tsv；NGS-UI worker 複製成 tertiary_output/{LIS_ID}/mito.annotated.tsv)
  → (HPO/panel: 用 /phenotype/ 工具產生 patient_phenotype/{...}_phenotype.txt)
  → (個案清單: 上傳 xlsx → patient_list/roster.json)

UI 流程：
  載入新個案 (modal) → POST /api/samples → 寫 sample_metadata.json + analyses/default/analysis.json
                                          + (write_version side-effect 寫 pheno_score.tsv；SNV in_panel 動態由 pheno_score.tsv 補值)
                                          + 記錄既有 vcf_from_tsv.vcf.gz；缺 VCF 時交給 Exomiser/LIRICAL worker 背景建立
  個案清單 (modal) → DELETE /api/samples/{id}?delete_pipeline_output={bool}
                  → 刪 NGS_UI_HOME/tertiary_output/{id}/
                  → 三選一確認：同時刪 /home/pipeline/tertiary_output/{id}/ + NGS-UI job logs（執行中拒絕刪除）、保留三級分析原始檔案、取消全部刪除
  選樣本 (combobox) → GET /api/samples/{id}  (核心 payload，aux_pending=true)
                    → 背景分別 GET /api/samples/{id}/cnv、/sv 和 /mito（CNV 可先顯示）
                    → GET /api/samples/{id}/report  (reviewer 編輯狀態)
  reviewer 在卡片標 1/2/C/0/X、寫 comment、勾 disease/gene checkbox、改 ACMG
  自動儲存 (1.5s debounce) → PUT /api/samples/{id}/report
  ▶ 開始分析 → POST /api/samples/{id}/phenotype (存 HPO/panels + 算 pheno_score)
            → POST /api/samples/{id}/jobs/exomiser_lirical (enqueue rerun worker)
  匯出診斷報告 → GET /api/samples/{id}/report.docx
  輸出 PDF → 瀏覽器列印視窗（報告區 causative / other / candidate 卡片摘要）
```

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

**「in_panel」概念**：`pheno_score.tsv` 裡 score>0 的基因 = 病人 HPO/panel 相關的基因（`phenotype_scorer.compute_pheno_match` 回 `{gene: matched_weight}` + total_weight；`compute_pheno_score` = 之後做 `100*matched/total` 正規化）。CNV/SV 的 Clinical tier、Mito 的排序 tie-breaker、CNV/SV 卡片基因表的 ⭐ 標記都用這個。CNV/SV 卡片「Pheno」欄顯示成 `matched/total`（乘 100 前的原始狀態）。`has_phenotype` = bool(hpo or panels)，前端用來在 CNV/SV Clinical 區空白時顯示「請先設定 phenotype」提示。`Disease-associated` 是另一層報告基因範圍，來自 `panel_loose.hgnc_canonical.txt`；SNV 主畫面預設用它過濾，DOCX §五.4 也只列這個範圍內的 HPO/panel genes。Mito disease 來自 `CLINVAR_DN`，前端用 `&` 拆成 checkbox，勾選者才進 DOCX 疾病句。`panel_deadzone.py` 會用 `HGNC_ID` 優先把 VEP symbol 轉成 HGNC-current 名稱，並提供 WES 20X / WGS DRAGEN 15X dead-zone exon 查詢。

---

## 5. 前端版面 / 卡片（`frontend/index.html` + `app.js` + `style.css`）

- **Topbar**（深色 `#24292f`，z-index 110 蓋過 login modal）：左 hamburger（toggle `#sidebar`），中標題「成大醫院基因醫學部 NGS 分析平台」，右 `登入/登出`（同一顆鈕 toggle：`data-loggedIn` 切 handler）·`輸入臨床表徵 (HPO/panel)`（`<a href="/phenotype/" target="_blank">`，**未登入也顯示**）·`上傳個案清單`（xlsx upload）。`.btn` 用 `text-decoration:none` + inline-flex，所以 `<a class="btn">` 跟 `<button class="btn">` 一樣。
- **首頁歡迎 / 版本紀錄**：未載入個案時，搜尋卡片下方顯示 `#welcome-card`，內容由 `frontend/VERSION.md` 讀取並用前端的簡易 Markdown renderer 呈現；載入任一 sample 後自動隱藏。之後若有影響判讀流程、報告輸出、資料載入或主要工具入口的更新，要評估是否同步更新 `frontend/VERSION.md`。
- **登入 modal**：「申請帳號：PYTHONPATH=backend python -m app create-user」那段提示字用 `.login-hint`/`.login-hint-code`（白底白字，反白才看得到）。
- **三個分析卡片**（`#card-snv`、`#card-cnv-sv`、`#card-mito`）各有 tier-tab bar + tier panels。`renderTierTabBar`/`renderCnvSvTabBar`/`renderMitoTabBar`；tab-click dispatch 統一處理三組（用 `data-tier` 判斷）。tier-panel 顏色：SNV 紅/黃系、CNV 藍 `#bfdbfe`、SV 紫 `#ddd6fe`、Mito teal `#99f6e4`。卡片要包在 `.block-body` 裡才有 inset 效果（`.tier-panel > .block-body { padding-top: 8px }`）。Mito/STR/ROH 卡片：STR 跟 ROH 還是「（無資料）」placeholder。
- **Staged loading**：`GET /samples/{id}` 只回核心（meta + reports + review TSV 的 SNV + analyses + has_phenotype，`aux_pending: true`，CNV/SV/Mito 是空 dict）。前端 `loadSample` render 完後背景分別打 `GET /samples/{id}/cnv`、`/sv` 和 `/mito`，CNV 載完會先顯示，SV 可繼續在背景載入；回來後 merge 進 `state.data` + re-render 那張卡。`state._auxLoadToken` 防 race（切樣本後晚到的回應丟掉）。等待時對應 CNV/SV、Mito panel 顯示「載入中…」、tab count「…」。SNV tier 只 render active tab，切 tab 時才建立該 tier 卡片 DOM。
- **SNV 顯示 filter**：分析區標題列依序是 `Disease-associated`（預設勾；限 `panel_loose`，只影響 SNV 主畫面，不限制 CNV/SV 或 gene search modal）、`In panel only`（預設勾）、`gnomAD_G_AF < 0.01`（預設勾；缺值保留）、`VAF ≥ 0.2`（預設勾；缺值隱藏）、`impact=MODIFIER`（預設不勾；勾後才顯示 MODIFIER）、OMIM 顯示 toggle。這些都是 UI-only，不刪原始 TSV；主畫面 review TSV 預載層只保留 `NGS_UI_CDS_CANDIDATE_BED` 內且 AF<0.01/AF缺值的點，ClinVar P/LP 與 reviewer 已標記點不受 BED/AF 限制，其他超過者仍可用 gene search 查回。`IMPACT=LOW` 一律保留。tier tab count 跟著 filter 重算。若 reviewer 標記不在 `panel_loose` 的 SNV，分析區會提示 DOCX §五.4 不會列該 gene，但不阻擋標記。
- **ESM1b 配色**：依 Bergquist et al. 2025 ClinGen SVI 校準區間壓成五色；ESM1b 與一般 predictor 方向相反，分數越低越偏致病（`≤-24` P、`≤-10.7` LP、`≤-6.2` VUS、`≤8.7` LB、其餘 B）。多 transcript TSV 用最低分作 worst case。
- **個案搜尋 filter**：主畫面搜尋框上方是可複選、可取消的 `WES` / `WGS` 圓形 chips；`matchSamples()` 與 Enter 精準查找都只在勾選的 test type 內搜尋。提示文字放在 input placeholder。
- **個案清單管理**：modal 內有獨立的 `WES` / `WGS` 圓形 chips 與全文搜尋框，會搜尋 ID、姓名、MRN、variant、疾病、comment 與日期等所有表格內容。首頁 `/api/samples` 只回輕量搜尋索引；打開個案清單才呼叫 `/api/samples/case-summary` 載入摘要。表格摘要 status `1` 的 causative SNV/CNV/SV、實際勾選的 SNV OMIM diseases 與 CNV/SV reviewer `Disease`、status `2` 的 other SNV/CNV/SV、主畫面 comment、簽收時間與載入時間；多個 variant / disease 逐行顯示，長 HGVS 可在任意字元換行。摘要依 metadata / TSV/index / OMIM signature 使用 process LRU + per-sample `case_summary.json` 持久化快取；SNV 摘要優先用 `snv_gene_index.sqlite` 依 variant id seek，CNV/SV 摘要只讀目標 id，避免 DRAGEN/WGS 大型 TSV 在首頁或重啟後被完整重掃。
- **Gene 搜尋**：SNV/Indel 與 CNV/SV 標題列右側各有一個 gene 搜尋框，支援以 `,`、`，`、`、` 或空白分隔多個基因；SNV 那邊還多 LIRICAL / Exomiser 兩顆按鈕 → 跳 `#gene-search-modal`（`max-width:1100px`，重用 `renderVariantCard`/`renderCnvSvCard`，所以卡片完整可互動）。SNV 會打 `GET /samples/{id}/snv-search?genes=...` 查完整 `snv_indel.annotated.tsv`，不受 review TSV 影響；搜尋 modal 預設勾選 `gnomAD_G_AF < 0.01`，取消後才顯示全部搜尋結果；搜尋 modal 內的 SNV 卡片也會用 `MANE_ALL` 把 ENST 顯示轉成 MANE SELECT `NM_`。純 CNV/SV 搜尋不顯示該 toggle。WGS raw TSV 不再於載入個案後自動背景預熱；三級分析結尾預建 `snv_gene_index.sqlite`（gene → raw TSV byte offset），搜尋時只從 raw TSV seek 指定 gene rows，舊樣本缺 index 才 fallback raw parse。
- **IGV.js modal**：SNV/Indel、Mitochondria 與 CNV/SV 卡片的 ext-links 最左側有 `IGV`；基本資料性別下方另有「☑ 已確認」核取方塊與「確認 SRY」，以 hg19/hg38 對應 gene region 開啟相同 modal。SRY 模式初始只帶當前 sample，仍可手動加入 sibling，alignment coverage 固定 `autoscale:false,max:100`；勾選狀態存 `reports.sry_confirmed`。一般 modal 點擊後標題顯示 sample、padded locus、variant HGVS/類型與原始 `[build]` 座標；alignment track 預設 `displayMode: "SQUISHED"`、`visibilityWindow: 5000000`。SNV/Indel 與 Mito 顯示前後 100bp，CNV/SV 原則上顯示原區間前後各 20% flanking area；若事件本身 ≤5 Mb，padding 會縮小以保持初始視窗 ≤5 Mb，直接載入 BAM coverage 而不先顯示 zoom-in 提示。browser-level `roi` 仍標出未加 padding 的實際 CNV/SV 區間；CNV/SV modal 內所有 BAM alignment coverage tracks 使用同一個 `autoscaleGroup`，讓 y-axis data range 一致以便比較 deletion/duplication；染色體先經 `_normalizeChrom()` 轉成 `chrN`。modal 先列 primary BAM + 同 batch 最多 2 個 sibling，可增刪 track，確認後才按「載入 IGV」。三級分析 job state 與 `pipeline_source.json` 保存原始 `source_sample_id` 與 `pipeline_type`；`GET /api/igv/bams` 會依來源分流，sample ID `-dragen` suffix 強制 DRAGEN roots、`-inhouse` 強制 in-house roots，無 suffix hint 才看 sidecar。DRAGEN 找 `NGS_UI_DRAGEN_BAM_ROOT/<run>/bam/{sample}.bam` 並排除 `{sample}.repeats.bam`，in-house 找 `NGS_UI_INHOUSE_BAM_ROOT/<batch>/{sample}/02_alignment/`。若 reviewer 自訂輸出 sample ID，會用 sidecar 回查原始 BAM ID；舊資料則移除已知 suffix 後 fallback 到原始 sample ID，並可容忍 in-house alignment 目錄中唯一的 `*.aligned.sorted.bam` 檔名不完全等於 sample ID。IGV modal 另有「其他路徑」下拉，呼叫 `/api/igv/bam-folders` 列出 DRAGEN `<run>/bam/` 與 in-house `<batch>/<sample>/02_alignment/`，選資料夾後再用 `/api/igv/bam-folder?dir=...` 列出直下一層 `.bam`（排除 `.repeats.bam`）；選 primary 後，同 batch 加入清單改用同資料夾其他 BAM。`routers/igv.py` 的 `/bams`、`/batch-samples`、`/bam-folders`、`/bam-folder`、`/genome`、`/file` 全部需登入；`/file` 支援 HTTP range，且只允許 BAM roots 下的 BAM/BAI/CRAM/CRAI 與 `NGS_UI_IGV_REF_DIR` 下的 reference 檔。本機 hg38 FASTA + `.fai` 避免 igv.js 連 hospital intranet 擋住的 AWS S3。
- **載入效能**：`snv_review.ensure_review_tsv()` 自動維護 `snv_indel.review.tsv`，主畫面只預載 `NGS_UI_CDS_CANDIDATE_BED` 內且 `GNOMAD_G_AF < 0.01` 或 AF 缺值的點，ClinVar P/LP rescue 不受 BED/AF 限制；reviewer 已標記但不在 review TSV 的 SNV 由 `snv_gene_index.sqlite` 依 variant id 補入，不會讓 review TSV 因 status 變動而重建。`run_stopgaps.sh` 在三級分析結尾先執行 `build_snv_review_tsv.py` 與 `build_snv_gene_index.py`；載入時 lazy rebuild 僅作為舊樣本 fallback。`sample_loader` 的 process-local bounded LRU 依 TSV / analysis / pheno / Exomiser / LIRICAL / OMIM signature 自動失效；完整 raw SNV TSV gene search 與個案清單 SNV 摘要優先走 `snv_gene_index.sqlite`，不在載入個案後自動預熱。`/api/samples` 不再計算 case-list summary；summary 延後到 `/api/samples/case-summary` 並以 `case_summary.json` 跨重啟保留。刪除個案只 invalidate 該 sample 的 SNV 與個案摘要 cache，不同步重掃全部清單。大型 JSON response 用 gzip。OMIM enrichment 每批只檢查一次 xlsx mtime。FastAPI startup 用 daemon thread 背景預熱 HPO、phenotype gene map、OMIM 與 mito ClinVar，避免大型 `phenotype_to_genes.txt` 擋住 HTTP port。
- **三級分析 modal / job log**：in-house 與 DRAGEN VCF 各自使用單一 typeahead 輸入格；點入輸入格即展開全部 VCF，輸入 sample / run / path 後即時縮小候選清單。In-house 選到 VCF 後會依檔案大小預設 `seq_type`（<100 MB = WES，>=100 MB = WGS），並顯示 WES/WGS segmented control 讓 reviewer 手動覆寫；CNV/SV/Mito sibling 路徑只收在可展開的「來源檔案」偵錯區，正式流程優先使用 pipeline v3.2 自動輸出的 `04_mito/` 與 `06_cnv_sv/`。Sample ID 可修改作為 UI 輸出資料夾與檔名前綴；選好後可按「加入批次」，同一批只能包含同一種來源（全 in-house 或全 DRAGEN），按「開始分析」後以多列 v3.x sample sheet 一次送出。worker 另外保存原始 `source_sample_id`，v3.x 預設為每個 job 產生 `samplesheet.csv`（`sample_id,pipeline_type,input_dir,seq_type,hpo`），執行 `nextflow ... --samplesheet ... --pipeline_type nckuh|dragen --out_dir /home/pipeline/tertiary_output -resume`；若 `NGS_UI_TERTIARY_ENV_SCRIPT`（預設 `/home/pipeline/pipeline_code/NGS2ndAnalysis_env.sh`）存在會先 source，不存在則沿用目前環境直接執行 `nextflow`。pipeline 以 source sample ID 輸出 `/home/pipeline/tertiary_output/{source_sample_id}/03_acmg/{source_sample_id}.snv_indel.acmg.tsv`、v3.2 `04_mito/{source_sample_id}.mito.tsv`，以及 `06_cnv_sv/{source_sample_id}.{cnv,sv}.annotated.tsv`；NGS-UI 逐筆複製成 `NGS_UI_HOME/tertiary_output/{Sample ID}/snv_indel.annotated.tsv`、`mito.annotated.tsv`、`cnv.annotated.tsv`、`sv.annotated.tsv` 並寫 `pipeline_source.json`（含 source VCF、source sample ID 與 `pipeline_type`，供載入新個案自動判斷 WGS）。若 pipeline CNV/SV 兩檔都複製成功，worker 會在 stop-gaps 加 `--skip-cnv`，避免再次跑本機 AnnotSV fallback；缺任一檔則保留 fallback。copy 前會先建立 `tertiary_output/{Sample ID}/`，避免 pipeline 已成功但 UI 目標資料夾不存在時失敗。只有設定 `NGS_UI_TERTIARY_LEGACY_STAGING=1` 時才回到舊 `stage_dragen_for_tertiary.sh` 單樣本 staging 流程。VCF 建立時間來自檔案 `mtime` Unix timestamp，前端固定以 `Asia/Taipei`（UTC+8）顯示。工具列內的 Extra VEP 預設不勾選，與不換行的 `↻ 更新索引`、`三級分析清單` 固定同高；tooltip 只描述補 SpliceAI 註解，只有勾選時才執行 Extra VEP。Extra VEP 和 GeneBe 一樣先送 `GNOMAD_G_AF ≤ 0.01` 或 AF 缺值候選點，再限於 `NGS_UI_CDS_CANDIDATE_BED`（預設 `$HOME/NGS_UI/biotools/cds_combined.bed`，MANE Select CDS±10bp + chrM / RefSeq supplement）內執行，最後把 annotation merge 回完整 TSV。執行進度顯示細分進度條：queued 1%、detect/samplesheet 2%、Nextflow 3-87%（VEP done 約 75%）、copy 87%、可見 stop-gaps / sample 預建 87-99%。Nextflow stdout、stop-gaps marker 與 sample-step marker 會即時更新 state step；黑底 log 預設收合，按右側三角形展開。ClinVar fallback 預設不跑且不顯示；`filter_snv_tsv.py` 仍在 GeneBe 前背景清理 TSV，但不顯示 log section；skip 的 Extra VEP / AnnotSV 不顯示空 section；review TSV / gene index 是 sample 預建步驟，不再以 stopgaps 命名，log 只保留 start/done/elapsed，隱藏中間 perf 行。worker 自己寫的 log 行、Nextflow process start/done（`[nextflow-step]`）、stop-gaps 與 sample 預建子步驟都帶台北時區 ISO timestamp；Nextflow 每個 process 的 start/done/elapsed 另存 `state.json.nextflow_step_history`，整體子程序結束時也記錄 elapsed seconds；`state.json` 另存 `step_started_at` + `step_history`。modal 內的「三級分析清單」合併 `PIPELINE_OUT_ROOT` 與 NGS-UI job state，因此失敗且沒有 output 目錄的 sample 仍會顯示 log；刪除時會一併刪除 `/home/pipeline/tertiary_output/{sample}/`、`NGS_UI_HOME/tertiary_output/{sample}/`；batch job log 會保留給同批其他 sample 查詢。新 job 狀態寫到 `data/jobs/tertiary/`；舊版 `data/jobs/dragen/` 僅作 read-only fallback。`run_stopgaps.sh` 不再建立 `.raw` snapshot；GeneBe 暫存 VCF 帶 primary contig header，AnnotSV 成功執行時只寫摘要 log，失敗才附尾端診斷。
- **三級分析 output reuse 保護**：重用既有 `/home/pipeline/tertiary_output/{sample}/` output 時只會剝除純 UI alias `-dragen` / `-inhouse`；`-WES` / `-WGS` 是 test-type ID，不可剝除，避免 `VAL-24-WGS` 誤用 `VAL-24` 舊 WES TSV。
- **三級分析進度 / 取消**：進度條依 6-sample DRAGEN WGS batch log 重新配重；Nextflow 佔主要時間，VEP 完成約 75%，Nextflow 結束約 87%，尾段依 `stopgap_sample_index/stopgap_sample_count` 與可見子步驟（GeneBe、可選 Extra VEP、可選 AnnotSV、review TSV、gene index）推進到完成。首頁每 15 秒查詢 `/api/dragen/jobs` 恢復 active job，所以其他瀏覽器或電腦登入後也能看到進度；進度面板「終止」按鈕會 `POST /api/dragen/jobs/{job_id}/cancel`，後端向 worker process group 送 SIGTERM 並標記 `cancelled`。
- **DOCX CNV/SV 表格**：「變異位置」欄使用 buffered wrap，內容寬度比欄寬少一格，避免座標字串貼到後面的拷貝數欄。
- Mito TSV 複製完成後會呼叫 `mitomap_mito.annotate_mito_tsv()` 靜默回補 `MITOMAP_*` 欄位；這不是進度 step，不會在 job log 額外顯示。回補欄位只給 `mito_tsv.py` 做 tier 分類，前端卡片不呈現 MITOMAP。
- **共用 class**：`.variant-head`（`#index` span + `.status-select` + ...）、`.btn-copy`（`COPY_ICON_SVG`）、`.ext-links`（Varsome 那種按鈕樣式）、`.cnv-sv-detail-box`（灰底兩行 + 折疊區）、`.cnv-sv-reasoning`、`.cnv-sv-comment-text`、`.acmg-class` + `.sig-p/.sig-lp/.sig-vus/.sig-lb/.sig-b` 五級色。**`.modal-card input[type=text]` 那條 catch-all（width:100%）被 `.variant-card`/`.cnv-sv-card` 裡的 input 排除**（不然 gene-search modal 裡的 ACMG points 那 3em 格會被撐爆）。
- **CNV/SV 卡片**（`renderCnvSvCard`）：每張左側無色條、inset 在 tier-panel 上；header = `#index` + status select（1/2/C/0）+ `CNV/SV` tag + SV-type pill（DEL紅/DUP藍/INV橙/INS紫/TRA灰）+ 座標 + 複製鈕 + `{chromN}{cytoband}`（如 `12p11.21-q24.33`）+ ext-links（UCSC/DECIPHER/dbVar/GeneCards）靠右；detail box 第一行 = ACMG 五級下拉（`.cnv-sv-acmg-select` + `sig-*` 色；reviewer override 存 `state.reports.edits[id].ACMG_class_sv`，跟 SNV ACMG 分開）+「涵蓋基因數: 1518（疾病相關：28）」+ 基因型 + Filter + Qual，第二行 = `AnnotSV 評分依據`（折疊）+ Score；基因表預設只顯示 phenotype 相關（`in_panel=true`）基因（首欄 checkbox → `state.reports.edits[id].report_genes`；Phenotype 與 Inheritance cell 點擊展開 `.gene-clip-cell`），按小三角形才 lazy-render phenotype score 為 0 的其餘完整列 / compact chips；「已知致病區域重疊」「已知良性區域重疊」兩段（DEL→只 P_loss/B_loss、DUP→只 gain、其他→全部；內容 CSS `-webkit-line-clamp:2` + 「展開全部」鈕；無資料顯示提示而非整段消失）；最後 Disease textarea（→ `state.reports.edits[id].disease`，DOCX 與個案清單使用）與 Comment textarea（→ `state.reports.edits[id].comment`）。同來源、同染色體、同 DEL/DUP 且相鄰 gap ≤250 kb 的片段預設自動整合；copy number 差異不再阻擋視覺整合，原始 segments 仍可展開查看各自 CN。UI / DOCX / 個案清單改用 parent CNV/SV。parent 取合併片段中最高 combined score 的代表資料，分析區與報告區皆依 `max_pheno_score + scaled AnnotSV score` 排序；前端依 sample payload 快取整合 view，render tiers 與 cards 時不會反覆排序及重建 parent。
- **Mito 卡片**（`renderMitoCard`）：header = `#index` + status select + locus pill（protein紅/tRNA橙/rRNA藍/control灰）+ mitochondrial `m.HGVS`（v3.2 由 `POS/REF/ALT` 產生，不使用 transcript `HGVS_C` 的 `c.`）+ 複製鈕 + gene + heteroplasmy 徽章（teal）+ gnomAD-MT / ClinVar 連結；detail box 第一行 = ACMG 五級下拉（`ACMG_classification_mito` reviewer override，報告區同步重 render）+ `ref→alt` / 類型 / `Heteroplasmy (AD·DP)` / genotype，只有來源 TSV 真的有 `FILTER` 時才顯示 Filter；第二行 = Consequence + Impact + Biotype + Protein change（v3.2 `HGVS_P` 會 URL decode 並移除 protein transcript，只留報告用 `p.xxx`）+ ClinVar + gnomAD-mito AF + TLOD（缺值顯示 `—`）。MITOMAP 欄位只供後端 tier 分類，MITOMAP 區塊與 MITOMAP 連結不顯示；`CLINVAR_DN` 以 `&` 拆成 ClinVar disease checkbox，勾選者才進 DOCX。
- **OMIM disease list**（SNV 卡片，`renderDiseaseList`）：`<details>` 列出完整 `Disease1..5`（跳過 NA），summary 有報告勾選 checkbox（→ `state.reports.edits[id].report_diseases`），展開的黃底框（`.disease-detail`）底部有「▴ 收合」鈕。只有「個案清單」摘要會呼叫 `omim_store.compact_disease_label()`，把 curator 長篇說明截到第一個可辨識的遺傳模式括號（如 `(AD)`、`(AR)`）結束；variant 卡片、原始 `OMIM.xlsx` 與 `OMIM_disease` 都保留完整內容。

---

## 6. pheno_score 自動寫入

`analyses_store.write_version()` 寫完 `analysis.json` 後 side-effect 算 `compute_pheno_score()` 並寫 `pheno_score.tsv`（HPO/panels 為空就刪掉舊檔）。所以 `register` 新個案、編輯 phenotype（`routers/phenotype.py`）、複製/重命名 version（`routers/analyses.py`）都會自動產生 pheno_score.tsv，不用等「▶ 開始分析」。`sample_loader` 還有 lazy backfill：載入時 pheno_score.tsv 缺失或比 analysis.json 舊就即時重算。SNV 的 `in_panel` 不再寫回 raw TSV 的 `IN_PANEL` 欄；主畫面載入與 gene search 都用當前 active analysis 的 `pheno_score.tsv` 動態補 `in_panel` / `pheno_score` / `total_score`。`write_pheno_table(sample_id, scores, target_dir=...)` 支援指定 version 目錄（不一定是 active 的）。

---

## 7. 轉換 / annotation scripts（`scripts/`）

| script | 用途 |
|---|---|
| `convert_anno_combined_to_tertiary_tsv.py` | 舊 R pipeline 的 `anno_combined.txt.gz` → `snv_indel.annotated.tsv`（去重 by 變異留最佳 transcript：MANE_SELECT > MANE_PLUS_CLINICAL > CANONICAL > any；缺的欄留空）。用法：`--in <file> --out tertiary_output/{LIS}/snv_indel.annotated.tsv` |
| `convert_old_json_to_tertiary_tsv.py` | 舊 webdata JSON → `snv_indel.annotated.tsv` |
| `filter_snv_tsv.py` | GUI TSV 的結構性清理：只移除 `REF/ALT=*` 與非 primary contig；AF、VAF、IMPACT 留在 TSV，由 UI 顯示 filter 控制。v3.1 正式 pipeline 自行處理 NCKUH/DRAGEN 前處理；舊 DRAGEN staging 只在 `NGS_UI_TERTIARY_LEGACY_STAGING=1` 時啟用 |
| `build_snv_review_tsv.py` | 從完整 `snv_indel.annotated.tsv` 預建主畫面用 `snv_indel.review.tsv`；`run_stopgaps.sh` 結尾自動執行 |
| `build_snv_gene_index.py` | 從完整 `snv_indel.annotated.tsv` 預建 `snv_gene_index.sqlite`，只保存 canonical gene、variant id 與 raw TSV byte offset/length，供 SNV gene search 快速 seek 指定 gene rows；`run_stopgaps.sh` 結尾自動執行 |
| `run_annotsv_cnv_sv.sh` | NGS-UI CNV/SV AnnotSV 統一入口：DRAGEN 用 hard-filtered VCF 尋找 sibling `{sample}.cnv.vcf.gz` / `{sample}.sv.vcf.gz`，in-house 用明確的 gCNV + Delly VCF，輸出 `cnv.annotated.tsv` / `sv.annotated.tsv`；目前用法摘要在 `docs/annotsv_current_usage.md` |
| `compare_genebe_spliceai_coverage.py` | 診斷 GeneBe 本機資料庫的 SpliceAI 完整性：從 `genebe_hg38.tsv.gz` 抽樣 variant，使用目前 extra-VEP 的 VEP SpliceAI plugin 路徑重算 `SPLICEAI_MAX`，輸出 coverage summary、mismatch 明細與 run metadata，供評估是否能用本機 GeneBe DB 取代 extra-VEP / GeneBe API；正式估計建議 `--max-sites 100000 --sample-mode chrom-balanced`。 |
| `import_fixed_panels.py` | `fixed_panel/WES-I.xlsx`、`WES-II.xlsx` 與 `other_panel/` → `phenotype_data/fixed_panels/index.json` + `gene_panels/*.txt`。WES Excel 只讀 `gene panel list` 標記列起始的基因區塊，不可把疾病名或資料來源列當成基因；顯示名稱只移除科別/級別前綴，需保留 `Non-syndromic` 這類 panel 名內的連字號。持久化 key 為相容既有 analysis metadata 仍沿用舊規則，和顯示名稱分開。更新來源 Excel 後要重新執行匯入。 |
| `annotate_dragen_mito_vcf.sh` | DRAGEN hard-filtered VCF 內含 chrM calls 時的 legacy mito wrapper。直接把 `{sample}.hard-filtered.vcf.gz` 交給 `parse_mito_vcf.py`，只取 `chrM/MT/chrMT` rows，`FORMAT/AF` 當 heteroplasmy、`FORMAT/SQ` 在缺 `INFO/TLOD` 時填入 TLOD 欄，輸出 `<outdir>/mito.annotated.tsv`。用法：`scripts/annotate_dragen_mito_vcf.sh --in /path/{sample}.hard-filtered.vcf.gz --sample {sample} --outdir tertiary_output/{sample}_legacy_mito`。 |
| `annotate_mito_vcf.sh` + `parse_mito_vcf.py` | Legacy/手動補跑用的 **MITOMAP-only（無 VEP）** mito 轉檔。正式三級分析目前優先使用 pipeline v3.2 `04_mito/{sample}.mito.tsv`，由 worker 複製成 `mito.annotated.tsv`。舊 script 可把 GATK Mutect2-mito VCF → `mito.annotated.tsv`：純 Python 讀 .vcf.gz、Python 端拆 multiallelic、HGVS_M 本地算（SNV `m.{pos}{ref}>{alt}`；indel 簡化 del/ins/dup）、gene/locus 用 rCRS 座標表（`_MT_GENES`/`_gene_at`，D-loop→`MT-CR`/control、OriL gap→`MT-OLR`、其他 gap→intergenic）、consequence + AA change 從 MITOMAP 的 `Amino Acid Change` 欄推、MITOMAP 只做精確 `(pos,ref,alt)` 比對（cc 用 `Nucleotide Change`、rna 用 `<ref><pos><alt>` 的 `Allele`；不做 POS-only fallback）、dedupe by `(pos,ref,alt)` 留 TLOD 最高。`MITOMAP_DIR` env（預設 `${REF_DIR:-/home/pipeline/reference/hg38}/tertiary/mitomap`）要有 `mitomap_mutations_coding_control.tsv`、`mitomap_mutations_rna.tsv`（**Latin-1 編碼**，loader 用 latin-1 讀）。輸出 `mito.annotated.tsv`（不帶 sample 前綴，直接放 `tertiary_output/{LIS}/`）。用法：`scripts/annotate_mito_vcf.sh --in <mito.vcf.gz> --sample {LIS} --outdir tertiary_output/{LIS}`。 |
| `migrate_to_versioned_layout.py` / `migrate_vcf_path.py` / `rewrite_vcf_paths.py` | 一次性的舊→新佈局遷移 |
| `probe_emr_api.py` | NCKU EMR API 診斷（urllib，dump 到 /tmp/emr_probe/；用內建 14 筆 MRN 清單，不讀 argv） |

`mito.annotated.tsv` adapter 同時支援舊欄位（22）：`CHROM POS REF ALT HGVS_M GENE LOCUS_TYPE CONSEQUENCE AA_CHANGE HETEROPLASMY AD DEPTH FILTER TLOD MITOMAP_DISEASE MITOMAP_STATUS MITOMAP_PLASMY MITOMAP_GB_FREQ MITOMAP_GB_SEQS MITOMAP_REFS MITOTIP_SCORE MITOMAP_ALLELE`，以及 v3.2 pipeline mito 欄位：`CHROM POS REF ALT GENE HGVS_C HGVS_P CONSEQUENCE IMPACT BIOTYPE GENOTYPE DP AF_SAMPLE GNOMAD_MITO_AF CLINVAR_SIG CLINVAR_DN GNOMAD_MITO_AF_HOM GNOMAD_MITO_AF_HET GNOMAD_MITO_AN CLINVAR_VARIATION_ID OMIM_IDS PIPELINE`。`AF_SAMPLE→heteroplasmy`、`DP→depth`、`HGVS_P→aa_change`（清成 `p.xxx`），`HGVS_C` 是 transcript `c.`，不拿來當卡片標題；adapter 會用 `POS/REF/ALT` 產生 `m.{pos}{ref}>{alt}`。若 TSV 含 `FILTER`，adapter 會照舊先保留 PASS；若沒有 `FILTER`，不做後端 FILTER 篩選，UI 也不顯示 Filter 欄。v3.2 TSV 複製後可由 `mitomap_mito.annotate_mito_tsv()` 精確 `(POS,REF,ALT)` 回補 MITOMAP 欄位；UI/報告卡片不呈現 MITOMAP，但 tier 1 會使用 MITOMAP confirmed/pathogenic 狀態。

`mito.vcf.gz` 是 GATK Mutect2 `--mitochondria-mode` 輸出（FilterMutectCalls + 黑名單 mask；chrM = GRCh38 chrM = rCRS/NC_012920.1；`FORMAT/AF` = heteroplasmy fraction，`FORMAT/DP` = depth，`FORMAT/AD`）。

`snv_indel.annotated.tsv` 主要來自 v3.1 `03_acmg/{SAMPLE_ID}.snv_indel.acmg.tsv`（65 欄），含 `CHROM POS REF ALT RS_ID GENE TRANSCRIPT TRANSCRIPT_TYPE HGVS_C HGVS_P CONSEQUENCE IMPACT EXON INTRON MANE_ALL CALLERS DP_DV AD_DV VAF_DV DP_HC AD_HC ZYGOSITY GT_DV GT_HC GNOMAD_G_AF GNOMAD_G_EAS_AF GNOMAD_E_AF GNOMAD_E_EAS_AF GNOMAD_E_AF_DBNSFP GNOMAD_E_EAS_AF_DBNSFP TG_EAS_AF CLINVAR_SIG CLINVAR_STARS CLINVAR_DN CLINVAR_SIGCONF CLINVAR_VARIATION_ID OMIM_IDS LOFTEE LOFTEE_FILTER LOFTEE_FLAGS LOFTOOL BAYESDEL_NOAF BAYESDEL_NOAF_PRED ALPHAMISSENSE ALPHAMISSENSE_PRED ESM1B ESM1B_PRED VARITY_R SIFT SIFT_PRED DANN PHACTBOOST PHYLOP100 GERP PKNN_LLR PKNN_EVIDENCE PANGOLIN_SCORE PANGOLIN_DETAIL DOMAINS SWISSPROT HGNC_ID ACMG_CRITERIA ACMG_SCORE ACMG_CLASS ACMG_NOTES`；NGS-UI stop-gaps 可能再補 `GENEBE_*`、`METARNN`、`SPLICEAI_MAX`、`IN_PANEL` 等 UI 欄位。
CNV/SV 是 AnnotSV 標準輸出（128 欄；`Annotation_mode` full=一個 SV 一列、split=每 gene 一列；adapter `annotsv_tsv.py` 用 index-based 解析只取 ~30 欄、聚合 full+split）。

`tertiary_code/` 裡有次級 pipeline 的 Nextflow（`main_tertiary.nf`、`modules/`、`scripts/parse_vep_csq.py` 等）+ config（`dgm` profile = NCKU 正式環境 `/home/pipeline/reference/hg38`、`/home/pipeline/nextflow_containers`、`--bind /home`；`local` profile = `/scratch/pylin1991/...`、`/data/pylin1991/nf-containers`）—— 純參考用。UI worker 使用正式環境 `/home/pipeline/tertiary_code/{main_tertiary.nf,scripts}`，Nextflow config 預設改用 `/home/n102968/NGS_UI/nextflow_tertiary_no_pgx.config`（可用 `NGS_UI_TERTIARY_CONFIG` 覆寫）以跳過 PGX；repo 內的 workflow/scripts 版本可能落後，不可整包覆寫正式 `scripts_dir`。

---

## 8. 輸入臨床表徵工具（`/phenotype/`）

`frontend/phenotype/`（從舊的 GitHub-backed hpo-docs 改寫，砍掉 GitHub OAuth/terminal/run-analysis）。**不需登入**，由 NGS-UI 伺服器靜態服務在 `/phenotype/`（`main.py` 加 `GET /phenotype` → 307 redirect `/phenotype/`，註冊在 StaticFiles mount 之前）。功能：HPO term 搜尋（Fuse.js + 本地 `hpo_data.json` 3.5MB，**在 repo 裡** `frontend/phenotype/`）；Gene Panels 搜尋（打 `GET /api/phenotype-tool/panels`）；**Custom panel**（名稱 + 基因清單 textarea + weight；按「產生 phenotype.txt」時 POST `/api/phenotype-tool/custom-panel` 建檔到 `gene_panels/`、即時更新 `phenotype_scorer` 記憶體、名稱自動清理成 `[A-Za-z0-9_-]{1,64}`、衝突 409、基因**不大寫**（`C7orf50` 保留小寫）、case-sensitive 去重）。「產生 phenotype.txt」一鍵：建 custom panel → 組 TSV → POST `/api/phenotype-tool/save` 寫到 `patient_phenotype/`。MRN 或 LIS_ID 至少填一個；檔名：兩個都填 `{code}_{mrn}_phenotype.txt`、只 LIS_ID `{code}_phenotype.txt`、只 MRN `{mrn}_phenotype.txt`。「載入既有資料」用 `GET /api/phenotype-tool/load?code=&mrn=`。phenotype.txt 格式：`phenotype\thpo_name\tweight` 表頭 + `HP:xxxxxxx\t<name>\t<weight>` / `<panel_name>\t\t<weight>` 列（`phenotype_io.parse` 讀這個）。

主畫面的 Patient phenotype card 與「載入新個案」modal 都有 `WES-I / WES-II / WGS / Other panel` tabs，預設展開 `Other panel`；固定 panel chip 由 `GET /api/phenotype-tool/fixed-panels` 載入，勾選後以 weight=1 寫入各自的 panels working copy，Other panel typeahead 則排除固定 panel key，避免同一個 panel 出現兩個入口。Patient phenotype、新個案 modal 與 `/phenotype/` 工具的 HPO / panel 下拉都支援鍵盤上下鍵與 Enter 選取；dropdown 開啟時 Enter 只會選項目，不會送出載入個案或開始分析。Patient phenotype 和 Comment 中間的 Dead zone card 依目前 HPO + panel gene set 顯示 cohort dead exons；WES 用 20X，WGS（in-house/DRAGEN 都一樣）用 DRAGEN 15X；主畫面預設只顯示前 10 列，可用小三角形展開全部。平台啟動讀取索引、載入個案核心資料或登錄新個案期間會顯示不可用背景點擊或 ESC 誤關閉的「資料載入中」遮罩；Exomiser/LIRICAL job 完成後的 sample refresh 使用 silent reload，不再用全頁遮罩阻擋判讀。

`routers/phenotype_tool.py`：`GET /api/phenotype-tool/panels`、`POST /api/phenotype-tool/save`、`GET /api/phenotype-tool/load`、`POST /api/phenotype-tool/custom-panel` —— **全公開無 auth**（intranet 信任 + 嚴格驗證：token 限 `[A-Za-z0-9_-]{1,32}`、檔名從驗證過的 token 拼、內容 ≤64KB、panel 基因 ≤5000）。

---

## 9. 上傳個案清單（roster）

`patient_list_store.py`：上傳 NCKU「未完成報告清單」xlsx（`POST /api/patient_list`，**需登入**），原始檔存到 `patient_list/{ts}_{name}.xlsx`，merge 進 `patient_list/roster.json`（**additive**，不刪舊的）。xlsx 格式：找 col 0 == `檢體編號` 的標題列，砍 `8BB1` 前綴得 LIS_ID（`8BB126WE0092`→`26WE0092`），`檢驗名稱`→WES/WGS，by LIS_ID 去重，欄位 `檢體編號|病歷號|姓名|檢驗名稱|...|科別|...`。`sample_loader.list_unregistered()` 用 roster 自動填「載入新個案」modal 的 MRN/姓名/Test type（科別只當提示文字）；若 `pipeline_source.json` 顯示 DRAGEN 來源，或 in-house source VCF 大於 100 MB，前端會把 Test type 預設為 WGS 並覆蓋 roster 預設。phenotype 檔查找順序：`{lis_id}_{roster_mrn}_phenotype.txt` → `{lis_id}_phenotype.txt` → `{lis_id}_*_phenotype.txt`（glob）→ `{roster_mrn}_phenotype.txt`。`GET /api/patient_list` = 看目前 roster（debug）。

> **`_index.json` 不要拿來放 roster** —— 它是「已登錄樣本清單快取」，`list_index()` 每次都重寫。roster 用獨立檔。

---

## 10. 其他

- **認證**：SQLite `data/users.db` + bcrypt（`users.py`）。建帳號：`PYTHONPATH=backend python -m app create-user [username]`（從 `backend/` 的 parent dir 跑；不用重啟服務）。`PYTHONPATH=backend python -m app list-users`。8h session cookie（SameSite=Lax）。沒有改密碼/刪帳號的指令。
- **OMIM annotation**：`omim_store.py` 啟動時讀 `OMIM.xlsx`（`_warm_caches` 預載；mtime 變了自動 reload；找不到檔就靜默 disable）。`sample_loader` 每個 SNV 變異 join `Disease1..5`/`OMIM_id`/`OMIM_disease`/`Inheritance`（OMIM_LINK 解析出的 OMIM_id 優先、gene_symbol fallback）。OMIM.xlsx 欄位：`OMIM_id | gene_symbol | OMIM_disease | Inheritance | Disease1..5 | Done`（17822 列，`OMIM_disease` 多行文字、每行 `<病名> (繼承碼)`）。`Disease1..5` 缺失時 `omim_store` 會從 `OMIM_disease` 的每一行合成。
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
- CNV/SV `genes` 陣列在 adapter 裡切到前 10（不然跨染色體的大 DEL 會塞 1500+ gene record 到 payload）；前端 overflow chip lazy-render（`<details>` 的 toggle 事件不 bubble，用 capture-phase listener）。
- `tier-tab` 的 click dispatch 用 `data-tier` 值判斷 SNV / CNV-SV / Mito 哪一組（`CNV_SV_TIER_ORDER`、`MITO_TIER_ORDER` 的 includes 檢查）。
- ACMG 五級色 class `.sig-p/.sig-lp/.sig-vus/.sig-lb/.sig-b` 是全域的；CNV/SV 的 ACMG 下拉用它們但**不要**加 `.acmg-class` class（會觸發 SNV 的 change handler 寫錯 state）。
- Mito tier 1/2/3 是「有 FILTER 則先 PASS 過濾，無 FILTER 則全收」；MITOMAP 只回補欄位與輔助 tier 分類，不做 disease-relevant 預篩，也不顯示在卡片。其他 PASS mtDNA variants 會進 `MITO-3`。
- 別 commit 沒被 `.gitignore` 排除的患者資料 / 大檔。

---

## 12. TODO / 還沒做

- **STR / ROH 卡片**（目前還是「（無資料）」placeholder；STRchive / ROH summary 還沒接）。
- mito **haplogroup**（Haplogrep2 sidecar）—— 沒做。
- PharmCAT / PGx 卡片（payload 裡 `pharmcat: {}` 是空的）。
