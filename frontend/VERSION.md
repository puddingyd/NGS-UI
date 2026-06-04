# 歡迎使用 NGS 分析平台

成大醫院基因醫學部 NGS 分析平台是院內三級分析判讀工具，用來整合次級 pipeline 產出的 SNV/Indel、CNV/SV 與 Mitochondria 變異資料，協助 reviewer 載入個案、檢視變異卡片、標記 causative / candidate / other、撰寫判讀意見，並輸出診斷報告。

請先在上方搜尋既有個案，或使用「載入新個案」登錄 pipeline 已輸出的樣本。平台功能若有影響判讀流程、報告輸出或資料載入方式的更新，後續都應評估是否寫入下方版本紀錄，讓使用者能追蹤重要變更。

## 版本紀錄

### v5.1 — 2026-06-02

- 整合相鄰 CNV/SV 片段，報告區、DOCX 與個案清單皆改用合併後的 parent variant。
- 調整 DOCX 報告細節，改善 CNV/SV 疾病描述、基因清單與版面呈現。
- 固定 panel 匯入邏輯更新，避免將疾病名稱或資料來源列誤當成基因。

### v5.0 — 2026-05-30

- 新增 IGV.js modal，可從 SNV/Indel、CNV/SV 卡片直接查看 BAM coverage。
- 新增 SRY 確認流程，支援性別確認與 sibling track 比對。
- 固定 panel 改為 WES-I / WES-II / WGS / Other panel 分頁，並支援匯入院內 panel Excel。

### v4.1 — 2026-05-29

- SNV/Indel adapter 改寫以支援新版 65 欄 pipeline TSV。
- GeneBe 結果改為 ACMG second opinion，並調整 tier / geno score 顯示邏輯。
- 加入舊格式提示與低 read depth 警示規則。

### v4.0 — 2026-05-16

- 新增「三級分析」功能，可由平台啟動 DRAGEN / in-house VCF 轉三級分析流程。
- 新增三級分析 modal、VCF typeahead、topbar progress 與 job 狀態追蹤。
- 串接 stop-gap annotation chain，包含 GeneBe、extra VEP、AnnotSV 與 review TSV 建立。

### v3.1 — 2026-05-28

- DOCX 診斷報告依院內模板重寫，包含字型、表格寬度、SNV/Mito/CNV/SV 描述格式。
- 報告區支援以正確卡片呈現 Mito 與 CNV/SV variant。
- Mito 加入 ClinVar runtime annotation 與手動 ACMG 欄位。

### v3.0 — 2026-05-12

- 新增 Mitochondria 變異分析卡片。
- 支援從 mitochondrial VCF 產生 `mito.annotated.tsv`。
- Mito 只顯示 `FILTER=PASS` 且具 MITOMAP 疾病關聯/致病性的位點。
- 個案載入改為 staged loading：先載入 SNV/Indel，再背景載入 CNV/SV 與 Mito。

### v2.1 — 2026-05-11

- 新增上傳個案清單功能，可匯入院內未完成報告清單 xlsx。
- 新增 SNV/Indel 與 CNV/SV gene search。
- 新增 Candidate 報告區，Secondary findings 改為折疊群組。
- 將 HPO/panel 輸入工具整合到 `/phenotype/`，支援自訂 panel。

### v2.0 — 2026-05-09

- 新增 CNV/SV 變異卡片，支援 AnnotSV 輸出。
- CNV/SV 依 phenotype gene relevance 與 ACMG class 分區顯示。
- 卡片加入涵蓋基因、疾病相關基因數、已知致病/良性區域重疊與 reviewer Disease 欄位。
- 自動寫入 `pheno_score.tsv`，讓 CNV/SV 與 SNV 都能使用 phenotype/panel 分數。

### v1.0 — 2026-05-04

- 第一版正式可用的 SNV/Indel 判讀平台。
- 建立 FastAPI 後端、原生 JS 前端、SNV TSV adapter 與 variant card。
- 支援 HPO / gene panel 編輯、pheno score、Exomiser/LIRICAL 重跑。
- 支援 reviewer 標記、comment、自動儲存、登入驗證與 DOCX 診斷報告匯出。
- SNV/Indel tier 改為 `1A / 1B / 1C / 2 / 3`，並依 genotype + phenotype score 排序。

### v0.3 — 2026-05-07

- 建立 analysis version 架構，支援多版本分析、重跑與版本切換。
- 「載入新個案」改為從 pipeline 已輸出的 sample directory 選取。
- 支援自動讀取 phenotype.txt、帶入 MRN、登錄後開始分析與自動儲存。

### v0.2 — 2026-05-05

- SNV/Indel 顯示加入 AF、in-panel、VAF 等篩選邏輯。
- SNV tier 改為水平 tab bar，改善大量 variant 的瀏覽方式。
- 將病人資料與 runtime data 移出 repo，改由 `NGS_UI_HOME` 管理。

### v0.1 — 2026-05-04

- 建立三級輸出格式規劃與第一批原型資料。
- 建立舊 webdata JSON 轉換成 `snv_indel.annotated.tsv` 的工具。
- 匯入 HPO、phenotype-to-gene 與 gene panel 參考資料。
