# 歡迎使用 NGS 分析平台

成大醫院基因醫學部 NGS 分析平台是院內 NGS 變異判讀工具，用來整合二級 pipeline 產出的 SNV/Indel、CNV/SV 與 Mitochondria 變異資料，進行三級分析 pipeline，並協助醫師及生物資訊工程師載入個案、檢視變異卡片、標記變異、撰寫判讀意見，並輸出診斷報告。

請先在上方搜尋既有個案，或使用「載入新個案」登錄 pipeline 已輸出的樣本，或使用右上方「三級分析」開始新的樣本分析。

## 版本紀錄

### v4.9 — 2026-06-12

- GeneBe 第二意見改用 lazy SQLite cache：上傳新的 `genebe_hg38.tsv.gz` 後，下一次三級分析會自動偵測來源變更並重建 `genebe_hg38.sqlite`，之後查詢不再每個 sample 重掃整顆 GeneBe TSV。
- 三級分析 VCF 搜尋選到候選後，輸入框只保留 sample name；run、大小與日期仍在下拉候選與提示中顯示，方便連續修改下一個 sample ID 搜尋。
- Disease-associated / reportable gene list、VEP alias map 與 cohort dead-zone tables 更新至 2026-06-12 版本；固定 WES-I / WES-II / WGS panel 內的舊 gene symbol 與誤植同步修正。
- 固定 WES-I / WES-II / WGS panel data 與 custom panel data 改為 repo 內版本化資料，server 更新程式碼時會一併更新 panel。
- HPO/panel gene set、既有 `pheno_score.tsv` 與 SNV/CNV/SV/Mito 變異 gene 統一以 HGNC canonical symbol 做 `pheno_score` / `in_panel` 比對，降低 VEP 舊名或 panel alias 漏算。
- Custom panel 匯入與儲存時會先套用安全的 HGNC-current alias 轉換；alias 來源改用 HGNC 官方 complete set、withdrawn table 與人工確認表合併產生，無法唯一確認的項目會保留原字串供人工確認。
- 輸入臨床表徵工具新增 HPO / panel gene-list drawer，可查看、篩選與複製各 HPO term 或 panel 的 canonical gene list；另新增 gene lookup，可反查某 gene 出現在哪些 HPO terms / panels，單一 HPO 查看也加速避免等待整份 phenotype cache 預熱。
- 主畫面 Dead zone 卡片新增臨床門檻下的 CDS dead percentage，並依比例由高到低排序；≥70% 深紅、50-70% 紅、30-50% 橘，其餘黑色。
- IGV alignment track height 調整為 SNV/Indel 與 Mito 300，CNV/SV 50。
- 三級分析 log 將尾段統一顯示為 post-processing；Nextflow process 完成時間改附在原 stdout 行尾並以分鐘顯示，進度百分比保持單調遞增。

### v4.8 — 2026-06-11

- 三級分析進度條改用 Nextflow process 權重，`queued` 事件只記錄 log、不推進百分比；前處理與 prepare 類快速步驟不再讓進度過早衝到 80%。
- Nextflow 結束點調整為約 82%，保留較多空間給 copy、GeneBe/Extra VEP/AnnotSV、review TSV 與 gene index 預建。

### v4.7 — 2026-06-10

- SNV/Indel 新增 GIAB stratification 標籤：變異若落在 homopolymer、tandem repeat、segmental duplication、low mappability、GC 極端或其他困難區域，會在標籤列提示，提醒該位點 short-read 判讀較不可靠。

### v4.6 — 2026-06-08

- 三級分析改接 pipeline v3.4，支援 Pharmacogenomics 輸出。
- 新增預設勾選的 PGx checkbox，可在需要時關閉 PGx。

### v4.5 — 2026-06-07

- 三級分析改接 pipeline v3.2，支援 Mitochondria 及 CNV/SV Nextflow 輸出。
- 優化首頁載入及樣本載入速度。

### v4.4 — 2026-06-05

- 三級分析改接 v3.1 sample sheet pipeline。
- 三級分析支援批次送出同一來源類型的多個 sample，產生多列 sample sheet。
- 主畫面 SNV review TSV 預載層改為 AF < 0.01 且限 CDS ±50 bp 候選區域，ClinVar P/LP 與已標記變異仍會保留。

### v4.3 — 2026-06-04

- 改用 HGNC gene symbol 顯示，降低 VEP 舊基因名稱造成的判讀落差。
- 新增 `Disease-associated` 標籤，預設只顯示 disease-associated gene list 內的基因。
- 新增 Dead zone 提醒，依目前 HPO / panel 基因列出 coverage 低於判讀門檻的 exon；WES 使用 20X，WGS 使用 15X。
- 診斷報告之基因清單限制在 disease-associated gene list 中；WGS 報告會在有 dead-zone 的基因後以括號標註，WES 報告僅在主畫面顯示 Dead zone。

### v4.2 — 2026-05-30

- 新增 IGV 視窗，可從分析平台直接查看 BAM 檔中的 SNV/Indel 及 CNV/SV 。
- 新增 SRY 確認流程，支援性別確認。
- 整合相鄰 CNV/SV 片段，改用整合後的 variant；相鄰 gap 門檻調整為 250 kb。
- CNV/SV 分析區、報告區與 DOCX 匯出改依 phenotype 與 AnnotSV score 綜合排序。
- 固定 panel 改為 WES-I / WES-II / WGS / Other panel 分頁。
- 調整匯出報告細節，改善 CNV/SV 疾病描述、基因清單與版面呈現。

### v4.1 — 2026-05-28

- 加入匯出診斷報告功能，診斷報告依院內報告模板撰寫，包含 SNV/Mito/CNV/SV 描述。
- 加入低 read depth 警示。

### v4.0 — 2026-05-16

- 新增「三級分析」功能，可由平台啟動 DRAGEN / in-house VCF 三級分析流程。

### v3.0 — 2026-05-12

- 新增 Mitochondria 變異分析。
- 個案載入改為 staged loading：先載入 SNV/Indel，再背景載入 CNV/SV 與 Mito，以減少載入時等待時間。

### v2.1 — 2026-05-11

- 新增上傳個案清單功能，可匯入院內未完成報告清單，供平台整合個案基本資料、基因資料與臨床資料。
- 新增 SNV/Indel 與 CNV/SV 搜尋基因功能，並在各點位新增「搜尋同基因」功能。
- 新增 Candidate 報告區，供生物資訊工程師或醫師挑選合適之 Candidate variant。
- 整合 HPO/panel 輸入工具，支援自訂 panel。

### v2.0 — 2026-05-09

- 新增 CNV/SV 區塊，支援 AnnotSV 註解與輸出。
- CNV/SV 依 phenotype gene relevance 與 ACMG class 分區顯示。
- CNV/SV 加入涵蓋基因、疾病相關基因數、已知致病/良性區域重疊與 Disease 欄位。

### v1.0 — 2026-05-04

- 分析平台完成架設，支援 SNV/Indel 分析。
- 建立 FastAPI 後端、原生 JS 前端、SNV TSV adapter 與 variant card。
- 支援臨床資料以 HPO / gene panel 輸入，結合 Pheno score、Exomiser/LIRICAL 分析。
- 支援標記、comment、自動儲存、登入驗證。
- SNV/Indel tier 設定為 `1A / 1B / 1C / 2 / 3`，並依 variant score + phenotype score 排序。
