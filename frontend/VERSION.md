# 歡迎使用 NGS 分析平台

成大醫院基因醫學部 NGS 分析平台是院內 NGS 變異判讀工具，用來整合二級 pipeline 產出的 SNV/Indel、CNV/SV 與 Mitochondria 變異資料，進行三級分析 pipeline，並協助醫師及生物資訊工程師載入個案、檢視變異卡片、標記變異、撰寫判讀意見，並輸出診斷報告。

請先在上方搜尋既有個案，或使用「載入新個案」登錄 pipeline 已輸出的樣本，或使用右上方「三級分析」開始新的樣本分析。

## 版本紀錄

### v5.0 — 2026-06-11

- GeneBe ACMG 第二意見改用院內本機 GeneBe 資料庫離線查詢，取代原本的 GeneBe 線上 API：三級分析不再需要 GeneBe 帳號或對外網路，也沒有 API 速率限制。
- 查不到的變異（資料庫未收錄，如罕見的 novel indel）就不顯示 GeneBe 第二意見，pipeline 自身的 ACMG 分類仍會顯示。
- 修正未設定 `NGS_UI_HOME` 的 standalone checkout 啟動路徑判斷，避免 restart 時誤把 checkout 上層目錄當成資料根而造成服務啟動失敗。

### v4.9 — 2026-06-10

- SNV/Indel 卡片新增 GIAB stratification 標籤：變異若落在 homopolymer、tandem repeat、segmental duplication、low mappability、GC 極端或其他困難區域，會在標籤列顯示對應的琥珀色 badge，提醒該位點 short-read 判讀較不可靠。
- 三級分析尾段新增 GIAB stratification 標註步驟（落在哪些困難區寫入 `GIAB_STRATA` 欄），純顯示用途，不影響 tier 排序或診斷報告。

### v4.8 — 2026-06-10

- In-house 三級分析預設 Sample ID 改為加上 `-nckuh` 後綴，並在 pipeline output 端產生對應的 `-nckuh` 目錄，和 DRAGEN 的 `-dragen` 分流一致。
- 後端也會強制補上新 job 的來源 suffix，避免舊前端快取或手動輸入造成 output 仍落在無 suffix 目錄。
- 三級分析 Nextflow 進度條改為單調遞增，避免 batch job 中較早 process 晚更新時造成百分比倒退。
- SNV/Indel 同基因搜尋與 More 內的 MANE transcript 顯示改為優先使用有實際 HGVS.c/p 的 MANE_SELECT `NM_` transcript，並容忍不同 TSV 引號格式的 `MANE_ALL` 欄位。
- 統一 SNV filter 計數字體，並讓 DRAGEN caller badge 使用和 DV/HC 一致的標籤樣式。

### v4.7 — 2026-06-09

- 三級分析 output 改以 `-dragen` / in-house suffix 目錄分流，避免同名 DRAGEN 與 in-house sample 混用。
- Nextflow reuse 改為檢查 SNV、Mito、CNV、SV 與已勾選時的 PGx/PharmCAT 輸出；缺檔會用 `-resume` 補跑。
- 三級分析 log 改用行尾時間戳，並記錄每個 Nextflow process 的完成比例與 elapsed time。
- SNV 同基因搜尋與 MANE 詳細資料優先顯示 RefSeq `NM_` transcript。

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
