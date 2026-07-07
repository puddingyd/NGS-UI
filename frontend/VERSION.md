# 歡迎使用 NGS 分析平台

成大醫院基因醫學部 NGS 分析平台是院內 NGS 變異判讀工具，用來整合二級 pipeline 產出的 SNV/Indel、CNV/SV、Mitochondria、STR 與 PGx 變異資料，進行三級分析 pipeline，並協助醫師及生物資訊工程師載入個案、檢視變異、標記變異、撰寫判讀意見，並輸出診斷報告。

請先在上方搜尋既有個案，或使用「載入新個案」登錄 pipeline 已輸出的樣本，或使用右上方「三級分析」開始新的樣本分析。

## 版本紀錄

### v6.2 — 2026-07-05

- SNV/Indel 新增「本院 AF（AF_nckuh）」：顯示該變異在院內 WGS cohort 的等位基因頻率。
- 個案清單新增 Phenotype 欄位，顯示該個案用於分析之 HPO/panel。
- 主畫面 SNV/Indel filter 改為預設隱藏 `VAF < 0.2` 或 `zygosity=ref` 的點位；勾選 `VAF < 0.2 / zygosity=ref` 後才顯示。
- 「載入新個案」的未登錄個案欄位可輸入 LIS ID、source sample、姓名或 MRN 搜尋；清單前端快取一天，可按「更新清單」手動重抓。
- 載入或重新載入個案時，SNV/CNV/Mito/STR 分頁與 gene search 輸入會回到預設狀態。
- SNV/Indel 卡片即使沒有 OMIM disease/ID，只要有 gene symbol 仍會顯示 OMIM 按鈕，連到 OMIM geneMap 搜尋頁。
- 服務端限制 SNV payload cache：大型完整 annotation TSV fallback 不再常駐 uvicorn 記憶體，降低長時間使用後的記憶體累積。

### v6.1 — 2026-06-29

- Secondary findings 新增血脂相關基因、腫瘤相關基因與中風相關基因區塊。
- 新增「匯出健檢報告」：匯出前可選擇報告項目，預設勾選 ACMG SF 與藥物基因體學。
- 健檢報告中的 ACMG SF 依血脂、腫瘤、心血管、代謝與內分泌、麻醉用藥風險及其它基因分組，疾病名稱使用 ClinGen / ACMG SF v3.3 表格。


### v6.0 — 2026-06-19

- 新增 STR 分析區，依 Pathogenic、Intermediate / Borderline、Normal / No threshold 分頁顯示。
- 新增 Pharmacogenomics 區塊，分成 Clinically actionable、Routine / negative screens 與 Additional PharmCAT genes。
- SNV/Indel variant card 支援同一 genomic variant 的多 transcript annotation；使用者可在卡片上切換 transcript。


### v5.1 — 2026-06-15

- 新增「二級分析」工具，可搜尋 WES/WGS FASTQ、建立 samplesheet，並產生可貼到 DGX2 執行的 tmux 指令。
- 二級分析支援批次加入與「加入同批全部樣本」；WGS 一律使用 lane FASTQ 並輸出 lane 欄；reanalysis 可改 samplesheet 的 Sample ID。
- 輸入臨床表徵工具新增 Clinical presentation 欄位，依病歷號自動儲存並帶入主畫面。

### v5.0 — 2026-06-14

- 新增 Secondary findings 顯示：ACMG SF、Proactive 與 Carrier screening 區塊，篩選各 panel 內符合 ClinVar P/LP 或 ACMG P/LP 之位點，Clinvar P/LP 之位點自動帶入報告，其他 ACMG P/LP 之位點需手動確認。
- 輸入臨床表徵工具，新增 `WES dead zone` 與 `WGS dead zone` 按鈕，可查詢 HPO / panel 之 Dead zone 及 CDS dead percentage。
- 修改主畫面 Dead zone 排序規則為先依 CDS dead percentage 分區（70-100%、50-70%、30-50%、<30%），同一區間內再依 pheno score 由高到低排序。

### v4.9 — 2026-06-12

- 基因套組的基因名稱套用 HGNC current name。
- 輸入臨床表徵工具新增 HPO / panel 之基因清單分頁，可查看、篩選各 HPO term 或 panel 的基因清單；另可查特定基因出現在哪些 HPO terms / panels。
- 主畫面 Dead zone 卡片新增臨床門檻下的 CDS dead percentage，並依比例由高到低排序。

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
- 新增 Dead zone 提醒，依目前 HPO / panel 基因列出 coverage 低於判讀門檻的 exon；WES 使用 20X，WGS 使用 DRAGEN 10X。
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
