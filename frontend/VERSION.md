# 歡迎使用 NGS 分析平台

成大醫院基因醫學部 NGS 分析平台是院內 NGS 變異判讀工具，用來整合二級 pipeline 產出的 SNV/Indel、CNV/SV、Mitochondria、STR 與 PGx 變異資料，進行三級分析 pipeline，並協助醫師及生物資訊工程師載入個案、檢視變異、標記變異、撰寫判讀意見，並輸出診斷報告。

請先在上方搜尋既有個案，或使用「載入新個案」登錄 pipeline 已輸出的樣本，或使用右上方「三級分析」開始新的樣本分析。

## 版本紀錄

### v6.11 — 2026-07-18

- 登錄新個案後會同步切換上方個案搜尋框；登錄狀態不再顯示 job ID。Exomiser/LIRICAL 改為至少有一個 HPO term 才排入，只有 panel 時仍會更新 pheno score，但不再產生 prepare-vcf 失敗訊息。
- 主畫面 Secondary findings 精簡為 ACMG SF、中風相關基因與 Carrier screening（PGx / PharmCAT 保留）；同一 variant 命中多個 panel 時仍會在每個分類分別顯示，但所有卡片共用同一筆勾選、ACMG、Comment、transcript 與疾病選取資料。既有勾選衝突以明確取消為準，健檢匯出時才合併去重。
- 健檢 DOCX 匯出選項精簡為 ACMG 疾病風險基因、中風相關基因、帶因者篩查與藥物基因體學；所選疾病 panel 點位先合併去重，再統一套用疾病風險／帶因者兩類規則。基因清單移除類別序號，附錄改為另起新頁、標題置中並留一空行。
- 三級分析畫面的黑底 log 可視高度調整為原本約 1.22 倍。

### v6.10 — 2026-07-17

- 修正三級分析清單刪除新 unified-path sample 時，檔案已刪除但 cache cleanup 因舊變數名稱回傳 500 的問題；確認視窗也會顯示新的 datalake 路徑與 legacy fallback 路徑。
- SNV/Indel 卡片的 Score、ClinVar、ACMG 與所有 in-silico tools 新增 `ⓘ` 註解；可查看計分公式、ClinVar 版本、ACMG 的 manual/GeneBe/in-house 來源比較，以及各工具的 PP3/BP4 threshold、文獻連結與 PMID。
- P-KNN 註解新增 LLR 的 ±1/±2/±4 門檻，空 evidence 改顯示 Uncertain；DANN 與 LOFTEE 不再套色，並精簡 Score、ClinVar、LOFTOOL 的註解文字。
- AlphaMissense、ESM1b、VARITY_R、BayesDel、REVEL、SpliceAI、PhyloP、GERP、SIFT 改依正式 calibration/ClinGen 建議套色；其他尚無通用 PP3/BP4 calibration 的工具會清楚標示為模型 cutoff 或 contextual evidence。SIFT 多 transcript 分數改取較低、較 deleterious 的值。
- 三級分析的 Extra VEP 現在會在 dbNSFP 實際提供 `REVEL_score` 時補入 REVEL，並與 MetaRNN、可選的 SpliceAI 一樣存入 sparse overlay；舊 dbNSFP 或未跑 Extra VEP 的個案保持缺值，不會被誤判為低分。
- SNV/Indel tier 簡化為 `1A / 1B / 1C / 2`；`1C — Predicted suspect` 除 ACMG points ≥4 外，會依 Core（AlphaMissense、BayesDel、Pangolin）與 Extra-VEP（REVEL、SpliceAI）門檻納入 reviewer 候選，但不顯示額外 trigger badge、不改 ACMG 分數。原 ClinVar P/LP 0★/conflicting tier 移除，其餘歸入 `2 — Other`。
- In-silico 顯示順序改為 P-KNN、AlphaMissense、Pangolin、REVEL、SpliceAI、ESM1b、VARITY_R、BayesDel、MetaRNN、DANN、PhactBoost、PhyloP、GERP、SIFT、LOFTOOL；前三個有值項目直接顯示。Extra VEP 預設 dbNSFP 改為 `biotools/dbnsfp/dbNSFP5.3.1a_grch38.gz`，SpliceAI 仍使用 server 原有 `biotools/spliceai/` 路徑。
- 三級輸出可在 `03_acmg/{source}.annotation_versions.json` 記錄 ClinVar release date；缺少 sidecar 的舊個案會明確顯示版本未提供，不會用檔案日期猜測。

### v6.9 — 2026-07-16

- 三級分析輸出統一到 `/home/datalake_Intermediate/pipeline/tertiary_output/{sample}/`：Nextflow 保留 03–07，UI 的判讀狀態與衍生檔集中在 `08_postprocessing/`。
- UI 直接讀取 03–07 的 SNV、CNV/SV、STR 與 PGx，不再永久複製大型 `snv_indel.annotated.tsv`。GeneBe、GIAB、院內 AF、MANE/extra-VEP 改存稀疏 SQLite overlay，再合併產生 review TSV 與 gene search 結果。
- 新增安全的舊個案遷移工具：支援 dry-run、逐案 canary、批次遷移與 marker-only rollback；舊資料會保留到人工確認完成。

### v6.8 — 2026-07-16

- SNV/Indel 主畫面移除重複的 gnomAD AF filter，直接顯示 review TSV 內通過其他顯示條件的候選點；`impact=MODIFIER` 仍預設隱藏一般 MODIFIER，但 ClinVar P/LP 點位會自動 rescue 顯示。Gene search 的 gnomAD AF filter 維持不變。

### v6.7 — 2026-07-14

- 健檢 DOCX 依最新 Word 註解調整段落與附錄：第一類標題後不留空行、第一類結果與第二類標題之間保留一行，兩個附錄之間保留一行；附錄標題改為「變異位點參考資料」。PGx 官方資訊改以 CPIC 為第一項並更新導讀文字，檢測限制只移除被標示的「ACMG SF 或」；ASCII 表格不再套用會在 Word 顯示黑色方塊的 keep-with-next 段落設定。
- 個案清單的「刪除」改為清除已載入個案的註冊/報告狀態，保留本機三級輸出檔案；刪除後 sample 會從個案清單消失並回到「載入新個案」的未登錄清單，不需要重新跑三級。
- 健檢 DOCX 副標題改為「基因醫學部基因檢測分析研究報告」；ACMG 主文分成「與疾病風險相關」與「符合帶因者狀態」兩類，並依遺傳模式、同基因變異數與基因型自動分類；六組 ACMG 基因分類與完整基因移到第五節基因清單。
- 健檢 PGx 主文改為「用藥建議概覽、藥物建議摘要、基因型與表現型」順序；摘要及完整建議只列實際驅動處置的基因，避免同一藥物把正常表型基因誤列為調整依據。CPIC 與 FDA 官方查詢網址及用藥警語移至 PGx 主文末端並顯示完整網址，完整用藥建議則集中到報告末端附錄。
- 健檢 ACMG 共用警語移到第一、第二類之前，變異表格的 ACMG/AMP 分級維持英文，表後說明與附錄參考資料使用中文。WGS 健檢報告平均深度更新為 27X，限制段落新增 CNV 未涵蓋範圍及 CYP2D6 藥物基因體學會專一納入 CNV 的例外說明；附錄依序列出變異位點參考資料及完整用藥建議，不再使用附錄一、二編號。

### v6.6 — 2026-07-12

- 輸入臨床表徵工具的 phenotype 預覽改為即時更新：載入既有病例即可顯示，新增或移除 HPO/panel 會同步重繪，切換到其他病例時不再殘留上一位病人的預覽。
- 健檢 DOCX 標題改為「基因醫學部基因檢測檢驗分析研究報告」，檢驗套組改寫為「ACMG疾病風險基因」，匯出健檢報告選單顯示「ACMG 疾病風險基因」；ACMG 警語使用次發現基因清單第 3.3 版，AR 單一 heterozygous 變異在第 1 點補上「符合報告條件之變異」帶因者狀態說明，AR 多變異改用相位未確認、建議家族成員檢測與雙等位基因致病型態句型，X-linked 女性單一變異改為受 X 染色體失活型態、家族史及疾病表現範圍影響的說明，建議句改為「遺傳諮詢或門診相關專科」，PGx 用藥建議分類後的預設句改為「其餘未列之藥物」。
- 健檢 ACMG SF 參考資料順序改為跟上方點位顯示順序一致，同一基因多位點會連續列出。

### v6.5 — 2026-07-12

- 二級分析的 WGS FASTQ 清單改為每個主 sample 只顯示一次，並標示 lane 與 FASTQ 檔案數；建立 samplesheet 時仍會展開成每 lane 一列，且排除 merged FASTQ。WES 行為不變。

### v6.4 — 2026-07-09

- 健檢報告匯出選單改用「可採取醫療處置之疾病風險基因」；ACMG SF 同基因多位點合併顯示，並修正 intron、中文分級與女性性聯遺傳說明。
- 健檢 PGx 改為完全使用 PharmCAT JSON：narrative 依 CPIC Strong/Moderate 非標準處置與 FDA Section 1 自動產生，完整表格同時保留 CPIC 與 FDA 原文。
- PGx 報告改為重點摘要、基因型與表現型、用藥建議分類、完整用藥建議四段式；完整建議改以藥物為主，CPIC/FDA 建議合併於同一藥物列，並依藥物英文名稱 A-Z 排序。
- 健檢報告調整 ACMG SF 參考資料位置；PGx 基因型表支援 MT-RNR1 由 TSV 補值，phenotype 空值會改用 allele function 或顯示 `No phenotype assigned`；PGx 重點摘要不列 FDA-label-only 項目，完整建議改為句尾括號標示 CPIC/FDA 來源，並清理原文引號與省略符號。
- PGx 完整用藥建議表格加寬「基因與表型」與「CPIC/FDA 建議」欄，讓總寬與基因型表格對齊。
- SNV 卡片與 DOCX 將 EXON 的 `. / - / NA` 視為缺值並改讀 INTRON，避免 intronic variant 顯示成 `exon.`。
- 健檢報告保留原報告名稱與大格式；檢驗套組依 UI 勾選動態組合。ACMG SF 維持原有變異表格，更新 v3.3/ClinVar 警語、AD/AR/X-linked 三點說明及集中式檢測限制。
- PGx 健檢結果新增中文 narrative summary，整合 CPIC Strong/Moderate 與 PharmCAT JSON 的 FDA therapeutic-management 資訊；完整英文建議仍保留於下方結果。
- 健檢報告 PGx 回報範圍更新為 21 個 CPIC Level A genes；IFNL3 從重點摘要、基因型與表現型、用藥建議分類、完整用藥建議及末尾基因清單排除，標題與參考依據不列 FDA；validation scope 待確效完成後再縮減。
- DRAGEN 三級分析新增 ploidy VCF sidecar 複製；基本資料性別欄會核對 M/XY、F/XX，一致顯示綠底，不一致或非典型性染色體 call 顯示紅底及 ploidy VCF 結果。
- 疾病清單固定將 OMIM disease 排在補充來源之前；OMIM 列改用稍深灰底，並以細分隔線區隔 GenCC / ClinGen 補充疾病。

### v6.3 — 2026-07-08

- SNV/Indel 疾病清單可合併 optional GenCC / ClinGen / MONDO gene-disease SQLite index，在 OMIM 未列 disease 或有額外非重複疾病時補充顯示；OMIM description / synopsis 仍即時取自最新版 OMIM.xlsx。
- OMIM disease 若只有單行、尚未補 description / synopsis，卡片 summary 後會以小 `*` 提醒補齊；此標記只影響 UI，不進報告或列印輸出。

### v6.2 — 2026-07-05

- WGS 固定套組新增血液科 Lymphoid Neoplasm Panel 與 Myeloid Neoplasm Panel，SNV/indel、CNV、STR、Mitochondria 欄位合併為單一基因清單。
- SNV/Indel 新增「本院 AF（AF_nckuh）」：顯示該變異在院內 WGS cohort 的等位基因頻率。
- 個案清單新增 Phenotype 欄位，顯示該個案用於分析之 HPO/panel。
- 主畫面 SNV/Indel filter 改為預設隱藏 `VAF < 0.2` 或 `zygosity=ref` 的點位；勾選 `VAF < 0.2 / zygosity=ref` 後才顯示。
- 「載入新個案」的未登錄個案欄位可輸入 LIS ID、source sample、姓名或 MRN 搜尋；清單前端快取一天，可按「更新清單」手動重抓。
- 載入或重新載入個案時，SNV/CNV/Mito/STR 分頁與 gene search 輸入會回到預設狀態。
- SNV/Indel 卡片即使沒有 OMIM disease/ID，只要有 gene symbol 仍會顯示 OMIM 按鈕，連到 OMIM geneMap 搜尋頁。
- 服務端限制 SNV payload cache；gene search 缺 index 時改為串流掃描大型完整 annotation TSV，不再整份物件化並常駐 uvicorn 記憶體。

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
