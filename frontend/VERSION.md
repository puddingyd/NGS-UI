# 歡迎使用 NGS 分析平台

成大醫院基因醫學部 NGS 分析平台是院內 NGS 變異判讀工具，用來整合二級 pipeline 產出的 SNV/Indel、CNV/SV、Mitochondria、STR 與 PGx 變異資料，進行三級分析 pipeline，並協助醫師及生物資訊工程師載入個案、檢視變異、標記變異、撰寫判讀意見，並輸出診斷報告。

請先在上方搜尋既有個案，或使用「載入新個案」登錄 pipeline 已輸出的樣本，或使用右上方「三級分析」開始新的樣本分析。

## 版本紀錄

### v9.4 — 2026-09-02

- 修正大量 case 執行三級分析時 VEP 進度長時間停在約 32%：Nextflow 在 batch channel 尚未收齊時可能先輸出動態 `1 of 1`、`2 of 2`，現在改以本批 sample 數作為進度分母，並以完成記號／完整批次判斷 process 結束，讓 VEP 權重隨完成 case 數逐步增加；modal 同時顯示目前 process 與 `完成數/總數`。

### v9.3 — 2026-09-01

- 新增唯讀 `scripts/audit_pgx_report_inputs.py`，可在正式機批次掃描指定 prefix 的 `07_pgx` PharmCAT JSON／PGx TSV，核對固定 21-gene 健檢表格、HLA risk allele copy 狀態與 MT-RNR1 TSV fallback，並輸出實際 allele 組合、No Result、所有 source diplotypes 及未定相 non-reference variant calls；不修改任何樣本輸出。
- 修正 SNV/Indel 卡片初次載入時 in-silico predictor 高度判斷偶爾過度展開：現在只以 LitVar2 實際可見內容（關閉的 details 只計 summary）決定可多顯示的列數，並在 LitVar2 尺寸穩定或改變後自動重新校正；保留至少六項直接顯示及按 More 展開其餘內容的行為。
- 修正「輸出 PDF」的 LitVar2 PMID 顯示：列印副本移除外部連結前，會先把最多五筆 `PMID:...` 與 `and N others` 轉成可換行的純文字，不再因連結被刪除而只剩一串逗號；`Merged from N LitVar2 records` 摘要仍會保留。
- 健檢 DOCX 依 TITAN v4 標註版調整 PGx 主文：移除用藥建議概覽、藥物建議摘要、表型欄與動態藥物數量文字，固定 21 genes 統一為「基因／檢測項目、Allele 1、Allele 2」三欄。一般 star allele／haplotype 與已分相 DPYD 各 allele 分欄；DPYD 多筆未定相 `sourceDiplotypes` 則跨兩個結果欄完整列出所有 named variants、het/hom 並明示「相位未定」，不再只取第一筆。CFTR Reference 不再重複長的 ivacaftor-response 功能名稱；ABCG2／VKORC1 以 rsID＋HGVS 搭配 transcript-oriented base 的 Reference／Variant 表示。MT-RNR1 與來源只有單一 allele 的 G6PD 不複製成兩份，第二欄顯示 `N/A`；完整用藥建議仍置於附錄。WGS 檢測方法另新增短讀長偵測限制與 phasing／haplotype 不確定性兩點，後續題號順延。
- 健檢 DOCX 的 HLA 結果改為逐 allele 呈現：保留 HLA-A／HLA-B 的 Allele 1／Allele 2 分型，並分別加列 `HLA-A*31:01` 與 `HLA-B*15:02`、`HLA-B*57:01`、`HLA-B*58:01`；四個篩檢子列的第一欄縮排兩格以呈現 HLA-A／HLA-B 展開階層。兩欄各自顯示 Positive／Negative，因此 homozygous risk allele 會顯示兩個 Positive；無可靠分型或明確狀態時顯示 `No Result`，不由 aggregate positive 推定兩個 copies。

### v9.2 — 2026-08-31

- 二級分析產生的 DGX2 tmux 啟動腳本會在載入環境後設定 `umask 0002`，讓共用 `dgm_nckuh` 群組的新 Nextflow work 目錄保持可共同清理；DGX2 清理指令也會檢查刪除狀態與目錄殘留，只有確實清空才顯示「已清理」。
- Manual ACMG/AMP modal 上方的 Manual、ERepo、GeneBe、In-house 來源摘要，所有可辨識 criteria 都會明確顯示 evidence strength，例如 `PM1_Moderate`、`PP3_Supporting`、`PVS1_Very_strong`；Manual 的 Global reusable 摘要也使用相同格式。這只調整畫面文字，不改原始 criteria 儲存格式、Apply 行為或 ACMG points 計算。
- 主畫面 Comment textarea 改為與 Clinical presentation／Genetic counseling 相同的自動增高行為：載入既有內容、展開卡片及輸入文字時，都會依完整內容高度展開，不再在欄位內顯示垂直捲軸。

### v9.1 — 2026-08-29

- SNV/Indel 卡片展開 `Transcripts` 後，可直接點選任一列切換卡片顯示的 transcript；目前選用列會顯示藍底與「目前顯示」，並與 HGVS 旁的小三角、DOCX、PDF 及個案清單摘要共用同一筆自動儲存選擇。`BEST_CONSEQUENCE` transcript type 也補上紫色 badge 底色，與其他 type 樣式一致。
- `1000G EAS` 移回 SNV/Indel 卡片最右側 AF 欄，排在 `AF`、`AF_eas`、`AF_nckuh` 後作為第四列；維持按 `More` 才顯示，不再出現在 info grid 下方。
- 修正切換 SNV transcript 後，疾病相關基因警示仍沿用原始 gene 的問題：警示現在依 reviewer 選定 transcript 的 gene 與其 disease-associated 狀態判斷。警示需要顯示時改用 DOCX「本次檢測基因包括」名稱，不再顯示內部段落編號 `§五.4`。
- 主畫面 SNV/Indel 卡片的 zygosity 統一使用 `het`／`hom`／`hemi` 縮寫；`hom` 與 `hemi` 會套用和 Likely pathogenic 相同的淺紅底，並相容 TSV 輸出的 `homozygous`／`hemizygous` 長寫法。此提示只影響顯示，不改 tier、filter 或報告判讀。

### v9.0 — 2026-08-27

- 主畫面 Patient phenotype 只有儲存 `default` 分析組合時，才會把目前 HPO、panel 與 weight 原子同步到病人層級的 `patient_phenotype/{MRN}_phenotype.txt`；後續新增的其他分析組合只寫入 sample-owned analysis version，不覆蓋病人主檔。載入新個案仍會保存最後採用的 default 快照。從個案清單移除 sample 時不刪此檔，因此相同 MRN 日後使用原檢體或新檢體重新登錄都能自動帶回。default 全部清空時會保留 header-only 明確空值，避免舊 LIS-specific phenotype 復活。
- phenotype lookup 全面改為 MRN-only 主檔優先，舊 `{LIS}_{MRN}_phenotype.txt`／`{LIS}_phenotype.txt` 只作非破壞 fallback；`/phenotype/` 有 MRN 時也固定寫 MRN-only，同一病人的不同檢體共用 phenotype。
- 「輸入臨床表徵」頁的「載入既有資料」右側新增同樣式 `EMR` 與 `GC紀錄`：EMR 直接開啟院內 HIS，GC紀錄經登入保護的 consultation API 即時抓取並以唯讀 modal 顯示看診日期、reason、diagnosis 與 counseling record。

### v8.9 — 2026-08-14

- 主畫面 SNV/Indel 的 in-silico predictor 欄改為每個工具名稱與結果同行顯示，結果依最長工具名稱統一對齊；固定顯示至少六個工具，若 ClinVar／ACMG／LitVar2 欄較高，會以固定行距自動多顯示可容納的 predictor，超出 LitVar2 最後一行後才收進 More。既有字體、大小、顏色與結果底色不變。

### v8.8 — 2026-08-13

- Documents modal 的上傳區新增檔案拖放操作，保留批次選檔與剪貼簿截圖；文件清單新增「下載全部 ZIP」，登入後依 MRN 將目前全部文件以顯示檔名即時串流打包，不在伺服器留下額外 ZIP。
- 「載入新個案」的 MRN 欄旁新增直接開啟院內 HIS EMR 的 `EMR` 按鈕，會隨手動輸入或未登錄個案清單帶入的 MRN 更新；原有 `EMR 同步` 功能維持不變。

### v8.7 — 2026-08-12

- OMIM 疾病 curator slots 從 `Disease1..5` 擴充為 `Disease1..16`：SNV 卡片、GenCC／ClinGen／MONDO evidence 合併、reviewer 疾病勾選、個案清單摘要、TXT/PDF 與診斷 DOCX 全程支援第 6–16 欄；OMIM ID 優先、gene symbol fallback 的變異命中規則不變，舊五欄 workbook 與既有 report state 可繼續使用。
- 只有 `Disease1..16` 全數為空時，才從多行 `OMIM_disease` 依序合成最多 16 個單行 fallback；只要已有任一 curator slot 就不自動補其餘空欄，避免 grouped rich description 與原始 disease 行錯位或重複。

### v8.6 — 2026-08-12

- 主畫面 Patient phenotype 與「載入新個案」的 HPO 搜尋改為比照「輸入臨床表徵」工具：支援名稱與 synonym 的拼字容錯，並在每個命中 term 下方顯示可直接選取的上一級 HPO term；滑鼠與鍵盤選取行為一致。

### v8.5 — 2026-08-11

- GPN-MSA 改為所有三級 case 都嘗試、但不再作為完成 gate：缺 score DB、`.tbi`、tabix 或查詢失敗時保留可用的 review TSV，在 job log 顯示 warning，並於 review manifest 記錄明確 `complete`／`skipped_*`／`failed` 狀態，不會讓整個三級分析失敗。
- 新增 `scripts/backfill_gpn_msa.py`，可指定完整 UI case ID、source sample ID（自動解析對應的 `-nckuh`／`-dragen` case）或以 `--all` 對舊 case 補 GPN-MSA；已有 review TSV 時只更新該衍生檔與 manifest，不修改或重掃完整 raw TSV。
- 首頁流程圖的 SNV/Indel Prioritization 改為「可編輯的 Filter → Category」左右分區，Category 直向列出 1A / 1B / 1C / 2 - Other，縮減 tier 的橫向空間。
- TITAN-WGS 個案的「報告」與「分析」標題旁新增 Secondary findings 小字說明，區分自動帶入的 ClinVar P/LP 與需勾選後才進報告的 Predicted suspect；WES/WGS 不顯示。

### v8.4 — 2026-08-11

- 二級分析產生的 Nextflow 指令改依定序類型選擇 DGX2 profile：WES 使用 `-profile dgx`，WGS 維持 `-profile dgx_single`；其餘分析參數不變。

### v8.3 — 2026-08-11

- 主畫面與「輸入臨床表徵」頁的 Clinical presentation 右側新增 `Documents`：登入後可依 MRN 共用上傳、下載、修改檔名與刪除病歷附件，也可直接貼上剪貼簿截圖。支援 PDF、JPG、PNG、TIF/TIFF；圖片可在 modal 預覽，多頁 TIFF 可翻頁且保留原始下載檔。附件串流寫入獨立 `patient_documents/`，不設固定單檔上限但保留磁碟安全空間；MRN 修改時會連同附件與 phenotype sidecar 安全搬移，遇到既有或共用病人資料時不自動合併。

### v8.2 — 2026-08-11

- 二級分析的 WES / WGS FASTQ 索引改為每天台北時間 02:00 自動更新；人工更新按鈕、首次建立與索引過期時的 fallback 仍保留。

### v8.1 — 2026-08-10

- 所有三級分析 case 的 compact SNV/Indel review TSV 新增固定 GRCh38 GPN-MSA pre-computed score；只以 tabix 查 final review SNVs，不改 `03_acmg` raw、不寫 sparse overlay，也不影響 tier、ACMG、filter、排序或報告。卡片將 GPN-MSA 放在 SpliceAI 後，低分代表較 deleterious，`≤ -7` 只作 contextual tail 提示，未校準成 PP3/BP4，且不產生 benign evidence。
- SpliceAI、GeneBe live API fallback 與 LitVar2 改為共用 review TSV 的最終保留規則：primary/non-mito、WES DP ≥20、baseline ClinVar P/LP rescue、rare/unknown AF + candidate BED，以及 weekly ClinVar `UP_TO_PLP` / `DOWN_FROM_PLP` rescue。GeneBe 本機主 DB 仍照舊查完整 working TSV；SpliceAI 仍只在勾選 Research-only 時執行。
- GPN-MSA reference 簡化為 `$NGS_UI_HOME/biotools/gpn_msa/scores.tsv.bgz` + `.tbi` 單一固定部署，不設 release/current 目錄或自動更新 timer；新三級 job 缺 score、index 或 tabix 會明確失敗，既有 case 則可在載入時依 review manifest 懶重建補分。

### v8.0 — 2026-08-04

- SNV/Indel 疾病清單改以 OMIM `Disease1..5` slot 為權威列：即使多筆疾病共用同一 phenotype MIM，仍依 Excel 原順序完整保留名稱與說明，不再被 GenCC / ClinGen / MONDO 合併覆寫或重複成同一列。同一 MIM 的 supplemental evidence 會附加到所有匹配的 OMIM slots，且不再額外重複顯示該 supplemental disease；只有 OMIM 完全沒有的 association 才列在下方補充區。
- LitVar2 外部連結圖示改為固定顯示：唯一命中時連到該 variant 結果頁，`No reference`、`NA` 等其他狀態連到 LitVar2 首頁。多筆命中改顯示可展開的 `Ambiguous match (N records)`，完整列出每筆 LitVar ID、rsID/gene/HGVS、PMID 數、前五篇 PubMed 與各自的 LitVar2 結果頁，不再只顯示無法追查的 ambiguous 文字。
- LitVar2 重複 records 會依 ClinGen CA 或 rsID + gene + normalized HGVS 聚合成 logical variant，完整 PMID 聯集去重後只顯示一次，來源收在 `Merged from N LitVar2 records`；只有仍不同的候選才顯示 `Ambiguous match (N variants)`。既有 v1 SQLite 可繼續讀取，手動／每月更新會直接用本地 bulk 升級 v2。
- 左側「基本資料／報告／分析」與子區導覽改為直接跳到目標高度，不再顯示快速滑過整頁的 smooth scrolling。
- 新版三級 v3.6 的 ClinGen Evidence Repository（ERepo）會在 SNV/Indel 卡片的 ClinVar 與 ACMG 之間顯示「分類（推導分數）」；點擊可在 structured ACMG modal 比較 VCEP criteria、Apply 至 manual 判讀，並標示這是 ClinGen experts 的評估。
- 三級分析選項補上 `Research-only`：勾選時 Nextflow 加入 `--academic_dbnsfp true` 使用 dbNSFP 5.3a，顯示 REVEL、MutPred2、VEST4、CADD；post-processing 只補 SpliceAI，不再重跑 dbNSFP。未勾選時維持 dbNSFP 4.9c 並略過 SpliceAI。
- In-silico 卡片順序更新為 P-KNN、AlphaMissense、Pangolin、REVEL、SpliceAI、ESM1b、VARITY_R、BayesDel、CADD、DANN、MutPred2、VEST4、PhactBoost、PhyloP、GERP、SIFT、LOFTOOL；新工具依文獻 PP3/BP4 calibration threshold 套色並附註。
- Nextflow ClinVar 基準固定為 2026-07-20；卡片標題固定顯示 `ClinVar (2026-07-20)`，內容也維持 baseline 分類，不再顯示 `ⓘ`。標題旁的外部連結圖示優先用 baseline Variation ID 開啟 NCBI ClinVar；baseline 沒有紀錄時 fallback 最新版 Variation ID，兩者都沒有就隱藏。每週三 04:00 更新的最新版只用來偵測箭頭；hover 箭頭才顯示最新版分類、review status 與日期。只有非 P/LP→P/LP 顯示紅色上箭頭、P/LP→非 P/LP 顯示綠色下箭頭；LP→P、LB→B、conflicting→VUS 不標。診斷與健檢 DOCX 的日期、分類及健檢候選集合仍固定使用 20260720 基準。

### v7.9 — 2026-08-03

- LitVar2 維持三級 filtered post-processing；有完成 marker 的個案在 gene search 時，會把尚未註解的 variants 批次查詢本地 SQLite，並依目前 LitVar2 DB 版本快取 hit 與 no-match，不需要網路也不逐張卡片查詢。
- 所有 SNV/Indel variant 卡片的 LitVar2 標題旁新增重新整理圖示，可強制用最新本地資料庫重查單一 variant；舊個案也可只更新選定點位，其餘仍顯示 `NA (請重跑三級)`。

### v7.8 — 2026-08-02

- 三級分析 modal 在「加入批次」左側新增「加入同目錄全部檢體」；選定一個 DRAGEN 或 NCKUH VCF 後，可將同一來源目錄的所有檢體依預設輸出 ID 一次加入批次，並自動排除重複 sample/path。

### v7.7 — 2026-07-31

- 三級分析的 DRAGEN / NCKUH VCF 索引改為每天台北時間 02:00 自動更新；人工更新按鈕與索引過期時的 fallback 仍保留。
- 分析區新增 `AF_NCKUH ≥ 0.05 & AC ≥ 50` 顯示開關，預設不勾選；SNV/Indel 與 mtDNA 中同時達到兩項門檻的 local common variants 預設隱藏，勾選後可展開查看。
- ClinVar non-conflicting P/LP、有效 ACMG P/LP、reviewer 已標記 `1/2/C` 的點位仍保留；mtDNA 另保留 MITOMAP pathogenic/reported 與手動 P/LP。Homozygous／hemizygous 及同一 recessive gene 的第二個 candidate 不單獨構成 rescue。
- 分析區 filter 名稱改為「疾病相關基因」與「臨床表現相關基因」；兩者後方數字會先套用其他啟用中的顯示條件，再呈現目前可顯示的數量。
- 統一分析區所有 filter 外框的最小高度，最右側 OMIM checkbox 現在與其他 filter 等高。

### v7.6 — 2026-07-30

- SNV/Indel 卡片在 ClinVar/ACMG 下新增 LitVar2 文獻列，顯示 bulk 版本日期；前五個 PMID 可直接開啟 PubMed，超過五篇以 `and N others` 連到該 variant 的 LitVar2 結果頁，標題旁的外部連結圖示也可固定開啟該結果頁，標題本身維持純文字。
- 三級 post-processing 只對 review filter 留下的 genomic variants 查本地 LitVar2 SQLite；先精確比對 rsID，缺 rsID 才用 gene + HGVS，結果只供 reviewer 參考，不影響 tier、ACMG、排序或報告。舊個案重跑三級後才補註。
- LitVar2 官方 bulk JSON 與 slim SQLite 放在 `NGS_UI_HOME/biotools/litvar2/`；三級 modal 新增手動更新按鈕，每月一號另由 systemd timer 背景更新。新資料完整建庫、驗證後才原子切換，更新失敗時保留舊 JSON/DB。
- LitVar2 更新進度移到三級 modal 標題列右上角；舊個案或尚未跑到 LitVar2 annotation 的卡片改顯示 `NA (請重跑三級)`，和真正查無文獻的 `No reference` 明確區分。
- 二級分析 modal 產生的 DGX2 Nextflow 指令，不論單樣本或多樣本都固定使用 `dgx_single` profile。

### v7.5 — 2026-07-29

- SNV/Indel ACMG 改為 structured modal：完整列出 28 個 ACMG/AMP 2015 criteria、原始規範與 ClinGen 後續 guidance，支援啟用／停用及 Supporting、Moderate、Strong、Very strong strength，PP5/BP6 保留可用並標示 ClinGen 停用建議。
- Modal 上方可比較並 Apply Manual、GeneBe、In-house 三個來源；儲存後即時計算 ACMG points、classification、variant/total score、tier 與排序。主卡片只顯示分類、分數及來源，不再直接編輯自由文字。
- Manual ACMG 以 hg38 + normalized chr-pos-ref-alt 跨使用者共用，並保存 append-only revision、reviewer 帳號名稱與來源 sample；case/family-specific criteria 只留在該 sample、不自動套到其他個案。每個 sample 同時保留自己的最終快照，診斷／健檢報告沿用該 sample 畫面結果。
- Manual ACMG 儲存後即由後端重算 tier、total score 與排序；若點位因此移動，畫面會自動切換到新 tier，並把同一張卡片維持在原本視野位置，不另外加醒目框。主卡片仍保留分類套色，並依文獻將 VUS 顯示細分為 low（0–1）、mid（2–3）、high（4–5 points）；黃色為共同底色，low 向右漸綠、high 向右漸紅，正式報告分類仍是 Uncertain significance。
- Manual ACMG 主卡片改成與 ClinVar 一致的單行字體與高度：分類和 points 套色，source 只以旁邊未套色小字顯示。Modal 標題改為 `Manual ACMG/AMP variant classification`；criteria 左側為淡紅 Pathogenic、右側為淡綠 Benign，並依 ACMG 2015 的八種 evidence type 分組排序。上下操作列都提供「關閉／儲存並套用」，計分與 VUS 子分級說明移至 criteria 下方。
- 目前仍為 Causative／Other 的其他個案會在 HGVS 下方顯示 `Observed (N)` badge；點擊可查看 sample ID、status、reviewer 與時間。Candidate、已取消 status 及目前 sample 不計入，Observed 與 Manual ACMG 完全分開。

### v7.4 — 2026-07-28

- Mitochondria 卡片第三行並列 `AF_gnomAD` 與本院 mtDNA carrier frequency（`AF_nckuh`）；後者顯示帶有該 ALT 的院內樣本數、可判讀粒線體樣本數，以及 homoplasmic / heteroplasmic carrier 數。
- 本院 mtDNA AF 直接從 indexed in-house AF VCF 的 chrM slice 載入，只供院內常見／重現變異參考，不影響 Mito tier 或 ACMG 分級。
- 三級分析改為 DRAGEN、NCKUH 各自共用一條 Nextflow cache lineage；同一 sample/input 即使前後放在不同 batch，也可由 Nextflow 沿用有效 cache，只有 hash 條件改變的 task 會重跑。
- 同模式另一個三級分析 job 會先顯示等待共享 cache，等目前的 Nextflow 結束後再啟動；DRAGEN 與 NCKUH 不互相阻塞，後續驗證與 post-processing 也可並行。

### v7.3 — 2026-07-28

- 重跑同一個三級分析 sample 時，不再因正式資料夾已有 SNV/Mito/STR/CNV/SV/PGx 就略過 Nextflow；每次都啟動 `-resume`，由 Nextflow 自行判斷沿用 cache 或重跑。
- Nextflow 先輸出到 job 專屬 staging，完整驗證 00–06、SNV/Mito/STR/CNV/SV 及勾選的 PGx 後才更新正式結果；Nextflow、驗證或 post-processing 失敗時，原本已載入個案與三級結果不受影響。
- 正式切換會保留 reviewer 的標記、comment、個案資料與 phenotype analyses；批次切換失敗可還原，同一 UI/source sample 也不能同時啟動兩個三級分析 job。

### v7.2 — 2026-07-27

- 三級分析的 GeneBe ACMG 第二意見改為主 DB 優先的 hybrid lookup：完整 TSV 先查本機 DB，再由中央 cache 補值；只有 DB/cache 都查不到且符合 SNV review 條件的點才送 live API。
- GeneBe live API 未設定或暫時失敗不會中止三級分析；成功結果會跨個案重用，明確 no-result 預設 30 天後才重查。
- API 新結果另存為正式 DB 相容的 7 欄去重 TSV 與 JSON sidecar，方便後續匯入；既有完整 working TSV 仍只在 post-processing 期間暫存並於結束時刪除，不增加永久大型 SNV 副本。
- PGx 藥物建議摘要直接沿用用藥建議概覽已選出的主要分類與依據，不再用摘要來源子集合重新分類，避免同一藥物在兩處顯示不同處置。

### v7.1 — 2026-07-27

- Ploidy status 改為 reviewer-first 畫面：先顯示 aneuploidy 結論、estimated karyotype／病歷性別核對及需要複核的染色體；完整染色體與原始測量值改為折疊顯示。
- 修正 DRAGEN ploidy 判讀：`FILTER=PASS` 不再被誤當成染色體正常，改以 `ALT=<DEL>/<DUP>` 判定 gain/loss，並分開顯示 PASS、LowQual 與 NCKUH SUSPECT confidence。
- DRAGEN 原始 VCF 沒有 RATIO 時，畫面會以 `DC / autosomeDepthOfCoverage` 顯示明確標註的 derived ratio；所有異常解讀維持 possible／疑似措辭。

### v7.0 — 2026-07-25

- SNV/Indel 卡片新增 `Strand bias PASS / WARN / MANUAL` 標籤；WARN 顯示 FS/SOR，缺少 FS/SOR 的 DeepVariant-only 位點會提示 IGV／人工複核。主畫面、報告區、Secondary findings 與 gene search 共用相同顯示。
- Read support 改為同時顯示 DP、AD、VAF；WGS/TITAN-WGS 總 DP<10 及 ALT AD<10 會分別標紅提醒，WES 仍維持 DP≥20 hard floor。
- DRAGEN 與 NCKUH ploidy VCF 都會帶入三級分析。性別欄固定顯示可點擊的 ploidy 結果，完整視窗列出所有 chromosome 的 FILTER/DC/NDC/RATIO；aneuploidy 或任一核染色體非 PASS 顯示紫色，未輸入／不符性別顯示紅色，符合顯示綠色。
- 新三級個案的 `08_postprocessing` 與 analysis sidecar 全面使用 `{LIS_ID}.<filename>`；舊版未加前綴檔仍可讀，並固定由 prefixed 檔優先。layout v3 marker 最後原子寫入，舊檔不會被自動刪除。

### v6.9 — 2026-07-22

- 首頁歡迎訊息與版本紀錄之間新增可由 `ANALYSIS_FLOW.md` 編輯的橫向 NGS 分析流程圖，呈現輸入/QC、五條二級分析、三級分析、Prioritization、Phenotype/clinical context、判讀與三種報告輸出。
- 上傳個案清單支援新版院內 xlsx：標題列與「檢體編號」可位於任一列／欄，可用「檢驗項目」取代「檢驗名稱」，並會合併解析所有工作表。
- 上傳個案清單可一次選擇多個 xlsx，完成後只顯示本次各檔案結果；歷次上傳紀錄預設收合，按小三角形才載入並展開。
- 二級分析產生的 WES/WGS DGX2 指令預設啟用 Manta、ExpansionHunter 與 AutoMap，直接保留三個模組的輸出；WES 仍同時啟用 gCNV。
- 修正二級分析 DGM→DGX2 FASTQ 路徑映射：`/home/datalake_Raw` 現在正確對應 DGX2 的 `/datalake_Raw/datalake_Raw`，避免 FASTP container 因 bind source 不存在而失敗。

### v6.8 — 2026-07-20

- 二級分析的 DGX2 清理指令「複製」改為一鍵寫入剪貼簿；plain-HTTP intranet 會自動使用隱藏文字框 fallback，不再顯示要求手動複製的視窗。
- PGx 分析區的基因名稱維持黑字、diplotype/allele 改用淡紫底，只有需注意 phenotype 顯示紅字。逐 gene 用藥建議簡化為藥物與 CPIC/FDA 建議兩欄，需調整項目優先；概覽把 generic fallback 改為「其他用藥建議（參考完整建議與最新藥品仿單）」，並確保相關 JSON 藥物與 MT-RNR1 TSV 補充藥物各自只進一個主要分類。star allele 說明改為欄名 hover tooltip。

### v6.7 — 2026-07-18

- Test type 新增 `TITAN-WGS`類別。主畫面個案搜尋與個案清單的三種類別旁新增 `only` 快速單選。
- TITAN-WGS 載入個案時預設隱藏診斷分析；基本資料下方可顯示／隱藏診斷分析區域。

### v6.6 — 2026-07-17

- SNV/Indel 卡片的 Score、ClinVar、ACMG 與所有 in-silico tools 新增 `ⓘ` 註解；可查看計分公式、ClinVar 版本、ACMG 的 manual/GeneBe/in-house 來源比較，以及各工具的 PP3/BP4 threshold、文獻連結與 PMID。

- AlphaMissense、ESM1b、VARITY_R、BayesDel、REVEL、SpliceAI、PhyloP、GERP、SIFT 改依正式 calibration/ClinGen 建議套色；其他尚無通用 PP3/BP4 calibration 的工具會清楚標示為模型 cutoff 或 contextual evidence。
- 三級分析的 Extra VEP 現在會在 dbNSFP 實際提供 `REVEL_score` 時補入 REVEL，並與 MetaRNN、SpliceAI 一同存入。
- SNV/Indel tier 簡化為 `1A / 1B / 1C / 2`；`1C — Predicted suspect` 除 ACMG points ≥4 外，會依 Core（P-KNN LLR、AlphaMissense、BayesDel、Pangolin）與 Extra-VEP（REVEL、SpliceAI）門檻納入 reviewer 候選。原 ClinVar P/LP 0★/conflicting tier 移除，其餘歸入 `2 — Other`。

### v6.5 — 2026-07-16

- 更改三級分析輸出（包涵 Nextflow 及 Post-processing 輸出）至 `/home/datalake_Intermediate/pipeline/tertiary_output/{sample}/`

### v6.4 — 2026-07-09

- 健檢報告 PGx 回報範圍更新為 21 個 CPIC Level A genes，報告改為重點摘要、基因型與表現型、用藥建議分類、完整用藥建議四段式。
- DRAGEN 三級分析新增 ploidy VCF 資訊；基本資料性別欄會核對 M/XY、F/XX，一致顯示綠底，不一致或非典型性染色體 call 顯示紅底及 ploidy VCF 結果。

### v6.2 — 2026-07-08

- SNV/Indel 疾病清單可合併 GenCC / ClinGen / MONDO gene-disease，在 OMIM 未列 disease 或有額外非重複疾病時補充顯示。
- SNV/Indel 新增「本院 AF（AF_nckuh）」：顯示該變異在院內 WGS cohort 的等位基因頻率。

- WGS 固定套組新增血液科 Lymphoid Neoplasm Panel 與 Myeloid Neoplasm Panel。

### v6.1 — 2026-06-29

- Secondary findings 新增中風相關基因區塊。
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
