# 台灣 CHD absolute risk：資料尋找與核對結果

核對日期：2026-11-02。**已取得台灣 GBD2019 基準、原論文驗算表與東亞候選關聯資料；尚未形成可輸出個案風險的正式參數包。** 本次沒有更動計分或風險引擎，也沒有部署。

## 取得狀態

| 所需資料 | 本次結果 | 用途與限制 |
|---|---|---|
| 台灣 IHD 發生、盛行、IHD 死亡、全死因死亡 | 已下載作者公開的三份 GBD2019 CSV，四類資料齊全 | 男女、年齡、2019 年與上下界；第一年齡組是 1–4 歲，沒有未滿 1 歲 |
| 原研究 CHD 預期輸出 | Supplementary Data 4：896 筆；Data 7：896 筆 | 四個地區 × 男女 × 七分組 × 16 年齡終點。Data 7 另含 meta-analysis HR 路徑結果 |
| CHD 模型選擇 | 已核對 Supplementary Table 3（PDF 第52頁）及 Data 4 | 男性 Age and Sex；女性 Sex。女性的 HR 不應強行套男性年齡斜率 |
| ICD／疾病定義 | Data 8、9 已取得 | CHD 與 GBD IHD 不能當完全相同；Data 9 在 UKB 的定義比對列為 82% |
| 歐洲原研究七分組 HR、SE、組內年齡中位數 | 尚未取得完整原始參數 | Data 2 為 per-SD HR；固定 Git tree 沒有腳本引用的原始 HR CSV |
| 台灣 TWB 關聯資料 | 作者 Drive 的 NBCTaiwan 目錄可讀；已下載男女分組表與 cohort 摘要、年齡四分位 | 檔名為 HR，但欄位是 **OR**，且區分 `prs`、`prs_norm`；不能直接當時間風險的 HR |
| 日本 BBJ 分組參數候選 | 已下載女性分組表與年齡性別分組表，均有 I9_CHD | 欄名 V1…，需核對欄位、score 版本、模型與 endpoint；不能據此宣稱台灣校正完成 |

## 1. 台灣基準：現在不需要先自行匯出 GBD2019

來源為[作者公開資料目錄](https://github.com/intervene-EU-H2020/flagship/tree/179986db2f9af3b6eb32e4e16dc9201a23d2a29c/AbsoluteRiskEstimation)，固定 commit `179986db2f9af3b6eb32e4e16dc9201a23d2a29c`：

- `Sex_Stratified_Incidence.csv`
- `Sex_Stratified_Prevalence.csv`
- `Sex_Stratified_Mortality.csv`（同時含 `All causes` 與 `Ischemic heart disease`）

保留來源地名 `Taiwan (Province of China)` 供精確篩選，畫面將來可顯示「台灣」。選 `metric=Rate`、`year=2019`、兩性與 IHD／All causes 後共160列；去掉 `All ages` 為152列，1–79歲範圍為128列。每個性別年齡有完整四類資料，鍵值無重複、數值有限且 `0 ≤ lower ≤ val ≤ upper`，全死因死亡點估計皆不小於 IHD 死亡點估計。

來源 Rate 為每十萬人；轉引擎時需除100000。發生／死亡為年度率，盛行則是比例，不能當年度 hazard。上下界有保留，但沒有 posterior draws 或跨年齡／性別／疾病的 covariance；不能聲稱已取得完整聯合不確定性分布。

原資料只有 1–4、5–9…75–79 歲，作者腳本卻將第一格乘5並標至5歲。原方法重現應明記這項約定；正式研究包不得默默把 1–4 改成 0–4。現有引擎要求 [0,80) 完整覆蓋，所以這份來源仍不能直接匯入。可後續取得 <1 資料，或另設有明確最低年齡的成人模型。2019 是版本對齊用的基準年，並非目前最新台灣風險；更新年度需整組替換四類資料及重新驗算。

## 2. 原論文可用的驗算證據

[Jermy et al., Nature Communications 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11169548/) 的12份 Supplementary Data、Supplementary Information PDF 與附件說明 PDF 已下載並解析。

- Data 4（`MOESM7`）：CHD 的 Finland、Estonia、United Kingdom、Massachusetts (USA)，5、10…80歲累積率及區間。單位為百分比；公開值有四捨五入，後續驗收容差需反映顯示精度。
- Data 7（`MOESM10`）：相同 CHD 維度，比較 biobank-specific 與 meta-analyzed HR 的累積率。輸出表不是 HR 參數表，不可反推後再當獨立驗證。
- Data 8（`MOESM11`）第9列：CHD ICD10 `I200|I21|I22|I23|I24|I25|I46`；ICD9 `410|4110|411|412|413|414|798`。保留作者格式，正式 ICD matching 要區分小數點正規化與前綴規則。
- Data 9（`MOESM12`）CHD 區塊的 J13 為0.82（Excel顯示82%）：這是原文在 UKB 中的定義重疊指標，不是台灣驗證率，也不是模型準確率。
- Data 11（`MOESM14`）：提供發病年齡四分位**切點**，並非男性每個 percentile × 年齡組的 `MedianAAO`。不能拿三個切點替代四組回歸所需的中位數。

Data 2（`MOESM5`）是 per-SD HR。另有需向作者核對的來源內部不一致：CHD Both / All ages 的 D7 Beta=0.28，G7 HR=1.66，與 `exp(0.28)` 不符；本次未自行更正或用此列生成參數。

## 3. 新找到的東亞資料：先審查再選模型

來源連結是作者 `AbsoluteRiskEstimation/scripts/multi_ancestry.R` 引用的[研究資料夾](https://drive.google.com/drive/folders/1W3tqDX6K51TfKZhHabi4Gn5xVcjLQZsB)。透過連接的 Google Drive 可讀取；這不等於確認匿名公開存取或再散布授權。檔案在2024原論文之後仍有更新，不能默認是該篇論文的凍結版本。

### 台灣 TWB

[NBCTaiwan 目錄](https://drive.google.com/drive/folders/1cOuc37civSfu1CY8fKaKSGPJY5Vi5_Ti) 已取得：

- `HR_FemaleSample_TWB_EAS.csv`
- `HR_MaleSample_TWB_EAS.csv`
- `TWB_EAS_baseline_summary_stats.csv`
- `TWB_EAS_age_quartiles.csv`

男女分組檔均有 I9_CHD，各32列，包含兩種 `PRS_TYPE` 與不同尾部分組。欄位有 `Beta, SE, OR, CIpos, CIneg`，case count 以 lower/upper 區間表示。尚不能從檔名確認 Cox 模型；不得把 OR 改名為 HR，亦不得把 count 區間中點當精確個案數。cohort 的盛行率／追蹤摘要不是一般台灣人口的年度發生率，不能代替 GBD baseline。

這組資料有價值，但仍需作者確認：OR／模型類型、prevalent 或 incident outcome、follow-up 定義、共變項、`prs_norm` 方法、PGS004879 對應與參考組、資料使用／引用方式，以及是否能提供 incident CHD Cox HR。只找到台灣關聯資料並不代表已有台灣 absolute-risk calibration。

### 日本 BBJ

[BiobankJapan 目錄](https://drive.google.com/drive/folders/13fWpt3iY3Ph6G15g9RqlYYo4Wmlgc0uS) 已取得 `HR_FemaleSample_BBJ.csv` 與 `HR_AgeandSexStratified_BBJ.csv`，I9_CHD 分別有16與128列；有 `<20%`、`>95%` 等合併尾組，不必用小尾組倍率平均假造。年齡性別檔有四年齡組與組內年齡欄。

作者程式給出原 positional schema，可作核對線索，但不能未審查便依位置直接套用。需要確認 V欄定義、HR 對應、PGS004879 與 A/B 子集相容性、cohort 選樣與病史定義、分組邊界及 score normalization。日本→台灣仍是跨族群／cohort 外推；東亞來源可評估為獨立敏感度路徑，不自動替代2024歐洲方法重現。

## 4. 最小剩餘資料需求與下一步

**原方法重現路徑**：需一個原研究 cohort 的完整 CHD 參數，或原研究 meta-analysis 的七分組參數。男性需六個非參考組 × 四年齡組的 log-HR、SE、年齡界線及 MedianAAO；女性需六組 log-HR／SE，40–60%作參考。還要 score 版本、covariances／bootstrap 規則（若可提供）與未四捨五入的輸出。作者程式引用 `R10_HR_AgeandSex_Stratified_FinnGen.csv`、`R10_HR_FemaleSample_FinnGen.csv`、`MetaAnalysisPercentilesAgeandSexStratified.csv` 等，但本次未在固定 Git tree／發表附件找到完整數據。沒有這些檔案時不能宣稱逐值重現原論文。

**台灣研究路徑**：先核對 TWB 檔案的 OR 與 score 定義；如無可用 Cox HR，再評估 BBJ 或歐洲 HR 外推，獨立版本並明列來源。接著處理年齡起點／疾病定義差異、建立 GBD adapter 與驗算報告，最後才產生 `unvalidated_external` 參數包。

本階段不需要使用者提供院內長期追蹤資料才能繼續整理與重現研究。要進一步宣稱台灣個人風險準確或做臨床用途，仍需要符合目標族群的外部／院內 outcome validation；資料取得與模型驗證是不同里程碑。

## 保存位置與重取

- 來源索引：`SOURCE_MANIFEST_20261102.json`（31個原檔及4個擷取／稽核檔的路徑、網址、SHA256）。
- 原始資料：`runtime/sources_20261102/{gbd2019,publication,author_code,candidate_associations}/`。
- 擷取：`runtime/sources_20261102/extracted/`，含台灣四類 rate、兩份 CHD 驗算表與 `SOURCE_AUDIT.json`；保留來源 row，未偽裝成風險包。
- `runtime/` 已 gitignore：**只有文件與 manifest 隨 git 同步，原檔不會隨 git pull 出現在伺服器**。重取時使用 manifest 中原始 URL，驗證 SHA256；Drive 檔案須具備讀取權限。原檔與其使用／再散布條款分開審查，不沿用程式 MIT 授權推定數據授權。

本次檢查為資料完整性、欄位與來源稽核；沒有重跑 WGS，沒有把 OR 或 per-SD HR 注入引擎，沒有輸出病人 absolute risk。
