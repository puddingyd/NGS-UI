# INTERVENE PRS → 台灣 CHD 研究流程

2026-11-02：第一階段程式納管與風險計算基礎實作。**目前沒有正式台灣風險參數包；沒有完成原研究數值重現或台灣 outcome validation。不可把合成資料輸出當作個案風險。**

## 目錄與既有版本

- `scoring_v211/`：從使用者 ZIP 選取的既有計分程式、測試與歷史文件，逐檔內容保持不變。
- `IMPORT_PROVENANCE.json`：來源 ZIP SHA256、實際匯入清單及各檔案 hash。它取代原 ZIP 的完整發行清單；未匯入正式 `config.json`、操作 manifest 或舊備份。
- `absolute_risk/`：Python 標準函式庫實作的研究風險引擎及 CLI。
- `examples/`：**全為合成參數與合成個案**，僅供離線驗證。沒有任何台灣發生率或發表 HR。
- `upstream/`：原作者 sex-stratified R 腳本、MIT 授權與固定 commit provenance，供公式核對；不是已選定的 CHD age/sex 模型。未執行該腳本。
- `tests/`：解析解、條件風險、七分組、輸入拒絕、CLI 與原始檔案完整性測試。

原 v2.1.1 已有三案 fresh、兩案並行等伺服器驗證紀錄，見 `scoring_v211/VALIDATION_REPORT_v211_2026-11-02.md`。這些是歷史證據，不代表本次重跑 WGS。大型 `canonical_v0` 模型、gVCF、Nextflow work、實測 outputs 均未納入 repo。

歷史 SOP 內的伺服器路徑、舊套件名稱及安裝指令保留作追溯，**不要整段照抄到 repo checkout 執行**。現行 scoring 入口是 `scoring_v211/pipeline_v21.py`。在院內另設 runtime config，把 `raw_script`、`method2_lock` 指到實際安裝位置；`raw_work` 仍為 PRS project root。model lock 內的絕對 model/source 路徑有嚴格核對，不能為換目錄任意重新 freeze。未匯入實際 manifest；新 manifest 用 `sample, mode, gvcf, raw_source` 四欄 TAB 格式，fresh 的 raw_source 必須空白。測試用 golden IDs 保留在原始程式，不能當作新個案設定。

## 離線驗證

在 repo 根目錄執行（Python 3.10+，沒有新增套件依賴）：

```bash
python3 -m unittest discover -s research/prs/tests -v
# macOS 的 /var → /private/var symlink 會影響原版路徑字串測試；使用實體 tmp 路徑。
TMPDIR=/private/tmp python3 -m unittest discover -s research/prs/scoring_v211/tests -q
python3 -m research.prs.absolute_risk validate --bundle research/prs/examples/synthetic_bundle.json
python3 -m research.prs.absolute_risk evaluate \
  --bundle research/prs/examples/synthetic_bundle.json \
  --summary research/prs/examples/synthetic_summary.json \
  --person research/prs/examples/synthetic_person.json \
  --allow-synthetic --output /tmp/prs_synthetic_result.json
```

輸出拒絕覆寫既有檔案；重新測試請用另一個輸出檔名。`validate` 只驗證結構與數值範圍，回傳 `scientific_validation=NOT_VERIFIED`。所有結果都包含 `research_only=true`、`clinical_validation=NOT_VALIDATED`、參數狀態與三份輸入的 SHA256。沒有實作信賴區間時明確輸出 `uncertainty_interval=null`，不產生假精確區間。

## 兩個明確分開的算法

七個分組採 `<20, 20-40, 40-60, 60-80, 80-90, 90-95, >=95`，界線值歸右側；不能先四捨五入 percentile。各組比例固定為 `.2,.2,.2,.2,.1,.05,.05`，這是原方法分位組的假設，不代表台灣各年齡疾病未發生者的實測組成。

每年齡／性別區間先計 `I0 = I_population / sum(proportion[g] * HR[g])`，`I_g = I0 * HR[g]`。用人口盛行率估 `hazard_g = I_g / (1-prevalence)`；其他原因死亡率為全死因減 CHD 死亡率。此處沿用各 PRS 組相同盛行率及其他死因死亡率的模型假設。

- `intervene_discrete_v1`：每五年增量 `S_start * (1-exp(-5*hazard))`，存活比例乘 `exp(-5*(hazard+mortality))`。只接受完整五年區間。用於重現作者的離散近似；僅公式經合成資料測試，**未與實際 CHD tables 對照**。
- `competing_exact_v1`：區間長度 dt，增量 `S_start * hazard/(hazard+mortality) * (1-exp(-dt*(hazard+mortality)))`。零總 hazard 時增量為 0；可切分非五歲整數的評估年齡。這是不同數值方法，不能宣稱與前者逐值相同。

未來風險從「目前仍存活且無 CHD」開始，初始 S=1，往指定年齡前進。不是出生累積率直接相減。10 年期間超過 80 歲時回傳 unavailable，不截短冒稱 10 年；仍可回傳目前至 80 歲風險。出生至 80 歲為另一種研究曲線，可使用引擎 start_age=0。

## 正式資料包規格

JSON schema 名稱為 `CHD_RESEARCH_BUNDLE_V1`，欄位形狀見 synthetic_bundle，但**不允許把 synthetic 改名當正式資料**。正式研究包狀態只能是 `unvalidated_external`，不得宣告 clinical validated。

必填 provenance：model_id、endpoint=`first_CHD`、population、baseline_source、baseline_year、baseline_release、hr_source、limitations、method、units、score_definition。來源欄應記 DOI/URL、資料檔 checksum、匯出 query、日期、HR 的祖源、年齡／性別模型與 phenotype mapping 審查。score_definition 固定 A 或 B、model_lock_sha256、reference_id、ancestry=EAS。Model A 作主分析、B 作敏感度；兩者必須使用獨立且經審查的風險包，不能平均或依個案結果選模型。

`bins` 每列包含 sex=`male/female`、age_start（含）／age_end（不含）、incidence、prevalence、all_cause_mortality、chd_mortality 及七組 hr。兩個性別都須連續覆蓋 [0,80)，禁止重疊或缺值。HR 必須正值，40-60 組 HR=1。率為**每人每年**，盛行率為 [0,1) 比例；不得直接匯入每十萬人數值、age-standardized rate 或把盛行率當年度發生率。原資料年齡起點、尤其 1–4 歲與五年 bins 的差異須先審查，不能按 CSV 列順序硬接。

輸入 summary 為既有 Method2 summary.json；必須有技術 PASS、相符的 model lock、樣本身分、正確 A/B 位點數與 percentile group。原預設 raw-only summary 的 ancestry=NOT_ASSESSED 會拒絕風險計算。既有 `--with-legacy-ancestry` 可帶入原 ancestry QC；該方法使用 PGS 位點交集，不是已驗證的獨立全基因組 ancestry panel。新 ancestry sidecar 的可追溯整合仍待實作，**勿手動把 NOT_ASSESSED 改成 EAS**。

person JSON 必填 sample、age（評估日期的歲數，可小數）、sex、assessment_date（ISO 日期）、prior_chd=false。病史未知、已有 CHD、祖源不明或超出年齡範圍均拒絕；不能從姓名推測性別或用 0 代替缺值。

## 缺件與下一個驗收點

1. **原方法重現資料**：原研究 CHD 分組 HR、SE/共變異、年齡切點、年齡插值／外推規則及選定模型；至少一個原研究國家的 baseline 與 published expected curves。作者 repo 腳本仍引用作者本機的 HR CSV；本次沒有取得完整參數，不能用 Catalog per-SD HR 代替七組年齡別 HR。補齊後逐性別、年齡、分組對照預期數值，保存誤差與來源。
2. **台灣基準**：同一 GBD release／年度的 IHD 年齡性別別發生、盛行、全死因與 IHD 死亡率及 uncertainty。重現先用 GBD2019；台灣研究版另固定年度，不能混用版本。確認 I9_CHD / IHD ICD 定義一致性。GBD 下載可能需使用者帳號，這次未下載或同意任何條款。
3. **外推與計分相容性**：2024 原研究 HR 為歐洲祖源；EAS735 的 ECDF 與台灣人口並不等價。查核 A/B 子集、全球3330頻率補值、來源PGS與published HR相容性；評估 A/B、高權重缺失、palindromic、EAS reference bootstrap 及跨組穩定性。現有 2pp 是 review flag，不是 CI。
4. **不確定性**：取得完整來源後再加入參數抽樣及固定 seed；區分 HR／baseline 區間、參考抽樣與 A/B sensitivity。歐洲→台灣的不可辨識外推偏差不能靠 narrow CI 消除。
5. **作業驗收**：mid-workflow partial-task recovery、5–10 WGS 資源測試需在院內執行；本機測試沒有代替它們。
6. **平台整合（尚未實作）**：獨立 PRS router/service/worker、全站資源限制、sample lock、取消與原子發布、歷史 snapshot；gVCF score 與年齡風險快照分開，換年齡不重跑 genotype。資料不要放進既有 SNV review TSV、ACMG tier 或單基因診斷報告。結果 schema 可作介面起點，但 API / UI 尚未啟用。資料庫、實際個案及輸出留院內 runtime，不進 Git。

## 來源

- [PGS004879](https://www.pgscatalog.org/score/PGS004879/)：原始 610,677 loci；A/B 為另行固定子集。
- [Jermy et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11169548/)：國家／年齡／性別基準與競爭死亡風險、歐洲祖源限制。
- [作者固定 commit](https://github.com/intervene-EU-H2020/flagship/tree/179986db2f9af3b6eb32e4e16dc9201a23d2a29c)：`upstream/PROVENANCE.json` 記錄保存範圍。
- [GBD Results](https://vizhub.healthdata.org/gbd-results/)：後續正式資料匯出入口。
