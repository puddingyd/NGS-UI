# INTERVENE PRS → 台灣多疾病研究流程

一般使用者請先看[PRS流程說明：從基因資料到疾病風險](../../docs/research/prs/PRS_WORKFLOW.md)（v10.32，2026-11-13），由概念、例子與選擇理由解釋整個流程；各步驟可展開生資方法、參數與程式對照，文末另列完整路徑與入口。平台PRS標題「研究用途」旁的「說明」可下載同一份文件。以下舊版本細節保留追溯，現行六疾病流程以該說明及操作文件為準。

後續方法審查請看[PRS審查回應與台灣在地化評估](../../docs/research/prs/PRS_REVIEW_RESPONSE_2026-11-14.md)（2026-11-14）：逐項回應外部 review，評估約1,200例成人健檢／社區WGS參考分布、台灣人口基準與BBJ HR。包含本次數值與來源核對、可做／待資料／不宜直接採用的建議；尚未變更現行模型或平台行為。

v9.91 新增未來20年風險（年齡≤60），與10年及至80歲使用相同參數／競爭風險算法。UI新病史表單預設否，已保存病史保留。舊評估按風險計算產生新版快照即可，不重跑基因計分；nf_work維持原路徑並透過既有bind mount使用SSD。

**六疾病基礎版 v9.90（2026-11-05）**：CHD＋T2D＋All cancers＋Breast cancer＋AF＋Epilepsy 共六疾病已接上共用計分／風險／UI／報告。準備與驗收請先看 [六疾病操作文件](../../docs/ops/PRS_MULTI_DISEASE.md)；來源與 1,232 筆重現見 [參數發行說明](releases/tw_multi_v1/README.md)。`multidisease/` 為共用計分程式，`absolute_risk/registry.py` 為模型對應，`releases/tw_multi_v1/` 為固定研究參數。大型參考模型須在院內準備；完整 WGS 尚未在本次驗收。以下 CHD V1 文件保留供舊版本／凍結計分流程追溯，新增疾病以六疾病文件為準。

2026-11-02 v9.72：接上三級自動計分與兩階段 PRS 表格。先由原始 gVCF 產生 score／percentile／祖源，再由使用者輸入臨床資料並計算風險。採 INTERVENE 原研究 EstBB HR＋台灣 GBD2019，輸出未來10年、目前至80歲風險與 PRS 40–60% 比較曲線。原文224筆 Estonia CHD 結果驗算通過；**台灣 outcome validation 尚未完成，結果僅供研究**。TWB／BBJ 未納入本模型。

設定計分環境、安裝參數包、操作 UI 與快照規格見 [操作文件](../../docs/ops/PRS_CHD_RESEARCH.md)。來源見 [DATA_SOURCES.md](DATA_SOURCES.md)，驗證範圍見 [VALIDATION.md](VALIDATION.md)。

## 目錄與既有版本

- `scoring_v211/`：從使用者 ZIP 選取的既有計分程式、測試與歷史文件，逐檔內容保持不變。
- `IMPORT_PROVENANCE.json`：來源 ZIP SHA256、實際匯入清單及各檔案 hash。它取代原 ZIP 的完整發行清單；未匯入正式 `config.json`、操作 manifest 或舊備份。
- `absolute_risk/`：Python 標準函式庫實作的研究風險引擎及 CLI。
- `examples/`：**全為合成參數與合成個案**，僅供離線驗證。沒有任何台灣發生率或發表 HR。
- `upstream/`：原作者 sex-stratified R 腳本、MIT 授權與固定 commit provenance，供公式核對；原檔未執行。實際選定模型與原資料來源見操作文件。
- `tests/`：解析解、條件風險、七分組、輸入拒絕、CLI 與原始檔案完整性測試。

原 v2.1.1 已有三案 fresh、兩案並行等伺服器驗證紀錄，見 `scoring_v211/VALIDATION_REPORT_v211_2026-11-02.md`。這些是歷史證據，不代表本次重跑 WGS。大型 `canonical_v0` 模型、gVCF、Nextflow work、實測 outputs 均未納入 repo。

歷史 SOP 內的伺服器路徑、舊套件名稱及安裝指令保留作追溯，**不要整段照抄到 repo checkout 執行**。現行 scoring 入口是 `scoring_v211/pipeline_v21.py`。在院內另設 runtime config，把 `raw_script`、`method2_lock` 指到實際安裝位置。**v9.80平台入口會覆寫 `project`／`raw_work`，個案產物預設放既有三級暫存根目錄的 `nf_work/prs/{個案ID}/`（跟隨三級work root設定），完整計分／QC另發布至08；原專案只供唯讀資源。**直接手動執行科學 CLI 仍遵循它收到的 config，須自行確認輸出位置。model lock 內的絕對 model/source 路徑有嚴格核對，不能為換目錄任意重新 freeze。未匯入實際 manifest；新 manifest 用 `sample, mode, gvcf, raw_source` 四欄 TAB 格式，fresh 的 raw_source 必須空白。測試用 golden IDs 保留在原始程式，不能當作新個案設定。

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

- `intervene_discrete_v1`：每五年增量 `S_start * (1-exp(-5*hazard))`，存活比例乘 `exp(-5*(hazard+mortality))`。只接受完整五年區間。用於重現作者的離散近似；已以 EstBB HR＋Estonia baseline 重現 Data 4 的224筆 CHD 結果（含男女各七組）。
- `competing_exact_v1`：區間長度 dt，增量 `S_start * hazard/(hazard+mortality) * (1-exp(-dt*(hazard+mortality)))`。零總 hazard 時增量為 0；可切分非五歲整數的評估年齡。這是不同數值方法，不能宣稱與前者逐值相同。

未來風險從「目前仍存活且無該疾病」開始，初始 S=1，往指定年齡前進。不是出生累積率直接相減。10／20 年期間超過 80 歲時回傳 unavailable，不截短冒稱完整期間；仍可回傳目前至 80 歲風險。出生至80歲為另一種研究曲線，僅適用含完整年齡範圍的參數；本版台灣包從20歲起，不提供出生累積率。

## 正式資料包規格

JSON schema 名稱為 `CHD_RESEARCH_BUNDLE_V1`，欄位形狀見 synthetic_bundle，但**不允許把 synthetic 改名當正式資料**。正式研究包狀態只能是 `unvalidated_external`，不得宣告 clinical validated。

必填 provenance：model_id、endpoint=`first_CHD`、population、baseline_source、baseline_year、baseline_release、hr_source、limitations、method、units、score_definition。來源欄應記 DOI/URL、資料檔 checksum、匯出 query、日期、HR 的祖源、年齡／性別模型與 phenotype mapping 審查。score_definition 固定 A 或 B、model_lock_sha256、reference_id、ancestry=EAS。Model A 作主分析、B 作敏感度；兩者必須使用獨立且經審查的風險包，不能平均或依個案結果選模型。

`bins` 每列包含 sex=`male/female`、age_start（含）／age_end（不含）、incidence、prevalence、all_cause_mortality、chd_mortality 及七組 hr。兩個性別都須連續覆蓋 [age_min,80)，age_min 預設0，本版成人包為20；禁止重疊或缺值。HR 必須正值，40-60 組 HR=1。率為**每人每年**，盛行率為 [0,1) 比例；不得直接匯入每十萬人數值、age-standardized rate 或把盛行率當年度發生率。原資料年齡起點、尤其 1–4 歲與五年 bins 的差異須先審查，不能按 CSV 列順序硬接。

輸入 summary 為既有 Method2 summary.json；必須有技術 PASS、相符的 model lock、樣本身分、正確 A/B 位點數與 percentile group。原預設 raw-only summary 的 ancestry=NOT_ASSESSED 會拒絕風險計算。既有 `--with-legacy-ancestry` 可帶入原 ancestry QC；該方法使用 PGS 位點交集，不是已驗證的獨立全基因組 ancestry panel。新 ancestry sidecar 的可追溯整合仍待實作，**勿手動把 NOT_ASSESSED 改成 EAS**。

person JSON 必填 sample、age（評估日期的歲數，可小數）、sex、assessment_date（ISO 日期）、prior_chd=false。病史未知、已有 CHD、祖源不明或超出年齡範圍均拒絕；不能從姓名推測性別或用 0 代替缺值。

## 現況與後續驗收

1. **已完成**：原研究 EstBB 男性 age/sex、女性 sex 參數核對；原文224筆 CHD 曲線重現；台灣成人參數包、API、UI、獨立快照、來源變更失效與缺件拒絕。程式隨 Git 更新，參數包須另安裝；沒有自動部署。
2. **院內操作驗收**：設定計分環境與資料包，從 UI 啟動一個真實 WGS 的 raw gVCF → score／percentile／祖源流程，再輸入臨床資訊計算風險，確認樣本對應及 UI 數值。raw-only 的 NOT_ASSESSED 仍不可估算，不可手改 EAS。
3. **研究限制**：GBD IHD 與 CHD 定義不完全一致；歐洲 HR、Model A 子集與 EAS735 ECDF 的台灣相容性未 outcome validated。不確定性區間尚未實作，明確輸出 null，不作臨床決策。
4. **已接上排程／worker**：三級發布後自動讀原 gVCF，完成score即發布，再補祖源；UI有獨立重試／取消與風險計算。院內fresh WGS、5–10 WGS負載與真實partial-task recovery仍待使用者部署後驗收，其他疾病之後再擴充。

## 來源

- [PGS004879](https://www.pgscatalog.org/score/PGS004879/)：原始 610,677 loci；A/B 為另行固定子集。
- [Jermy et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11169548/)：國家／年齡／性別基準與競爭死亡風險、歐洲祖源限制。
- [作者固定 commit](https://github.com/intervene-EU-H2020/flagship/tree/179986db2f9af3b6eb32e4e16dc9201a23d2a29c)：`upstream/PROVENANCE.json` 記錄保存範圍。
- [GBD Results](https://vizhub.healthdata.org/gbd-results/)：後續正式資料匯出入口。
