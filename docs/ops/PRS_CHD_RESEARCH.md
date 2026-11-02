# CHD PRS 研究風險 UI — v9.71

## 本版範圍

讀取既有 PGS004879 Method2 Model A 結果，輸出：

1. 評估日起未來10年 CHD 風險。
2. 評估時年齡至80歲 CHD 累積風險。
3. 兩條比較曲線：此個案與相同年齡／性別的 PRS 40–60% 參考組。

僅20至79歲、確認尚未有 CHD、具備合格 EAS 祖源證據及相符 model lock 的個案可估算。年齡大於70歲時，10年結果為 unavailable，仍可輸出至80歲；不截短期間。所有風險都以「評估時仍存活且無 CHD」為條件，不是出生累積率直接相減。既有 CHD 或病史未知不輸出。

本版沒有 gVCF 排程／重跑按鈕，不從 SNV TSV 計分，也不將研究結果加入診斷／健檢 DOCX 或列印報告。計分仍用既有 v2.1.1 pipeline，UI 接收完成的 summary。換年齡只重算風險，不重跑 genotype。TWB、BBJ 未納入本模型。

## 固定參數與科學驗算

- HR：INTERVENE 原研究 Estonian Biobank（EstBB）歐洲祖源，非跨 cohort meta-analysis。男性使用六個非參考 PRS 組各四個年齡組的 HR/SE/MedianAAO，按作者 `1/SE²` 加權線性回歸擬合 log-HR；外側年齡固定於最近組的**擬合值**。女性使用六組 sex-specific HR；參考組 HR=1。
- 基準：作者固定 commit 的台灣 GBD2019 IHD incidence/prevalence/all-cause mortality/IHD mortality。
- 原方法驗算：以 Estonia baseline 與作者五年離散算法，逐一比對 Supplementary Data 4 的 **224筆 CHD 結果，全部通過**。最大絕對誤差 `0.004995100317904511` 百分點，小於表格保留兩位小數所容許的 `0.0050001` 百分點。沒有重現 bootstrap CI，也沒有聲稱原 R 執行環境完全重建。
- UI 算法：`competing_exact_v1`，同一組 HR／基準率，採每個五年區間內固定疾病／死亡 hazard 的連續競爭風險積分，支援50.5歲等非五年邊界。**它不是與論文五年離散近似逐值相同的數值方法**。台灣兩性七組、20–75歲每5歲起點的10年／至80歲測試網格，兩方法最大差異為 `1.4436119070748987` 百分點；不同方法必須保留名稱，不得混稱。
- 成人參數包只含 [20,80)，未將來源的1–4歲改寫成0–4歲。驗算出生曲線時，才明記沿用作者第一格乘5年的約定。
- Model A、EAS735 參考分布及原始 model lock 固定；分位組跨族群與 A 子集的 outcome 相容性仍未驗證。原研究歐洲 HR 套用台灣率是研究外推，並非完成台灣臨床校正。

模型 ID：`intervene_estbb_tw2019_A_exact_v1`。`clinical_validation=NOT_VALIDATED`；`uncertainty_interval=null`。不生成吸菸／血壓／血脂或治療情境，也不套臨床高低風險門檻。

## 安裝資料（由操作者部署，不由 Codex SSH）

程式與文件隨 Git 更新；參數資料不放 Git。此次本機已產生：

```
research/prs/runtime/tw_chd_v1/bundle.json
research/prs/runtime/tw_chd_v1/reproduction.json
research/prs/runtime/tw_chd_v1.zip
```

將資料包解壓到伺服器 `$NGS_UI_HOME/biotools/prs/tw_chd_v1/`，或以 `NGS_UI_PRS_BUNDLE` 指定 `bundle.json` 的完整路徑。正式 production checkout 的 NGS_UI_HOME 通常是 `~/NGS_UI`。包內 `reproduction.json` 保留逐筆 expected/calculated/error 與來源 hash。

本版 bundle SHA256：

```
322d4e68d49a2047d7326172a609b3898b8754f658f800f50806f4e2276a0a27
```

後端固定檢查這個 digest 與模型版本；synthetic、被改寫或僅自行標記 PASS 的參數不能啟用。原檔 bytes 必須保留，不可重新排版 JSON。後續更新模型需另做驗算、版本與 release digest 更新。

讀取資料包檢查（在 repo 根目錄）：

```bash
PYTHONPATH=backend python3 -c 'from app.services.prs import load_bundle; b,h=load_bundle(); print(b["model_id"],h)'
```

## 匯入一個已完成的計分結果

```bash
PYTHONPATH=backend python3 scripts/import_prs_summary.py \
  --sample-id '實際已登錄個案ID' \
  --summary '/伺服器上既有PRS結果/summary.json'
```

`summary.sample` 必須符合 UI 個案 ID，或該個案 `pipeline_source.json` 內的 `source_sample_id`；不接受任意指定 alias。僅接受 `PGS004879_METHOD2_V0_20261031`、`METHOD2_TECHNICAL_PASS`。匯入保留原始 summary、原始 SHA256、來源識別與匯入時間，連結到當時 pipeline source marker 的 hash。

「成功匯入」不等於「適合風險計算」。UI 還會檢查 Model A、locus counts、percentile/group 一致性、model lock 與合格 EAS 祖源證據。原 pipeline 預設 raw-only 的 `NOT_ASSESSED` 無法估算；應取得既有流程的真實 ancestry QC，不得手改為 EAS 或偽造 PASS。

使用者自行更新／重啟平台後，載入個案，在「基本資料」下方開啟「CHD 多基因風險」，填入評估日期、當時年齡、模型使用性別與 CHD 病史後按「估算風險」。未安裝參數／未匯入 score 時，UI 明示缺件，不顯示合成數值。

## API 與保存

- `GET /api/samples/{sample_id}/prs`：登入保護；回傳 available、issues、score 摘要與有效的上次評估。
- `POST /api/samples/{sample_id}/prs/assessments`：只接受 `age, sex, assessment_date, prior_chd`。percentile、HR、參數路徑皆由伺服器讀取，不能由 request 指定。
- 若輸入缺失／不合格，不會寫評估快照。future date、NaN、bool age、未知性別、既往 CHD 及超出20–79歲皆拒絕。

依既有 state resolver 選擇 v3 的 `08_postprocessing/` 或 legacy state directory；新 PRS 檔名一律加個案前綴：

```
{LIS_ID}.prs.summary.json
{LIS_ID}.prs/assessments/{UUID}.json
```

每次評估使用新 UUID 原子保存：完整 bundle／summary／person、輸入 hash、建立時間、操作者、兩個風險與兩條曲線。未寫入 SNV review／ACMG／report state。summary 或模型、個案 metadata 改變時，不再顯示舊快照為目前結果；pipeline source marker 改變則要求重新確認並匯入 score。舊快照留作歷史稽核。前端改輸入即隱藏舊數值，切換個案會中止／丟棄舊 request，防止交叉顯示。

## 來源與重建

原 HR 來自作者 GitHub `MetaAnalysis/FullSample_MetaAnalysisandHetero_GoogleDrive.R` 引用的[原研究資料根目錄](https://drive.google.com/drive/folders/1bwedgU4lb4Y4i1pLsHXzjVaLBusjfcOQ)，其中 EstBB Full_cohort → Update_v3_211122：

- [男性年齡／性別 HR](https://drive.google.com/file/d/18ylwVsjMsaaQ-FASwHat_XuhTyItmI3i/view)：`HR_AgeandSexStratified_EstBB_211122.csv`。
- [女性 HR](https://drive.google.com/file/d/14IQVWEdxTvp2OloAji2637VDmI8VaCRC/view)：`HR_FemaleSample_EstBB.csv`。

同一資料根目錄的 MGB 候選檔與發表曲線未吻合，因此未採用；FinnGen 候選的年齡性別檔缺 Sex 欄，不猜測男女順序。EstBB 是本版唯一採用的 HR。

本機原始資料在 `research/prs/runtime/original_hr_20261102/`、GBD 與論文表在 `runtime/sources_20261102/`。builder 先檢查所有 input SHA256，再依精確 schema 解析；錯檔、缺組、驗算失敗皆不建立輸出包：

```bash
python3 -m research.prs.absolute_risk.build_taiwan \
  --sources research/prs/runtime/sources_20261102 \
  --hr-dir research/prs/runtime/original_hr_20261102 \
  --output-dir /tmp/tw_chd_v1_new
```

輸出目錄必須尚不存在。資料使用與再散布授權未因程式 MIT 授權而自動成立；原 HR／GBD 檔未提交 Git。

## 本機驗證

- 原文 Estonia CHD：224/224符合發表表格精度。
- Python：17個 API／輸入／snapshot 測試案例及12個風險引擎測試，共29項通過。
- 既有 frontend Node 回歸：93項通過。
- 獨立瀏覽器：桌面、390px手機、兩曲線、API數值一致、修改輸入清除舊值、>70歲、既往CHD、切換個案的晚到回應，全部通過。預覽使用明示的合成個案，不是實際病人結果。
- 未執行院內 fresh WGS、未測試院內部署、未進行台灣 clinical validation。
