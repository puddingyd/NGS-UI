# CHD PRS 研究風險 UI — v9.72

## 本版範圍

平台採兩階段；第一個疾病為 CHD，未建立 DM／cancer 的空白模型。

1. **三級分析**勾選「PRS」：00–08 發布完成後，自動排入原始 gVCF → 固定 v2.1.1 Model A/B score／percentile → 祖源 QC 的背景任務。不需要病人年齡或病史，不計算 absolute risk。PRS 啟動或執行失敗不回滾已完成三級結果。
2. **Secondary findings → PRS**：上方單行輸入「評估時年齡、性別、已有 CHD？」後按「風險計算」。日期自動記為當日；每次保存獨立評估快照。改臨床輸入只重算 risk，不重跑 gVCF。

主表每列一個疾病，固定欄位為「疾病／原始 PRS score／PRS percentile／未來10年風險／目前至80歲累積風險」。未計算的兩個風險欄顯示「請輸入臨床資訊後計算」；曲線放表格下方。未安裝 risk bundle 或祖源尚未通過時，合格 score／percentile 仍可獨立顯示，風險按鈕停用並說明原因。

表格的原始 score 指 **Model A total_sum（含固定模型缺失補值）**，與該列 EAS735 ECDF percentile 使用同一定義；不是先前610,677 loci的 observed-only SUM，也不是 pgsc_calc ancestry 的 normalized score。A/B完整結果仍保留於計分 output，B 暫不作UI風險主分析。

風險適用20至79歲、確認尚未有 CHD、具備合格 EAS 祖源證據與相符 model lock 的個案。>70歲不將截短期間冒稱10年風險，仍可計算至80歲。以「評估時存活且無 CHD」為條件，不是出生累積率相減。TWB／BBJ未納入；目前不收血壓、血脂、吸菸或家族史，因為本模型不使用這些欄位。結果不加入診斷／健檢 DOCX 或列印報告。

## 執行環境與首次啟用

程式本身已接上 queue／worker；**Git pull 不會安裝大型計分 reference，也不會自動複製原伺服器的 config**。使用原 v2.1.1 已驗證的設定檔，讓服務讀取：

- `NGS_UI_PRS_CONFIG`：既有完整 config.json 的絕對路徑；未設定時讀 `$NGS_UI_HOME/data/prs_config.json`。欄位見 `research/prs/scoring_v211/config.example.json`，勿直接沿用不存在的示範路徑。
- `project`：PRS 原始工作／結果 root。保留 `raw_work`、scorefile、model lock、Model A/B reference、FASTA、dbSNP、containers、兩份 Nextflow config 與 ancestry reference 的有效路徑。model lock 的絕對路徑鎖不可為換位置重新 freeze。
- 平台 worker 的 PATH 需可執行 `nextflow, apptainer, bcftools, samtools`，並具有原工作目錄讀寫權限。`raw_script` bytes 必須等於 repo 中已驗證的 vendor 腳本；既有十檔核心 checksum 仍由 pipeline_v21 驗證。
- `NGS_UI_PRS_BUNDLE`：第二階段的風險參數包，安裝方式見下節；第一階段不依賴它。

在服務相同使用者／環境中檢查（只讀，不執行個案分析）：

```bash
PYTHONPATH=backend python3 -c 'from app.services.prs_jobs import settings; c,h=settings(); print("PRS settings OK",h)'
```

三級視窗「PRS」預設未勾選；勾選時先檢查基本執行環境，缺設定立即回報，不會讓操作者等完整三級後才知道缺設定。每個樣本在三級發布後核對 gVCF；WES、找不到來源或多個候選均記錄到三級 Log／`prs_jobs` 子任務摘要，其他三級結果保留。載入個案後可在 PRS 區補跑或確認 gVCF。

### 原始 gVCF 定位

綁定 `pipeline_source.json` 中的 source sample、原 VCF anchor、來源類型與新保存的 `seq_type`，不從 UI 顯示姓名或 LIS 字串猜檔案。僅WGS；舊 marker 缺類型時可讀已登錄 metadata 的 test_type，仍不明則拒絕。

DRAGEN 預設精確搜尋同一 anchor 目錄或同 run 的 `gvcf.gz/` 中 `{source}.hard-filtered.gvcf.gz`／`{source}.gvcf.gz`。不同實體檔案不自動選最新；多候選由操作者在 UI 下拉確認。單候選自動使用。

In-house 沒有猜測的預設 caller／目錄；院內管理者確認原始gVCF後，在同一 config 加入 `gvcf_templates`。這個欄位會取代預設 mapping，需保留所需的 DRAGEN templates。只支援精確模板，不接受 glob 或 browser 提供任意檔案路徑：

```json
{
  "gvcf_templates": {
    "dragen": [
      "{anchor_dir}/{source}.hard-filtered.gvcf.gz",
      "{run_dir}/gvcf.gz/{source}.hard-filtered.gvcf.gz"
    ],
    "inhouse": []
  }
}
```

`anchor_dir` 是原 SNV VCF 的目錄，`run_dir` 是它的上一層，`source` 是原 sequencing ID。In-house 只有核對 caller、組裝及實際命名後才能填入模板。每個選定輸入都要求單樣本 header 精確符合 source ID、GRCh38 chr1長度248956422、gVCF reference-block allele宣告，以及 tbi/csi index。正式執行另由原 pipeline／bcftools 驗證資料。

### 排程、重試與保存

- 工作狀態／Log：`$NGS_UI_HOME/data/jobs/prs/{UUID}/`；保存該次 config snapshot、source／gVCF／resource fingerprints、操作者與 parent tertiary job。
- raw／Method2：原 `project/integrated_runs/ui_prs_{UUID}/{source}/`；大型中間檔保留在原 PRS project，不複製到08。
- ancestry：相同 raw 完成後，以原 pipeline.py 的 reuse 模式執行於 `ui_prs_{UUID}_ancestry/`。**不使用它的 normalized percentile 覆蓋 Method2 percentile**。原始 Method2 summary 不改寫；平台 envelope 另記 QC hash 並附上同源 ancestry evidence。
- 每次先原子發布已完成 score／percentile，再跑祖源。祖源失敗時顯示「計分完成，祖源檢查未完成」；score 保留，risk 不啟用。非 EAS／低信心亦不啟用風險。
- 全站同時只跑一個 PRS 計分工作，其他排隊。與三級／Somatic共用 UI及source sample locks，防止同個案來源替換。三級完成不代表其獨立 PRS 子工作已完成；在 PRS 區看進度／Log。
- 「重新分析／重試」在輸入 fingerprints 完全相同時沿用原 run ID，交給 v2.1.1 guarded resume；不刪除 raw audit。原 pipeline 不能恢復的 early partial run 可用「重新開始」建立新 run ID，原檔保留。
- 取消由 worker 收到 cancel flag 後終止自己啟動的工具 process group；原有 score／快照保留。完整刪除個案 pipeline 時，active PRS 阻擋刪除，完成工作狀態與08結果一起移除；原 project 的大型 raw audit 不在此刪除範圍。

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

## 相容入口：匯入既有計分結果（新流程不需要手動匯入）

```bash
PYTHONPATH=backend python3 scripts/import_prs_summary.py \
  --sample-id '實際已登錄個案ID' \
  --summary '/伺服器上既有PRS結果/summary.json'
```

`summary.sample` 必須符合 UI 個案 ID，或該個案 `pipeline_source.json` 內的 `source_sample_id`；不接受任意指定 alias。僅接受 `PGS004879_METHOD2_V0_20261031`、`METHOD2_TECHNICAL_PASS`。匯入保留原始 summary、原始 SHA256、來源識別與匯入時間，連結到當時 pipeline source marker 的 hash。

「成功匯入」不等於「適合風險計算」。UI 還會檢查 Model A、locus counts、percentile/group 一致性、model lock 與合格 EAS 祖源證據。原 pipeline 預設 raw-only 的 `NOT_ASSESSED` 無法估算；應取得既有流程的真實 ancestry QC，不得手改為 EAS 或偽造 PASS。

使用者自行更新／重啟平台後，載入個案，在「Secondary findings → PRS」填入年齡、性別與 CHD 病史後按「風險計算」。未安裝參數／未匯入 score 時，UI 明示缺件，不顯示合成數值。

## API 與保存

- `GET /api/samples/{sample_id}/prs`：登入保護；回傳 available、issues、score 摘要與有效的上次評估。
- `POST /api/samples/{sample_id}/prs/jobs`：啟動／重試；只接受伺服器列出的 candidate_id 與布林 new_run。
- `POST /api/samples/{sample_id}/prs/jobs/{job_id}/cancel`、`GET .../log`：個案綁定，禁止跨個案讀取／取消。
- `POST /api/samples/{sample_id}/prs/assessments`：只接受 `age, sex, assessment_date, prior_chd`。percentile、HR、參數路徑皆由伺服器讀取，不能由 request 指定。
- 若輸入缺失／不合格，不會寫評估快照。future date、NaN、bool age、未知性別、既往 CHD 及超出20–79歲皆拒絕。

依既有 state resolver 選擇 v3 的 `08_postprocessing/` 或 legacy state directory；新 PRS 檔名一律加個案前綴：

```
{LIS_ID}.prs.summary.json
{LIS_ID}.prs/assessments/{UUID}.json
```

每次評估使用新 UUID 原子保存：完整 bundle／summary／person、輸入 hash、建立時間、操作者、兩個風險與兩條曲線。未寫入 SNV review／ACMG／report state。summary 或模型、個案 metadata 改變時，不再顯示舊快照為目前結果；pipeline source marker 或 worker 記錄的 gVCF／index signature 改變則要求重新分析。舊快照留作歷史稽核。前端改輸入即隱藏舊數值，切換個案會中止／丟棄舊 request，防止交叉顯示。

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
- v9.72 Python：PRS API／queue／worker／引擎及相關三級回歸，共66項通過；外部WGS工具以隔離假工具驗證發布流程，非真實WGS驗收。原版計分27項離線測試亦通過。
- 既有 frontend Node 回歸：93項通過。
- 獨立瀏覽器：桌面單行輸入／疾病表格、390px手機表格捲動、待計算文案、兩曲線、API數值、修改輸入清除舊值、>70歲、祖源不合格仍顯示score、開始／Log／取消、切換個案晚到回應，全部通過。預覽使用明示的合成個案，不是實際病人結果。
- 未執行院內 fresh WGS、未測試院內部署、未進行台灣 clinical validation。

可重複執行的測試入口：

```bash
PYTHONPATH=backend python -m pytest tests/test_prs.py tests/test_prs_jobs.py research/prs/tests/test_absolute_risk.py tests/test_dragen_same_directory_batch.py tests/test_dragen_rerun.py tests/test_dragen_delete_output.py -q
TMPDIR=/private/tmp python3 -m unittest discover -s research/prs/scoring_v211/tests -q
node --test tests/frontend_*.test.cjs
# 需可載入 Playwright；CHROME_EXECUTABLE_PATH 可指定已安裝 Chrome。
# 測試自建本機靜態預覽，所有API為合成資料，不連正式平台。
node tests/browser_prs_ui.cjs
```
