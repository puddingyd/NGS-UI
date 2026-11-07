# 六疾病 PRS：準備、啟用與驗收（v10.02）

2026-11-05。CHD 與新增五疾病共用一次 gVCF 重建、一次六欄計分及一次祖源分析。五個新疾病一併啟用；沒有逐疾病開關。程式、台灣風險參數及 UI 已納入 repo。依操作者11/5貼回Log，院內共用模型準備與CHD 3,330人回歸已通過，26T00970完成六疾病計分；祖源因輔助檔表頭缺欄位失敗並已取消。**v9.92修復此格式問題；完整院內祖源／風險／報告驗收仍待重試確認，不能將本機測試視為院內流程已驗收。**

## v9.92 修復後如何重試

已完成`setup_prs_multi.py --activate`者不需重下準備命令，也不用改已安裝的`ancestry_score.txt.gz`或manifest。由操作者更新程式、重啟平台後，在原個案PRS卡片按「重新分析／重試」。本次程式版本改變會自動選新run並重做個案計分；舊正式結果及audit保留，不需重跑整個三級，也不重新下載或做3,330人校準。同版本且同輸入的後續重試仍沿用checksum checkpoint。

祖源開始前，程式在個案`multi/ancestry/`產生格式完整的輔助檔，變異列與已核對模型完全相同，只補metadata；留來源／輸出SHA256，通過`pgscatalog-format`預檢才啟動Nextflow。預設使用`plink_img`同層既有`ghcr.io-pgscatalog-pygscatalog-pgscatalog-utils-2.0.4.img`；若實際檔名／位置不同，在PRS設定加入絕對路徑`pgscatalog_utils_img`。容器與格式檢查暫存均沿用個案工作區。

若Nextflow回報確定的workflow失敗，平台會顯示失敗原因及收尾狀態並停止該工作process group，再保留score／percentile、標記祖源未完成；不允許算risk。一般Nextflow task retry warning不會觸發停止。摘要Log保留ValidationError的欄位與原因；完整工具Log仍可下載。新祖源checkpoint另核對輔助檔、來源sidecar與格式檢查紀錄。

## 操作者先做什麼

1. 由操作者依既有程序更新程式、重啟平台；Codex 不執行部署。既有 CHD 設定及結果仍可使用。
2. 在 repo 根目錄、使用與平台相同的 PRS 設定及工具環境執行：

   ```bash
   PYTHONPATH=backend .venv/bin/python scripts/setup_prs_multi.py --activate
   ```

   設定預設讀 `NGS_UI_PRS_CONFIG`，未設時讀 `NGS_UI_HOME/data/prs_config.json`。若 systemd 中設定的環境變數未帶入登入 shell，請加 `--config /實際路徑/prs_config.json`。需使用原來成功執行 CHD 的 `bcftools`、`samtools`、`apptainer`、Nextflow 環境，另需 `zstd` 解開參考資料。

   工具會下載並核對五份 PGS Catalog GRCh38 scorefiles（壓縮檔合計約 84 MB），讀原 CHD lock 及原 HGDP+1kGP archive，建立六疾病固定模型。完整參考 PGEN、解壓副本、SQLite 及計分輸出遠大於 scorefile 本身，準備期間需保留工作磁碟空間。若已有解開的完整 GRCh38 PGEN，可加 `--reference-pfile /完整前綴`（不含 `.pgen`），省去解壓；不能使用僅有 CHD 位點的子集。

   只有六個模型及 CHD 3,330 人參考分數回歸檢查全部通過，才寫 `manifest.json`。`--activate` 另檢查無執行中的三級／PRS 工作，備份現有設定，再原子寫入 `multi_model_dir`；不啟動任何個案、不重啟服務。若模型已完成但當時有工作執行，工作結束後重下同一指令即可啟用，不重建模型。
3. 載入既有 WGS 個案，在 PRS 卡片按重新分析；不必重跑整個三級。之後新的三級工作勾 PRS 時會使用同一個六疾病流程，仍等 PRS 結束才顯示 100%。先依下列三案完成回歸驗收，再用於研究結果輸出。

如果準備中斷，已下載的 scorefiles／已完成的參考解壓可沿用；尚未完成的模型目錄不覆寫。保留失敗稽核，以新的 `--output` 目錄重試。不得直接修改 manifest 把失敗改成 PASS。需要退回舊 CHD 計分時，在無 active job 時由操作者還原 `.before-multi-*` 設定備份；既有個案的六疾病結果保留，但設定不再相符時不允許重新計算／匯出六疾病風險。

## 疾病、模型與基準資料

| 疾病 | PGS | GBD 2019 cause | EstBB 相對風險分層 |
|---|---|---|---|
| 冠心病 | PGS004879 | Ischemic heart disease | 男性年齡＋性別；女性性別 |
| 第二型糖尿病 | PGS004870 | Diabetes mellitus type 2 | 年齡＋性別 |
| 整體癌症 | PGS004875 | Total cancers | 年齡 |
| 乳癌 | PGS004869 | Breast cancer | 完整樣本 HR；風險僅適用女性 |
| 心房顫動 | PGS004878 | Atrial fibrillation and flutter | 年齡＋性別 |
| 癲癇 | PGS004881 | Idiopathic epilepsy | 年齡 |

五個新增疾病都有作者使用的 GBD 對應資料，本版讀其台灣、2019、年齡別、性別別資料，不使用年齡標準化率。各疾病均用發生率、盛行率、疾病死亡率及全死因死亡率，按與 CHD 相同的競爭風險公式計算。疾病代碼／命名與 GBD 並不完全等同；具體限制由 registry 隨 UI 與報告輸出。

所有模型都使用固定 3,330 名參考者的平均基因型補值，分別與相同 735 名東亞參考者的該疾病分數比較。**同一批參考者不代表共用一個分數分布**。CHD 固定 610,305 位點與原參考分布不變；新增模型在完整參考 PGEN 上核對等位基因，方向轉換須有 dbSNP 身分支持，不唯一的匹配拒絕。每個模型至少 99.5% 原始位點須匹配，實際固定位點數寫入 manifest，不在驗收前猜測數量。

本版平均補值與東亞百分位屬平台固定參考校準，不是原研究歐洲 cohort 的個人分數。歐洲 EstBB HR 外推台灣及這套分位組對應，**尚未完成台灣實際疾病追蹤驗證**；不整合 TWB，不需要院內追蹤資料才可執行，但也不因此變成臨床驗證模型。

## 共用流程與路徑

```text
原始 gVCF
  → 六模型位點聯集；每條染色體只掃一次，DP≥10／GQ≥20
  → 一個 union PGEN，逐位點核對匯入／匯出基因型
  → 一次 PLINK 六欄計分，另以 Python 權重加總獨立核對
  → 六個 score／各自 EAS735 percentile，發布至 08
  → 一次祖源分析（固定 CHD canonical 位點，避免新增模型改變 PCA 位點政策）
  → 每疾病品質及共用祖源 gate；三級後處理結束
  → UI 填年齡／性別／各疾病病史 → 風險 → 保存快照 → 健檢報告
```

| 內容 | 預設位置 |
|---|---|
| 準備完成的共用模型 | `NGS_UI_HOME/biotools/prs/multi_v1/` |
| 一次性下載／參考解壓工作區 | `NGS_UI_PRS_WORK_ROOT/_reference_prepare_v1/`；未設 override 時為 `NGS_UI_HOME/nf_work/prs/_reference_prepare_v1/` |
| 個案重建、audit、Nextflow | `nf_work/prs/{個案ID}/integrated_runs/{run}/{source}/multi/` |
| 個案容器暫存 | `nf_work/prs/{個案ID}/.runtime/tmp/` |
| 正式 score／QC／風險快照 | 原個案 `08_postprocessing/`，沿用既有 prefixed resolver 與不可變 artifacts |
| 台灣風險參數 | repo `research/prs/releases/tw_multi_v1/`，發布檔 SHA256 固定核對 |

原 PRS 專案僅讀取既有 CHD／FASTA／dbSNP／container／reference 資源，所有新個案產物不寫回原專案。新增流程沿用既有容器 temp bind、sample lock、PRS slot、work lease 與三級取消機制。

2026-11-05 核對「規劃磁碟重新配置方案 A」（對話 `01a1a057-3959-7473-8c83-96bf63bcbfc8`）最後的操作者紀錄：11/3 已將 `/scratch/n102968/nf_work` bind mount 回 `/home/n102968/NGS_UI/nf_work`，`findmnt` 指向 NVMe 的 `ubuntu--vg-ubuntu--lv`，服務的 mountpoint／device 檢查通過。因此上表的個案、容器暫存與 `_reference_prepare_v1` 都隨既有掛載使用 SSD；`biotools/prs/multi_v1` 是共用模型資源，不屬於 nf_work。**維持原 work root 設定，不另外改成 `/scratch`，避免改變續跑 fingerprint。** 本次未 SSH 或重新掛載；[搬遷與回復文件](TERTIARY_SCRATCH_MIGRATION.md)保留原操作細節。

染色體重建及計分完成各有 checksum checkpoint。輸入／程式／模型／設定變動不沿用舊 checkpoint。計分成功、祖源失敗時仍可顯示六個 score／percentile，但禁止風險；重試沿用已核對的計分結果。Log 中標示疾病的 score／percentile 取正式表格結果；Nextflow 輔助計分不取代表格值。

## UI 與報告

- PRS 維持 PGx 前方，一列一疾病、六列固定順序；舊 CHD 個案保留原列，其餘顯示尚未分析。
- 最上方輸入評估年齡、性別、各疾病病史（是／否／尚未確認）。依使用者要求，v9.91 新表單各疾病預設「否」，已儲存的「是」或「尚未確認」會保留；仍須按「風險計算」才產生評估。
- 表格順序：疾病／原始 score／percentile／未來10年／未來20年／目前至80歲。20年與10年從同一評估年齡起算，採同一競爭風險積分；20年不是10年乘二。20至60歲可提供三種期間，>60至70歲提供10年及至80歲，>70至未滿80歲只提供至80歲。不能把超出80歲終點的期間截短冒稱完整10／20年。
- 各疾病獨立判斷：已有疾病不估復發、未知病史不計算、計分品質不足不計算。男性不提供乳癌 absolute risk；乳癌已患病與整體癌症未患病的矛盾輸入會被拒絕。
- 年齡、性別或病史變更即清除畫面舊風險及全部曲線，需重新按風險計算。v9.95起下方六疾病曲線同時呈現，桌機3欄2排（順序同表格），≤1050px為2欄、≤640px為單欄；不再用下拉選單。每圖保留個案／40–60參考組及展開數值表；不適用或未完成分析者保留原位置並說明原因，不產生曲線。各圖縱軸依疾病風險調整，須讀取百分比刻度，不用圖形高度直接比較疾病。本次只調整呈現，已有評估重新載入即可，不需重算。
- 健檢報告仍用一個「多基因風險評估」選項，位於 PGx 前。至少一疾病有目前資料對應的完整風險快照才預設勾選；v10.02起表格保留全部已計分疾病的百分位，不適用原因直接寫在該列風險欄，患病者顯示「已有此疾病，不估算風險」。採與其他章節相同的`====`等寬文字表，模型清單包含全部已計分疾病。文字與Word段落設定依20則註解調整，詳見[逐則對照](PRS_HEALTH_REPORT_COPY_REVIEW.md)。已有有效快照重新匯出即可，不需重算。
- 新健檢表為疾病／百分位／10年／20年／至80歲五欄。v9.90以前快照不被改寫；UI保留其原結果並提示補算20年，報告須重新按一次「風險計算」產生新版快照才可勾選。新評估在>60歲明確保存`risk_20_year=null`，屬合法不適用，與舊快照缺欄位不同。不需重跑gVCF或重新準備模型／風險參數包。
- 整體癌症與乳癌事件重疊，風險不可相加。各疾病模型也不計算同時罹患多病的聯合機率。

## 參數重現與資料限制

使用原 HR／GBD，以原研究五年離散算法重現 Supplementary Data 4 的 Estonia 結果；發表數字僅作驗算，不用來擬合 HR。UI 則沿用 CHD 的區間內連續競爭風險積分，支援非五歲整數年齡。兩種數值方法不可宣稱逐值相同。

| 疾病 | 驗算筆數 | 最大誤差（百分點） |
|---|---:|---:|
| CHD | 224 | 0.0049951003 |
| T2D | 224 | 0.0049983042 |
| All cancers | 224 | 0.0049869211 |
| Breast cancer | 112 | 0.0049325866 |
| AF | 224 | 0.0049995434 |
| Epilepsy | 224 | 0.0049904313 |

全部 1,232 筆落在發表數值四捨五入至小數兩位的容差（0.0050001 百分點）內。未重現信賴區間，個人不確定性區間仍為 null。

**整體癌症的分組方向**：公開 EstBB age-stratified CSV 的分組方向與發表風險曲線相反。本版以 `<20→>80%`、`20–40→60–80%`、`40–60→40–60%`、`60–80→20–40%`、`80–90→10–20%`、`90–95→5–10%`、`≥95→<5%` 對應原始 HR 標籤，重現全部 224 筆結果。這是從公開資料及發表結果推得的來源對應，不宣稱作者明文記載此轉換；沒有倒數 HR、調整係數或放寬容差。來源對應與證據寫入 bundle，也在 UI／報告揭露研究性重建限制。

**AF 年輕區間**：固定 commit 的原始 GBD 六個 CSV 中，台灣男性／女性 20–24、25–29 歲發生率及盛行率明列為零；這四格疾病死亡率則缺列，由作者處理流程補零。版本內保存原始資料的台灣 500 筆相關列與六檔 SHA256，每個基準率逐格核對，只有上述四格死亡率允許作者補零。此期間模型的零風險不代表個人不會發病；不把任意缺列當零。

來源、原始 URL、checksum、逐列驗算見 [發行資料說明](../../research/prs/releases/tw_multi_v1/README.md)。

## 本機驗證與院內待驗收

本機已測：條件風險及七分組、疾病個別 gate、性別／病史／年齡、快照防竄改及過期、舊 CHD 相容；使用真正 bcftools／samtools／PLINK 的 22 染色體小型合成 gVCF 六欄計分，包含多等位、低 MIN_DP、平均補值、PGEN 往返及 checkpoint 竄改拒絕；worker 先發布計分再發布祖源、祖源失敗保留計分。桌機／手機 Chrome 測試與合成男性／女性健檢 DOCX 逐頁檢查完成。合成資料不是病患結果。

v9.90 驗證紀錄：相關 Python 測試 111 passed／5 skipped（五項需真實 Nextflow 環境）、前端單元測試 6 passed，另有完整 Chrome 操作測試。五份官方 scorefile 已全量解析及核對 checksum／同模型等位基因唯一性：乳癌 872060、T2D 930506、整體癌症 818872、AF 785779、癲癇 605432 列；這些是原 scorefile 列數，**不是**尚未準備完成的固定參考匹配數。

v9.91 驗證紀錄：20年期間、60／70歲邊界、舊快照與報告相容性等相關 Python 測試 56 passed、前端單元測試 4 passed；Chrome 操作測試確認六疾病病史預設「否」、已儲存病史保留、新舊評估顯示及手機橫向捲動。男性／女性合成健檢 DOCX 共六頁逐頁檢查完成。

v9.92 驗證紀錄：相關Python測試117 passed（含官方pgscatalog.core 1.0.2格式工具重現舊錯誤／新檔通過，以及真實Nextflow 24.10.4失敗／pending task清理）、前端單元4 passed及Chrome桌機／手機整合通過。驗證舊模型壓縮檔與變異列不變、來源／列數不符拒絕、格式錯誤不啟動Nextflow、程式變更自動新run與過往Log錯誤不誤殺新步驟；本機未執行院內Linux容器或完整WGS，院內啟動時仍須通過實際容器的格式預檢。

11/5操作者紀錄：完整參考模型準備及3,330人CHD回歸已通過；26T00970 CHD表格0.366424／35.10百分位、609032可靠位點，與下列基準相符。尚待院內確認：修復後祖源與容器預檢、其餘真實WGS六模型匹配與完整耗時／記憶體／磁碟、Nextflow中斷恢復、風險／報告及三級整合。部署後以原三案核對 **CHD 表格正式值**，不是 Nextflow 輔助分數：

| 個案 source ID | CHD score（含固定缺失補值） | EAS percentile | observed／missing |
|---|---:|---:|---:|
| 26T00970 | 0.366424148909 | 35.102041 | 609032／1273 |
| 26T00969 | 0.474838562769 | 60.544218 | 609162／1143 |
| 26T00968 | 0.556804458045 | 79.455782 | 609272／1033 |

來源是凍結版 `VALIDATION_REPORT_v211_2026-11-02.md` 的固定 610305 位點主分析；score 顯示六位及 percentile 顯示兩位須與舊 UI 一致，未四捨五入值則與已保存的原 summary 比對。若不同，保留 audit 追查，不自動覆寫基準或放寬 gate。新增五模型需確認數值有限、匹配率、正反股／高權重缺失及分組邊界。最後驗收三級等 PRS 完成才 done，及 UI 重算風險、報告前後一致。
