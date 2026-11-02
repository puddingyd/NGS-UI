# 個案載入效能與排查

「資料載入中」涵蓋多個步驟，不能從這個視窗判定三級或 PRS 還在計算。載入新個案依序執行 POST `/api/samples`（登錄與院內 EMR、phenotype／panel設定）、GET `/api/samples`（索引）、GET `/api/samples/{id}`（核心 SNV資料）及 `/report`。核心載入完成後，Secondary SNV／CNV／SV／Mito／STR／PGx／ROH才各自背景讀取。

PRS 在個案載入時只重設顯示；展開 PRS 才讀取其 API。更新前啟動、未勾選 PRS 的三級工作不會自動補跑 PRS。v9.72 的 active-job 查詢包含 PRS job state，但不會因此啟動計分，也不會載入 PRS reference。

## 當下如何處理

若已載入完成，直接使用現有個案，不需刪除或重跑三級。若仍卡住，先確認三級清單狀態；另開分頁從既有個案清單搜尋同一個 LIS ID，已登錄者從既有清單開啟，避免重複送登錄請求。重新整理瀏覽器不會終止獨立三級 worker；重新整理也不保證原本的 HTTP request 會立刻停止於伺服器。

在瀏覽器開發者工具 Network 看仍 pending 的請求，或在伺服器讀取效能 Log。部署／SSH仍由操作者執行，以下只讀：

```bash
journalctl -u ngs-ui --since '1 hour ago' --no-pager | grep -E '實際LIS_ID|snv_review|sample.review_tsv' | tail -n 80
```

把 `實際LIS_ID` 換成目錄使用的個案ID，可能含 `-dragen`／`-nckuh` 後綴。不要貼姓名、病歷號、EMR內容或 secrets。

| Log | 意義 |
|---|---|
| `patient_store.register.phenotype_emr` | 登錄等待院內 EMR |
| `patient_store.register.write_version`／`compute_pheno_score` | 建立 phenotype／panel分數 |
| `router.samples.post.total` | 登錄 request總耗時 |
| `snv_review.rebuild_start`（v9.73起） | 準備掃完整raw；reason／changed_fields指出觸發條件 |
| `snv_review.ensure action=rebuild` | 完整重建review完成，含scanned／kept／耗時 |
| `snv_review.ensure action=hit` | 直接沿用既有review |
| `snv_review.ensure action=gpn-refresh` | 獨立補分或更新manifest；是否查GPN需看各annotation狀態 |
| `sample.review_tsv`／`snv.enriched_cache`／`sample.load` | review選擇、解析註解或個案組裝階段 |

只有取得該次 Log 才能歸因個別延遲。沒有現場證據時，不可把所有等待都歸因於 PRS，也不可假設某個修正已消除全部延遲。

## v9.73 修正範圍

WGS與TITAN-WGS使用相同的SNV review篩選規則。舊manifest的 `test_type=WGS` 在以 `TITAN-WGS` 開啟時，先前仍會因名稱不同而重掃完整raw，再重新做GPN annotation。小型重現測試已確認此行為；大型WGS的成本取決於檔案大小、磁碟與annotation資料。

新版在快取比對時視兩個名稱為同一政策，既有review bytes與mtime維持不變，manifest可更新實際標籤。WES不合併；raw大小／mtime、BED、overlay等仍嚴格核對，實際輸入變動仍重建。需要重掃時，先記錄原因，讓執行中也能定位。

目前的預建腳本本身已會將年＋T的個案記成TITAN-WGS，因此此修正主要保護舊／混合版本輸出及標籤轉換；不能僅憑個案名稱推定它一定遇到這個問題。

本機驗證：`tests/test_review_test_type.py`、`test_gpn_msa.py`、`test_gpn_star.py`、`test_snv_overlay.py` 共46項通過，涵蓋標籤互換不重掃、不重新GPN補分，以及WES深度政策／raw／BED真實變更仍重建。

## v9.74歷史方案：GPN同步刷新與重疊請求（UI背景補分已由v9.75取代）

現場效能紀錄確認一類延遲是GPN刷新：登錄約2.3秒，但核心資料request最久約742秒，其中review整理約738秒；另一重疊request的`action=gpn-refresh`約417秒，其後另一次約171秒。這些是重疊請求，不能相加當成單一耗時。完整5MB review的解析約2.3秒、之後cache命中約0.2秒；PRS讀取在核心資料完成後才發生。這組紀錄支持GPN刷新阻擋載入，並不支持v9.73的檢測名稱差異是此次主要原因。

舊程式沒有同案寫入鎖，且manifest使用共同`.tmp`路徑；多個request可能重複補分、互搶暫存或覆蓋彼此的provider欄位。紀錄出現OSError後沿用既有review，但當時未記例外，無法單憑該行斷定具體errno。舊`error=fallback_raw`標籤亦不可靠，實際selected仍是review；新版記錄`fallback_review`及exception。

本版處理：

- UI的core／secondary SNV讀取若base signature有效，立即使用已發布review；可選GPN更新改排背景。同process每案去重、同時一件背景工作，同案跨process由持久flock保護；其他同步reader取得鎖後重看manifest，避免重算。
- 更新MSA／Star M／V都在working copy，完成後一次發布；讀者可見舊完整或新完整註解。工作失敗／服務重啟會保留已發布結果；下次載入可再排程。背景完成後重新載入即可顯示新值。
- 只有成功complete的註解可忽略runtime-only差異：Star的Polars版本、MSA的tabix路徑。資料庫、index、manifest與shard signature仍嚴格比對；失敗／缺套件狀態仍可因環境修復觸發重試。`*_values`保留成功查詢的原始provenance，不以新的runtime覆寫。
- raw／BED／overlay／WES depth政策失效、review缺失或strict／force診斷仍需同步處理，不可用stale generation冒充有效結果；本版不保證所有首次載入永遠不需重建。
- 新增`snv_review.gpn_start`／`gpn_done`記各provider的觸發原因、改變欄位、狀態與耗時。`gpn-background`表示已排背景，`busy-reuse`表示使用有效的舊generation；`background status=complete`指工作結束，provider是否成功仍以各自annotation status為準。

當時為何首次判定GPN需要更新，仍需原manifest／執行環境或新版`gpn_start`紀錄判斷；runtime差異是已確認可觸發的不必要更新機制，不能憑現場耗時推定它就是當時的唯一觸發原因。

本機52項相關測試通過：原有lookup／版本刷新／overlay／WGS policy回歸，加上runtime-only不重查、背景排程、重疊請求只查一次、讀者不見M/V中間狀態、background失敗釋放名額、raw變更時不得盲用舊review及strict／force不被延後。未SSH、未部署；院內實際載入耗時需更新服務後驗證。

## v9.75歷史方案：載入只讀（v9.76已恢復版本刷新）

工作保存的Research-only=true能證明系統收到該設定，不能推斷使用者是否誤按。缺少分數本身不應啟用Star。舊loader依review manifest的持久opt-in允許lazy refresh，因此即使資料庫沒更新，也可能因CLI／web環境差異、先前skipped狀態或其他signature差異觸發查詢；必須看三級GPN輸出／當時manifest，才能區分三級未成功與載入誤判過期。原有耗時紀錄只能確定載入時做了GPN刷新，無法單獨判定最初觸發欄位。

v9.75當時行為（已由下節取代）：

- 載入core／secondary SNV以`compute_gpn=False`取得review，完全不檢查GPN資料庫或查分，不會因缺分、套件／版本差異或新資料庫自動補跑。v9.74的UI背景排程已移除。
- 未勾Research-only：Star不計算；review manifest未明記opt-in=true時，payload不顯示殘留的M/V分數，MSA仍可顯示。磁碟原始分數不因此改寫或刪除。
- 真正raw／BED／篩選政策改變時，仍可在載入重建篩選列，GPN標為`pending_postprocessing`（未啟用Star維持skipped）；由三級或明確補建補齊。不能拿舊generation掩蓋來源變更。
- 三級使用與web相同的Python直譯器，並將該bin加入PATH，讓`env python3`的post-processing腳本能取得相同Polars環境；原工具PATH其餘部分保留。
- 新pipeline_source記`research_only`、`job_id`。三級post-processing結束時驗證旗標相符並保存`gpn_annotations`；所需provider未complete則記`annotation_warnings`，三級視窗顯示未完成原因。原best-effort保留，但使用者不再靠載入隱性補跑。
- 同案寫入鎖、工作副本發布、成功lookup忽略runtime-only差異、逐provider觸發原因與耗時仍保留，適用CLI與明確計算路徑。

本機89項相關Python測試與93項既有前端回歸通過，涵蓋資料庫變動仍不允許UI查GPN、缺review重建不跑predictor、非research隱藏殘留分數、三級旗標追溯與未完成提醒。未SSH／未部署；既有已完成個案無須為此重跑整套三級。

## v9.76：根因證據與恢復正確版本刷新

補齊原三級Log後，已確認原工作在review建立時M／V皆為`skipped_missing_polars: No module named polars`、annotated_rows=0，MSA則complete。這代表三級根本沒有成功產生Star分數；載入web環境具備Polars，signature的polars由空字串變成已安裝版本，才觸發原本允許的補分。當時不是資料庫無更新卻無故重算complete結果。Research-only工作設定true只能確認系統收到該設定，不能據此判斷是否誤按。

此次修正對準執行環境，並撤回v9.75的正常載入只讀限制：

1. worker沿用web的`sys.executable`，並設定`NGS_UI_PYTHON_BIN`；shell review步驟明確用此Python執行builder。即使shell PATH的python3沒有Polars，review仍使用正確環境。builder先印`[gpn-runtime]`的直譯器、Polars版本與Research-only旗標。獨立手動呼叫shell時可設定同一環境變數，未設定則沿用其python3。
2. core／secondary正常載入恢復計算與刷新：同一份raw、篩選、DB與runtime signature，三級完成後立即載入必須hit；GPN資料庫／index／manifest／shard／Polars套件版本變動時更新對應provider，不再忽略套件版本差異。compute_gpn=False僅給明確只讀操作，正常UI不使用。
3. 缺分本身不會開啟Research-only；Star仍只在持久opt-in=true下刷新。MSA獨立，不隨Star開關消失。非research殘留Star值不顯示。
4. 同案鎖、完整工作副本發布與三級未完成提醒保留；不引入UI背景補分。每次provider查詢前記`signature_diff`實際葉欄位的before／after，最多12項並標示截斷，之後可直接追出polars、path或mtime差異。

重現測試採無Polars的isolated Python執行原review CLI，證實M／V被略過；再以有Polars的測試環境載入時，只刷新缺的M／V，所有DB檔案bytes與mtime皆未變，signature差異只有polars。修正回歸則把shell PATH指向無Polars環境、明確指定平台Python，實際跑run_stopgaps後立即載入，禁止任何再次查詢並驗證review bytes／mtime不變。另以套件版本改變測試保證仍會更新M／V。

這些都是合成資料／小型Parquet測試，沒有執行院內病人WGS，也未SSH／部署。已成功補齊的既有個案不需為本修正重跑三級；更新服務後，後續三級工作會使用一致的執行環境。

原三級Log沒有記錄直譯器路徑，因此它能確定的是「三級當時無法import Polars，平台補分時可以」，無法再細分當時是否為兩個不同Python，或期間曾安裝套件。舊啟動方式以python3／PATH尋找直譯器確有混用風險；本版明確指定同一Python並記錄版本，消除這條風險並留下可驗證證據。

本版90項相關Python測試通過（含PRS任務與三級回歸）；未執行院內真實WGS或部署。
