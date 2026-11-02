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

## v9.74：GPN同步刷新與重疊請求

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
