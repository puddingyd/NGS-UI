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
