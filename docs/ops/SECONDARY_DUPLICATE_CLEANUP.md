# DGX／DGM 二級重複資料：盤點與確認清理

v10.34，2026-11-13。這是獨立操作者工具，不改平台服務／分析設定，也不代替使用者到 server 執行。

## DGM：本次只搜尋

開啟`scripts/paste_cleanup_secondary_duplicates_dgm.sh`，**全選整份內容**，貼到DGM互動終端。完整程式已包含在檔案中，不需要先複製repo、沿用DGX的`/tmp`檔案或安裝套件。明確指定：

```text
二級：/home/datalake_Intermediate/pipeline/nextflow_output
三級：/home/datalake_Intermediate/pipeline/tertiary_output
```

入口先印執行帳號，將程式與`plan.json`／`duplicates.tsv`放在DGM新建的`/tmp/ngs-secondary-review-dgm.XXXXXX/`，**不傳`--review`，所以只有盤點，不出現刪除確認、也不執行刪除**。清單的DELETE僅表示通過本次檢查的待確認候選。修改日期／11月5日界線／最新輸出與相依保護沿用下列規則。終端同時印出實際二級及三級根，方便核對。

不使用sudo、setfacl或chmod。DGX回報`setfacl: Operation not supported`表示該掛載不支援ACL；換到DGM後能否讀取700的目錄，取決於執行帳號是否為擁有者或有其他有效讀取權限，不能只因資料在DGM建立就認定可讀。仍遇到讀取錯誤時，保留所有待刪候選，附出錯路徑；不要略過PRS相依子目錄。

本次先將DGM產生的完整報告交回確認。DGM入口沒有自動進入清理的後續步驟；既有DGX確認碼／報告也不能跨主機沿用。

## DGX：盤點後可確認清理

開啟 `scripts/paste_cleanup_secondary_duplicates.sh`，**全選整份內容**，貼到 DGX 互動終端。DGX 不必有這個 repo、不需連 GitHub、不需安裝 Python 套件。使用 Python 3.8+ 標準函式庫及 Linux 的 `ps`。

程式固定傳入使用者指定的 root：

```text
/datalake_Intermediate/pipeline/nextflow_output
```

先在 `/tmp/ngs-secondary-review.XXXXXX/` 放程式和報告；資料根先唯讀。列出每個重複 sample 的每次sample資料夾修改時間（台北，判斷依據）、完整路徑、KEEP／DELETE 與原因。完整紀錄在 `duplicates.tsv` 和 `plan.json`。

有可刪項目時最後詢問：

```text
輸入 DELETE <本次報告的12碼確認碼> 才永久刪除；Enter或其他輸入一律不刪。
```

如果只是要先盤點、稍後將報告交回討論，**直接按 Enter**。不要輸入 DELETE。程式與報告留在印出的 `/tmp` 目錄；不會自動刪除資料。若要稍後重看及確認，可在 DGX 執行印出目錄內的程式，重新盤點：

```bash
python3 /tmp/ngs-secondary-review.實際代碼/cleanup_secondary_duplicates.py --review
```

也可直接使用 repo 內 `scripts/cleanup_secondary_duplicates.py`；**省略 `--review` 時一律只有盤點**。`--root`預設仍為DGX二級根，DGM必須明確指定/home路徑，或直接使用上述完整DGM可貼版。`--report-dir` 指定報告位置（必須在資料根以外），`--reference-root` 指定實際本機三級資料根，`--work-root` 指定共用 Nextflow work root；預設分別為同層的 `tertiary_output` 與 `/raid/DGM/work`。

## 本次規則

- 只辨識既有 `<root>/<batch>/<sample>/02_alignment` 或 `04_snv_indel` 的兩層 sample 目錄；sample ID採完整字串相同，不合併 S1 與 S10，也不刪整個混有其他個案的 batch。
- 依使用者確認的新規則，每個重複sample保留**sample資料夾修改時間（mtime）最新者**；相同mtime並列最新者全部保留。這是`<root>/<batch>/<sample>`目錄本身的mtime，不用batch名稱日期，也不以子檔案的最大mtime取代。
- 日期界線保守設為 **2026-11-05 00:00:00 +08:00**：**sample資料夾修改日期為11/5當天及以後的版本全部保留**。因此不保證最後每個sample只有一份；這是使用者指定保護例外。
- 不查birthtime，也不再因建立時間未知而整組保留；不把修改日期稱為建立日期。
- mtime是本次經確認的排序／日期界線；ctime與整個樹的metadata仍存入確認用指紋，僅用於偵測「盤點後又有變更」，不再因ctime或子檔案時間晚於11/5就額外改變日期篩選規則。
- 同ID出現明確 WES 與 WGS batch標籤時，整組先保留，避免清掉不同檢測類型。標籤判斷只看batch名稱中的完整WES／WGS詞。
- 最新版本必須具備非空的對應BAM、VCF及BAI／TBI或CSI，才把更舊版本列為候選；這是最低檔案檢查，不代表全部工具或臨床品質已驗證。

前次操作者貼回145個重複sample、298個目錄全部birthtime未知，故v10.31只有盤點，沒有刪除；使用者隨後明確確認改採「資料夾修改日期」。新版JSON記`date_basis=sample_directory_mtime`，終端與TSV也標示同一依據。先前保存在`/tmp/ngs-secondary-review.w6rLy9/`的程式不會自動更新，請重新複製新版可貼程式，重新產生清單與確認碼。

## 相依資料與執行中工作

候選目錄會核對：

- 二級root及三級root內的symlink是否指回該目錄（不遞迴跟隨symlink）。
- 三級`08_postprocessing/*pipeline_source.json`的`source_path`／`source_vcf_path`是否仍引用該版本。只針對已知的共享掛載，依選定root將DGM的`/home/datalake_Intermediate/...`及DGX的`/datalake_Intermediate/...`對應為當前主機寫法再比較；支援兩個方向，亦適用另一主機路徑在本機不存在的symlink。採完整路徑元件界線，不改寫相似名稱或任意自訂root。

被引用者保留。相依root沒掛載或無法完整核對時，仍輸出重複清單，但所有待刪項目改為保留。

**執行刪除前，操作者需確認DGX及DGM都沒有分析工作，清理期間不啟動新工作。** DGX程式會再檢查可見Nextflow／`.command.run/.sh`程序，並取得現有`/raid/DGM/.work.ngs-ui-locks/admission.lock`、確認其他scope leases皆閒置。無法讀取／取得鎖就停止；不自行啟停任何服務或程序。DGX無法從本機ps證明DGM或其他主機沒有工作，因此終端確認文字也包含操作者對跨主機閒置的確認；手動繞過平台鎖的啟動不得並行。

## 刪除與紀錄

確認後，在共用admission鎖內重新盤點。日期、輸出、檔案樹metadata或相依關係任一變動，都要求重新檢視，不沿用舊確認。每筆刪除前再核對舊目錄及最新保留版本。metadata digest不讀取整個BAM／VCF內容；核心輸出亦追蹤symlink實際目標的metadata。

只刪已標記DELETE的batch/sample，不刪FASTQ原始資料、不刪三級、不清Nextflow work、不刪batch samplesheet或pipeline_info。rmtree不跟隨symlink；內含掛載點拒絕刪除。實際刪除前先將該sample原子移到同batch下的私人暫存目錄並核對inode，再刪除；這是操作期間的暫存，**不是可復原垃圾桶**。

`actions.jsonl`在移動前就記原路徑／暫存路徑，逐階段flush及fsync。成功記deleted；失敗或中斷保留剩餘暫存路徑並停止，不繼續刪下一筆，不聲稱可以回復已刪內容。

清理後可在平台的檔案清單按一次「更新」，或等原本每日更新；這個工具不直接寫平台索引。批次samplesheet與Nextflow work保留，因此之後刻意重跑舊batch仍可能重新產生那些sample輸出。

## 驗證與界線

合成資料測試涵蓋資料夾mtime最新／日期並列、11/5界線、完全不呼叫birthtime查詢、ctime與子檔日期不取代sample目錄mtime、最新核心檔不齊、WES/WGS分開保留、symlink／三級marker引用、DGM及DGX兩方向掛載引用保護／相似前綴不誤判、相依子目錄權限錯誤全部保留、DGM唯讀入口不詢問不刪除、最新版本變動時拒刪、實際rmtree只刪核准sample、外部link目標與同batch其他sample保留、Enter不刪、busy lease阻擋、相依root缺失不刪、兩份可貼版與受測Python原碼一致。測試只在本機合成暫存目錄；沒有到DGX或DGM搜尋或刪除真實資料。
