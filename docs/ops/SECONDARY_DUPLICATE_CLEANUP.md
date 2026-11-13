# DGX 二級重複資料：盤點與確認清理

v10.31，2026-11-13。這是獨立操作者工具，不改平台服務／分析設定，也不代替使用者到 server 執行。

## 執行方式

開啟 `scripts/paste_cleanup_secondary_duplicates.sh`，**全選整份內容**，貼到 DGX 互動終端。DGX 不必有這個 repo、不需連 GitHub、不需安裝 Python 套件。使用 Python 3.8+ 標準函式庫及 Linux 的 `stat`／`ps`。

程式固定傳入使用者指定的 root：

```text
/datalake_Intermediate/pipeline/nextflow_output
```

先在 `/tmp/ngs-secondary-review.XXXXXX/` 放程式和報告；資料根先唯讀。列出每個重複 sample 的每次建立時間、修改時間（另行標示，不冒稱建立時間）、完整路徑、KEEP／DELETE 與原因。完整紀錄在 `duplicates.tsv` 和 `plan.json`。

有可刪項目時最後詢問：

```text
輸入 DELETE <本次報告的12碼確認碼> 才永久刪除；Enter或其他輸入一律不刪。
```

如果只是要先盤點、稍後將報告交回討論，**直接按 Enter**。不要輸入 DELETE。程式與報告留在印出的 `/tmp` 目錄；不會自動刪除資料。若要稍後重看及確認，可在 DGX 執行印出目錄內的程式，重新盤點：

```bash
python3 /tmp/ngs-secondary-review.實際代碼/cleanup_secondary_duplicates.py --review
```

也可直接使用 repo 內 `scripts/cleanup_secondary_duplicates.py`；**省略 `--review` 時一律只有盤點**。`--report-dir` 指定報告位置（必須在資料根以外），`--reference-root` 指定實際 DGX 三級資料根，`--work-root` 指定共用 Nextflow work root；預設分別為同層的 `tertiary_output` 與 `/raid/DGM/work`。

## 本次規則

- 只辨識既有 `<root>/<batch>/<sample>/02_alignment` 或 `04_snv_indel` 的兩層 sample 目錄；sample ID採完整字串相同，不合併 S1 與 S10，也不刪整個混有其他個案的 batch。
- 每個重複 sample 保留建立時間最新者；建立時間並列的最新版本全部保留。
- 日期界線保守設為 **2026-11-05 00:00:00 +08:00**：11/5當天及以後建立的版本全部保留。因此不保證最後每個sample只有一份；這是使用者指定保護例外。
- 建立時間使用真正的檔案系統 birthtime：macOS的`st_birthtime`或Linux GNU `stat %W`。任何版本取不到建立時間，就保留整組，**不拿批次名稱／mtime／ctime猜建立時間**。修改時間另外列出，方便核對。
- 舊目錄即使建立於界線以前，只要樹內任一項目的mtime／ctime仍在11/5起有變動，也保留，避免清掉舊目錄內新產生／覆寫的內容。
- 同ID出現明確 WES 與 WGS batch標籤時，整組先保留，避免清掉不同檢測類型。標籤判斷只看batch名稱中的完整WES／WGS詞。
- 最新版本必須具備非空的對應BAM、VCF及BAI／TBI或CSI，才把更舊版本列為候選；這是最低檔案檢查，不代表全部工具或臨床品質已驗證。

如DGX檔案系統沒有birthtime，報告會列「建立時間無法取得」並保留，不會因為人已按確認就猜測日期強刪。這時先檢視盤點報告，再由使用者決定其他日期證據與新的清理規則。

## 相依資料與執行中工作

候選目錄會核對：

- 二級root及三級root內的symlink是否指回該目錄（不遞迴跟隨symlink）。
- 三級`08_postprocessing/*pipeline_source.json`的`source_path`／`source_vcf_path`是否仍引用該版本。DGM記錄的`/home/datalake_Intermediate/...`會對應成DGX的`/datalake_Intermediate/...`再比較。

被引用者保留。相依root沒掛載或無法完整核對時，仍輸出重複清單，但所有待刪項目改為保留。

**執行刪除前，操作者需確認DGX及DGM都沒有分析工作，清理期間不啟動新工作。** DGX程式會再檢查可見Nextflow／`.command.run/.sh`程序，並取得現有`/raid/DGM/.work.ngs-ui-locks/admission.lock`、確認其他scope leases皆閒置。無法讀取／取得鎖就停止；不自行啟停任何服務或程序。DGX無法從本機ps證明DGM或其他主機沒有工作，因此終端確認文字也包含操作者對跨主機閒置的確認；手動繞過平台鎖的啟動不得並行。

## 刪除與紀錄

確認後，在共用admission鎖內重新盤點。日期、輸出、檔案樹metadata或相依關係任一變動，都要求重新檢視，不沿用舊確認。每筆刪除前再核對舊目錄及最新保留版本。metadata digest不讀取整個BAM／VCF內容；核心輸出亦追蹤symlink實際目標的metadata。

只刪已標記DELETE的batch/sample，不刪FASTQ原始資料、不刪三級、不清Nextflow work、不刪batch samplesheet或pipeline_info。rmtree不跟隨symlink；內含掛載點拒絕刪除。實際刪除前先將該sample原子移到同batch下的私人暫存目錄並核對inode，再刪除；這是操作期間的暫存，**不是可復原垃圾桶**。

`actions.jsonl`在移動前就記原路徑／暫存路徑，逐階段flush及fsync。成功記deleted；失敗或中斷保留剩餘暫存路徑並停止，不繼續刪下一筆，不聲稱可以回復已刪內容。

清理後可在平台的檔案清單按一次「更新」，或等原本每日更新；這個工具不直接寫平台索引。批次samplesheet與Nextflow work保留，因此之後刻意重跑舊batch仍可能重新產生那些sample輸出。

## 驗證與界線

合成資料測試涵蓋最新／日期並列、11/5界線、未知建立時間、近期變動、最新核心檔不齊、WES/WGS分開保留、symlink／三級marker引用、最新版本變動時拒刪、實際rmtree只刪核准sample、外部link目標與同batch其他sample保留、Enter不刪、busy lease阻擋、相依root缺失不刪、可貼版與受測Python原碼一致。測試只在本機合成暫存目錄；沒有到DGX搜尋或刪除真實資料。
