# Nextflow 暫存安全清理（v9.77）

v9.82三級staging、post-processing及工具runtime已移入同一scope的`.runtime/`，既有lifetime lease一併保護計算、全批驗證複製、發布與清理。整包work可bind mount至NVMe，保留原路徑與父目錄guard鎖；不可用symlink替代根目錄。見 [SSD搬遷操作](TERTIARY_SCRATCH_MIGRATION.md)。

二級、三級都採「跳過使用中的範圍，清理其他確認安全的暫存」。清理不終止任何工作。

## 二級：在 DGX2 執行

更新 UI 後，重新產生啟動指令及「顯示 DGX2 清理指令」；先前複製的舊指令不會自動更新。
guard 直接內嵌在指令中，DGX2 只需 Python 3.9+ 標準函式庫與原有環境，不需另行部署 helper。
啟動指令自動寫出 `/tmp/ngs2_<batch>.guard.py`，讓傳給 tmux 的命令維持短小。

- 每個 `NGS_UI_SECONDARY_DGX_WORK_ROOT/<batch>` 是一個清理單位。
- 新版啟動先取得 batch 的 shared lifetime lease，再執行 setup、Nextflow、WES QC、Log 複製。
- tmux 的互動 shell 在 lease 外，完成後即使 pane 未關閉也可清理該批次。
- 清理先顯示候選批次並詢問，再重新檢查當下狀態；遇到 busy lease 就跳過。
- `NGS_UI_SECONDARY_DGX_OUTPUT_ROOT` 內仍以符號連結引用的 work 保留。

## 三級：從三級分析清單清理

- v9.80起PRS重建／raw audit／Nextflow暫存集中於同一work root下的`prs/{個案ID}/`（一般為`~/NGS_UI/nf_work/prs/`）。清理明確保留整個`prs/`及目前自訂PRS root所在單位，即使PRS已完成也保留稽核與接續分析資料；不會把它當閒置DRAGEN／NCKUH task清除。完整刪除個案才依ownership marker清除該案PRS work。PRS worker與直接runner使用同一guard root、`prs` scope的lease；沒有DRAGEN／NCKUH的work-dir pointer不表示PRS尚未初始化。
- DRAGEN、NCKUH 各自以整組共用 lineage 為單位；不在執行中的同一組內挑選已完成 task 刪除。
- `start_job()` 在寫 queued 前取得 lifetime lease，透過檔案描述元交給 detached worker。
  Worker 及其子程序保留 lease，直到 post-processing、promotion、staging 清除完成。
- 原有 `.nextflow_locks/{dragen,nckuh}.lock` 繼續保護 Nextflow LevelDB；清理也檢查舊鎖。
- 執行中、等待中及後處理中的 scope 均保留。`work-dir.txt` 指向的舊目錄與相依 context 一起保留；
  尚在初始化、沒有有效 pointer 時保留所有可能承接的舊暫存。
- 檢查正式及 legacy pipeline output、受保護 work 的符號連結，保留中間目錄與最終目標。
  若輸出依賴 work 連結，即使分析已完成也不刪除；本功能不自動搬移或改寫輸出。
- UI 顯示已清理、保留、失敗的數量與逐項原因。部分刪除失敗不計為成功或釋放容量。

## 併發保護與相容性

鎖放在 work root 同層的 `.<root名>.ngs-ui-locks/`，清理永遠不刪除此目錄。
新建鎖目錄採 work root 的 group 與 mode `2775`；DGX2 指令使用 `umask 0002` 建立共用鎖檔。
執行帳號必須能在 work root 的父目錄建立此目錄，並能共同讀寫鎖檔；權限不足會停止，不能繞過鎖清理。

啟動及清理共用 admission lock：清理的檢查、引用分析、刪除期間不能註冊新的受管理工作。
既有工作持有各自 lifetime lease 繼續跑；清理以 nonblocking exclusive probe 跳過使用中的 scope。
Lease 關閉時不明確 `LOCK_UN`，因此 launcher 異常結束後，已繼承描述元的子程序仍受到保護。

清理用 `ps` 檢查 Nextflow shell／Java 程序，並從 Linux `/proc/<pid>/environ` 核對 root、scope 與 live lease。
舊版、手動啟動、orphan 或其他無法核對的 Nextflow，會使本次所有候選保留。
若程序清單無法讀取則停止；輸出／引用無法完整檢查則全部保留。
服務升級不會替已啟動的舊工作補上 lifetime lease，仍使用工作狀態、舊 cache 鎖及保守程序檢查。

完整的啟動／清理併發保護適用於新版受管理啟動流程。
手動繞過 guard 的 Nextflow 不應與清理同時啟動；無法用一次程序掃描保證阻擋掃描後才出現的任意外部命令。
不依檔案年齡、已完成標記或是否有開啟檔案來推斷可刪除。

清理會使已刪除的暫存失去原本的 `-resume` 效益；跳過的範圍保持原狀。
本次不部署或重啟伺服器，由使用者自行更新 UI 服務。

## 驗證

`PYTHONPATH=backend python -m pytest tests/test_nextflow_work.py tests/test_secondary_analysis.py tests/test_dragen*.py tests/test_secondary_qc_report.py tests/test_secondary_index_schedule.py -q`

測試包含真實跨程序 lease 繼承、啟動／清理互斥、使用中 A 批與閒置 B 批、queued 初始化、舊 cache 鎖、
共用舊 work、多層輸出連結、未知程序、無法讀取引用及部分刪除失敗；使用合成暫存，不清理正式資料。
