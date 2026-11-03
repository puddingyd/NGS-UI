# 三級工作區整包移至 NVMe — v9.82

此流程由操作者在 `n102968@server` 執行；Codex 不連入伺服器部署。前置條件為已完成方案 A：`ubuntu-vg/ubuntu-lv` 與 ext4 擴充至約 951 GiB，`/scratch` 使用根檔案系統。2026-11-03 實際查詢：根目錄 `df` 936 GiB、可用856 GiB，`nf_work` 35 GiB，來源仍為 `/dev/sda` HDD。

## 資料配置

實際資料：`/scratch/n102968/nf_work/`。以 bind mount 掛回 `/home/n102968/NGS_UI/nf_work/`。程式、既有 session、work-dir pointer、PRS config snapshot 與工作檔仍使用原絕對路徑；不是兩套運算資料，也不建立符號連結。

保留服務現有 `NGS_UI_TERTIARY_NF_WORK_ROOT=/home/n102968/NGS_UI/nf_work`（未設定時預設即此路徑）。PRS未覆寫時自然使用其下 `prs/`；若另設 `NGS_UI_PRS_WORK_ROOT`，須確認它也在這棵樹中。**不要改成 `/scratch/...` 後又期待舊的路徑 fingerprint 不變。** 搬遷工具在服務執行時只讀取並核對相關路徑環境，不輸出其他環境或密碼。

```text
nf_work/                              # 邏輯路徑不變；實際為 NVMe
  contexts/shared-{dragen,nckuh}/
    launch/.nextflow/{cache,history}   # 既有 session／續跑紀錄
    work/                             # 或 work-dir.txt 指定的既有子目錄
    .runtime/
      nextflow/{framework,plugins,assets}
      containers/, apptainer/, conda/
      jobs/{jobID}/
        nextflow.config               # include 原 config，只覆寫執行路徑
        tmp/                          # host / Nextflow / container 暫存
        staging/{sourceID}/00–08       # 完整三級＋post-processing 運算
  prs/{個案ID}/                       # 原有 PRS 工作、runtime、稽核一起搬移
  ...                                 # 所有舊式 work 目錄完整保留
```

`work-dir.txt` 沒有列出時，不據此刪除任何資料；舊式目錄仍整包搬移。舊 repo `.nextflow` 仍可作共享 lineage 的首次 seed 來源。

正式結果仍在 `/home/datalake_Intermediate/pipeline/tertiary_output/`。原始 VCF/BAM、大型 reference、模型、程式與平台 jobs/判讀狀態使用既有路徑。

## 程式行為

- 三級 staging 及 SNV working TSV／overlay／review／索引／CNV rescue 改在上述 scope 下運算；同一 scope 的 lifetime lease 保護整個計算、發布與清理過程。
- 三級 Nextflow/Apptainer/Singularity/Conda 可寫入快取與 TMPDIR/TMP/TEMP/NXF_TEMP 置於工作區；既有 framework/plugins/assets 首次複製至該 scope，避免只因切換路徑重新下載既有 Nextflow 安裝。
- worker source 院內環境腳本後，再套用 runtime 路徑；保留腳本加上的工具 PATH 與 reference bind。衍生 config 保留原參數及 `dgm` profile，追加 work root／tmp bind，將容器 `/tmp`、`/var/tmp` 與環境變數導向 job tmp；`process.scratch` 亦在此處。
- 全批運算成功後，先複製到正式輸出檔案系統的 `.staging/{jobID}/{LIS_ID}/`，逐檔比較 SHA256、保留mtime，將 symlink 實體化成一般檔案，不依賴即將刪除的SSD staging。
- 全批複製完成才開始原有 promotion。SQLite/JSON 的原始資料路徑與 overlay signature 以目的地實際檔案重定位；marker保持相對 raw 路徑。`.rollback` 與 `.staging` 都在正式磁碟，以同檔案系統 rename 替換，layout marker 最後寫入。判讀、metadata、analysis version 沿用既有保留規則。
- 複製損壞、空間不足或取消不發布；後案發布失敗會回復前案。若回復本身失敗，Log明示並保留 `.rollback` 和 publication staging 供復原，不刪掉舊資料唯一副本。
- 正常完成／失敗會清除自己的SSD job tmp/staging；保留Nextflow work/cache、PRS audit及既有共享鎖。部分發布回復失敗的保留目錄須人工核對。

## 操作順序

### 1. 先更新程式並檢查（此時服務維持執行）

```bash
cd /home/n102968/NGS_UI/NGS-UI
git pull --ff-only
sudo python3 scripts/migrate_nf_work_to_scratch.py
```

此指令預設唯讀，確認來源、目的地檔案系統與服務實際三級／PRS／DATA root。若顯示 DATA root 不同，後續兩次指令都加上提示的 `--data-root /實際路徑`。若 PRS 使用工作樹外的目錄，工具會停止，先釐清並遷移該配置，不能把它誤當成已隨 nf_work 搬移。

目的地不得已存在、為 symlink、互相包含，或位於 HDD；工作樹若有其他掛載點亦停止。根目錄若已是掛載點，代表可能已搬遷，不重複搬。確認已安裝 Ubuntu `rsync`（支援 `-aHAXS`）。

### 2. 等所有三級／PRS／Somatic 工作結束，再停止網頁服務

先在 UI 確認沒有 queued/running/cancelling；需要中止的工作從 UI 取消並等結束。

```bash
sudo systemctl stop ngs-ui
sudo python3 scripts/migrate_nf_work_to_scratch.py --apply
```

停服務不等於 detached worker 已停止。工具會再次檢查狀態、Nextflow程序與工作保護鎖；有舊式／手動／無法歸屬程序時停止，不使用 kill 或修改 state.json 繞過。整個搬移期間持有 admission 與 DRAGEN/NCKUH/PRS exclusive locks，父目錄鎖位置不變。

工具依序：

1. 核對容量，至少保留10 GiB可用空間；建立個人scratch父目錄。
2. 用 `rsync -aHAXS --numeric-ids` 複製完整工作樹，保留隱藏檔、symlink本身、hardlink、mtime、ACL、xattr及sparse files，**不使用 `-L` 展開大型輸入/reference連結**。
3. 用 checksum＋dry-run核對內容、額外檔、權限及連結；不一致時不切換來源。
4. 將HDD來源改名為 `nf_work.hdd-backup-{時間}`，建立空掛載點，bind mount SSD目錄回原路徑；驗證samefile。
5. 備份 `/etc/fstab`，新增：

   ```fstab
   /scratch/n102968/nf_work /home/n102968/NGS_UI/nf_work none bind,nofail 0 0
   ```

6. 建立 `ngs-ui.service.d/50-ngs-scratch.conf`，加入 `RequiresMountsFor`、`mountpoint` 及根檔案系統device核對，執行daemon-reload。`nofail`讓工作目錄掛載失敗時作業系統仍可開機，服務本身則必須通過掛載檢查，不能回落HDD。**工具不啟動服務。**

若複製／核對失敗，來源不切換，目的地副本保留待查；不要直接重跑或刪除來源。若掛載失敗，工具恢復HDD原目錄。若掛載成功但永久設定失敗，先勿啟動服務，依Log所列fstab及HDD備份復原或修復。

### 3. 核對掛載，再啟動服務

```bash
findmnt -T /home/n102968/NGS_UI/nf_work
df -hT / /home /home/n102968/NGS_UI/nf_work
sudo systemctl start ngs-ui
sudo systemctl status ngs-ui --no-pager
```

`nf_work` 應顯示 NVMe 的 root LV（可能帶 `[/scratch/n102968/nf_work]` 子路徑），`/home` 仍是 `/dev/sda`。先驗收 DRAGEN、NCKUH各一案及PRS，再恢復批次量；相同輸入／版本／參數應看到cached tasks，實際是否重跑由Nextflow判定。新版runtime設定可能使受影響task失效，不保證每一task都cached。

三級Log應顯示 `[workspace]` 的staging/tmp及 `[publish-copy] ... verified ...`，工作檔實際使用NVMe，完成後正式03/08仍可讀。PRS應維持原個案work路徑，並通過 `CONTAINER_TMP_PASS`；已完成結果與風險快照可讀。

檢查至少一個新執行task的`.command.run`，確認含work root與tmp binds；在院內以真實容器驗證讀寫、TMPDIR及 `/tmp`、`/var/tmp`。本機測試使用真實Nextflow及合成資料／模擬容器，未在院內執行WGS。

### 4. 空間與備份

```bash
df -h /
du -sh /scratch/n102968/nf_work
systemctl status fstrim.timer --no-pager
```

以整個根檔案系統使用量管理：60%提醒、70%整理、80%暫停新大型工作，並估計每批峰值。這是維運規則，沒有自動刪檔或自動終止工作的新增功能。既有「三級Nextflow暫存清理」仍保留active scope及PRS audit，不能對整包work設定按天自動刪除。fstrim.timer未啟用時由管理者確認SSD discard支援後設定。

HDD舊副本在驗收前保留，非另一套工作區；保留期間會佔用HDD約35 GiB。不列入一般分析或續跑，也不修改副本。確認成功後由操作者清理。

## 回復

服務已開始寫SSD後，不能卸載就直接啟動舊HDD副本，否則會遺失搬遷後的工作狀態。先等所有工作結束、停止服務並確認沒有lease；用與搬遷相同的rsync保留選項，將 **SSD目前內容** 同步至選定的HDD回復目錄並checksum核對。需刪除同步多餘檔案時，先檢視dry-run，不對未確認目錄使用delete。

完成後才卸載bind、移除對應fstab項與`50-ngs-scratch.conf`、將核對過的HDD目錄恢復為`nf_work`、daemon-reload並驗證HDD掛載。原始資料、reference、正式結果不因工作區回復而搬動。此節是故障復原指引，不由工具自動執行。

## 本機驗證

164項相關測試通過，涵蓋三級rerun/CNV rescue/logs/batch、Nextflow lease清理與PRS工作。使用Nextflow24.10.4實際解析一般／dgm profile、執行合成workflow並在更換job tmp後以`-resume`取用cache；以模擬container engine核對真實`.command.run`的bind/環境傳遞。另驗證跨磁碟rename限制、全批先複製、checksum錯誤／空間不足／取消不碰live、後案失敗回復前案、回復失敗保留備份，以及複製後overlay/CNV來源signature有效。搬遷工具測試沒有真的掛載或修改`/etc`；院內Linux掛載與真實容器WGS由操作者驗收。
