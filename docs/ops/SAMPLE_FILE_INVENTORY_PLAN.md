# 檔案清單、三級資料瀏覽與共用每日索引

v9.96，2026-11-06。已完成本機實作與合成驗證；院內部署、全量檔名與實際下載由操作者驗收。Codex 不 SSH、不重啟服務、不安裝 server timer。

## 使用方式

主畫面登入後，右上角排列為「三級分析 → 檔案清單 → 上傳個案清單」。清單支援 sample ID 搜尋、每頁50筆、上一頁／下一頁及手動更新。

六欄固定為 sample ID、fastq、bam、SNV/indel VCF、其他檔案、三級資料。檔案欄只列含副檔名的檔名，點擊下載；旁邊沿用 SNV 卡片的 Copy 圖示，複製 server 完整路徑。沒有提供的檔案留白，不輸出缺檔或不適用字樣。

- FASTQ：每個 sample 先顯示兩檔，超過時以 `more` 展開、`less` 收合。保留全部 run／lane，依資料夾及檔名排序，使開頭兩檔維持同資料夾的 R1／R2 順序。lane 數不寫死；11/5使用者貼回的26T01345實際為L001–L008、16檔。
- BAM：主 BAM、`.bam.bai` 或 `.bai`；排除 repeats BAM 及其 index。
- SNV/indel VCF：原始 DRAGEN hard-filtered 或 in-house ensemble fixed VCF。
- 其他檔案：只列 CNV VCF、SV VCF、合併 CNV/SV VCF、DRAGEN gVCF 與 targeted JSON。In-house 不列 gVCF／合併 CNV/SV／targeted JSON；不列 `.tbi`／`.csi`、CNVkit、mito、ploidy、三級衍生 TSV 等其他項目。
- 三級資料：沿用 SNV ClinVar/LitVar 外部連結圖示，另開資料夾瀏覽頁。多個三級個案目錄各有一個圖示，title 顯示目錄名。

同 source ID 的不同 run／batch 都保留，檔名相同者可由 hover title 或複製路徑區分。這是檔案候選彙整，不能據此宣稱某 BAM/VCF 必然來自哪一組 FASTQ；分析仍沿用既有來源驗證規則。

## 搜尋來源

`S` 為原始 sequencing ID，`C` 為 UI 個案 ID。所有 root 沿用執行中服務的既有環境覆寫。

| 類型 | 預設搜尋位置 |
|---|---|
| WES FASTQ | `/home/datalake_Raw/NextSeq2000/`，遞迴辨識 `_S數字_R1/2_001.fastq.gz` 等配對檔名 |
| 重分析 FASTQ | `/home/datalake_Raw/Other/Reanalysis/`；兼容 `/datalake_Raw/Other/Reanalysis/`。`{source}.R1/2.clean.fastq.gz` 採父資料夾為 sample，保留原始 source 別名 |
| WGS FASTQ | `/home/datalake_Raw/Novaseq/{run}/fastq.gz/{S}_S數字_L數字_R1/2_001.fastq.gz` |
| DRAGEN BAM | `/home/datalake_Raw/Novaseq/{run}/bam/{S}.bam`；另查 `/home/datalake_Intermediate/n102968/` |
| In-house BAM | `/home/datalake_Intermediate/pipeline/nextflow_output/{batch}/{S}/02_alignment/{S}.aligned.sorted.bam` |
| DRAGEN VCF | 上述 DRAGEN roots 的 `{run}/vcf.gz/{S}.hard-filtered.vcf.gz`、`.cnv.vcf.gz`、`.sv.vcf.gz`、`.cnv_sv.vcf.gz` |
| In-house VCF | `nextflow_output/**/{S}/04_snv_indel/{S}.ensemble.fixed.vcf.gz`；`05_cnv_sv/{S}.gcnv.vcf.gz`、`{S}.delly.vcf.gz` |
| DRAGEN gVCF | 同DRAGEN roots 下原有PRS路徑：`vcf.gz/` 或 `gvcf.gz/` 的 `{S}.hard-filtered.gvcf.gz`／`{S}.gvcf.gz`，以及 `other/{S}/germline_seq/{S}.hard-filtered.gvcf.gz` |
| targeted JSON | 同DRAGEN roots 下 `{S}.targeted.json`，包含 `other/{S}/germline_seq/` |
| 三級資料 | `NGS_UI_TERTIARY_ROOT`，並保留 legacy pipeline／UI root；依精確 source marker 連回 sample |

索引的 FASTQ 名稱解析另外兼容 `_R1/2`、`.R1/2`、`.fq.gz`、有無 `_S數字`／lane 及 `_merged`。**二級分析的可執行輸入仍由原 parser 決定**：WGS維持lane成對規則、忽略merged；不完整R1/R2只出現在檔案清單，不能因此變成可執行samplesheet。沒有符合命名的FASTQ不猜sample；數量與前30個範例存於索引meta供院內補驗收。

`n102968` 的 PGx GetRM 驗證 FASTQ 不建立臨床清單 sample，因為 FASTQ 僅從二級設定的 raw roots 收集。11/5貼回的140筆壞FASTQ link均位於該驗證目錄，正式掃描不將它們當臨床缺檔。其他相關檔案的權限／來源錯誤仍會阻止新generation發布。

只有 `/datalake_Raw/Other/Reanalysis` 這個已驗證不存在的替代路徑允許首次跳過；若先前已成功掃描，後來消失仍報錯。必要root離線／權限不足會保留舊索引。目錄symlink不遞迴跟隨，需將核對過的實際目標設定為對應root；跳過的目錄link與未辨識檔名均留meta。自訂DRAGEN gVCF若在既有roots之外，需把其資料根加入既有DRAGEN roots後重新更新。檔案清單不修改PRS候選模板或分析設定。

## 共用索引與更新

`backend/app/services/file_index.py` 是唯一正式掃描協調器，資料位於 `NGS_UI_DATA_ROOT/sample_file_index.sqlite`。

1. 正規化既有 FASTQ、VCF、BAM roots，重疊的實體目錄只走訪一次，將同一份檔案發現結果餵給既有二／三級 parser。
2. SQLite 存同一generation的sample清單、可下載檔案registry、二級FASTQ視圖及三級VCF視圖。檔案清單以FASTQ建sample；三級視圖保留沒有FASTQ的VCF候選。
3. 共用跨程序flock，鎖內重查generation、請求時間、root設定及過期時間；重疊手動更新／冷啟動／背景更新合併。同一每日02:00窗口的舊worker也去重，不因舊timer暫未遷移而再掃一次。
4. 先建立完整新SQLite，commit／fsync後原子replace；更新期間讀舊完整generation。掃描失敗不發布部分資料、不抹除既有索引。root設定變動不等24小時，下一次查詢／更新即視為過期。
5. 二級及三級 `load_index()` 優先讀共用視圖；尚無共用DB時可先讀舊JSON並排背景建立。舊JSON不再是新generation的寫入來源。
6. 三級資料夾只枚舉sample目錄與小型source marker，prefixed exact name優先、legacy exact name fallback。任意子目錄內容由資料夾頁即時讀取一層，不做每日深掃。

既有 `/api/secondary/index/refresh`、`/api/dragen/index/refresh` 與新 `/api/files/refresh` 都更新同一generation。舊 `secondary_index_update`／`tertiary_index_update` worker保留相容入口，皆轉至共用worker；舊installer亦轉至新installer。

### Server 排程遷移（操作者執行）

先更新程式並重啟web，再執行一次安裝腳本：

```bash
cd ~/NGS_UI/NGS-UI
git pull
sudo systemctl restart ngs-ui
bash scripts/install_file_index_timer.sh
```

腳本安裝 `ngs-ui-file-index-update.service/.timer`，先成功核對／建立當日完整generation，啟用每日 **02:00 Asia/Taipei** 的新timer，再停用兩個舊timer。若首次更新失敗，腳本停止，舊timer不會被關掉。腳本只啟動索引更新unit，不重啟web。

新unit的 `scripts/update_file_index.py` 讀執行中 `ngs-ui.service` 的MainPID及 `/proc/{PID}/environ`，只擷取程式列出的路徑／TTL設定，不輸出或複製其他環境值。這讓每日worker與web使用相同roots；無法讀取、服務未執行時明確失敗，不偷偷使用另一組預設。以相同帳號執行，仍需在院內確認 `/proc` 讀取政策。

檢查：

```bash
systemctl list-timers --all --no-pager 'ngs-ui-*index*'
journalctl -u ngs-ui-file-index-update.service -n 30 --no-pager
```

ClinVar／LitVar／HPO等資料庫更新排程維持原本頻率，不混入檔案索引。

## 下載與資料夾頁

新API均沿用平台登入保護：

- `GET /api/files`：sample ID搜尋及分頁，只讀已發布索引；冷啟動才建立，過期背景更新。
- `POST /api/files/refresh`：強制共用更新。
- `GET/HEAD /api/files/download/{id}`：只接受索引registry中的檔案ID，下載時重新核對配置root。
- `GET /api/files/browse?root=...&path=...`：三級root內逐層瀏覽，含路徑導覽及上一層；上層到三級root為止。
- `GET/HEAD /api/files/tertiary-download?root=...&path=...`：下載三級root內檔案。

下載逐元件以no-follow開啟，確認regular file後串流，不將整個BAM/FASTQ載進記憶體；支援HEAD、單段Range及suffix range。越界／外部symlink／特殊檔案不讀取；關閉或中斷時關閉descriptor。資料夾頁以HTML escape處理檔名，任意檔案都以attachment下載，不直接執行其中HTML；不提供上傳、改名或刪除。dot目錄（包含`.staging`／`.rollback`）不開放瀏覽。院內的特殊symlink layout需個別驗收。

## 驗證

- Python合成檔案樹涵蓋8lane16檔、BAI兩種命名、限定其他類型、空欄、in-house排除、raw-only／單邊配對、多run、重分析source、VCF-only三級相容、root變更、重疊／舊timer去重、失敗保留generation。
- API驗證登入、搜尋／分頁、下載內容、HEAD／Range、刪除後404、越界／外部symlink、資料夾路徑導覽與惡意檔名escaping。
- Chrome合成前端驗證按鈕位置、六欄、16檔展開／收合、Copy、下載、三級連結、空格、搜尋／分頁／更新、Escape，以及桌機／手機版面。手機可在表格內橫向捲動，不擴張整個頁面。
- 現有二／三級與前端回歸一併執行；實際院內全量命名、各root權限及大型檔下載仍待部署驗收，不能以合成測試宣稱已完成。

## 歷史唯讀探測工具

`scripts/probe_sample_files.py`／`scripts/paste_probe_sample_files.sh` 是11/5前置盤點版本，為保留可重現性不改寫其歷史行為。它列的類型比正式表格多，且可選讀PRS templates；不維護正式索引，也不代表目前六欄篩選規則。已取得的 `/tmp/ngs-file-probe.2hnBhH/inventory.json` 在server上；附件／病患明細不入Git。
