# PRS 研究批次搬到本機資料碟

更新：2026-11-20，v10.43。此工具只搬既有 `nf_work/prs/_validation_*/study`，不搬整個 nf_work、不啟停服務。由操作者在院內 Linux 主機執行；本機只有合成資料測試，沒有連入院內搬檔。

## 一、為什麼搬，以及保留什麼

操作者回報根磁碟只剩23 GiB，nf_work約831 GiB，其中PRS約736 GiB；本次cohort_v1約581 GiB，排除ancestry約170 GiB，祖源約411 GiB。批次在第504筆途中正常Ctrl-C；迴圈序號不是重新核對後的完成檢體總數，先前試跑也可能包含排序較後的檢體。

這次把完整study（含已完成、未完成、checkpoint、基因型、audit、祖源工作檔及既有報告）複製到本機HDD，核對後bind掛回原路徑，再釋放記錄中的SSD副本。**不修改plan、模型、設定或科學程式，也不重新prepare。**不要直接mv後改`--study`，也不要用symlink取代此流程。

```text
程式繼續使用：
/home/n102968/NGS_UI/nf_work/prs/_validation_20261114/cohort_v1

搬移後的實體資料：
/home/n102968/NGS_UI/prs_validation_work/_validation_20261114/cohort_v1
```

其他三級contexts及其他PRS工作目錄維持原位置。工具只在PRS執行名額及本study加鎖，使用既有PRS cleanup lease，允許其他NCKUH／DRAGEN scope繼續運作；搬移期間其他PRS計算會等待或拒絕啟動。仍須監看根磁碟：既有三級可能繼續消耗空間，不能因PRS停止就認為用量不再增加。

## 二、先只取工具，避免更新正在執行的三級程式

不必git pull或重啟服務。以下只fetch程式物件，再取一個獨立工具到私人的工具目錄，工作目錄中的程式保持原版。工具從目前repo載入既有的Nextflow鎖模組，不載入新的科學流程。

```bash
bash <<'BASH'
set -euo pipefail
cd ~/NGS_UI/NGS-UI
git fetch origin claude/plan-ngs-ui-RQW8J
PRS_TOOL_DIR="$HOME/NGS_UI/.prs_storage_tool"
install -d -m 700 "$PRS_TOOL_DIR"
git show FETCH_HEAD:scripts/migrate_prs_study_storage.py > "$PRS_TOOL_DIR/migrate_prs_study_storage.py"
sudo python3 -B "$PRS_TOOL_DIR/migrate_prs_study_storage.py"
BASH
```

預設為檢查：確認study plan、凍結參考、來源／目標範圍、無symlink／子掛載、無程序存取、目的地確為另一個本機/home檔案系統及可用容量。預留至少100 GiB或study大小10%的餘裕；現有部分副本僅接受相符搬移紀錄。既有目的地沒有紀錄時不覆寫。這不是重新跑PRS，也沒有複製／刪除study；全面逐檔checksum在下一步。

## 三、複製、核對，再掛回原路徑

PRS批次保持停止。建議在tmux執行，維持終端連線；581 GiB的複製與逐檔核對可能需要數小時。請在原repo目錄執行：

```bash
cd ~/NGS_UI/NGS-UI
sudo python3 -B "$HOME/NGS_UI/.prs_storage_tool/migrate_prs_study_storage.py" --apply
```

工具會依序：

1. 持有原study鎖、全站PRS執行鎖與PRS清理lease，禁止相同study寫入與cleanup；不長時間持有全站admission鎖。
2. 以rsync `-aHAXS --numeric-ids --protect-args`複製，保留隱藏檔、權限、mtime、ACL、xattr、hardlink及sparse；不展開symlink。新資料目錄沿用來源uid/gid，非一律留成root擁有。
3. 以`--checksum --dry-run --itemize-changes --delete`核對內容、metadata及多餘檔案。此處delete只在dry-run，沒有刪除。比對非空即停，兩份資料均保留；另核對plan原始hash與凍結模型／參考hash。
4. 同時鎖住來源與目的地的execution.lock，避免切換inode後新工作繞過鎖。將原study改名為同層`cohort_v1.ssd-backup-時間戳`，把HDD副本bind掛回原study路徑。
5. 加入專屬`/etc/fstab`的`bind,nofail`設定並保留fstab備份；只做systemd daemon-reload，不重啟服務。既有父層nf_work掛載不改。底下空掛載點設為000，若日後bind未掛上，研究流程會因權限／缺plan停止，不默默在SSD重新開study。
6. 顯示掛載後磁碟用量；**此時SSD副本仍在，根磁碟尚未釋放約581 GiB。**

控制紀錄在HDD study同層的`.cohort_v1.migration.json`，不在科學study內容內；紀錄phase、來源inode/device、原plan hash與核對後目錄metadata簽章，不把任何病人內容提交Git。

## 四、釋放已核對的SSD副本

看到第三節「HDD掛載與永久設定完成」後，先不要report或恢復分析，接著執行：

```bash
cd ~/NGS_UI/NGS-UI
sudo python3 -B "$HOME/NGS_UI/.prs_storage_tool/migrate_prs_study_storage.py" --release-ssd
```

這一步會再次逐檔核對SSD副本與HDD資料，所以也需要時間。只在bind身分、永久設定、備份名稱／inode、無使用程序／子掛載及目標未變更均通過後，刪除**本次紀錄中那個SSD副本**。HDD study完整保留，原始gVCF、其他PRS目錄及三級contexts不動。不接受通配符刪除、不執行Nextflow clean。

釋放操作會先記錄releasing；途中Ctrl-C可重貼相同`--release-ssd`，在目的地metadata未變且掛載／備份身分仍相符時，另核對剩餘SSD檔案仍與HDD相同，再接續清除；備份新添或變更的檔案會使操作停止。若掛載後先啟動分析或產生新報告，目標簽章改變就拒絕清除；先回傳錯誤，不手改紀錄。

成功顯示`搬移狀態：released`後，檢查：

```bash
findmnt -T /home/n102968/NGS_UI/nf_work/prs/_validation_20261114/cohort_v1
df -hT / /home /home/n102968/NGS_UI/nf_work/prs/_validation_20261114/cohort_v1
```

study應顯示/home資料碟（目前為/dev/sda）的子目錄掛載，根磁碟可用空間應增加約study原佔用量，實際受其他三級寫入及仍開啟檔案影響。工具同時掃描logical與原SSD實體alias的程序引用，只輸出PID、不列命令／環境秘密。其他與study無關的NCKUH／容器程序不應被誤判或終止。

## 五、驗收後沿原命令接續

先產生一份報告，讓原流程讀取並核對已完成checkpoint。成功後才接續計分；下面整段在第一個錯誤即停：

```bash
bash <<'BASH'
set -euo pipefail
cd ~/NGS_UI/NGS-UI
PRS_STUDY=/home/n102968/NGS_UI/nf_work/prs/_validation_20261114/cohort_v1
.venv/bin/python scripts/run_prs_cohort_validation.py report --study "$PRS_STUDY"
.venv/bin/python scripts/run_prs_cohort_validation.py run \
  --study "$PRS_STUDY" --group healthy --all-healthy --stage score
BASH
```

仍使用原study與設定，不重新prepare。已完成且校驗通過者沿用；中斷檢體依原checkpoint重做尚未完成的部分。report與接續核對本身也會讀取大型檔案，不等於重新從gVCF計分。若來源、資源或程式真的變更，原有拒混版規則仍生效，不能靠搬移繞過。

## 六、失敗、後續保存與測試範圍

- 複製／checksum階段中斷：原study維持SSD原樣，有相符紀錄時可重貼`--apply`接續rsync，不覆寫科學plan。
- mount失敗且沒有成功掛上：恢復原study名稱，HDD副本保留。若切換期間主機斷電或被SIGKILL，保留紀錄及兩份資料，先回傳錯誤，不自行mv／umount／改plan。
- bind已成功但fstab／daemon-reload失敗：SSD備份保留，未允許release；以`--apply`完成永久設定，衝突時停止，不覆寫其他管理者設定。
- 若日後需移除研究study，必須先解除其掛載並移除這次fstab專屬列；不可直接刪除/home目標目錄或對掛載點rmtree。此工具不提供研究資料最終刪除功能。
- 搬到HDD不等於完成封存或自動清理。大型祖源參考重複檔仍保留；共享參考、暫存清理與長期datalake證據封存是後續工作。全量祖源維持暫緩，不把每案約15 GiB的祖源工作複製擴到1516份。
- 本輪相關測試66項通過、5項選配工具檢查未執行。合成測試涵蓋實際rsync的隱藏／hardlink／symlink／sparse複製、同大小同mtime內容損壞、額外檔案、鎖衝突／其他scope可執行、mount失敗回復、fstab衝突及失敗保留、精確備份刪除與中斷接續、目標改變拒刪、流程冪等、程序引用與bind alias。沒有在本機執行真實Linux mount、fstab或院內WGS搬移；需以上院內驗收。
