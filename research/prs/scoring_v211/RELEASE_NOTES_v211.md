# PGS004879 integrated v2.1.1｜Release notes（2026-11-02）

**性質：validated operational candidate / documentation update。** 此更新未更改固定 Model A/B 或任何已鎖定的科學計分程式；新的操作測試有成功紀錄，但仍保留完整 mid-workflow interrupted-task recovery 的待辦。

## 版本組成

| 元件 | 版本／狀態 |
|---|---|
| `pipeline.py`, `pipeline_v2.py`, `method2_score.py`, vendor | v2.0 golden baseline，10/10 checksum PASS；**不修改** |
| `pipeline_v21.py`, `v21_raw_resume.py` | v2.1 orchestration，sample lock、fingerprint、batch、post-roundtrip resume |
| `pipeline_v21.py` batch hotfix | v2.1.1：`expected={r['sample'] for r in jobs}` → `expected={r.sample for r in jobs}` |
| `tests/test_v21_batch_success.py` | v2.1.1 新增成功／失敗 batch 路徑測試（+2，總數 27） |
| `config.json`（伺服器既有） | `raw_work` 修正為 project root；實測 command line 正確 |
| 本版更新內容 | 主手冊、Model 2 使用指引、v1/v2 歷史文件標註、驗證紀錄、next steps、example config、版本／checksum metadata |

## 確認的伺服器實測

1. 27/27 unit tests，golden baseline 10/10，Model 2 lock PASS。
2. `method2_v21_reuse_validation_02`：三案 reuse 3/3、Model A/B 6/6 PASS。
3. `method2_v211_fresh_validation_01`：三案 fresh 3/3、Model A/B 6/6 PASS。
4. 上述 fresh run 再次 `--resume`：三案均 `V21_ALREADY_COMPLETE`，9 個 `method2` 核心檔案 SHA256 不變。
5. OS 層跨程序 sample lock：相同 sample blocked、不同 sample allowed、釋放後再次取得 PASS。
6. `method2_v211_recovery_test_01`：第一次 Nextflow 啟動前 fault injection；保留 609,346-locus / 22-chr upstream；`--resume` 後 A/B 2/2 PASS，上游 checksum 不變。
7. 同一 recovery run completed Nextflow `-resume` probe：cached tasks 出現（最多一個 process 顯示 22 cached），3 項核心輸出 SHA256 不變。
8. **`method2_v211_parallel_fresh_01`**：兩個不同樣本同一 UTC 秒啟動、工作時間重疊；兩案 raw roundtrip／reconciliation／Method 2 全 PASS，A/B 四組 regression **4/4 PASS**。

以上為使用者提供的伺服器終端紀錄；此 ZIP 並無內含 raw gVCF、完整 run 結果或 Nextflow task logs。詳 [`VALIDATION_REPORT_v211_2026-11-02.md`](VALIDATION_REPORT_v211_2026-11-02.md)。

## 未涵蓋／release gate

- **Mid-workflow task interruption → `nextflow -resume` 實際重用先前已完成的一部分 tasks：尚未實測。** 本次 start-failure recovery 不能替代它。
- `--jobs 2` correctness 已驗證，但 `--jobs 3..8`、大量 WGS throughput、RAM/CPU/I/O 峰值、磁碟成長、Nextflow image cache 可離線運作尚未驗證。
- 台灣族群 EAS/ancestry calibration、INTERVENE 年齡／性別 CHD cumulative incidence/absolute-risk mapping、CHD outcome validation 尚未完成；目前研究 PRS percentile 不可直接視為臨床 absolute risk。
- 若 `pipeline_v21.py` 後續有任何程式修改，fingerprint 會變；不要用新程式直接 resume 舊 run，應單獨升版與回歸驗證。

## 發行一致性

`SHA256SUMS.txt` 為本次文件更新後的快照清單；`SHA256SUMS_v2_ORIGINAL.txt` 保留 v2.0 原清單以供歷史查核；`BASELINE_V2_SHA256SUMS.txt` 仍固定十個 v2 核心 hash。bytecode `__pycache__` 不屬於發行必要檔案，已自更新包排除。
