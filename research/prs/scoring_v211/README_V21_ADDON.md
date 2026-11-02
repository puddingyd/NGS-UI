# PGS004879 integrated v2.1.1｜Orchestration layer

這個目錄中的 `pipeline_v21.py` 與 `v21_raw_resume.py` 由 **v2.1 add-on** 安裝到原 v2 目錄，後續完成 **v2.1.1 batch `Namespace` hotfix**（`expected={r.sample for r in jobs}`），並新增 `tests/test_v21_batch_success.py` 兩項測試。上傳的原始伺服器 ZIP **已包含**這些檔案，不需要再次執行舊 add-on installer。

- 科學計算（`pipeline.py`、`pipeline_v2.py`、`method2_score.py`、vendor raw script）**未更動**；`BASELINE_V2_SHA256SUMS.txt` 十項核對均通過。
- 專案內 `config.json` 已將 `raw_work` 改成 `/home/n102968/2026_INTERVENE_PRS`；新版 `config.example.json` 也對齊此值。
- Server evidence（2026-11-02）：unit 27/27、reuse 3/3 + 6/6、fresh 3/3 + 6/6、completed-resume SHA256 9/9、cross-process lock、post-roundtrip failure recovery、Nextflow cached task probe、**2-sample parallel fresh 2/2 + 4/4 全數 PASS**。
- **未驗證**：Nextflow 已完成部分 tasks 時發生中斷，重新啟動實際重用部分 tasks 的 recovery；5–10+ samples 壓力測試、量產資源／container cache 管理。

目前主手冊 [`V21_SOP.md`](V21_SOP.md)；查實測證據 [`VALIDATION_REPORT_v211_2026-11-02.md`](VALIDATION_REPORT_v211_2026-11-02.md)；後續工作 [`NEXT_STEPS_v211.md`](NEXT_STEPS_v211.md)。

**研究用途／無 CHD outcome-calibrated absolute risk。**
