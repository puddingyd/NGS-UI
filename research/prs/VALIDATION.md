# 本次驗證紀錄 — 2026-11-02

本機檢查；未登入院內主機，未重新執行 WGS，未部署。

| 檢查 | 結果 | 證據範圍 |
|---|---|---|
| 上傳 ZIP 發行清單 | 45/45 checksum 一致 | 原附件完整性 |
| 選取匯入程式 | IMPORT_PROVENANCE 全部一致 | 匯入後未改寫既有核心 |
| 原版 unit tests | 27/27 通過 | 離線與 mock 測試，不是 fresh WGS |
| 新增風險模組測試 | 11/11 通過 | 常數 hazard 解析解、條件機率、離散公式、分組、資料拒絕、CLI、來源 hash |
| CLI validate | PASS / scientific_validation=NOT_VERIFIED | schema 檢查，不是流行病學驗證 |
| Published CHD curves reproduction | 尚未完成 | 已取得 baseline 與 expected outputs；缺原版歐洲七分組 HR 及完整數值驗算 |
| Taiwan calibration / outcome validation | 尚未完成 | GBD2019 baseline 已取得但未建包；院內追蹤資料可用性未確定 |
| CI / uncertainty propagation | 尚未實作 | 結果明確標 null / NOT_ESTIMATED |
| 原 pipeline 中途 task recovery／負載 | 未在本次執行 | 待院內 isolated run |

原版27項測試在 macOS 預設暫存路徑下初次為26/27：`/var` 與 resolve 後 `/private/var` 的路徑字串不同。以實體 `TMPDIR=/private/tmp` 重跑27/27通過；未為此變更原版核心或測試。固定算法的兩種數值模式有意分開，不能以合成測試宣稱重現發表結果。

資料尋找稽核：31個原檔已保存；台灣四類 Rate 160列（含 All ages），鍵值、上下界與競爭死亡點估計檢查通過。Data 4／7 各擷取896筆 CHD 結果。TWB 的 OR 欄位與 BBJ positional schema 標為待審查，未進引擎。詳見 DATA_SOURCES.md 及 SOURCE_MANIFEST_20261102.json。本次僅更新來源與文件，未改動數值程式。
