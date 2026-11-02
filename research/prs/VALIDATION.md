# 本次驗證紀錄 — 2026-11-02

本機檢查；未登入院內主機，未重新執行 WGS，未部署。

| 檢查 | 結果 | 證據範圍 |
|---|---|---|
| 上傳 ZIP 發行清單 | 45/45 checksum 一致 | 原附件完整性 |
| 選取匯入程式 | IMPORT_PROVENANCE 全部一致 | 匯入後未改寫既有核心 |
| 原版 unit tests | 27/27 通過 | 離線與 mock 測試，不是 fresh WGS |
| 風險模組測試 | 12/12 通過 | 常數 hazard 解析解、條件機率、離散公式、分組、資料拒絕、CLI、來源 hash |
| CLI validate | PASS / scientific_validation=NOT_VERIFIED | schema 檢查，不是流行病學驗證 |
| Published CHD curves reproduction | **224/224 通過** | 原 EstBB HR＋Estonia baseline／五年離散算法；最大誤差0.004995100317904511百分點，符合原表兩位小數精度 |
| API／快照／拒絕條件 | 17/17 通過 | 身分與來源綁定、認證、快照、參數 checksum、QC／祖源、病史、年齡與併發更新 |
| 既有前端 Node 回歸 | 93/93 通過 | 非院內部署驗收 |
| 獨立瀏覽器 | 通過 | 桌面／390px手機、兩曲線、API數值、改輸入清舊值、>70歲、既有CHD拒絕、切個案晚到回應、無JS錯誤；合成個案 |
| Taiwan calibration / outcome validation | 尚未完成 | GBD2019 成人研究包已建立；院內追蹤資料可用性未確定 |
| CI / uncertainty propagation | 尚未實作 | 結果明確標 null / NOT_ESTIMATED |
| 原 pipeline 中途 task recovery／負載 | 未在本次執行 | 待院內 isolated run |

原版27項測試在 macOS 預設暫存路徑下初次為26/27：`/var` 與 resolve 後 `/private/var` 的路徑字串不同。以實體 `TMPDIR=/private/tmp` 重跑27/27通過；未為此變更原版核心或測試。固定算法的兩種數值模式有意分開，不能以合成測試宣稱重現發表結果。

資料尋找稽核：31個原檔已保存；台灣四類 Rate 160列（含 All ages），鍵值、上下界與競爭死亡點估計檢查通過。Data 4／7 各擷取896筆 CHD 結果。TWB 的 OR 欄位與 BBJ positional schema 標為待審查，未進引擎。詳見 DATA_SOURCES.md 及 SOURCE_MANIFEST_20261102.json。以上31檔為先前來源稽核；後續新增 EstBB 原始 HR，完成重現、成人參數包及平台實作，見操作文件。

本版 UI 採 `competing_exact_v1` 支援任意評估年齡；發表曲線驗算採 `intervene_discrete_v1`。同一組台灣參數在兩性七組、20–75歲每5歲起點的10年／至80歲測試網格，最大差異1.4436119070748987百分點；不能宣稱兩算法逐值一致。原研究 CI 未重現，也未完成歐洲 HR→Model A/EAS735→台灣 outcome 的外推驗證。

參數 release SHA256：`322d4e68d49a2047d7326172a609b3898b8754f658f800f50806f4e2276a0a27`。完整224筆對照與來源 hash 位於另行提供參數包中的 `reproduction.json`；API僅接受固定 checksum。原27項計分測試及ZIP檢查為本次功能前階段紀錄，原版計分檔本次未改寫。
