# PGx DRAGEN 比對與人工採用（v9.94；規則 v2）

2026-11-05。畫面與報告沿用固定 21 個 PGx genes；既有個案更新服務、重新載入即可，不需重跑 PGx。未連線院內 server，真實路徑權限與院內 Word 仍由操作者驗收。

## 操作

- 相同結果不顯示額外註記。順序不同但組成相同、已確認的新舊命名相容，皆視為相同。
- CYP2B6 結果落在多個 DRAGEN 候選之一時只顯示「相位未定」。
- DRAGEN 明列 LowGQ 位點時只加「LowGQ」；真的不同才另加「結果不一致」。不另標 MT-RNR1 TSV 補值來源。
- 不一致的基因在分析區顯示兩邊基因型及「人工確認後採用」下拉選單。選擇平台或某個 DRAGEN 候選即立即保存，不需再按主畫面的儲存。
- 畫面、健檢 Word 的 allele 表、報告預覽／PDF 都使用已保存選擇。可隨時切回平台；原始檔不修改。
- 採用 DRAGEN 後，只用其明列表型；沒有表型就顯示未指派，不套原 PharmCAT 的舊表型／activity score。
- 現有 PharmCAT 用藥建議是以原基因型產生。採用 DRAGEN 後暫停受影響的原 annotation，連同包含該 gene 的多基因建議，並顯示「用藥建議待重新註解」。這不代表無須調整用藥。此版不重新執行 PharmCAT、不自行生成新指引建議；切回平台即恢復原有建議。

## 來源

平台結果由 `pgx_tsv.load_pgx()` + 共用 `docx_export` presenter 提供，JSON 優先，只有 MT-RNR1 Unknown 可以依既有規則以 TSV 補值。TSV 的 NOTES 作為原 pipeline 紀錄保留，不作新判定唯一依據。

DRAGEN 優先讀已通過 specialized manifest／generation／checksum 校驗的完整本地 source 快照。既有快照損壞時顯示來源問題，不繞過校驗改讀別的 run。

沒有 bundle 的舊案，僅由 `pipeline_source.json.source_vcf_path`、`source_sample_id` 定位：

```text
<VCF所在目錄>/<source>.targeted.json
<run>/other/<source>/germline_seq/<source>.targeted.json
```

兩處同時存在而內容不同時拒絕比較。使用 source ID，不用可能帶 `-dragen` 的平台個案 ID。sampleId 與 build 不符也拒絕比較，不全磁碟掃描。原始 JSON 的 locusAnnotations 與頂層 cyp2d6/cyp2b6 都讀取。

## 判定規則

- `/` 兩側及同一染色體 `+` copy 成分排序後比較，保留拷貝數與染色體分組；`x2` 與同一 allele 重複兩次可對照，bracket combination 不等同多 copy 結構。HGVS 內的 `+` 不拆成 copy。
- 不靠表型相同推定基因型相同。不同命名系統無可靠對照時顯示「無法比對」，不硬判不一致。
- CFTR 的 Reference 與 non-responsive sequence、G6PD 的 B(reference)/B(wildtype) 使用 gene-specific 等價規則。
- NAT2 僅對 DRAGEN `PharmGKB Database Version: Snapshot-2022.01.01` 且新版（2024-03 起）PharmCAT definition 套用限定對照：舊 *13/*13A→*4，*7A/*7B→*7。其他字尾不任意刪除。
- DPYD／RYR1 多筆未定相 source variants 不當作多個完整 diplotype，也不自行推斷 cis/trans。
- 明確 `(heterozygous)/(het)`、`(homozygous)/(hom)` 可展開成 variant/Reference 或 variant/variant；HGVS 後的 star 名稱註解（如 `c.1627A>G (*5)`）不改變同一 HGVS 位點的身分。裸變異不預設為 het。
- G6PD 單一 allele 與相同 allele 的同型合子表示視為相同；不從這項比較推定個案核型，不修改原始 ploidy／allele 顯示，也不把真正異型合子壓成單倍體。
- LowGQ 取自原始 missingVariantSites 的明確標記，不新增任意 GQ 門檻。

NAT2 對照來源：[PharmVar change log v1.7](https://a.storyblok.com/f/70677/x/253c992c9f/change-log_nat2_v1-7.pdf)、[PharmVar GeneFocus](https://pmc.ncbi.nlm.nih.gov/articles/PMC12997514/)。`+` 順序忽略是使用者指定的本平台比對政策，不宣稱兩個 caller 已驗證實體排列相同。原始表示法仍保留。

### 既有 Phase-2 多檢體報告的適用範圍

依使用者提供的 `PHASE2_DIPLOTYPE_CONCORDANCE_REPORT.md` 檢視；附件未放入 repo。該附件 SHA256：`47a51de89199ffac3ddef0aca4ec09c5ebb2ca2dbbc6d96215485b9a67621e3a`。本次重建規則情境作合成回歸，沒有原始126案及報告引用的逐案裁決／gVCF檔，未獨立重算其一致率，也不將報告提出的pipeline修復或批次重產建議當成本次執行指令。

| 報告情況 | 本版處理 |
|---|---|
| CYP2D6 `*36+*10` 順序 | 既有排序規則已處理，仍保留 copy number 及染色體分組 |
| G6PD B/B vs B | 基因限定的表示正規化，相同時不顯示註記；B/Radlowo vs B 仍不同 |
| DPYD 單一變異 het/hom 及 star 註解 | 明確表示可正規化，一致不加註記 |
| DPYD／RYR1 多筆元件 | 從全部 source＋完整 comparison_variants 取變異劑量，吻合多變異未定相只標相位未定。呈現完整清單，不能拿第一筆或 recommendationDiplotypes 當作確定相位 |
| NUDT15 `*1/*3` vs `*1/*2;*3/*6` | 只在下列定義／位點證據齊備時視為相位相容，僅標相位未定 |
| RYR1 c.11266C>G vs Reference | 舊報告指出定義範圍差異，但沒有標明 DRAGEN 軟體／panel 版本，targeted.json 也未列完整陰性 panel。維持差異或無法比對，不把所有 Reference 自動當一致 |
| CYP2D6 真實 CNV／分組不同、UGT1A1 *6、G6PD 真正 het | 不由舊群體報告覆蓋當前結果，保留警示與既有人工選擇 |

NUDT15 限定平台 `CLINPGX` definition `2026-02-09-10-28`、DRAGEN `PharmGKB Database Version: Snapshot-2022.01.01`；平台 `*1/*3` 或 `*3/*6`、DRAGEN 列完整兩解 `*1/*2;*3/*6`；PharmCAT 必須保留 hg38 chr13:48037782 rs746071566 的 `AGGAGTC/AGGAGTCGGAGTC` 與 chr13:48045719 rs116855232 的 `C/T`，兩位點均無 phased／PS／warning，gene 亦未定相。缺資料、其他 allele 長度、不同版本、其他候選解均不套用特例。

[ClinPGx／CPIC 2026更新](https://blog.clinpgx.org/cpic-guideline-for/)確認舊 *2 合併為 *3 子型，以及兩處 het 未定相可有 *1/*3、*3/*6。此二解表型可不同，故內部狀態為 `phase_compatible`，不是功能／表型等價；不因此改寫原始結果或產生新用藥建議。

元件比對只核對變異與劑量 multiset，不能驗證已知 cis/trans；若來源宣稱已定相但元件標籤未保留分組，列無法比對。comparison_variants 保留所有位點（含ref），不使用已篩非ref／截短80筆的畫面清單作證據。畫面／健檢表共用 DPYD／RYR1 的完整元件呈現。完整 diplotype 形式的多個替代解則維持替代解性質，不串成假基因型。

規則由v1改為v2，既有人工採用的來源fingerprint失效，須重新確認或恢復平台；不靜默沿用舊裁決。

## 人工選擇與一致性

新增登入保護 `PUT /api/samples/{sid}/pgx/choice`。傳 gene、source、當前候選 genotype、fingerprint、confirmed=true；後端不接受自填 allele／phenotype，也不接受過期來源。執行中三級／PRS持有sample/source鎖時拒絕修改，避免跨代發布。

reviewer-owned `{ID}.pgx_review.json` 不在 pipeline managed 清單內，正常重跑保留。選擇與 audit 以獨立 flock + 原子 replace 寫入，記確認者、時間、兩邊結果、品質及來源 fingerprint。主畫面一般 report PUT 不會覆蓋它；回復平台清除 active choice、保留 audit。

fingerprint 包含比較規則版本、pipeline source、PharmCAT 內容hash／身分／版本／該gene資料及 DRAGEN 內容。來源改變後舊選擇失效，提示重新確認，PGx 匯出暫停；可明確重新採用 DRAGEN 或恢復平台。遺失來源不會靜默把舊選擇當作平台結果匯出。

選擇保存中不允許前端匯出。切個案後遲到的回應不套用。Word 建立前後核對 fingerprint；來源或選擇在建立期間變更則不下載、不歸檔。未選 PGx 的健檢報告不受這項檢查影響。

## 驗證

使用合成個案驗證順序／dosage／分組、候選解、NAT2版本映射、LowGQ、MT fallback、來源衝突與快照損壞、API權限、執行中鎖、保存／回復／audit、來源更新失效、Word實際下載及匯出中變動拒絕。桌機／390px手機瀏覽器測試選擇、回復、保存失敗與安靜的相同結果。

```sh
PYTHONPATH=backend python -m pytest tests/test_pgx_phase2.py tests/test_pgx_concordance.py tests/test_pgx_tsv.py tests/test_pgx_main_screen.py tests/test_health_report.py tests/test_sample_layout.py tests/test_specialized.py tests/test_prs.py -q
node --test tests/frontend_*.test.cjs
PYTHONPATH=backend node tests/browser_pgx_ui.cjs
```

瀏覽器測試需可解析 Playwright 及本機 Chrome；可設定 PYTHON、CHROME_EXECUTABLE_PATH、PGX_SCREENSHOT_DIR。使用者附件只在本機唯讀核對，未放入 repo。
