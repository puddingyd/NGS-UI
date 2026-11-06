# PGx DRAGEN 比對與人工採用（v9.98；規則 v2）

2026-11-06。畫面與報告沿用固定 21 個 PGx genes；不需重跑 PGx pipeline。未連線院內 server，工具環境、真實路徑權限與院內 Word 仍由操作者驗收。

## 操作

- v9.97（2026-11-06）僅改呈現：警示及人工採用選單在基因標題同排靠右，窄螢幕局部橫向捲動；不一致紅色，其他警示維持黃色。兩邊基因型與確認紀錄收進原本「詳細」，報告預覽基因欄警示亦同排靠右。不改規則或保存結果，不使既有選擇失效。
- v9.98 順序固定為「結果不一致 → 人工確認後採用 → LowGQ 等其他警示」。
- 相同結果不顯示額外註記。順序不同但組成相同、已確認的新舊命名相容，皆視為相同。
- CYP2B6 結果落在多個 DRAGEN 候選之一時只顯示「相位未定」。
- DRAGEN 明列 LowGQ 位點時只加「LowGQ」；真的不同才另加「結果不一致」。不另標 MT-RNR1 TSV 補值來源。
- 不一致的基因在分析區顯示兩邊基因型及「人工確認後採用」下拉選單。選擇平台或某個 DRAGEN 候選即立即保存，不需再按主畫面的儲存。
- 人工採用保存後在背景重新執行固定版本 PharmCAT 的 phenotyper／reporter；基因型、表型、activity score、共同用藥建議與報告一起切換。表型與分數由目前固定版本 PharmCAT 計算，不採用舊 DRAGEN snapshot 的表型覆寫。
- 更新中保留上一份**完整**結果，暫停 PGx 匯出；成功後畫面自動更新。失敗顯示原因，可「重試更新」或「恢復原結果」；不把新基因型配上舊建議，也不把缺失建議說成無需調整用藥。健檢選項停用 PGx 時，其他項目仍可匯出。
- 更新中「恢復原結果」會取消背景工具並恢復全部原平台結果；更新完成後可逐基因切回平台。相同來源及選取組合可重用校驗通過的已完成版本。原始資料不修改。
- 舊版已保存的人工選擇仍保留，但先顯示原平台完整結果並暫停 PGx 匯出，按「重試更新」完成後才能採用。不改比對規則 v2，不要求重跑 pipeline。

## 固定工具準備（由部署者執行一次）

使用 PharmCAT **3.2.0**，資料版本 **2026-02-09-10-28**。官方 JAR SHA256 固定為 `a1ceda4c4f8e2c53ed17ff4d104df82360c04d4b15226d09cb7f854bc585a0e7`。不從最新資料庫臨時抓指引、不在個案工作中下載程式。

預設依序使用 `biotools/pharmcat/pharmcat-3.2.0-all.jar`（host Java 17以上），否則使用既有 `/home/datalake_Intermediate/pipeline/nextflow_containers/pharmcat_3.2.0.sif`，舊 `/home/pipeline/nextflow_containers` 為 fallback。SIF 內 `/pharmcat/pharmcat.jar` 亦須通過相同 SHA256，不僅看檔名。若要讓既有容器提供 Java、JAR 另外固定，可同時設定下列兩個環境變數。

| 設定 | 用途 |
|---|---|
| `NGS_UI_PGX_PHARMCAT_JAR` | 管理者指定官方固定 JAR 的絕對路徑 |
| `NGS_UI_PGX_PHARMCAT_IMAGE` | 管理者指定提供 Java 的既有 Apptainer SIF；指定時 JAR 唯讀掛入 |
| `JAVA_BIN` | 未使用 SIF 時的 Java 可執行檔；沿用平台設定 |

```sh
# 使用與服務相同的帳號及路徑環境；預設只讀，不改樣本或服務
python scripts/check_pgx_pharmcat.py

# 僅當需要準備官方 JAR 時，由部署者明示下載；先驗 SHA 才寫入，不覆寫不同版本
python scripts/check_pgx_pharmcat.py --install
```

工具缺失、Java／容器無法啟動、版本不符、原建議無法零修改重現或選取 allele 無法辨識時，保留完整舊結果並拒絕 PGx 匯出，需修復環境重試或恢復原結果。安裝工具不會啟停服務。Codex 未在院內執行上述準備／部署。

## 重新註解與發布

原 full report 的 `genes` 是 PharmCAT 原生 `GeneReport`。原樣放入 native phenotype snapshot 的 `geneReports`，連同 `matcherMetadata`、`unannotatedGeneCalls` 跑 reporter；先逐一核對完整 gene 資料與所有正規化用藥註解 multiset。不是從簡化 TSV 重建全部 outside calls，也不偽造 matcher JSON，因此保留 CYP2D6 組合／原 activity score、HLA、DPYD／RYR1 多筆元件、缺失位點與相位資訊。

僅把人工修改基因的候選 genotype 寫入 outside-call TSV，phenotype／activity score 留空，由 phenotyper 依固定資料重新註解。把這些新的 GeneReports 合併回原 snapshot，對整份跑 reporter（CPIC／DPWG／FDA，須先通過原來源重播）。包含 CYP2C9 的 phenytoin／fosphenytoin、fluvastatin 等多基因 annotation 整筆重算；其他 gene 的 call、變異、相位、表型、AS 不得改動。reporter 的 `relatedDrugs` 由新建議重建。未辨識 allele／Indeterminate 不套舊表型猜測。

執行紀錄、輸入快照與中間結果保留在 `data/jobs/pgx_reannotation/{job}/`，不使用 Nextflow 工作區，也不重跑 WGS／StellarPGx／HLA。每工作獨立 tmp，工具 process group 有逾時與取消收尾，全站同時只跑一個 PharmCAT。worker 先持 sample/source locks 再取 PGx slot，與三級／PRS 保持相同鎖順序。queued 也擋取消登錄／刪除；刪除操作持相同 sample/source locks 避免競爭。

完成 report 與 manifest 存在 `08_postprocessing/{ID}.pgx_annotations/{job}/`。manifest 綁定原 report hash、全部比對 fingerprint、人工選取組合、引擎 SHA、零改重播驗證與新 report SHA。發布前重驗來源、當前請求與工具，僅以 `{ID}.pgx_review.json.published` 原子切換整份臨床結果；原比較仍以原平台來源計算。衍生目錄、review／append audit 都是 reviewer-owned，正常 pipeline 重跑保留但來源改變即失效。暫存／歷史紀錄含檢體資料，不放 Git。

GET PGx 不啟動工作。背景 job 用獨立 lease 判斷存活，不把被重用的 PID 當成原 worker；服務重新載入可從 job/review 找回狀態，worker 死亡顯示失敗可重試。登入保護 `POST /api/samples/{sid}/pgx/reannotation` 只接受 `action=retry` 或 `restore`，不接受路徑、命令或自填判型。舊的 choice API 同樣只接受當前不一致候選。前端輪詢與遲到回應均核對切案 token。

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

reviewer-owned `{ID}.pgx_review.json` 不在 pipeline managed 清單內，正常重跑保留。選擇與 audit 以獨立 flock + 原子 replace 寫入，記確認者、時間、兩邊結果、品質及來源 fingerprint；另保存請求、背景 job、published revision 與歷史 manifest pointer。主畫面一般 report PUT 不會覆蓋它；回復平台清除請求及 published pointer、保留 audit／歷史結果。

fingerprint 包含比較規則版本、pipeline source、PharmCAT 內容hash／身分／版本／該gene資料及 DRAGEN 內容。來源改變後舊選擇失效，提示重新確認，PGx 匯出暫停；可明確重新採用 DRAGEN 或恢復平台。遺失來源不會靜默把舊選擇當作平台結果匯出。

選擇保存或重新註解未完成時不允許 PGx 匯出。切個案後遲到的回應不套用。Word 建立前後核對包含 published revision 的 fingerprint；來源、選擇或版本在建立期間變更則不下載、不歸檔。未選 PGx 的健檢報告不受這項檢查影響。

## 驗證

使用合成個案驗證順序／dosage／分組、候選解、NAT2版本映射、LowGQ、MT fallback、來源衝突與快照損壞、API權限、執行中鎖、保存／回復／audit、來源更新失效、Word實際下載及匯出中變動拒絕。桌機／390px手機瀏覽器測試選擇、回復、保存失敗與安靜的相同結果。

v9.98 另驗 pending 保留完整基線、後續修改失敗保留已發布完整版本、重試、取消 owned 工具與子程序、死亡 worker、queued 刪除保護、cache hash 與來源改變拒絕發布。使用官方固定 JAR 與合成 CYP2C9／HLA-B／SLCO1B1 實跑重算，驗共同建議改變及未知 allele 拒絕。本機對使用者提供的 26T01086 實跑 detached worker：零改 184 筆註解一致，修改 CYP2C9 後其餘 20 個顯示基因結果維持、共同建議重算、Word 使用新結果、恢復原結果、原始三檔 byte 不變；檢體資料／輸出均不入 Git。未重跑 Phase2 的126案，院內 Apptainer、完整服務部署與正式 Word 仍待操作者驗收。

```sh
PYTHONPATH=backend python -m pytest tests/test_pgx_phase2.py tests/test_pgx_concordance.py tests/test_pgx_tsv.py tests/test_pgx_main_screen.py tests/test_health_report.py tests/test_sample_layout.py tests/test_specialized.py tests/test_prs.py -q
node --test tests/frontend_*.test.cjs
PYTHONPATH=backend node tests/browser_pgx_ui.cjs
# 選用：官方固定 JAR 實跑合成資料（其餘測試不需 Java）
NGS_UI_TEST_PHARMCAT_JAR=/absolute/path/pharmcat-3.2.0-all.jar PYTHONPATH=backend python -m pytest tests/test_pgx_reannotation_engine.py -q
```

瀏覽器測試需可解析 Playwright 及本機 Chrome；可設定 PYTHON、CHROME_EXECUTABLE_PATH、PGX_SCREENSHOT_DIR。使用者附件只在本機唯讀核對，未放入 repo。
