# NGS 分析平台 (NGS-UI)

成大醫院基因醫學部的 NGS 三級分析判讀工具。次級 pipeline（Nextflow，跑在另一台 compute cluster）產出 per-sample 的註解 TSV，本平台讓 reviewer 載入個案、檢視 SNV/Indel + CNV/SV + Mitochondria 變異、標記 causative / candidate / other、整理 ACMG SF / proactive / carrier secondary findings、撰寫判讀意見，並匯出診斷報告 (docx)。另附一個獨立的「臨床表徵輸入 (HPO / gene panel)」工具掛在 `/phenotype/`。

- 後端：FastAPI + uvicorn（Python 3.10+）
- 前端：原生 HTML/CSS/JS，**無 build step**（直接 serve `frontend/`）
- 背景工作：Redis + RQ（Exomiser / LIRICAL 重跑）
- 帳號：SQLite (`data/users.db`) + bcrypt
- 部署：內網 `192.168.84.91:8765`，systemd unit `ngs-ui`

> 開發者 / 接手者請另看 `CLAUDE.md`（架構、資料流、各模組細節、踩雷紀錄）與 `docs/`。

---

## 1. 目錄佈局

所有路徑都從 `NGS_UI_HOME` 推導（每個子路徑也都可以用各自的環境變數覆寫）：

```
NGS_UI/                    ← NGS_UI_HOME
├── NGS-UI/                ← 這個 git checkout（REPO_ROOT）
├── biotools/              ← Exomiser / LIRICAL CLI + data
├── vcf/                   ← per-sample VCF
├── tertiary_output/       ← per-sample TSV + sidecar（不進 git）
│   ├── _index.json        ← 個案清單快取（在 tertiary_output 旁，往上一層）
│   └── {LIS_ID}/          ← snv_indel.annotated.tsv, snv_indel.review.tsv,
│                             snv_gene_index.sqlite,
│                             cnv.annotated.tsv,
│                             sv.annotated.tsv, mito.annotated.tsv,
│                             sample_metadata.json, case_summary.json,
│                             qc_summary.json,
│                             roh_summary.json, analyses/{ver}/...
├── patient_phenotype/     ← {LIS_ID}_{MRN}_phenotype.txt（自動帶入 HPO）
├── patient_list/          ← 上傳的「未完成報告清單」xlsx + 衍生 roster.json
├── phenotype_data/        ← git-tracked fixed/custom panels; large HPO refs live in NGS_UI_HOME
├── ngs_panel_deadzone/    ← expanded reportable gene list、HGNC alias map、dead-zone tables
├── OMIM/OMIM.xlsx         ← OMIM 疾病註解表（缺檔則 Disease 欄留空）
└── data/                  ← server runtime state（users.db, jobs/, ...）
```

開發 checkout 不需要這整棵樹：`NGS_UI_HOME` 未設時，只有典型正式部署路徑 `NGS_UI/NGS-UI` 會把 parent 當作資料根；一般 standalone checkout（例如桌面上的 `NGS-UI`）會 fallback 成 repo 自己，所有路徑都落在 repo 內。

---

## 2. 安裝與啟動

```bash
# 1. 取得程式
git clone <repo> NGS_UI/NGS-UI && cd NGS_UI/NGS-UI

# 2. Python 套件
python3 -m pip install -r backend/requirements.txt

# 3. Redis（背景工作佇列需要；沒有的話 Exomiser/LIRICAL 重跑會無法 enqueue）
#    sudo apt install redis-server && sudo systemctl enable --now redis

# 4. （正式環境）把 patient data / runtime state 移出 repo
#    參考 scripts/migrate_layout.sh。大型 phenotype reference 仍放在
#    NGS_UI_HOME/phenotype_data；固定 panel data 則保留在 repo 內隨 git 更新：
#      mv NGS_UI/NGS-UI/tertiary_output NGS_UI/tertiary_output
#      mv NGS_UI/NGS-UI/data            NGS_UI/data

# 5. 建第一個帳號
PYTHONPATH=backend python3 -m app create-user <username>     # 互動輸入密碼
PYTHONPATH=backend python3 -m app list-users

# 6. 啟動
PYTHONPATH=backend NGS_UI_HOME=/path/to/NGS_UI \
  python3 -m uvicorn app.main:app --host 0.0.0.0 --port 8765
```

正式環境用 systemd（`scripts/migrate_layout.sh` 會幫忙產生 unit）：

```ini
[Service]
WorkingDirectory=/path/to/NGS_UI/NGS-UI
Environment=PYTHONPATH=/path/to/NGS_UI/NGS-UI/backend
Environment=NGS_UI_HOME=/path/to/NGS_UI
ExecStart=/usr/bin/env python3 -m uvicorn app.main:app --host 0.0.0.0 --port 8765
Restart=on-failure
```

```
sudo systemctl daemon-reload && sudo systemctl restart ngs-ui
```

背景 worker（Exomiser / LIRICAL 重跑）另外跑：

```bash
PYTHONPATH=backend NGS_UI_HOME=/path/to/NGS_UI python3 -m app.workers.run   # 視實作而定，見 backend/app/workers/
```

---

## 3. 主要環境變數

| 變數 | 預設 | 用途 |
|---|---|---|
| `NGS_UI_HOME` | repo 的上層（找不到則 repo 自己） | 整棵資料樹的根 |
| `TERTIARY_OUTPUT_ROOT` | `$NGS_UI_HOME/tertiary_output` | per-sample TSV |
| `NGS_UI_DATA_ROOT` | `$NGS_UI_HOME/data` | users.db, jobs/ |
| `NGS_UI_VCF_DIR` | `$NGS_UI_HOME/vcf` | per-sample VCF |
| `NGS_UI_PHENOTYPE_DIR` | `$NGS_UI_HOME/patient_phenotype` | `{LIS}_{MRN}_phenotype.txt` |
| `NGS_UI_PATIENT_LIST_DIR` | `$NGS_UI_HOME/patient_list` | 上傳清單 + roster.json |
| `NGS_UI_PHENO_DATA_DIR` | `$NGS_UI_HOME/phenotype_data` | hp.obo / phenotype_to_genes 等大型 HPO reference |
| `NGS_UI_GENE_PANELS_DIR` | `REPO_ROOT/phenotype_data/gene_panels` | git-tracked fixed panel gene lists |
| `NGS_UI_FIXED_PANELS_DIR` | `REPO_ROOT/phenotype_data/fixed_panels` | git-tracked fixed panel UI index |
| `NGS_UI_CUSTOM_GENE_PANELS_DIR` | `REPO_ROOT/phenotype_data/custom_panels` | git-tracked custom panel gene lists |
| `NGS_UI_PANEL_DEADZONE_DIR` | `REPO_ROOT/ngs_panel_deadzone` | expanded reportable gene list、HGNC alias map、WES/WGS dead-zone tables |
| `NGS_UI_OMIM_XLSX` | `$NGS_UI_HOME/OMIM/OMIM.xlsx` | OMIM 疾病註解（缺檔 = 停用） |
| `NGS_UI_BIOTOOLS_DIR` | `$NGS_UI_HOME/biotools` | Exomiser / LIRICAL |
| `NGS_UI_INHOUSE_BAM_ROOT` | `/home/datalake_Intermediate/pipeline/nextflow_output` | IGV 搜尋 in-house / Nextflow BAM 的根目錄；可用 `:` 分隔多個 root；舊 `NGS_UI_BAM_ROOT` 仍作為 fallback |
| `NGS_UI_DRAGEN_BAM_ROOT` | `/home/datalake_Raw/Novaseq` | IGV 搜尋 DRAGEN raw BAM 的根目錄；尋找 `<run>/bam/{sample}.bam`，排除 `{sample}.repeats.bam` |
| `NGS_UI_IGV_REF_DIR` | `/home/pipeline/reference/hg38` | IGV 本機 hg38 FASTA + `.fai` |
| `EXOMISER_HOME` / `LIRICAL_HOME` / `EXOMISER_DATA_HG38` ... | `biotools/...` | 工具與 data 路徑 |
| `JAVA_BIN` / `JAVA_OPTS` | `java` / `-Xms4g -Xmx16g` | 跑 Exomiser/LIRICAL 用 |
| `REDIS_URL` | `redis://127.0.0.1:6379/0` | RQ 佇列 |
| `NGS_UI_EMR_CLIENT_ID` | `""`（空 = 停用所有 EMR 路徑） | NCKU 內網 HIS / APIM |

固定 WES-I / WES-II / WGS panel 與 custom panel 檔現在都保留在 git 內，server `git pull` 後會直接更新。

---

## 4. Reviewer 操作流程

1. **登入** — 右上角登入（帳號由管理者用 `create-user` 建立）。

### 載入效能

- 首頁登入後的 `/api/samples` 只載入搜尋用的輕量樣本索引，不同步計算個案清單摘要。
- 「個案清單」開啟時才呼叫 `/api/samples/case-summary` 載入 causative / disease / other variant 摘要。
- 個案摘要會寫入每個 sample 目錄的 `case_summary.json`，依 metadata、TSV/index 與 OMIM 簽章自動失效；已標記 SNV 優先用 `snv_gene_index.sqlite` 依 variant id 查找，CNV/SV 只讀目標 id，避免大型 DRAGEN TSV 每次重掃。
2. **載入新個案** — 點「載入新個案」：
   - LIS_ID 下拉會列出 pipeline 已丟進 `tertiary_output/` 但尚未登錄的目錄；
   - 若先用「上傳個案清單」匯入過「未完成報告清單」xlsx，MRN / 姓名 / Test type 會自動帶入（來自 `patient_list/roster.json`）；
   - 若來源為 DRAGEN，或 in-house 來源 VCF 大於 100 MB，Test type 會預設為 `WGS`；送出前仍可手動改回 `WES`；
   - HPO / gene panel 可在這裡選；gene panel 與主畫面同樣使用 `WES-I / WES-II / WGS / Other panel` tabs，預設展開 `Other panel`，固定 panel chip 與搜尋下拉都會顯示基因數量；HPO / panel 下拉可用上下鍵選取並以 Enter 加入，避免 Enter 誤送出載入個案；若存在 `patient_phenotype/{LIS}_{MRN}_phenotype.txt` 會自動讀入；
   - 登錄新個案不再同步掃完整 TSV 產生 `vcf_from_tsv.vcf.gz`；若 VCF 尚不存在，Exomiser/LIRICAL 背景 job 開始前會自動建立或刷新。
   - HPO/panel 的 in-panel 狀態來自 `pheno_score.tsv` 動態補值，不再寫回大型 `snv_indel.annotated.tsv` 的 `IN_PANEL` 欄。
   - 勾「登錄後開始分析」會順便把 Exomiser/LIRICAL 排入佇列。
   - 旁邊的「個案清單」可查看已載入個案、依 `WES` / `WGS` 篩選並全文搜尋；表格會摘要 causative / other SNV、CNV、SV，已勾選 OMIM disease、CNV/SV reviewer 輸入的 Disease、主畫面 comment、簽收與載入時間。多個 variant / disease 會逐行顯示，長 HGVS 可自動折行。也可刪除 `NGS_UI_HOME/tertiary_output/{LIS_ID}/`；刪除時可選擇同步刪除或保留 `/home/pipeline/tertiary_output/{LIS_ID}/`，按取消則完全不刪除。
   - 主畫面搜尋框上方有可複選的 `WES` / `WGS` 圓形 filter，取消勾選後對應 test type 不出現在搜尋清單。
3. **看變異卡片** — 個案載入後先顯示 SNV/Indel（分段載入），CNV/SV 與 Mitochondria 在背景載完後補上：
   - 平台剛開啟讀取索引、個案核心資料載入與新個案登錄期間都會顯示不可誤關閉的「資料載入中」遮罩，避免重複點擊。
   - SNV/Indel tier：`1A / 1B / 1C / 2 / 3`（互斥）
   - Patient phenotype 與 Comment 之間的 Dead zone 卡片會依目前 HPO + panel gene set 顯示 cohort-level dead exons；WES 用 20X，WGS（含 in-house / DRAGEN）用 DRAGEN 15X。主畫面依臨床門檻下的 CDS dead percentage 由高到低排序並顯示比例，預設只顯示前 10 列，可用小三角形展開全部。
   - SNV/Indel 顯示 filter 預設啟用 `Disease-associated`（優先使用 `ngs_panel_deadzone/panel/panel_loose_plus_clinical.hgnc_canonical.txt`，缺檔時 fallback 到 `panel_loose.hgnc_canonical.txt`）、`In panel only`、`gnomAD_G_AF < 0.01`、`VAF ≥ 0.2`；`impact=MODIFIER` 預設不顯示，可手動勾選展開。`IMPACT=LOW` 仍會顯示。CNV/SV 與 gene search modal 不受 `Disease-associated` filter 限制。
   - TSV post-processing 只移除 `REF/ALT=*` 與非 primary contig；正式 v3.1 pipeline 會自行處理 NCKUH/DRAGEN 前處理與 PASS/chrM 分流，舊版 DRAGEN staging 僅在 `NGS_UI_TERTIARY_LEGACY_STAGING=1` 時啟用。
   - 主畫面讀取自動衍生的 `snv_indel.review.tsv`：保留 `NGS_UI_CDS_CANDIDATE_BED` 內且 `GNOMAD_G_AF < 0.01` 或 AF 缺值的位點，並 rescue ClinVar P/LP。reviewer 已標記但不在 review TSV 的 SNV 會用 `snv_gene_index.sqlite` 依 variant id 補入，不會因標記狀態變動而重建 review TSV。三級分析的 post-processing 階段會先建立它；舊樣本載入時仍可自動補建。原始 `snv_indel.annotated.tsv` 不會被覆寫。
   - SNV/Indel 卡片標籤列會依 GIAB genome-stratification 標出困難區（homopolymer / tandem repeat / segdup / low mappability / GC extreme / other difficult）的琥珀色 badge，提醒 reviewer 該位點 short-read calling 較不可靠。資料由 `scripts/annotate_giab_strata.py` 在三級分析尾段寫入 `snv_indel.annotated.tsv` 的 `GIAB_STRATA` 欄，BED 與 `strata_manifest.json` 放在 `NGS_UI_GIAB_STRAT_DIR`（部署時用 `scripts/download_giab_strata.sh` 下載；不在 git 內），缺 BED 目錄則自動略過。純 UI 提示，不影響 tier 排序或診斷報告。
   - ACMG SF / Proactive / Carrier screening 只列 custom panel `ACMG_SF_v3.3`、`proactive`、`carrier_mackenzie_1300+` 內的 SNV/Indel，且需符合 ClinVar P/LP 或卡片最終 ACMG P/LP（GeneBe 優先、再用 pipeline ACMG）。ClinVar P/LP 預設 ✓ 並進 Secondary findings 報告區；其他 ACMG P/LP 預設留在分析區，reviewer 可手動勾 ✓。手動取消會寫入 `secondary_findings.{section}.dismissed`，不會下次載入又被預設勾回。
   - GeneBe ACMG 第二意見（`GENEBE_ACMG_*` 欄）改由本機 GeneBe 資料庫離線提供，取代原本的 GeneBe API（不需帳號/網路）：post-processing 對整張 TSV 的變異查 `NGS_UI_GENEBE_DB`（預設 `biotools/genebe/genebe_hg38.tsv.gz`）的 `acmg_score`/`acmg_criteria`。預設會在同目錄 lazy 建立 `genebe_hg38.sqlite` key-value cache；若上傳/替換新的 `genebe_hg38.tsv.gz`，下一次三級分析會依來源檔 size/mtime/ctime 偵測變更、用 file lock 重建 SQLite，之後 sample annotation 直接查 SQLite，不再每個 sample 串流掃整顆 GeneBe TSV。SQLite 建置或查詢失敗時會 fallback 到舊的 single-pass streaming。DB 為預算好的快取，查不到的點（如 novel coding indel）就沒有 GeneBe 第二意見，pipeline 自己的 ACMG 仍會顯示。DB 部署/更新用 `scripts/deploy_genebe_db.sh`（濾壞行 → bgzip → tabix → 原子換檔），建置端要求見 `docs/genebe_db_requirements.md`。
   - SNV tier 只在點開時建立該 tier 的卡片 DOM，避免一次 render 全部卡片。
   - SNV/Indel 與 CNV/SV gene 搜尋支援多個基因，以 `,` 或 `、` 分隔；SNV 搜尋由 `/api/samples/{id}/snv-search` 查完整原始 TSV，不受 review TSV 限制。三級分析結尾會預建 `snv_gene_index.sqlite`（gene → raw TSV byte offsets），所以 WGS gene search 不需在載入個案時掃 1–2GB raw TSV；舊樣本若缺 index 才 fallback 到 raw TSV parse。modal 預設勾選 `gnomAD_G_AF < 0.01`，取消後才顯示全部搜尋結果。搜尋 modal 內的 SNV 卡片會依 `MANE_ALL` 優先用有 HGVS.c/p、非 MODIFIER 的 MANE_SELECT `NM_` transcript 顯示，且可容忍不同 TSV 引號格式。
   - Variant 狀態用 `1 / 2 / C / 0` 圓形按鈕；`C` 與 `0` 可並存，`1` 或 `2` 會反選其他狀態；再次點擊已選項目可清空狀態。同一個 variant 在分析區、報告區與搜尋 modal 的按鈕會同步上色。
   - Variant 卡片保留 OMIM `Disease1..5` 完整內容；「個案清單」中的疾病摘要才截到第一個可辨識的遺傳模式括號（例如 `(AD)`、`(AR)`）為止。
   - SNV/Indel、Mitochondria 與 CNV/SV 卡片有 `IGV` 按鈕；基本資料的性別下方另有「☑ 已確認」核取方塊與「確認 SRY」按鈕。SRY 會沿用同一個 modal 開啟 hg19/hg38 對應區域，初始只載入當前 sample，仍可手動加入 sibling，coverage data range 固定 `0-100`；勾選確認狀態會寫進 reviewer metadata。一般 modal 標題會顯示 sample、padded locus、variant 註解與原始座標；alignment 預設用 squished 模式，`visibilityWindow=5 Mb`，SNV/Indel 與 Mito track height 固定 300，CNV/SV track height 固定 50。SNV/Indel 與 Mito 顯示前後 100bp，CNV/SV 原則上顯示前後各 20% flanking area；若事件本身 ≤5 Mb，padding 會自動縮小以維持初始視窗 ≤5 Mb，直接載入 BAM coverage，不需先 zoom in。IGV 另以 ROI 標出實際 CNV/SV 區間。CNV/SV modal 會把所有 BAM coverage tracks 放進同一個 autoscale group，讓 y-axis data range 一致，方便和 sibling 比較 deletion/duplication。先確認 primary BAM 與同 batch sibling tracks，再按「載入 IGV」；若 UI sample ID 帶 `-dragen` / `-nckuh` / legacy `-inhouse` / `-WES` / `-WGS`，BAM 查詢可 fallback 到去 suffix 的原始 sample ID。DRAGEN 來源會在 raw Novaseq root 找 `<run>/bam/{sample}.bam` 並排除 `{sample}.repeats.bam`；in-house 來源維持找 Nextflow output 的 `02_alignment`。若自動搜尋不適用，可按「其他路徑」從下拉選擇 DRAGEN run 或 in-house batch，列出該 run/batch 可用 BAM 後選 primary；同 batch 加入清單會改用同一資料夾的其他 BAM。BAM range request 與 hg38 FASTA 都由後端在內網 proxy。
   - CNV：`CNV-1A`（Clinical）、`CNV-1B`（Pathogenic）；SV：`SV-2A / SV-2B`。CNV/SV 分析區與報告區皆依 `max_pheno_score + scaled AnnotSV ranking score` 由高到低排序；CNV/SV 卡片的基因表預設只顯示 phenotype 相關基因，按小三角形才展開 phenotype score 為 0 的其餘基因。
   - 同來源、同染色體、同為 deletion 或同為 duplication，且相鄰 gap ≤ `250 kb` 的 CNV/SV 會預設自動整合；copy number 差異不再阻擋視覺整合，原始片段仍可展開查看各自 CN。UI、DOCX 與個案清單改用合併後 parent；parent 會取代最佳原始 segment 的位置，不會掉到 tier 最後。前端會依 sample payload 快取整合結果，避免卡片 render 時反覆重建 parent。
   - Mitochondria：`MITO-1`（ClinVar P/LP 或 MITOMAP confirmed/pathogenic）、`MITO-2`（rare / reported mtDNA variant）、`MITO-3`（other variant）— 若來源 TSV 有 `FILTER` 欄，只列 `FILTER=PASS`；v3.2 `04_mito/{sample}.mito.tsv` 沒有 `FILTER` 時不顯示 Filter 欄。MITOMAP 欄位只在後端分類使用，卡片不顯示 MITOMAP；Disease 改用 `CLINVAR_DN` 依 `&` 拆成 checkbox，勾選者才進 DOCX。
4. **標記與判讀** — 在每個變異上標 causative / candidate / other，編輯 ACMG/分類、寫 comment；SNV/Indel ACMG 優先序為 reviewer override → GeneBe → pipeline `ACMG_CLASS`。Mito ACMG 下拉會同步更新分析區與報告區卡片，CNV/SV 卡片另有 Disease 欄供 DOCX 與個案清單使用。變更會自動存到 `tertiary_output/{LIS}/sample_metadata.json`；開啟個案清單或匯出 DOCX 前會先 flush 尚未完成的自動儲存。

SNV/Indel 卡片的 ESM1b 依 ClinGen SVI 校準區間上色；ESM1b 分數越低越偏致病，多 transcript TSV 會取最低分作為 worst case。
5. **匯出報告** — 「匯出診斷報告」下載 `GET /api/samples/{LIS}/report.docx`；DOCX 依序輸出第一類、第二類、固定建議文字，再集中列出各 variant 的參考資料。CNV/SV Disease 欄會優先覆蓋單基因預設疾病，也可為片段型 CNV/SV 指定發報告疾病；多基因或無 OMIM gene 的片段會把 Disease 直接接在第 1 點位置描述，不另列一點。未涵蓋 OMIM 疾病相關基因時不再列出一般基因清單。SNV/Indel 會帶 TSV 的 `RS_ID`，RS ID 欄會預留尾端空格；SNV/Indel 與 Mito 核苷酸欄每行最多 13 字元且蛋白質括號另起一行。§五.4 grouped gene list 的 HPO/panel 區塊之間會留空行，且只列 expanded reportable disease-associated genes；WGS 報告若基因有 dead-zone exon，會以 `GENE（exon 2, 4-6 <15X）` 這類短括號標註，基因清單後空一行並以「註：括號中標示之 exon 為 cohort dead-zone，代表該 exon coverage 低於本檢測判讀門檻」說明；WES 報告不輸出 dead-zone 標註。旁邊的「輸出 PDF」會開啟列印視窗，輸出報告區的 causative / other / candidate 卡片摘要，可由瀏覽器另存為 PDF。列印版會略過 comment、More、Secondary findings 與 CNV/SV overlap 明細。

三級分析 modal 的 in-house 與 DRAGEN VCF 各自使用單一 typeahead 輸入格；點入輸入格即展開全部 VCF，輸入 sample / run / path 後即時縮小候選清單，候選列顯示 sample、run、檔案大小與建立時間；選定後輸入框只保留 sample name，完整資訊留在 title。worker 會保留原始 `source_sample_id`，並為 v3.x pipeline 產生 sample sheet 後執行 `nextflow ... --samplesheet ... --pipeline_type nckuh|dragen --out_dir /home/pipeline/tertiary_output -resume`。pipeline 完成後，NGS-UI 逐筆複製 SNV/Mito/CNV/SV TSV 並寫 `pipeline_source.json`；若 pipeline CNV/SV 兩檔都複製成功，post-processing 會跳過本機 AnnotSV fallback，缺任一檔才保留 fallback。執行進度以 step-based 進度條顯示，Nextflow 會從 stdout 追蹤 PREPARE_VCF / VEP / Pangolin / CSQ parse / ACMG 等內部 process；完成事件不再插入獨立 `[nextflow-step]` 或 `[step] nextflow:*` 行，而是在原本的 Nextflow process 行尾補上 `elapsed=...m [YYYY-MM-DD HH:MM:SS]`。post-processing 只追蹤實際執行的 GeneBe、extra VEP、AnnotSV，review TSV 與 gene index 則列為 sample 預建步驟。詳細 log 預設收合，跳過的 Extra VEP / AnnotSV 不留下空 section；Nextflow process history 仍另存於 `state.json.nextflow_step_history`，`step_history` 則保留 UI 進度 milestone。「三級分析清單」會合併 `/home/pipeline/tertiary_output/` 與 NGS-UI job state，因此失敗且尚未建立 output 的 sample 仍會顯示 log；刪除時會一併刪除 pipeline output、`$NGS_UI_HOME/tertiary_output/{sample}/` 與該 sample 的 NGS-UI job log 目錄，執行中的 sample 不可刪除。新 job 狀態寫在 `data/jobs/tertiary/`；舊版 `data/jobs/dragen/` 紀錄仍可讀取。底層相容腳本仍名為 `run_stopgaps.sh`，但 log / UI 一律顯示為 post-processing。AnnotSV 的目前用法另見 `docs/annotsv_current_usage.md`。

三級分析重用既有 `/home/pipeline/tertiary_output/{sample}/` output 時，優先使用 UI sample ID（例如 `VAL-57-dragen` / `VAL-57-nckuh`）作為 pipeline output 目錄，且後端 job API / worker 會為新 job 強制補上來源 suffix，避免 DRAGEN 與 in-house 同名 sample 互相覆蓋；舊的 source-ID-only output 只作 legacy fallback。reuse 也不再只看 `03_acmg/*.snv_indel.acmg.tsv`，而是檢查 SNV、Mito、CNV、SV；PGx checkbox 有勾時也會要求 PGx/PharmCAT 輸出存在，缺任一項就讓 Nextflow 用 `-resume` 補齊。已跑完但 pipeline output 還沒有 suffix 的 case，可先 dry-run `python scripts/repair_pipeline_output_suffixes.py`，確認後加 `--apply` 複製成 suffixed 目錄；工具會依 `pipeline_source.json` 修復，也會預設把 `/home/pipeline/tertiary_output/VAL-數字/` 視為 legacy DRAGEN 並複製成 `VAL-數字-dragen/`（可用 `--legacy-dragen-pattern` 調整或設空字串關閉）。worker-owned log timestamp 會放在行尾 `[YYYY-MM-DD HH:MM:SS]`；Nextflow process 完成時只在原 stdout 行尾追加 elapsed 分鐘與時間戳，batch 進度條會依同一 process 的完成比例推進。

三級分析進度條依實測 DRAGEN / in-house batch log 配重：Nextflow 佔主要時間，但 `queued` 事件只寫 timing log，不推進 UI 百分比；實際 `start/done` 才依 process 權重更新。前處理、mito、STR、prepare CNV/SV 等快速步驟只佔少量進度，`ANNOTSV_SV`、VEP、Pangolin、parse CSQ、ACMG/PGx 等耗時步驟佔主要權重。Nextflow 結束約落在 82%，post-processing 依目前第幾個 sample 與子步驟推進到完成。Nextflow batch 會同時更新多個 process，因此後端寫入 `nextflow_progress_pct` 時會保留目前最大值，避免較早 process 晚更新造成進度條倒退。首頁會定期查詢後端 active job，因此其他瀏覽器或電腦登入後也會看到正在跑的三級分析；進度面板提供「終止」按鈕，會呼叫後端取消該 job 並向 worker process group 送出終止訊號。

DOCX CNV/SV 表格的「變異位置」欄使用 buffered wrap，內容寬度比欄寬少一格，避免座標字串貼到後面的欄位。

可用 `scripts/compare_genebe_spliceai_coverage.py` 抽樣 GeneBe 本機資料庫（例如 `/home/n102968/NGS_UI/biotools/genebe/genebe_hg38.tsv.gz`），再走目前 extra-VEP 的 VEP SpliceAI plugin 路徑產生對照，輸出 `summary.tsv`、`summary_by_kind.tsv`、`mismatches.tsv` 與 `run_metadata.tsv`，評估 GeneBe DB 的 SpliceAI 覆蓋率是否足以取代 extra-VEP / GeneBe API。正式估計建議使用 `--max-sites 100000 --sample-mode chrom-balanced`，避免只取 TSV 前段或讓大型染色體主導結果；未指定 `--seed` 時會由系統亂數產生並記錄在 metadata。

服務啟動時會以 daemon thread 在背景預熱 HPO、phenotype gene map、OMIM 與 mito ClinVar cache，避免解析大型 `phenotype_to_genes.txt` 期間擋住 HTTP port。完整 SNV TSV 的 gene search 走 per-sample `snv_gene_index.sqlite`；index 由 tertiary job 的 post-processing 預建，載入個案後不再自動掃描 WGS raw TSV。Secondary findings 會用同一個 gene index 只補齊 ACMG SF / proactive / carrier panel 內且符合 ClinVar P/LP 或 ACMG P/LP 的 SNV/Indel，不把整個 carrier panel 的所有 variant 帶進前端。`in_panel` / `pheno_score` 由 active analysis 的 `pheno_score.tsv` 在載入與搜尋時補入，不再為 HPO/panel 變更重寫完整 raw TSV。刪除個案時只 invalidate 該 sample 的 SNV 與個案摘要 cache，不會同步重掃全部清單。

三級分析另保存 `source_sample_id`（原始 sequencing sample ID）到 job state 與 `pipeline_source.json`。IGV 先依 sample ID suffix 判斷來源，`-dragen` 強制 DRAGEN roots、`-nckuh` / legacy `-inhouse` 強制 in-house roots；沒有 suffix hint 時才看 `pipeline_source.json` 的 `pipeline_type`。DRAGEN 找 `NGS_UI_DRAGEN_BAM_ROOT/<run>/bam/{sample}.bam`，in-house 找 `NGS_UI_INHOUSE_BAM_ROOT/<batch>/{sample}/02_alignment/`。自訂輸出 ID 找不到時，會用 sidecar 回查原始 BAM ID；舊個案缺少 sidecar 欄位時，會移除已知 suffix 後用原始 sample ID 搜尋 BAM。IGV 的「其他路徑」資料夾下拉由 `/api/igv/bam-folders` 列出 configured BAM roots 內的 DRAGEN `<run>/bam/` 與 in-house `<batch>`，選 in-house batch 後會列出該 batch 底下各 sample `02_alignment` 的 BAM，排除 `.repeats.bam`。

### 臨床表徵工具 `/phenotype/`

獨立頁面（內網信任、無需登入）：搜尋 HPO term、套用 / 自訂 gene panel、把結果存成 token 之後在「載入新個案」帶入。Token 限 `[A-Za-z0-9_-]{1,32}`，內容 ≤64KB，panel ≤5000 個基因；自訂 panel 的基因 symbol 不會被轉大寫（`C7orf50` 保持原樣），但會先套用安全 alias 轉成 HGNC-current symbol。

主畫面的 Patient phenotype card 也提供 `WES-I / WES-II / WGS / Other panel` tabs，預設展開 `Other panel`；固定 panel 可直接點 chip 選取，其他 panel 仍可用 typeahead 搜尋。

固定 panel 的來源是 `reference/fixed_panel_sources/WES-I.xlsx`、`reference/fixed_panel_sources/WES-II.xlsx` 與 `reference/fixed_panel_sources/other_panel/`。更新 Excel 後執行 `PYTHONPATH=backend python scripts/import_fixed_panels.py`，會同步重建三個入口共用、且會進 git 的 `phenotype_data/fixed_panels/index.json` 與 `phenotype_data/gene_panels/*.txt`。WES Excel 只會匯入 `gene panel list` 標記列起始的基因區塊，避免把疾病名或資料來源列誤算成基因。

HPO reference、固定 panel、custom panel 與既有 `pheno_score.tsv` 讀入時都會先透過 `ngs_panel_deadzone` 的 HGNC alias map 轉成 canonical gene symbol；從 `/phenotype/` 建立 custom panel 時，後端也會先套用 `ngs_panel_deadzone/panel/panel_gene_aliases.tsv` 的安全 alias 再寫入 repo 內 panel 檔。`panel_gene_aliases.tsv` 由 HGNC 官方 `reference/hgnc/hgnc_complete_set.txt`（`prev_symbol` / 唯一 `alias_symbol`）與 `reference/hgnc/withdrawn.txt`（唯一 merged/split replacement）加上 `reference/hgnc/manual_panel_aliases.tsv` 產生，重建指令為 `python scripts/build_hgnc_panel_aliases.py`；衝突項輸出到 `docs/ops/hgnc_alias_conflicts.tsv`，既有 custom panel 仍非 current HGNC 的項目列在 `docs/ops/custom_panel_hgnc_review_20260613.tsv`。Custom panel 檔案第一行可用 `#source:` 記錄來源（空白也可），loader 會略過註解行；`/phenotype/` 的 gene-list drawer 會顯示這個 source。SNV/CNV/SV/Mito 變異端也用同一套 canonicalization，再做 `pheno_score` / `in_panel` join。

`/phenotype/` 的 HPO term、fixed panel chip 與 panel 搜尋列都有「查看」按鈕，會呼叫 `GET /api/phenotype-tool/gene-list?kind=hpo|panel&key=...`，在右側 drawer 顯示 canonical gene list、來源、清單內篩選與複製功能，不把完整基因清單塞進主畫面。Topbar 的「搜尋基因」會呼叫 `GET /api/phenotype-tool/gene-memberships?gene=...`，反查某個 canonical gene 出現在哪些 HPO terms / panels。單一 HPO gene-list 查詢在 full phenotype scorer cache 尚未預熱完成時會走 fast path，只掃該 HPO term，避免第一次點「查看」被整份 `phenotype_to_genes.txt` 載入卡住。

---

## 5. 新增一個個案 / mitochondrial annotation

次級 pipeline 通常會直接把 `snv_indel.annotated.tsv` 等放進 `tertiary_output/{LIS_ID}/`。目前三級 pipeline v3.2 會輸出 `04_mito/{SAMPLE_ID}.mito.tsv`，NGS-UI worker 會在三級分析完成後複製成 `tertiary_output/{LIS_ID}/mito.annotated.tsv`，並靜默用本地 MITOMAP 表回補 `MITOMAP_*` 欄位（不額外寫進 job log）。前端 adapter 同時相容舊欄位與 v3.2 欄位。v3.2 的 `HGVS_C` 是 transcript `c.`，所以 UI 的 mito 卡片標題改由 `POS/REF/ALT` 產生 mitochondrial `m.` 寫法；`HGVS_P` 會清成報告用的 `p.xxx`。Mito 判讀以 ClinVar P/LP、MITOMAP confirmed/pathogenic、gnomAD-mito rare 或缺值/已報告變異分成三層；`CLINVAR_DN` 會拆成可勾選 disease，MITOMAP 不顯示在卡片中。

舊樣本或手動補跑時，仍可用本 repo 的 legacy script 從 GATK Mutect2 `--mitochondria-mode` 的 VCF 產生 mito TSV（純 Python，只需 VCF + 本地 MITOMAP 表，不需 VEP/bcftools）：

```bash
# MITOMAP_DIR 預設 ${REF_DIR}/tertiary/mitomap，內含
#   mitomap_mutations_coding_control.tsv  與  mitomap_mutations_rna.tsv
scripts/annotate_mito_vcf.sh \
  --in   /path/to/{LIS_ID}.mito.vcf.gz \
  --sample {LIS_ID} \
  --outdir tertiary_output/{LIS_ID}/
# → tertiary_output/{LIS_ID}/mito.annotated.tsv
```

DRAGEN run 若 mtDNA calls 直接包在 `{sample}.hard-filtered.vcf.gz` 的 `chrM/MT` rows，可在 dev 機用專用 wrapper 跑同一套 MITOMAP-only 流程：

```bash
scripts/annotate_dragen_mito_vcf.sh \
  --in /home/datalake_Raw/Novaseq/20260428_LH00873_0015_B23NG3WLT4/vcf.gz/VAL-33.hard-filtered.vcf.gz \
  --sample VAL-33 \
  --outdir /home/n102968/NGS_UI/tertiary_output/VAL-33_legacy_mito
```

批次：

```bash
for s in 26WE0043 26WE0044 26WE0045 26WE0046 26WE0047 26WE0048 26WE0074; do
  scripts/annotate_mito_vcf.sh --in vcf/$s.mito.vcf.gz --sample $s \
    --outdir tertiary_output/$s/
done
```

> 注意：MITOMAP 那兩個 TSV 是 **Latin-1** 編碼（含 0xa0 byte），不是 UTF-8。
> 其他轉檔/遷移/診斷 script 見 `scripts/`（`convert_anno_combined_to_tertiary_tsv.py`、`compare_genebe_spliceai_coverage.py`、`migrate_layout.sh`、`migrate_to_versioned_layout.py` 等）與 `CLAUDE.md`。

---

## 6. 帳號管理

```bash
PYTHONPATH=backend python3 -m app create-user [USERNAME]   # 互動輸入密碼（bcrypt，上限 72 bytes）
PYTHONPATH=backend python3 -m app list-users
```

Session cookie 8 小時、`SameSite=Lax`、`https_only=False`（內網可能還沒 HTTPS）。

---

## 7. 注意事項

- **不要把病人資料 / 大檔 commit 進 git**：`.gitignore` 已排除 `tertiary_output/`、`data/`、`patient_list/`、`phenotype_data/` 內的大型 HPO reference、`_index.json`；但固定與自訂 panel 的 `phenotype_data/fixed_panels/`、`phenotype_data/gene_panels/` 與 `phenotype_data/custom_panels/` 例外追蹤。
- 首頁歡迎文字與版本紀錄放在 `frontend/VERSION.md`，前端啟動時會讀取並顯示在尚未載入個案的首頁。之後若有影響判讀流程、報告輸出、資料載入或主要工具入口的更新，需評估是否同步更新這份版本紀錄。
- 大型 HPO reference 必須放在 `NGS_UI_PHENO_DATA_DIR`（`hp.obo`、`phenotype_to_genes.txt` 等）；固定與自訂 panel data 則在 repo 的 `phenotype_data/fixed_panels/`、`phenotype_data/gene_panels/` 與 `phenotype_data/custom_panels/`，會跟著 git 更新。
- EMR 相關功能預設停用，需設 `NGS_UI_EMR_CLIENT_ID` 才會啟用，且只在內網可達。
- `/api/phenotype-tool/*` 與 `/api/healthz` 是刻意公開無認證；`/api/patient_list` 與其餘 `/api/*` 需登入。
- 大型 JSON response 會在瀏覽器支援時自動 gzip；SNV parse + phenotype / Exomiser / LIRICAL / OMIM join 使用有上限的 process-local LRU cache，輸入 TSV 或 sidecar 更新後自動失效。
