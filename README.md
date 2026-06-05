# NGS 分析平台 (NGS-UI)

成大醫院基因醫學部的 NGS 三級分析判讀工具。次級 pipeline（Nextflow，跑在另一台 compute cluster）產出 per-sample 的註解 TSV，本平台讓 reviewer 載入個案、檢視 SNV/Indel + CNV/SV + Mitochondria 變異、標記 causative / candidate / other、撰寫判讀意見，並匯出診斷報告 (docx)。另附一個獨立的「臨床表徵輸入 (HPO / gene panel)」工具掛在 `/phenotype/`。

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
│                             cnv.annotated.tsv,
│                             sv.annotated.tsv, mito.annotated.tsv,
│                             sample_metadata.json, qc_summary.json,
│                             roh_summary.json, analyses/{ver}/...
├── patient_phenotype/     ← {LIS_ID}_{MRN}_phenotype.txt（自動帶入 HPO）
├── patient_list/          ← 上傳的「未完成報告清單」xlsx + 衍生 roster.json
├── phenotype_data/        ← hp.obo, phenotype_to_genes.txt, gene_panels/*.txt
├── ngs_panel_deadzone/    ← panel_loose、HGNC alias map、dead-zone tables
├── OMIM/OMIM.xlsx         ← OMIM 疾病註解表（缺檔則 Disease 欄留空）
└── data/                  ← server runtime state（users.db, jobs/, ...）
```

開發 checkout 不需要這整棵樹：`NGS_UI_HOME` 未設且找不到上層 `NGS-UI/` 時，會 fallback 成 repo 自己，所有路徑都落在 repo 內。

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
#    參考 scripts/migrate_layout.sh，重點是讓 phenotype_data 等資料放在
#    NGS_UI_HOME 底下而非 repo 內：
#      mv NGS_UI/NGS-UI/phenotype_data NGS_UI/phenotype_data
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
| `NGS_UI_PHENO_DATA_DIR` | `$NGS_UI_HOME/phenotype_data` | hp.obo / phenotype_to_genes / gene_panels |
| `NGS_UI_PANEL_DEADZONE_DIR` | `REPO_ROOT/ngs_panel_deadzone` | panel_loose、HGNC alias map、WES/WGS dead-zone tables |
| `NGS_UI_OMIM_XLSX` | `$NGS_UI_HOME/OMIM/OMIM.xlsx` | OMIM 疾病註解（缺檔 = 停用） |
| `NGS_UI_BIOTOOLS_DIR` | `$NGS_UI_HOME/biotools` | Exomiser / LIRICAL |
| `NGS_UI_BAM_ROOT` | `/home/datalake_Intermediate/pipeline/nextflow_output` | IGV 搜尋 BAM 的根目錄；可用 `:` 分隔多個 root |
| `NGS_UI_IGV_REF_DIR` | `/home/pipeline/reference/hg38` | IGV 本機 hg38 FASTA + `.fai` |
| `EXOMISER_HOME` / `LIRICAL_HOME` / `EXOMISER_DATA_HG38` ... | `biotools/...` | 工具與 data 路徑 |
| `JAVA_BIN` / `JAVA_OPTS` | `java` / `-Xms4g -Xmx16g` | 跑 Exomiser/LIRICAL 用 |
| `REDIS_URL` | `redis://127.0.0.1:6379/0` | RQ 佇列 |
| `NGS_UI_EMR_CLIENT_ID` | `""`（空 = 停用所有 EMR 路徑） | NCKU 內網 HIS / APIM |

---

## 4. Reviewer 操作流程

1. **登入** — 右上角登入（帳號由管理者用 `create-user` 建立）。
2. **載入新個案** — 點「載入新個案」：
   - LIS_ID 下拉會列出 pipeline 已丟進 `tertiary_output/` 但尚未登錄的目錄；
   - 若先用「上傳個案清單」匯入過「未完成報告清單」xlsx，MRN / 姓名 / Test type 會自動帶入（來自 `patient_list/roster.json`）；
   - HPO / gene panel 可在這裡選；gene panel 與主畫面同樣使用 `WES-I / WES-II / WGS / Other panel` tabs，預設展開 `Other panel`，固定 panel chip 與搜尋下拉都會顯示基因數量；若存在 `patient_phenotype/{LIS}_{MRN}_phenotype.txt` 會自動讀入；
   - 勾「登錄後開始分析」會順便把 Exomiser/LIRICAL 排入佇列。
   - 旁邊的「個案清單」可查看已載入個案、依 `WES` / `WGS` 篩選並全文搜尋；表格會摘要 causative / other SNV、CNV、SV，已勾選 OMIM disease、CNV/SV reviewer 輸入的 Disease、主畫面 comment、簽收與載入時間。多個 variant / disease 會逐行顯示，長 HGVS 可自動折行。也可刪除 `NGS_UI_HOME/tertiary_output/{LIS_ID}/`；刪除時可選擇同步刪除或保留 `/home/pipeline/tertiary_output/{LIS_ID}/`，按取消則完全不刪除。
   - 主畫面搜尋框上方有可複選的 `WES` / `WGS` 圓形 filter，取消勾選後對應 test type 不出現在搜尋清單。
3. **看變異卡片** — 個案載入後先顯示 SNV/Indel（分段載入），CNV/SV 與 Mitochondria 在背景載完後補上：
   - 平台剛開啟讀取索引、個案核心資料載入與新個案登錄期間都會顯示不可誤關閉的「資料載入中」遮罩，避免重複點擊。
   - SNV/Indel tier：`1A / 1B / 1C / 2 / 3`（互斥）
   - Patient phenotype 與 Comment 之間的 Dead zone 卡片會依目前 HPO + panel gene set 顯示 cohort-level dead exons；WES 用 20X，WGS（含 in-house / DRAGEN）用 DRAGEN 15X。主畫面預設只顯示前 10 列，可用小三角形展開全部。
   - SNV/Indel 顯示 filter 預設啟用 `Disease-associated`（限 `ngs_panel_deadzone/panel/panel_loose.hgnc_canonical.txt`）、`In panel only`、`gnomAD_G_AF < 0.01`、`VAF ≥ 0.2`；`impact=MODIFIER` 預設不顯示，可手動勾選展開。`IMPACT=LOW` 仍會顯示。CNV/SV 與 gene search modal 不受 `Disease-associated` filter 限制。
   - TSV stop-gap 只移除 `REF/ALT=*` 與非 primary contig；正式 v3.1 pipeline 會自行處理 NCKUH/DRAGEN 前處理與 PASS/chrM 分流，舊版 DRAGEN staging 僅在 `NGS_UI_TERTIARY_LEGACY_STAGING=1` 時啟用。
   - 主畫面讀取自動衍生的 `snv_indel.review.tsv`：保留 `NGS_UI_CDS_CANDIDATE_BED` 內且 `GNOMAD_G_AF < 0.01` 或 AF 缺值的位點，並 rescue ClinVar P/LP 與 reviewer 已標記點。`run_stopgaps.sh` 在三級分析結尾先建立它；舊樣本載入時仍可自動補建。原始 `snv_indel.annotated.tsv` 不會被覆寫。
   - SNV tier 只在點開時建立該 tier 的卡片 DOM，避免一次 render 全部卡片。
   - SNV/Indel 與 CNV/SV gene 搜尋支援多個基因，以 `,` 或 `、` 分隔；SNV 搜尋由 `/api/samples/{id}/snv-search` 查完整原始 TSV。登錄新個案或載入既有個案後都會在背景預熱 raw TSV cache；modal 預設勾選 `gnomAD_G_AF < 0.01`，取消後才顯示全部搜尋結果。
   - Variant 狀態用 `1 / 2 / C / 0` 圓形按鈕；`C` 與 `0` 可並存，`1` 或 `2` 會反選其他狀態；再次點擊已選項目可清空狀態。同一個 variant 在分析區、報告區與搜尋 modal 的按鈕會同步上色。
   - Variant 卡片保留 OMIM `Disease1..5` 完整內容；「個案清單」中的疾病摘要才截到第一個可辨識的遺傳模式括號（例如 `(AD)`、`(AR)`）為止。
   - SNV/Indel 與 CNV/SV 卡片有 `IGV` 按鈕；基本資料的性別下方另有「☑ 已確認」核取方塊與「確認 SRY」按鈕。SRY 會沿用同一個 modal 開啟 hg19/hg38 對應區域，初始只載入當前 sample，仍可手動加入 sibling，coverage data range 固定 `0-100`；勾選確認狀態會寫進 reviewer metadata。一般 modal 標題會顯示 sample、padded locus、variant 註解與原始座標；alignment 預設用 squished 模式，`visibilityWindow=5 Mb`。SNV/Indel 顯示前後 100bp，CNV/SV 原則上顯示前後各 20% flanking area；若事件本身 ≤5 Mb，padding 會自動縮小以維持初始視窗 ≤5 Mb，直接載入 BAM coverage，不需先 zoom in。IGV 另以 ROI 標出實際 CNV/SV 區間。CNV/SV modal 會把所有 BAM coverage tracks 放進同一個 autoscale group，讓 y-axis data range 一致，方便和 sibling 比較 deletion/duplication。先確認 primary BAM 與同 batch sibling tracks，再按「載入 IGV」；BAM range request 與 hg38 FASTA 都由後端在內網 proxy。
   - CNV：`CNV-1A`（Clinical）、`CNV-1B`（Pathogenic）；SV：`SV-2A / SV-2B`。CNV/SV 卡片的基因表預設只顯示 phenotype 相關基因，按小三角形才展開 phenotype score 為 0 的其餘基因。
   - 同來源、同染色體、同為 deletion 或同為 duplication、copy number 相容且相鄰 gap ≤ `250 kb` 的 CNV/SV 會預設自動整合。UI、DOCX 與個案清單改用合併後 parent，原始片段仍可展開查看；parent 會取代最佳原始 segment 的位置，不會掉到 tier 最後。前端會依 sample payload 快取整合結果，避免卡片 render 時反覆重建 parent。
   - Mitochondria：`MITO-1`（Pathogenic）、`MITO-2`（Disease-associated）— 只列 `FILTER=PASS` 且具 MITOMAP 疾病關聯/致病性的位點
4. **標記與判讀** — 在每個變異上標 causative / candidate / other，編輯 ACMG/分類、寫 comment；SNV/Indel ACMG 優先序為 reviewer override → GeneBe → pipeline `ACMG_CLASS`。CNV/SV 卡片另有 Disease 欄供 DOCX 與個案清單使用。變更會自動存到 `tertiary_output/{LIS}/sample_metadata.json`；開啟個案清單或匯出 DOCX 前會先 flush 尚未完成的自動儲存。

SNV/Indel 卡片的 ESM1b 依 ClinGen SVI 校準區間上色；ESM1b 分數越低越偏致病，多 transcript TSV 會取最低分作為 worst case。
5. **匯出報告** — 「匯出診斷報告」下載 `GET /api/samples/{LIS}/report.docx`；DOCX 依序輸出第一類、第二類、固定建議文字，再集中列出各 variant 的參考資料。CNV/SV Disease 欄會優先覆蓋單基因預設疾病，也可為片段型 CNV/SV 指定發報告疾病；多基因或無 OMIM gene 的片段會把 Disease 直接接在第 1 點位置描述，不另列一點。未涵蓋 OMIM 疾病相關基因時不再列出一般基因清單。SNV/Indel 會帶 TSV 的 `RS_ID`，RS ID 欄會預留尾端空格；SNV/Indel 與 Mito 核苷酸欄每行最多 13 字元且蛋白質括號另起一行。§五.4 grouped gene list 的 HPO/panel 區塊之間會留空行，且只列 `panel_loose` disease-associated genes；WGS 報告若基因有 dead-zone exon，會以 `GENE（exon 2, 4-6 <15X）` 這類短括號標註，基因清單後空一行並以「註：括號中標示之 exon 為 cohort dead-zone，代表該 exon coverage 低於本檢測判讀門檻」說明；WES 報告不輸出 dead-zone 標註。旁邊的「輸出 PDF」會開啟列印視窗，輸出報告區的 causative / other / candidate 卡片摘要，可由瀏覽器另存為 PDF。列印版會略過 comment、More、Secondary findings 與 CNV/SV overlap 明細。

三級分析 modal 的 in-house 與 DRAGEN VCF 各自使用單一 typeahead 輸入格；點入輸入格即展開全部 VCF，輸入 sample / run / path 後即時縮小候選清單。Sample ID 可修改作為 UI 輸出資料夾與檔名前綴；選好後可按「加入批次」，同一批只能包含同一種來源（全 in-house 或全 DRAGEN），按「開始分析」後會用多列 v3.1 sample sheet 一次送出。worker 會另外保留原始 `source_sample_id`，並為 v3.1 pipeline 產生 sample sheet（`sample_id,pipeline_type,input_dir,seq_type,hpo`），再執行 `nextflow ... --samplesheet ... --pipeline_type nckuh|dragen --out_dir /home/pipeline/tertiary_output -resume`；若 `NGS_UI_TERTIARY_ENV_SCRIPT`（預設 `/home/pipeline/pipeline_code/NGS2ndAnalysis_env.sh`）存在會先 source，不存在則沿用目前環境直接執行 `nextflow`。v3.1 pipeline 輸出以 source sample ID 落在 `/home/pipeline/tertiary_output/{source_sample_id}/03_acmg/{source_sample_id}.snv_indel.acmg.tsv`，NGS-UI 再逐筆複製到 `$NGS_UI_HOME/tertiary_output/{Sample ID}/snv_indel.annotated.tsv` 並寫 `pipeline_source.json`；只有設定 `NGS_UI_TERTIARY_LEGACY_STAGING=1` 時才回到舊的 `stage_dragen_for_tertiary.sh` 單樣本 staging 流程。VCF 建立時間來自檔案 `mtime` Unix timestamp，前端固定以台北時區（UTC+8）顯示。工具列的 Extra VEP 僅顯示 checkbox，與不換行的 `↻ 更新索引`、`三級分析清單` 使用固定高度；tooltip 只描述補 SpliceAI 註解，pipeline 仍照常處理既有 annotation。Extra VEP 與 GeneBe 一樣先送 `GNOMAD_G_AF ≤ 0.01` 或 AF 缺值的候選點，再限於 `NGS_UI_CDS_CANDIDATE_BED`（預設 `$HOME/NGS_UI/biotools/cds_combined.bed`，MANE Select CDS±10bp + chrM / RefSeq supplement）內執行，最後把結果 merge 回完整 TSV。執行進度以 step-based 進度條顯示，Nextflow 會從 stdout 追蹤 PREPARE_VCF / PREPARE_VCF_DRAGEN、VEP、Pangolin、CSQ parse、ACMG 等內部 process，stop-gaps 會追蹤 ClinVar fallback / filter / GeneBe / extra VEP / AnnotSV / review TSV；詳細 log 預設收合。worker-owned log 行與 stop-gaps 子步驟帶 ISO timestamp，子程序完成時另記錄 elapsed seconds；`state.json` 的 `step_history` 可供後續依實測耗時調整百分比。「三級分析清單」會合併 `/home/pipeline/tertiary_output/` 與 NGS-UI job state，因此失敗且尚未建立 output 的 sample 仍會顯示 log；刪除時會一併刪除 pipeline output、`$NGS_UI_HOME/tertiary_output/{sample}/` 與該 sample 的 NGS-UI job log 目錄，執行中的 sample 不可刪除。新 job 狀態寫在 `data/jobs/tertiary/`；舊版 `data/jobs/dragen/` 紀錄仍可讀取。`run_stopgaps.sh` 不再建立 `.raw` snapshot；GeneBe VCF 會帶 contig header，AnnotSV 成功執行時只保留摘要 log。

服務啟動時會以 daemon thread 在背景預熱 HPO、phenotype gene map、OMIM 與 mito ClinVar cache，避免解析大型 `phenotype_to_genes.txt` 期間擋住 HTTP port。完整 SNV TSV 的 gene-search cache 也由單一 daemon queue 預熱；重啟服務時不會等待尚未完成的預熱工作。刪除個案時只 invalidate 該 sample 的 SNV 與個案摘要 cache，不會同步重掃全部清單。

三級分析另保存 `source_sample_id`（原始 sequencing sample ID）到 job state 與 `pipeline_source.json`。IGV 先以目前個案 ID 找 BAM；自訂輸出 ID 找不到時，會用 sidecar 回查原始 BAM ID。舊個案缺少 sidecar 欄位時，僅在移除已知 suffix 後得到唯一 BAM 命中才採用 fallback。

### 臨床表徵工具 `/phenotype/`

獨立頁面（內網信任、無需登入）：搜尋 HPO term、套用 / 自訂 gene panel、把結果存成 token 之後在「載入新個案」帶入。Token 限 `[A-Za-z0-9_-]{1,32}`，內容 ≤64KB，panel ≤5000 個基因；自訂 panel 的基因 symbol 不會被轉大寫（`C7orf50` 保持原樣）。

主畫面的 Patient phenotype card 也提供 `WES-I / WES-II / WGS / Other panel` tabs，預設展開 `Other panel`；固定 panel 可直接點 chip 選取，其他 panel 仍可用 typeahead 搜尋。

固定 panel 的來源是 `fixed_panel/WES-I.xlsx`、`fixed_panel/WES-II.xlsx` 與 `fixed_panel/other_panel/`。更新 Excel 後執行 `PYTHONPATH=backend python scripts/import_fixed_panels.py`，會同步重建三個入口共用的 `phenotype_data/fixed_panels/index.json` 與 `phenotype_data/gene_panels/*.txt`。WES Excel 只會匯入 `gene panel list` 標記列起始的基因區塊，避免把疾病名或資料來源列誤算成基因。

---

## 5. 新增一個個案 / 跑 mitochondrial annotation

次級 pipeline 通常會直接把 `snv_indel.annotated.tsv` 等放進 `tertiary_output/{LIS_ID}/`。Mitochondria 的 TSV 用本 repo 的 script 從 GATK Mutect2 `--mitochondria-mode` 的 VCF 產生（純 Python，只需 VCF + 本地 MITOMAP 表，不需 VEP/bcftools）：

```bash
# MITOMAP_DIR 預設 ${REF_DIR}/tertiary/mitomap，內含
#   mitomap_mutations_coding_control.tsv  與  mitomap_mutations_rna.tsv
scripts/annotate_mito_vcf.sh \
  --in   /path/to/{LIS_ID}.mito.vcf.gz \
  --sample {LIS_ID} \
  --outdir tertiary_output/{LIS_ID}/
# → tertiary_output/{LIS_ID}/mito.annotated.tsv
```

批次：

```bash
for s in 26WE0043 26WE0044 26WE0045 26WE0046 26WE0047 26WE0048 26WE0074; do
  scripts/annotate_mito_vcf.sh --in vcf/$s.mito.vcf.gz --sample $s \
    --outdir tertiary_output/$s/
done
```

> 注意：MITOMAP 那兩個 TSV 是 **Latin-1** 編碼（含 0xa0 byte），不是 UTF-8。
> 其他轉檔/遷移 script 見 `scripts/`（`convert_anno_combined_to_tertiary_tsv.py`、`migrate_layout.sh`、`migrate_to_versioned_layout.py` 等）與 `CLAUDE.md`。

---

## 6. 帳號管理

```bash
PYTHONPATH=backend python3 -m app create-user [USERNAME]   # 互動輸入密碼（bcrypt，上限 72 bytes）
PYTHONPATH=backend python3 -m app list-users
```

Session cookie 8 小時、`SameSite=Lax`、`https_only=False`（內網可能還沒 HTTPS）。

---

## 7. 注意事項

- **不要把病人資料 / 大檔 commit 進 git**：`.gitignore` 已排除 `tertiary_output/`、`data/`、`patient_list/`、`phenotype_data/`、`_index.json`。
- 首頁歡迎文字與版本紀錄放在 `frontend/VERSION.md`，前端啟動時會讀取並顯示在尚未載入個案的首頁。之後若有影響判讀流程、報告輸出、資料載入或主要工具入口的更新，需評估是否同步更新這份版本紀錄。
- `phenotype_data/` 必須放對位置（`NGS_UI_PHENO_DATA_DIR` 沒有 fallback）；正式環境記得 `mv NGS_UI/NGS-UI/phenotype_data NGS_UI/phenotype_data`。
- EMR 相關功能預設停用，需設 `NGS_UI_EMR_CLIENT_ID` 才會啟用，且只在內網可達。
- `/api/phenotype-tool/*` 與 `/api/healthz` 是刻意公開無認證；`/api/patient_list` 與其餘 `/api/*` 需登入。
- 大型 JSON response 會在瀏覽器支援時自動 gzip；SNV parse + phenotype / Exomiser / LIRICAL / OMIM join 使用有上限的 process-local LRU cache，輸入 TSV 或 sidecar 更新後自動失效。
