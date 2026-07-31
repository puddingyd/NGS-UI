# NGS 分析平台 (NGS-UI)

成大醫院基因醫學部的 NGS 三級分析判讀工具。次級 pipeline（Nextflow，跑在另一台 compute cluster）產出 per-sample 的註解 TSV，本平台讓 reviewer 載入個案、檢視 SNV/Indel + CNV/SV + Mitochondria + STR + PGx 變異、標記 causative / candidate / other、整理 ACMG SF / 中風 / Carrier screening secondary findings、撰寫判讀意見，並匯出診斷報告或健檢報告 (docx)。另附一個獨立的「臨床表徵輸入 (HPO / gene panel)」工具掛在 `/phenotype/`。

- 後端：FastAPI + uvicorn（Python 3.10+）
- 前端：原生 HTML/CSS/JS，**無 build step**（直接 serve `frontend/`）
- 背景工作：Redis + RQ（Exomiser / LIRICAL 重跑）
- 帳號：SQLite (`data/users.db`) + bcrypt
- 部署：內網 `192.168.84.91:8765`，systemd unit `ngs-ui`

> 開發者 / 接手者請另看 `CLAUDE.md`（架構、資料流、各模組細節、踩雷紀錄）與 `docs/`。

---

## 1. 目錄佈局

程式與一般 runtime 路徑由 `NGS_UI_HOME` 推導；三級輸出統一放在 datalake：

```
NGS_UI/                    ← NGS_UI_HOME
├── NGS-UI/                ← 這個 git checkout（REPO_ROOT）
├── biotools/              ← Exomiser / LIRICAL CLI + data
│   └── litvar2/           ← 官方 bulk JSON、slim SQLite、版本/更新狀態
├── vcf/                   ← per-sample VCF
├── tertiary_output/       ← 舊 UI sample root；只在舊個案遷移期間讀取
├── patient_phenotype/     ← {LIS_ID}_{MRN}_phenotype.txt（自動帶入 HPO）
│                             + {MRN}_clinical_presentation.txt
├── patient_list/          ← 上傳的「未完成報告清單」xlsx + 衍生 roster.json
├── phenotype_data/        ← git-tracked fixed/custom panels; large HPO refs live in NGS_UI_HOME
│   └── gene_disease/       ← optional GenCC / ClinGen / MONDO raw + SQLite disease補充
├── ngs_panel_deadzone/    ← expanded reportable gene list、HGNC alias map、dead-zone tables
├── OMIM/OMIM.xlsx         ← OMIM 疾病註解表（缺檔則 Disease 欄留空）
└── data/                  ← server runtime state（users.db, manual_acmg.sqlite, jobs/, ...）

/home/datalake_Intermediate/pipeline/tertiary_output/
├── _index.json, _case_table.json
└── {LIS_ID}/
    ├── 03_acmg/{source}.snv_indel.acmg.tsv   ← 唯一完整 SNV source of truth
    │           {source}.annotation_versions.json ← ClinVar 等資料庫版本 sidecar
    ├── 04_mito/, 05_str/, 06_cnv_sv/, 07_pgx/ ← UI 直接讀，不再複製
    └── 08_postprocessing/
        ├── {LIS_ID}.layout.json, {LIS_ID}.pipeline_source.json
        ├── {LIS_ID}.snv_annotations.sqlite   ← post-processing 稀疏欄位 overlay
        ├── {LIS_ID}.snv_indel.review.tsv, {LIS_ID}.snv_gene_index.sqlite
        ├── {LIS_ID}.sample_metadata.json, {LIS_ID}.case_summary.json
        ├── analyses/{ver}/{LIS_ID}.analysis.json、{LIS_ID}.pheno_score.tsv、...
        └── {LIS_ID}.ploidy.vcf.gz、{LIS_ID}.ploidy_qc.txt、{LIS_ID}.vcf_from_tsv.vcf.gz 等
```

開發 checkout 不需要這整棵樹：`NGS_UI_HOME` 未設時，只有典型正式部署路徑 `NGS_UI/NGS-UI` 會把 parent 當作資料根；一般 standalone checkout（例如桌面上的 `NGS-UI`）會 fallback 成 repo 自己，所有路徑都落在 repo 內。

`08_postprocessing` layout v3 的 sample-owned 檔案一律加 `{LIS_ID}.` 前綴，讀取順序固定為 prefixed exact name → 舊版 unprefixed exact name，不使用 wildcard。舊 layout v2 的 `layout.json` 與未加前綴檔仍可直接讀寫；完整重跑或 migration 會非破壞性複製成 prefixed 檔，最後才寫 `{LIS_ID}.layout.json`，不會自動刪除舊檔。

---

## 2. 安裝與啟動

```bash
# 1. 取得程式
git clone <repo> NGS_UI/NGS-UI && cd NGS_UI/NGS-UI

# 2. Python 套件
python3 -m pip install -r backend/requirements.txt

# 3. Redis（背景工作佇列需要；沒有的話 Exomiser/LIRICAL 重跑會無法 enqueue）
#    sudo apt install redis-server && sudo systemctl enable --now redis

# 4. （正式環境）把 runtime state 移出 repo；舊三級個案另用
#    scripts/migrate_tertiary_output_layout.py 做 canary 後批次遷移。
#    大型 phenotype reference 仍放在 NGS_UI_HOME/phenotype_data。

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
MemoryMax=16G
KillMode=process
```

```
sudo systemctl daemon-reload && sudo systemctl restart ngs-ui
```

LitVar2 每月一號背景更新另有 systemd timer（正式機路徑已寫在 unit；首次可先從三級 modal 手動更新）：

```bash
scripts/install_litvar2_timer.sh
systemctl list-timers ngs-ui-litvar2-update.timer
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
| `NGS_UI_TERTIARY_ROOT` | `/home/datalake_Intermediate/pipeline/tertiary_output` | 新三級分析統一 root（Nextflow 00-07 + UI `08_postprocessing`） |
| `NGS_UI_PIPELINE_OUT_ROOT` | 同上 | `NGS_UI_TERTIARY_ROOT` 的相容 alias |
| `NGS_UI_LEGACY_TERTIARY_OUTPUT_ROOT` | `$NGS_UI_HOME/tertiary_output` | 遷移前 UI sample root；`TERTIARY_OUTPUT_ROOT` 仍可作舊 alias |
| `NGS_UI_LEGACY_PIPELINE_OUT_ROOT` | `/home/pipeline/tertiary_output` | 遷移前 Nextflow output root（read-only fallback） |
| `NGS_UI_DATA_ROOT` | `$NGS_UI_HOME/data` | users.db, jobs/ |
| `NGS_UI_MANUAL_ACMG_DB` | `$NGS_UI_HOME/data/manual_acmg.sqlite` | hg38 + normalized chr-pos-ref-alt 全域 manual ACMG revisions/current 與 active Causative/Other observations |
| `NGS_UI_SNV_CACHE_MAX` | `2` | uvicorn process 內保留的 SNV review payload LRU 筆數；設 `0` 可關閉 |
| `NGS_UI_SNV_CACHE_MAX_RAW_MB` | `100` | 只有小於此大小的完整 `03_acmg` raw TSV fallback 解析結果才可進 SNV cache；gene search 缺 index 時串流掃描，不會整份物件化 |
| `NGS_UI_TERTIARY_NF_WORK_ROOT` | `$NGS_UI_HOME/nf_work` | 三級分析 Nextflow cache/work root；DRAGEN 與 NCKUH 分別共用 `contexts/shared-dragen/`、`contexts/shared-nckuh/` 的 launch/session lineage；新 lineage 的 work 預設在其 `work/`，若承接舊 session 則持續使用該 session 原 work root |
| `NGS_UI_VCF_DIR` | `$NGS_UI_HOME/vcf` | per-sample VCF |
| `NGS_UI_PHENOTYPE_DIR` | `$NGS_UI_HOME/patient_phenotype` | `{LIS}_{MRN}_phenotype.txt`、`{MRN}_clinical_presentation.txt` |
| `NGS_UI_PATIENT_LIST_DIR` | `$NGS_UI_HOME/patient_list` | 上傳清單 + roster.json |
| `NGS_UI_PHENO_DATA_DIR` | `$NGS_UI_HOME/phenotype_data` | hp.obo / phenotype_to_genes 等大型 HPO reference |
| `NGS_UI_GENE_DISEASE_DB` | `$NGS_UI_HOME/phenotype_data/gene_disease/gene_disease.sqlite` | optional GenCC / ClinGen / MONDO gene-disease SQLite index；缺檔時靜默停用 |
| `NGS_UI_GENE_DISEASE_RAW_DIR` | `$NGS_UI_HOME/phenotype_data/gene_disease/raw` | GenCC / ClinGen / MONDO raw downloads |
| `NGS_UI_GENE_DISEASE_TSV` | `$NGS_UI_HOME/phenotype_data/gene_disease/gene_disease.tsv` | audit / legacy fallback TSV；runtime 優先讀 SQLite |
| `NGS_UI_GENE_PANELS_DIR` | `REPO_ROOT/phenotype_data/gene_panels` | git-tracked fixed panel gene lists |
| `NGS_UI_FIXED_PANELS_DIR` | `REPO_ROOT/phenotype_data/fixed_panels` | git-tracked fixed panel UI index |
| `NGS_UI_CUSTOM_GENE_PANELS_DIR` | `REPO_ROOT/phenotype_data/custom_panels` | git-tracked custom panel gene lists |
| `NGS_UI_PANEL_DEADZONE_DIR` | `REPO_ROOT/ngs_panel_deadzone` | expanded reportable gene list、HGNC alias map、WES/WGS dead-zone tables |
| `NGS_UI_OMIM_XLSX` | `$NGS_UI_HOME/OMIM/OMIM.xlsx` | OMIM 疾病註解（缺檔 = 停用） |
| `NGS_UI_BIOTOOLS_DIR` | `$NGS_UI_HOME/biotools` | Exomiser / LIRICAL |
| `NGS_UI_GENEBE_DB` | `$NGS_UI_HOME/biotools/genebe/genebe_hg38.tsv.gz` | GeneBe 主 DB；永遠優先於 API cache/live API |
| `NGS_UI_GENEBE_API_CACHE` | 主 DB 同目錄 `genebe_api_cache.sqlite` | live API 成功結果與 30-day no-result cache |
| `NGS_UI_GENEBE_API_PENDING_DIR` | 主 DB 同目錄 `api_pending/` | 可匯入主 DB 的精確 7 欄去重 TSV + JSON sidecar |
| `NGS_UI_GENEBE_API_ENABLED` | `true` | 設 `false` 關閉 live fallback；主 DB 查詢不受影響 |
| `GENEBE_USER` / `GENEBE_API_KEY` / `GENEBE_SIF` | 空／`$NGS_UI_HOME/biotools/genebe.sif` | optional GeneBe live API 帳密與預先部署的 pygenebe Apptainer image；可放 mode 0600、git-ignored 的 `$NGS_UI_HOME/secrets.env` |
| `NGS_UI_LITVAR2_DIR` | `$NGS_UI_HOME/biotools/litvar2` | LitVar2 bulk、SQLite、manifest、更新狀態與 log |
| `NGS_UI_LITVAR2_BULK_PATH` / `NGS_UI_LITVAR2_DB` | 上述目錄內 `litvar2_variants.json.gz` / `litvar2.sqlite` | 可個別覆寫 bulk 與 slim index 路徑 |
| `NGS_UI_LITVAR2_BULK_URL` | NCBI LitVar 官方 bulk URL | 每月／手動更新來源；三級分析本身不呼叫 LitVar2 API |
| `phenotype_reference/` | `REPO_ROOT/phenotype_reference` | git-tracked Exomiser / LIRICAL rerun YAML templates |
| `NGS_UI_INHOUSE_BAM_ROOT` | `/home/datalake_Intermediate/pipeline/nextflow_output` | IGV 搜尋 in-house / Nextflow BAM 的根目錄；可用 `:` 分隔多個 root；舊 `NGS_UI_BAM_ROOT` 仍作為 fallback |
| `NGS_UI_DRAGEN_BAM_ROOT` | `/home/datalake_Raw/Novaseq:/home/datalake_Intermediate/n102968` | IGV 搜尋 DRAGEN raw / test BAM 的根目錄；尋找 `<run>/bam/{sample}.bam`，排除 `{sample}.repeats.bam` |
| `NGS_UI_IGV_REF_DIR` | `/home/pipeline/reference/hg38` | IGV 本機 hg38 FASTA + `.fai` |
| `NGS_UI_SECONDARY_WES_FASTQ_ROOTS` | `/home/datalake_Raw/NextSeq2000:/home/datalake_Raw/Other/Reanalysis:/datalake_Raw/Other/Reanalysis` | 二級分析 FASTQ 搜尋 WES roots；可用 `:` 分隔 |
| `NGS_UI_SECONDARY_WGS_FASTQ_ROOTS` | `/home/datalake_Raw/Novaseq` | 二級分析 FASTQ 搜尋 WGS roots；WGS 一律使用 lane FASTQ，不使用 merged FASTQ |
| `NGS_UI_SECONDARY_OUTPUT_ROOT` | `/home/datalake_Intermediate/pipeline/nextflow_output` | 二級分析 samplesheet 寫入位置 |
| `NGS_UI_SECONDARY_SAMPLESHEET_STAGING_ROOT` | `/home/datalake_Intermediate/pipeline/nextflow_samplesheet_staging` | NGS-UI 建立 samplesheet 的 staging 位置；DGX2 指令會複製到 output dir |
| `NGS_UI_SECONDARY_DGX_RAW_ROOT` | `/datalake_Raw/datalake_Raw` | DGX2 上對應 DGM `/home/datalake_Raw` 的 raw data root |
| `NGS_UI_SECONDARY_DGX_OUTPUT_ROOT` | `/datalake_Intermediate/pipeline/nextflow_output` | 產生 DGX2 指令時使用的 output root |
| `NGS_UI_SECONDARY_DGX_WORK_ROOT` | `/raid/DGM/work` | 二級分析 DGX2 Nextflow work/cache root；二級 modal 依此產生需貼到 DGX2 執行的清理指令 |
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
- SNV 載入會把畫面中的 variant IDs 一次送進 `manual_current` / `observations` indexed bulk query，再疊加 per-sample snapshot；不會逐張卡片查 SQLite。Observed badge 只先取 count，實際 case list 點擊後才 lazy 查詢。
- 「個案清單」開啟時呼叫 `/api/samples/case-summary`，後端直接讀統一三級 root 的 `_case_table.json` 總表；若總表缺失或 sample 集合不一致才重建。
- 個案摘要仍會寫入每個 sample 目錄的 `case_summary.json`，依 metadata、TSV/index、Mito TSV 與 OMIM 簽章自動失效；已標記 Mito 會先用 `mito.annotated.tsv` adapter 抓 mito 卡片資料，已標記 SNV 再用 `snv_gene_index.sqlite` 依 variant id 查找，CNV/SV 只讀目標 id。註冊個案、metadata/report/phenotype/analysis 變更與刪除個案時會同步更新 `_case_table.json` 的單一 row；個案清單的刪除只清除註冊/報告狀態並讓 sample 回到未登錄清單，避免開啟個案清單時重掃大型 DRAGEN TSV。
2. **載入新個案** — 點「載入新個案」：
   - LIS_ID 下拉會列出 `08_postprocessing/{LIS_ID}.layout.json`（v3）或舊 `layout.json`（v2）已啟用、但尚未登錄的三級個案；尚在 post-processing 或失敗而沒有 marker 的新目錄不會提早出現。未遷移舊個案仍從 legacy UI root 列出；
   - 若先用「上傳個案清單」匯入過「未完成報告清單」xlsx，MRN / 姓名 / Test type 會自動帶入（來自 `patient_list/roster.json`）。檔案選擇器可一次多選 xlsx，前端會依序上傳並只顯示當次各檔案的成功／失敗與解析統計；歷次紀錄預設收合，按小三角形後才載入。匯入器會掃描所有工作表，標題列可位於任一列／欄；至少需有 `檢體編號`，並同時支援 `檢驗名稱` 與 `檢驗項目` 欄名。三級分析輸出若使用 `{LIS_ID}-dragen` / `{LIS_ID}-nckuh` / `{LIS_ID}-inhouse` 這類 UI 後綴，會先保留後綴作為 sample ID，再回查未加後綴的 roster LIS_ID；
   - Test type 提供 `WES`、`WGS`、`TITAN-WGS`；LIS/sample ID 以兩位數年份加 `T` 開頭（例如 `25T...`、`26T...`、`27T...`）時會自動歸為 `TITAN-WGS`。其他來源若為 DRAGEN，或 in-house 來源 VCF 大於 100 MB，仍預設為 `WGS`；`TITAN-WGS` 的分析門檻、dead zone 與報告方法文字皆沿用 WGS 規則。TITAN-WGS 載入時會預設隱藏診斷分析，可從基本資料下方的小三角切換；Secondary findings、PGx/PharmCAT、健檢報告及儲存仍顯示，而且報告區、分析區內的所有 Secondary findings 子區域都預設展開。WES/WGS 不顯示此切換且維持既有開合狀態；
   - Clinical presentation 可在檢體編號 / 病歷號下方輸入，會依病歷號 debounce 自動儲存為 `patient_phenotype/{MRN}_clinical_presentation.txt`，載入新個案與主畫面 Clinical presentation 會自動帶入；主畫面 reviewer 修改後也會同步寫回此檔，供 `/phenotype/` 後續載入。若沒有 MRN，才 fallback 使用 LIS_ID 暫存。
   - HPO / gene panel 可在這裡選；gene panel 與主畫面同樣使用 `WES-I / WES-II / WGS / Other panel` tabs，預設展開 `Other panel`，固定 panel chip 與搜尋下拉都會顯示基因數量；HPO / panel 下拉可用上下鍵選取並以 Enter 加入，避免 Enter 誤送出載入個案；若存在 `patient_phenotype/{LIS}_{MRN}_phenotype.txt` 會自動讀入；
   - 登錄新個案的未登錄個案欄位可直接輸入 LIS ID、source sample、姓名或 MRN 搜尋，也可從下拉清單選擇；清單在前端快取一天，需要看到最新 pipeline output 時可按「更新清單」手動重抓。登錄完成會同步切換主畫面與上方個案搜尋框。登錄時不再同步掃完整 TSV 產生 `vcf_from_tsv.vcf.gz`；至少有一個 HPO term 時才會排入 Exomiser/LIRICAL，若 VCF 尚不存在，背景 job 開始前會自動建立或刷新。
   - HPO/panel 的 in-panel 狀態來自 `pheno_score.tsv` 動態補值，不再寫回大型 `snv_indel.annotated.tsv` 的 `IN_PANEL` 欄。
   - 旁邊的「個案清單」可查看已載入個案、依 `WES` / `WGS` / `TITAN-WGS` 篩選並全文搜尋；表格會摘要目前 active analysis 的 HPO/panel、causative / other SNV、CNV、SV、Mito、疾病與 comment。表格內「刪除」會移除 prefixed 與 legacy metadata/case summary 及 `analyses/` 註冊/報告狀態，保留 00-07 pipeline output 與其他衍生檔；完成後該 sample 會回到「載入新個案」清單。
   - 主畫面搜尋框與個案清單上方有可複選的 `WES` / `WGS` / `TITAN-WGS` 圓形 filter；每類旁邊的 `only` 可一次只保留該類、反選另外兩類。
3. **看變異卡片** — 個案載入後先顯示 SNV/Indel（分段載入），CNV/SV、Mitochondria、STR 與 PGx 在背景載完後補上：
   - 平台剛開啟讀取索引、個案核心資料載入與新個案登錄期間都會顯示不可誤關閉的「資料載入中」遮罩，避免重複點擊。
   - 載入或重新載入個案時會重設 sample-scoped UI 狀態：SNV/CNV/Mito/STR 分頁回到預設頁籤，主畫面 gene search 輸入與 gene search modal 舊結果會清空。
   - SNV/Indel tier：`1A / 1B / 1C / 2`（互斥）。`1C — Predicted suspect` 保留 ACMG points ≥4，並納入 Core predictors（P-KNN LLR ≥1、AlphaMissense、BayesDel_noAF、Pangolin）及 Extra-VEP predictors（REVEL、SpliceAI）的規劃門檻；predictor 只用於收入 reviewer 候選，不顯示額外 trigger badge，也不修改 ACMG points/class。原本 ClinVar P/LP 0★/conflicting 專屬 tier 已移除，未符合其他 1C 條件者歸入 `2 — Other`。
   - SNV/Indel 卡片的 Score、ClinVar 與各個 in-silico predictor 旁有可 hover/focus 的 `ⓘ`：Score 說明 variant/phenotype 計分公式；ClinVar 顯示三級輸出的 release date；predictor 顯示本位點 evidence、PP3/BP4 calibrated threshold 或未校準警語，Reference URL 可直接點開。ACMG 改成可點擊的「classification (points) + manual/GeneBe/in-house source」摘要，不再顯示自由文字 criteria 輸入格。顯示順序固定為 P-KNN、AlphaMissense、Pangolin、REVEL、SpliceAI、ESM1b、VARITY_R、BayesDel、MetaRNN、DANN、PhactBoost、PhyloP、GERP、SIFT、LOFTOOL；有值的前三項直接顯示，其餘收在 More。AlphaMissense、ESM1b、VARITY_R、BayesDel、REVEL、SpliceAI、PhyloP、GERP、SIFT 使用文獻/ClinGen SVI threshold；P-KNN 依三級輸出的 `PKNN_EVIDENCE` 套色，註解列出 LLR 的 ±1/±2/±4 門檻。沒有可參考 cutoff 的 DANN 與非 score 型的 LOFTEE 不套色。
   - ClinVar/ACMG 下方顯示 `LitVar2 (版本日期)`：標題維持純文字，旁邊的 `SquareArrowOutUpRight` 外部連結圖示可開啟該 variant 的 LitVar2 結果頁；前五個 PMID 分別直連 PubMed，超過五篇顯示 `and N others` 並連回 LitVar2。只讀三級 post-processing 從本地 bulk SQLite 寫入的結果，先 exact rsID、再 gene + HGVS；不在開卡片時查 API，也不影響 tier、ACMG、排序或報告。舊個案重跑三級後才會補上。
   - ClinVar release date 不從檔案 mtime、今天日期或報告常數猜測。三級 pipeline 應在 `03_acmg/{source}.annotation_versions.json` 寫入 `databases.clinvar.release_date`；舊資料也接受 `*.tsv.meta.json`、同目錄 `annotation_versions.json` 或 `08_postprocessing/pipeline_source.json.annotation_versions`。缺 sidecar 時 UI 明確顯示「三級輸出未提供」。
   - Patient phenotype 與 Comment 之間的 Dead zone 卡片會依目前 HPO + panel gene set 顯示 cohort-level dead exons；WES 用 20X，WGS（含 in-house / DRAGEN）用 DRAGEN 10X。主畫面先依 CDS dead percentage 分成 70-100%、50-70%、30-50%、<30% 四個區間，區間內再依 gene 的 pheno score 由高到低排序，最後用 CDS percentage 與 gene name 當 tie-breaker；預設只顯示 CDS ≥50% 的列，可用標題列或區塊底部的小三角形展開全部。Dead zone 列改用淡背景警示色：≥70% rose、50-70% orange、30-50% amber、<30% yellow。
   - SNV/Indel 主畫面直接顯示 `snv_indel.review.tsv` 內的候選點，不再另套 gnomAD AF filter；顯示 filter 預設啟用「疾病相關」（內部 `disease_associated`，優先使用 `ngs_panel_deadzone/panel/panel_loose_plus_clinical.hgnc_canonical.txt`，缺檔時 fallback 到 `panel_loose.hgnc_canonical.txt`）與「臨床相關」（內部 `in_panel`）。兩者後方數字會先套用其餘啟用中的 SNV 顯示條件，再計算目前候選點。其後的 `AF_NCKUH ≥ 0.05 & AC ≥ 50` 預設不勾選：SNV/Indel 與 mtDNA 同時達到兩項門檻時預設隱藏，勾選後才顯示；ClinVar non-conflicting P/LP、有效 ACMG P/LP、reviewer `1/2/C` 會 rescue，mtDNA 另 rescue MITOMAP pathogenic/reported 與手動 P/LP。Homozygous／hemizygous 與同一 recessive gene 第二個 candidate 不 rescue。`VAF < 0.2 / zygosity=ref` 預設不勾選，主畫面預設隱藏低 VAF 或 ref genotype 點位，勾選後才顯示；`impact=MODIFIER` 預設不顯示，可手動勾選展開，但 ClinVar P/LP 點位不受 MODIFIER filter 限制。`IMPACT=LOW` 仍會顯示。Gene search modal 保留自己的 `gnomAD_G_AF < 0.01` filter；CNV/SV 與 gene search modal 不受「疾病相關」filter 限制。
   - `03_acmg/*.snv_indel.acmg.tsv` 保持 immutable；`REF/ALT=*` 與非 primary contig 改在 review/index/gene search/VCF consumer 端排除，不再整份改寫 raw TSV。worker 會接受 v3.1-v3.4 的 65 欄 TSV；v3.5 移除 `MANE_ALL`、新增 `STRAND_BIAS`，仍為 65 欄且同一 genomic variant 可有多列 transcript。
   - 主畫面讀取 `08_postprocessing/{LIS_ID}.snv_indel.review.tsv`（舊檔 fallback `snv_indel.review.tsv`）：SNV/Indel 會排除 chrM/MT，WES 先濾 DP < 20，再保留 candidate BED 內且 `GNOMAD_G_AF < 0.01` 或 AF 缺值的位點，並 rescue ClinVar P/LP。GeneBe、GIAB、院內 AF、MANE/extra-VEP 等 post-processing 差異存於 prefixed `snv_annotations.sqlite`，產生 review、gene search 與 reviewer-marked 補點時才合併；完整來源固定為 `03_acmg/*.snv_indel.acmg.tsv`，不再永久保存第二份數 GB 的 `snv_indel.annotated.tsv`。gene index offsets 也指向 03_acmg raw；缺 index 時仍可 bounded-memory 串流搜尋。
   - v3.5 `STRAND_BIAS` 在既有 badge 列顯示：`PASS` 為綠色，`WARN(FS=...,SOR=...)` 為紅色，`.` 為紅色 `Strand bias MANUAL`；tooltip 說明 SNV（FS>60/SOR>3.0）或 indel（FS>200/SOR>10.0）門檻。舊 TSV 沒有此欄時不顯示 badge。它只作 QC 提示，不改 tier、排序、ACMG、filter 或報告。
   - Structured ACMG modal 完整列出 28 個 ACMG/AMP 2015 criteria、說明、ACMG 2015 來源及可用的 ClinGen guidance/PubMed；每項可啟用／停用並調 Supporting、Moderate、Strong、Very strong（BA1 另有 Stand-alone）。PP5/BP6 保留可用，但顯示 ClinGen 建議停用的警語。上方並列 Manual、GeneBe、In-house 的 classification + criteria 與 Apply；Apply 只改 editor，按「儲存並套用」才寫入。
   - Manual ACMG 全域 match 只用 `hg38 + normalized chr-pos-ref-alt`，不帶 gene/transcript/condition。每次儲存 append 一筆全域 revision（含 users table 內部 id、username snapshot、來源 sample、時間），同時在該 sample 的 `reports.edits[id].manual_acmg` 留最終快照。`PS2/PM3/PM6/PP1/PP4/BS4/BP2/BP5` 標為 case/family-specific：完整保留在來源 sample/revision，但 Global Manual Apply 與其他 sample 的自動 overlay 會排除。畫面／報告有效優先序為 per-sample snapshot → global reusable manual current → GeneBe → in-house；其他 sample 已有自己的 snapshot 時不會被新的 global revision 覆蓋。儲存後後端權威重算 ACMG points、classification、variant/total score、1C eligibility、tier 與各區排序；若點位換 tier 或排序位置，前端會切到新 tier，並把同一張卡片捲回原本視野位置。正式分類仍維持五級，主卡片以和 ClinVar 相同的單行字體／高度顯示分類與 points 色塊，旁邊 source 為不套色小字；0–5 points 的 VUS 顯示細分成 `VUS-low`（0–1）、`VUS-mid`（2–3）、`VUS-high`（4–5），low 向右漸綠、high 向右漸紅。Modal 左右固定為 Pathogenic／Benign criteria 淡紅／淡綠欄，兩欄皆依 ACMG 2015 evidence type 順序分組；上下各有同一組「關閉／儲存並套用」按鈕。
   - 與 manual ACMG 分開的 observation registry 只保存各 sample 目前 status 恰為 `1`（Causative）或 `2`（Other）的 SNV；Candidate、`0`、已取消與已刪除 sample 不保留。卡片在 HGVS 下一行的既有 badge 列顯示其他 sample 的 `Observed (N)`，點擊獨立 modal 才列 sample ID、Causative/Other、當時 reviewer 與更新時間，目前 sample 自己不計數。
   - Read support 固定顯示 `DP · AD · VAF`。WGS/TITAN-WGS 總 DP<10 會標紅；WES DP<20 仍直接排除。另以第一個 ALT 的 AD<10 獨立顯示紅色 `Low ALT support`，建議 IGV、必要時 Sanger；AD 缺值不會被誤判。這些警示不改 tier 或 filter。
   - SNV/Indel 卡片標籤列會依 GIAB genome-stratification 標出困難區 badge。`scripts/annotate_giab_strata.py` 在暫存 working TSV 產生 `GIAB_STRATA`，最後壓入 `snv_annotations.sqlite`；BED 與 manifest 放在 `NGS_UI_GIAB_STRAT_DIR`。純 UI 提示，不影響 tier 或報告。
   - SNV/Indel 卡片的 OMIM 外部連結會優先使用 TSV 內 OMIM link / ID；若該變異沒有 OMIM disease/ID 但有 gene symbol，仍會顯示 OMIM 按鈕並連到該 gene 的 OMIM geneMap 搜尋頁，方便確認 OMIM 是否已有更新。
   - SNV/Indel 卡片的疾病清單仍以最新版 `OMIM.xlsx` 的 `Disease1..5` 作為 rich description / synopsis 來源；若部署了 `NGS_UI_GENE_DISEASE_DB`，後端會把 GenCC / ClinGen / MONDO association 依 gene 合併進 `disease_associations`，用 phenotype MIM、MONDO ID、source disease ID 與 gene-scoped normalized disease name 去重，補上 OMIM 沒列到或額外的非重複疾病。DB 由 `python scripts/update_gene_disease_db.py` 建置：下載 GenCC submissions CSV、ClinGen gene validity CSV、MONDO JSON，建立 `gene_disease.sqlite`、audit `gene_disease.tsv`、`source_manifest.json` 與 `build_report.tsv`；GenCC/ClinGen 只納入 `Limited` 以上，低於 `Limited`、disputed/refuted 等會略過。UI 固定先列 OMIM disease，以稍深灰底顯示，第一筆補充疾病前另有細分隔線。補充列目前只在 UI 顯示，不改 DOCX 的 `report_diseases` 選取；OMIM disease 只有單行、沒有下方 description/synopsis 時，summary 後會以小 `*` 提醒 curator 補描述，星號不寫回資料也不進列印輸出。
   - Secondary findings 從 `snv_indel.review.tsv` 產生，只保留 ACMG SF、中風相關基因（`WGS__神經科__Stroke`）與 Carrier screening；需符合 panel gene、ClinVar P/LP 或 SNV tier `1A/1B/1C`、VAF ≥ 0.2、zygosity 非 ref，因此會納入 LOFTEE HC、ACMG points ≥4、P-KNN LLR ≥1 及與主 SNV 卡片相同的其他 predictor 候選。分析區保留全部候選點位；同一 variant 命中多個 panel 時會依各 panel 基因內容分別顯示在每個分類，但所有卡片共用同一筆 variant review state，勾選／取消、ACMG、Comment、transcript 與疾病選取會同步。只有 ClinVar P/LP 預設 ✓；其他 1B/1C 候選預設未勾，報告區只顯示 ✓ 點位。任一所屬 panel 有明確取消即以取消為準，健檢匯出時才將不同 panel 的相同 variant 合併去重。
   - 「匯出健檢報告」只提供 ACMG 疾病風險基因、中風相關基因、帶因者篩查與藥物基因體學四個選項；ACMG 與 PGx 預設勾選，中風與帶因者預設不勾。所選疾病 panel 的點位會先合併去重，再統一分成「第一類：與疾病風險相關之致病性或疑似致病性變異位點」與「第二類：符合帶因者狀態之致病性或疑似致病性變異位點」；AD、X-linked、純 AR 同基因兩個以上變異及純 AR homozygous 變異列第一類，純 AR 單一 heterozygous 變異列第二類。§五「本次檢測基因包括」的 panel 標題不再帶「第一類／第二類」等前綴。附錄另起新頁，「附錄」置中且下方空一行，再依序列出變異位點參考資料與完整用藥建議。PGx 主文與附錄仍使用 `pharmcat.report.json` 的 CPIC/FDA 建議，固定列 21 個 CPIC Level A genes，不含 CYP3A4 與 IFNL3。
   - 報告區 Causative / Other / Candidate 卡片依 SNV/Indel、CNV/SV、Mitochondria 排列；Secondary findings 報告區預設只展開 ACMG SF，其餘收合。分析區 Secondary findings 預設只展開 ACMG SF，PGx / PharmCAT 預設展開，其餘 secondary panels 預設收合。
   - GeneBe ACMG 第二意見（`GENEBE_ACMG_*` 欄）採 **本機主 DB → 本機 API cache → live API**。post-processing 仍先對整張 TSV 查 `NGS_UI_GENEBE_DB`（預設 `biotools/genebe/genebe_hg38.tsv.gz`）的 `acmg_score/acmg_criteria`；主 DB 的 lazy `genebe_hg38.sqlite` 會依來源檔 signature 自動重建，SQLite 失敗才 fallback single-pass streaming。只有主 DB/API cache 都 miss 且符合 review TSV 條件的點才送 API：primary/non-mito、WES DP≥20、ClinVar P/LP rescue，否則 `GNOMAD_G_AF` 缺值或 <0.01 且位於 candidate BED；WGS 不設 DP hard floor。API 使用 `GENEBE_USER`、`GENEBE_API_KEY` 與預先部署的 `GENEBE_SIF`；缺設定、timeout 或 batch 失敗只警告，三級分析仍以主 DB／pipeline ACMG 完成。成功結果寫到 `NGS_UI_GENEBE_API_CACHE`，並在 `NGS_UI_GENEBE_API_PENDING_DIR` 產生可匯入正式 DB 的精確 7 欄去重 TSV（`#chr,pos,ref,alt,acmg_classification,acmg_score,acmg_criteria`）及 JSON sidecar；明確 no-result 預設 negative-cache 30 天。完整 working TSV 仍只是既有暫存檔，overlay 建好後由 worker `finally` 刪除，不新增永久大型 SNV 副本。主 DB 部署/更新用 `scripts/deploy_genebe_db.sh`，格式要求見 `docs/genebe_db_requirements.md`。
   - 需要跨批次整併 API success cache 時執行 `python scripts/export_genebe_api_cache.py --out /path/genebe_api_rows.tsv`；輸出仍是相同 7 欄、全域 variant-key 去重及自然座標排序，另附 `.tsv.json` sidecar。正式主 DB 的 merge/publish 維持人工、原子部署流程。
   - SNV tier 只在點開時建立該 tier 的卡片 DOM，避免一次 render 全部卡片。
   - SNV/Indel 與 CNV/SV gene 搜尋支援多個基因，以 `,` 或 `、` 分隔；SNV 搜尋由 `/api/samples/{id}/snv-search` 查完整原始 TSV，不受 review TSV 限制。三級分析結尾會預建 `snv_gene_index.sqlite`（gene → raw TSV byte offsets），所以 WGS gene search 不需在載入個案時掃 1–2GB raw TSV；舊樣本若缺 index 才 fallback 到 raw TSV parse。modal 預設勾選 `gnomAD_G_AF < 0.01`，取消後才顯示全部搜尋結果。SNV adapter 會把同一 genomic variant 的多 transcript TSV rows 合併成一張卡，預設顯示 consequence 較嚴重者（同嚴重度時 MANE 優先）；卡片 HGVS 優先用 post-processing 補上的 RefSeq transcript，沒有 RefSeq 時才用 Ensembl。HGVS 旁的小三角可切換 transcript，選擇會寫入 reviewer edits，DOCX、PDF 卡片與個案清單摘要使用同一個 transcript。
   - Variant 狀態用 `1 / 2 / C / 0` 圓形按鈕；`C` 與 `0` 可並存，`1` 或 `2` 會反選其他狀態；再次點擊已選項目可清空狀態。同一個 variant 在分析區、報告區與搜尋 modal 的按鈕會同步上色。
   - Variant 卡片保留 OMIM `Disease1..5` 完整內容；「個案清單」中的疾病摘要才截到第一個可辨識的遺傳模式括號（例如 `(AD)`、`(AR)`）為止。
   - SNV/Indel、Mitochondria 與 CNV/SV 卡片有 `IGV` 按鈕；基本資料的性別下方另有「☑ 已確認」核取方塊與「確認 SRY」按鈕。SRY 會沿用同一個 modal 開啟 hg19/hg38 對應區域，初始只載入當前 sample，仍可手動加入 sibling，coverage data range 固定 `0-100`；勾選確認狀態會寫進 reviewer metadata。一般 modal 標題會顯示 sample、padded locus、variant 註解與原始座標；alignment 預設用 squished 模式，`visibilityWindow=5 Mb`，SNV/Indel 與 Mito track height 固定 300，CNV/SV track height 固定 50。SNV/Indel 與 Mito 顯示前後 100bp，CNV/SV 原則上顯示前後各 20% flanking area；若事件本身 ≤5 Mb，padding 會自動縮小以維持初始視窗 ≤5 Mb，直接載入 BAM coverage，不需先 zoom in。IGV 另以 ROI 標出實際 CNV/SV 區間。CNV/SV modal 會把所有 BAM coverage tracks 放進同一個 autoscale group，讓 y-axis data range 一致，方便和 sibling 比較 deletion/duplication。先確認 primary BAM 與同 batch sibling tracks，再按「載入 IGV」；若 UI sample ID 帶 `-dragen` / `-nckuh` / legacy `-inhouse` / `-WES` / `-WGS`，BAM 查詢可 fallback 到去 suffix 的原始 sample ID。DRAGEN 來源會在 raw Novaseq root 或 `/home/datalake_Intermediate/n102968` 測試 root 找 `<run>/bam/{sample}.bam` 並排除 `{sample}.repeats.bam`；in-house 來源維持找 Nextflow output 的 `02_alignment`。若自動搜尋不適用，可按「其他路徑」從下拉選擇 DRAGEN run 或 in-house batch，列出該 run/batch 可用 BAM 後選 primary；同 batch 加入清單會改用同一資料夾的其他 BAM。BAM range request 與 hg38 FASTA 都由後端在內網 proxy。
   - CNV：`CNV-1A`（Clinical）、`CNV-1B`（Pathogenic）；SV：`SV-2A / SV-2B`。CNV/SV 分析區與報告區皆依 `max_pheno_score + scaled AnnotSV ranking score` 由高到低排序；CNV/SV 卡片的基因表預設只顯示 phenotype 相關基因，按小三角形才展開 phenotype score 為 0 的其餘基因。
   - 同來源、同染色體、同為 deletion 或同為 duplication，且相鄰 gap ≤ `250 kb` 的 CNV/SV 會預設自動整合；copy number 差異不再阻擋視覺整合，原始片段仍可展開查看各自 CN。UI、DOCX 與個案清單改用合併後 parent；parent 會取代最佳原始 segment 的位置，不會掉到 tier 最後。前端會依 sample payload 快取整合結果，避免卡片 render 時反覆重建 parent。
   - Mitochondria：`MITO-1`（ClinVar P/LP 或 MITOMAP confirmed/pathogenic）、`MITO-2`（rare / reported mtDNA variant）、`MITO-3`（other variant）— 若來源 TSV 有 `FILTER` 欄，只列 `FILTER=PASS`；v3.2 `04_mito/{sample}.mito.tsv` 沒有 `FILTER` 時不顯示 Filter 欄。MITOMAP 欄位只在後端分類使用，卡片不顯示 MITOMAP；Disease 改用 `CLINVAR_DN` 依 `&` 拆成 checkbox，勾選者才進 DOCX。卡片第三行並列 `AF_gnomAD` 與本院 `AF_nckuh`；後者從 sites VCF 的 indexed chrM slice 顯示 `AF (carrier AC / callable mito AN; homoplasmic/heteroplasmic carrier 數)`，只供辨識院內常見/重現變異，不影響 tier 或 ACMG。分析區共用 `AF_NCKUH ≥ 0.05 & AC ≥ 50` 開關，符合門檻且未被已知致病／reviewer 狀態 rescue 的 mtDNA 卡片預設隱藏。
   - STR：`STR-P`（pathogenic）、`STR-I`（intermediate / borderline）、`STR-N`（normal / no_threshold）— 讀取 pipeline `05_str/{sample}.str.tsv` 複製出的 `str.tsv`，顯示 STRchive locus、疾病、遺傳模式、repeat count、confidence interval、depth 與 benign/intermediate/pathogenic threshold。
   - PGx / PharmCAT：以精簡後的 `pharmcat.report.json` 為主，後端直接共用健檢 DOCX 規則產生固定 CPIC Level A genes 的畫面資料。報告區顯示用藥建議概覽、藥物摘要及基因型/表型表；分析區顯示六類完整概覽、完整 CPIC/FDA 建議，以及每個 gene 的相關藥物、建議與 allele/variant evidence。generic fallback 以較清楚的「其他用藥建議（參考完整建議與最新藥品仿單）」呈現；相關 JSON 藥物與 MT-RNR1 TSV 補充藥物依藥名去重並只進一個主要分類，明確處置優先，其餘才進「不需調整／無明確調整建議」。藥物建議摘要直接沿用概覽 drug group 已選出的 primary action/source，不會針對摘要來源重新分類。基因名稱維持黑字、diplotype/allele 用淡紫底膠囊，只有需注意 phenotype 顯示紅字；逐 gene 建議拿掉重複的基因／表型欄，需調整藥物排在不需調整之前。normal/uncertain 等仍列相關藥物但標示不需調整。PharmCAT object message 不顯示，MT-RNR1 raw `Unknown` allele placeholder 也隱藏；variant evidence 的 star allele 長清單預設收合，欄名 hover tooltip 說明這是可參與定義的清單而不是病人同時帶有全部 allele。pipeline `07_pgx/{sample}.pgx.tsv` 只在 MT-RNR1 的 JSON 為 Unknown/No Result 時補基因型、風險與藥物範圍，不再決定其他 gene 的主畫面範圍。
   - DRAGEN 依原始 hard-filtered VCF basename 從同層精確尋找 `*.ploidy.vcf.gz`；NCKUH 則從 `04_snv_indel/{source}.ensemble.fixed.vcf.gz` 的 sibling `03_alignment_qc/{source}.ploidy.vcf.gz` 尋找。VCF 複製到 `08_postprocessing/{LIS_ID}.ploidy.vcf.gz`；`ploidy_qc.txt` 也複製成 `{LIS_ID}.ploidy_qc.txt` 留待未來使用，但目前 UI/報告不讀。性別欄一律顯示可點擊的 `ploidy VCF: ...`。modal 先顯示 chromosome dosage 結論、estimated karyotype／病歷性別核對與需要複核的染色體，完整染色體與 ALT/FILTER/QUAL/DC/NDC/raw ratio/END 收在折疊區。DRAGEN 以 `<DEL>/<DUP>` 作 dosage call，PASS 與 LowQual 分開顯示；NCKUH 以 `SUSPECT` 作 signal，方向依 NDC 判斷。DRAGEN 未原生提供 RATIO 時，UI 以 `DC / autosomeDepthOfCoverage` 顯示明確標註的 derived ratio。非 XX/XY 或任一核染色體 dosage signal 優先顯示紫色；沒有 signal 時，病歷性別空白/U 或與 XX/XY 不符顯示紅色，M↔XY、F↔XX 顯示綠色。
   - 健檢 DOCX：ACMG SF 維持既有 ASCII 變異表格，表格內 ACMG/AMP 分級保留英文；表後依 AD、AR、X-linked、zygosity、同基因回報變異數及性染色體組成自動產生中文結果意義與建議。性染色體組成優先讀取已複製的 `ploidy.vcf.gz` 之 `estimatedSexKaryotype`，其次使用 EMR 性別，無資料時採中性文字。AR 單一 heterozygous 變異會寫明「本次僅檢出一個符合報告條件之變異，檢測結果符合帶因者狀態」；AR 多變異會寫明兩變異相位尚未確認，建議進行家族成員檢測以釐清是否位於不同等位基因，方能判斷是否符合體染色體隱性疾病之雙等位基因致病型態。建議句使用「遺傳諮詢或門診相關專科」。X-linked 女性單一變異採可能受 X 染色體失活型態及疾病表現範圍影響的保守文字。WGS 健檢報告平均深度為 27X，限制段落列出 CNV 未涵蓋範圍，並另註明 CYP2D6 藥物基因體學判讀會專一納入該基因 CNV；診斷報告的 WGS 深度仍為 30X。PGx 報告與主畫面共用同一套後端 projection；基因型、表型、用藥建議概覽、藥物建議摘要與完整 CPIC/FDA 建議由 PharmCAT JSON 產生，只有 MT-RNR1 在 JSON 無有效結果時可由 TSV 補值。藥物建議摘要與附錄完整用藥建議皆以藥物 A-Z 呈現，摘要不列 FDA-label-only 項目，完整建議同一藥物只列一次且 CPIC/FDA 來源放在建議句尾括號。末尾固定列 21 個 CPIC Level A genes，IFNL3 不會出現在 PGx 主文、附錄或基因清單，標題與參考依據不列 FDA。
   - 健檢 ACMG SF 的附錄參考資料順序跟上方點位顯示順序一致；同一基因多位點在表格合併時，參考資料也連續列出。
4. **標記與判讀** — 在每個變異上標 causative / candidate / other、用 structured modal 編輯 ACMG criteria/strength、寫 comment。SNV/Indel ACMG 優先序為 per-sample manual snapshot → global manual current → GeneBe → pipeline `ACMG_CLASS`；sample reviewer state 存在 `{LIS}/08_postprocessing/{LIS}.sample_metadata.json`（舊個案仍可寫 legacy 名稱），global revision/observation registry 存在 `NGS_UI_MANUAL_ACMG_DB`。開啟個案清單或匯出 DOCX 前會先 flush 尚未完成的自動儲存。

SNV/Indel 卡片的 ESM1b 依 ClinGen SVI 校準區間上色；ESM1b 分數越低越偏致病，多 transcript TSV 會取最低分作為 worst case。
5. **匯出報告** — 「匯出診斷報告」下載 `GET /api/samples/{LIS}/report.docx`；DOCX 依序輸出第一類、第二類、固定建議文字，再集中列出各 variant 的參考資料。CNV/SV Disease 欄會優先覆蓋單基因預設疾病，也可為片段型 CNV/SV 指定發報告疾病；多基因或無 OMIM gene 的片段會把 Disease 直接接在第 1 點位置描述，不另列一點。未涵蓋 OMIM 疾病相關基因時不再列出一般基因清單。SNV/Indel 會帶 TSV 的 `RS_ID`，RS ID 欄會預留尾端空格；ClinVar 欄也預留尾端空格，避免 `conflicting classifications` 與 ACMG 欄黏在一起。SNV transcript header 若有 RefSeq 對到 Ensembl，會顯示 `GENE (ENST...; NM_...)`，沒有 RefSeq 則只顯示 Ensembl；SNV/Indel 與 Mito 核苷酸欄每行最多 13 字元且蛋白質括號另起一行。§五.4 grouped gene list 的 HPO/panel 區塊之間會留空行，且只列 expanded reportable disease-associated genes；WES/WGS DOCX 都不輸出 dead-zone exon 標註，dead-zone 僅保留在主畫面與 phenotype 工具查詢。旁邊的「輸出 PDF」會開啟列印視窗，輸出報告區的 causative / other / candidate 卡片摘要與本次檢測基因清單；PDF 與 DOCX 一樣會先詢問基因清單要按 HPO/panel 分組或全部合併去重，並在按下輸出時才呼叫 `/api/samples/{id}/report-gene-list` 載入基因清單。列印版會略過 comment、More、Secondary findings 與 CNV/SV overlap 明細。

二級分析 modal 由右上角「二級分析」開啟，後端掃描 WES roots（NextSeq2000 與 Reanalysis）與 WGS Novaseq FASTQ，提供 WES/WGS 兩個 typeahead 搜尋框、更新索引、批次清單與「加入同批全部樣本」。WGS 下拉選單依主 sample 聚合顯示，每個 sample 註明 lane 數與 FASTQ 檔案數；只收 `*_S*_L00*_R[12]_001.fastq.gz`，不使用同資料夾內的 merged FASTQ。建立 samplesheet 時，後端會把每個 WGS sample 展開成每 lane 一列並寫入 `lane` 欄，交由 FASTP 分 lane QC、FQ2BAM 合併；WES 的掃描、顯示與 samplesheet 行為不變。按「建立 sample sheet」後會先在 `NGS_UI_SECONDARY_SAMPLESHEET_STAGING_ROOT/<batch_name>/samplesheet.csv` 寫出 CSV；DGM `/home/datalake_Raw/...` 會映射到 DGX2 `NGS_UI_SECONDARY_DGX_RAW_ROOT/...`（預設 `/datalake_Raw/datalake_Raw/...`），已經位於 `/datalake_Raw/...` 的路徑不重寫。產生的 tmux launch block 會由 DGX2 runner 建立 `OUT_DIR`，再把 staged samplesheet 複製到 `OUT_DIR/samplesheet.csv` 後啟動 Nextflow，避免 NGS-UI 建立的 output dir owner 造成權限問題。Nextflow 結束或失敗時會複製 `${LAUNCH_DIR}/.nextflow.log` 到 `${OUT_DIR}/nextflow.log`，並保留 tmux pane 不自動關閉，方便確認最後狀態。不論單一或多個主 sample，產生的指令都固定使用 `dgx_single` profile；WES 另加 `--run_gcnv true`，WES/WGS 都預設加入 `--run_manta`、`--run_expansionhunter`、`--run_automap`，讓二級分析直接產出三個模組的結果。Reanalysis FASTQ 也直接指向原始 `*.R1/R2.clean.fastq.gz`，不再複製或建立 symlink；只改 samplesheet 的 `sample` 欄。因 NGS-UI 跑在 DGM，不能直接清 DGX2；modal 的「顯示 DGX2 清理指令」會依 `NGS_UI_SECONDARY_DGX_WORK_ROOT`（預設 `/raid/DGM/work`）產生 guarded shell script，必須複製到 DGX2 terminal 執行。script 會先阻擋仍有 Nextflow process 的狀況、列出待刪項目並要求確認；「複製」會直接寫入剪貼簿，plain-HTTP intranet 自動使用隱藏文字框 fallback，不再跳出手動複製視窗。清理後既有批次不能再用 `-resume`。

三級分析 modal 的 in-house 與 DRAGEN VCF 各自使用單一 typeahead。PGx 和「更新索引」之間的「更新 LitVar2」會啟動 detached worker，背景下載官方 bulk（壓縮檔約 1.8 GB）並重建本地 SQLite；新資料先驗證再原子切換，失敗繼續使用舊版，modal 會顯示進度、版本與錯誤。每次送出三級分析都會啟動 Nextflow `-resume`，不再由 UI 依正式資料夾中的檔案決定是否略過；DRAGEN 與 NCKUH 各自固定使用一條共享 launch/session/work lineage，因此同一 sample 即使前後放在不同 batch，只要 Nextflow 判定 task 的 input、script、container 與其他 hash 條件相同，就能沿用 cache，條件改變的 task 則自行重跑。同一 pipeline 類型一次只允許一個 Nextflow process 開啟共享 cache；後送 job 會停在 `waiting-nextflow-cache`，DRAGEN 與 NCKUH 彼此不互鎖，Nextflow 結束後的驗證與 post-processing 也可繼續並行。共享 lineage 第一次建立時會從舊 launch contexts 選擇與本次 sample 重疊最多的可用 session 承接；不同 session ID 的舊 cache 無法合併，未承接到的 task 會由 Nextflow 正常重算並逐步累積到共享 lineage。相同 UI/source sample 正在執行時，重複送出仍會回 HTTP 409。

Nextflow 的 `--out_dir` 指向 job-private `.staging/{job_id}/`，不直接覆寫正式個案。結束後先驗證該 source sample 的 `00_prepare`–`06_cnv_sv` 目錄以及 SNV、Mito、STR、CNV、SV 產物；有勾 PGx 時再要求 `07_pgx` 與 PGx/PharmCAT 產物。驗證通過後，SNV post-processing 也在 staging 中完成 GeneBe/extra-VEP/GIAB/院內 AF/MANE/LitVar2、prefixed overlay/review/index、MITOMAP 與 ploidy 衍生檔；LitVar2 僅查與 review TSV 相同篩選範圍的 genomic variants，缺 DB 或查詢失敗為 warning/no-op。hidden working TSV 仍在 `finally` 刪除。最後進入 `promote-output`，以可 rollback 的 rename 切換正式 00–07 與 worker-owned 08 衍生檔，保留當下最新的 reviewer metadata、case summary 與 analyses，並最後才切換 `{LIS_ID}.layout.json`。任何 Nextflow、完整性驗證或 post-processing 失敗都只清理 staging，原本正式結果與已載入個案維持可用；批次切換中途失敗則還原同批已切換的 sample。若取消 PGx，舊 `07_pgx` 會在成功切換時移除，避免顯示上一次的 PGx 結果。三級分析清單會忽略 `.staging` / `.rollback`；完整刪除仍會清掉對應 pipeline/UI tree 與 job log，執行中的 sample 不可刪。

Extra VEP 預設讀取 `NGS_UI_HOME/biotools/dbnsfp/dbNSFP5.3.1a_grch38.gz`（dev 機即 `/home/n102968/NGS_UI/biotools/dbnsfp/dbNSFP5.3.1a_grch38.gz`）補 MetaRNN 與 REVEL；必須有同名 `.tbi`，啟動時會檢查實際 header，仍可用 `NGS_UI_EXTRA_VEP_DBNSFP` 或 `--extra-vep-dbnsfp` 覆寫。SpliceAI 維持讀取 `/home/n102968/NGS_UI/biotools/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz` 與 `spliceai_scores.raw.indel.hg38.vcf.gz`，兩個 score VCF 都存在時才補 `SPLICEAI_MAX`，不使用桌面 `/Volumes` 路徑。

正式 output 仍以 UI sample ID（例如 `VAL-57-dragen` / `VAL-57-nckuh`）管理；Nextflow staging 先依 samplesheet 的原始 `source_sample_id` 接收產物，通過驗證後才 promotion 到 UI sample 目錄。legacy `/home/pipeline/tertiary_output` 僅作過渡讀取，不參與新的 cache／切換判斷。

舊個案請先 dry-run，再挑 1–2 個 canary；工具永遠保留舊資料，`{LIS_ID}.layout.json` 也只在 prefixed overlay/review/index 與 metadata checksum 驗證通過後才寫入：

```bash
# 只看計畫，不寫入
PYTHONPATH=backend python scripts/migrate_tertiary_output_layout.py --sample SAMPLE_ID

# canary 實際遷移，確認 UI / gene search / 報告後再跑第二個
PYTHONPATH=backend python scripts/migrate_tertiary_output_layout.py --sample SAMPLE_ID --apply

# marker-only rollback：立即回到 legacy UI root，不刪新舊任何檔案
PYTHONPATH=backend python scripts/migrate_tertiary_output_layout.py --sample SAMPLE_ID --rollback --apply

# canary 驗證完成後批次搬移
PYTHONPATH=backend python scripts/migrate_tertiary_output_layout.py --all --apply
```

批次完成後先不要刪 `/home/n102968/NGS_UI/tertiary_output` 或 `/home/pipeline/tertiary_output`；待所有個案抽查完成再另排清理。遷移程式會略過舊 `snv_indel.annotated.tsv`、CNV/SV/STR/PGx 純複本，以 03 raw + 舊 enriched TSV 建 sparse overlay，因此 target 不會永久多留一份完整 SNV TSV。

三級分析進度條依實測 DRAGEN / in-house batch log 配重：等待同模式共享 cache 鎖時停在 `waiting-nextflow-cache`（2%）；取得鎖後 Nextflow 佔主要時間，但 `queued` 事件只寫 timing log，不推進 UI 百分比；實際 `start/done` 才依 process 權重更新。前處理、mito、STR、prepare CNV/SV 等快速步驟只佔少量進度，`ANNOTSV_SV`、VEP、Pangolin、parse CSQ、ACMG/PGx 等耗時步驟佔主要權重。Nextflow 結束約落在 82%，post-processing 依目前第幾個 sample 與子步驟推進到 99%，`promote-output` 維持 99%，切換完成才到 100%。Nextflow batch 會同時更新多個 process，因此後端寫入 `nextflow_progress_pct` 時會保留目前最大值，避免較早 process 晚更新造成進度條倒退。首頁會定期查詢後端 active job，因此其他瀏覽器或電腦登入後也會看到正在跑的三級分析；systemd unit 需設定 `KillMode=process`，避免重啟 web service 時把同一個 cgroup 裡的 Nextflow worker 一起 SIGTERM。若 persisted job state 仍是 `queued`/`running`，但 worker PID 已不存在或在 Linux 上已成 zombie，後端讀取 job 時會自動標記為 failed，避免 UI 狀態列永久停在最後百分比。進度面板提供「終止」按鈕，會呼叫後端取消該 job 並向 worker process group 送出終止訊號。三級分析清單右上角提供「清理 Nextflow 暫存」，可清空 `NGS_UI_TERTIARY_NF_WORK_ROOT`（預設 `$NGS_UI_HOME/nf_work`，dev 機為 `/home/n102968/NGS_UI/nf_work`）底下的共享 work 與 `contexts/` cache；若仍有 queued/running job，後端會拒絕清理。

DOCX CNV/SV 表格的「變異位置」欄使用 buffered wrap，內容寬度比欄寬少一格，避免座標字串貼到後面的欄位。

可用 `scripts/compare_genebe_spliceai_coverage.py` 抽樣 GeneBe 本機資料庫（例如 `/home/n102968/NGS_UI/biotools/genebe/genebe_hg38.tsv.gz`），再走目前 extra-VEP 的 VEP SpliceAI plugin 路徑產生對照，輸出 `summary.tsv`、`summary_by_kind.tsv`、`mismatches.tsv` 與 `run_metadata.tsv`，評估 GeneBe DB 的 SpliceAI 覆蓋率是否足以取代 extra-VEP / GeneBe API。正式估計建議使用 `--max-sites 100000 --sample-mode chrom-balanced`，避免只取 TSV 前段或讓大型染色體主導結果；未指定 `--seed` 時會由系統亂數產生並記錄在 metadata。

服務啟動時會以 daemon thread 在背景預熱 HPO、phenotype gene map、OMIM 與 mito ClinVar cache，避免解析大型 `phenotype_to_genes.txt` 期間擋住 HTTP port。完整 SNV TSV 的 gene search 走 per-sample `snv_gene_index.sqlite`；index 由 tertiary job 的 post-processing 預建，載入個案後不再自動掃描 WGS raw TSV。`snv_indel.review.tsv` 的 manifest 會記錄 test type 與 WES read-depth hard floor，SNV/Indel 會排除 chrM/MT，WES 樣本會在 review TSV 產生階段直接排除 DP < 20；SNV/Indel 主畫面、gene search 與 reviewer-marked 補點也都套用同一個 WES DP ≥ 20 gate。若舊樣本缺少或尚未更新 `snv_gene_index.sqlite`，報告區會對少量已標記 SNV id 做 raw TSV fallback scan，把 gene search 標記進 1/2/C 的點位補回 payload。Exomiser/LIRICAL 只有 active analysis 至少有一個 HPO term 才可排入；背景 job 會由完整 TSV 產生 `vcf_from_tsv.vcf.gz`，v3.5 多 transcript rows 會先依 `(CHROM,POS,REF,ALT)` 去重，避免同一 variant 在 VCF 重複出現。Secondary findings 透過 `/api/samples/{id}/secondary-snv` 背景載入，來源是 `snv_indel.review.tsv`，只計算 ACMG SF、中風與 Carrier panel，再套 ClinVar/ACMG P/LP、VAF ≥ 0.2、zygosity 非 ref 篩選；不掃完整 raw TSV，也不帶入 BED 外、低 DP、低 VAF/ref 點位。重疊 variant 仍會在每個命中的 panel 顯示，但 review state 跨 panel 共用；健檢 DOCX 才會合併去重。健檢 PGx 概覽會列出摘要表中的全部藥物並交代僅列於附錄的項數；HLA 多 allele 結果只要有 positive 就保留其 actionable 建議，不會因同列其他 negative 被排除。PDF gene list 由 `/api/samples/{id}/report-gene-list` lazy 載入。`in_panel` / `pheno_score` 由 active analysis 的 `pheno_score.tsv` 動態補入，不重寫完整 raw TSV。

三級分析另保存 `source_sample_id`（原始 sequencing sample ID）到 job state 與 `pipeline_source.json`。IGV 先依 sample ID suffix 判斷來源，`-dragen` 強制 DRAGEN roots、`-nckuh` / legacy `-inhouse` 強制 in-house roots；沒有 suffix hint 時才看 `pipeline_source.json` 的 `pipeline_type`。DRAGEN 找 `NGS_UI_DRAGEN_BAM_ROOT/<run>/bam/{sample}.bam`，預設 roots 包含 `/home/datalake_Raw/Novaseq` 與 `/home/datalake_Intermediate/n102968`；in-house 找 `NGS_UI_INHOUSE_BAM_ROOT/<batch>/{sample}/02_alignment/`。自訂輸出 ID 找不到時，會用 sidecar 回查原始 BAM ID；舊個案缺少 sidecar 欄位時，會移除已知 suffix 後用原始 sample ID 搜尋 BAM。IGV 的「其他路徑」資料夾下拉由 `/api/igv/bam-folders` 列出 configured BAM roots 內的 DRAGEN `<run>/bam/` 與 in-house `<batch>`，選 in-house batch 後會列出該 batch 底下各 sample `02_alignment` 的 BAM，排除 `.repeats.bam`。若 systemd 已設定 `NGS_UI_DRAGEN_BAM_ROOT`，需把 `/home/datalake_Intermediate/n102968` 併入該 env，否則會覆蓋程式預設。

### 臨床表徵工具 `/phenotype/`

獨立頁面（內網信任、無需登入）：搜尋 HPO term、套用 / 自訂 gene panel、把結果存成 token 之後在「載入新個案」帶入。Token 限 `[A-Za-z0-9_-]{1,32}`，內容 ≤64KB，panel ≤5000 個基因；自訂 panel 的基因 symbol 不會被轉大寫（`C7orf50` 保持原樣），但會先套用安全 alias 轉成 HGNC-current symbol。

「載入既有資料」會先清空目前頁面的 HPO terms、fixed panel chips、free panel rows 與舊預覽，再依新查詢結果填回；既有 phenotype 載入後會立即顯示預覽，之後新增、移除或修改 HPO/panel 也會即時重繪，不必先儲存。若下一位病人查無 phenotype，頁面與預覽會維持空白，不會沿用上一位病人的 HPO/panel。

主畫面的 Patient phenotype card 也提供 `WES-I / WES-II / WGS / Other panel` tabs，預設展開 `Other panel`；固定 panel 可直接點 chip 選取，其他 panel 仍可用 typeahead 搜尋。

固定 panel 的來源是 `reference/fixed_panel_sources/WES-I.xlsx`、`reference/fixed_panel_sources/WES-II.xlsx` 與 `reference/fixed_panel_sources/other_panel/`。更新 Excel 後執行 `PYTHONPATH=backend python scripts/import_fixed_panels.py`，會同步重建三個入口共用、且會進 git 的 `phenotype_data/fixed_panels/index.json` 與 `phenotype_data/gene_panels/*.txt`。WES Excel 只會匯入 `gene panel list` 標記列起始的基因區塊，避免把疾病名或資料來源列誤算成基因。WGS 固定套組另包含血液科 Lymphoid Neoplasm Panel 與 Myeloid Neoplasm Panel；來源 CSV 的 SNV/indel、CNV、STR、Mitochondria 欄位已合併成單一 gene list，不在 phenotype panel 層分 variant type。

HPO reference、固定 panel、custom panel 與既有 `pheno_score.tsv` 讀入時都會先透過 `ngs_panel_deadzone` 的 HGNC alias map 轉成 canonical gene symbol；從 `/phenotype/` 建立 custom panel 時，後端也會先套用 `ngs_panel_deadzone/panel/panel_gene_aliases.tsv` 的安全 alias 再寫入 repo 內 panel 檔。`panel_gene_aliases.tsv` 由 HGNC 官方 `reference/hgnc/hgnc_complete_set.txt`（`prev_symbol` / 唯一 `alias_symbol`）與 `reference/hgnc/withdrawn.txt`（唯一 merged/split replacement）加上 `reference/hgnc/manual_panel_aliases.tsv` 產生，重建指令為 `python scripts/build_hgnc_panel_aliases.py`；衝突項輸出到 `docs/ops/hgnc_alias_conflicts.tsv`，既有 custom panel 仍非 current HGNC 的項目列在 `docs/ops/custom_panel_hgnc_review_20260613.tsv`。Custom panel 檔案第一行可用 `#source:` 記錄來源（空白也可），loader 會略過註解行；`/phenotype/` 的 gene-list drawer 會顯示這個 source。SNV/CNV/SV/Mito 變異端也用同一套 canonicalization，再做 `pheno_score` / `in_panel` join。

`/phenotype/` 的 HPO term、fixed panel chip 與 panel 搜尋列都有「查看」按鈕，會呼叫 `GET /api/phenotype-tool/gene-list?kind=hpo|panel&key=...`，在右側 drawer 顯示 canonical gene list、來源、清單內篩選與複製功能，不把完整基因清單塞進主畫面。Topbar 的「搜尋基因」會呼叫 `GET /api/phenotype-tool/gene-memberships?gene=...`，反查某個 canonical gene 出現在哪些 HPO terms / panels。單一 HPO gene-list 查詢在 full phenotype scorer cache 尚未預熱完成時會走 fast path，只掃該 HPO term，避免第一次點「查看」被整份 `phenotype_to_genes.txt` 載入卡住。

---

## 5. 新增一個個案 / mitochondrial annotation

三級 pipeline 的 Mito source 是 `04_mito/{SAMPLE_ID}.mito.tsv`。UI 平常直接讀此檔；若本機 MITOMAP 表確實補到 `MITOMAP_*` 欄位，worker 才在 `08_postprocessing/mito.annotated.tsv` 留一份小型衍生檔並優先讀取。前端 adapter 同時相容舊欄位與 v3.2 欄位；Mito 判讀以 ClinVar P/LP、MITOMAP confirmed/pathogenic、gnomAD-mito rare 或缺值/已報告變異分成三層。本院 AF 不寫回 Mito TSV：FastAPI 啟動背景 warm-up 會透過 `.tbi` 與 `bcftools`/`tabix` 只讀 `NGS_UI_INHOUSE_AF_DB` 的 chrM records 並快取；DB、index 或查詢工具缺失時該欄靜默顯示 `—`。

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
- 首頁歡迎文字與版本紀錄放在 `frontend/VERSION.md`；首頁中間的橫向 NGS 分析流程由 `frontend/ANALYSIS_FLOW.md` 管理。流程檔使用固定的 `##` section 與 `|` 前 key 決定版面位置，一般文字修改只需編輯 `|` 後內容。前端在未載入個案時依序顯示歡迎訊息、流程圖與版本紀錄；流程檔缺失或格式不完整時只隱藏流程區。之後若有影響判讀流程、報告輸出、資料載入或主要工具入口的更新，需評估是否同步更新這兩份檔案。
- 大型 HPO reference 必須放在 `NGS_UI_PHENO_DATA_DIR`（`hp.obo`、`phenotype_to_genes.txt` 等）；固定與自訂 panel data 則在 repo 的 `phenotype_data/fixed_panels/`、`phenotype_data/gene_panels/` 與 `phenotype_data/custom_panels/`，會跟著 git 更新。
- GenCC / ClinGen / MONDO 補充疾病資料用 `python scripts/update_gene_disease_db.py` 更新；若正式機不能出網，先手動下載 `gencc_submissions.csv`、`clingen_gene_validity.csv`、`mondo.json` 到 `NGS_UI_GENE_DISEASE_RAW_DIR`，再跑 `python scripts/update_gene_disease_db.py --skip-download`。
- EMR 相關功能預設停用，需設 `NGS_UI_EMR_CLIENT_ID` 才會啟用，且只在內網可達。
- `/api/phenotype-tool/*` 與 `/api/healthz` 是刻意公開無認證；`/api/patient_list` 與其餘 `/api/*` 需登入。
- 大型 JSON response 會在瀏覽器支援時自動 gzip；SNV parse + phenotype / Exomiser / LIRICAL / OMIM join 使用有上限的 process-local LRU cache，輸入 TSV 或 sidecar 更新後自動失效。
