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
│   └── {LIS_ID}/          ← snv_indel.annotated.tsv, cnv.annotated.tsv,
│                             sv.annotated.tsv, mito.annotated.tsv,
│                             sample_metadata.json, qc_summary.json,
│                             roh_summary.json, analyses/{ver}/...
├── patient_phenotype/     ← {LIS_ID}_{MRN}_phenotype.txt（自動帶入 HPO）
├── patient_list/          ← 上傳的「未完成報告清單」xlsx + 衍生 roster.json
├── phenotype_data/        ← hp.obo, phenotype_to_genes.txt, gene_panels/*.txt
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
| `NGS_UI_OMIM_XLSX` | `$NGS_UI_HOME/OMIM/OMIM.xlsx` | OMIM 疾病註解（缺檔 = 停用） |
| `NGS_UI_BIOTOOLS_DIR` | `$NGS_UI_HOME/biotools` | Exomiser / LIRICAL |
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
   - HPO / gene panel 可在這裡選；若存在 `patient_phenotype/{LIS}_{MRN}_phenotype.txt` 會自動讀入；
   - 勾「登錄後開始分析」會順便把 Exomiser/LIRICAL 排入佇列。
3. **看變異卡片** — 個案載入後先顯示 SNV/Indel（分段載入），CNV/SV 與 Mitochondria 在背景載完後補上：
   - SNV/Indel tier：`1A / 1B / 1C / 2 / 3`（互斥）
   - CNV：`CNV-1A`（Clinical）、`CNV-1B`（Pathogenic）；SV：`SV-2A / SV-2B`
   - Mitochondria：`MITO-1`（Pathogenic）、`MITO-2`（Disease-associated）— 只列 `FILTER=PASS` 且具 MITOMAP 疾病關聯/致病性的位點
4. **標記與判讀** — 在每個變異上標 causative / candidate / other，編輯 ACMG/分類、寫 comment；變更會自動存到 `tertiary_output/{LIS}/analyses/{ver}/analysis.json`。
5. **匯出報告** — `GET /api/samples/{LIS}/report.docx`（UI 上的「匯出」按鈕），目前涵蓋 SNV/Indel；CNV/SV/Mito 的 docx 匯出在 TODO（見 `CLAUDE.md`）。

### 臨床表徵工具 `/phenotype/`

獨立頁面（內網信任、無需登入）：搜尋 HPO term、套用 / 自訂 gene panel、把結果存成 token 之後在「載入新個案」帶入。Token 限 `[A-Za-z0-9_-]{1,32}`，內容 ≤64KB，panel ≤5000 個基因；自訂 panel 的基因 symbol 不會被轉大寫（`C7orf50` 保持原樣）。

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
- `phenotype_data/` 必須放對位置（`NGS_UI_PHENO_DATA_DIR` 沒有 fallback）；正式環境記得 `mv NGS_UI/NGS-UI/phenotype_data NGS_UI/phenotype_data`。
- EMR 相關功能預設停用，需設 `NGS_UI_EMR_CLIENT_ID` 才會啟用，且只在內網可達。
- `/api/phenotype-tool/*` 與 `/api/healthz` 是刻意公開無認證；`/api/patient_list` 與其餘 `/api/*` 需登入。
