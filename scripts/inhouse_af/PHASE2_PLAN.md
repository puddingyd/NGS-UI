# Phase 2 — 把 in-house AF 接進 NGS-UI（實作 plan）

> Phase 0/1 已完成：全量 677 樣本增量版 `inhouse_af.hg38.vcf.gz` 建好並對 GLnexus
> 驗證通過（AN≥1300：SNP r=0.9999 / indel r=0.9960，且為 GLnexus 超集合）。
> 這份是 **Phase 2（UI 整合）** 的實作清單。整合要在 **UI 主 branch 的最新版**
> 上開新 branch 做（見文末），不是在 inhouse_af 這個分支。

模式完全比照現有的 **GeneBe / GIAB** annotation：一支 TSV 加欄 script → 塞進
`run_stopgaps.sh` → adapter 帶出 → 前端顯示。**不是 `bcftools annotate`**
（`snv_indel.annotated.tsv` 是 TSV 不是 VCF）。

## ⚠️ 命名
`config.py:211 INHOUSE_VCF_ROOTS` 的 "in-house" 已經是「在院（NCKUH）定序來源」的意思
（對比 DRAGEN），前端大量 `inhouse` 也是這個。AF 這個新功能一律用
**`INHOUSE_AF_*`（後端欄位）/ `inhouse_af`（payload）/ `AF_nckuh`（畫面標籤）**，
不要跟「定序來源」混淆。

## A. 資料 / 設定層
- [ ] 把 DGX 的 `inhouse_af.hg38.vcf.gz`(+`.tbi`) 放到部署機
      `NGS_UI_HOME/biotools/inhouse_af/`（不進 git，同 gnomAD/GeneBe DB）。
- [ ] `config.py`：新增 `INHOUSE_AF_DB`（env `NGS_UI_INHOUSE_AF_DB`，預設
      `BIOTOOLS_DIR/inhouse_af/inhouse_af.hg38.vcf.gz`），**缺檔靜默 disable**
      —— 比照 `GENEBE_DB`（`config.py:120`）。

## B. 註記層（新 script + 接 stopgaps）
- [ ] 新 `scripts/annotate_inhouse_af.py`（比照 `annotate_giab_strata.py`）：
  - 讀 `snv_indel.annotated.tsv`，依 `(CHROM,POS,REF,ALT)` 對 `INHOUSE_AF_DB`
    **merge-join**（TSV 排序後 vs `bgzip -dc` 串流單次掃描；低記憶體，避免對
    WGS 數百萬列做隨機 tabix seek）。**決定：用 merge-join。**
  - 寫入 `INHOUSE_AC / INHOUSE_AN / INHOUSE_AF` 三欄，fill-or-augment +
    atomic replace，缺 DB → no-op `exit 0`。
  - 正規化：lookup key 要和 sites VCF 一致（都 left-align / `norm -m-`），
    否則複雜 indel 會漏配（同 GeneBe 的已知侷限，可接受）。
- [ ] `run_stopgaps.sh` 插一步 `inhouse-af`，位置在 **giab-strata 之後
      （`run_stopgaps.sh:187` 之後）、review-tsv/gene-index 之前**。
      加欄會整檔改寫、位移 byte offset，而 gene-index 本就排在其後重建，順序天然正確。
      同步加 `--skip-inhouse-af` 開關（比照 `--skip-giab`）。
- [ ] 新 `scripts/backfill_inhouse_af.sh`（比照 `backfill_giab_strata.sh`）：
      舊樣本補註記 → **重建 review TSV → 重建 gene index**（三個一起，因 offset
      位移——GIAB backfill 踩過的雷）。每批增量更新完 `inhouse_af.hg38.vcf.gz`
      後，用它把既有樣本補上。

## C. 後端 payload
- [ ] `adapters/snv_tsv.py`：在 `AF` 那段（`~514`）旁邊加
  - `"inhouse_af": _to_num(row.get("INHOUSE_AF"))`
  - `"inhouse_ac": _to_int(row.get("INHOUSE_AC"))`
  - `"inhouse_an": _to_int(row.get("INHOUSE_AN"))`
- review TSV 保留全欄（`snv_review.py:187` 是列過濾、不砍欄，已確認）；gene search
  走 raw TSV。兩邊都會帶到，主畫面即可顯示。

## D. 前端（`app.js` / `index.html` / `style.css`）—— 只做主畫面
- [ ] **AF 那排（`app.js:3094-3100`）改版：**
  - 在 `AF_eas` 底下新增一列 **`AF_nckuh`**，格式**比照 AF / AF_eas**：顯示實際
    數字（`fmtNum(v.inhouse_af,5)`，**不是 %**），後面括號 `(AC/AN)`。
    例：`0.009 (12/1354)`。`inhouse_af` 缺值顯示 `—`（AC/AN 也一併省略）。
  - 把現有的 **`1000G EAS`（`TG_eas_af`）移進 More**（例如挪到 `app.js:3103`
    的 `more-extras`），讓主排只留 AF / AF_eas / AF_nckuh。
- [ ] **不做** in-house AF 的顯示 filter（原 plan 第 8 點取消）。
- [ ] tooltip / 說明：本院 AF 是**疾病轉介族群**，用來標「本院常見 / 重現 artifact」，
      **不作 ACMG BA1/BS1**（那仍看 gnomAD）。N=677 → 最小可解析 AF ≈ 1/1354，
      所以要顯示 AC/AN 讓樣本數看得見。
- [ ] **個案清單 / DOCX：先不做**（決定：只做主畫面），之後再評估。

## E. 文件
- [ ] 更新 `AGENTS.md`（in-house AF 資料流 + `INHOUSE_AF_*` 欄位 + stopgaps 新步驟）、
      `scripts/inhouse_af/README.md`（把「Phase 2 not done」改成實作說明）、
      `frontend/VERSION.md`（首頁版本紀錄）。

## 已敲定的決定
| 決定 | 選擇 |
|---|---|
| annotate 方式 | **merge-join**（非 tabix 逐點） |
| 顯示位置/格式 | 主排新增 `AF_nckuh` 列，`0.009 (12/1354)` 樣式；1000G 移進 More |
| in-house AF filter | **不做** |
| 個案清單 / DOCX | 只做主畫面，暫不進清單/報告 |

## 踩雷備忘
- annotate 改寫整檔 → gene index 一定要一起重建（backfill 尤其）。
- 正規化不一致 → 複雜 indel 漏配（GeneBe-like，可接受）。
- 命名別撞 `INHOUSE_VCF_ROOTS`（= 定序來源）。
- in-house AF 不作族群頻率 ACMG 證據。
