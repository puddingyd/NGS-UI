# 表型引擎文獻驗證（計畫第 11 章、第 15 章階段 6）

| 檔案 | 內容 |
|---|---|
| `dev.txt` | 開發資料 400 例（第 12.1 節先導實測）。表型引擎參數在此選定 |
| `fresh.txt` | 第二組 481 例。2026-11-24 起與開發資料合併估計名次點數，不再是未接觸的驗證組 |
| `test.txt` | 鎖定測試組 152 例，2026-11-24 抽一次；之後不得重抽或依結果調整參數 |
| `results_20261124/` | 第一次完整報告（`report.md`、`summary.json`、每例名次 `rows_*.tsv`） |

病例 ID 是公開的 phenopacket-store ID，不含院內資料。

## 測試組怎麼抽的

time-split 合格（HPO 2023-10-09 知識庫內有該疾病、病例文獻沒有被該疾病的舊註記引用、≥3 個 HPO、有致病基因）的 4,722 例中，排除開發與第二組的病例，以及和它們同一篇文獻的病例，剩 1,060 例；每個 cohort（基因資料夾）最多 3 例，seed 2026，共 152 例。開發組與第二組已用掉所有合格 cohort 的前 3 例，所以無法再抽到「基因完全不重疊」的測試組（計畫第 11.3 節原本的設計）；這組是「相同基因、不同病人與文獻」。

## 重跑

需要 phenopacket-store v0.1.27（`all_phenopackets.zip` 解壓）、HPO 2023-10-09 與目前版本的 `hp.obo`／`phenotype.hpoa`／`genes_to_disease.txt`、`mondo.obo`。先各建一個知識庫快照（放在獨立資料夾，不要啟用），再跑：

```bash
.venv/bin/python scripts/benchmark_pheno.py run \
  --cases phenopacket-store/0.1.27 \
  --old-snapshot bench_kb/2023/2023-10-09_2023-10-09 --old-gene-map hpo_2023-10-09/genes_to_disease.txt \
  --current-snapshot bench_kb/cur/2026-09-01_2026-09-02 --current-gene-map genes_to_disease.txt \
  --out pheno_benchmark_out
```

疾病–基因對應用 HPO 官方 `genes_to_disease.txt`（MENDELIAN／UNKNOWN），time-split 用 2023 版，避免新知識洩漏；約 5 分鐘。

LIRICAL／Exomiser 對照：把每例結果整理成 `{case}.tsv`（欄位 `gene_symbol`、`disease_id`，依名次排列），以 `--lirical DIR`、`--exomiser DIR` 讀入；第一次報告沒有這兩項（見計畫第 12.3 節）。
