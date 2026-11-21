# Pangolin 分數快取（PANGOLIN_SCORE 加速）

## 為什麼要做

DGM 的 RTX 2000 Ada 在跑 Pangolin 時，GPU 已經跑滿（功耗約 66 W，上限 70 W），每例 WES 約 3,400 個 splice 候選，要花 8 分鐘左右。36 例的批次光 Pangolin 就要 4–5 小時。換 GPU 的跑法（同時開多個、一次送多個變異）已經沒有什麼空間。

但 Pangolin 分數只跟變異本身（CHROM/POS/REF/ALT）、參考序列、gene model DB、Pangolin 容器版本和參數有關，**跟是哪個 sample 無關**。統計現有 803 份 `02_pangolin` 輸出：總共 927 萬筆，不重複的變異只有 15 萬個，重複率 98.4%。

所以只要把算過的分數存起來，每例只把沒算過的變異送給 GPU，其他直接查表。

## 改了什麼

| 檔案 | 內容 |
|---|---|
| `pangolin_cache.py` | 快取工具（只用 Python 標準庫加 SQLite，在 Pangolin 容器裡或主機上都能跑） |
| `PANGOLIN_SCORE.nf` | 取代 `modules/snv_annotation.nf` 裡的 `PANGOLIN_SCORE` |

`PANGOLIN_SCORE.nf` 和正式版相比有三點不同：

1. **修正 `params.containsKey(...)`**：對 `params` 呼叫方法，會讓 Nextflow 把整個 params 算進 task hash。每個 NGS-UI job 的 `out_dir`、`samplesheet` 都不一樣，所以 PANGOLIN_SCORE 從來沒有被 `-resume` 沿用過。改成直接讀 `params.use_gpu_lock ?: false`、`params.pangolin_num_gpus ?: 1`。
2. **分數快取**：`params.pangolin_cache_db` 是非空字串時啟用：
   - `split`：查表，只把沒算過的變異寫成 `pangolin_miss.vcf`。
   - Pangolin 只跑 `pangolin_miss.vcf`；全部都查得到時就完全不啟動模型。
   - `merge`：依原本的列順序，把查到的值和新算的值組回 `{sample}.pangolin.vcf`，再把新算的值寫進快取。

   `pangolin_cache_db = ''` 時，執行的就是原本的 Pangolin 指令，行為和現在完全相同。
3. **Pangolin 參數只寫一次**：`--distance 50`、`-m True` 定義成變數，命令列和快取分區共用，兩邊不會不一致。

**輸出不變**：下游 `parse_vep_csq.load_pangolin_scores` 只讀 CHROM/POS/REF/ALT 和 `INFO/Pangolin`，這兩樣和直接跑 Pangolin 完全相同，所以 `PANGOLIN_SCORE`/`PANGOLIN_DETAIL` 欄位和 ACMG BP7 都不變。唯一的差別是 VCF 其他欄位保留原始格式。原本 Pangolin 用 PyVCF 重新輸出，例如 QUAL 的 `50` 會變成 `50.0`。

### 快取的設計

- **分區（namespace）**：快取依 `容器檔名 | gene model DB（檔名+大小）| 參考序列（檔名+大小）| distance | mask` 分區。換容器、換 DB、換參考序列或改參數，都會自動開一個新的空白分區，不會拿到舊條件下算出的分數。
- **沒有分數也記**：Pangolin 有跑、但沒給分數的變異（例如不在任何 gene model 上），會存成 NULL，下次不會重算。
- **只增不改**：已經存在的值不會被覆寫。每一列都記錄來源 sample、主機和時間，供評鑑追溯。
- **安全檢查**：Pangolin 輸出缺了任何一個送去計算的變異，`merge` 就會失敗，這個 task 的結果也不會寫進快取。
- **同時寫入**：SQLite 使用 WAL 模式，等待上限 10 分鐘。DGX 最多 6 個 Pangolin 同時寫入也不會衝突。快取檔要放在本機磁碟，不要放在 NFS 上。
- **大小**：15 萬個變異約 15 MB。

## 部署步驟

以下以 DGM 為例。`REF=/home/pipeline/reference/hg38`，請依實際的 `params.ref_fasta` 和 `params.pangolin_db` 調整。

### 0. 確認容器裡的 Python 有 sqlite3

```bash
apptainer exec /home/pipeline/nextflow_containers/pangolin_cu121_1.0.0.sif python3 -c "import sqlite3, sys; print(sys.version.split()[0], sqlite3.sqlite_version)"
```

### 1. 放入工具

```bash
cp docs/ops/pangolin_cache/pangolin_cache.py /home/pipeline/tertiary_code/scripts/
```

### 2. 在 config 宣告快取位置

在 `nextflow_tertiary.config` 從 `ref_dir` 推導路徑的全域 params 區塊（`pangolin_db` 那一帶）加上：

```groovy
pangolin_cache_db = "${params.ref_dir}/tertiary/pangolin/cache/pangolin_scores.sqlite"
```

用執行 Nextflow 的帳號（也就是跑 NGS-UI 服務的帳號）確認有寫入權限：

```bash
mkdir -p /home/pipeline/reference/hg38/tertiary/pangolin/cache && touch /home/pipeline/reference/hg38/tertiary/pangolin/cache/.write_test && rm /home/pipeline/reference/hg38/tertiary/pangolin/cache/.write_test && echo ok
```

### 3. 替換 PANGOLIN_SCORE

```bash
cp /home/pipeline/tertiary_code/modules/snv_annotation.nf /home/pipeline/tertiary_code/modules/snv_annotation.nf.bak
```

把 `snv_annotation.nf` 裡從 `process PANGOLIN_SCORE {` 到對應的 `}` 換成 `PANGOLIN_SCORE.nf` 裡的 process（檔頭說明可以不貼）。

### 4.（建議）用現有結果預先建立快取

先查出正式設定的分區字串，確認路徑正確：

```bash
python3 /home/pipeline/tertiary_code/scripts/pangolin_cache.py namespace --sif-name pangolin_cu121_1.0.0.sif --pangolin-db /home/pipeline/reference/hg38/tertiary/pangolin/gencode.v47.annotation.db --ref-fasta <params.ref_fasta> --distance 50 --mask True
```

`db=` 和 `ref=` 後面要顯示檔案大小，不能是 `missing`。

先試跑（不寫入），並產生衝突清單：

```bash
python3 /home/pipeline/tertiary_code/scripts/pangolin_cache.py seed --dry-run --cache /home/pipeline/reference/hg38/tertiary/pangolin/cache/pangolin_scores.sqlite --sif-name pangolin_cu121_1.0.0.sif --pangolin-db /home/pipeline/reference/hg38/tertiary/pangolin/gencode.v47.annotation.db --ref-fasta <params.ref_fasta> --distance 50 --mask True --conflicts ~/pangolin_seed_conflicts.tsv /home/pipeline/tertiary_output/*/02_pangolin/*.pangolin.vcf.gz
```

- 同一個變異在不同檔案裡分數不同，就是「衝突」。這些變異**不會**寫入快取，之後遇到會重新計算。目前約有 111 個，可以用 `head ~/pangolin_seed_conflicts.tsv` 看看是哪種差異：小數末位不同，通常是 DGX 和 DGM 不同 GPU 的浮點差異；gene ID 不同，則是不同 DB 版本。
- 如果較舊的輸出是用別的容器或 DB 算的，加上 `--since YYYY-MM-DD`，只匯入目前設定上線後產生的檔案。
- 沒有問題就把 `--dry-run` 拿掉，再跑一次正式寫入。

預先建立快取還有一個好處：第 3 步改了 PANGOLIN_SCORE 的程式，Nextflow 的 cache 會失效一次，所有 sample 的 Pangolin 都要重跑。有了快取，這次重跑幾乎全部查表，很快就結束。

### 5. 驗證

1. 先把 2–3 個已經跑過的 sample 的 `02_pangolin/*.pangolin.vcf.gz` 複製到別處。重跑時 `publishDir` 會覆蓋原檔。
2. 用 NGS-UI 重跑這幾個 sample。
3. 在 Pangolin task 的 `.command.err` 會看到：
   ```
   [PANGOLIN-CACHE] <sample>: candidates=3400 cached=3390 to_score=10 namespace=...
   ```
4. 比對新舊輸出，`differ=0` 加上 `only_left=0 only_right=0` 就代表完全一致：
   ```bash
   python3 /home/pipeline/tertiary_code/scripts/pangolin_cache.py compare <舊的>.pangolin.vcf.gz /home/pipeline/tertiary_output/<sample>/02_pangolin/<sample>.pangolin.vcf.gz
   ```
5. 再原樣跑一次同一批，Nextflow 狀態裡的 PANGOLIN_SCORE 應該顯示 `cached`。這表示 `containsKey` 的修正生效了。
6. 查看快取狀態：
   ```bash
   python3 /home/pipeline/tertiary_code/scripts/pangolin_cache.py stats --cache /home/pipeline/reference/hg38/tertiary/pangolin/cache/pangolin_scores.sqlite
   ```
   正常只會有一個分區。出現兩個，代表第 4 步預建時用的路徑和 pipeline 實際用的不一樣。

### 停用 / 回復

- 暫時停用：config 改成 `pangolin_cache_db = ''`，或在命令列加 `--pangolin_cache_db ''`。
- 完全回復：還原 `snv_annotation.nf.bak`。

## 建議一併設定（非必要）

上次 24WE0166 是在負載高時，Apptainer 啟動 GPU 容器逾時失敗（`--nv` 的 `/bin/sed` 一分鐘 timeout）。可以在 dgm profile 的 `withLabel: process_gpu` 加上自動重試：

```groovy
errorStrategy = { task.attempt <= 2 ? 'retry' : 'terminate' }
maxRetries = 2
```

`errorStrategy`/`maxRetries` 不屬於 task hash，加了不會讓既有的 Nextflow cache 失效。
