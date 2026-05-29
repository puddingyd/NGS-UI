# DGX-2 NGS Pipeline 預設執行問題回報：Nextflow work directory 與 PyTensor cache 權限

## 問題摘要

目前 DGX-2 NGS pipeline 在多使用者環境下，若照現行 guide 的預設方式執行 WES + gCNV，可能會遇到兩類問題：

1. **Nextflow launch/work directory 共用導致 cache DB lock 或 file lock 問題**
2. **GATK gCNV / PyTensor cache directory 權限問題**

這兩個問題會使一般使用者即使正確準備 `samplesheet.csv`，仍無法順利照預設指令完成分析。

---

## 一、Nextflow launch/work directory 問題

### 現行 guide 的執行方式

目前 guide 建議使用者：

```bash
cd /raid/DGM/work
```

然後直接執行：

```bash
nextflow -c ${PIPELINE_CONFIG} run ${PIPELINE_CODE}/main.nf \
    -profile dgx \
    --input_csv ${OUT_DIR}/samplesheet.csv \
    --seq_type WES \
    --run_gcnv true \
    --out_dir ${OUT_DIR} \
    -resume
```

### 實際遇到的錯誤

實際執行時曾出現：

```text
ERROR ~ Can't open cache DB: /raid/DGM/work/.nextflow/cache/.../db

Nextflow needs to be executed in a shared file system that supports file locks.
Alternatively, you can run it in a local directory and specify the shared work
directory by using the -w command line option.
```

### 問題分析

如果所有 batch 都直接在 `/raid/DGM/work` 執行 Nextflow，會導致不同使用者或不同 batch 共用：

```bash
/raid/DGM/work/.nextflow/cache
```

這會增加 Nextflow cache DB lock、metadata 混雜或 resume 行為不穩定的風險。

Nextflow 比較合理的使用方式應是：

- 使用獨立的 launch directory 儲存 `.nextflow` metadata
- 使用 `-w` 指定 process work directory

### 建議修正方式

建議 guide 或 pipeline wrapper 預設改為每個 batch 使用獨立 launch directory 與獨立 work directory，例如：

```bash
BATCH_NAME="260603_ped_WES"

OUT_DIR="/datalake_Intermediate/pipeline/nextflow_output/${BATCH_NAME}"
LAUNCH_DIR="/datalake_Intermediate/pipeline/nextflow_launch/${BATCH_NAME}"
WORK_DIR="/raid/DGM/work/${BATCH_NAME}"

mkdir -p "${LAUNCH_DIR}" "${WORK_DIR}"

cd "${LAUNCH_DIR}"

nextflow -c ${PIPELINE_CONFIG} run ${PIPELINE_CODE}/main.nf \
    -profile dgx \
    --input_csv "${OUT_DIR}/samplesheet.csv" \
    --seq_type WES \
    --run_gcnv true \
    --out_dir "${OUT_DIR}" \
    -w "${WORK_DIR}" \
    -resume
```

這樣可以避免所有使用者或所有 batch 共用 `/raid/DGM/work/.nextflow/cache`。

---

## 二、GATK gCNV / PyTensor cache 權限問題

### 執行情境

在執行 WES + gCNV pipeline 時，例如：

```bash
BATCH_NAME="260603_ped_WES"

nextflow -c ${PIPELINE_CONFIG} run ${PIPELINE_CODE}/main.nf \
    -profile dgx \
    --input_csv "${SAMPLESHEET}" \
    --seq_type WES \
    --run_gcnv true \
    --out_dir "${OUT_DIR}" \
    -w "${WORK_DIR}" \
    -resume
```

pipeline 在 GATK gCNV `DetermineGermlineContigPloidy` 步驟失敗。

### 錯誤訊息

錯誤訊息顯示：

```text
ValueError: compiledir '/raid/DGM/pytensor_cache/n102968' exists but you don't have read, write or listing permissions.

java.lang.RuntimeException: A required Python package ("gcnvkernel") could not be imported into the Python environment.
```

外層錯誤看起來像是 `gcnvkernel` 無法 import，但實際原因是 `gcnvkernel` import 時會載入 `pymc` 與 `pytensor`，而 PyTensor 嘗試使用一個目前使用者沒有寫入權限的 compiledir。

### 權限檢查結果

執行：

```bash
ls -ld /raid/DGM/pytensor_cache
ls -ld /raid/DGM/pytensor_cache/n102968
```

結果：

```text
drwxrwxrwx 6 n101569 dgm_nckuh 20480 May 22 03:22 /raid/DGM/pytensor_cache
drwxr-xr-x 2 n101569 dgm_nckuh 4096 May 22 01:51 /raid/DGM/pytensor_cache/n102968
```

也就是：

```text
/raid/DGM/pytensor_cache/n102968
```

雖然目錄名稱是目前 user ID `n102968`，但 owner 是 `n101569`，所以 `n102968` 無法寫入。

### 進一步檢查 `.command.run`

檢查 Nextflow work directory 中的 `.command.run`：

```bash
grep -n "apptainer\|singularity\|PYTENSOR\|THEANO" .command.run .command.sh .command.err
```

發現 Singularity 執行時 pipeline config 硬塞了：

```bash
--env PYTENSOR_FLAGS=compiledir=/raid/DGM/pytensor_cache/${currentUser}
```

實際展開後變成：

```bash
--env PYTENSOR_FLAGS=compiledir=/raid/DGM/pytensor_cache/n102968
```

因此即使使用者在 shell 中另外 export：

```bash
export PYTENSOR_FLAGS="compiledir=/raid/DGM/work/260603_ped_WES/pytensor_cache"
```

仍會被 Nextflow/Singularity `runOptions` 裡的：

```bash
--env PYTENSOR_FLAGS=compiledir=/raid/DGM/pytensor_cache/${currentUser}
```

覆蓋。

---

## 目前 workaround

目前可行的 workaround 是：

1. 建立 batch 專用 `LAUNCH_DIR`
2. 建立 batch 專用 `WORK_DIR`
3. 在 `WORK_DIR` 下建立 batch 專用 `pytensor_cache`
4. 複製一份 local config 到 batch 專用 launch directory
5. 將 local config 內的 PyTensor compiledir 改成該 batch 自己可寫入的路徑
6. 使用這份 local config 執行 Nextflow

範例：

```bash
BATCH_NAME="260603_ped_WES"

OUT_DIR="/datalake_Intermediate/pipeline/nextflow_output/${BATCH_NAME}"
SAMPLESHEET="${OUT_DIR}/samplesheet.csv"
LAUNCH_DIR="/datalake_Intermediate/pipeline/nextflow_launch/${BATCH_NAME}"
WORK_DIR="/raid/DGM/work/${BATCH_NAME}"
PYTENSOR_CACHE="${WORK_DIR}/pytensor_cache"

mkdir -p "${LAUNCH_DIR}" "${WORK_DIR}" "${PYTENSOR_CACHE}"
chmod -R u+rwx "${PYTENSOR_CACHE}"

source /datalake_Intermediate/pipeline/pipeline_code/NGS2ndAnalysis_env.sh

LOCAL_CONFIG="${LAUNCH_DIR}/nextflow_main.${BATCH_NAME}.config"

cp "${PIPELINE_CONFIG}" "${LOCAL_CONFIG}"

sed -i "s#compiledir=/raid/DGM/pytensor_cache/\${currentUser}#compiledir=${PYTENSOR_CACHE}#g" "${LOCAL_CONFIG}"
sed -i "s#compiledir=/raid/DGM/pytensor_cache/n102968#compiledir=${PYTENSOR_CACHE}#g" "${LOCAL_CONFIG}"
sed -i "s#compiledir=/home/pipeline/pytensor_cache/\${currentUser}#compiledir=${PYTENSOR_CACHE}#g" "${LOCAL_CONFIG}"

cd "${LAUNCH_DIR}"

nextflow -c "${LOCAL_CONFIG}" run ${PIPELINE_CODE}/main.nf \
    -profile dgx \
    --input_csv "${SAMPLESHEET}" \
    --seq_type WES \
    --run_gcnv true \
    --out_dir "${OUT_DIR}" \
    -w "${WORK_DIR}" \
    -resume
```

這個 workaround 不會修改共用的主 config，只影響該 batch。

---

## 建議 pipeline 端修正方向

### 1. Guide 預設改為每個 batch 獨立 launch/work directory

建議不要再建議使用者直接：

```bash
cd /raid/DGM/work
```

而是預設：

```bash
LAUNCH_DIR="/datalake_Intermediate/pipeline/nextflow_launch/${BATCH_NAME}"
WORK_DIR="/raid/DGM/work/${BATCH_NAME}"

mkdir -p "${LAUNCH_DIR}" "${WORK_DIR}"
cd "${LAUNCH_DIR}"

nextflow ... -w "${WORK_DIR}" -resume
```

這樣可以避免不同 batch 共用 `/raid/DGM/work/.nextflow/cache`。

---

### 2. PyTensor cache 不要固定使用 `/raid/DGM/pytensor_cache/${currentUser}`

目前 config 中類似以下設定在多使用者環境下容易出錯：

```bash
--env PYTENSOR_FLAGS=compiledir=/raid/DGM/pytensor_cache/${currentUser}
```

建議改成以下其中一種方式。

#### 方案 A：每個 batch 使用自己的 PyTensor cache

```bash
PYTENSOR_CACHE="${WORK_DIR}/pytensor_cache"
```

並在 Singularity runOptions 中指定：

```bash
--env PYTENSOR_FLAGS=compiledir=${PYTENSOR_CACHE}
```

#### 方案 B：保留 user-level cache，但執行前檢查權限

例如：

```bash
PYTENSOR_CACHE="/raid/DGM/pytensor_cache/${USER}"

if [[ -e "${PYTENSOR_CACHE}" && ! -w "${PYTENSOR_CACHE}" ]]; then
    echo "ERROR: ${PYTENSOR_CACHE} exists but is not writable by ${USER}"
    exit 1
fi

mkdir -p "${PYTENSOR_CACHE}"
chmod u+rwx "${PYTENSOR_CACHE}"
```

#### 方案 C：提供 pipeline 參數讓使用者指定

例如：

```bash
--pytensor_cache_dir /path/to/cache
```

再由 config 帶入 Singularity runOptions。

---

### 3. 避免 hard-code 或覆蓋使用者 shell 中的 `PYTENSOR_FLAGS`

目前即使使用者在 shell 中 export `PYTENSOR_FLAGS`，仍會被 config 內的 `--env PYTENSOR_FLAGS=...` 覆蓋。

建議允許使用者透過參數或環境變數覆寫，例如：

```groovy
params.pytensor_cache_dir = params.pytensor_cache_dir ?: "${workDir}/pytensor_cache"
```

或至少在 guide 中說明如何指定。

---

### 4. 檢查 Nextflow 25.10.2 config 相容性

執行時也曾遇到：

```text
ERROR ~ Unknown config attribute `singularity.params.sif_dir`
```

建議檢查 `nextflow_main.config` 中是否有不被 Nextflow 25.10.2 接受的舊寫法，例如：

```groovy
singularity.params.sif_dir
```

可能需要改成目前 Nextflow 支援的設定方式。

---

## 總結

目前 WES + gCNV 預設流程在多使用者環境下會遇到兩個主要問題：

1. 所有 batch 若都從 `/raid/DGM/work` launch Nextflow，會共用 `/raid/DGM/work/.nextflow/cache`，可能導致 cache DB lock 或 file lock 問題。
2. GATK gCNV 會透過 `PYTENSOR_FLAGS` 使用 `/raid/DGM/pytensor_cache/${currentUser}`，但該目錄可能已由其他使用者建立，造成目前使用者無法寫入，導致 PyTensor import 失敗，進而讓 GATK 報 `gcnvkernel` 無法 import。

建議 pipeline 預設改成：

```bash
LAUNCH_DIR="/datalake_Intermediate/pipeline/nextflow_launch/${BATCH_NAME}"
WORK_DIR="/raid/DGM/work/${BATCH_NAME}"
PYTENSOR_CACHE="${WORK_DIR}/pytensor_cache"
```

並在 Nextflow 指令中使用：

```bash
-w "${WORK_DIR}"
```

同時讓 Singularity/Apptainer 使用該 batch 自己可寫入的 PyTensor cache。

這樣其他使用者才能照 guide 預設方式執行，不需要手動複製 local config 或修改 cache 路徑。
