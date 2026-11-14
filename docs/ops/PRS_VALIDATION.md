# PRS 離線來源與計分驗證

v10.35，2026-11-14。這些工具只讀輸入、寫入指定的新稽核檔，不會啟用模型、更新個案或修改風險參數。結果範圍與待資料項目見[本次執行紀錄](../research/prs/PRS_VALIDATION_2026-11-14.md)。

## 一、準備資料與輸出位置

在 NGS-UI repo 根目錄執行 `python3 -m research.prs.validation`；院內也可用 `.venv/bin/python`。輸出放既有 PRS work root 之下的獨立稽核目錄，例如 `NGS_UI_PRS_WORK_ROOT/_validation_v1035/`；未自訂時是三級 work root 的 `prs/_validation_v1035/`。不要把稽核產物放回原 PRS 專案、biotools 模型或 reference。

下列變數由操作者按現有 runtime config 填入，都是既有資源；範例不是新增平台環境設定。先確認磁碟／使用者可寫入輸出目錄，使用全新檔名；CLI 拒絕覆寫。

```bash
PRS_AUDIT=/home/n102968/NGS_UI/nf_work/prs/_validation_v1035
PRS_SCORES=/path/to/pinned/scorefiles
PRS_CHD=/path/to/PGS004879_hmPOS_GRCh38.txt.gz
PRS_MODELS=/path/to/biotools/prs/multi_v1
mkdir -p "$PRS_AUDIT"
```

`PRS_SCORES` 放五個固定的 `PGSxxxxxx_hmPOS_GRCh38.txt.gz`，CHD 另指定；模型路徑以 `NGS_UI_PRS_CONFIG` 中的 `multi_model_dir` 為準。CHD 原檔須符合 repo lock 的解壓後 hash，不能拿新下載但內容不同的檔案冒充舊來源。

## 二、來源與保留模型核對

```bash
.venv/bin/python -m research.prs.validation sources \
  --score-dir "$PRS_SCORES" --chd "$PRS_CHD" \
  --output "$PRS_AUDIT/sources.json"

.venv/bin/python -m research.prs.validation model \
  --score-dir "$PRS_SCORES" --chd "$PRS_CHD" \
  --model-dir "$PRS_MODELS" --output "$PRS_AUDIT/models.json"
```

如需重做全部癌症原始／harmonized 對照，在 `sources` 加 `--cancer-original /path/to/PGS004875.txt.gz`。來源檢查對五病比 compressed SHA256、CHD 比 decompressed SHA256；不同內容直接拒絕。

`model` 先驗證既有 manifest 與 artifacts hash，再以原始 PGS 的列身分、座標、allele、beta 核對保留模型。反向互補另計數，不把它當成整個分數取負號。通過表示來源對應與數值定義相符，不是獨立基因型或疾病效能通過。

## 三、獨立基因型基準與缺失敏感度

輸入必須是另一條來源的 **GRCh38、單樣本、雙等位 SNV VCF**，有明確 diploid GT；可先從 reference PGEN 匯出指定公開樣本與保留位點，並確認樣本真正在參考中。參考基因型與院內 DRAGEN gVCF 的來源關係、caller、已知可靠範圍仍要記錄。

只有變異列的普通 VCF 通常沒有所有 0/0，不能因此宣稱完整覆蓋。工具對缺席／缺失／非 PASS 或 `.` 以外 FILTER 的位點視為未知；不讀 gVCF reference blocks，也不把未出現的列當 0/0。部分真值只輸出部分加總及涵蓋數，不給全模型百分位。

```bash
.venv/bin/python -m research.prs.validation benchmark \
  --model-dir "$PRS_MODELS" --score-dir "$PRS_SCORES" --chd "$PRS_CHD" \
  --vcf /path/to/independent_single_sample_full_sites.vcf.gz \
  --mask /path/to/case/missing_loci.tsv.gz \
  --eas-frequency "$PRS_AUDIT/eas.afreq" \
  --expected /path/to/case/multi/summary.json \
  --output "$PRS_AUDIT/benchmark_01.json"
```

`--mask`、`--eas-frequency`、`--expected` 均可省略。遮罩採既有 `missing_loci.tsv.gz` 的 `disease`／`id`；`--expected` 必須是對應**同一個人、同一保留模型與相同遮罩**的六病 summary。若用某個院內個案的遮罩模擬另一名參考者，不應提供該院內個案的 summary 作數值真值。

EAS 頻率以同一份 reference PGEN、固定 EAS735 keep 另算，不能用已按個案裁切的頻率。可用已設定的 PLINK 容器執行以下參數；輸出在稽核目錄，不覆寫原頻率：

```text
plink2 --pfile <multi_model_dir>/reference
       --keep <multi_model_dir>/reference_calibration_eas_735.keep
       --freq --out <PRS_AUDIT>/eas
```

工具核對頻率檔的 REF／ALT、頻率範圍與有效 allele count 不超過 1,470；**這些欄位檢查不能獨自證明它真由固定 735 人產生**，需保存上述 keep 的 hash、PLINK Log／命令與參考版本。缺任一保留 allele 的頻率會拒絕，不以零補缺值。

有完整真值時輸出：未遮罩、目前 3,330 頻率補值、EAS 頻率補值的分數；以及相對既有 EAS735 的排名位移、分數標準差位移與跨組情況。這是在固定 ECDF 上隔離補值影響，並未建立可啟用的新補值模型。

提供 expected 且完整覆蓋時，絕對分數差 `≤1e−5` 僅作數值對帳；不相符會標 `FAIL` 並非零退出。部分真值標 `PARTIAL`；沒有 expected 則明記只計算、未對照平台。產生檔案的成功退出，不表示台灣臨床驗證通過。

## 四、多次比較的彙整

```bash
.venv/bin/python -m research.prs.validation aggregate \
  --inputs "$PRS_AUDIT/benchmark_01.json" "$PRS_AUDIT/benchmark_02.json" \
  --output "$PRS_AUDIT/aggregate.json"
```

彙整同一模型版本，拒絕相同 VCF／遮罩／頻率 hash 的重複比較。只把完整真值納入百分位統計，另列部分覆蓋及數值不相符的數量。每病輸出平均有方向位移、絕對位移中位數／最大值及跨組次數。

比較檔數不等於獨立人數：同一個人套多個遮罩會得到多次比較。正式族群結論仍需標明抽樣、親緣、重複使用參考者、頻率估計相依性、尾部與批次。這一輪不由這些數字自動更換 QC gate。

## 五、平台工具追溯及相容性

新六病 score／ancestry 工作會各產生一份 `provenance.*.json`，由 worker 發布到既有 `08_postprocessing/{sid}.prs/runs/{job}/score` 或 ancestry artifacts。內容包括：

- gVCF 可取得的 source／版本／reference header；版本不明保持 unknown，不從檔名推測。
- 執行環境中工具路徑、版本查詢輸出與成功／缺失狀態。
- 計分 image／Nextflow config hash、科學程式 hash、可取得的 repo commit。
- Nextflow 請求版本、cached checkout commit、當次 launch revision 與 Log 中可取得的容器 hash；Log 沒列出的容器不宣稱已全數盤點。

有完整 checkpoint 的續跑不覆寫原 provenance；失敗後新的嘗試會記錄當次檢查。來源 fingerprint、完成檔 checksum、祖源 gate 與風險包鎖定維持；追溯蒐集失敗不假造版本。這次 runner 內容變動會由既有 fingerprint 安排新的分析 run，舊已發布分數及風險快照保留。

新 UI／Word 提示可由既有已驗證快照衍生，不需重跑基因分析；舊快照不因顯示改變而重寫。此處沒有替操作者更新伺服器或啟停服務。
