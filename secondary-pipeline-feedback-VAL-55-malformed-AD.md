# 二級分析 Pipeline 問題回報：multiallelic VCF 的 FORMAT/AD 欄位未正確對齊

- 日期：2026-07-17
- 樣本：`VAL-55`
- 二級分析輸出：`VAL-55.ensemble.fixed.vcf.gz`
- 失敗階段：三級分析 `PREPARE_VCF:ADD_CALLERS_TAG`
- 失敗位置：`chr1:83829`

## 1. 問題摘要

二級分析產生的 ensemble VCF 在 multiallelic 位點中，`FORMAT/AD` 的欄位數量與 `REF + ALT` 數量不一致。

三級 pipeline 使用 `bcftools norm -m -any` 拆分 multiallelic variant 時，因無法判斷 AD 數值分別屬於哪一個 ALT，因此直接中止。

這不是 Nextflow、Java 或 Apptainer 的問題，而是輸入 VCF 不符合 `FORMAT/AD Number=R` 的格式要求。

## 2. 錯誤 Log

```text
ERROR ~ Error executing process > 'PREPARE_VCF:ADD_CALLERS_TAG (1)'

Caused by:
  Process `PREPARE_VCF:ADD_CALLERS_TAG (1)` terminated with an error exit status (255)

Command executed:

  bcftools norm -m -any \
      -f /home/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta \
      -c w \
      VAL-55.ensemble.fixed.vcf.gz \
      -Oz -o VAL-55.norm.vcf.gz

  python3 /home/pipeline/tertiary_code/scripts/add_callers_tag.py \
      --input VAL-55.norm.vcf.gz \
      --sample VAL-55 \
      --output VAL-55.callers_tagged.vcf

Command exit status:
  255

Command error:
  Error: wrong number of fields in FMT/AD at chr1:83829,
  expected 8, found 6. Use --force to proceed anyway.

Work dir:
  /home/n102968/NGS_UI/nf_work/VAL-55-WGS-nckuh/9e/520f9da39b00a6978c19ef3e604376
```

三級分析因此無法產生 SNV、Mito、CNV/SV、PGx 等後續輸出。

## 3. 問題 VCF Record

```text
#CHROM  POS    ID  REF                        ALT                                                                  QUAL    FILTER  INFO  FORMAT                VAL-55_DV                 VAL-55_HC
chr1   83829   .   GAGAAAGAAAGAAAGAAAGAAAGAA G,GAGAAAGAAAGAAAGAAAGAA,GAGAAAGAAAGAAAGAAAGAAAGAAAGAA 312.02  PASS    ...   GT:PS:AD:DP:GQ:PL    1|2:83829:.:.:.:.        1/3:.:0,6,2:8:84:329,84,104,.,.,.,245,0,.,239
```

此 record 有三個 ALT：

```text
ALT1 = G
ALT2 = GAGAAAGAAAGAAAGAAAGAA
ALT3 = GAGAAAGAAAGAAAGAAAGAAAGAAAGAA
```

因此總共有四個 allele：

```text
0 = REF
1 = ALT1
2 = ALT2
3 = ALT3
```

## 4. FORMAT/AD 的格式要求

VCF header 中的 AD 通常定義為：

```text
##FORMAT=<ID=AD,Number=R,Type=Integer,...>
```

`Number=R` 表示每個 sample 的 AD 必須按照以下順序，為 REF 與每個 ALT 各提供一個位置：

```text
REF, ALT1, ALT2, ALT3
```

這筆 variant 有三個 ALT，所以每個 sample 應有四個 AD 位置。

### VAL-55_DV

```text
GT = 1|2
AD = .
```

整個 AD 設為 `.` 是合法的完整缺值。

### VAL-55_HC

目前內容為：

```text
GT = 1/3
AD = 0,6,2
```

AD 只有三個值。依照 VCF 的位置規則，會被解讀成：

```text
REF  = 0
ALT1 = 6
ALT2 = 2
ALT3 = 缺少
```

但 genotype 是 `1/3`，表示 HC 實際呼叫的是 ALT1 與 ALT3，而不是 ALT2。

PL 欄位也支持 `1/3`：

```text
PL = 329,84,104,.,.,.,245,0,.,239
                              ^
                         1/3 的 PL 為 0
```

這表示 ensemble 合併前，HC 很可能只有：

```text
REF, ALT1, ALT3 = 0,6,2
```

在 DV 與 HC 的 ALT 合併後，程式雖然已把 ALT2 插入 ALT list，也正確重排了 GT 與 PL，卻沒有在 HC 的 AD 中為 ALT2 插入缺值位置。

在確認原始 HC VCF 的 allele 順序後，預期正確內容應為：

```text
GT = 1/3
AD = 0,6,.,2
```

也就是：

```text
REF  = 0
ALT1 = 6
ALT2 = .
ALT3 = 2
```

請勿直接補成：

```text
0,6,0,2
```

`0` 表示該 allele 經過評估且讀數為零；`.` 才表示 HC caller 沒有提供該 allele 的 AD。

## 5. 為什麼錯誤訊息是 expected 8, found 6

此 record 有：

```text
3 ALT + 1 REF = 每個 sample 應有 4 個 AD 位置
```

又有兩個 sample：

```text
VAL-55_DV
VAL-55_HC
```

所以 bcftools 預期的內部欄位總數為：

```text
2 samples × 4 AD positions = expected 8
```

目前 record 的 AD vector 寬度只有三個位置，因此為：

```text
2 samples × 3 AD positions = found 6
```

所以錯誤訊息顯示：

```text
expected 8, found 6
```

## 6. 推測的根本原因

二級 pipeline 在合併 DeepVariant 與 HaplotypeCaller 的 multiallelic records 時，建立了所有 caller 的 ALT union，但沒有同步依新的 ALT index 重排 `FORMAT/AD`。

目前看起來：

- `GT` 已正確轉換成新的 allele index。
- `PL` 已依新的 genotype 組合插入 `.`。
- `AD` 沒有依 ALT union 插入缺少的 allele 位置。
- 可能還需要檢查其他 `Number=A`、`Number=R`、`Number=G` 欄位是否有相同問題。

## 7. 建議的正式修復方式

### 7.1 修正 ensemble 合併流程

合併 DV 與 HC allele 時，先建立完整 ALT union，例如：

```text
Combined ALT:
ALT1, ALT2, ALT3
```

接著針對每個 sample 建立「原始 allele index → combined allele index」對照，並同步重排：

- `GT`
- `AD`（`Number=R`）
- `PL`（`Number=G`）
- 其他 `Number=A` 或 `Number=R` 的 FORMAT/INFO fields

對於 caller 未提供的 allele：

```text
使用 .
```

如果無法可靠判斷原始 AD 與 ALT 的對應關係，應將該 sample 的完整 AD 設為：

```text
AD=.
```

不應保留長度不足的 AD vector，也不應猜測或任意補零。

### 7.2 優先使用理解 VCF cardinality 的工具

如果目前是以文字處理方式合併 VCF，建議改用能依 VCF header 處理 `Number=A/R/G` 的工具，或在自訂合併程式中明確實作 allele remapping。

不要只更新：

```text
ALT
GT
PL
```

而漏掉 `AD`。

### 7.3 重新產生輸出

修正後應重新產生：

```text
VAL-55.ensemble.fixed.vcf.gz
VAL-55.ensemble.fixed.vcf.gz.tbi
```

建議從原始 DV/HC VCF 重新執行 ensemble 合併，不要只手動修改單一 record，因為其他 multiallelic 位點可能也有相同問題。

## 8. 建議增加發布前檢查

在二級 pipeline 發布 `ensemble.fixed.vcf.gz` 前，增加完整 VCF normalization preflight：

```bash
bcftools norm \
    -m -any \
    -f /home/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta \
    -c w \
    -Ou \
    -o /dev/null \
    VAL-55.ensemble.fixed.vcf.gz
```

此檢查必須在不使用 `--force` 的情況下成功。

建議讓二級 pipeline 在發現以下錯誤時停止發布：

```text
wrong number of fields in FMT/AD
Incorrect number of fields
```

如此可在二級分析階段提早攔截，避免三級 pipeline 啟動多個工作後才失敗。

## 9. 不建議作為正式修復的方式

三級 pipeline 可以暫時加入：

```bash
bcftools norm -m -any --force ...
```

但 `bcftools norm --force` 對 malformed `Number=R` field 的處理方式是丟棄無法解析的 tag，而不是修復 allele mapping。

這可能造成：

- 該位點失去 AD。
- 無法可靠計算 VAF。
- `add_callers_tag.py` 可能無法正確判斷 caller。
- UI 或後續報告缺少深度資訊。
- 問題被靜默掩蓋，其他 malformed records 繼續流入。

因此 `--force` 最多只能作為緊急繞過方案，不應取代二級 pipeline 的正式修復。

## 10. 驗收條件

修復完成後應確認：

1. 原始 HC VCF 的 allele 順序證實該位點 AD 應為：

   ```text
   0,6,.,2
   ```

2. `VAL-55.ensemble.fixed.vcf.gz` 可在不使用 `--force` 的情況下通過：

   ```bash
   bcftools norm -m -any
   ```

3. `chr1:83829` 拆分後，每個 ALT 取得正確的 AD。

4. GT、AD 與 PL 的 allele index 一致。

5. 重新掃描整份 WGS VCF，確認沒有其他 multiallelic record 發生相同的 `Number=A/R/G` 欄位數錯誤。

6. 重新建立 tabix index。

7. 三級 pipeline 使用 `-resume` 後可順利完成 `PREPARE_VCF:ADD_CALLERS_TAG`。
