# Specialized caller 資料包與操作

2026-10-16；DRAGEN v4.3 targeted schema。實作位於 `backend/app/services/specialized.py`，由三級 worker 的 staging post-processing 呼叫。部署由使用者執行，不需要修改外部 Nextflow pipeline。

## 資料流與檔案契約

```text
同一 source sample 的 targeted JSON + targeted / hard-filtered VCF
  → 三級 post-processing：樣本驗證、索引區域擷取、同源區 evidence、allele projection
  → staging/08_postprocessing/{LIS}.specialized.*
  → 既有 managed promotion / rollback
  → 平台 / gene search / secondary / 摘要 / 報告：只讀三級結果
```

| 檔案 | 用途 |
| --- | --- |
| `{LIS}.specialized.manifest.json` | schema/parser/caller 版本、source sample、assembly、generation、來源追溯、本地各檔 SHA256、locus 狀態；最後發布 |
| `{LIS}.specialized.source.json` | 完整小型 targeted JSON，包括尚未投影的其他 caller；保留 null、原始欄位與版本 |
| `{LIS}.specialized.evidence.jsonl.gz` | 各來源 VCF header 與四 locus 區域原始紀錄，保留 ref/no-call、同源區、multiallelic、FILTER、GT、ALLELE_ID、JIDS、joint FORMAT |
| `{LIS}.specialized.calls.json` | stable findings、組成、copy number、review warnings/token、locus summary、被整合的 component IDs |

JSON 不裁掉其他 caller，方便日後擴充；本版 VCF 區域快照只涵蓋 GBA1、CYP21A2、SMN、HBA，其他 caller 的 VCF 若未納入，日後需另做擷取，不能宣稱已有完整證據。

只在三級 post-processing 或明確指定 `--from-source` 時讀來源資料夾。原始路徑僅作 provenance，loader 不對其 stat/open；個案清單的原始 VCF 大小改讀 `pipeline_source.json` 中保存的 `source_vcf_size`，舊 metadata 無該欄則顯示未知。IGV 仍依原機制讀 BAM。

大型 hard-filtered VCF（≥100,000,000 bytes）須有 bcftools 與 tbi/csi，僅查目標及同源區。小型 VCF 可串流讀取，過程有內容上限。來源 JSON、單一輸出、資料包總量及 evidence 解壓內容均限制在 100 MB 以下；超限報錯，不截斷、不複製完整大型 VCF。來源錯樣本、JSON 衝突、格式錯誤或不支援 assembly 時不發布；來源本來缺 caller 則以 `not_available` 發布，不能當陰性。checksum／generation 不符時拒絕載入，不 fallback 過時結果。

VCF 兩來源對同 ALT 的 evidence 優先取較完整 joint metrics，原始兩筆仍保存；不同 GT／copy count 保留衝突。一般 SNV 的完整註解仍在 03，資料包只保存對應少量 allele 的既有註解；複合 allele 不繼承單一 constituent 的分類。一般單點沿用 genomic ID；合成 findings 的 ID 依 assembly 與生物學身分計算，與路徑、時間、evidence 深度無關。GBA1 原 ID 不變。

## 判讀與呈現

所有 findings 沿用卡片、Comment、人工 ACMG、1/2/C 標記及報告選取。額外摘要保留正常／未判定／未提供的 locus；「無已命名重組」不等同沒有其他變異。未知 allele 不自動改為 P/LP，聯合 reads 不轉成 gene-specific DP/VAF。

| Locus | 處理方式與限制 |
| --- | --- |
| GBA1 | 沿用既有七種 haplotype 與一般 variants、RecNciI whole-allele curation；保存 joint metrics，取消載入時讀原始 VCF補證據 |
| CYP21A2 | 每條 `+` 串接的 caller haplotype 為一筆，保留 constituent HGVS；非重組小變異獨立處理。combined CN 包含 CYP21A1P，不能直接當 CYP21A2 CN；異常劑量、缺 component evidence、FILTER、定位／來源衝突需複核 |
| SMN | 正常也顯示 SMN1、SMN2、SMN2Δ7–8 CN。SMN1 CN0/1 建立拷貝數卡，不偽造 SNV、斷點或 het/hom；marker 單獨不能確診 silent carrier，SMN1 CN≥2 不排除 2+0。SMN2 CN不當 SMN1 疾病 allele |
| HBA | 一張完整兩 haplotype genotype 卡，小變異另存；`--` 不推定 SEA，保留所有已支援 3.7/4.2 缺失、triplication 組合與原 FILTER/GQ。未知 genotype 不以總 CN 猜組合；小變異 copy number 不直接轉成 het/hom |

review warnings 可經獨立驗證釐清時，需勾選卡片確認並在 Comment 記錄依據；確認綁定 evidence token，重解析證據改變後失效。未知定義、no-call 等 hard-blocked findings 及單獨 SMN marker 不能透過 checkbox 放行為疾病結果。

## 報告

保留既有診斷八欄／健檢七欄表格、表下列點與附錄。組成說明放第 2 點並於附錄重複；RS ID／結構缺值為全形 `－`。HBA1/HBA2 長 gene label 換行，避免與 RS ID 黏在一起。CN 卡用實際拷貝數或 genotype，不使用 heterozygous SNV 句型。健檢只有實際回報專一性結果時才補充檢測範圍例外，不宣稱涵蓋該基因所有變異。

- GBA1、CYP21A2：一條 het P/LP allele 第二類；兩個獨立 P/LP、相位未明仍第一類並保留相位說明。同一 haplotype 的成分只算一個 allele。
- SMN1：CN1 第二類、CN0 第一類；劑量結果須先完成人工 P/LP 分類。單獨 marker 不列為確診 carrier。
- HBA：常見一／二個功能拷貝缺失組合第二類，HbH／Hb Bart's 對應 genotype 第一類；`--/aaa3.7`、`--/aaa4.2` 仍依完整組合為 carrier。小變異、triplication／複合結果需人工指定健檢第一／第二類、當前 evidence token、Comment 與 P/LP 分類。同一 locus 多筆選取的分類必須一致，不能只數 variant rows 判定。
- 帶因篩檢既有 `carrier_mackenzie_1300+` 清單補入 HBA1、HBA2，供 HBA 結果選取與報告檢測範圍一致呈現；CYP21A2、SMN1 原已在清單。其他疾病 panel 按既有 gene membership，不把四個 locus 加入 ACMG SF 清單。
- caller 未完整判定且在該健檢 panel 範圍內，報告加「不能據此視為陰性」的限制。
- 過時 finding、已整合 component 或 haplotype 與 constituent 同時被選取時拒絕匯出，要求重新核對；不自動搬移舊選取。

## 舊個案補建與離線重解析

部署後，尚無資料包的個案依原本 `pipeline_source.json` 精確定位來源。請使用與服務相同的 Python 環境與 NGS_UI_* 設定，在 repo 執行：

```bash
# 舊個案：先驗證，讀原始 targeted JSON / VCF，但不發布
python scripts/rebuild_specialized.py --sample 26T00054 --from-source
# 確認摘要後發布；不重跑 Nextflow / VEP
python scripts/rebuild_specialized.py --sample 26T00054 --from-source --apply

# 已有資料包：只讀三級快照，原始來源離線也可執行
python scripts/rebuild_specialized.py --sample 26T00054
python scripts/rebuild_specialized.py --sample 26T00054 --apply
```

來源搬遷時可明確給 `--from-source --source-vcf /new/run/vcf.gz/SOURCE.hard-filtered.vcf.gz`；仍驗證 JSON sampleId、VCF header 與既存 source sample 相同。缺少原始 JSON 的舊個案不能憑空補出不存在的證據；舊 `gba.json` 仍可讀，但不會再向原始 VCF 補數字。

補建共用三級 sample lock，遇到 active job 拒絕競爭；先在 08 隱藏暫存目錄完成，再逐檔 promotion、manifest 最後。失敗 rollback 所有已移動檔案；程序硬中斷造成不完整檔案時 checksum/generation 會拒絕讀取，重新補建即可。其餘 00–07 與 metadata 的 review fields 不變；完整三級重跑則採既有 worker promotion。載入快取依本地四檔簽章失效。

## 驗證與適用範圍

`tests/test_specialized.py` 包含來源消失／資料包搬遷後讀取與離線重解析、sample mismatch、checksum、大小上限、multiallelic／來源衝突、CYP cis／component 防重複、SMN CN／marker、HBA genotype／人工分類、過時 token、報告模板與 publish rollback。前端測試涵蓋 joint evidence、CN0、HTML escaping、未知結果與 IGV 使用本地 evidence。

以提供的 26T00054（DRAGEN v4.3.16）實際 JSON 與 indexed hard-filtered VCF 驗證：四 locus 共約 62 KB，保留 RecNciI joint metrics。CYP21A2／SMN／HBA 異常組合使用合成 fixtures；這是軟體流程驗證，不取代實驗室對各 caller／變異類型的分析效能驗證。DOCX 使用既有模板渲染檢查；本機預覽以替代中文字型呈現，正式報告仍保留 MingLiU。

主要 schema 與判讀來源：

- [DRAGEN v4.3 CYP21A2 Caller](https://help.dragen.illumina.com/dragen-v4.3/product-guide/dragen-v4.3/dragen-dna-pipeline/targeted-caller/cyp21a2-calling)
- [DRAGEN v4.3 SMN Caller](https://help.dragen.illumina.com/dragen-v4.3/product-guide/dragen-v4.3/dragen-dna-pipeline/targeted-caller/smn-calling)
- [DRAGEN v4.3 HBA Caller](https://help.dragen.illumina.com/dragen-v4.3/product-guide/dragen-v4.3/dragen-dna-pipeline/targeted-caller/hba-calling)
- [GBA1 現有 curation 與規則](GBA1_SPECIALIZED_CALLER.md)
