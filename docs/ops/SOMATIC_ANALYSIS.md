# 指定範圍 Somatic SNV/Indel

## 操作與合併規則

載入個案後，SNV/Indel 標題旁「Somatic 分析」開啟 modal。支援多基因（空白、逗號、分號或換行分隔）、`chr1:100000`、`chr1:100000-101000` 或指定 allele `chr1:100000:A>G`。以上僅為格式範例。座標固定 GRCh38、1-based inclusive；基因／座標取聯集。預設所有 transcript 的 exon 聯集加兩側各 20 bp，也可選完整 gene span。可預覽合併後的確切區域。

BAM 只允許目前個案的 IGV resolver 找到且有 index 的檔案。卡片 IGV 使用本次實際選取的 BAM。Modal 可關閉，背景工作繼續；重新開啟可看進度、log、取消、被過濾候選和重新執行。

完整 germline raw TSV 和 Mutect2 結果以相同 FASTA 做 bcftools normalization，再依 CHROM/POS/REF/ALT 排除既有 allele。比對不套 germline DP/VAF/BED/AF/主畫面篩選；`CALLERS=NONE`、reference rows 不算 germline ALT。同座標不同 ALT 保留。不能先從 calling interval 排除 germline 座標，以免漏掉不同 ALT。原始 Mutect2 VCF 保留稽核，送往 annotation 的 `novel.vcf` 已排除 germline 重複點位。

原 germline 卡片不加 somatic 標籤，不覆寫 AD/DP/VAF。僅新增卡片顯示「Somatic pipeline」、Mutect2 read support、FILTER、品質資訊及人工驗證狀態，沿用 SNV tiers、reviewer 標記、comment、transcript、搜尋與報告。來源本身不提高 tier，也不自動進 secondary findings 或選入報告。

「☑ Somatic」只在成功完成的 run（含零候選）後出現，預設勾選。它只控制新增卡片，不放寬 germline 的篩選。Somatic 候選不受既有 VAF、疾病基因、HPO/panel、MODIFIER、local-common 顯示限制；使用自己的 caller/filter 結果。PASS 預設納入，其他 FILTER 須在 modal 明確「加入判讀」，不改成 PASS。`AS_FilterStatus` 逐 allele 處理，避免 site PASS 掩蓋個別 ALT 的失敗。

## 執行端設定

本版本由 UI server 啟動 detached Python worker，在**同一主機**執行工具；不需要 RQ、SSH 或重跑完整 Nextflow 00–07。原生工具或管理員配置的 Apptainer command prefix 均可。若工具/reference 僅存在 DGX，須先讓 worker 主機可執行並存取相同檔案；本版本不含跨主機排程。

必備：GATK 4（Mutect2、LearnReadOrientationModel、FilterMutectCalls）、bcftools、Samtools（支援 `depth -s`）、VEP offline cache、GRCh38 FASTA/FAI/dictionary、germline population-AF resource/index、版本固定的 ClinVar VCF、gene/exon BED。工具／容器版本應在院內固定並驗證；worker 保存版本輸出及 config/resource signatures。

1. 以相容的 GRCh38 GTF 建立五欄 BED（`chrom start0 end gene exon|gene`），保留所有 transcript：

   ```bash
   python scripts/build_somatic_gene_regions.py \
     --gtf /path/to/reference.GRCh38.gtf.gz \
     --release 'your-reference-release' \
     --out /path/to/somatic_gene_regions.bed
   ```

   Reference/BAM 使用 `chr1` 命名；gene alias 由現有 HGNC/panel canonicalizer 處理。

2. 複製 `deploy/somatic_config.example.json` 到 `NGS_UI_HOME/data/somatic_config.json`，填入實際路徑與 release；也可透過 `NGS_UI_SOMATIC_CONFIG` 指定檔案。`vep_cache_version` 是實際 numeric cache version。設定檔由管理員維護，request 不得提供命令。

   原生工具範例：`"gatk_command": ["/path/to/gatk"]`。

   容器範例：`"gatk_command": ["apptainer", "exec", "--bind", "/home:/home", "/path/to/gatk.sif", "gatk"]`。其他工具亦使用 JSON argv 陣列，不經 shell。所有 input/reference/staging/output 路徑須在容器內以相同路徑可見。

3. `pon`、`contamination_sites` 可選；缺少會在完成紀錄提示。PoN 需與平台／建庫方式相容，不用院內 germline AF 取代。污染估計使用背景 resource loci，不限於指定的小基因區域。預設同時一個 somatic job、最多 200 genes／200 positions／20 Mb，concurrency 與總長度可由設定檔調整。

4. 部署程式並重啟 UI server。缺設定／資源／工具時 modal 明示原因並停用啟動；既有 germline 功能可正常使用。

5. 在院內用已確認陽性、陰性／artifact 及稀釋資料驗證 BAM/reference、SNV/indel、低 VAF、重跑與報告。軟體及合成資料測試不代表生物學或臨床檢出極限已完成驗證。

## Pipeline

`workers/somatic_run.py`：

1. 檢查 BAM、單一 read-group sample、reference contig lengths、資源與輸入穩定性。
2. Mutect2 單檢體模式，指定 interval 加 100 bp assembly padding，可 force-call 指定 allele。
3. LearnReadOrientationModel；選配 GetPileupSummaries/CalculateContamination；FilterMutectCalls。
4. bcftools norm 拆 ALT/left-align；完整 germline 同法正規化後去重，結果依原始指定區域收錄。
5. VEP offline JSON（固定 cache version、所有 transcript）、固定 ClinVar、本地 GeneBe（存在時；不呼叫 live API）。UI 共用 HPO/panel、OMIM/gene-disease、有效 ACMG overlay；未配置的額外 germline predictor/plugin 不假造分數或 tier evidence。
6. Samtools depth（BQ/MQ ≥20、paired overlap 不重複）與指定座標覆蓋摘要。零候選仍完成，無 call 不表示排除變異，coverage 不是 validated LOD。
7. 檢查 annotation 候選完整、input/resource 未變，複製至目的 filesystem 隱藏目錄後 rename，最後原子發布 index。

## 儲存、重跑與報告

```text
data/jobs/somatic/{run_id}/
  state.json, config.json, log.txt, cancel
tertiary_output/{LIS_ID}/09_somatic/{run_id}/
  mutect2.vcf.gz, filtered.vcf.gz, novel.vcf
  annotations.tsv, candidates.json, coverage.json, manifest.json
  targets.bed, tool-version files, orientation/contamination evidence
08_postprocessing/{LIS_ID}.somatic.json
```

Legacy 個案的 index 由現有 state resolver 放在原 UI state 目錄。每 run 獨立，不修改 germline 00–07。暫存的完整 germline VCF/normalization 檔最後刪除，不保留第二份完整 raw。取消／失敗不切換 index，共用三級 sample lock；active job 阻擋取消登錄和刪 output。

不同範圍累積，重新執行某次 run 成功後 archive 前次；已標記 1/2/C 的歷史點位保留顯示。相同 variant 使用較新有效 observation，不累加 DP/AD。Germline raw signature 改變後停止載入 stale somatic、提示重跑；已標記點位因此缺失時阻擋診斷 DOCX，避免靜默漏報。

診斷 DOCX 對人工選取的 somatic-only 變異附來源、VAF、FILTER、驗證狀態、該次 ClinVar release；germline 報告不變。完整刪 pipeline 清除 somatic jobs/results，單純取消登錄保留結果。

## 軟體驗證

`tests/test_somatic.py`：區域、別名、多 ALT、低 VAF/DP、權限、FILTER 明確納入、stale、歷史、DOCX，以及使用合成工具輸出的 worker 全流程。`tests/frontend_somatic.test.cjs`：checkbox 顯示條件與 germline/somatic 篩選隔離。另與既有 layout、case summary、secondary finding、三級刪除及報告測試一起執行。

參考：[Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/21905083931035-Mutect2)、[GATK somatic workflow](https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect2.wdl)、[VEP 格式](https://www.ensembl.org/info/docs/tools/vep/vep_formats.html)。
