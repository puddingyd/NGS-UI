# NGS 二級／三級分析 Pipeline — 改動總覽（2026-07）

這一輪從評鑑回饋、Broad 對齊、depth／phasing 修復，到 sub-workflow 重構與 DRAGEN 交叉檢核的完整改動清單，依主題分組。

- **Repo**：`NGSSecondaryAnalysis`（FASTQ → 變異）、`NGSTertiaryAnalysis`（變異 → 註釋／分類）
- **開發分支**：`claude/ngs-quality-compression-analysis-r7j5wp`
- **標記**：每項標示涉及的 repo（二級／三級／共用）與類型（新增／修復／重構／對齊·調整／文件／驗證）

---

## A. 主流程整合與旗標

1. 退役 `main_research.nf`，研究工具整合進單一 `main.nf`，改用 flag 啟動：`--run_manta`、`--run_expansionhunter`、`--run_automap`（非商用，預設 OFF）與 `--run_roh`（bcftools roh，可商用，預設 ON）。　`二級 · 重構`
2. `--run_phasing` 與 WES 的 `--run_gcnv` 改為**預設 ON**（二級＋三級）；`--run_gcnv false` 保留為「PON 尚未建好」的逃生口。　`二級／三級 · 調整`

## B. CNV／SV 與 Broad germline-CNV 對齊（需重建 PON）

3. Delly 只發布 `FILTER=PASS`（單樣本不能用 `delly filter -f germline`），砍掉 LowQual 噪音；臨床相關事件確認仍在 PASS 集內。　`二級 · 對齊`
4. CNVkit 移除 `--filter cn`（對 germline 過度激進、會漏真實 CNV），改由 `params.cnvkit_filter_cn` 控制、預設 OFF。　`二級 · 對齊`
5. gCNV 超參數對齊 Broad germline-CNV WDL 預設（`p-alt 5e-4`、`coherence 10000`、`p-active 1e-1`），並參數化到 `nextflow_pon.config`；改動後**必須重建 PON**（`gcnv_model/`、`cnvkit_reference/`）。　`二級 · 對齊`
6. 其他 Broad 對齊：VQSR SNP 加 `-an DP`；CPU HaplotypeCaller 改先 `ApplyBQSR`（GATK4 已移除 on-the-fly BQSR）；mtDNA 用 `--mitochondria-mode` ＋ blacklist mask；bwa 加 `-K 100000000`（thread-deterministic）。　`二級 · 對齊`

## C. SNV／indel：phasing、compound 合併、depth 修復

7. 新增 `--run_phasing`：各 caller 先 `whatshap phase` → `combine_phased.py` 把相鄰／重疊的 cis 變異合成單一 MNV（如 SUZ12 `c.2168_2170delAAAinsTT`），再進 ensemble，讓三級 VEP 報出正確的 combined `p.`。　`二級 · 新增`
8. **修復合併記錄 depth 消失**：combined MNV 改為繼承 anchor（cluster 內最寬 biallelic 記錄）的完整 FORMAT，保留 `AD`／`DP`／`VAF`，只覆寫 `GT`／`PS`；四種無法安全重建的情況直接 passthrough。修掉先前只剩 `GT:PS`、depth 掉成 `.` 的 bug（DRAGEN VAL-10 145k+ 筆）。　`共用 · 修復`
9. ensemble 前調和 `FORMAT/AD`→`Number=R`、`PL`→`Number=G` header，並加發布前 preflight，避免 merge／三級因欄位數不符崩潰；`combine_phased.py` 二級／三級逐位元組同步（md5 一致）。　`共用 · 修復`
10. 新增 haploid combine 模式（男性非 PAR `chrX`／`chrY` hemizygous compound）；`combine_phased.py` 改為 staged path input，修掉 `-resume` 讀到舊快取的陷阱。　`共用 · 新增`

## D. 性別／倍體 QC

11. 新增 `ploidy_check.py`（mosdepth-based）：輸出 `ploidy.vcf.gz`（`DC:NDC:RATIO`，header key 對齊 DRAGEN `##estimatedSexKaryotype` 等）＋人可讀 `ploidy_qc.txt`；warn-only 的 sex 比對與非整倍體提示，不改 ploidy、不讓流程失敗。　`二級 · 新增`
12. **修復 WGS 性染色體 QC**：mosdepth 對 `--by` BED 外的 contig 會吐 `*_region=0`，使男生被誤判成 `X0?`／MISMATCH；改為 `*_region` 有值且 > 0 才採用、否則回退整條 contig mean。輸出限縮到主要 contig（chr1-22, X, Y, M）。　`二級 · 修復`
13. 三級接上 DRAGEN 原生 ploidy（`PLOIDY_REPORT_DRAGEN` ＋ `parse_dragen_ploidy.py`）→ 產生與二級統一格式的 `ploidy_qc.txt`，NDC 正規化語意一致。　`三級 · 新增`
14. `sex_ploidy_GRCh38.txt` 收為單一真相來源、移到 `ref_dir`（`bcftools +fixploidy` 的 ploidy map）。　`二級 · 調整`

## E. 三級 annotation 強化

15. SNV／indel 表新增 `STRAND_BIAS` 欄（`parse_vep_csq.py`）：依 FisherStrand／SOR 標 `PASS`／`WARN(FS=..,SOR=..)`；DeepVariant-only 位點無 FS/SOR → `.`（人工複核）。germline 只標記不硬刪。　`三級 · 新增`
16. 新增 DRAGEN PGx 交叉註記（`compare_dragen_pgx.py` ＋ `PGX_DRAGEN_CONCORDANCE`）：DRAGEN 原生 `targeted.json` 對照我們的 `pgx.tsv`，在 `NOTES` 標 `DRAGEN 一致／不一致／未比對` ＋原始 genotype（reference 跨寫法正規化、模糊多重解落在候選集即一致；欄位不變；找不到 json 則跳過不報錯）。　`三級 · 新增`

## F. 架構重構（sub-workflow 化）

17. 二級 `main.nf` 全面 sub-workflow 化：`ALIGNMENT_QC`／`CALL_SNV`（內含 DEEPVARIANT、HAPLOTYPECALLER、VQSR）／`CALL_CNV_SV`（內含 DELLY、GCNV）／`CALL_STR`／`CALL_MITO`／`CALL_ROH`；workflow body 變成純組合。　`二級 · 重構`
18. 三級 `main_tertiary.nf` sub-workflow 化：`ANNOTATE_SNV`、`ANNOTATE_CNV_SV_{NCKUH,DRAGEN}`、`ANNOTATE_STR_{NCKUH,DRAGEN}`；DRAGEN PGx 交叉註記收進 `PGX_ANNOTATE` 內部。　`三級 · 重構`
19. ExpansionHunter 輸出改名 `*.expansionhunter.*`，避免與 GangSTR 的 `*.str.vcf` 在 `06_repeat` 互相覆蓋（同名 clobber）。　`二級 · 修復`

## G. 文件與驗證

20. 兩 repo 的公開 `README.md`（架構圖＋每個輸出資料夾檔案／欄位說明）、三級開發筆記 `readme.md`（v3.5）與同事使用說明 `DGM_tertiary_pipeline_guide`（v3.5）全部同步更新。　`共用 · 文件`
21. 端到端驗證：二級 NA12878 WES ＋ VAL55 WGS；三級 NA12878 WES／WGS（NCKUH）＋ VAL-10（DRAGEN）。確認 DRAGEN AD 保留（`5,940,465 / 5,940,563`）、性別 QC 正確（WGS 男→XY、WES 女→XX）、WES→gCNV／WGS→CNVkit 分流、PGx WES／WGS 分流、Delly PASS-only。　`共用 · 驗證`

---

共 21 項改動，橫跨 `NGSSecondaryAnalysis`（FASTQ→變異）與 `NGSTertiaryAnalysis`（變異→註釋／分類）。預設路徑僅使用可商用工具；非商用工具（Manta／ExpansionHunter／AutoMap）需以 flag 明確啟用。
