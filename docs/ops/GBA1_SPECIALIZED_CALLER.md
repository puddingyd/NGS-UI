# GBA1 specialized caller 整合

三級 DRAGEN worker 在 staging 的 `08_postprocessing/{LIS_ID}.gba.json` 保存 caller 證據與 allele projection，隨既有 promotion／rollback 原子發布；不改寫 00–07，不另建 SNV 完整註解副本。NCKUH 與缺輸入的個案寫明 `not_available`，重跑不會殘留前一次 GBA 結果。舊個案需重跑三級（Nextflow 仍用 resume），才會產生此檔。

## 輸入與樣本身分

- 以原始 hard-filtered VCF 的 source sample 為準，精確找同層 `{source}.targeted.json`，或同 run 的 `other/{source}/germline_seq/{source}.targeted.json`。
- 同層 `{source}.targeted.vcf.gz` 與原始 hard-filtered VCF 均讀取 GBA 區域；相同 ALT 優先保留有完整 joint metrics／allele ID 的紀錄，其他只見於其中一份的 GBA1 allele 也保留。VCF 必須為單一相同 source sample；JSON sampleId 亦須相符。不同內容的多份候選 JSON、錯誤格式或樣本不符會阻止整批發布。
- 不使用 LIS ID 猜測來源、不跨 run 搜尋、不把其他樣本的 targeted VCF 合併。JSON SHA256、軟體／組裝版本、VCF 路徑與 size/mtime、原始 gba object、targeted records、FILTER、GT、ALLELE_ID、JIDS 與 joint FORMAT 留存。原始三級 TSV signature 改變時須重新產生衍生結果。
- 有 bcftools 與 `.tbi` 時只讀 GBA 區域；否則串流 VCF。JSON 缺失不代表陰性。`deletionBreakpointInGene=null` 顯示未知；CN 是 GBA1/GBAP1 combined CN。
- 舊 sidecar 若重組 allele 缺 JDP／JAD／JAF，loader 僅在既有原始 VCF 的 size／mtime 簽章相符、可用 bcftools 及索引時，唯讀補回對應 evidence 與 IGV 座標；不改 sidecar、不重算 allele 身分、分類或報告選取。若原始檔不在、已改變或無索引，須重做該樣本的 GBA1 post-processing。

## Allele 與判讀

支援 DRAGEN 的 A495P、L483P、D448H、c.1263del（1263del alias）、RecNciI、RecTL、c.1263del+RecTL，以及 `variants[]` 中的 alleleId／alleleCopyNumber／genotypeQuality／filter。未知 haplotype、no-call、非 PASS、無法唯一對應或 copy number 衝突會保留證據／待複核提示，不自動推定為正常或致病。沒有可確認組成的未知 allele 不允許正式匯出。

每個 recombinant haplotype 為一筆可標記的變異，同名兩份為 homozygous，兩個不同 caller haplotype 為兩筆。單一 recombinant-free 空字串只表示無已命名 recombinant haplotype，不能據此排除其他 SNV。非重組 variant 不限定白名單，以 VCF 的 allele-index 對應 ALLELE_ID 與 GT；ref-only 與 GBAP1 counterpart 不建立卡片。只取同 allele 的既有三級註解，不把一個 constituent 的 pathogenicity 套給整個多位點 allele。

複合 allele 的 stable ID 為 `GBA1-<24 hex>`，由 assembly／類型／名稱／deletion state 決定，可沿用既有 manual ACMG、review state 與 observation 機制。對應確定的單點保留 genomic ID。RecNciI 的全 allele 判讀參考 ClinVar VCV000004297.9（2026-10-15 核對，Pathogenic），來源標示為 GBA1 allele curation，不偽裝成七月 SNV baseline 或由 ACMG points 推算。其餘 allele 使用相符的單點既有註解或人工判讀，**不因 specialized caller 有 call 就全部指定 Pathogenic**。擴充全 allele 分類須有獨立可追溯證據。

主畫面、gene search、case summary、secondary findings 與 DOCX 共用投影；組成 marker 不重複計數。原先已選的組成位點不自動轉成父 allele，匯出時要求重新選取，避免把既有 component 註解轉移到完整 allele。Joint JDP／JAD／JAF 不代替 gene-specific 深度與 VAF；secondary 候選對合格 GBA caller 結果使用專屬品質判定，不捏造 0.5 VAF。GBA1 仍依原 panel gene membership 進 carrier screening，不加入 ACMG SF 清單。

卡片沿用 SNV 模板，以紫色組成說明取代標題的 opaque stable ID；下一行直接顯示 combined CN、deletion breakpoint、allele copies 與 FILTER。Read support 按組成位點順序列出 JDP／JAD／JAF 原值，游標提示其為 GBA1/GBAP1 合併證據，不能視為 GBA1 的 DP／AD／VAF；完整原始證據仍保存在 sidecar。複合 allele 不產生指向單個 marker 的 Varsome／Franklin／GeneBe 連結。有 JIDS 時 IGV 同時開啟 GBA1 與對應 GBAP1 區間，作為支持性複核，不能單憑三個 marker 證明 locus specificity。

## 報告

健檢仍用原有七欄表格、表下列點與附錄參考資料。RecNciI 的核苷酸欄顯示名稱，缺少 RS ID／結構時填全形 `－` 並保持表格對齊；全文組成在表下第 2 點及附錄各列一次，原本致病性與建議依序後移。所有健檢報告共用的陰性限制段落以「若」起句。診斷 DOCX 同樣保留原表格與增加組成列點。

健檢 GBA1 使用「高雪氏症（Gaucher disease）」及 AR 模式，單一 het allele 列第二類；兩個獨立 P/LP 變異即使相位未明仍列第一類，沿用原本相位分析說明。RecNciI 內部已知 cis 的三個成分只算一筆。沒有把 Parkinson susceptibility 改成 AD 診斷或 carrier 疾病名稱。完整 allele 的族群頻率沒有提供時明示缺值，不使用 component frequency，也不宣稱缺值表示罕見。

## 驗證與來源

`tests/test_gba.py` 涵蓋七種 haplotype、aliases、hom／het、multiallelic index、no-call、sample mismatch、component 去重、carrier／第一類、列點順序與舊選取；`tests/frontend_gba.test.cjs` 核對證據與雙區間 IGV。

- [DRAGEN v4.3 GBA Caller](https://help.dragen.illumina.com/dragen-v4.3/product-guide/dragen-v4.3/dragen-dna-pipeline/targeted-caller/gba-calling)
- [RecNciI whole haplotype — ClinVar 4297](https://www.ncbi.nlm.nih.gov/clinvar/variation/4297/)

Nextflow pipeline 本身位於伺服器外部 checkout；本次改動在本平台的 staging post-processing 層，無須修改外部 pipeline。部署由使用者執行，之後重跑個案三級以匯入原始 caller 結果。
