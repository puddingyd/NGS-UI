# Panel 內容修改彙整（2026-06-12）

本文件彙整「要把 panel 內容做更改的部分」，供各科回饋與來源 panel 檔（DGM 上之 `phenotype_data/gene_panels/*.txt`）更新之依據。

- 最終 reportable 基因清單：**7,363** 個 HGNC 現行符號
  （= curated 疾病關聯清單 6,240 ∪ 各科臨床 panel 對齊後新增 1,123）。
- 清單檔：`results/panel_alignment/panel_loose_plus_clinical.hgnc_canonical.txt`
- 各科「真正新增」基因（逐科）：`results/panel_alignment/clinical_panel_genuinely_new_by_panel.tsv`

---

## A. 符號誤植：請各科把來源 panel 改為現行 HGNC 名稱

下列為各科 panel 之誤植，已於 reportable 清單更正；**請各科一併更新來源 .txt**，
否則下次比對會再被標出。多數更正後之基因原本即在清單內（純命名問題）。

| 科別 / panel | 原符號 | 更正為 | HGNC | 備註 |
|---|---|---|---|---|
| 神經科 ALS-HSP | FTPS | PTS | HGNC:9689 | 已在清單 |
| 神經科 ALS-HSP | NIPAI | NIPA1 | HGNC:17043 | 已在清單 |
| 神經科 ALS-HSP | MT-T1 | MT-TI | HGNC:7488 | 已在清單（mt-tRNA）|
| 神經科 Dementia | EIF2B1-5 | EIF2B1 | HGNC:3257 | 已在清單 |
| 神經科 Dementia | LMNB | LMNB1 | HGNC:6637 | 已在清單 |
| 神經科 Stroke | LOC100505841 | ZNF475 | HGNC:53564 | **新增** |
| 神經科 Stroke | ND1 / ND4 / ND5 / ND6 | MT-ND1 / MT-ND4 / MT-ND5 / MT-ND6 | HGNC:7455… | 已在清單（均為粒線體基因）|
| 神經科 Stroke | PROC4 | PROC | HGNC:9451 | 已在清單 |
| 神經科 Epilepsy | DMN1 | DNM1 | HGNC:2972 | 已在清單 |
| 神經科 Epilepsy | PI3KCA | PIK3CA | HGNC:8975 | 已在清單 |
| 神經科 Epilepsy | SCL35A2 | SLC35A2 | HGNC:11022 | 已在清單 |
| 神經科 Epilepsy | SCN8 | SRP54 | HGNC:11301 | 已在清單 |
| 兒科 Growth_panel | HMAG2 | HMGA2 | — | 已在清單 |
| 兒科 Growth_panel | SMARCB2 | SMARCA2 | — | 已在清單 |

### A-2. 原為「染色體區域 / 基因簇」符號 → 改用具體致病基因

| 科別 / panel | 原符號 | 改為 | 說明 |
|---|---|---|---|
| 神經科 Epilepsy | DGS / VCFS | TBX1 | DiGeorge / 22q11；TBX1 已在清單，DGS、VCFS 移除 |
| 神經科 Epilepsy | HOXD | HOXD13 | HOXD 基因簇 → 具體致病基因 HOXD13（已在清單）|

---

## B. 請各科刪除（誤植，reportable 清單未納入）

| 科別 / panel | 刪除 | 原因 |
|---|---|---|
| 神經科 Neuropathy_NMJ_Myopathy / Dementia / Stroke | MT-CR | mitochondrial control region 誤植，非基因 |
| 神經科 Stroke | NR2D2 | 誤植 |
| 神經科 Epilepsy | PCDHG | 誤植（PCDHGA4/GA5/GA8/GC4 已在清單）|

---

## C. 待處理：染色體連鎖基因座（非 HGNC 基因，先擱置）

文獻把這些 band 列入 panel，但非 HGNC 基因、無 transcript／CDS，VEP 不會註解；
目前 CNV/SV 端未區分 dead zone，故**先不納入、僅記錄**。

| locus | NCBI Gene | 位置 (GRCh38) | 大小 |
|---|---|---|---|
| GFND1 | 100689213 | 1:198,700,001-214,400,000 (1q32) | ~15.7 Mb |
| HHT4 | 791087 | 7:28,800,001-43,300,000 (7p14) | ~14.5 Mb |
| MYMY4 | 100653379 | X:148,000,001-156,040,895 (Xq28) | ~8.0 Mb |

---

## D. 各科 panel 使用舊名／別名（建議來源更新為現行 HGNC，不影響 reportable 清單）

這些在比對時被標為「缺漏」，但其現行 HGNC 名稱其實已在清單內——只是來源 panel 寫了舊名。
建議各科一併更新（下次比對才會乾淨）。共 106 個。

（建議回饋各科把舊名改為現行名）

**`WES-II__兒科__Inborn_error_of_metabolism.txt`** （17）：
ODAD2 → ARMC4; POPDC1 → BVES; MTRFR → C12ORF65; CFAP418 → C8ORF37; VMA22 → CCDC115; CPAP → CENPJ; HYCC1 → FAM126A; G6PC1 → G6PC; GBA1 → GBA; BLTP1 → KIAA1109; ASPNAT → NAT8L; COXFA4 → NDUFA4; AFG2A → SPATA5; TAFAZZIN → TAZ; CRIPTO → TDGF1; VMA12 → TMEM199; IFT54 → TRAF3IP1

**`WES-II__兒科__遺傳性心臟疾病.txt`** （1）：
TAFAZZIN → TAZ

**`WES-II__神經科__ALS_HSP.txt`** （2）：
GBA1 → GBA; COXFA4 → NDUFA4

**`WES-II__神經科__Myopathy_CMS.txt`** （3）：
POPDC1 → BVES; G6PC1 → G6PC; TAFAZZIN → TAZ

**`WES-II__神經科__PD_dystonia.txt`** （4）：
FERRY3 → C12ORF4; GBA1 → GBA; COXFA4 → NDUFA4; AFG2A → SPATA5

**`WES-II__耳鼻喉科__syndromic__syndromic_hearing_loss.txt`** （1）：
HARS1 → HARS

**`WES-II__腫瘤醫學__心臟毒性.txt`** （1）：
TAFAZZIN → TAZ

**`WES-I__兒科__先天新陳代謝疾病.txt`** （1）：
GBA1 → GBA

**`WES-I__兒科__先天神經肌肉疾病.txt`** （1）：
POPDC1 → BVES

**`WES-I__兒科__遺傳性腎臟疾病.txt`** （2）：
PHB1 → PHB; NHERF1 → SLC9A3R1

**`WES-I__皮膚科__Ichthyosis.txt`** （1）：
GBA1 → GBA

**`WES-I__皮膚科__Inflammatory_genodermatoses.txt`** （1）：
RIGI → DDX58

**`WES-I__皮膚科__PPK_and_PC.txt`** （1）：
SACK1G → FAM83G

**`WES-I__眼科__Lens_disease.txt`** （3）：
HYCC1 → FAM126A; BLTP1 → KIAA1109; SLC9D1 → TMCO3

**`WES-I__眼科__Retina_disease.txt`** （3）：
MTRFR → C12ORF65; CFAP418 → C8ORF37; IFT38 → CLUAP1

**`WES-I__腫瘤醫學__遺傳癌症.txt`** （1）：
GBA1 → GBA

**`WGS__兒科__Early_onset_syndromic_epilepsy.txt`** （45）：
ADPRS → ADPRHL2; COA8 → APOPT1; ARSL → ARSE; ATP5F1A → ATP5A1; ATP5F1B → ATP5B; ATP5F1D → ATP5D; ATP5F1E → ATP5E; ATP5MC3 → ATP5G3; ATP5PO → ATP5O; MTRFR → C12ORF65; KICS2 → C12ORF66; MICOS13 → C19ORF70; VMA22 → CCDC115; CERT1 → COL4A3BP; DARS1 → DARS; FCSK → FUK; G6PC1 → G6PC; GARS1 → GARS; GBA1 → GBA; CBLIF → GIF; H3-3A → H3F3A; H3-3B → H3F3B; HJV → HFE2; CRPPA → ISPD; KARS1 → KARS; PRORP → KIAA0391; BLTP1 → KIAA1109; KIFBP → KIF1BP; LARS1 → LARS; MMUT → MUT; NARS1 → NARS; ASPNAT → NAT8L; COXFA4 → NDUFA4; EPRS1 → QARS; RARS1 → RARS; RNU2-2 → RNU2-2P; SARS1 → SARS; SKIC2 → SKIV2L; AFG2A → SPATA5; AFG2B → SPATA5L1; TAFAZZIN → TAZ; VMA12 → TMEM199; RXYLT1 → TMEM5; SKIC3 → TTC37; VARS1 → VARS

**`WGS__兒科__Growth_panel.txt`** （4）：
CPAP → CENPJ; G6PC1 → G6PC; GNAS → GNAS1; BLM → RECQL3

**`WGS__皮膚科__Inflammatory_genodermatoses.txt`** （1）：
RIGI → DDX58

**`WGS__皮膚科__PPK_and_PC.txt`** （1）：
SACK1G → FAM83G

**`WGS__神經科__ALS_HSP.txt`** （10）：
SPG21 → ACP33; ADAR → ADAR1; COQ8A → CABC1; DARS1 → DARS; RETREG1 → FAM134B; GBA1 → GBA; SPG11 → KIAA1840; COXFA4 → NDUFA4; ATP2B4 → PMCA4; RAB3GAP1 → RAB3GAP

**`WGS__神經科__Dementia.txt`** （7）：
COA8 → APOPT1; DARS1 → DARS; EPRS1 → EPRS; HYCC1 → FAM126A; FDX2 → FDX1L; GBA1 → GBA; RARS1 → RARS

**`WGS__神經科__Epilepsy.txt`** （18）：
COA8 → APOPT1; FAM194C → C3ORF20; CPAP → CENPJ; CYP27A1 → CTX; CLXN → EFCAB1; EPRS1 → EPRS; FITM2 → FIT2; GBA1 → GBA; GUCY1A1 → GUCY1A3; HADH → HADHSC; HNRNPH1 → HNRPH1; HNRNPR → HNRPR; CILK1 → ICK; BLTP1 → KIAA1109; ASPNAT → NAT8L; ATXN2 → SCA2; AFG2A → SPATA5; AFG2B → SPATA5L1

**`WGS__神經科__Leukoencephalopathy.txt`** （6）：
COA8 → APOPT1; DARS1 → DARS; EPRS1 → EPRS; HYCC1 → FAM126A; FDX2 → FDX1L; RARS1 → RARS

**`WGS__神經科__Neuropathy_NMJ_Myopathy.txt`** （11）：
POPDC1 → BVES; G6PC1 → G6PC; GAN → GAN1; HARS1 → HARS; KARS1 → KARS; LARGE1 → LARGE; MYL11 → MYLPF; RAB7A → RAB7; NHERF1 → SLC9A3R1; TAFAZZIN → TAZ; YARS1 → YARS

**`WGS__神經科__PD_Dystonia_Chorea_Ataxia_Myoclonus.txt`** （7）：
ADPRS → ADPRHL2; FERRY3 → C12ORF4; GBA1 → GBA; MRE11 → MRE11A; COXFA4 → NDUFA4; AFG2A → SPATA5; SPART → SPG20

**`WGS__神經科__Stroke.txt`** （21）：
COQ8A → ADCK3; MT-CYB → CYTB; GUCY1A1 → GUCY1A3; FERMT3 → KIND3; MMUT → MUT; MT-ND4 → ND4; MT-ND5 → ND5; MT-ND6 → ND6; STN1 → OBFC1; WRN → RECQL2; STING1 → STING; MT-TC → TRNC; MT-TF → TRNF; MT-TH → TRNH; MT-TK → TRNK; MT-TL1 → TRNL1; MT-TQ → TRNQ; MT-TS1 → TRNS1; MT-TS2 → TRNS2; MT-TV → TRNV; MT-TW → TRNW

**`WGS__耳鼻喉科__Hearing_loss.txt`** （1）：
RIPOR2 → FAM65B
