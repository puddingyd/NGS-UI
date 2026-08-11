<!--
Edit the text after each "|" to update the homepage flowchart.
Keep the ## section names and the keys before "|" unchanged because the UI
uses them to place content in the fixed horizontal layout.
-->

# NGS 分析流程

從定序資料、二級與三級分析，到臨床判讀及報告輸出。

## Input & QC

- input | WES / WGS FASTQ or DRAGEN VCF / BAM
- qc | Alignment & QC

## Phenotype & clinical context

- items | HPO · Gene panels × Clinical presentation × Exomiser / LIRICAL / pheno_score

## SNV / Indel

- secondary | Variant calling
- secondary_tools | DeepVariant + HaplotypeCaller or DRAGEN
- tertiary | VEP annotation
- tertiary_tools | ClinVar 20260720 + weekly UI comparison · ClinGen ERepo · LitVar2 · gnomAD · dbNSFP · LOFTEE · Pangolin
- review_filter | AF < 0.01 · VAF > 0.2 · ClinVar rescue
- 1A | ClinVar P/LP ≥ 1★
- 1B | LOFTEE HC
- 1C | ACMG score · AlphaMissense · P-KNN · Pangolin · BayesDel
- Other

## CNV / SV

- secondary | CNV / SV calling
- secondary_tools | gCNV (WES) · CNVkit (WGS) · Delly or DRAGEN
- tertiary | AnnotSV
- tertiary_tools | Gene content · ACMG class · region evidence
- prioritization | Clinical
- prioritization | Pathogenic
- context | Panel-gene overlap · pheno_score

## mtDNA

- secondary | mtDNA calling
- secondary_tools | NCKUH / DRAGEN mitochondrial calls
- tertiary | mtDNA annotation
- tertiary_tools | VEP · ClinVar · gnomAD-mito · NCKUH AF
- prioritization | Pathogenic
- prioritization | Rare / reported
- prioritization | Other
- context | Disease / gene context

## STR

- secondary | Repeat calling
- secondary_tools | GangSTR
- tertiary | STRchive
- tertiary_tools | Repeat count · thresholds · inheritance
- prioritization | Pathogenic
- prioritization | Intermediate
- prioritization | Normal
- context | Disease / gene context

## PGx

- secondary | PGx typing
- secondary_tools | VCF / BAM-based genotyping
- tertiary | PGx annotation
- tertiary_tools | PharmCAT · StellarPGx · OptiType
- prioritization | CPIC / FDA recommendations
- context | Medication × Gene × CPIC / FDA recommendations

## Review

- summary | Evidence review
- detail | IGV genome viewer · Genotype-phenotype correlation · Structured manual ACMG · Observed cases
- item | Diagnostic findings
- item | ACMG SF · Stroke · Carrier screening
- item | PGx recommendations

## Reports

- 診斷報告 | Causative · Other
- 健檢報告 | Secondary findings · PGx
- PDF 摘要 | Selected variants
