<!--
Edit the text after each "|" to update the homepage flowchart.
Keep the ## section names and the keys before "|" unchanged because the UI
uses them to place content in the fixed horizontal layout.
-->

# NGS 分析流程

從定序資料、二級與三級分析，到臨床判讀及報告輸出。

## Input & QC

- input | WES / WGS FASTQ or DRAGEN VCF / BAM
- qc | QC & Alignment
- tools | fastp · Parabricks / BWA · Mosdepth
- split | 分流為五條分析線

## Phenotype & clinical context

- items | HPO · Gene panels × Clinical presentation × Exomiser / LIRICAL / pheno_score

## SNV / Indel

- secondary | Variant calling
- secondary_tools | DeepVariant + HaplotypeCaller or DRAGEN
- tertiary | VEP annotation
- tertiary_tools | ClinVar · gnomAD · dbNSFP · LOFTEE · Pangolin
- 1A | ClinVar P/LP ≥ 1★
- 1B | LOFTEE HC
- 1C | ACMG score · P-KNN · AlphaMissense · Pangolin · BayesDel
- Other

## CNV / SV

- secondary | CNV / SV calling
- secondary_tools | gCNV (WES) · CNVkit (WGS) · Delly or DRAGEN
- tertiary | AnnotSV annotation
- tertiary_tools | Gene content · ACMG class · region evidence
- prioritization | Clinical · Pathogenic
- context | Panel-gene overlap · pheno_score

## mtDNA

- secondary | mtDNA calling
- secondary_tools | NCKUH / DRAGEN mitochondrial calls
- tertiary | Mito annotation
- tertiary_tools | VEP · ClinVar · gnomAD-mito
- prioritization | Pathogenic · Rare / reported · Other
- context | Disease / gene context

## STR

- secondary | Repeat calling
- secondary_tools | GangSTR · ExpansionHunter
- tertiary | STRchive
- tertiary_tools | Repeat count · thresholds · inheritance
- prioritization | Pathogenic · Intermediate · Normal
- context | Disease / gene context

## PGx

- secondary | PGx typing
- secondary_tools | VCF / BAM-based genotyping
- tertiary | PGx annotation
- tertiary_tools | PharmCAT · StellarPGx · OptiType
- prioritization | CPIC / FDA recommendations
- context | Medication / clinical context

## Review

- summary | Evidence review
- detail | IGV · Disease / phenotype concordance · Reviewer selection
- item | Diagnostic findings
- item | ACMG SF · Stroke · Carrier
- item | PGx recommendations

## Reports

- 診斷報告 | Causative · Other
- 健檢報告 | Secondary findings · PGx
- PDF 摘要 | Selected variants
