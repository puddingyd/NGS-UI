"""Path / env config.

Layout (production):
    NGS_UI/                    ← NGS_UI_HOME
    ├── NGS-UI/                ← REPO_ROOT (this git checkout)
    ├── biotools/              ← Exomiser + LIRICAL CLIs
    ├── vcf/                   ← per-sample VCFs
    ├── tertiary_output/       ← per-sample TSV + sidecars (NOT in git)
    │   └── _index.json        ← optional sample-list cache (lives next to samples)
    └── data/                  ← server runtime state (users.db, jobs/, …)

Every path is derived from NGS_UI_HOME so the whole tree can be moved
by changing one env var (or the default below). Each piece can still be
overridden individually with its own env var when the layout differs.
"""
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def _env_int(name: str, default: int, *, minimum: int | None = None) -> int:
    raw = os.environ.get(name)
    if raw in (None, ""):
        return default
    try:
        value = int(raw)
    except ValueError:
        return default
    if minimum is not None and value < minimum:
        return minimum
    return value

# Default: parent of the repo only for the production-style
# NGS_UI/NGS-UI checkout. A standalone checkout named NGS-UI (for
# example ~/Desktop/NGS-UI) must fall back to the repo itself; otherwise
# startup tries to create runtime dirs beside the checkout.
if "NGS_UI_HOME" in os.environ:
    NGS_UI_HOME = Path(os.environ["NGS_UI_HOME"])
elif REPO_ROOT.name == "NGS-UI" and REPO_ROOT.parent.name == "NGS_UI":
    NGS_UI_HOME = REPO_ROOT.parent
else:
    NGS_UI_HOME = REPO_ROOT

TERTIARY_OUTPUT_ROOT = Path(os.environ.get(
    "TERTIARY_OUTPUT_ROOT",
    NGS_UI_HOME / "tertiary_output",
))

TERTIARY_NF_WORK_ROOT = Path(os.environ.get(
    "NGS_UI_TERTIARY_NF_WORK_ROOT",
    NGS_UI_HOME / "nf_work",
))

# Generated docx reports — every export saves a copy here so reviewers
# can re-download or audit past versions without re-running the report.
REPORT_OUTPUT_DIR = Path(os.environ.get(
    "NGS_UI_REPORT_DIR",
    NGS_UI_HOME / "report",
))
REPORT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Pipeline (林伯昱's in-house Nextflow) production output root. The
# new pipeline always writes here; the NGS-UI worker reads it
# (one-time copy per sample) and otherwise keeps its own state
# in TERTIARY_OUTPUT_ROOT. Override via env when running on a
# non-DGM machine.
PIPELINE_OUT_ROOT = Path(os.environ.get(
    "NGS_UI_PIPELINE_OUT_ROOT",
    "/home/pipeline/tertiary_output",
))

DATA_ROOT = Path(os.environ.get(
    "NGS_UI_DATA_ROOT",
    NGS_UI_HOME / "data",
))

INDEX_PATH = Path(os.environ.get(
    "NGS_UI_INDEX_PATH",
    NGS_UI_HOME / "tertiary_output" / "_index.json",
))

# Parsed SNV payloads are Python dict-heavy and can be far larger than
# their TSV source. Cache compact review TSV payloads, but avoid retaining
# multi-GB complete annotation fallbacks in a long-lived uvicorn process.
SNV_CACHE_MAX = _env_int("NGS_UI_SNV_CACHE_MAX", 2, minimum=0)
SNV_CACHE_MAX_RAW_MB = _env_int("NGS_UI_SNV_CACHE_MAX_RAW_MB", 100, minimum=0)

VCF_DIR = Path(os.environ.get(
    "NGS_UI_VCF_DIR",
    NGS_UI_HOME / "vcf",
))

PHENOTYPE_DIR = Path(os.environ.get(
    "NGS_UI_PHENOTYPE_DIR",
    NGS_UI_HOME / "patient_phenotype",
))

# Uploaded "未完成報告清單" xlsx files + the derived roster.json that
# maps LIS_ID → {mrn, name, test_name, department}. The 載入新個案
# modal reads the roster to auto-fill MRN / 姓名 / Test type.
PATIENT_LIST_DIR = Path(os.environ.get(
    "NGS_UI_PATIENT_LIST_DIR",
    NGS_UI_HOME / "patient_list",
))

# Reference data for phenotype scoring + HPO search: hp.obo and
# phenotype_to_genes.txt stay under NGS_UI_HOME because they are large
# deploy-time data. Fixed WES-I / WES-II / WGS panel files are small,
# curated repo data so git pull can update them together with code.
PHENO_DATA_DIR = Path(os.environ.get(
    "NGS_UI_PHENO_DATA_DIR",
    NGS_UI_HOME / "phenotype_data",
))
GENE_PANELS_DIR = Path(os.environ.get(
    "NGS_UI_GENE_PANELS_DIR",
    REPO_ROOT / "phenotype_data" / "gene_panels",
))
FIXED_PANELS_DIR = Path(os.environ.get(
    "NGS_UI_FIXED_PANELS_DIR",
    REPO_ROOT / "phenotype_data" / "fixed_panels",
))
CUSTOM_GENE_PANELS_DIR = Path(os.environ.get(
    "NGS_UI_CUSTOM_GENE_PANELS_DIR",
    REPO_ROOT / "phenotype_data" / "custom_panels",
))

# Delivered panel/HGNC/dead-zone package. It is intentionally small enough
# to live in the repo, but can be swapped at deploy time after panel/HGNC
# refreshes without touching application code.
NGS_PANEL_DEADZONE_DIR = Path(os.environ.get(
    "NGS_UI_PANEL_DEADZONE_DIR",
    REPO_ROOT / "ngs_panel_deadzone",
))

BIOTOOLS_DIR = Path(os.environ.get(
    "NGS_UI_BIOTOOLS_DIR",
    NGS_UI_HOME / "biotools",
))

# GIAB / GA4GH genome-stratification BEDs (homopolymers, tandem repeats,
# segdups, low mappability, GC extremes, other difficult regions). Used by
# scripts/annotate_giab_strata.py to flag variants in difficult regions and
# by the variant card to badge them. BEDs are large and not committed; place
# them (plus strata_manifest.json) under this dir via download_giab_strata.sh.
GIAB_STRAT_DIR = Path(os.environ.get(
    "NGS_UI_GIAB_STRAT_DIR",
    BIOTOOLS_DIR / "giab_stratification",
))

# Local GeneBe ACMG database (bgzip TSV + lazy SQLite cache) used by
# scripts/annotate_acmg_genebe.py as an OFFLINE replacement for the live
# GeneBe API. Large, not committed; placed under biotools.
GENEBE_DB = Path(os.environ.get(
    "NGS_UI_GENEBE_DB",
    BIOTOOLS_DIR / "genebe" / "genebe_hg38.tsv.gz",
))

# In-house allele frequency sites VCF (built from our own WGS cohort by
# scripts/inhouse_af/publish_af.py). scripts/annotate_inhouse_af.py joins its
# INHOUSE_AC/AN/AF into snv_indel.annotated.tsv and the SNV card shows an
# AF_nckuh row. Large + patient-derived, not committed; deploy under biotools
# via scripts/inhouse_af/deploy_inhouse_af_db.sh. Missing file = silently off.
INHOUSE_AF_DB = Path(os.environ.get(
    "NGS_UI_INHOUSE_AF_DB",
    BIOTOOLS_DIR / "inhouse_af" / "inhouse_af.hg38.vcf.gz",
))

# OMIM annotation table (xlsx). Loaded once at first use and lazily
# reloaded when the file mtime changes. Empty value or missing file
# disables OMIM annotation; variants render with empty Disease lists.
OMIM_XLSX = Path(os.environ.get(
    "NGS_UI_OMIM_XLSX",
    NGS_UI_HOME / "OMIM" / "OMIM.xlsx",
))

JOBS_DIR = DATA_ROOT / "jobs"
JOBS_DIR.mkdir(parents=True, exist_ok=True)

FRONTEND_DIR = Path(os.environ.get(
    "FRONTEND_DIR",
    REPO_ROOT / "frontend",
))

# ---- Bioinformatics tool paths (override via env on the server) ----

EXOMISER_HOME = Path(os.environ.get(
    "EXOMISER_HOME",
    BIOTOOLS_DIR / "exomiser-cli-14.1.0",
))
EXOMISER_JAR  = Path(os.environ.get(
    "EXOMISER_JAR",
    EXOMISER_HOME / "exomiser-cli-14.1.0.jar",
))
EXOMISER_PROPS = Path(os.environ.get(
    "EXOMISER_PROPS",
    EXOMISER_HOME / "application.properties",
))
EXOMISER_DATA_HG38 = Path(os.environ.get(
    "EXOMISER_DATA_HG38",
    EXOMISER_HOME / "data" / "2508_hg38",
))
EXOMISER_DATA_HG19 = Path(os.environ.get(
    "EXOMISER_DATA_HG19",
    EXOMISER_HOME / "data" / "2508_hg19",
))

LIRICAL_HOME = Path(os.environ.get(
    "LIRICAL_HOME",
    BIOTOOLS_DIR / "lirical-cli-2.2.1",
))
LIRICAL_JAR  = Path(os.environ.get(
    "LIRICAL_JAR",
    LIRICAL_HOME / "lirical-cli-2.2.1.jar",
))

JAVA_BIN  = os.environ.get("JAVA_BIN",  "java")
JAVA_OPTS = os.environ.get("JAVA_OPTS", "-Xms4g -Xmx16g")

REDIS_URL = os.environ.get("REDIS_URL", "redis://127.0.0.1:6379/0")

# NCKU intranet EMR (HIS phenotype list + APIM consultation gateway).
# Empty value disables every EMR-related code path so the UI can run
# off-network without raising. The legacy VIP_API.sh hard-coded this
# id; we now read it from systemd's Environment= so it lives outside
# the repo. EMR endpoints are intranet-only so this is currently
# secondary defense, not auth.
EMR_CLIENT_ID = os.environ.get("NGS_UI_EMR_CLIENT_ID", "")

# DRAGEN VCF roots scanned by /api/dragen/vcfs. Each entry is a path
# (with shell globs) — files matching `*hard-filtered.vcf.gz` under any
# of them show up as candidates in the 三級分析 modal. Override via env
# (`:`-separated list of paths/globs) when DRAGEN deposits land
# elsewhere.
DRAGEN_VCF_ROOTS = [
    Path(p) for p in os.environ.get(
        "NGS_UI_DRAGEN_VCF_ROOTS",
        "/home/datalake_Raw/Novaseq:/home/datalake_Intermediate/n102968",
    ).split(":") if p
]
TERTIARY_JOBS_DIR = DATA_ROOT / "jobs" / "tertiary"
TERTIARY_JOBS_DIR.mkdir(parents=True, exist_ok=True)
# Existing deployments may still have completed jobs here. New runs
# write to TERTIARY_JOBS_DIR; the service keeps this path read-only.
LEGACY_DRAGEN_JOBS_DIR = DATA_ROOT / "jobs" / "dragen"

# In-house Nextflow ensemble pipeline outputs scanned by
# /api/dragen/vcfs under mode=inhouse. Per-sample layout:
#   <root>/.../<SID>/04_snv_indel/<SID>.ensemble.fixed.vcf.gz
#   <root>/.../<SID>/05_cnv_sv/<SID>.gcnv.vcf.gz
#   <root>/.../<SID>/05_cnv_sv/<SID>.delly.vcf.gz
#   <root>/.../<SID>/07_mitochondria/<SID>.mito.vcf.gz
# Override via env (`:`-separated). The 三級分析 modal anchors on the
# SNV/Indel VCF; the three siblings are discovered relative to it.
INHOUSE_VCF_ROOTS = [
    Path(p) for p in os.environ.get(
        "NGS_UI_INHOUSE_VCF_ROOTS",
        "/home/datalake_Intermediate/pipeline/nextflow_output",
    ).split(":") if p
]

# Cached scan results for both DRAGEN + in-house VCF discovery. find(1)
# across the datalake can take 1–30 s; the modal reads this file
# directly so it opens instantly. A 🔄 button POSTs to
# /api/dragen/index/refresh to rescan on demand. The file is treated as
# stale after PIPELINE_VCF_INDEX_TTL_HOURS so the modal kicks off a
# background refresh on the next open.
PIPELINE_VCF_INDEX_PATH = DATA_ROOT / "pipeline_vcf_index.json"
PIPELINE_VCF_INDEX_TTL_HOURS = 24

# Secondary-analysis FASTQ discovery and DGX-2 samplesheet creation.
# The UI server scans the server-mounted raw datalake paths and writes
# samplesheets into a DGX-readable staging tree. The generated launch
# command creates the final output dir as the DGX runner and copies the
# samplesheet there, avoiding owner/permission conflicts.
SECONDARY_WES_FASTQ_ROOTS = [
    Path(p) for p in os.environ.get(
        "NGS_UI_SECONDARY_WES_FASTQ_ROOTS",
        "/home/datalake_Raw/NextSeq2000:/home/datalake_Raw/Other/Reanalysis:/datalake_Raw/Other/Reanalysis",
    ).split(":") if p
]
SECONDARY_WGS_FASTQ_ROOTS = [
    Path(p) for p in os.environ.get(
        "NGS_UI_SECONDARY_WGS_FASTQ_ROOTS",
        "/home/datalake_Raw/Novaseq",
    ).split(":") if p
]
SECONDARY_OUTPUT_ROOT = Path(os.environ.get(
    "NGS_UI_SECONDARY_OUTPUT_ROOT",
    "/home/datalake_Intermediate/pipeline/nextflow_output",
))
SECONDARY_SAMPLESHEET_STAGING_ROOT = Path(os.environ.get(
    "NGS_UI_SECONDARY_SAMPLESHEET_STAGING_ROOT",
    "/home/datalake_Intermediate/pipeline/nextflow_samplesheet_staging",
))
SECONDARY_DGX_OUTPUT_ROOT = Path(os.environ.get(
    "NGS_UI_SECONDARY_DGX_OUTPUT_ROOT",
    "/datalake_Intermediate/pipeline/nextflow_output",
))
SECONDARY_DGX_SAMPLESHEET_STAGING_ROOT = Path(os.environ.get(
    "NGS_UI_SECONDARY_DGX_SAMPLESHEET_STAGING_ROOT",
    "/datalake_Intermediate/pipeline/nextflow_samplesheet_staging",
))
SECONDARY_DGX_LAUNCH_ROOT = Path(os.environ.get(
    "NGS_UI_SECONDARY_DGX_LAUNCH_ROOT",
    "/datalake_Intermediate/pipeline/nextflow_launch",
))
SECONDARY_DGX_WORK_ROOT = Path(os.environ.get(
    "NGS_UI_SECONDARY_DGX_WORK_ROOT",
    "/raid/DGM/work",
))
SECONDARY_DGX_ENV_SCRIPT = os.environ.get(
    "NGS_UI_SECONDARY_DGX_ENV_SCRIPT",
    "/datalake_Intermediate/pipeline/pipeline_code/NGS2ndAnalysis_env.sh",
)
SECONDARY_FASTQ_INDEX_PATH = DATA_ROOT / "secondary_fastq_index.json"
SECONDARY_FASTQ_INDEX_TTL_HOURS = int(os.environ.get(
    "NGS_UI_SECONDARY_FASTQ_INDEX_TTL_HOURS",
    "24",
))
