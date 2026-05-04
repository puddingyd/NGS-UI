# Load packages
pkgs <- c("dplyr", "readxl", "stringr", "tidyr", "DT", "htmlwidgets", "parallel", "purrr", "data.table", "readr", "reactable", "htmltools", "beepr", "emayili", "magrittr", "jsonlite")
invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(
  library(p, character.only = TRUE)
)))

# Make CLI tools reachable when R is launched from RStudio/Finder (GUI sessions
# inherit the minimal launchd PATH, so tools in /usr/local/bin are otherwise
# invisible to system() calls such as bgzip or git).
Sys.setenv(PATH = paste("/usr/local/bin", Sys.getenv("PATH"), sep = ":"))

NGS_pipeline <- "/Users/changym/Desktop/NGS_pipeline/"
REPO         <- normalizePath("/Users/changym/Desktop/NGS_pipeline/vcf-analysis-hg38-R", mustWork = FALSE)  # Same repo as hg38; webdata is shared
GENOME_BUILD <- "hg19"

### Load ClinVar — auto-refresh from NCBI FTP when local copy is older than
### 14 days. Successful refresh stores the txt + a .rds cache and deletes
### older clinvar_*.{txt,rds} so only the active version stays on disk.
{
  clinvar_dir     <- "/Volumes/genetics/humandb/hg19/"
  clinvar_url     <- "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"
  clinvar_suffix  <- ".txt"
  clinvar_max_age <- 14L

  setwd(clinvar_dir)
  suffix_re <- gsub("\\.", "\\\\.", clinvar_suffix)
  txt_re    <- paste0("^clinvar_\\d{8}", suffix_re, "$")
  existing  <- list.files(clinvar_dir, pattern = txt_re)
  local_date <- if (length(existing)) {
    max(regmatches(existing, regexpr("\\d{8}", existing)))
  } else NA_character_
  age_days <- if (!is.na(local_date)) {
    as.integer(Sys.Date() - as.Date(local_date, "%Y%m%d"))
  } else NA_integer_

  needs_update <- is.na(local_date) || age_days > clinvar_max_age
  result <- list(ok = FALSE)
  if (needs_update) {
    message("ClinVar local copy is ",
            if (is.na(local_date)) "missing" else paste0(local_date, " (", age_days, "d old)"),
            "; fetching latest from NCBI FTP...")
    result <- tryCatch({
      download.file(clinvar_url, "clinvar.vcf.gz", quiet = FALSE, method = "libcurl")
      system("bgzip -df clinvar.vcf.gz")
      hdr <- readLines("clinvar.vcf")
      new_date <- gsub("-", "", sub("^##fileDate=", "",
                                    hdr[grepl("^##fileDate=", hdr)][1]))
      txt_path <- paste0("clinvar_", new_date, clinvar_suffix)
      writeLines(hdr[grep("^#CHROM", hdr):length(hdr)], txt_path)

      cv <- data.table::fread(txt_path, sep = "\t", quote = "", header = FALSE,
                              col.names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"))
      cv[, CLNHGVSn    := stringr::str_extract(INFO, "(?<=CLNHGVS=)[^;]+")]
      cv[, CLNSIGn     := stringr::str_extract(INFO, "(?<=CLNSIG=)[^;]+")]
      cv[, CLNSIGCONFn := stringr::str_extract(INFO, "(?<=CLNSIGCONF=)[^;]+")]
      cv[, CLNREVSTATn := stringr::str_extract(INFO, "(?<=CLNREVSTAT=)[^;]+")]
      cv[, ID := paste0("chr", CHROM, "-", POS, "-", REF, "-", ALT)]
      data.table::fwrite(cv, txt_path, sep = "\t", quote = FALSE)
      saveRDS(as.data.frame(cv), sub("\\.txt$", ".rds", txt_path))
      file.remove("clinvar.vcf")
      list(ok = TRUE, Clinvar = as.data.frame(cv), clinvar_date = new_date)
    }, error = function(e) {
      warning("ClinVar update failed (", conditionMessage(e),
              "); falling back to local copy.")
      list(ok = FALSE)
    })
  }

  if (result$ok) {
    Clinvar      <- result$Clinvar
    clinvar_date <- result$clinvar_date
  } else {
    if (is.na(local_date)) stop("ClinVar: no local copy and update failed.")
    clinvar_date <- local_date
    txt_path <- paste0("clinvar_", clinvar_date, clinvar_suffix)
    rds_path <- sub("\\.txt$", ".rds", txt_path)
    if (file.exists(rds_path)) {
      Clinvar <- readRDS(rds_path)
    } else {
      Clinvar <- read.table(txt_path, header = TRUE, sep = "\t", quote = "", fill = TRUE)
      saveRDS(Clinvar, rds_path)
    }
  }

  # Auto-delete stale clinvar_<otherdate>.{txt,rds}
  rds_suffix_re <- sub("\\\\\\.txt", "\\\\.rds", suffix_re)
  all_files <- list.files(clinvar_dir,
                          pattern = paste0("^clinvar_\\d{8}", suffix_re, "$|^clinvar_\\d{8}", rds_suffix_re, "$"))
  to_remove <- all_files[!grepl(paste0("^clinvar_", clinvar_date, "[._]"), all_files)]
  if (length(to_remove)) {
    message("Removing stale ClinVar files: ", paste(to_remove, collapse = ", "))
    file.remove(file.path(clinvar_dir, to_remove))
  }
}

### Create directory for each file
{
  setwd("/Users/changym/Desktop/VCF")
  files <- list.files(pattern = "_ann\\.txt$")

  for (file in files) {
    dir_name <- sub("_ann\\.txt$", "", file)
    if (!dir.exists(dir_name)) {
      dir.create(dir_name)
    }
    file.rename(file, file.path(dir_name, file))
    #file.copy("patient_phenotype.xlsx", file.path(dir_name, "patient_phenotype.xlsx"))
    file.copy("genes_interested.xlsx", file.path(dir_name, "genes_interested.xlsx"))
  }
}

### Start here — resolve which IDs to process, plus any per-sample gene overrides.
### CLI args form an interleaved sequence of (ID [, genes])* — any arg that
### contains a comma is treated as the gene list for the preceding ID. Example:
###   Rscript vcf-analysis-hg19.R 26WE0060 'SLC25A13, ASS1' 26WE0061
### gives IDs = c("26WE0060", "26WE0061") and writes SLC25A13, ASS1 into
### 26WE0060/genes_interested.xlsx (overwriting any prior content). Samples
### without an attached gene arg keep whatever genes_interested.xlsx is
### already in their dir. Falls back to analysis_quere.txt when run without
### CLI args so RStudio workflows still work.
cli_args <- commandArgs(trailingOnly = TRUE)
cli_genes_map <- list()
if (length(cli_args)) {
  cli_ids <- character()
  cur_id  <- NULL
  for (arg in cli_args) {
    if (grepl(",", arg, fixed = TRUE)) {
      if (is.null(cur_id)) {
        warning("Gene list '", arg, "' appeared before any sample ID; ignored.")
      } else {
        cli_genes_map[[cur_id]] <- trimws(arg)
      }
    } else {
      cur_id <- trimws(arg)
      if (nzchar(cur_id)) cli_ids <- c(cli_ids, cur_id)
    }
  }
  IDs <- unique(cli_ids)
  if (!length(IDs)) stop("No valid sample IDs in CLI args.")
  message("Using CLI-provided IDs: ", paste(IDs, collapse = ", "))
  if (length(cli_genes_map)) {
    message("CLI-provided genes: ",
            paste(sprintf("%s → %s", names(cli_genes_map),
                          unlist(cli_genes_map)), collapse = " | "))
  }
} else {
  id_file <- "/Users/changym/Desktop/VCF/analysis_quere.txt"
  if (!file.exists(id_file)) stop("ID list not found: ", id_file)
  IDs <- unique(trimws(readLines(id_file, warn = FALSE)))
  IDs <- IDs[IDs != ""]
  if (!length(IDs)) stop("ID list is empty: ", id_file)
}

### Refresh webdata/index.json from NGS_list.xlsx so the web tool (and the
### phenotype seeder below) can look up samples by LIS_ID / Name / MRN.
### Also syncs the `yield` column back. Identical to the hg38 script —
### both builds share the same REPO and index.json, so running either one
### keeps the index up to date.
{
  ngs_xlsx <- file.path(REPO, "NGS_list.xlsx")
  if (file.exists(ngs_xlsx)) {
    wes <- suppressMessages(readxl::read_excel(ngs_xlsx, sheet = "WES list"))

    # Collect yield per LIS_ID from all saved reports
    reports_dir <- file.path(REPO, "webdata", "reports")
    yield_map <- list()
    if (dir.exists(reports_dir)) {
      for (rf in list.files(reports_dir, pattern = "\\.json$", full.names = TRUE)) {
        lid <- sub("\\.json$", "", basename(rf))
        rep <- tryCatch(jsonlite::read_json(rf), error = function(e) NULL)
        if (!is.null(rep)) {
          if (!is.null(rep$yield)) {
            yield_map[[lid]] <- as.integer(rep$yield)
          } else if (!is.null(rep$status)) {
            yield_map[[lid]] <- as.integer(any(unlist(rep$status) == "1"))
          }
        }
      }
    }

    if (!"yield" %in% names(wes)) wes$yield <- NA_integer_
    if (length(yield_map)) {
      m <- match(as.character(wes$LIS_ID), names(yield_map))
      has <- !is.na(m)
      wes$yield[has] <- unlist(yield_map)[m[has]]
    }

    idx <- wes %>%
      transmute(
        LIS_ID   = as.character(LIS_ID),
        Name     = as.character(Name),
        MRN      = as.character(ID),                # xlsx "ID" column is the medical record number
        Test     = as.character(Test),
        Category = as.character(Category),
        Yield    = suppressWarnings(as.integer(yield)),
        Tag      = if ("Tag" %in% names(.)) as.character(Tag) else NA_character_
      ) %>%
      dplyr::filter(!is.na(LIS_ID) & LIS_ID != "")
    cat_options <- sort(unique(na.omit(idx$Category)))
    cat_options <- cat_options[cat_options != ""]
    dir.create(file.path(REPO, "webdata"), showWarnings = FALSE, recursive = TRUE)
    jsonlite::write_json(idx,
                         file.path(REPO, "webdata", "index.json"),
                         auto_unbox = TRUE, na = "null", pretty = TRUE)
    jsonlite::write_json(list(category_options = cat_options),
                         file.path(REPO, "webdata", "options.json"),
                         auto_unbox = TRUE, na = "null", pretty = TRUE)
  } else {
    warning("NGS_list.xlsx not found at ", ngs_xlsx, " — webdata/index.json not updated")
  }
}

# ---------------------------------------------------------------------
# Remote VCF fetch + ANNOVAR fallback
#
# When a sample's <ID>_ann.txt isn't already present locally, the wrapper
# below SSHes to the institutional server, locates the GPU-called VCF
# (path pattern */NextSeq2000/<run>/<ID>_<idx>/nv_result/<ID>_<idx>_gpu.vcf),
# scp's it down, and runs table_annovar.pl with the same protocol set the
# downstream R code reads. A flat path index in ~/.vcf_pipeline_cache/
# avoids re-walking the remote tree on every miss; refreshed on cache miss
# or when the file is older than INDEX_MAX_AGE_DAYS.
#
# Requires SSH key auth (BatchMode=yes; will fail loudly if not set up).
# Cache: ~/.vcf_pipeline_cache/server_vcf_index.tsv. Built on first miss
# (no file yet) and refreshed on a subsequent miss; deliberately never
# stale-by-time so a known sample stays a fast lookup forever — new
# samples that aren't in the cache trigger one auto-rebuild and are
# served from the refreshed index.
REMOTE_SERVER      <- "n102968@192.168.84.91"
REMOTE_ROOT        <- "/home/datalake_Intermediate/NextSeq2000"
HUMANDB_HG19       <- "/Volumes/genetics/humandb/hg19"
ANNOVAR_BIN        <- "~/bin/annovar"
INDEX_CACHE_FP     <- path.expand("~/.vcf_pipeline_cache/server_vcf_index.tsv")

refresh_remote_index <- function() {
  dir.create(dirname(INDEX_CACHE_FP), showWarnings = FALSE, recursive = TRUE)
  message(sprintf("Refreshing remote VCF index from %s …", REMOTE_SERVER))
  cmd <- sprintf(
    "ssh -o BatchMode=yes %s 'find %s -path \"*/nv_result/*_gpu.vcf\" 2>/dev/null'",
    REMOTE_SERVER, REMOTE_ROOT
  )
  paths <- suppressWarnings(system(cmd, intern = TRUE))
  if (!length(paths)) {
    warning("Remote find returned no paths; check SSH key + server availability")
    return(invisible(FALSE))
  }
  writeLines(paths, INDEX_CACHE_FP)
  message(sprintf("Indexed %d VCFs", length(paths)))
  invisible(TRUE)
}
find_remote_vcf <- function(LIS_ID) {
  if (!file.exists(INDEX_CACHE_FP)) refresh_remote_index()
  if (!file.exists(INDEX_CACHE_FP)) return(NA_character_)
  paths <- readLines(INDEX_CACHE_FP, warn = FALSE)
  pat <- sprintf("/%s_[^/]+/nv_result/[^/]+_gpu\\.vcf$", LIS_ID)
  hit <- grep(pat, paths, value = TRUE)
  if (length(hit)) return(hit[1])
  # Cache miss — refresh once and retry before giving up.
  message(sprintf("%s not in cached index; refreshing once …", LIS_ID))
  refresh_remote_index()
  paths <- readLines(INDEX_CACHE_FP, warn = FALSE)
  hit <- grep(pat, paths, value = TRUE)
  if (length(hit)) hit[1] else NA_character_
}
ensure_ann_txt <- function(LIS_ID) {
  ann_fp <- paste0(LIS_ID, "_ann.txt")
  if (file.exists(ann_fp)) return(invisible())

  remote_path <- find_remote_vcf(LIS_ID)
  if (is.na(remote_path) || !nzchar(remote_path)) {
    stop(sprintf("Cannot find VCF for %s under %s on %s",
                 LIS_ID, REMOTE_ROOT, REMOTE_SERVER))
  }

  dir.create("analysis_files", showWarnings = FALSE, recursive = TRUE)
  local_vcf <- file.path("analysis_files", paste0(LIS_ID, "_remote.vcf"))
  message(sprintf("Downloading %s → %s", remote_path, local_vcf))
  rc <- system(sprintf("scp -q -o BatchMode=yes %s:%s %s",
                       REMOTE_SERVER, shQuote(remote_path), shQuote(local_vcf)))
  if (rc != 0 || !file.exists(local_vcf)) stop("scp failed for ", remote_path)

  message("Running ANNOVAR table_annovar.pl …")
  out_prefix <- file.path("analysis_files", LIS_ID)
  cmd <- sprintf(
    paste("perl %s/table_annovar.pl %s %s -buildver hg19",
          "-out %s -remove",
          "-protocol refGeneWithVer,gnomad211_exome,twnaf_annovarin,clinvar_20250623,dbscsnv11",
          "-operation g,f,f,f,f",
          "-nastring . --vcfinput -polish"),
    ANNOVAR_BIN, shQuote(local_vcf), shQuote(HUMANDB_HG19), shQuote(out_prefix)
  )
  rc <- system(cmd)
  if (rc != 0) stop("ANNOVAR failed for ", LIS_ID)

  # ANNOVAR writes <prefix>.hg19_multianno.txt; rename to the
  # <ID>_ann.txt the downstream read.table expects, in the cwd. Drop
  # the matching multianno.vcf — nothing reads it. The downloaded
  # _remote.vcf stays put for any future re-annotation.
  multianno_fp <- paste0(out_prefix, ".hg19_multianno.txt")
  if (!file.exists(multianno_fp)) stop("ANNOVAR output missing: ", multianno_fp)
  file.rename(multianno_fp, ann_fp)
  unlink(paste0(out_prefix, ".hg19_multianno.vcf"))
  message("Wrote ", ann_fp)
  invisible()
}

for (ID in IDs) {
message("\n=== Processing sample: ", ID, " ===")
VCF <- paste0("/Users/changym/Desktop/VCF/", ID, "/")

# CLI gene override: if run-vcf / Rscript was called with a gene list for
# this sample, write it to <VCF>/genes_interested.xlsx (cell A1, overwrite).
# Happens after the startup block has already copied the template xlsx into
# each sample dir, so the file path is guaranteed to exist for new samples
# too. The webdata block later reads cell A1 back, so we don't need to
# touch any other code path.
if (!is.null(cli_genes_map[[ID]])) {
  dir.create(VCF, showWarnings = FALSE, recursive = TRUE)
  genes_xlsx <- file.path(VCF, "genes_interested.xlsx")
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("openxlsx required for CLI gene override; install.packages('openxlsx')")
  }
  openxlsx::write.xlsx(
    data.frame(V1 = cli_genes_map[[ID]], stringsAsFactors = FALSE),
    genes_xlsx,
    colNames = FALSE
  )
  message("Wrote CLI genes for ", ID, " → ", genes_xlsx,
          " (", cli_genes_map[[ID]], ")")
}

### Genotype prioritization
{ # Set up
  WES_pipeline <- "/Users/changym/Desktop/NGS_pipeline/"
  setwd(VCF)
  # Remove regenerated outputs from previous runs of THIS sample so a fresh
  # run doesn't sit alongside stale copies. ID-prefixed so other samples
  # sharing this project root aren't touched. Includes legacy filenames the
  # current script no longer produces (so existing checkouts get cleaned up
  # on their next run). The previous list.files() form had no path = arg
  # and only scanned the cwd, so it never actually reached analysis_files/
  # and was effectively a no-op.
  unlink(file.path("analysis_files", paste0(ID, c(
    "_annotated.txt",         # legacy, no longer written
    "_anno_combined.txt",     # legacy uncompressed form
    "_anno_combined.txt.gz",  # current
    "_rank.txt"
  ))))
  
  dir.create("analysis_files", showWarnings = FALSE)
  # One-time migration: gzip pre-existing uncompressed ANNOVAR / VEP caches
  # in place so a sample re-run after the format change still hits the cache
  # instead of re-running the slow steps. Safe no-op once the .gz file exists.
  for (suffix in c(".hg19_multianno.txt", ".vep.vcf")) {
    old <- paste0("analysis_files/", ID, suffix)
    new <- paste0(old, ".gz")
    if (file.exists(old) && !file.exists(new)) {
      system(sprintf("gzip %s", shQuote(old)))
    }
  }
  
  if (file.exists(paste0(ID, "_ann.txt.gz")) && !file.exists(paste0(ID, "_ann.txt"))) {
    system(paste("bgzip -d", paste0(ID, "_ann.txt.gz")))
  }

  # If neither <ID>_ann.txt nor _ann.txt.gz is here, fetch the upstream
  # _gpu.vcf from the institutional server and run ANNOVAR locally to
  # rebuild the same _ann.txt the downstream code expects. No-op when
  # the file already exists.
  ensure_ann_txt(ID)

  # Make a vcf file
  vcf_ann <- read.table(paste0(ID, "_ann.txt"), header = TRUE, sep = "\t", quote = "", fill = TRUE)
  message(sprintf("[%s] read %d rows from _ann.txt", ID, nrow(vcf_ann)))
  if (nrow(vcf_ann)) {
    message(sprintf("[%s] sample row Otherinfo10 (FILTER) = %s", ID, vcf_ann$Otherinfo10[1]))
    message(sprintf("[%s] sample row Otherinfo13 (SAMPLE) = %s", ID, vcf_ann$Otherinfo13[1]))
    message(sprintf("[%s] FILTER value distribution: %s", ID,
                    paste(names(table(vcf_ann$Otherinfo10)),
                          unname(table(vcf_ann$Otherinfo10)),
                          sep="=", collapse=", ")))
  }

  # Normalise chromosome notation. Source VCFs vary: ANNOVAR-against-
  # GRCh37 strips "chr" so the multianno output ends up "1"/"2"/…/"X",
  # while older _ann.txt copies kept "chr1"/"chr2"/…/"chrX". Downstream
  # code (chr_all filter at L281, ID = chrN-pos-ref-alt at L324)
  # consistently expects the "chr" prefix, so coerce both Chr and
  # Otherinfo4 (= the renamed #CHROM) here.
  if ("Chr" %in% names(vcf_ann)) {
    vcf_ann$Chr <- paste0("chr", sub("^chr", "", as.character(vcf_ann$Chr)))
  }
  if ("Otherinfo4" %in% names(vcf_ann)) {
    vcf_ann$Otherinfo4 <- paste0("chr", sub("^chr", "", as.character(vcf_ann$Otherinfo4)))
  }

  # Some upstream callers (e.g. NVIDIA Parabricks DeepVariant) leave
  # FILTER as "." instead of "PASS" because they don't run a hard-filter
  # pass. Accept both — we still rely on AF / DP / VAF cutoffs below to
  # drop low-quality calls.
  vcf_ann <- vcf_ann %>% filter(Otherinfo10 == "PASS" | Otherinfo10 == ".")
  message(sprintf("[%s] after FILTER==PASS|.: %d rows", ID, nrow(vcf_ann)))
  
  vcf_ann$AF <- as.numeric(vcf_ann$AF)
  vcf_ann$AF[is.na(vcf_ann$AF)] <- 0
  vcf_ann <- vcf_ann %>% filter(!(AF > 0.01))
  message(sprintf("[%s] after AF<=0.01: %d rows", ID, nrow(vcf_ann)))

  chr_all <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
  vcf_ann <- vcf_ann[vcf_ann$Chr %in% chr_all, ]
  message(sprintf("[%s] after chr filter: %d rows", ID, nrow(vcf_ann)))

  vcf_ann <- vcf_ann %>% mutate(DP = sub('^([^:]*:){2}([^:]*):.*', '\\2', vcf_ann$Otherinfo13))
  vcf_ann$DP <- as.numeric(vcf_ann$DP)
  vcf_ann <- vcf_ann %>% mutate(Alt_DP = sub('^[^,]*,([^:]*):[^:]*:.*', '\\1', vcf_ann$Otherinfo13))
  vcf_ann$Alt_DP <- as.numeric(vcf_ann$Alt_DP)
  vcf_ann <- vcf_ann %>% mutate(Alt_ratio = vcf_ann$Alt_DP / vcf_ann$DP)
  vcf_ann <- vcf_ann %>% filter(Alt_ratio >= 0.2)
  message(sprintf("[%s] after Alt_ratio>=0.2: %d rows (sample DP=%s, Alt_ratio=%s)",
                  ID, nrow(vcf_ann),
                  if (nrow(vcf_ann)) vcf_ann$DP[1] else "—",
                  if (nrow(vcf_ann)) round(vcf_ann$Alt_ratio[1], 3) else "—"))
  vcf_ann <- vcf_ann %>% filter(DP >= 20)
  message(sprintf("[%s] after DP>=20: %d rows", ID, nrow(vcf_ann)))

  # Surface the read-depth / VAF fields the webdata block + web UI expect.
  # FORMAT layout is GT:AD:DP:GQ:PL:VAF, so AD is the 2nd colon-delimited
  # field of Otherinfo13. DP / Alt_ratio are already parsed above.
  vcf_ann <- vcf_ann %>% mutate(
    AD          = sub("^[^:]*:([^:]*):.*", "\\1", Otherinfo13),
    total_depth = DP,
    alt_af      = Alt_ratio
  )
  
  vcf_ann <- select(vcf_ann, Otherinfo4, Otherinfo5, Otherinfo6, Otherinfo7, Otherinfo8, Otherinfo9, Otherinfo10, Otherinfo11, Otherinfo12, Otherinfo13, everything())
  {
    colnames(vcf_ann)[colnames(vcf_ann) == "Otherinfo1"] <- "zygosity"
    colnames(vcf_ann)[colnames(vcf_ann) == "Otherinfo4"] <- "#CHROM"
    colnames(vcf_ann)[colnames(vcf_ann) == "Otherinfo5"] <- "POS"
    colnames(vcf_ann)[colnames(vcf_ann) == "Otherinfo6"] <- "ID"
    colnames(vcf_ann)[colnames(vcf_ann) == "Otherinfo7"] <- "REF"
    colnames(vcf_ann)[colnames(vcf_ann) == "Otherinfo8"] <- "ALT"
    colnames(vcf_ann)[colnames(vcf_ann) == "Otherinfo9"] <- "QUAL"
    colnames(vcf_ann)[colnames(vcf_ann) == "Otherinfo10"] <- "FILTER"
    colnames(vcf_ann)[colnames(vcf_ann) == "Otherinfo11"] <- "INFO"
    colnames(vcf_ann)[colnames(vcf_ann) == "Otherinfo12"] <- "FORMAT"
    colnames(vcf_ann)[colnames(vcf_ann) == "Otherinfo13"] <- "sample"
    vcf_ann <- subset(vcf_ann, ALT != "0")
  }
  vcf_clean <- vcf_ann[, 1:10]
  write.table(vcf_clean, "processing.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  variant <- readLines("processing.txt")
  header <- readLines(paste0(WES_pipeline, "header.txt"))
  writeLines(c(header, variant), paste("analysis_files/", ID, "_hg19.vcf", sep = ""))
  file.remove("processing.txt")
    
  ### Filter
  if (nrow(vcf_ann) == 0) {
    stop(sprintf(
      "vcf_ann is empty after upstream filters for %s — check FILTER (Otherinfo10) values, chromosome notation, and the AF/DP cutoffs in %s_ann.txt",
      ID, ID))
  }
  vcf_ann$ID <- paste0(vcf_ann$'#CHROM', "-", vcf_ann$POS, "-", vcf_ann$REF, "-", vcf_ann$ALT)
  vcf_input <- vcf_ann
  
  ## Filter function and Clinvar
  # Incorporate Clinvar
  vcf_input <- left_join(vcf_input, Clinvar[, c("ID", "CLNHGVSn", "CLNSIGn", "CLNSIGCONFn", "CLNREVSTATn")], by = "ID")
  
  # ADA score > 0.957813, RF score > 0.584
  vcf_input$dbscSNV_ADA_SCORE <- as.numeric(vcf_input$dbscSNV_ADA_SCORE)
  vcf_input$dbscSNV_ADA_SCORE[is.na(vcf_input$dbscSNV_ADA_SCORE)] <- 0
  vcf_input$dbscSNV_RF_SCORE <- as.numeric(vcf_input$dbscSNV_RF_SCORE)
  vcf_input$dbscSNV_RF_SCORE[is.na(vcf_input$dbscSNV_RF_SCORE)] <- 0
  
  # Filter
  vcf_input_func <- vcf_input %>% filter((Func.refGeneWithVer == "exonic" | Func.refGeneWithVer == "splicing" | Func.refGeneWithVer == "exonic;splicing") & !(ExonicFunc.refGeneWithVer == "synonymous SNV"))
  vcf_input_splicing <- vcf_input %>% filter(dbscSNV_ADA_SCORE > 0.957813 | dbscSNV_RF_SCORE > 0.584)
  vcf_input_clinvar <- vcf_input %>% filter(CLNSIGn == "Pathogenic" | CLNSIGn == "Likely_pathogenic" | CLNSIGn == "Pathogenic/Likely_pathogenic" | CLNSIGn == "Uncertain_significance" | CLNSIGn == "Conflicting_classifications_of_pathogenicity")
  #vcf_input <- rbind(vcf_input_func, vcf_input_splicing, vcf_input_clinvar)
  #vcf_input <- vcf_input %>% distinct()
  
  ### Generate a filtered vcf file for further annotation
  # Extract header of the vcf
  text <- readLines(paste0("analysis_files/", ID, "_hg19.vcf"))
  header <- text[1:(grep("#CHROM", text))]
  
  # Combine the variants and header
  vcf_filtered <- select(vcf_input, '#CHROM', POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample)
  write.table(vcf_filtered, "processing.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  writeLines(c(header, readLines("processing.txt")), "processing.vcf")
  file.remove("processing.txt")
  
  # Annovar annotation
  setwd(VCF)
  input_vcf <- "processing.vcf"
  output_vcf <- paste0(ID)
  
  annovarCommand <- sprintf('perl ~/bin/annovar/table_annovar.pl %s /Volumes/genetics/humandb/hg19 -buildver hg19 -out %s -remove -protocol dbnsfp47a -operation f -nastring . --vcfinput -polish', input_vcf, output_vcf)
  
  if (!file.exists(paste0("analysis_files/", ID, ".hg19_multianno.txt.gz"))) {
    system(annovarCommand)
    file.rename(paste0(ID, ".hg19_multianno.txt"), paste0("analysis_files/", ID, ".hg19_multianno.txt"))
    system(sprintf("gzip %s", shQuote(paste0("analysis_files/", ID, ".hg19_multianno.txt"))))
  }

  annovar <- fread(cmd = paste("gunzip -c", shQuote(paste0("analysis_files/", ID, ".hg19_multianno.txt.gz"))), sep   = "\t", quote = "")
  colnames(annovar)[colnames(annovar) == "Otherinfo6"] <- "ID"
  
  
  ### VEP - HGVS
  vep_cmd <- paste(
    "docker run --rm",
    paste0("-v ", shQuote(VCF),  ":/vcf"),
    paste0("-v ", shQuote("/Volumes/genetics/humandb/hg19/vep_cache"), ":/opt/vep/.vep"),
    paste0("-v ", shQuote("/Volumes/genetics/humandb/hg19/dbNSFP"), ":/dbNSFP"),
    paste0("-v ", shQuote("/Volumes/genetics/humandb/hg19/spliceai"), ":/spliceai"),
    paste0("-v ", shQuote("/Volumes/genetics/humandb/hg19/maxentscan"), ":/maxentscan"),
    "ensemblorg/ensembl-vep:latest",
    "vep",
    "-i /vcf/processing.vcf",
    paste0("-o /vcf/analysis_files/", ID, ".vep.vcf"),
    "--offline --cache --dir_cache /opt/vep/.vep",
    "--cache_version 114",
    "--assembly GRCh37",
    "--fasta /opt/vep/.vep/Homo_sapiens.GRCh37.dna.toplevel.fa",
    "--hgvs --symbol --canonical --pick --refseq --mane",
    "--vcf",
    #"--fields", shQuote("Uploaded_variation,Location,Allele,SYMBOL,Gene,Feature,Consequence,HGVSc,HGVSp,CANONICAL,MANE_SELECT,MANE_PLUS_CLINICAL,Existing_variation,IMPACT"),
    "--plugin", shQuote("dbNSFP,/dbNSFP/dbNSFP5.2a_grch37.gz,VARITY_R_score,ESM1b_score,MutPred2_score"), # gnomAD4.1_joint_AF,gnomAD4.1_joint_EAS_AF,SIFT_score,SIFT4G_score,Polyphen2_HDIV_score,Polyphen2_HVAR_score,MutationTaster_score,MutationAssessor_score,PROVEAN_score,VEST4_score,MetaSVM_score,MetaLR_score,M-CAP_score,MVP_score,MPC_score,PrimateAI_score,DEOGEN2_score,BayesDel_addAF_score,BayesDel_noAF_score,ClinPred_score,LIST-S2_score,CADD_phred,DANN_score,fathmm-XF_coding_score,Eigen-raw_coding,Eigen-PC-raw_coding,GERP++_RS,phyloP100way_vertebrate,phastCons100way_vertebrate,MutPred2_score,AlphaMissense_score,MetaRNN_score,REVEL_score,VARITY_R_score,ESM1b_score
    '--plugin "SpliceAI,snv=/spliceai/spliceai_scores.raw.snv.hg19.vcf.gz,indel=/spliceai/spliceai_scores.raw.indel.hg19.vcf.gz,cutoff=0.5"',
    "--plugin MaxEntScan,/maxentscan/fordownload",
    "--no_stats --force_overwrite --verbose"
  )
  
  if (!file.exists(paste0("analysis_files/", ID, ".vep.vcf.gz"))) {
    system(vep_cmd)
    system(sprintf("gzip %s", shQuote(paste0("analysis_files/", ID, ".vep.vcf"))))
  }

  ## Make annotation into datatable. read_lines (readr) and fread both
  ## auto-detect gzip from the .gz extension, so the rest of the block works
  ## unchanged.
  vep <- paste0("analysis_files/", ID, ".vep.vcf.gz")
  hdr <- read_lines(vep, n_max = 500) 
  csq_line <- hdr[grepl("^##INFO=<ID=CSQ", hdr)]
  csq_names <- strsplit(str_match(csq_line, "Format: (.+?)\"")[,2], "\\|", fixed = FALSE)[[1]]
  
  # Read variant rows
  # Pipe through gunzip rather than letting fread auto-detect, which would
  # require R.utils. read_lines() above handles gz transparently via readr.
  vep <- fread(cmd = paste("gunzip -c", shQuote(vep)), sep="\t", skip = "#CHROM", col.names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE"), showProgress = FALSE)
  vep[, variant_row := .I]
  
  # Pull just the CSQ payload from INFO (everything between 'CSQ=' and next ';')
  csq_payload <- str_match(vep$INFO, "CSQ=([^;]+)")[,2]
  
  # Expand to one row per (variant × transcript), then split by '|'
  csq_lst <- lapply(seq_along(csq_payload), function(i) {
    s <- csq_payload[i]
    if (is.na(s)) return(NULL)
    entries <- strsplit(s, ",", fixed = TRUE)[[1]]        # multiple transcripts
    cols <- tstrsplit(entries, "\\|", fixed = FALSE, fill = "")
    x <- as.data.table(cols)
    n <- ncol(x); nm <- csq_names
    if (n > length(nm)) nm <- c(nm, paste0("EXTRA", seq_len(n - length(nm))))
    setnames(x, nm[seq_len(n)])
    x[, variant_row := i][]
  })
  
  csq_table <- rbindlist(csq_lst, use.names = TRUE, fill = TRUE)
  
  # Join back any variant-level columns you want
  csq_keep <- csq_table[ is.na(CANONICAL) | CANONICAL %chin% c("YES","1") ]
  vep <- merge(
    vep[, .(variant_row, CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = ".", FORMAT, SAMPLE)],
    csq_keep,
    by = "variant_row",
    all.x = TRUE
  )
  
  # Extract FORMAT into columns
  vep <- vep %>%
    mutate(
      fmt_keys  = str_split(FORMAT, ":", simplify = FALSE),
      fmt_vals  = str_split(SAMPLE, ":", simplify = FALSE),
      fmt_named = map2(fmt_keys, fmt_vals, ~ rlang::set_names(.y, .x))
    ) %>%
    bind_cols(dplyr::bind_rows(.$fmt_named)) %>%
    select(-fmt_keys, -fmt_vals, -fmt_named)
  
  ## Adjust column content
  # ID
  vep$ID <- paste0(vep$CHROM, "-", vep$POS, "-", vep$REF, "-", vep$ALT)
  
  # HGVS
  vep$HGVSp <- sub("^[^:]*:", "", vep$HGVSp)
  vep$HGVSp <- utils::URLdecode(vep$HGVSp) # decode %3D to = (for synonymous variant)
  vep <- vep %>% mutate(
    HGVS = case_when(
      is.na(SYMBOL) | SYMBOL == "" | is.na(HGVSc) | HGVSc == "" ~ NA,
      is.na(HGVSp) | HGVSp == "" ~ paste0(SYMBOL, ":", HGVSc),
      TRUE ~ paste0(SYMBOL, ":", HGVSc, ":", HGVSp)
    )
  )
  
  # Consequence -> pick most severe
  sev_order <- c(
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "feature_elongation", "feature_truncation", # High impact
    "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant", # Moderate impact
    "splice_donor_5th_base_variant", "splice_region_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant", "start_retained_variant", "stop_retained_variant", "synonymous_variant", # Low impact
    "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant", "coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification", "regulatory_region_variant", "intergenic_variant", "sequence_variant"  # Modifier impact
  )
  sev_rank <- setNames(seq_along(sev_order), sev_order)
  
  pick_most_severe <- function(x, sep="&") {
    sapply(strsplit(ifelse(is.na(x), "", x), sep, fixed = TRUE), function(tokens) {
      tokens <- trimws(tokens)
      tokens <- tokens[tokens != "" & tokens != "."]
      if (!length(tokens)) return(NA_character_)
      r <- sev_rank[tokens]
      r[is.na(r)] <- length(sev_order) + 1L
      tokens[ which.min(r) ]
    }, USE.NAMES = FALSE)
  }
  
  vep$Consequence <- pick_most_severe(vep$Consequence)
  
  # SpliceAI
  vep$SpliceAI_score <- pmax(
    vep$SpliceAI_pred_DS_AG,
    vep$SpliceAI_pred_DS_AL,
    vep$SpliceAI_pred_DS_DG,
    vep$SpliceAI_pred_DS_DL,
    na.rm = TRUE
  )
  vep$SpliceAI_score[is.infinite(vep$SpliceAI_score)] <- NA_real_
  vep$SpliceAI_score <- as.numeric(vep$SpliceAI_score)
  
  ## In-silico predictors
  cell_max <- function(v) sapply(str_split(replace_na(v, ""), fixed("&")), function(tokens) {
    nums <- suppressWarnings(as.numeric(str_trim(tokens[tokens != "" & tokens != "."])))
    if (!length(nums) || all(is.na(nums))) NA_real_ else max(nums, na.rm = TRUE)
  })
  
  cell_negmax <- function(v) sapply(str_split(replace_na(v, ""), fixed("&")), function(tokens) {
    nums <- suppressWarnings(as.numeric(str_trim(tokens[tokens != "" & tokens != "."])))
    if (!length(nums) || all(is.na(nums))) NA_real_ else max(-nums, na.rm = TRUE)
  })
  
  insilico_cols <- c(
    #"SIFT_score", "SIFT4G_score", "Polyphen2_HDIV_score", "Polyphen2_HVAR_score",
    #"MutationTaster_score", "MutationAssessor_score", "PROVEAN_score", "VEST4_score", "MetaSVM_score", "MetaLR_score", "M-CAP_score", "MVP_score", "MPC_score", "PrimateAI_score", "DEOGEN2_score",
    #"BayesDel_addAF_score", "BayesDel_noAF_score", "ClinPred_score", "LIST-S2_score", "CADD_phred", "DANN_score", "fathmm-XF_coding_score", "Eigen-raw_coding", "Eigen-PC-raw_coding",
    #"GERP++_RS", "phyloP100way_vertebrate", "phastCons100way_vertebrate",
    #"AlphaMissense_score", "MetaRNN_score", "REVEL_score",
    "MutPred2_score",
    "VARITY_R_score",
    "ESM1b_score"
  )
  
  negmax_cols <- c(
    #"SIFT_score", "SIFT4G_score", "PROVEAN_score",
    "ESM1b_score"
  )
  
  max_cols <- setdiff(insilico_cols, negmax_cols)
  
  vep <- vep %>%
    mutate(
      # higher = worse
      across(all_of(max_cols),    cell_max),
      
      # lower = worse (SIFT, SIFT4G, PROVEAN, ESM1b)
      across(all_of(negmax_cols), cell_negmax)
    )
  
  
  ### GeneBe
  setwd(VCF)
  
  genebe_command <- sprintf('docker run -v $(pwd)/processing.vcf:/tmp/input.vcf --rm genebe/pygenebe:0.0.18 genebe annotate --genome hg19 --input /tmp/input.vcf --output /dev/stdout %s > analysis_files/%s.genebe.vcf',
                            #"--username puddingyd@gmail.com --api_key ak-f16ru4f3pewnR5VwRdXQ4ztA6"
                            #"--username i54996195@gs.ncku.edu.tw --api_key ak-XgBffTBcx0dx1Uk9KN43SLb2s1VMZUYI"
                            #"--username pikachucym@gmail.com --api_key ak-bhdTsgS4EG47iTULa6MWn9eJu87YkWL9"
                            #"--username changym@myyahoo.com --api_key ak-DWkr5ubcPGEwNamJXzoWBi8aTlT3tqoQ"
                            #"--username ped2021cr@gmail.com --api_key ak-HUbYOUXMHmoklafqpbUSEiQCTuvH4LTr"
                            #"--username n102968@mail.hosp.ncku.edu.tw --api_key ak-ZyntCcxNi6Zd8trukRuI935SPnBG4xpd"
                            #"--username nckuhgenetics@gmail.com --api_key ak-I246LEEYMnA6Vk0X2WnleYW6EYw7nylN"
                            #"--username melissachuang@gmail.com --api_key ak-zqB2bLAtXtf6EXn0F3pwKmve247d7GOm"
                            #"--username complicated.buzzard.bntf@protectsmail.net --api_key ak-DaCueZPihwuEd5CDggJGqvgs74M08TH1"
                            "--username verbal.penguin.qrkd@protectsmail.net --api_key ak-cUhrZSD2SSQdEfrzOpXAoAm3nXCpgfC5"
                            , ID
                            )
  
  if (!file.exists(paste0("analysis_files/", ID, ".genebe.vcf"))) {
    system(genebe_command)
  } 
  
  genebe <- read.table(paste0("analysis_files/", ID, ".genebe.vcf"), header = FALSE, sep = "\t")
  colnames(genebe) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample")
  genebe$ID <- paste0(genebe$`#CHROM`, "-", genebe$POS, "-", genebe$REF, "-", genebe$ALT)
  
  genebe <- genebe %>%
    mutate(ACMG_score = str_extract(INFO, "(?<=acmg_score=)[^;]+")) %>%
    mutate(ACMG_criteria = str_extract(INFO, "(?<=acmg_criteria=)[^;]+")) %>%
    mutate(gnomad_exomes_AF = str_extract(INFO, "(?<=gnomad_exomes_af=)[^;]+"))
  
  
  ### OMIM
  git_pull_repo <- function(repo_dir) {
    if (!dir.exists(repo_dir)) {
      warning("Git repo directory does not exist: ", repo_dir)
      return(FALSE)
    }
    
    cmd <- sprintf("cd %s && git pull", shQuote(repo_dir))
    
    res <- tryCatch(
      {
        out <- system(cmd, intern = TRUE)
        message(paste(out, collapse = "\n"))
        TRUE
      },
      error = function(e) {
        warning("git pull failed, using local files instead: ", e$message)
        FALSE
      }
    )
    
    res
  }
  git_pull_repo("/Users/changym/Desktop/NGS_pipeline/omim")
  
  gene_to_mim <- suppressMessages(read_excel(paste0(WES_pipeline, "omim/OMIM.xlsx")))
  colnames(gene_to_mim)[colnames(gene_to_mim) == "gene_symbol"] <- "Gene.refGeneWithVer"
  gene_to_mim <- gene_to_mim %>%
  mutate(
    split_omim = str_split(OMIM_disease, "\n"),
    split_omim = map(split_omim, ~ .x[1:5]),      # keep at most 5
    Disease1_was_na = is.na(Disease1)             # flag ORIGINAL NA status
  ) %>%
  mutate(
    Disease1 = if_else(Disease1_was_na, map_chr(split_omim, ~ .x[1] %||% NA_character_), Disease1), 
    Disease2 = if_else(Disease1_was_na, map_chr(split_omim, ~ .x[2] %||% NA_character_), Disease2),
    Disease3 = if_else(Disease1_was_na, map_chr(split_omim, ~ .x[3] %||% NA_character_), Disease3),
    Disease4 = if_else(Disease1_was_na, map_chr(split_omim, ~ .x[4] %||% NA_character_), Disease4),
    Disease5 = if_else(Disease1_was_na, map_chr(split_omim, ~ .x[5] %||% NA_character_), Disease5)
  ) %>%
  select(-split_omim, -Disease1_was_na)
  
  ### Annotation
  setwd(VCF)
  vcf_to_annotate <- select(vcf_input, ID, zygosity,
                            AD, total_depth, alt_af,
                            Func.refGeneWithVer, Gene.refGeneWithVer, GeneDetail.refGeneWithVer, ExonicFunc.refGeneWithVer, AAChange.refGeneWithVer,
                            TaiwanBioBank, AF, AF_eas,
                            CLNHGVSn, CLNSIGn, CLNSIGCONFn, CLNREVSTATn,
                            # Old CLNSIG from the input _ann.txt — kept around so the
                            # ClinVar-upgrade detection can compare it to the freshly
                            # joined CLNSIGn. Drops out cleanly if _ann.txt didn't
                            # carry a CLNSIG column on this build.
                            tidyselect::any_of(c("CLNSIG", "CLNSIGCONF", "CLNREVSTAT")),
                            dbscSNV_ADA_SCORE, dbscSNV_RF_SCORE)
  vcf_to_annotate <- left_join(vcf_to_annotate, vep[, c("ID", "SYMBOL", "Consequence", "HGVS", "IMPACT",
                                                        "EXON", "INTRON",
                                                        "SpliceAI_score", "MaxEntScan_ref", "MaxEntScan_alt", "MaxEntScan_diff",
                                                        #"gnomAD4.1_joint_AF","gnomAD4.1_joint_EAS_AF",
                                                        #"SIFT_score","SIFT4G_score","Polyphen2_HDIV_score","Polyphen2_HVAR_score","MutationTaster_score","MutationAssessor_score","PROVEAN_score","VEST4_score","MetaSVM_score","MetaLR_score","M-CAP_score",
                                                        #"MVP_score","MPC_score","PrimateAI_score","DEOGEN2_score","BayesDel_addAF_score","BayesDel_noAF_score","ClinPred_score","LIST-S2_score","CADD_phred","DANN_score","fathmm-XF_coding_score",
                                                        #"Eigen-raw_coding","Eigen-PC-raw_coding","GERP++_RS","phyloP100way_vertebrate","phastCons100way_vertebrate",
                                                        #"AlphaMissense_score","MetaRNN_score","REVEL_score",
                                                        "VARITY_R_score","ESM1b_score", "MutPred2_score")], by = "ID")
  vcf_to_annotate <- left_join(vcf_to_annotate, annovar[, c("ID", 
                                                          #  "AF", "AF_eas", 
                                                            "SIFT_score", "SIFT4G_score", "Polyphen2_HDIV_score", "Polyphen2_HVAR_score", "MutationTaster_score", "MutationAssessor_score", "PROVEAN_score", "VEST4_score", "MetaSVM_score", "MetaLR_score", "M-CAP_score", 
                                                            "MVP_score", "MPC_score", "PrimateAI_score", "DEOGEN2_score", "BayesDel_addAF_score", "BayesDel_noAF_score", "ClinPred_score", "LIST-S2_score", "CADD_phred", "DANN_score", "fathmm-XF_coding_score",
                                                            "Eigen-raw_coding", "Eigen-PC-raw_coding", "GERP++_RS", "phyloP100way_vertebrate", "phastCons100way_vertebrate",
                                                            "AlphaMissense_score", "MetaRNN_score", "REVEL_score")], by = "ID")
  vcf_to_annotate <- left_join(vcf_to_annotate, genebe[, c("ID", "ACMG_score", "ACMG_criteria")], by = "ID") #, relationship = "many-to-many")
  vcf_to_annotate <- left_join(vcf_to_annotate, gene_to_mim[, c("Gene.refGeneWithVer", "OMIM_id", "OMIM_disease", "Disease1", "Disease2", "Disease3", "Disease4", "Disease5")], by = "Gene.refGeneWithVer", relationship = "many-to-many")
  
  vcf_to_annotate <- vcf_to_annotate %>% mutate(HGVS = case_when(is.na(HGVS) & AAChange.refGeneWithVer != "." ~ AAChange.refGeneWithVer, TRUE ~ HGVS))
  vcf_to_annotate <- vcf_to_annotate %>% mutate(Consequence = case_when(is.na(Consequence) & ExonicFunc.refGeneWithVer != "." ~ ExonicFunc.refGeneWithVer, TRUE ~ Consequence))
  
  ## Filter
  vcf <- vcf_to_annotate

  
  
  ## Calculate genotype score
  # BayesDel_addAF score transformation
  vcf <- vcf %>%
    mutate(BayesDel_addAF_score = as.numeric(BayesDel_addAF_score)) %>%
    mutate(BayesDel_addAF_transform = case_when(
      BayesDel_addAF_score <= 0.0692655 ~ 0.4215 * BayesDel_addAF_score + 0.4708,
      BayesDel_addAF_score > 0.0692655 ~ 0.7335 * BayesDel_addAF_score + 0.4492,
      TRUE ~ 0
    ))
  # quadratic equation: y = 0.1670x^2 + 0.5965x + 0.4579
  # piecewise function: For −1.11707 ≤ x < 0.0692655 => y = 0.4215x + 0.4708; For 0.0692655 ≤ x ≤ 0.750927 => y = 0.7335x + 0.4492
  
  # Calcalate AF score
  vcf$AF <- as.numeric(vcf$AF)
  vcf <- vcf %>%
    mutate(AF_score = case_when(
      AF > 0.001 & AF <= 0.01 ~ -0.1,
      AF > 0.0001 & AF <= 0.001 ~ -0.05,
      TRUE ~ 0
    ))
  
  vcf <- vcf %>%
    mutate(TaiwanBioBank = as.numeric(TaiwanBioBank)) %>%
    mutate(TWB_score = case_when(
      TaiwanBioBank > 0.001 & TaiwanBioBank <= 0.01 ~ -0.1,
      TaiwanBioBank > 0.0001 & TaiwanBioBank <= 0.001 ~ -0.05,
      TRUE ~ 0
    ))
  
  # Calculate clinvar score
  vcf <- vcf %>%
    mutate(clinvar_score = case_when(
      CLNSIGn %in% c("Pathogenic", "Pathogenic/Likely_pathogenic") ~ 1,
      CLNSIGn %in% c("Likely_pathogenic") ~ 0.7,
      CLNSIGn %in% c("Likely_benign") ~ -0.2,
      CLNSIGn %in% c("Benign", "Benign/Likely_benign") ~ -0.4,
      grepl("athogenic", CLNSIGCONFn) ~ 0.5,
      TRUE ~ 0 
    ))
  
  # Calculate consequence score -> most severe (https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences)
  vcf <- vcf %>%
    mutate(impact_score = case_when(
      IMPACT %in% c("HIGH") ~ 0.8,
      ExonicFunc.refGeneWithVer %in% c("frameshift insertion", "frameshift deletion", "stopgain", "stoploss", "startloss") ~ 0.8,
      TRUE ~ 0
    ))
  
  # Calculate PTV score
  vcf <- vcf %>%
    mutate(PTV_score = case_when(
      Consequence %in% c("stop_gained", "frameshift_variant", "stop_lost", "start_lost", "splice_acceptor_variant", "splice_donor_variant", "transcript_ablation", "transcript_amplification") ~ 0.8,
      TRUE ~ 0
    ))
  
  # Calculate metascore
  vcf$BayesDel_addAF_transform[vcf$BayesDel_addAF_transform == 0 | vcf$BayesDel_addAF_transform > 1] <- NA
  vcf$MetaRNN_score <- as.numeric(vcf$MetaRNN_score)
  vcf$MetaRNN_score[vcf$MetaRNN_score > 1] <- NA
  vcf$ClinPred_score <- as.numeric(vcf$ClinPred_score)
  vcf$ClinPred_score[vcf$ClinPred_score > 1] <- NA
  vcf$AlphaMissense_score <- as.numeric(vcf$AlphaMissense_score)
  vcf$REVEL_score <- as.numeric(vcf$REVEL_score)
  vcf$BayesDel_noAF_score <- as.numeric(vcf$BayesDel_noAF_score)
  vcf <- vcf %>% mutate(Meta_score = coalesce(MetaRNN_score, BayesDel_addAF_transform, ClinPred_score, 0))
  #vcf <- vcf %>% mutate(Meta_score = 100 * rowMeans(select(., BayesDel_addAF_transform, MetaRNN_score, ClinPred_score), na.rm = TRUE))
  
  
  # Calculate in silico predictors
  vcf$CADD <- vcf$CADD_phred
  vcf$PrimateAI <- vcf$PrimateAI_score
  data.table::setnames(
    vcf,
    old = c("M-CAP_score", "LIST-S2_score",
            "fathmm-XF_coding_score", "Eigen-raw_coding", "Eigen-PC-raw_coding"),
    new = c("M_CAP_score",  "LIST_S2_score",
            "fathmm_XF_coding_score", "Eigen_raw_coding", "Eigen_PC_raw_coding")
  )
  tool_names <- c("SIFT_score", "SIFT4G_score", "PROVEAN_score", 
                  "MutationTaster_score", "MutationAssessor_score",
                  "M_CAP_score", "MutPred2_score", "MVP_score",
                  "PrimateAI_score", "DEOGEN2_score", "Polyphen2_HVAR_score", 
                  "LIST_S2_score", "CADD_phred", "DANN_score", 
                  "fathmm_XF_coding_score", "Eigen_raw_coding", "Eigen_PC_raw_coding")
  vcf[tool_names] <- lapply(vcf[tool_names], as.numeric)
  
  vcf$SIFT_score <- -vcf$SIFT_score
  vcf$SIFT4G_score <- -vcf$SIFT4G_score
  vcf$PROVEAN_score <- -vcf$PROVEAN_score
  
  thresholds_df <- read.table(paste0(WES_pipeline, "tools_score.txt"), header = TRUE, sep = "\t")  # Load the thresholds table
  
  convert_score_to_points <- function(raw_score, tool_name) {  # Define a function to convert raw scores to points for a single tool
    tool_thresholds <- thresholds_df %>% filter(tool == tool_name)
    if (is.na(raw_score)) {
      return(NA)
    } else if (raw_score <= tool_thresholds$BP4_SS) {
      return(-4)
    } else if (raw_score <= tool_thresholds$BP4_M) {
      return(-2)
    } else if (raw_score <= tool_thresholds$BP4_S) {
      return(-1)
    } else if (raw_score < tool_thresholds$VUS) {
      return(0)
    } else if (raw_score < tool_thresholds$PP3_S) {
      return(1)
    } else if (raw_score < tool_thresholds$PP3_M) {
      return(2)
    } else if (raw_score < tool_thresholds$PP3_SS) {
      return(4)
    } else {
      return(NA)
    }
  }
  
  for (tool_name in tool_names) {
    vcf[[tool_name]] <- sapply(vcf[[tool_name]], convert_score_to_points, tool_name = tool_name)  # Iterate over each tool name and convert raw scores to points
  }
  
  patho_subset <- vcf[, tool_names]
  pathogenicity_count <- function(row) { # count the numbers of P/VUS/B
    P_VS <- sum(row == 8, na.rm = TRUE)
    P_SS <- sum(row == 4, na.rm = TRUE)
    P_M <- sum(row == 2, na.rm = TRUE)
    P_S <- sum(row == 1, na.rm = TRUE)
    vus <- sum(row == 0, na.rm = TRUE)
    B_S <- sum(row == -1, na.rm = TRUE)
    B_M <- sum(row == -2, na.rm = TRUE)
    B_SS <- sum(row == -4, na.rm = TRUE)
    B_VS <- sum(row == -8, na.rm = TRUE)
    
    pathogenic <- P_VS * 8 + P_SS * 4 + P_M * 2 + P_S * 1
    benign <- B_VS * 8 + B_SS * 4 + B_M * 2 + B_S * 1
    
    return(paste(pathogenic, "-", vus, "-", benign))
  }
  vcf$in_silico_prediction <- apply(patho_subset, 1, pathogenicity_count)
  vcf$in_silico_score <- rowMeans(vcf[, tool_names], na.rm = TRUE)
  
  # Split and duplicate rows based on Gene.refGene
  #vcf <- vcf %>% separate_rows(Gene.refGeneWithVer, sep = ";")
  
  # Calculate ACMG score and assign classification
  vcf$ACMG_score <- as.numeric(vcf$ACMG_score)
  vcf <- vcf %>%
    mutate(ACMG_classification = case_when(
      ACMG_score >= 10 ~ "Pathogenic",
      ACMG_score <= 9 & ACMG_score >= 6 ~ "Likely pathogenic",
      ACMG_score <= 5 & ACMG_score >= 0 ~ "VUS",
      ACMG_score <= -1 & ACMG_score >= -6 ~ "Likely benign",
      ACMG_score <= -7 ~ "Benign",
      TRUE ~ "VUS"
    ))
  
  vcf <- vcf %>%
    mutate(geno_score_ACMG = 5 * (10 + case_when(ACMG_score >= 10 ~ 10, ACMG_score <= -10 ~ -10, TRUE ~ ACMG_score)))
  
  # Low impact -> minus score
  vcf <- vcf %>% mutate(geno_score_ACMG = case_when(IMPACT %in% c("LOW", "MODIFIER") ~ geno_score_ACMG - 30, TRUE ~ geno_score_ACMG))
  vcf <- vcf %>% mutate(geno_score_ACMG = case_when(geno_score_ACMG < 0 ~ 0, TRUE ~ geno_score_ACMG))
  
  #vcf <- vcf %>% mutate(geno_score = case_when(is.na(geno_score) ~ geno_score_ACMG, TRUE ~ geno_score))
  
  
  # Calculate In-silico predictor score
  # Transform raw score to pathogenicity score (ref: varsome and https://pubmed.ncbi.nlm.nih.gov/40084623/)
  # Benchmark of 65 in-silico tools: https://www.sciencedirect.com/science/article/pii/S0888754325000527
  cols <- c("MetaRNN_score","AlphaMissense_score","REVEL_score","ESM1b_score","VARITY_R_score","BayesDel_noAF_score")
  vcf[cols] <- lapply(vcf[cols], as.numeric)
  
  vcf <- vcf %>%
    mutate(MetaRNN_ps = case_when(
      MetaRNN_score >= 0.939 ~ 4,
      MetaRNN_score >= 0.900 & MetaRNN_score <  0.939 ~ 3,
      MetaRNN_score >= 0.841 & MetaRNN_score <  0.900 ~ 2,
      MetaRNN_score >= 0.748 & MetaRNN_score <  0.841 ~ 1,
      MetaRNN_score >  0.430 & MetaRNN_score <  0.748 ~ 0,
      MetaRNN_score >  0.267 & MetaRNN_score <= 0.430 ~ -1,
      MetaRNN_score >  0.108 & MetaRNN_score <= 0.267 ~ -2,
      MetaRNN_score <= 0.108 ~ -4,
      TRUE ~ NA
    ))
  
  vcf <- vcf %>%
    mutate(AlphaMissense_ps = case_when(
      AlphaMissense_score >= 0.990 ~ 4,
      AlphaMissense_score >= 0.972  & AlphaMissense_score <  0.990 ~ 3,
      AlphaMissense_score >= 0.906  & AlphaMissense_score <  0.972 ~ 2,
      AlphaMissense_score >= 0.792  & AlphaMissense_score <  0.906 ~ 1,
      AlphaMissense_score >  0.170  & AlphaMissense_score <  0.792 ~ 0,
      AlphaMissense_score >  0.100  & AlphaMissense_score <= 0.170 ~ -1,
      AlphaMissense_score >  0.071  & AlphaMissense_score <= 0.100 ~ -2,
      AlphaMissense_score <= 0.071 ~ -3,
      TRUE ~ NA
    ))
  
  vcf <- vcf %>%
    mutate(REVEL_ps = case_when(
      REVEL_score >= 0.932 ~ 4,
      REVEL_score >= 0.895 & REVEL_score <  0.932 ~ 3,
      REVEL_score >= 0.829 & REVEL_score <  0.895 ~ 2,
      REVEL_score >= 0.737 & REVEL_score <  0.829 ~ 1,
      REVEL_score >  0.392 & REVEL_score <  0.737 ~ 0,
      REVEL_score >  0.198 & REVEL_score <= 0.392 ~ -1,
      REVEL_score >  0.032 & REVEL_score <= 0.198 ~ -2,
      REVEL_score >  0.011 & REVEL_score <= 0.032 ~ -2,
      REVEL_score <= 0.011 ~ -4,
      TRUE ~ NA
    ))
  
  vcf <- vcf %>%
    mutate(ESM1b_ps = case_when(
      ESM1b_score >=  24 ~ 4,
      ESM1b_score >=  14.0 & ESM1b_score <   24.0 ~  3,
      ESM1b_score >=  12.1 & ESM1b_score <   14.0 ~  2,
      ESM1b_score >=  10.6 & ESM1b_score <   12.1 ~  1,
      ESM1b_score >   6.3  & ESM1b_score <   10.6 ~  0,
      ESM1b_score >   3.1  & ESM1b_score <=  6.3  ~ -1,
      ESM1b_score >  -8.8  & ESM1b_score <=  3.1  ~ -2,
      ESM1b_score <  -8.8  ~ -3,
      TRUE ~ NA
    ))
  
  vcf <- vcf %>%
    mutate(VARITY_R_ps = case_when(
      VARITY_R_score >= 0.965 ~ 4,
      VARITY_R_score >= 0.915 & VARITY_R_score <  0.965 ~  3,
      VARITY_R_score >= 0.842 & VARITY_R_score <  0.915 ~  2,
      VARITY_R_score >= 0.675 & VARITY_R_score <  0.842 ~  1,
      VARITY_R_score >  0.252 & VARITY_R_score <  0.675 ~  0,
      VARITY_R_score >  0.117 & VARITY_R_score <= 0.252 ~ -1,
      VARITY_R_score >  0.064 & VARITY_R_score <= 0.117 ~ -2,
      VARITY_R_score >  0.037 & VARITY_R_score <= 0.064 ~ -3,
      VARITY_R_score <= 0.037 ~ -4,
      TRUE ~ NA
    ))
  
  vcf <- vcf %>%
    mutate(BayesDel_noAF_ps = case_when(
      BayesDel_noAF_score >=  0.500 ~ 4,
      BayesDel_noAF_score >=  0.410 & BayesDel_noAF_score <   0.500 ~  3,
      BayesDel_noAF_score >=  0.270 & BayesDel_noAF_score <   0.410 ~  2,
      BayesDel_noAF_score >=  0.130 & BayesDel_noAF_score <   0.270 ~  1,
      BayesDel_noAF_score >  -0.179 & BayesDel_noAF_score <   0.130 ~  0,
      BayesDel_noAF_score >  -0.359 & BayesDel_noAF_score <= -0.179 ~ -1,
      BayesDel_noAF_score >  -0.519 & BayesDel_noAF_score <= -0.359 ~ -2,
      BayesDel_noAF_score <= -0.519 ~ -3,
      TRUE ~ NA
    ))
  
  # Calculate pathogenicity score (average of the above 6 predictors)
  cols <- c("MetaRNN_ps","AlphaMissense_ps","REVEL_ps","ESM1b_ps","VARITY_R_ps","BayesDel_noAF_ps")
  
  vcf <- vcf %>%
    rowwise() %>%
    mutate(
      patho_raw = if (all(is.na(c_across(all_of(cols))))) NA_real_
      else mean(c_across(all_of(cols)), na.rm = TRUE),
      patho_n   = sum(!is.na(c_across(all_of(cols)))),
      patho_raw_pen_sqrt = patho_raw * sqrt(patho_n / length(cols)),  # penalty: mean * sqrt ( n of tools with score / n of total tools) 
      patho_score = (pmin(pmax(patho_raw_pen_sqrt, -4), 4) + 4) / 8  # rescale to 0~1
    ) %>%
    ungroup()
  
  
  # Calculate Splice score
  # Transform raw score to pathogenicity score (ref: varsome)
  cols <- c("SpliceAI_score","dbscSNV_ADA_SCORE","dbscSNV_RF_SCORE","MaxEntScan_diff")
  vcf[cols] <- lapply(vcf[cols], as.numeric)
  
  vcf <- vcf %>%
    mutate(SpliceAI_ps = case_when(
      SpliceAI_score >= 0.8 ~ 4,
      SpliceAI_score >= 0.5  & SpliceAI_score <  0.8 ~ 2,
      SpliceAI_score >= 0.2  & SpliceAI_score <  0.8 ~ 1,
      SpliceAI_score <  0.2 ~ 0,
      TRUE ~ NA
    ))
  
  vcf <- vcf %>%
    mutate(dbscSNV_ADA_ps = case_when(
      dbscSNV_ADA_SCORE >= 0.999925 ~ 4,
      dbscSNV_ADA_SCORE >= 0.999322  & dbscSNV_ADA_SCORE <  0.999925 ~ 2,
      dbscSNV_ADA_SCORE >= 0.957813  & dbscSNV_ADA_SCORE <  0.999322 ~ 1,
      dbscSNV_ADA_SCORE <  0.957813 ~ 0,
      TRUE ~ NA
    ))
  
  vcf <- vcf %>%
    mutate(dbscSNV_RF_ps = case_when(
      dbscSNV_RF_SCORE >= 0.994 ~ 4,
      dbscSNV_RF_SCORE >= 0.832  & dbscSNV_RF_SCORE <  0.994 ~ 2,
      dbscSNV_RF_SCORE >= 0.584  & dbscSNV_RF_SCORE <  0.832 ~ 1,
      dbscSNV_RF_SCORE <  0.584 ~ 0,
      TRUE ~ NA
    ))
  
  vcf <- vcf %>%
    mutate(MaxEntScan_ps = case_when(
      abs(MaxEntScan_diff) >= 7.65 ~ 4,
      abs(MaxEntScan_diff) >= 5.96  & abs(MaxEntScan_diff) <  7.65 ~ 2,
      abs(MaxEntScan_diff) >= 4.24  & abs(MaxEntScan_diff) <  5.96 ~ 1,
      abs(MaxEntScan_diff) <  4.24 ~ 0,
      TRUE ~ NA
    ))
  
  # Calculate pathogenicity score (average of the above 6 predictors)
  cols <- c("SpliceAI_ps","dbscSNV_ADA_ps","dbscSNV_RF_ps","MaxEntScan_ps")
  
  vcf <- vcf %>%
    rowwise() %>%
    mutate(
      splice_raw = if (all(is.na(c_across(all_of(cols))))) NA_real_
      else mean(c_across(all_of(cols)), na.rm = TRUE),
      splice_n   = sum(!is.na(c_across(all_of(cols)))),
      splice_raw_pen_sqrt = splice_raw * sqrt(splice_n / length(cols)),  # penalty: mean * sqrt ( n of tools with score / n of total tools) 
      splice_score = splice_raw_pen_sqrt / 4  # rescale to 0~1
    ) %>%
    ungroup()
  
  
  # Calculate genotype score
  cols <- c("clinvar_score","impact_score","Meta_score", "splice_score","AF_score", "TWB_score")
  vcf <- vcf %>% mutate(impact_score = case_when(Consequence %in% c("inframe_deletion", "inframe_insertion") ~ 0.5, TRUE ~ impact_score))
  vcf$geno_score <- rowSums(vcf[ , cols, drop = FALSE], na.rm = TRUE) * 100
  vcf$geno_score <- pmax(pmin(vcf$geno_score, 100), 0)
  
  # Calculate genotype score 2  
  cols <- c("clinvar_score","impact_score","patho_score", "splice_score","AF_score", "TWB_score")
  vcf$geno_score_2 <- rowSums(vcf[ , cols, drop = FALSE], na.rm = TRUE) * 100
  vcf$geno_score_2 <- pmax(pmin(vcf$geno_score_2, 100), 0)
  
  # `vcf` flows directly into the post-PharmCAT block as `variants`; we used to
  # spill it to analysis_files/<ID>_annotated.txt and re-read it there, but the
  # data is unchanged in between so the round-trip is pure disk waste.
  file.remove("processing.vcf", paste0("analysis_files/", ID, ".vep.vcf_warnings.txt"), paste0(ID, ".avinput"), paste0(ID, ".hg19_multianno.vcf"))
}

### Phenotype prioritization - Start here if changing phenotype
{
  setwd(VCF)
  ## Git pull
  git_pull_repo <- function(repo_dir) {
    if (!dir.exists(repo_dir)) {
      warning("Git repo directory does not exist: ", repo_dir)
      return(FALSE)
    }
    
    cmd <- sprintf("cd %s && git pull", shQuote(repo_dir))
    
    res <- tryCatch(
      {
        out <- system(cmd, intern = TRUE)
        message(paste(out, collapse = "\n"))
        TRUE
      },
      error = function(e) {
        warning("git pull failed, using local files instead: ", e$message)
        FALSE
      }
    )
    
    res
  }
  git_pull_repo(paste0(WES_pipeline, "hpo-translator"))

  ## Commit + push a single path in a given git repo; robust against
  ## remote-ahead races (fetch + rebase + retry once). Used by the
  ## phenotype seeder below so a freshly-created _phenotype.txt lands
  ## in the hpo-translator repo automatically.
  git_commit_push_path <- function(repo_dir, path_spec, commit_msg) {
    if (!dir.exists(file.path(repo_dir, ".git"))) {
      warning("push skipped: not a git repo at ", repo_dir)
      return(invisible(FALSE))
    }
    q   <- shQuote(repo_dir)
    pq  <- shQuote(path_spec)
    run <- function(cmd) system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    tryCatch({
      run(sprintf("git -C %s add %s", q, pq))
      changed <- run(sprintf("git -C %s diff --cached --quiet -- %s", q, pq)) != 0
      if (!changed) {
        message("No changes to commit at ", repo_dir, ": ", path_spec)
        return(invisible(TRUE))
      }

      if (run(sprintf("git -C %s commit -m %s", q, shQuote(commit_msg))) != 0) {
        warning("commit failed at ", repo_dir); return(invisible(FALSE))
      }

      br <- suppressWarnings(trimws(system(
        sprintf("git -C %s rev-parse --abbrev-ref HEAD", q), intern = TRUE)))
      if (!length(br) || br == "HEAD") {
        warning("push skipped: detached HEAD in ", repo_dir)
        return(invisible(FALSE))
      }
      br_q <- shQuote(br)

      if (run(sprintf("git -C %s push origin %s", q, br_q)) == 0) {
        message("Pushed to ", repo_dir, ": ", commit_msg); return(invisible(TRUE))
      }

      run(sprintf("git -C %s fetch origin %s", q, br_q))
      if (run(sprintf("git -C %s pull --rebase --autostash origin %s", q, br_q)) != 0) {
        run(sprintf("git -C %s rebase --abort", q))
        warning("rebase failed in ", repo_dir, " (committed locally, will retry next run)")
        return(invisible(FALSE))
      }
      if (run(sprintf("git -C %s push origin %s", q, br_q)) == 0) {
        message("Pushed to ", repo_dir, ": ", commit_msg, " (after rebase)")
        return(invisible(TRUE))
      }
      warning("push failed in ", repo_dir, " (committed locally)")
      invisible(FALSE)
    }, error = function(e) {
      warning("push errored in ", repo_dir, ": ", e$message)
      invisible(FALSE)
    })
  }

  ## Load patient phenotype, with a three-step fallback:
  ##   1) hpo-translator/patient_phenotype/<ID>*_phenotype.txt (canonical)
  ##   2) <VCF>/patient_phenotype.xlsx — if present, seed it back into
  ##      hpo-translator/patient_phenotype/ so future runs hit step 1
  ##   3) Nothing → empty phenotype; Exomiser / LIRICAL get skipped and
  ##      pheno_score stays at 0, but the pipeline still produces webdata.
  pheno_dir  <- file.path(WES_pipeline, "hpo-translator", "patient_phenotype")
  pheno_file <- list.files(
    path = pheno_dir,
    pattern = paste0("^", ID, ".*\\_phenotype\\.txt$"),
    full.names = TRUE
  )

  has_phenotype <- TRUE
  if (length(pheno_file) == 1) {
    patient_phenotype <- read.delim(
      pheno_file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
  } else if (length(pheno_file) > 1) {
    stop(
      "Multiple phenotype files found for ID: ",
      ID, "\n",
      paste(pheno_file, collapse = "\n")
    )
  } else {
    xlsx_fp <- file.path(VCF, "patient_phenotype.xlsx")
    if (file.exists(xlsx_fp)) {
      patient_phenotype <- as.data.frame(read_excel(xlsx_fp))
      # xlsx commonly carries trailing phantom columns (read_excel names
      # them "...4", "...5", etc., or shoves the sample-id header into an
      # empty column) and trailing blank rows. Keep only the meaningful
      # fields in a known order, and drop rows that have no phenotype.
      keep_cols <- intersect(c("phenotype", "hpo_name", "weight"),
                             names(patient_phenotype))
      patient_phenotype <- patient_phenotype[, keep_cols, drop = FALSE]
      patient_phenotype <- patient_phenotype[!is.na(patient_phenotype$phenotype) &
                                               trimws(patient_phenotype$phenotype) != "", ]
      # Panel rows (phenotype doesn't start with "HP:") should leave hpo_name
      # blank, matching the HPO Term Manager layout.
      if ("hpo_name" %in% names(patient_phenotype)) {
        non_hp <- !grepl("^HP:", patient_phenotype$phenotype)
        patient_phenotype$hpo_name[non_hp] <- ""
      }
      dir.create(pheno_dir, showWarnings = FALSE, recursive = TRUE)

      # Look up MRN from webdata/index.json so the seeded filename carries it
      # (easier to browse the hpo-translator directory by patient).
      mrn <- NA_character_
      idx_lookup <- tryCatch(
        jsonlite::fromJSON(file.path(REPO, "webdata", "index.json")),
        error = function(e) NULL
      )
      if (is.data.frame(idx_lookup) && nrow(idx_lookup)) {
        hit <- idx_lookup[as.character(idx_lookup$LIS_ID) == as.character(ID), , drop = FALSE]
        if (nrow(hit)) {
          m <- as.character(hit$MRN[1])
          if (!is.na(m) && nzchar(m)) mrn <- m
        }
      }
      fn_stem   <- if (!is.na(mrn)) paste0(ID, "_", mrn) else ID
      seeded_fp <- file.path(pheno_dir, paste0(fn_stem, "_phenotype.txt"))
      write.table(patient_phenotype, seeded_fp,
                  sep = "\t", row.names = FALSE, quote = FALSE)
      message("Seeded phenotype for ", ID, " from xlsx → ", seeded_fp)
      # Auto-commit + push the seeded file to the hpo-translator repo
      git_commit_push_path(
        repo_dir   = file.path(WES_pipeline, "hpo-translator"),
        path_spec  = file.path("patient_phenotype", paste0(fn_stem, "_phenotype.txt")),
        commit_msg = paste0("patient_phenotype: seed ", ID,
                            if (!is.na(mrn)) paste0(" (MRN ", mrn, ")") else "",
                            " from xlsx")
      )
    } else {
      warning("No phenotype file found for ", ID,
              " (neither ", file.path(pheno_dir, paste0(ID, "_phenotype.txt")),
              " nor ", xlsx_fp, "). Continuing with empty phenotype.")
      patient_phenotype <- data.frame(
        phenotype = character(), weight = numeric(), hpo_name = character(),
        stringsAsFactors = FALSE
      )
      has_phenotype <- FALSE
    }
  }

  patient_phenotype <- patient_phenotype[!is.na(patient_phenotype$phenotype), ]
  if (nrow(patient_phenotype) == 0) has_phenotype <- FALSE
  
  ## Exomiser
  exomiser_tsv <- paste0(VCF, "analysis_files/exomiser_result.variants.tsv")
  if (has_phenotype && !file.exists(exomiser_tsv)) {

  # Make yml file for exomiser
  hp_for_exomiser <- patient_phenotype$phenotype
  hp_for_exomiser <- hp_for_exomiser[grepl("^HP:", hp_for_exomiser)]
  hp_for_exomiser <- paste0("[", paste0("'", hp_for_exomiser, "'", collapse = ", "), "]")

  file.copy("/Users/changym/Desktop/VCF/exomiser_input.yml", paste0(VCF, "analysis_files/"), overwrite = TRUE)

  yml_path <- "analysis_files/exomiser_input.yml"
  x <- readLines(yml_path, warn = FALSE)

  i <- grep("^\\s*hpoIds\\s*:", x)
  indent <- sub("^(\\s*).*$", "\\1", x[i])  # 保留原縮排
  x[i] <- paste0(indent, "hpoIds: ", hp_for_exomiser) # 替換成 inline list

  i <- grep("^\\s*vcf\\s*:", x)
  indent <- sub("^(\\s*).*$", "\\1", x[i])
  x[i] <- paste0(indent, "vcf: ", paste0("analysis_files/", ID, "_hg19.vcf"))

  writeLines(x, yml_path)

  # Run exomiser
  system2(
    "java",
    args = paste0(
      "-Xmx4g ",
      "-Dspring.config.location=/Users/changym/biotools/exomiser/exomiser-cli-14.1.0/application.properties ",
      "-jar /Users/changym/biotools/exomiser/exomiser-cli-14.1.0/exomiser-cli-14.1.0.jar ",
      "--analysis analysis_files/exomiser_input.yml"
    )
  )

  }

  # Load exomiser result (or stub if we had to skip because of missing phenotype)
  if (file.exists(exomiser_tsv)) {
    exomiser_var <- read.delim(exomiser_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    exomiser_var$ID <- paste0("chr", sub("_.*$", "", exomiser_var$ID))
    exomiser_var <- exomiser_var[order(exomiser_var$EXOMISER_GENE_PHENO_SCORE, decreasing = TRUE), ]
    exomiser_var <- exomiser_var[!duplicated(exomiser_var$ID), ]
  } else {
    message("Exomiser skipped for ", ID, " (no phenotype); stub result used.")
    exomiser_var <- data.frame(
      ID = character(),
      EXOMISER_GENE_PHENO_SCORE = numeric(),
      EXOMISER_GENE_COMBINED_SCORE = numeric(),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  
  
  ## LIRICAL
  setwd(VCF)
  lirical_gene_tsv <- "analysis_files/lirical_result.gene.tsv"
  lirical_var_tsv  <- "analysis_files/lirical_result.variant.tsv"
  if (has_phenotype && !file.exists(lirical_var_tsv)) {

  # Make yaml file for LIRICAL
  hp_for_exomiser <- patient_phenotype$phenotype
  hp_for_exomiser <- hp_for_exomiser[grepl("^HP:", hp_for_exomiser)]
  hp_for_exomiser <- paste0("[", paste0("'", hp_for_exomiser, "'", collapse = ", "), "]")

  file.copy("/Users/changym/Desktop/VCF/lirical_input.yaml", paste0(VCF, "analysis_files/"), overwrite = TRUE)

  yml_path <- "analysis_files/lirical_input.yaml"
  x <- readLines(yml_path, warn = FALSE)

  i <- grep("^\\s*hpoIds\\s*:", x)
  indent <- sub("^(\\s*).*$", "\\1", x[i])  # 保留原縮排
  x[i] <- paste0(indent, "hpoIds: ", hp_for_exomiser) # 替換成 inline list

  i <- grep("^\\s*vcf\\s*:", x)
  indent <- sub("^(\\s*).*$", "\\1", x[i])
  x[i] <- paste0(indent, "vcf: ", paste0("analysis_files/", ID, "_hg19.vcf"))

  writeLines(x, yml_path)

  # Run LIRICAL for gene
  system(paste(
    'lirical yaml',
    '-y analysis_files/lirical_input.yaml',
    '-x lirical_result.gene -o analysis_files/ -f html -f tsv'
    ))

  # Run LIRICAL for variant
  system(paste(
    'lirical yaml',
    '-y analysis_files/lirical_input.yaml',
    '--assembly hg19 -ed19 /Users/changym/biotools/exomiser/exomiser-cli-14.1.0/data/2508_hg19',
    '-x lirical_result.variant -o analysis_files/ -f html -f tsv -f'
  ))

  }

  # Load LIRICAL gene result (or stub if no phenotype → no run)
  if (file.exists(lirical_gene_tsv)) {
  x <- readLines(lirical_gene_tsv, warn = FALSE)
  header_idx <- grep("^rank\\t", x)[1]
  lirical_gene <- read.delim(lirical_gene_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, skip = header_idx - 1)

  genemap2 <- read.table(paste0(WES_pipeline, "genemap2.txt"), header = TRUE, sep = "\t", quote = "", fill = TRUE, colClasses = "character")

  # Add gene name and mim number to lirical table
  lirical_gene$pheno_OMIM <- as.character(sub("^OMIM:", "", lirical_gene$diseaseCurie))

  escape_regex <- function(x) gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", x)

  # ---- 同時回填 gene_OMIM 與 gene_symbol ----
  res <- lapply(lirical_gene$pheno_OMIM, function(id) {
    id <- trimws(as.character(id))
    if (is.na(id) || id == "") {
      return(c(gene_OMIM = NA_character_, gene_symbol = NA_character_))
    }

    pat <- paste0("(?<!\\d)", escape_regex(id), "(?!\\d)")
    hit <- grepl(pat, genemap2$Phenotypes, perl = TRUE)

    if (!any(hit, na.rm = TRUE)) {
      return(c(gene_OMIM = NA_character_, gene_symbol = NA_character_))
    }

    mim_vals  <- unique(genemap2$MIM.Number[hit & !is.na(genemap2$MIM.Number) & genemap2$MIM.Number != ""])
    gene_vals <- unique(genemap2$Approved.Gene.Symbol[hit & !is.na(genemap2$Approved.Gene.Symbol) & genemap2$Approved.Gene.Symbol != ""])

    c(
      gene_OMIM   = if (length(mim_vals)  == 0) NA_character_ else paste(mim_vals,  collapse = ";"),
      gene_symbol = if (length(gene_vals) == 0) NA_character_ else paste(gene_vals, collapse = ";")
    )
  })

  res_mat <- do.call(rbind, res)

  lirical_gene$gene_OMIM   <- res_mat[, "gene_OMIM"]
  lirical_gene$gene_symbol <- res_mat[, "gene_symbol"]

  # compositeLR -> clamp to [-10, 10] -> rescale to [0, 100]
  x <- suppressWarnings(as.numeric(lirical_gene$compositeLR))
  x <- pmax(pmin(x, 10), -10)
  lirical_gene$pheno_score_lirical <- as.integer(round((x + 10) / 20 * 100, 0))

  lirical_gene <- lirical_gene %>%
    group_by(gene_symbol) %>%
    summarise(
      pheno_score_lirical = if (all(is.na(pheno_score_lirical))) NA_real_ else max(pheno_score_lirical, na.rm = TRUE),
      .groups = "drop"
    )
  } else {
    message("LIRICAL gene skipped for ", ID, " (no phenotype); stub result used.")
    lirical_gene <- data.frame(
      gene_symbol = character(),
      pheno_score_lirical = numeric(),
      stringsAsFactors = FALSE
    )
  }


  # Load LIRICAL variant result (or stub)
  if (file.exists(lirical_var_tsv)) {
  x <- readLines(lirical_var_tsv, warn = FALSE)
  header_idx <- grep("^rank\\t", x)[1]
  lirical_variant <- read.delim(lirical_var_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, skip = header_idx - 1)

  lirical_variant <- lirical_variant %>%
    mutate(variants = as.character(variants)) %>%
    separate_rows(variants, sep = "\\s*;\\s*") %>%   # 用 ; 展開成多列（自動 trim）
    mutate(
      pathogenicity = str_match(
        variants,
        "pathogenicity\\s*:\\s*([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"
      )[, 2] %>% as.numeric()
    ) %>%
    filter(!is.na(variants), variants != "") %>%      # 去掉空段
    filter(!is.na(pathogenicity)) %>%                # 若你想保留「沒有 pathogenicity」的段落，把這行拿掉
    filter(pathogenicity >= 0.5)                        # 去掉 pathogenicity = 0

  lirical_variant$ID <- sub("\\s.*$", "", as.character(lirical_variant$variants))
  lirical_variant$ID <- gsub("[:>]", "-", lirical_variant$ID) # 1) Replace ":" and ">" with "-"
  lirical_variant$ID <- gsub("([0-9])([A-Za-z])", "\\1-\\2", lirical_variant$ID)   # 2) Add "-" between number and alphabet (both directions, for robustness)
  lirical_variant$ID <- paste0("chr", lirical_variant$ID) # Add "chr"

  lirical_variant <- lirical_variant %>%
    group_by(ID) %>%
    slice_min(order_by = rank, n = 1, with_ties = FALSE) %>%  # 同 ID 取 rank 最小
    ungroup()

  names(lirical_variant)[names(lirical_variant) == "rank"] <- "rank_lirical_variant"
  # Rescale per-variant compositeLR to 0-100 the same way as the gene-
  # level pheno_score_lirical: clamp to [-10, 10] then linear map to
  # [0, 100]. lirical_variant_pathogenicity from the variants TSV often
  # clamps at 1.0 for ClinVar P/LP and isn't useful as a gradient.
  if ("compositeLR" %in% names(lirical_variant)) {
    x <- suppressWarnings(as.numeric(lirical_variant$compositeLR))
    x <- pmax(pmin(x, 10), -10)
    lirical_variant$lirical_variant_score <- as.integer(round((x + 10) / 20 * 100, 0))
  }
  } else {
    message("LIRICAL variant skipped for ", ID, " (no phenotype); stub result used.")
    lirical_variant <- data.frame(
      ID = character(),
      rank_lirical_variant = integer(),
      lirical_variant_score = integer(),
      stringsAsFactors = FALSE
    )
  }

  
  ## HPO score weighted
  if (!file.exists(paste0("analysis_files/", ID, ".pheno.txt"))) {
  
  setwd(WES_pipeline)
  
  ## Load HPO
  hp_db <- read.table("hpo-translator/data/phenotype_to_genes.txt", header = TRUE, sep = "\t", quote = "", fill = TRUE, colClasses = "character")
  hp_db <- select(hp_db, gene_symbol, hpo_id)
  hp_db <- hp_db %>% distinct()
  
  # Load custom panel and combine to HPO
  custom_panels_files <- list.files("hpo-translator/data/gene_panels", pattern = "\\.txt$")
  custom_panels_df <- list()
  
  for (custom_panel in custom_panels_files) {
    panel_name <- sub("\\.txt$", "", custom_panel)
    custom_panels_df[[panel_name]] <- read.table(paste0("hpo-translator/data/gene_panels/", custom_panel), header = FALSE, sep = "\t")
    colnames(custom_panels_df[[panel_name]]) <- c("gene_symbol")
    custom_panels_df[[panel_name]]$hpo_id <- paste(panel_name)
    hp_db <- rbind(custom_panels_df[[panel_name]], hp_db)
    hp_db <- hp_db %>% distinct()
  }
  
  # Create a new table showing the association of each gene with hp/panel
  genes_to_HP <- hp_db %>%
    group_by(gene_symbol) %>%
    summarise(HP = paste(hpo_id, collapse = ", "))
  genes_to_HP <- subset(genes_to_HP, gene_symbol != "-")
  
  
  ## Calculate phenotype score (vectorized via data.table join)
  total_weight <- sum(patient_phenotype$weight)

  if (has_phenotype && total_weight > 0) {
    hp_dt <- as.data.table(hp_db)[gene_symbol != "-"]
    pt <- data.table(hpo_id = patient_phenotype$phenotype,
                     weight = patient_phenotype$weight)

    pheno_score <- hp_dt[pt, on = "hpo_id", nomatch = NULL,
                         allow.cartesian = TRUE
                        ][, .(pheno_score = 100 * sum(weight) / total_weight),
                          by = gene_symbol]
    pheno_score <- pheno_score[pheno_score > 0][order(-pheno_score)]
  } else {
    pheno_score <- data.table(gene_symbol = character(), pheno_score = numeric())
  }

  write.table(pheno_score, paste0(VCF, "analysis_files/", ID, ".pheno.txt"), sep = "\t", row.names = FALSE)
  }
  
  ## Incorporate phenotype score
  setwd(VCF)
  # In-memory hand-off from the annotation block above (was previously
  # round-tripped through _annotated.txt). as.data.frame so subsequent
  # `colnames(variants)[…] <- …` doesn't mutate the upstream object by ref.
  variants <- as.data.frame(vcf)
  pheno_score <- read.table(paste0("analysis_files/", ID, ".pheno.txt"), sep = "\t", header = TRUE,
                            colClasses = c(gene_symbol = "character", pheno_score = "numeric"))
  
  colnames(variants)[colnames(variants) == "Gene.refGeneWithVer"] <- "gene_symbol"
  variants <- left_join(variants, pheno_score[, c("gene_symbol", "pheno_score")], by = "gene_symbol")
  variants$pheno_score[is.na(variants$pheno_score)] <- 0
  
  variants <- left_join(variants, exomiser_var[, c("ID", "EXOMISER_GENE_PHENO_SCORE", "EXOMISER_GENE_COMBINED_SCORE")], by = "ID") # incorporate exomiser result
  variants$EXOMISER_GENE_PHENO_SCORE <- as.numeric(variants$EXOMISER_GENE_PHENO_SCORE)
  variants$EXOMISER_GENE_PHENO_SCORE[is.na(variants$EXOMISER_GENE_PHENO_SCORE)] <- 0
  
  colnames(variants)[colnames(variants) == "EXOMISER_GENE_PHENO_SCORE"] <- "pheno_score_exomiser"
  variants$pheno_score_exomiser <- 100 * variants$pheno_score_exomiser
  
  variants$EXOMISER_GENE_COMBINED_SCORE <- as.numeric(variants$EXOMISER_GENE_COMBINED_SCORE)
  colnames(variants)[colnames(variants) == "EXOMISER_GENE_COMBINED_SCORE"] <- "total_score_exomiser_variant"
  
  variants <- left_join(variants, lirical_gene[, c("gene_symbol", "pheno_score_lirical")], by = "gene_symbol")

  # LIRICAL per-variant join — pathogenicity (0-1) + variant rank — so the
  # webdata can show each variant's LIRICAL result alongside the in-house
  # scores. Rename `pathogenicity` so it doesn't shadow the column added
  # earlier in the in-silico aggregator.
  lirical_variant_keep <- intersect(c("ID", "rank_lirical_variant", "lirical_variant_score"),
                                    names(lirical_variant))
  variants <- left_join(variants, lirical_variant[, lirical_variant_keep], by = "ID")
  
  write.table(variants, gzfile(paste0("analysis_files/", ID, "_anno_combined.txt.gz")), sep = "\t", row.names = FALSE)

  # Rank-comparison table: refilter to the high-impact / disease-relevant
  # subset, compute total scores under the in-house, Exomiser and LIRICAL
  # weighting schemes, then convert each total to a rank. Useful for
  # benchmarking the in-house ranking against Exomiser / LIRICAL on the
  # same sample. (Restored from pre-refactor versions of this script.)
  rank_test <- variants %>%
    filter(
      (
        IMPACT %in% c("HIGH", "MODERATE") |
          (Func.refGeneWithVer %in% c("exonic", "splicing") & ExonicFunc.refGeneWithVer != "synonymous SNV") |
          Consequence %in% c("splice_donor_5th_base_variant", "splice_region_variant", "splice_donor_region_variant") |
          CLNSIGn %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic") |
          grepl("P", CLNSIGCONFn)
      ) &
        !gene_symbol %in% c("ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "HTT", "TBP", "ZFHX3", "THAP11", "JPH3", "PABPN1", "HLA-A", "HLA-B") &
        !is.na(OMIM_disease)
    )

  # rank_lirical_variant is already on `variants` (joined upstream into
  # webdata). Pull only compositeLR here — joining rank_lirical_variant
  # again would collide and leave the column as `.x` / `.y`-suffixed.
  rank_test <- left_join(rank_test, lirical_variant[, c("ID", "compositeLR")], by = "ID")
  colnames(rank_test)[colnames(rank_test) == "compositeLR"] <- "total_score_lirical_variant"

  add_na0 <- function(a, b) rowSums(cbind(a, b), na.rm = TRUE)
  rank_test$total_score_single_weighted    <- add_na0(rank_test$geno_score,   rank_test$pheno_score)
  rank_test$total_score_multiple_weighted  <- add_na0(rank_test$geno_score_2, rank_test$pheno_score)
  rank_test$total_score_single_exomiser    <- add_na0(rank_test$geno_score,   rank_test$pheno_score_exomiser)
  rank_test$total_score_multiple_exomiser  <- add_na0(rank_test$geno_score_2, rank_test$pheno_score_exomiser)
  rank_test$total_score_single_lirical     <- add_na0(rank_test$geno_score,   rank_test$pheno_score_lirical)
  rank_test$total_score_multiple_lirical   <- add_na0(rank_test$geno_score_2, rank_test$pheno_score_lirical)

  score_cols <- c(
    "total_score_single_weighted",
    "total_score_multiple_weighted",
    "total_score_single_exomiser",
    "total_score_multiple_exomiser",
    "total_score_single_lirical",
    "total_score_multiple_lirical",
    "total_score_exomiser_variant"
  )

  for (col in score_cols) {
    x <- rank_test[[col]]
    if (!is.numeric(x)) {
      x <- suppressWarnings(as.numeric(as.character(x)))
      rank_test[[col]] <- x
    }
    rank_col <- sub("^total_score_", "rank_", col)
    rank_test[[rank_col]] <- rank(-rank_test[[col]], ties.method = "min", na.last = "keep")
  }

  rank_test <- select(rank_test, "ID", "HGVS",
                      "geno_score", "geno_score_2", "pheno_score", "pheno_score_exomiser", "total_score_exomiser_variant", "pheno_score_lirical", "total_score_lirical_variant",
                      "total_score_single_weighted", "total_score_multiple_weighted", "total_score_single_exomiser", "total_score_multiple_exomiser", "total_score_single_lirical", "total_score_multiple_lirical",
                      "rank_single_weighted", "rank_multiple_weighted", "rank_single_exomiser", "rank_multiple_exomiser", "rank_single_lirical", "rank_multiple_lirical", "rank_exomiser_variant", "rank_lirical_variant")

  write.table(rank_test, paste0("analysis_files/", ID, "_rank.txt"), sep = "\t", row.names = FALSE)

  ### Write webdata JSON for the online analysis tool (hg19 build)
  {
    LIS_ID <- ID

    # Column aliases so the same scoring/filter code works for hg19 (CLNSIGn etc.)
    # and matches the key names the web tool expects (CLNSIG / CLNSIGCONF / CLNREVSTAT).
    col_clnsig     <- "CLNSIGn"
    col_clnsigconf <- "CLNSIGCONFn"
    col_clnrevstat <- "CLNREVSTATn"

    # ClinVar review stars from the review-status string. Computed for
    # both the new (CLNREVSTATn) and — if present — the old (CLNREVSTAT
    # from _ann.txt) review status, so the web tool can render the
    # original ClinVar with its own star count next to the upgrade arrow.
    stars_from_revstat <- function(s) {
      dplyr::case_when(
        s == "practice_guideline"                                 ~ 4L,
        s == "reviewed_by_expert_panel"                           ~ 3L,
        grepl("multiple_submitters,_no_conflicts", s)             ~ 2L,
        grepl("single_submitter|conflicting_classifications", s)  ~ 1L,
        TRUE                                                       ~ 0L
      )
    }
    variants <- variants %>%
      mutate(
        clinvar_stars     = stars_from_revstat(.data[[col_clnrevstat]]),
        clinvar_stars_old = if ("CLNREVSTAT" %in% names(.)) stars_from_revstat(CLNREVSTAT) else NA_integer_
      )

    # Use the in-silico-ensemble variant score (patho_score-based) to match
    # hg38, which folds the MetaRNN/REVEL/AlphaMissense/ESM1b/VARITY_R/
    # BayesDel_noAF consensus (via patho_score) instead of MetaRNN alone.
    # Both hg19 geno_score_2 and hg38 geno_score emit the same scoring
    # philosophy — overwrite here so the web receives the ensemble value
    # under the same "geno_score" JSON key.
    variants$geno_score <- variants$geno_score_2

    # total_score: simple sum used by the web's ranking (same as hg38)
    variants$total_score_web <- rowSums(cbind(coalesce(variants$geno_score, 0),
                                              coalesce(variants$pheno_score, 0)), na.rm = TRUE)

    gene_panels_dir <- paste0(NGS_pipeline, "hpo-translator/data/gene_panels/")
    read_panel <- function(fn) {
      p <- file.path(gene_panels_dir, fn)
      if (!file.exists(p)) { warning("Gene panel not found: ", p); return(character(0)) }
      unique(trimws(readLines(p, warn = FALSE)))
    }
    acmg_sf_panel   <- read_panel("ACMG_SF_v3.3.txt")
    proactive_panel <- read_panel("proactive.txt")
    carrier_panel   <- read_panel("carrier_mackenzie_1300+.txt")

    # Genes flagged as "interested" in genes_interested.xlsx (cell A1, comma-separated)
    candidate_genes <- tryCatch({
      f <- file.path(VCF, "genes_interested.xlsx")
      if (!file.exists(f)) character(0) else {
        df <- suppressMessages(read_excel(f, sheet = 1, col_names = FALSE))
        if (is.null(df) || ncol(df) == 0 || nrow(df) == 0) character(0) else {
          val <- df[[1]][1]
          if (is.na(val)) character(0) else {
            g <- strsplit(as.character(val), ",\\s*|，\\s*")[[1]]
            unique(trimws(g[trimws(g) != ""]))
          }
        }
      }
    }, error = function(e) character(0))

    # Repeat-expansion / HLA genes excluded from Diagnostic & Pathogenic
    repeat_genes <- c(
      "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7",
      "HTT", "TBP", "ZFHX3", "THAP11", "JPH3", "PABPN1",
      "HLA-A", "HLA-B"
    )

    # Parse CLNSIGCONF like "Uncertain_significance(1)|Likely_benign(1)|Pathogenic(2)"
    clnsigconf_sum <- function(s, labels) {
      if (is.na(s) || s == "") return(0L)
      parts <- strsplit(s, "\\|")[[1]]
      tot <- 0L
      for (p in parts) {
        m <- regmatches(p, regexec("^(.+?)\\((\\d+)\\)$", p))[[1]]
        if (length(m) == 3 && m[2] %in% labels) tot <- tot + as.integer(m[3])
      }
      tot
    }
    plp_labels <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")
    blb_labels <- c("Benign", "Likely_benign", "Benign/Likely_benign")

    # ClinVar significance → 4-tier ranking, where lower index = more
    # pathogenic. Used by the upgrade detector below.
    #   T1: Pathogenic / Likely_pathogenic / Pathogenic/Likely_pathogenic
    #   T2: Uncertain_significance / Conflicting_*
    #   T3: Benign / Likely_benign / Benign/Likely_benign
    #   T4: not_provided / no_classification_* / "" / "." / NA  (and any other)
    tier_clnsig <- function(s) {
      s <- ifelse(is.na(s) | s == "" | s == ".", "missing", as.character(s))
      result <- rep(4L, length(s))
      result[grepl("Benign", s) | grepl("Likely_benign", s)] <- 3L
      result[grepl("Uncertain_significance", s) | grepl("Conflicting", s)] <- 2L
      result[grepl("Pathogenic", s) | grepl("Likely_pathogenic", s)] <- 1L
      result
    }

    web_variants <- variants %>%
      mutate(
        is_clinvar_plp = .data[[col_clnsig]] %in% plp_labels,
        clnsigconf_plp_n = vapply(.data[[col_clnsigconf]], clnsigconf_sum, integer(1), labels = plp_labels),
        clnsigconf_blb_n = vapply(.data[[col_clnsigconf]], clnsigconf_sum, integer(1), labels = blb_labels),
        is_clinvar_conflict = !is.na(.data[[col_clnsig]]) &
          .data[[col_clnsig]] == "Conflicting_classifications_of_pathogenicity",
        is_clinvar_conflict_plp = is_clinvar_conflict & clnsigconf_plp_n > 0,
        is_clinvar_conflict_plp_majority =
          is_clinvar_conflict & clnsigconf_plp_n > 0 & clnsigconf_plp_n >= clnsigconf_blb_n,
        is_acmg_plp = ACMG_classification %in% c("Pathogenic", "Likely pathogenic"),
        is_acmg_blb = ACMG_classification %in% c("Benign", "Likely benign"),
        # ClinVar upgrade detection. Skipped (NA) if _ann.txt didn't carry
        # an old CLNSIG column. Two rules now (narrowed per user spec):
        #   ↑↑ : new = T1 (LP/P/P+LP) and old was T2/T3/T4
        #   ↑  : old was {B, LB, B/LB, plain VUS, or Conflicting NOT
        #        majority pathogenic} AND new is Conflicting with
        #        majority pathogenic (the existing
        #        is_clinvar_conflict_plp_majority flag on the new side)
        clinvar_old_tier = if ("CLNSIG" %in% names(.)) tier_clnsig(CLNSIG) else NA_integer_,
        clinvar_new_tier = tier_clnsig(.data[[col_clnsig]]),
        old_clnsigconf_plp_n = if (all(c("CLNSIG", "CLNSIGCONF") %in% names(.))) {
          vapply(CLNSIGCONF, clnsigconf_sum, integer(1), labels = plp_labels)
        } else NA_integer_,
        old_clnsigconf_blb_n = if (all(c("CLNSIG", "CLNSIGCONF") %in% names(.))) {
          vapply(CLNSIGCONF, clnsigconf_sum, integer(1), labels = blb_labels)
        } else NA_integer_,
        old_is_conflict = "CLNSIG" %in% names(.) & !is.na(CLNSIG) &
                          CLNSIG == "Conflicting_classifications_of_pathogenicity",
        old_is_conflict_plp_majority = old_is_conflict &
          !is.na(old_clnsigconf_plp_n) & old_clnsigconf_plp_n > 0 &
          !is.na(old_clnsigconf_blb_n) & old_clnsigconf_plp_n >= old_clnsigconf_blb_n,
        old_qualifies_for_arrow_up =
          ("CLNSIG" %in% names(.)) & (
            (clinvar_old_tier == 3L) |
            (!is.na(CLNSIG) & CLNSIG == "Uncertain_significance") |
            (old_is_conflict & !old_is_conflict_plp_majority)
          ),
        clinvar_upgrade = dplyr::case_when(
          is.na(clinvar_old_tier) ~ NA_character_,
          clinvar_new_tier == 1L & clinvar_old_tier >= 2L ~ "↑↑",
          is_clinvar_conflict_plp_majority &
            old_qualifies_for_arrow_up ~ "↑",
          TRUE ~ NA_character_
        )
      )

    cand_df <- web_variants %>%
      dplyr::filter(
        gene_symbol %in% candidate_genes & (
          IMPACT %in% c("HIGH", "MODERATE") |
          (Func.refGeneWithVer %in% c("exonic", "splicing") &
             !(ExonicFunc.refGeneWithVer %in% c("synonymous SNV"))) |
          Consequence %in% c("splice_donor_5th_base_variant",
                             "splice_region_variant",
                             "splice_donor_region_variant") |
          is_clinvar_plp | is_clinvar_conflict_plp
        )
      ) %>%
      arrange(desc(total_score_web))

    # Diagnostic (hg19): pass the impact / ClinVar pre-filter (excluding
    # STR / HLA genes), then keep only variants with pheno_score > 0,
    # sort by total_score, take top 30, and cluster same-gene variants
    # so the Report reads row by gene.
    diag_df <- web_variants %>%
      dplyr::filter(
        (
          IMPACT %in% c("HIGH", "MODERATE") |
          (Func.refGeneWithVer %in% c("exonic", "splicing") &
             !(ExonicFunc.refGeneWithVer %in% c("synonymous SNV"))) |
          Consequence %in% c("splice_donor_5th_base_variant",
                             "splice_region_variant",
                             "splice_donor_region_variant") |
          is_clinvar_plp |
          is_clinvar_conflict_plp
        ) &
        !(gene_symbol %in% repeat_genes)
      ) %>%
      dplyr::filter(!is.na(pheno_score) & pheno_score > 0) %>%
      arrange(desc(total_score_web)) %>%
      slice_head(n = 30) %>%
      mutate(gene_first = match(gene_symbol, unique(gene_symbol))) %>%
      arrange(gene_first, desc(total_score_web)) %>%
      select(-gene_first)

    patho_df <- web_variants %>%
      dplyr::filter(
        !(gene_symbol %in% repeat_genes) & !is_acmg_blb & (
          (!is.na(ACMG_score) & ACMG_score >= 4) |
          is_clinvar_plp |
          is_clinvar_conflict_plp_majority
        )
      ) %>%
      arrange(desc(total_score_web))

    # ClinVar upgraded — variants whose ClinVar classification has been
    # re-graded since the original _ann.txt. Top half (now LP/P, ↑↑)
    # ahead of bottom half (now VUS/Conflict, ↑); within each half,
    # sort by total_score_web desc.
    upgraded_df <- web_variants %>%
      dplyr::filter(!is.na(clinvar_upgrade) & !(gene_symbol %in% repeat_genes)) %>%
      arrange(clinvar_upgrade == "↑", desc(total_score_web))

    panel_hits <- function(panel) {
      web_variants %>%
        dplyr::filter(gene_symbol %in% panel & (is_clinvar_plp | is_acmg_plp)) %>%
        arrange(desc(total_score_web))
    }
    sf_df        <- panel_hits(acmg_sf_panel)
    proactive_df <- panel_hits(proactive_panel)
    carrier_df   <- panel_hits(carrier_panel)

    # ID (e.g. "chr1-1234-C-A") → per-variant CHROM/POS/REF/ALT for the web.
    # hg19 drops CHROM/POS/REF/ALT earlier in the pipeline, so recover from ID.
    split_id_parts <- function(ids) {
      m <- stringr::str_match(ids, "^([^-]+)-([0-9]+)-([^-]+)-([^-]+)$")
      list(CHROM = m[, 2], POS = suppressWarnings(as.integer(m[, 3])),
           REF = m[, 4], ALT = m[, 5])
    }

    slim_fields <- function(df) {
      parts <- split_id_parts(df$ID)
      df %>% transmute(
        id                  = ID,
        CHROM               = parts$CHROM,
        POS                 = parts$POS,
        REF                 = parts$REF,
        ALT                 = parts$ALT,
        HGVS                = HGVS,
        gene_symbol         = gene_symbol,
        Consequence         = Consequence,
        total_score         = total_score_web,
        geno_score          = geno_score,
        pheno_score         = pheno_score,
        zygosity            = zygosity,
        AD                  = AD,
        alt_af              = suppressWarnings(as.numeric(alt_af)),
        total_depth         = suppressWarnings(as.integer(total_depth)),
        exon_or_intron      = dplyr::case_when(
          !is.na(EXON)   & EXON   != "" ~ paste0("Exon ",   sub("/.*$", "", EXON)),
          !is.na(INTRON) & INTRON != "" ~ paste0("Intron ", sub("/.*$", "", INTRON)),
          TRUE                          ~ NA_character_
        ),
        CLNSIG              = .data[[col_clnsig]],
        CLNSIGCONF          = .data[[col_clnsigconf]],
        CLNREVSTAT          = .data[[col_clnrevstat]],
        clinvar_stars       = clinvar_stars,
        ACMG_classification = ACMG_classification,
        ACMG_criteria       = ACMG_criteria,
        ACMG_score          = ACMG_score,
        AlphaMissense_score = AlphaMissense_score,
        MetaRNN_score       = MetaRNN_score,
        SpliceAI_score      = SpliceAI_score,
        in_silico_prediction = in_silico_prediction,
        MaxEntScan_diff      = suppressWarnings(as.numeric(MaxEntScan_diff)),
        AF                  = suppressWarnings(as.numeric(AF)),
        AF_eas              = suppressWarnings(as.numeric(AF_eas)),
        TaiwanBioBank       = suppressWarnings(as.numeric(TaiwanBioBank)),
        OMIM_id             = OMIM_id,
        OMIM_disease        = OMIM_disease,
        Disease1 = Disease1, Disease2 = Disease2, Disease3 = Disease3,
        Disease4 = Disease4, Disease5 = Disease5,
        # ClinVar upgrade marker — "↑↑" when the variant moved up to the
        # P/LP tier since the original _ann.txt, "↑" for the narrower
        # VUS/Conflict-not-majority-P → Conflict-majority-P transition,
        # NA otherwise. Drives both the filter for the new ClinVar
        # upgraded section and the inline arrow prefix on the HGVS line.
        clinvar_upgrade = if ("clinvar_upgrade" %in% names(.)) clinvar_upgrade else NA_character_,
        # Old CLNSIG / CLNSIGCONF / stars from _ann.txt, kept around so
        # the web tool can render the previous ClinVar entry next to the
        # current one with the same compact "LP(1)|VUS(2)(★)" formatter.
        CLNSIG_old        = if ("CLNSIG"     %in% names(.)) CLNSIG     else NA_character_,
        CLNSIGCONF_old    = if ("CLNSIGCONF" %in% names(.)) CLNSIGCONF else NA_character_,
        clinvar_stars_old = if ("clinvar_stars_old" %in% names(.)) clinvar_stars_old else NA_integer_,
        # Exomiser / LIRICAL per-variant results, exposed for the More
        # extras panel. rank_exomiser_variant is computed in-house from
        # total_score_exomiser_variant; rank_lirical_variant comes from
        # LIRICAL's own TSV. Variants Exomiser / LIRICAL didn't see end up
        # NA and the web UI just shows "—".
        pheno_score_exomiser           = pheno_score_exomiser,
        # Exomiser combined score is 0-1 in the source TSV; rescale to
        # 0-100 integer for parity with the in-house pheno_score and
        # the LIRICAL variant score (also 0-100 below).
        total_score_exomiser_variant   = as.integer(round(suppressWarnings(as.numeric(total_score_exomiser_variant)) * 100, 0)),
        rank_exomiser_variant          = rank(-suppressWarnings(as.numeric(total_score_exomiser_variant)),
                                              ties.method = "min", na.last = "keep"),
        pheno_score_lirical            = pheno_score_lirical,
        # Per-variant compositeLR rescaled to 0-100 (computed upstream).
        lirical_variant_score          = if ("lirical_variant_score" %in% names(.)) as.integer(lirical_variant_score) else NA_integer_,
        rank_lirical_variant           = if ("rank_lirical_variant" %in% names(.)) as.integer(rank_lirical_variant) else NA_integer_
      )
    }

    keep_ids <- unique(c(cand_df$ID, diag_df$ID, patho_df$ID, upgraded_df$ID,
                         sf_df$ID, proactive_df$ID, carrier_df$ID))
    combined_slim <- slim_fields(web_variants %>% dplyr::filter(ID %in% keep_ids))
    variants_map <- setNames(
      lapply(seq_len(nrow(combined_slim)), function(i) as.list(combined_slim[i, ])),
      combined_slim$id
    )

    idx_df <- tryCatch(
      jsonlite::fromJSON(file.path(REPO, "webdata", "index.json")),
      error = function(e) data.frame()
    )
    meta_row <- if (is.data.frame(idx_df) && nrow(idx_df)) idx_df[idx_df$LIS_ID == LIS_ID, ] else data.frame()
    meta <- if (nrow(meta_row)) as.list(meta_row[1, ])
            else list(LIS_ID = LIS_ID, Name = NA, MRN = NA, Test = NA, Category = NA)

    pheno_list <- tryCatch({
      df <- patient_phenotype %>%
        as.data.frame() %>%
        mutate(across(everything(), as.character)) %>%
        mutate(weight = suppressWarnings(as.numeric(weight)))
      if (!"hpo_name" %in% names(df)) df$hpo_name <- NA_character_
      df %>% mutate(
        label = dplyr::case_when(
          !is.na(hpo_name) & hpo_name != "" ~ hpo_name,
          grepl("^HP:", phenotype)          ~ NA_character_,
          TRUE                              ~ phenotype
        )
      )
    }, error = function(e) data.frame(phenotype = character(), weight = numeric(),
                                      hpo_name = character(), label = character()))

    web_obj <- list(
      meta              = meta,
      genome_build      = GENOME_BUILD,   # "hg19" — web uses this for the [hg19] tag
      clinvar_date      = clinvar_date,
      patient_phenotype = pheno_list,
      # Gene panel for the diagnostic report's "本次檢測基因包括" section.
      # Sourced from the pheno-score table (every gene with score > 0),
      # not from the variant set, so the listing covers genes the patient
      # has no variant in too.
      pheno_genes       = tryCatch({
        ps <- read.table(paste0("analysis_files/", ID, ".pheno.txt"),
                         sep = "\t", header = TRUE,
                         colClasses = c(gene_symbol = "character", pheno_score = "numeric"))
        as.list(sort(unique(ps$gene_symbol[ps$pheno_score > 0])))
      }, error = function(e) list()),
      variants          = variants_map,
      categories        = list(
        candidate         = as.list(cand_df$ID),
        diagnostic        = as.list(diag_df$ID),
        pathogenic        = as.list(patho_df$ID),
        clinvar_upgraded  = as.list(upgraded_df$ID),
        acmg_sf           = as.list(sf_df$ID),
        proactive         = as.list(proactive_df$ID),
        carrier           = as.list(carrier_df$ID)
      ),
      generated_at      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    )

    dir.create(file.path(REPO, "webdata", "variants"), showWarnings = FALSE, recursive = TRUE)
    jsonlite::write_json(
      web_obj,
      file.path(REPO, "webdata", "variants", paste0(LIS_ID, ".json")),
      auto_unbox = TRUE, na = "null", pretty = TRUE
    )

    # Push the updated webdata to GitHub (same resilient helper as hg38).
    git_push_webdata <- function(repo_dir, lis_id) {
      if (!dir.exists(file.path(repo_dir, ".git"))) {
        warning("webdata push skipped: not a git repo at ", repo_dir)
        return(invisible(FALSE))
      }
      q   <- shQuote(repo_dir)
      run <- function(cmd) system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

      tryCatch({
        run(sprintf("git -C %s add webdata", q))
        changed <- run(sprintf("git -C %s diff --cached --quiet -- webdata", q)) != 0
        if (changed) {
          if (run(sprintf("git -C %s commit -m %s", q,
                          shQuote(paste0("webdata: update ", lis_id, " (hg19)")))) != 0) {
            warning("webdata commit failed for ", lis_id); return(invisible(FALSE))
          }
        } else {
          message("webdata unchanged for ", lis_id)
        }

        br <- suppressWarnings(trimws(system(
          sprintf("git -C %s rev-parse --abbrev-ref HEAD", q), intern = TRUE)))
        if (!length(br) || br == "HEAD") {
          warning("webdata push skipped: detached HEAD in ", repo_dir)
          return(invisible(FALSE))
        }
        br_q <- shQuote(br)

        if (run(sprintf("git -C %s push origin %s", q, br_q)) == 0) {
          message("Pushed webdata for ", lis_id); return(invisible(TRUE))
        }

        run(sprintf("git -C %s fetch origin %s", q, br_q))
        if (run(sprintf("git -C %s pull --rebase --autostash origin %s", q, br_q)) != 0) {
          run(sprintf("git -C %s rebase --abort", q))
          warning("webdata rebase failed for ", lis_id,
                  " (committed locally, will retry next run)")
          return(invisible(FALSE))
        }
        if (run(sprintf("git -C %s push origin %s", q, br_q)) == 0) {
          message("Pushed webdata for ", lis_id, " (after rebase)")
          return(invisible(TRUE))
        }
        warning("webdata push failed for ", lis_id, " (committed locally)")
        invisible(FALSE)
      }, error = function(e) {
        warning("webdata push errored for ", lis_id, ": ", e$message)
        invisible(FALSE)
      })
    }
    git_push_webdata(REPO, LIS_ID)
  }
}

system(paste("bgzip", paste0(ID, "_ann.txt")))
unlink(Filter(dir.exists, Sys.glob("*_variants_files")), recursive = TRUE)

}  # end for (ID in IDs)
