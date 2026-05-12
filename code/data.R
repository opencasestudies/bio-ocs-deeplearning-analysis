#===============================================================================
#
#  PROGRAM: data.R
#
#  AUTHOR:  Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: Full preprocessing pipeline used to construct the training and 
#           testing datasets for the deep learning epigenetic clock case study.
#
#           This script orchestrates the following end-to-end workflow:
#
#             1. Define GEO training and testing studies
#             2. Retrieve phenotype metadata from GEO
#             3. Download supplementary files (raw IDATs when available)
#             4. Process raw methylation arrays using sesame
#             5. Apply probe filtering and platform harmonization
#             6. Parse and standardize age metadata
#             7. Align phenotype records with methylation matrices
#             8. Construct shared CpG feature space across studies
#             9. Combine study-level datasets into train/test matrices
#            10. Produce scorcher-ready training and testing objects
#
#  DEPENDS: This script sources helper functions defined in:
#
#             code/data_helpers.R
#
#           These helpers provide reusable functionality for:
#
#             - GEO metadata retrieval
#             - Raw IDAT discovery and preprocessing
#             - Probe filtering and CpG harmonization
#             - Phenotype-methylation alignment
#             - Age extraction and rule-based parsing
#             - Matrix combination and imputation
#
#           The helper file must be available in the project directory
#           before running this script.
#
#  INPUT:   GEO DNA methylation studies (Illumina 27K and 450K arrays)
#
#           Training Studies:
#
#             GSE106648, GSE125105, GSE128235, GSE19711, GSE27044, GSE30870,
#             GSE40279,  GSE41037,  GSE52588,  GSE53740, GSE58119, GSE67530,
#             GSE77445,  GSE77696,  GSE81961,  GSE84624, GSE97362
#
#           Testing Studies:
#
#             GSE102177, GSE103911, GSE105123, GSE107459, GSE107737, GSE112696,
#             GSE34639,  GSE37008,  GSE59065,  GSE61496,  GSE79329,  GSE87582,
#             GSE87640,  GSE98876,  GSE99624
#
#           Optional reference files located in:
#
#               data/reference/
#
#             - cross_reactive_probes.csv
#             - age_rules.csv
#
#  OUTPUT:  Preprocessed phenotype and methylation datasets stored in:
#
#               data/processed/
#
#           Key outputs include:
#
#             - pheno_train_final / pheno_test_final
#             - beta_train_cpg_x_sample / beta_test_cpg_x_sample
#             - scorcher_train_test.qs
#
#           The final object contains:
#
#             x_train  (samples x CpGs)
#             y_train  (chronological age)
#             x_test
#             y_test
#
#           along with metadata and preprocessing manifests.
#
#  NOTES:   Please note the following:
#
#             1. Raw array processing
#
#                 - When IDAT files are available, methylation beta values are
#                   computed using `sesame::openSesame()` with default settings.
#                 - When IDAT files are unavailable, GEO series matrices are
#                   used as a fallback and transformed to beta scale when
#                   necessary.
#
#             2. Probe harmonization
#
#                 - The feature space is restricted to CpGs shared across
#                   Illumina 27K and 450K arrays.
#                 - Sex chromosome probes and cross-reactive probes are removed.
#
#             3. Age parsing
#
#                 - Age is extracted from GEO phenotype metadata using curated
#                   study-specific rules when available.
#                 - A heuristic parser provides a fallback when no rule is
#                   defined.
#                 - Integer ages are shifted by +0.5 years to match the
#                   convention used in the reference preprocessing pipeline.
#
#             4. Reproducibility
#
#                 - Intermediate processing manifests and QC summaries are
#                   written to disk so each step of the pipeline can be
#                   inspected or reproduced.
#
#  UPDATED: 2026-04-27
#
#===============================================================================

#=== PRE-SETUP PROJECT CHECKS ==================================================

message("\nPRE-SETUP PROJECT CHECKS\n")

#--- 1. CHECK FOR CORRECT WORKING DIRECTORY ------------------------------------ 

wd_found <- FALSE

expected_marker <- "code/data.R"

if (file.exists(expected_marker)) {
  wd_found <- TRUE
  message(" + Working directory correct (found ", expected_marker, ")")
}

if (!wd_found) {
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    tryCatch({
      wd_try <- rstudioapi::getActiveProject()
      if (!is.null(wd_try) && file.exists(file.path(wd_try, expected_marker))) {
        setwd(wd_try)
        wd_found <- TRUE
        message(" + RStudio project detected: ", wd_try)
      }
    }, error = function(e) NULL)
  }
}

if (!wd_found) {
  current <- getwd()
  for (i in 1:3) {
    if (file.exists(file.path(current, expected_marker))) {
      setwd(current)
      wd_found <- TRUE
      message(" + Project root found at: ", current)
      break
    }
    current <- dirname(current)
  }
}

if (!wd_found) {
  stop("\nx Cannot locate project root.\n\n",
       "Expected to find: code/data.R\n",
       "Current directory: ", getwd(), "\n\n",
       "Please do one of:\n",
       "  1. Run this script from the project root directory\n",
       "  2. In RStudio, open the .Rproj file first\n",
       "  3. Manually: setwd('/path/to/bio-ocs-deeplearning-analysis')\n",
       "               source('code/data.R')")
}

#--- 2. VERIFY REQUIRED DIRECTORIES --------------------------------------------

message("\n + Checking directory structure...")

required_dirs <- c("code", "data")

for (dir in required_dirs) {
  if (!dir.exists(dir)) {
    stop("x Required directory not found: ", dir)
  }
}

#--- 3. SET UP DATA DIRECTORY PATHS --------------------------------------------

#  data/geoquery/           - GEO series matrix files (getGEO cache)
#  data/raw/downloads/      - Downloaded supplemental archives (unchanged)
#  data/raw/extracted/      - First-pass extraction (tar/tar.gz/tgz only)
#  data/raw/tar_extracted/  - Secondary extraction (_RAW.tar files)
#  data/raw/idat/           - Final .idat files (from .idat.gz or nested tars)
#  data/intermediate/       - Processing checkpoints
#  data/processed/          - Final outputs
#  data/reference/          - Reference files (cross_reactive_probes, age_rules)

DIR_GEO_PHENO      <- "data/geoquery"
DIR_RAW_DOWNLOADS  <- "data/raw/downloads"
DIR_RAW_EXTRACTED  <- "data/raw/extracted"
DIR_RAW_TAR        <- "data/raw/tar_extracted"
DIR_RAW_IDAT       <- "data/raw/idat"
DIR_INTER          <- "data/intermediate"
DIR_OUT            <- "data/processed"
DIR_REF            <- "data/reference"

fs::dir_create(DIR_GEO_PHENO,     recurse = TRUE)
fs::dir_create(DIR_RAW_DOWNLOADS, recurse = TRUE)
fs::dir_create(DIR_RAW_EXTRACTED, recurse = TRUE)
fs::dir_create(DIR_RAW_TAR,       recurse = TRUE)
fs::dir_create(DIR_RAW_IDAT,      recurse = TRUE)
fs::dir_create(DIR_INTER,         recurse = TRUE)
fs::dir_create(DIR_OUT,           recurse = TRUE)
fs::dir_create(DIR_REF,           recurse = TRUE)

CROSS_REACTIVE_FILE <- fs::path(DIR_REF, "cross_reactive_probes.csv")
AGE_RULES_FILE      <- fs::path(DIR_REF, "age_rules.csv")

#--- 4. REFERENCE FILE VALIDATION ---------------------------------------------- ### NEED TO IMPLEMENT

message(" + Checking reference files...")

ref_files <- list(
  age_rules = "data/reference/age_rules.csv",
  cross_reactive = "data/reference/cross_reactive_probes.csv"
)

ref_status <- list()

for (name in names(ref_files)) {
  path <- ref_files[[name]]
  exists <- file.exists(path)
  ref_status[[name]] <- exists

  if (exists) {
    size <- file.size(path) / 1024
    message("   + ", basename(path), " (", round(size, 0), " KB)")
  } else {
    message("   x ", basename(path), " (MISSING - Optional but recommended)")
  }
}

#--- 5. DISK SPACE CHECK -------------------------------------------------------

message(" + Checking available disk space...")

data_dir <- "data"

disk_check <- tryCatch({
  df <- fs::dir_info(data_dir, recurse = FALSE)
  used_gb <- sum(df$size, na.rm = TRUE) / 1024^3
  used_gb < 500  # Warn if > 500GB already used in data dir
}, error = function(e) TRUE)

if (!disk_check) {
  warning(" x Warning: Over 500GB already used in data/ directory.\n",
          "   Consider cleaning up old runs or expanding storage.")
}

#--- 6. REFERENCE FILE VALIDATION FOR DATA QUALITY -----------------------------

message(" + Validating reference files required for accuracy...")

missing_refs <- c()

if (!file.exists(CROSS_REACTIVE_FILE)) {
  missing_refs <- c(missing_refs, "cross_reactive_probes.csv")
  message("   x", CROSS_REACTIVE_FILE, " (will use all probes)")
} else {
  message("   + cross_reactive_probes.csv found")
}

if (!file.exists(AGE_RULES_FILE)) {
  missing_refs <- c(missing_refs, "age_rules.csv")
  message("   x ", AGE_RULES_FILE, " (will use heuristic age parsing)")
} else {
  message("   + age_rules.csv found")
}

if (length(missing_refs) > 0) {
  warning("\n x WARNING: Missing reference files:\n",
          "  - ", paste(missing_refs, collapse="\n  - "), "\n",
          "\n  For better results, provide:\n",
          "    data/reference/age_rules.csv\n",
          "    data/reference/cross_reactive_probes.csv\n")
}

message("\n + Project setup checks complete\n")

#=== SETUP ===================================================================== 

#--- SESSION OPTIONS -----------------------------------------------------------

options(timeout = 7200,
        repos = c(CRAN = "https://cloud.r-project.org"),
        BioC_mirror = "https://bioconductor.org")

#--- SET SEED FOR REPRODUCIBILITY ----------------------------------------------

set.seed(1)

#--- PACKAGE MANAGEMENT --------------------------------------------------------

cran_pkgs <- c(
  "data.table",
  "qs2",
  "fs",
  "Matrix",
  "tidyverse",
  "R.utils"
)

missing_cran <- cran_pkgs[!vapply(
  cran_pkgs,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

invisible(lapply(cran_pkgs, library, character.only = TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "GEOquery",
  "Biobase",
  "sesame",
  "sesameData",
  "minfi",
  "BiocParallel",
  "IlluminaHumanMethylation27kanno.ilmn12.hg19",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)

missing_bioc <- bioc_pkgs[!vapply(
  bioc_pkgs,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]

if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

invisible(lapply(bioc_pkgs, library, character.only = TRUE))

#--- SESAME DATA CACHING -------------------------------------------------------

# NOTE: The sesame package requires cached reference data to process IDAT files.
#       Download and cache this data once at the start.
#       This must complete successfully before IDAT processing can work.

message("\nPreparing sesame data cache (one-time setup)...\n")

sesame_cache_ok <- FALSE

tryCatch({
  ExperimentHub::ExperimentHub()
  sesameData::sesameDataCache()
  message("   + sesame data cache ready")
  sesame_cache_ok <- TRUE
}, error = function(e) {
  message("   x sesame data cache setup failed:")
  message("    Error: ", conditionMessage(e))
})

if (!sesame_cache_ok) {
  stop("\n x CRITICAL: sesame data cache unavailable.\n\n",
       "This is required for processing raw IDAT methylation files.\n\n",
       "Try the following:\n",
       "  1. In R console, run:\n",
       "       install.packages('BiocManager')\n",
       "       BiocManager::install('sesameData')\n",
       "       sesameData::sesameDataCache()\n",
       "  2. Restart R\n",
       "  3. Re-run this script\n\n",
       "Note: If this remains unavailable, IDAT files will be skipped\n",
       "and the pipeline will fall back to GEO series matrix data.")
}

#--- HELPER FUNCTIONS ----------------------------------------------------------

source("code/data_helpers.R")

#=== DEFINE DATA ===============================================================

#--- DEFINE STUDIES ------------------------------------------------------------

trn_gses <- c(
  "GSE106648", "GSE125105", "GSE128235", "GSE19711", "GSE27044", "GSE30870",
  "GSE40279",  "GSE41037",  "GSE52588",  "GSE53740", "GSE58119", "GSE67530",
  "GSE77445",  "GSE77696",  "GSE81961",  "GSE84624", "GSE97362"
)

tst_gses <- c(
  "GSE102177", "GSE103911", "GSE105123", "GSE107459", "GSE107737", "GSE112696",
  "GSE34639",  "GSE37008",  "GSE59065",  "GSE61496",  "GSE79329",  "GSE87582",
  "GSE87640",  "GSE98876",  "GSE99624"
)

all_gses <- c(trn_gses, tst_gses)

#--- LOAD REFERENCE OBJECTS ---------------------------------------------------- ### NEED TO IMPLEMENT

AGE_RULES <- load_age_rules(AGE_RULES_FILE)

COMMON_CPGS <- common_cpgs_27k_450k(
  drop_sex = TRUE,
  drop_cross_reactive = TRUE,
  cross_reactive_file = CROSS_REACTIVE_FILE,
  save_manifest = TRUE
)

message("Fixed common CpG panel size: ", length(COMMON_CPGS))

#--- FETCH RAW PHENOTYPE TABLES ------------------------------------------------

message("\nFetching raw phenotype tables...\n")

ph_trn_raw <- dplyr::bind_rows(lapply(trn_gses, fetch_gse_pheno))
ph_tst_raw <- dplyr::bind_rows(lapply(tst_gses, fetch_gse_pheno))

safe_checkpoint_save(ph_trn_raw, fs::path(DIR_INTER, "ph_trn_raw.qs"))
safe_checkpoint_save(ph_tst_raw, fs::path(DIR_INTER, "ph_tst_raw.qs"))

#--- DOWNLOAD AND EXTRACT SUPPLEMENTARY FILES ----------------------------------

#  STEP 1: Download all supplementary archives from GEO
#           -> data/raw/downloads/{GSE_ID}/ (archives keep original names)
#
#  STEP 2: Extract top-level archives (tar/tar.gz/tgz)
#           -> data/raw/extracted/{GSE_ID}/ (archives stay in downloads/)
#
#  STEP 3: Extract nested TAR files (e.g., _RAW.tar)
#           -> data/raw/tar_extracted/{GSE_ID}/ (nested extraction)
#
#  STEP 4: Extract .idat.gz files
#           -> data/raw/idat/{GSE_ID}/ (final decompressed IDATs)

message("\nStarting supplementary file data extraction...\n")

message("\nStep 1: Downloading supplementary files...")

invisible(lapply(all_gses, download_geo_raw)) 

message("\nStep 2: Extracting top-level archives...")

invisible(lapply(all_gses, extract_geo_raw))

message("\nStep 3: Extracting nested TAR files...")

invisible(lapply(all_gses, extract_nested_tar_files))

message("\nStep 4: Extracting .idat.gz files...")

invisible(lapply(all_gses, extract_idat_gz_files))

cross_reactive <- load_cross_reactive_probes(CROSS_REACTIVE_FILE)                ### NEED TO IMPLEMENT

if (length(cross_reactive) > 0) {
  message(length(cross_reactive), " cross-reactive probes to remove.")
} else {
  message("No cross-reactive probe file found; skipping that filter.")
}

#=== PREPROCESS ALL STUDIES ====================================================

# NOTE: Each step processes all studies independently and defines checkpoints.
#
#         - First run: Run the script normally. All steps will execute and save.
#         - Subsequent runs: Modify only the steps you need to change, e.g.,
#
#            SKIP_STEP_1 <- TRUE  # Load phenotype checkpoint (default)
#            SKIP_STEP_2 <- FALSE # Rerun this step
#
# NOTE: You can modify, rerun, and debug each step without affecting the others.

message("\nProcessing all studies...\n")

studies <- setNames(vector("list", length(all_gses)), all_gses)

#--- STEP 1: FETCH PHENOTYPE DATA ----------------------------------------------

# Retrieve raw phenotype metadata from GEO for each study.
# Extract and parse age information using heuristic fallback.
#
# NOTE: Phenotype data is loaded from previously-fetched checkpoints (ph_trn_raw
#       and ph_tst_raw) to avoid redundant GEO API calls. This step can be 
#       skipped if a full checkpoint exists. Set SKIP_STEP_1 <- FALSE to rerun,
#       or TRUE to load saved results.

SKIP_STEP_1 <- FALSE

STEP_1_CHECKPOINT <- fs::path(DIR_INTER, "step_1_phenotype_checkpoint.qs")

if (SKIP_STEP_1 && file.exists(STEP_1_CHECKPOINT)) {                             ### QC CODE CHECKS UP TO HERE

  message("\n[STEP 1/6] Loading phenotype data from checkpoint...")

  studies_loaded <- safe_checkpoint_read(STEP_1_CHECKPOINT)

  if (!is.null(studies_loaded)) {
    studies <- studies_loaded
    message("   + Step 1 checkpoint loaded: ", length(all_gses), " studies")
  } else {
    message("   x Checkpoint corrupted, rerunning STEP 1...")
    SKIP_STEP_1 <- FALSE
  }
}

if (!SKIP_STEP_1 || is.null(studies_loaded)) {

  message("\n[STEP 1/6] Processing phenotype data for all studies...")
  
  # Combine previously-fetched phenotype data (ph_trn_raw and ph_tst_raw)
  # to avoid redundant GEO API calls
  message("   Building phenotype lookup table from cached data...")
  
  ph_all <- dplyr::bind_rows(ph_trn_raw, ph_tst_raw)
  ph_lookup <- split(ph_all, f = ph_all$gse, drop = TRUE)
  
  message("   Phenotype data available for ", 
    length(unique(ph_all$gse)), " unique studies")

  for (gse_id in all_gses) {

    message("  [", match(gse_id, all_gses), "/", length(all_gses), "] ", gse_id, 
      " - Processing phenotype data...", appendLF = FALSE)

    # Use cached phenotype data instead of refetching from GEO
    if (gse_id %in% names(ph_lookup)) {
      ph <- ph_lookup[[gse_id]]
      message(" (", nrow(ph), " samples)", appendLF = FALSE)
    } else {
      message(" (phenotype data not found in cache)", appendLF = FALSE)
      ph <- data.frame()  # Empty dataframe if not found
    }
    
    # Extract age information with progress indication
    message(" -> extracting age...", appendLF = FALSE)
    
    ph$age <- extract_age_robust(ph)
    
    n_age_success <- sum(!is.na(ph$age))
    pct_age_success <- if (nrow(ph) > 0) 
      round(100 * n_age_success / nrow(ph), 1) else 0
    
    message(" (", n_age_success, "/", nrow(ph), " = ", pct_age_success, "% parsed)")
    
    # Debug: Show available columns if age parsing failed
    if (pct_age_success == 0 && nrow(ph) > 0) {
      char_cols <- names(ph)[vapply(ph, function(z) is.character(z) || is.factor(z), logical(1))]
      num_cols <- names(ph)[vapply(ph, is.numeric, logical(1))]
      age_like_char <- char_cols[stringr::str_detect(tolower(char_cols), "age|gestational|donor_age|chronological")]
      age_like_num <- num_cols[stringr::str_detect(tolower(num_cols), "age|gestational|donor_age|chronological")]
      
      if (length(age_like_char) > 0 || length(age_like_num) > 0) {
        message("    [DEBUG] Available age-like columns: ",
          paste(c(age_like_char, age_like_num), collapse = ", "))
      } else {
        message("    [DEBUG] No age-like columns found. Available columns: ",
          paste(head(names(ph), 10), collapse = ", "),
          if (ncol(ph) > 10) " ... and more" else "")
      }
    }

    studies[[gse_id]]$pheno <- ph

    studies[[gse_id]]$n_pheno_total <- nrow(ph)

    studies[[gse_id]]$n_age_parsed <- n_age_success
  }

  message("   + Step 1 complete: Phenotype data processed for ",

    length(all_gses), " studies")

  message("  Saving checkpoint for future runs...")

  safe_checkpoint_save(studies, STEP_1_CHECKPOINT)

  message("  Checkpoint saved to: ", STEP_1_CHECKPOINT)
}

#--- STEP 2: READ METHYLATION DATA ---------------------------------------------

# Load methylation beta values from IDAT files (primary method).
# Fallback to GEO series matrix if IDATs unavailable.
# Normalize all values to beta scale [0, 1].
#
# NOTE: IDAT files are extracted in the supplementary file extraction steps
#       (STEPS 1-4 above) and stored in data/raw/idat/{GSE_ID}/.

SKIP_STEP_2 <- FALSE                                                             

STEP_2_CHECKPOINT <- fs::path(DIR_INTER, "step_2_methylation_checkpoint.qs")

if (SKIP_STEP_2 && file.exists(STEP_2_CHECKPOINT)) {

  message("\n[STEP 2/6] Loading methylation data from checkpoint...")

  studies_loaded <- safe_checkpoint_read(STEP_2_CHECKPOINT)

  if (!is.null(studies_loaded)) {
    studies <- studies_loaded
    message("   + Step 2 checkpoint loaded")
  } else {
    message("   x Checkpoint corrupted, rerunning STEP 2...")
    SKIP_STEP_2 <- FALSE
  }
}

if (!SKIP_STEP_2 || is.null(studies_loaded)) {

  message("\n[STEP 2/6] Reading methylation data for all studies...")

  for (gse_id in all_gses) {

    message("  [", match(gse_id, all_gses), "/", length(all_gses), "] ", gse_id, 
      " - Reading methylation data...", appendLF = FALSE)

    gse_dir <- fs::path(DIR_RAW_IDAT, gse_id)

    # Try IDAT files first
    beta <- read_beta_from_idats(gse_dir)

    method <- "idat_sesame_openSesame_beta"

    # Fallback to series matrix if IDATs unavailable
    if (is.null(beta)) {
      message(" (no IDATs found, fallback to series matrix)", appendLF = FALSE)
      expr <- read_expr_from_series_matrix(gse_id)
      beta <- normalize_expr_fallback_to_beta(expr)
      method <- "series_matrix_exprs_fallback"
    }

    if (is.null(beta)) {
      message(" [FAILED]")
      message("    [!] No methylation data available; marking as failed")
      studies[[gse_id]]$meth <- NULL
      studies[[gse_id]]$method <- NA_character_

    } else {
      message(" (OK) (", nrow(beta), " probes x ", ncol(beta), " samples)")
      studies[[gse_id]]$meth <- beta
      studies[[gse_id]]$method <- method
    }
  }

  message("   + Step 2 complete: Methylation data loaded")

  # Save checkpoint for future runs

  message("  Saving checkpoint for future runs...")

  safe_checkpoint_save(studies, STEP_2_CHECKPOINT)
}

#--- STEP 3: FILTER PROBES -----------------------------------------------------

# Remove problematic probes:
#   3a. Deduplicate probe rows (average if same probe appears multiple times)
#   3b. Remove non-autosomal probes                                              ### NEED TO IMPLEMENT
#   3c. Remove cross-reactive / multi-mapping probes                             ### NEED TO IMPLEMENT

SKIP_STEP_3 <- FALSE

STEP_3_CHECKPOINT <- fs::path(DIR_INTER, "step_3_probe_filtering_checkpoint.qs")

if (SKIP_STEP_3 && file.exists(STEP_3_CHECKPOINT)) {
  message("\n[STEP 3/6] Loading probe-filtered data from checkpoint...")
  studies_loaded <- safe_checkpoint_read(STEP_3_CHECKPOINT)
  if (!is.null(studies_loaded)) {
    studies <- studies_loaded
    message("   + Step 3 checkpoint loaded")
  } else {
    message("   x Checkpoint corrupted, rerunning STEP 3...")
    SKIP_STEP_3 <- FALSE
  }
}

if (!SKIP_STEP_3 || is.null(studies_loaded)) {
  message("\n[STEP 3/6] Filtering probes for all studies...")

  for (gse_id in all_gses) {
    message("  ", gse_id)

    # Skip studies with no methylation data
    
    if (is.null(studies[[gse_id]]$meth)) {
      message("    [!] Skipping (no methylation data)")
      next
    }

    beta <- studies[[gse_id]]$meth

    # 3a. Deduplicate probe rows
    
    beta <- dedup_probe_rows(beta)
    message("    - After deduplication: ", nrow(beta), " probes")

    # 3b. Remove non-autosomal probes
    
    beta <- drop_non_autosomal_probes(beta)
    message("    - After removing non-autosomal: ", nrow(beta), " probes")

    # 3c. Remove cross-reactive probes
    
    beta_before <- nrow(beta)
    beta <- drop_cross_reactive_probes(beta, cross_reactive = cross_reactive)
    beta_after <- nrow(beta)
    removed <- beta_before - beta_after
    message("    - After removing cross-reactive: ", 
            nrow(beta), " probes (removed ", removed, ")")

    studies[[gse_id]]$meth <- beta
  }

  message("   + Step 3 complete: Probes filtered")

  # Save checkpoint for future runs
  
  message("  Saving checkpoint for future runs...")
  safe_checkpoint_save(studies, STEP_3_CHECKPOINT)
}

#--- STEP 4: ALIGN PHENOTYPE TO METHYLATION MATRIX -----------------------------

# Match phenotypes to methylation matrix columns by normalized sample keys.
# Remove unmatched samples.
# Standardize column names to GEO accession IDs.

SKIP_STEP_4 <- FALSE

STEP_4_CHECKPOINT <- fs::path(DIR_INTER, "step_4_alignment_checkpoint.qs")

if (SKIP_STEP_4 && file.exists(STEP_4_CHECKPOINT)) {
  message("\n[STEP 4/6] Loading aligned data from checkpoint...")
  studies_loaded <- safe_checkpoint_read(STEP_4_CHECKPOINT)
  if (!is.null(studies_loaded)) {
    studies <- studies_loaded
    message("   + Step 4 checkpoint loaded")
  } else {
    message("   x Checkpoint corrupted, rerunning STEP 4...")
    SKIP_STEP_4 <- FALSE
  }
}

if (!SKIP_STEP_4 || is.null(studies_loaded)) {
  message("\n[STEP 4/6] Aligning phenotype to methylation matrix...")

  for (gse_id in all_gses) {
    message("  [", match(gse_id, all_gses), "/", length(all_gses), "] ", gse_id, 
      " - Aligning samples...", appendLF = FALSE)

    # Skip studies with no methylation data
    
    if (is.null(studies[[gse_id]]$meth)) {
      message(" [SKIPPED - no methylation data]")
      next
    }

    ph <- studies[[gse_id]]$pheno
    mat <- studies[[gse_id]]$meth
    
    n_pheno_start <- nrow(ph)
    n_meth_start <- ncol(mat)

    aligned <- align_pheno_to_matrix(ph, mat)
    ph2 <- aligned$ph
    mat2 <- aligned$mat

    n_matched <- ncol(mat2)
    
    if (n_matched == 0) {
      message(" [WARNING: 0 samples matched - check sample identifiers]")
    } else {
      message(" (OK) (", n_matched, "/", n_meth_start, " methylation samples matched)")
    }

    studies[[gse_id]]$pheno <- ph2
    studies[[gse_id]]$meth <- mat2
  }

  message("   + Step 4 complete: Phenotype and methylation aligned")

  # Save checkpoint for future runs
  
  message("  Saving checkpoint for future runs...")
  safe_checkpoint_save(studies, STEP_4_CHECKPOINT)
}

#--- STEP 5: SUBSET BY AGE AVAILABILITY ----------------------------------------

# Retain only samples with successfully parsed age values.
# Ensure 1:1 correspondence between matrix columns and phenotype rows.

SKIP_STEP_5 <- FALSE

STEP_5_CHECKPOINT <- fs::path(DIR_INTER, "step_5_age_subset_checkpoint.qs")

if (SKIP_STEP_5 && file.exists(STEP_5_CHECKPOINT)) {
  message("\n[STEP 5/6] Loading age-subset data from checkpoint...")
  studies_loaded <- safe_checkpoint_read(STEP_5_CHECKPOINT)
  if (!is.null(studies_loaded)) {
    studies <- studies_loaded
    message("   + Step 5 checkpoint loaded")
  } else {
    message("   x Checkpoint corrupted, rerunning STEP 5...")
    SKIP_STEP_5 <- FALSE
  }
}

if (!SKIP_STEP_5 || is.null(studies_loaded)) {
  message("\n[STEP 5/6] Subsetting to samples with age data...")

  for (gse_id in all_gses) {
    message("  ", gse_id)

    # Skip studies with no methylation data
    
    if (is.null(studies[[gse_id]]$meth)) {
      message("    [!] Skipping (no methylation data)")
      next
    }

    ph <- studies[[gse_id]]$pheno
    mat <- studies[[gse_id]]$meth

    # Keep only samples with age
    
    keep <- !is.na(ph$age)
    n_before <- nrow(ph)
    n_after <- sum(keep)

    ph <- ph[keep, , drop = FALSE]
    mat <- mat[, keep, drop = FALSE]

    message("    - Samples with age: ", n_after, " / ", n_before)

    studies[[gse_id]]$pheno <- ph
    studies[[gse_id]]$meth <- mat
  }

  message("   + Step 5 complete: Age subsetting done")

  # Save checkpoint for future runs
  
  message("  Saving checkpoint for future runs...")
  safe_checkpoint_save(studies, STEP_5_CHECKPOINT)
}

#--- STEP 6: CREATE GLOBALLY UNIQUE SAMPLE IDS ---------------------------------

# Prefix sample IDs with GSE accession to guarantee uniqueness across studies.
# Validate alignment between phenotype and matrix.
# Store final processed results.

message("\n[STEP 6/6] Creating globally unique sample identifiers...")

proc_all <- setNames(vector("list", length(all_gses)), all_gses)

for (gse_id in all_gses) {
  message("  ", gse_id)

  ph <- studies[[gse_id]]$pheno
  mat <- studies[[gse_id]]$meth

  # Handle case with no methylation data
  
  if (is.null(mat)) {
    message("    [!] Skipping (no methylation data)")
    proc_all[[gse_id]] <- list(
      gse = gse_id,
      method = studies[[gse_id]]$method,
      ph = ph[0, , drop = FALSE],
      mat = NULL,
      n_parsed_age = studies[[gse_id]]$n_age_parsed,
      n_total_pheno = studies[[gse_id]]$n_pheno_total
    )
    next
  }

  # Create globally unique sample IDs: GSE:GSM
  
  ph$sample_id <- paste0(gse_id, ":", ph$geo_accession)
  colnames(mat) <- ph$sample_id

  # Validate alignment
  
  if (!identical(colnames(mat), ph$sample_id)) {
    stop("Study ", gse_id, " (STEP 6): Sample ID mismatch!\n",
         "  Matrix has ", ncol(mat), " samples\n",
         "  Phenotype has ", nrow(ph), " samples\n",
         "  Column names don't match. Check STEP 4 alignment.")
  }
  if (ncol(mat) != nrow(ph)) {
    stop("Study ", gse_id, " (STEP 6): Dimension mismatch!\n",
         "  Matrix columns: ", ncol(mat), "\n",
         "  Phenotype rows: ", nrow(ph), "\n",
         "  These must be equal. Check STEP 4 alignment.")
  }

  message("    - Final: ", nrow(ph), " samples x ", nrow(mat), " probes")

  # Store final processed results
  
  proc_all[[gse_id]] <- list(
    gse = gse_id,
    method = studies[[gse_id]]$method,
    ph = ph,
    mat = mat,
    n_parsed_age = studies[[gse_id]]$n_age_parsed,
    n_total_pheno = studies[[gse_id]]$n_pheno_total
  )
}

message("   + Step 6 complete: Unique IDs assigned\n")

message("All studies preprocessed!")

#--- STEP 7: CREATE TRAINING AND TESTING DATA SPLITS ---------------------------

proc_trn <- proc_all[trn_gses]
proc_tst <- proc_all[tst_gses]

# Per-split manifests before harmonization

train_manifest <- tibble::tibble(
  gse = names(proc_trn),
  method = vapply(proc_trn, `[[`, character(1), "method"),
  n_pheno_total = vapply(proc_trn, `[[`, numeric(1), "n_total_pheno"),
  n_age_parsed = vapply(proc_trn, `[[`, numeric(1), "n_parsed_age"),
  n_cpg_raw = vapply(proc_trn, 
    function(x) if (is.null(x$mat)) NA_integer_ else nrow(x$mat), integer(1)),
  n_samples_final = vapply(proc_trn, 
    function(x) if (is.null(x$mat)) NA_integer_ else ncol(x$mat), integer(1))
)

test_manifest <- tibble::tibble(
  gse = names(proc_tst),
  method = vapply(proc_tst, `[[`, character(1), "method"),
  n_pheno_total = vapply(proc_tst, `[[`, numeric(1), "n_total_pheno"),
  n_age_parsed = vapply(proc_tst, `[[`, numeric(1), "n_parsed_age"),
  n_cpg_raw = vapply(proc_tst, 
    function(x) if (is.null(x$mat)) NA_integer_ else nrow(x$mat), integer(1)),
  n_samples_final = vapply(proc_tst, 
    function(x) if (is.null(x$mat)) NA_integer_ else ncol(x$mat), integer(1))
)

save_qs(train_manifest, fs::path(DIR_OUT, "train_manifest.qs"))
save_qs(test_manifest,  fs::path(DIR_OUT, "test_manifest.qs"))

readr::write_csv(train_manifest, fs::path(DIR_OUT, "train_manifest.csv"))
readr::write_csv(test_manifest,  fs::path(DIR_OUT, "test_manifest.csv"))

# Build one global CpG intersection across all usable studies

COMMON_CPGS <- global_common_cpgs(proc_all)

message("Global shared CpGs across all processed studies: ", 
        length(COMMON_CPGS))

# Check that we have usable data

if (length(COMMON_CPGS) == 0) {
  stop(" x CRITICAL: No shared CpGs found across studies!\n\n",
       "This indicates a serious problem:\n",
       "  1. All studies failed to load methylation data, OR\n",
       "  2. No valid 27K/450K probes overlap, OR\n",
       "  3. Probe harmonization filtered everything out\n\n",
       "Check STEP 2 output above for which studies loaded successfully.\n",
       "Ensure GEOquery successfully downloaded phenotypes and IDAT files.")
}

message(" + Data validation passed")

proc_trn <- subset_proc_list_to_cpgs(proc_trn, COMMON_CPGS)
proc_tst <- subset_proc_list_to_cpgs(proc_tst, COMMON_CPGS)

beta_train <- combine_mats(proc_trn)
beta_test  <- combine_mats(proc_tst)

pheno_train <- combine_pheno(proc_trn)
pheno_test  <- combine_pheno(proc_tst)

# Order pheno exactly to matrix columns

pheno_train <- pheno_train[
  match(colnames(beta_train), pheno_train$sample_id), , drop = FALSE]

pheno_test  <- pheno_test[
  match(colnames(beta_test),   pheno_test$sample_id),  , drop = FALSE]

# Validate final alignment

if (!identical(colnames(beta_train), pheno_train$sample_id)) {
  stop(" x Training set alignment failed!\n",
       "  Matrix columns: ", ncol(beta_train), "\n",
       "  Phenotype rows: ", nrow(pheno_train), "\n",
       "  Sample IDs don't match between matrix and phenotype table.")
}

if (!identical(colnames(beta_test), pheno_test$sample_id)) {
  stop(" x Test set alignment failed!\n",
       "  Matrix columns: ", ncol(beta_test), "\n",
       "  Phenotype rows: ", nrow(pheno_test), "\n",
       "  Sample IDs don't match between matrix and phenotype table.")
}

if (!identical(rownames(beta_train), rownames(beta_test))) {
  stop(" x Feature (CpG) mismatch between train and test!\n",
       "  Train features: ", nrow(beta_train), "\n",
       "  Test features: ", nrow(beta_test), "\n",
       "  Both sets must have identical CpGs after harmonization.")
}

message(" + All alignment validations passed")

# Impute missing values using training probe medians only

imp <- impute_by_train_probe_median(beta_train, beta_test)

beta_train_imp <- imp$train
beta_test_imp  <- imp$test

# Final scorcher-ready objects: rows = samples, cols = CpGs

sc_train <- to_scorcher_xy(beta_train_imp, pheno_train)
sc_test  <- to_scorcher_xy(beta_test_imp, pheno_test)

scorcher_obj <- list(
  x_train = sc_train$x,
  y_train = sc_train$y,
  x_test = sc_test$x,
  y_test = sc_test$y,
  feature_names = sc_train$feature_names,
  train_sample_ids = sc_train$sample_ids,
  test_sample_ids = sc_test$sample_ids,
  pheno_train = pheno_train,
  pheno_test = pheno_test,
  beta_train_cpg_x_sample = beta_train_imp,
  beta_test_cpg_x_sample = beta_test_imp,
  train_probe_medians = imp$train_probe_medians,
  common_cpgs = COMMON_CPGS,
  train_manifest = train_manifest,
  test_manifest = test_manifest
)

final_summary <- tibble::tibble(
  dataset = c("train", "test"),
  n_samples = c(nrow(sc_train$x), nrow(sc_test$x)),
  n_cpg = c(ncol(sc_train$x), ncol(sc_test$x)),
  age_min = c(min(sc_train$y, na.rm = TRUE), min(sc_test$y, na.rm = TRUE)),
  age_median = c(stats::median(sc_train$y, na.rm = TRUE), 
                 stats::median(sc_test$y,  na.rm = TRUE)),
  age_max = c(max(sc_train$y, na.rm = TRUE), max(sc_test$y, na.rm = TRUE))
)

save_qs(pheno_train, fs::path(DIR_OUT, "pheno_train_final.qs"))
save_qs(pheno_test,  fs::path(DIR_OUT, "pheno_test_final.qs"))

readr::write_csv(pheno_train, fs::path(DIR_OUT, "pheno_train_final.csv"))
readr::write_csv(pheno_test,  fs::path(DIR_OUT, "pheno_test_final.csv"))

save_qs(beta_train_imp, fs::path(DIR_OUT, "beta_train_cpg_x_sample.qs"))
save_qs(beta_test_imp,  fs::path(DIR_OUT, "beta_test_cpg_x_sample.qs"))
save_qs(scorcher_obj,   fs::path(DIR_OUT, "scorcher_train_test.qs"))

save_qs(final_summary, fs::path(DIR_OUT, "final_summary.qs"))
readr::write_csv(final_summary, fs::path(DIR_OUT, "final_summary.csv"))

#=== SUMMARIZE DATA PROCESSING =================================================

message("\nOUTPUT FILES GENERATED:\n")

required_outputs <- c(
  "scorcher_train_test.qs",
  "pheno_train_final.qs",
  "pheno_test_final.qs",
  "beta_train_cpg_x_sample.qs",
  "beta_test_cpg_x_sample.qs",
  "pheno_train_final.csv",
  "pheno_test_final.csv",
  "final_summary.csv"
)

all_exist <- TRUE

for (output in required_outputs) {
  path <- fs::path(DIR_OUT, output)
  if (file.exists(path)) {
    size <- file.size(path) / 1024^2
    message("   + ", basename(path), " (", round(size, 1), " MB)")
  } else {
    message("   x MISSING: ", basename(path))
    all_exist <- FALSE
  }
}

if (!all_exist) {
  warning("\n x Some output files are missing. This may indicate an issue.")
}

# Data quality summary

message("\nDATA QUALITY SUMMARY:\n")
message("  Training samples:  ", nrow(sc_train$x))
message("  Test samples:      ", nrow(sc_test$x))
message("  Total samples:     ", nrow(sc_train$x) + nrow(sc_test$x))
message("  Features (CpGs):   ", ncol(sc_train$x))
message("  ")
message("  Age range (train): ", round(min(sc_train$y, na.rm=T), 1),
        " - ", round(max(sc_train$y, na.rm=T), 1), " years")
message("  Age range (test):  ", round(min(sc_test$y, na.rm=T), 1),
        " - ", round(max(sc_test$y, na.rm=T), 1), " years")
message("  ")
message("  Missing values imputed using training set medians")
message("\n")

message("NEXT STEPS:\n")
message("  1. Review outputs in: ", DIR_OUT, "/")
message("  2. Load main object: scorcher_obj ",
        "<- qs2::qs_read('data/processed/scorcher_train_test.qs')")
message("  3. Access training data: ",
        "x_train = scorcher_obj$x_train, y_train = scorcher_obj$y_train")
message("  4. Access test data: ", 
        "x_test = scorcher_obj$x_test, y_test = scorcher_obj$y_test")
message("\n")

message("Done.")

print(final_summary)

#=== END =======================================================================
sessionInfo()
