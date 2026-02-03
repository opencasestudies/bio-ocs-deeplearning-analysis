#===============================================================================
#
#  PROGRAM: data.R
#
#  AUTHOR:  Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: Download and preprocess publicly available GEO data to train a 
#           deep learning-based epigenetic clock. 
#
#  INPUT:   Select GEO DNA methylation studies (27K and 450K arrays):
#
#           Training Studies: 
#
#           GSE106648, GSE125105, GSE128235, GSE19711, GSE27044, GSE30870,
#           GSE40279,  GSE41037,  GSE52588,  GSE53740, GSE58119, GSE67530,  
#           GSE77445,  GSE77696,  GSE81961,  GSE84624, GSE97362
#
#           Testing Studies:  
#
#           GSE102177, GSE103911, GSE105123, GSE107459, GSE107737, GSE112696, 
#           GSE34639,  GSE37008,  GSE59065,  GSE61496,  GSE79329,  GSE87582,  
#           GSE87640,  GSE98876,  GSE99624
#
#  OUTPUT:  Preprocessed phenotype data frames and DNA methylation matrices:
#
#           - phenotype data frames (train/test)                                 ### FIX/FINISH
#           - methylation matrices aligned to a common CpG set                   ### FIX/FINISH
#
#  NOTES:
#
#    1. Age parsing:
#
#         - The reference paper's preprocessing effectively requires a per-GSE 
#           mapping from GEO fields (characteristics_ch1, title, etc.) to age 
#           (and units).
#
#         - `age_guess` is only a placeholder for now. Need to go through and 
#           make sure I have a clear/reproducible pipeline for doing this. 
#
#    2. Probe filtering:
#
#         - The reference paper removes sex chromosomes (done here) AND 
#           multi-mappers/cross-reactives.
#
#         - Need to implement a published cross-reactive probe removal list and 
#           cite it.
#
#    3. Platform normalization:
#
#         - Using `preprocessNoob` as a modern default. If needed, I can switch
#           back to `lumi` as in the original paper for the preprocessing step.
#
#    4. Fallback matrices:
#
#         - series_matrix_exprs_fallback may not be Beta values (e.g., M-values,
#           normalized intensities, etc.)
#
#         - Part of this preprocessing code creates a manifest that records 
#           which method was used so I can go back and manually check each 
#           study.
#
#  UPDATED: 2026-02-03
#
#===============================================================================

#=== SETUP =====================================================================

options(timeout = 1000)

#--- PACKAGE MANAGEMENT --------------------------------------------------------

#-- CRAN Packages

# NOTE: Using `pacman` for CRAN package management

if (!requireNamespace("pacman", quietly=TRUE)) {
  
  install.packages("pacman")
}

library(pacman)

p_load(data.table, qs2, fs, Matrix, tidyverse, update=FALSE)                     ### TODO: Use the `here` package 
                                                                                 ### throughout for relative filepaths?
#-- BioC Packages

# NOTE: Using `BiocManager` for Bioconductor package management

if (!requireNamespace("BiocManager", quietly=TRUE)) {
  
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  
  "GEOquery",
  "minfi",
  "limma",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylation27kmanifest",
  "IlluminaHumanMethylation27kanno.ilmn12.hg19"
)

missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, 
                                  
  requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_bioc) > 0) {
  
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

# Load BioC Libraries Explicitly

library(GEOquery)
library(Biobase)
library(minfi)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation27kanno.ilmn12.hg19)

# NOTE: The original paper used `lumi` for color correction and normalization, 
# but minfi seems to be better supported. I will see how this looks and if not, 
# I can go back to `lumi`-based processing.

#--- PATHS ---------------------------------------------------------------------

DIR_GEO_PHENO <- "data/geoquery"
DIR_RAW_IDAT  <- "data/raw/idat"
DIR_INTER     <- "data/intermediate"
DIR_OUT       <- "data/processed"

dir_create(DIR_GEO_PHENO, recurse = TRUE)
dir_create(DIR_RAW_IDAT,  recurse = TRUE)
dir_create(DIR_INTER,     recurse = TRUE)
dir_create(DIR_OUT,       recurse = TRUE)

#--- HELPER FUNCTIONS ----------------------------------------------------------

#-- Fetch GEO Phenotype (Series Matrix) as a data.frame

fetch_gse_pheno <- function(gse_id, dest = "data/geoquery") {
  
  dir.create(dest, recursive=TRUE, showWarnings=FALSE)
  
  gse <- GEOquery::getGEO(gse_id, GSEMatrix=TRUE, getGPL=FALSE, destdir=dest)
  
  if (is.list(gse)) gse <- gse[[1]]
  
  ph <- Biobase::pData(gse) |>
    
    as.data.frame()
  
  ph$gse <- gse_id
  
  ph
}

#-- Download GEO Supplementary Files

download_geo_raw <- function(gse_id, dest = "data/raw/idat") {
  
  gse_dir <- fs::path(dest, gse_id)
  
  dir.create(gse_dir, recursive=TRUE, showWarnings=FALSE)
  
  message("Downloading supplementary files for ", gse_id)
  
  GEOquery::getGEOSuppFiles(gse_id, makeDirectory = FALSE, baseDir = gse_dir)
  
  archives <- fs::dir_ls(gse_dir, regexp = "\\.(tar|tar.gz|tgz)$")
  
  for (a in archives) {
    
    message("Extracting ", fs::path_file(a))
    
    utils::untar(a, exdir = gse_dir)
  }
  
  invisible(gse_dir)
}

#-- Get Potential Age Columns from Study

extract_age_guess <- function(ph) {
  
  cn <- names(ph)
  
  idx <- which(stringr::str_detect(tolower(cn), "age"))
  
  if (!length(idx)) return(rep(NA_real_, nrow(ph)))
  
  col <- cn[idx[1]]
  
  suppressWarnings(as.numeric(stringr::str_extract(ph[[col]], "\\d+\\.?\\d*")))
}

#-- Add 0.5 to Integer-Only Ages (Matches Paper's Convention)

add_half_if_integer <- function(age) {
  
  is_int <- is.finite(age) & abs(age - round(age)) < 1e-8
  
  age2 <- age
  
  age2[is_int] <- age2[is_int] + 0.5
  
  age2
}

#-- Identify Common CpGs Shared between 27K and 450K, Drop Chr X/Y

# NOTE: The paper also removes multi-mappers / cross-reactive probes.

# TODO: Apply a published probe-filter list (e.g., Chen et al., Pidsley et al.)
# and document it in the case study.

common_cpgs_27k_450k <- function(drop_sex = TRUE) {
  
  ann450 <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
  ann27  <- minfi::getAnnotation(IlluminaHumanMethylation27kanno.ilmn12.hg19)
  
  cpgs <- intersect(rownames(ann450), rownames(ann27))
  
  if (drop_sex) {
    
    chr450 <- ann450[cpgs, "chr"]
    
    cpgs <- cpgs[!(chr450 %in% c("chrX", "chrY"))]
  }

  cpgs
}

#-- Process a Single GSE into a CpG x Sample Matrix

# NOTE: The preference order for raw data files is:
#   1. raw IDATs -> preprocessNoob -> Beta
#   2. fallback: Series Matrix exprs (Need to document! May not be Beta)

process_one_gse <- function(gse_id) {
  
  gse_dir <- fs::path(DIR_RAW_IDAT, gse_id)
  
  # Try IDAT
  
  beta <- read_beta_from_idats(gse_dir)
  
  method <- "idat_minfi_noob_beta"
  
  if (is.null(beta)) {
    
    # Fallback
    
    beta <- read_expr_from_series_matrix(gse_id)
    
    method <- "series_matrix_exprs_fallback"
  }
  
  if (is.null(beta)) {
    
    warning("No methylation matrix could be created for ", gse_id)
    
    return(list(gse = gse_id, method = NA_character_, mat = NULL))
  }
  
  # Restrict to Common CpGs
  
  beta2 <- restrict_to_common_cpgs(beta, COMMON_CPGS)
  
  list(gse = gse_id, method = method, mat = beta2)
}

#-- Find IDAT Files for a GSE

find_idats <- function(gse_dir) {
  
  fs::dir_ls(gse_dir, recurse = TRUE, 
             
    regexp = "\\.idat(\\.gz)?$", type = "file")
}

#-- Try to Construct a minfi "targets" Sheet from standard IDAT Naming

#  Naming: <SentrixID>_<SentrixPosition>_Red.idat / _Grn.idat
#  Returns a data.frame with columns required by read.metharray.exp:
#   - Basename: Full path without _Red/_Grn suffix

make_targets_from_idats <- function(idat_files) {
  
  # Normalize gz (minfi can read gz in many cases; if trouble, gunzip upstream)
  
  idat_files <- as.character(idat_files)
  
  # Strip _Red.idat, _Grn.idat, plus optional .gz
  
  base <- idat_files |>
    
    stringr::str_replace("(_Red|_Grn)\\.idat(\\.gz)?$", "")
  
  # Keep only basenames that have both channels
  
  tab <- tibble::tibble(file = idat_files, base = base) |>
    
    dplyr::group_by(base) |>
    
    dplyr::summarise(
      
      has_red = any(stringr::str_detect(file, "_Red\\.idat")),
      has_grn = any(stringr::str_detect(file, "_Grn\\.idat")),
      
      .groups = "drop") |>
    
    dplyr::filter(has_red & has_grn)
  
  if (nrow(tab) == 0) return(NULL)
  
  targets <- tibble::tibble(Basename = tab$base)
  
  # Extract Sentrix Identifiers if Possible 
  # Base Typically Ends with ".../<SentrixID>_<SentrixPosition>"
  
  fname <- fs::path_file(targets$Basename)
  
  parts <- stringr::str_split(fname, "_", simplify = TRUE)
  
  if (ncol(parts) >= 2) {
    
    targets$Sentrix_ID <- parts[, 1]
    
    targets$Sentrix_Position <- parts[, 2]
  }
  
  as.data.frame(targets)
}

#-- Read + preprocess IDATs into Beta values (minfi)

# Uses preprocessNoob (default)
# Returns a Matrix: CpGs x Samples (beta values in [0,1])

read_beta_from_idats <- function(gse_dir) {
  
  idats <- find_idats(gse_dir)
  
  if (length(idats) == 0) return(NULL)
  
  targets <- make_targets_from_idats(idats)
  
  if (is.null(targets) || nrow(targets) == 0) return(NULL)
  
  # minfi Reads by Directory (uses Basename)
  # Set the Working Directory Approach by Passing Targets Explicitly
  # read.metharray.exp will locate IDATs via Basename
  
  rg <- tryCatch(
    
    minfi::read.metharray.exp(targets = targets, verbose = TRUE),
    
    error = function(e) {
      
      message("IDAT read failed in ", gse_dir, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(rg)) return(NULL)
  
  mset <- tryCatch(
    
    minfi::preprocessNoob(rg),
    
    error = function(e) {
      message("preprocessNoob failed in ", gse_dir, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(mset)) return(NULL)
  
  beta <- minfi::getBeta(mset)
  
  beta
}

#-- Grab a Methylation-Like Matrix from a GEO Series Matrix

# NOTE: This is a fallback when raw IDATs are not available.

# NOTE: This may be M-values, Beta, or something else depending on submitter.
# Need to document which GSEs are using this fallback in the case study.

read_expr_from_series_matrix <- function(gse_id, dest = DIR_GEO_PHENO) {
  
  gse <- GEOquery::getGEO(gse_id, 
                          
    GSEMatrix = TRUE, getGPL = FALSE, destdir = dest)
  
  if (is.list(gse)) gse <- gse[[1]]
  
  expr <- Biobase::exprs(gse)
  
  if (is.null(expr) || nrow(expr) == 0) return(NULL)
  
  storage.mode(expr) <- "numeric"
  
  expr
}

#-- Restrict to Common CpGs and Align Column Names / Sample IDs

restrict_to_common_cpgs <- function(mat, common_cpgs) {
  
  keep <- intersect(rownames(mat), common_cpgs)
  
  mat[keep, , drop = FALSE]
}

#-- Save a Methylation Matrix in a Compact Format

save_matrix <- function(mat, out_path) {
  
  # NOTE: qs is fast, and matrices can be large but manageable for these GSEs.
  
  qs2::qs_save(mat, out_path)
}

#-- Combine study-level preprocessed data into training and testing sets

combine_mats <- function(mat_list) {
  
  mats <- purrr::keep(mat_list, ~ !is.null(.x$mat))
  
  if (length(mats) == 0) return(NULL)
  
  # Ensure same CpG ordering (use intersection across all processed mats)
  
  common_rows <- Reduce(intersect, lapply(mats, function(x) rownames(x$mat)))
  
  common_rows <- sort(common_rows)
  
  mats2 <- lapply(mats, function(x) {
    
    m <- x$mat[common_rows, , drop = FALSE]
    
    colnames(m) <- paste0(x$gse, ":", colnames(m))
    
    m
  })
  
  do.call(cbind, mats2)
}

#=== QUERY GEO DATA ============================================================

#--- 0. DEFINE STUDIES ---------------------------------------------------------

# NOTE: These are just the studies used in the original paper. We can augment 
# them with other studies for potential side quests/homework assignments.

#-- Training Set

trn_gses <- c(
  
  "GSE106648", "GSE125105", "GSE128235", "GSE19711", "GSE27044", "GSE30870",
  "GSE40279",  "GSE41037",  "GSE52588",  "GSE53740", "GSE58119", "GSE67530",
  "GSE77445",  "GSE77696",  "GSE81961",  "GSE84624", "GSE97362"
)

#-- Testing Set

tst_gses <- c(
  
  "GSE102177", "GSE103911", "GSE105123", "GSE107459", "GSE107737", "GSE112696",
  "GSE34639",  "GSE37008",  "GSE59065",  "GSE61496",  "GSE79329",  "GSE87582",
  "GSE87640",  "GSE98876",  "GSE99624"
)

#--- 1. FETCH + SAVE PHENOTYPES ------------------------------------------------

message("Fetching phenotype data (Series Matrix) ...")

ph_trn <- dplyr::bind_rows(lapply(trn_gses, fetch_gse_pheno))
ph_tst <- dplyr::bind_rows(lapply(tst_gses, fetch_gse_pheno))

ph_trn$age_guess <- add_half_if_integer(extract_age_guess(ph_trn))
ph_tst$age_guess <- add_half_if_integer(extract_age_guess(ph_tst))

qs2::qs_save(ph_trn, fs::path(DIR_INTER, "ph_trn.qs"))
qs2::qs_save(ph_tst, fs::path(DIR_INTER, "ph_tst.qs"))

#--- 2. DOWNLOAD SUPPLEMENTARY FILES -------------------------------------------

message("Downloading supplementary (raw) files for training GSEs ...")

invisible(sapply(trn_gses, download_geo_raw))

message("Downloading supplementary (raw) files for testing GSEs ...")

invisible(sapply(tst_gses, download_geo_raw))

#--- 3. BUILD DNA METHYLATION MATRICES -----------------------------------------

# NOTE: The paper preprocesses the methylation to a common CpG set (24,538) and 
#       used only CpGs shared between 450K and 27K and removed sex chromosomes
#       and multi-mappers.

# NOTE: Orthologous sequences on multiple chromosomes removal can be 
# implemented more strictly using published probe filtering lists (cross-
# reactive probes). I need to implement this and cite it in the case study.

#-- Common CpGs Used Across Platforms

COMMON_CPGS <- common_cpgs_27k_450k(drop_sex = TRUE)

message("Common CpGs between 27K and 450K, dropping chr X/Y): ", 
        
  length(COMMON_CPGS))

#-- TRAINING DATA

message("Processing training methylation matrices ...")

trn_list <- lapply(trn_gses, process_one_gse)

#- Save Per-GSE Outputs + A Processing Manifest

trn_manifest <- tibble::tibble(
  
  gse    = vapply(trn_list, `[[`, character(1), "gse"),
  
  method = vapply(trn_list, `[[`, character(1), "method"),
  
  n_cpg  = vapply(trn_list, 
                  
    function(x) if (is.null(x$mat)) NA_integer_ else nrow(x$mat), integer(1)),
  
  n_samp = vapply(trn_list, 
                  
    function(x) if (is.null(x$mat)) NA_integer_ else ncol(x$mat), integer(1))
)

qs2::qs_save(trn_manifest, fs::path(DIR_OUT, "train_manifest.qs"))

for (x in trn_list) {
  
  if (!is.null(x$mat)) {
    
    save_matrix(x$mat, fs::path(DIR_OUT, paste0(x$gse, "_train_mat.qs")))
  }
}

#-- TESTING DATA

message("Processing testing methylation matrices ...")

tst_list <- lapply(tst_gses, process_one_gse)

tst_manifest <- tibble::tibble(
  
  gse    = vapply(tst_list, `[[`, character(1), "gse"),
  
  method = vapply(tst_list, `[[`, character(1), "method"),
  
  n_cpg  = vapply(tst_list, 
                  
    function(x) if (is.null(x$mat)) NA_integer_ else nrow(x$mat), integer(1)),
  
  n_samp = vapply(tst_list, 
                  
    function(x) if (is.null(x$mat)) NA_integer_ else ncol(x$mat), integer(1))
)

qs2::qs_save(tst_manifest, fs::path(DIR_OUT, "test_manifest.qs"))

for (x in tst_list) {
  
  if (!is.null(x$mat)) {
    
    save_matrix(x$mat, fs::path(DIR_OUT, paste0(x$gse, "_test_mat.qs")))
  }
}

#--- 4. COMBINE INTO SINGLE TRAIN/TEST MATRICES --------------------------------

# NOTE: This requires column-binding across GSEs after restricting CpGs.

# NOTE: Sample IDs can collide across studies; prefixing with GSE for safety.

train_big <- combine_mats(trn_list)

test_big  <- combine_mats(tst_list)

if (!is.null(train_big)) {
  
  save_matrix(train_big, fs::path(DIR_OUT, "train_big.qs"))
} 

if (!is.null(test_big)) {
  
  save_matrix(test_big,  fs::path(DIR_OUT, "test_big.qs"))
}

#=== END =======================================================================

sessionInfo()
