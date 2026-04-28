#===============================================================================
#
#  PROGRAM: data_helpers.R
#
#  AUTHOR:  Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: Define helper functions used by the GEO methylation preprocessing
#           pipeline for the epigenetic clock case study:
#
#             - GEO phenotype retrieval
#             - Supplementary file download and IDAT discovery
#             - Raw IDAT preprocessing with sesame
#             - Fallback series matrix handling
#             - Probe filtering and CpG harmonization
#             - Sample ID alignment between phenotype and methylation data
#             - Age parsing and curated study-specific age extraction
#             - Missing data imputation and scorcher-ready output formatting
#
#  INPUT:   This file does not execute the preprocessing pipeline on its own.
#           Instead, it provides reusable functions that are called by the main
#           data preprocessing script. Expected external inputs used by some
#           helper functions include:
#
#             - GEO study accession IDs (e.g., GSE40279)
#             - Downloaded GEO series matrix files
#             - Downloaded supplementary raw files (e.g., IDATs)
#             - Optional reference files such as:
#
#                 * cross_reactive_probes.csv
#                 * age_rules.csv
#
#  OUTPUT:  This file does not directly write final processed datasets. 
#           Its functions return intermediate objects used by the main 
#           preprocessing workflow, including:
#
#             - Phenotype data frames
#             - Methylation beta-value matrices
#             - Shared CpG vectors
#             - Aligned sample metadata
#             - Parsed age vectors
#             - Scorcher-ready x/y objects
#
#  NOTES:   Please note the following:
#
#             1. Separation of Concerns:
#
#                 - This file is intended to keep reusable helper logic 
#                   separate from the main execution script so that the 
#                   primary workflow in data.R remains easier to read, 
#                   debug, and teach from.
#
#             2. Reproducibility:
#
#                 - Some helper functions implement pragmatic fallback 
#                   behavior (e.g., age parsing, series matrix processing).
#
#                 - For publication-quality reproduction, study-specific age 
#                   rules and a documented cross-reactive probe list should 
#                   be supplied and reviewed explicitly.
#
#             3. Platform harmonization:
#
#                 - The intended preprocessing uses CpGs shared across the 
#                   Illumina 27K and 450K arrays, excluding sex chromosome 
#                   probes and probes with problematic mappings, so these 
#                   helpers include functionality for building and applying 
#                   a fixed common CpG panel.
#
#             4. Scope:
#
#                   - These functions are written for the epigenetic clock 
#                     case study, but many are general enough to be reused 
#                     in related public methylation preprocessing pipelines.
#
#  UPDATED: 2026-03-24
#
#===============================================================================

#=== HELPER FUNCTIONS ==========================================================

#--- NULL-COALESCING HELPER ----------------------------------------------------

# PURPOSE:
#   Return `y` when `x` is NULL or length-0; otherwise return `x`.
#
# INPUT:
#   x - Any R object.
#   y - Fallback value to use when `x` is NULL or empty.
#
# OUTPUT:
#   Either `x` or `y`.

`%||%` <- function(x, y) {
  
  if (is.null(x) || length(x) == 0) y else x
}

#--- NORMALIZE SAMPLE IDENTIFIERS ----------------------------------------------

# PURPOSE:
#   Standardize sample identifiers so phenotype records and methylation matrix
#   columns can be matched more reliably across GEO metadata sources.
#
# INPUT:
#   x - Character vector of sample names, file paths, GSM IDs, or IDAT prefixes.
#
# OUTPUT:
#   Character vector of normalized sample keys.

sample_key_normalize <- function(x) {
  
  x |>
    as.character() |>
    basename() |>
    stringr::str_replace("\\.gz$", "") |>
    stringr::str_replace("\\.idat$", "") |>
    stringr::str_replace("(_Red|_Grn)$", "") |>
    stringr::str_replace("\\s+", "_") |>
    stringr::str_replace_all("[^A-Za-z0-9._-]", "_") |>
    toupper()
}

#--- CLIP BETA VALUES TO [0,1] -------------------------------------------------

clip_beta <- function(x) {
  
  x[x < 0] <- 0
  x[x > 1] <- 1
  x
}

#--- INTEGER-ISH CHECK ---------------------------------------------------------

is_integerish <- function(x, tol = 1e-8) {
  
  is.finite(x) & abs(x - round(x)) < tol
}

#--- ADD 0.5 TO INTEGER AGES ---------------------------------------------------

add_half_if_integer <- function(age) {

  out <- age
  idx <- is_integerish(out)
  out[idx] <- out[idx] + 0.5
  out
}

#--- COALESCE NUMERIC VALUES ---------------------------------------------------

coalesce_numeric <- function(...) {
  xs <- list(...)
  if (length(xs) == 0) return(numeric())
  out <- xs[[1]]
  if (length(xs) > 1) {
    for (j in 2:length(xs)) {
      idx <- is.na(out) & !is.na(xs[[j]])
      out[idx] <- xs[[j]][idx]
    }
  }
  out
}

#--- FETCH A GEO SERIES OBJECT (WITH RETRY LOGIC) -----

# PURPOSE:
#   Fetch GEO object with robust error handling and retry logic.
#   Network calls to GEO servers can be unreliable, so we retry on failure.
#
# INPUT:
#   gse_id - GEO Series accession (e.g., "GSE40279")
#   dest   - Directory to cache downloaded files
#   max_retries - Number of times to retry on network failure
#
# OUTPUT:
#   GEO ExpressionSet object, or NULL on failure
#
# NOTES:
#   - Retries up to 3 times on network/API errors
#   - Returns NULL on final failure (graceful degradation)
#   - Logs which attempts succeeded/failed

fetch_gse_object <- function(gse_id, dest = DIR_GEO_PHENO, max_retries = 3) {

  for (attempt in 1:max_retries) {
    tryCatch({
      gse <- GEOquery::getGEO(
        gse_id,
        GSEMatrix = TRUE,
        getGPL = FALSE,
        destdir = dest
      )

      if (is.list(gse)) gse <- gse[[1]]

      return(gse)

    }, error = function(e) {
      error_msg <- conditionMessage(e)

      if (attempt < max_retries) {
        message("    [Retry ", attempt, "/", max_retries, "] Failed to fetch ", gse_id, ": ",
                substr(error_msg, 1, 80), "...")
      } else {
        message("    [FAILED] Could not fetch ", gse_id, " after ", max_retries, " attempts:")
        message("    Error: ", error_msg)
      }
    }, warning = function(w) NULL)
  }

  # All retries failed
  return(NULL)
}

#--- FETCH GEO PHENOTYPE DATA (WITH ERROR HANDLING) ---

# PURPOSE:
#   Retrieve raw phenotype metadata for a GEO study with error handling.
#   If GEO fetch fails, returns an empty but valid data frame.
#
# INPUT:
#   gse_id - GEO Series accession (e.g., "GSE40279")
#   dest   - Directory to cache downloaded GEO objects
#
# OUTPUT:
#   Data frame with phenotype variables as columns, where each row is a sample.
#   Always includes geo_accession (GSM IDs) and gse (study ID) columns.
#   Returns empty data frame (0 rows) if fetch fails.
#
# NOTES:
#   - Uses fetch_gse_object() which includes retry logic
#   - Returns consistent data frame shape even on failure
#   - Caller should check nrow() to detect failures

fetch_gse_pheno <- function(gse_id, dest = DIR_GEO_PHENO) {

  gse <- fetch_gse_object(gse_id, dest = dest)

  if (is.null(gse)) {
    # Return empty data frame with required columns for consistency
    return(data.frame(
      geo_accession = character(),
      gse = character(),
      stringsAsFactors = FALSE
    ))
  }

  tryCatch({
    ph <- Biobase::pData(gse) |>
      as.data.frame()

    ph$geo_accession <- rownames(ph)
    ph$gse <- gse_id

    ph

  }, error = function(e) {
    message("    [ERROR] Could not extract phenotype for ", gse_id, ": ",
            conditionMessage(e))
    # Return empty frame on error
    return(data.frame(
      geo_accession = character(),
      gse = character(),
      stringsAsFactors = FALSE
    ))
  })
}

#--- DOWNLOAD GEO SUPPLEMENTARY FILES ------------------------------------------
#
# PURPOSE:
#   Download supplementary files for a GEO study.
#   This is STEP 1: Downloads are cached as-is (no extraction).
#
# INPUT:
#   gse_id - GEO Series accession (e.g., "GSE40279")
#   dest   - Directory to store downloaded files (defaults to DIR_RAW_DOWNLOADS)
#
# OUTPUT:
#   Invisibly returns gse_dir path.
#
# NOTES:
#   - Creates gse_dir if it doesn't exist
#   - Downloads are cached; safe to call multiple times
#   - Archives are NOT extracted here; all zipped files stay intact
#   - STEP 2: Call extract_geo_supplemental() to extract archives
#   - STEP 3: Call extract_nested_tar_files() for secondary _RAW.tar files
#   - STEP 4: Call extract_idat_gz_files() for methylation IDAT data

download_geo_raw <- function(gse_id, dest = DIR_RAW_DOWNLOADS) {

  gse_dir <- fs::path(dest, gse_id)
  fs::dir_create(gse_dir, recurse = TRUE)

  message("  [STEP 1/4] Downloading supplementary files for ", gse_id)

  tryCatch({
    GEOquery::getGEOSuppFiles(
      GEO = gse_id,
      makeDirectory = FALSE,
      baseDir = gse_dir
    )
    message("    + Downloaded to: ", gse_dir)
  }, error = function(e) {
    message("    [!] Download failed: ", conditionMessage(e))
  })

  invisible(gse_dir)
}

#--- ExTRACT GEO SUPPLEMENTARY FILES -------------------------------------------
#
# PURPOSE:
#   Recursively extract tar archives in a GEO study directory.
#   Tracks which files have been extracted to prevent infinite loops.
#
# INPUT:
#   gse_id - GEO Series accession (e.g., "GSE40279")
#   dest   - Directory containing downloaded files
#
# OUTPUT:
#   Invisibly returns gse_dir path.
#
# NOTES:
#   - Extracts all tar/tar.gz/tgz files recursively
#   - Tracks extracted files to prevent re-extraction of the same archive
#   - Keeps new tar files that appear after extraction to handle nested archives
#   - Safe from infinite loops: each original tar is only processed once

# extract_geo_raw <- function(gse_id, dest = DIR_RAW_IDAT) {
# 
#   gse_dir <- fs::path(dest, gse_id)
#   fs::dir_create(gse_dir, recurse = TRUE)
# 
#   message("Extracting supplementary files for ", gse_id)
# 
#   # Track archives we've already processed to prevent infinite loops
#   extracted_set <- character()
# 
#   repeat {
#     archives <- fs::dir_ls(
#       gse_dir,
#       recurse = TRUE,
#       regexp = "\\.(tar|tar\\.gz|tgz)$",
#       type = "file"
#     )
# 
#     # Convert to character for comparison
#     archives <- as.character(archives)
# 
#     # Find archives not yet extracted
#     new_archives <- setdiff(archives, extracted_set)
# 
#     if (length(new_archives) == 0) break
# 
#     for (a in new_archives) {
#       message("  Extracting ", fs::path_file(a))
#       try(utils::untar(a, exdir = gse_dir, extras = "-o"), silent = TRUE)
#       extracted_set <- c(extracted_set, a)
#     }
#   }
# 
#   invisible(gse_dir)
# }

extract_geo_raw <- function(gse_id, src = DIR_RAW_DOWNLOADS, dest = DIR_RAW_ExTRACTED) {

  #--- PURPOSE -----------------------------------------------------------------
  # STEP 2: Extract top-level archives from downloaded files.
  #
  # This function extracts tar/tar.gz/tgz files from the download directory.
  # Extracted files go to a SEPARATE directory to keep original archives intact.
  #
  # INPUT:
  #   gse_id - GEO Series accession
  #   src    - Source directory where archives were downloaded (DIR_RAW_DOWNLOADS)
  #   dest   - Destination directory for extracted files (DIR_RAW_ExTRACTED)
  #
  # OUTPUT:
  #   Invisibly returns dest_dir path.
  #
  # NOTES:
  #   - Keeps original archives safe in src directory
  #   - Extracts all tar/tar.gz/tgz files recursively
  #   - Tracks extracted files to prevent infinite loops
  #   - After this, call extract_nested_tar_files() for _RAW.tar files
  #   - Then call extract_idat_gz_files() for .idat.gz files

  src_dir <- fs::path(src, gse_id)
  dest_dir <- fs::path(dest, gse_id)

  fs::dir_create(dest_dir, recurse = TRUE)

  message("  [STEP 2/4] Extracting supplemental archives for ", gse_id)

  # Check if source directory has any archives
  archives <- fs::dir_ls(
    src_dir,
    recurse = TRUE,
    regexp = "\\.(tar|tar\\.gz|tgz)$",
    type = "file"
  ) |> as.character()

  if (length(archives) == 0) {
    message("    x No archives found in download directory")
    return(invisible(dest_dir))
  }

  extracted_set <- character()

  repeat {
    current_archives <- fs::dir_ls(
      src_dir,
      recurse = TRUE,
      regexp = "\\.(tar|tar\\.gz|tgz)$",
      type = "file"
    ) |> as.character()

    new_archives <- setdiff(current_archives, extracted_set)

    if (length(new_archives) == 0) break

    for (a in new_archives) {
      archive_name <- fs::path_file(a)
      message("    Extracting: ", archive_name)

      ok <- FALSE

      # Try system tar with overwrite (-o flag)
      res <- try(
        utils::untar(a, exdir = dest_dir, extras = "-o"),
        silent = TRUE
      )

      if (!inherits(res, "try-error")) {
        ok <- TRUE
        message("      + Extracted successfully")
      } else {
        # Fallback to internal tar
        res2 <- try(
          utils::untar(a, exdir = dest_dir, tar = "internal"),
          silent = TRUE
        )
        if (!inherits(res2, "try-error")) {
          ok <- TRUE
          message("      + Extracted (internal tar)")
        }
      }

      if (!ok) {
        warning("Failed to extract archive: ", a, call. = FALSE)
      }

      extracted_set <- c(extracted_set, a)
    }
  }

  message("    + Extraction complete for ", gse_id)
  invisible(dest_dir)
}


#--- ExTRACT NESTED TAR FILES ---------------------------------------------------

extract_nested_tar_files <- function(gse_id, src = DIR_RAW_ExTRACTED, dest = DIR_RAW_TAR) {

  #--- PURPOSE -----------------------------------------------------------------
  # STEP 3: Extract secondary/nested TAR archives (e.g., _RAW.tar).
  #
  # Some GEO studies have nested tar files that need secondary extraction.
  # These typically contain raw platform files (_RAW.tar) or additional archives.
  #
  # INPUT:
  #   gse_id - GEO Series accession
  #   src    - Source directory with already-extracted files (DIR_RAW_ExTRACTED)
  #   dest   - Destination directory for nested tar extraction (DIR_RAW_TAR)
  #
  # OUTPUT:
  #   Invisibly returns dest_dir path.
  #
  # NOTES:
  #   - Only processes tar/tar.gz/tgz files found as files (not archives)
  #   - Keeps source extraction directory intact
  #   - After this, you can look for .idat.gz files in either:
  #     * dest directory (from nested extraction)
  #     * src directory (if .idat.gz was in top-level files)

  src_dir <- fs::path(src, gse_id)
  dest_dir <- fs::path(dest, gse_id)

  fs::dir_create(dest_dir, recurse = TRUE)

  message("  [STEP 3/4] Extracting nested TAR files for ", gse_id)

  # Look for tar files in extracted directory
  nested_tars <- fs::dir_ls(
    src_dir,
    recurse = TRUE,
    regexp = "\\.(tar|tar\\.gz|tgz)$",
    type = "file"
  ) |> as.character()

  if (length(nested_tars) == 0) {
    message("    x No nested TAR files found")
    return(invisible(dest_dir))
  }

  extracted_set <- character()

  repeat {
    current_tars <- fs::dir_ls(
      src_dir,
      recurse = TRUE,
      regexp = "\\.(tar|tar\\.gz|tgz)$",
      type = "file"
    ) |> as.character()

    new_tars <- setdiff(current_tars, extracted_set)

    if (length(new_tars) == 0) break

    for (tar_file in new_tars) {
      tar_name <- fs::path_file(tar_file)
      message("    Extracting: ", tar_name)

      ok <- FALSE

      # Try system tar
      res <- try(
        utils::untar(tar_file, exdir = dest_dir, extras = "-o"),
        silent = TRUE
      )

      if (!inherits(res, "try-error")) {
        ok <- TRUE
        message("      + Extracted successfully")
      } else {
        # Fallback to internal tar
        res2 <- try(
          utils::untar(tar_file, exdir = dest_dir, tar = "internal"),
          silent = TRUE
        )
        if (!inherits(res2, "try-error")) {
          ok <- TRUE
          message("      + Extracted (internal tar)")
        }
      }

      if (!ok) {
        warning("Failed to extract nested tar: ", tar_file, call. = FALSE)
      }

      extracted_set <- c(extracted_set, tar_file)
    }
  }

  message("    + Nested TAR extraction complete")
  invisible(dest_dir)
}

#--- ExTRACT .IDAT.GZ FILES ---------------------------------------------------

extract_idat_gz_files <- function(gse_id,
                                  src_dirs = c(DIR_RAW_ExTRACTED, DIR_RAW_TAR),
                                  dest = DIR_RAW_IDAT) {

  #--- PURPOSE -----------------------------------------------------------------
  # STEP 4: Extract .idat.gz files to final IDAT directory.
  #
  # This function finds all .idat.gz files (compressed IDAT methylation data)
  # in the source directories and decompresses them to the final IDAT directory.
  #
  # INPUT:
  #   gse_id   - GEO Series accession
  #   src_dirs - Vector of directories to search for .idat.gz files
  #              (default: both extracted and tar-extracted directories)
  #   dest     - Destination directory for decompressed IDAT files (DIR_RAW_IDAT)
  #
  # OUTPUT:
  #   Invisibly returns dest_dir path.
  #
  # NOTES:
  #   - Searches recursively through all src_dirs
  #   - Only decompresses if uncompressed file doesn't already exist
  #   - Uses R.utils::gunzip() for decompression
  #   - Skips files that are already decompressed (.idat without .gz)

  dest_dir <- fs::path(dest, gse_id)
  fs::dir_create(dest_dir, recurse = TRUE)

  message("  [STEP 4/4] Extracting .idat.gz files for ", gse_id)

  gz_idats <- character()

  # Search all source directories
  for (src_dir in src_dirs) {
    full_src <- fs::path(src_dir, gse_id)

    if (!dir.exists(full_src)) next

    found <- tryCatch({
      fs::dir_ls(
        full_src,
        recurse = TRUE,
        regexp = "\\.idat\\.gz$",
        type = "file"
      ) |> as.character()
    }, error = function(e) character())

    gz_idats <- c(gz_idats, found)
  }

  if (length(gz_idats) == 0) {
    message("    x No .idat.gz files found")
    return(invisible(dest_dir))
  }

  gz_idats <- unique(gz_idats)
  message("    Found ", length(gz_idats), " .idat.gz file(s)")

  success_count <- 0
  fail_count <- 0

  for (gz_file in gz_idats) {
    uncompressed_file <- stringr::str_replace(gz_file, "\\.gz$", "")
    uncompressed_name <- fs::path_file(uncompressed_file)

    # Copy to destination and decompress in one step
    dest_gz <- fs::path(dest_dir, basename(gz_file))
    dest_uncompressed <- stringr::str_replace(dest_gz, "\\.gz$", "")

    tryCatch({
      # Only decompress if uncompressed file doesn't already exist
      if (!file.exists(dest_uncompressed)) {
        R.utils::gunzip(gz_file, destname = dest_uncompressed, remove = FALSE)
        message("      + ", uncompressed_name)
        success_count <- success_count + 1
      } else {
        message("      ~ ", uncompressed_name, " (already exists)")
      }
    }, error = function(e) {
      message("      [!] Failed to decompress ", basename(gz_file), ": ", conditionMessage(e))
      fail_count <<- fail_count + 1
    })
  }

  message("    + Decompressed ", success_count, " file(s)",
          if (fail_count > 0) paste0(" (", fail_count, " failed)") else "")

  invisible(dest_dir)
}

#--- LOCATE IDAT FILES ---------------------------------------------------------

find_idats <- function(gse_dir) {

  fs::dir_ls(
    gse_dir,
    recurse = TRUE,
    regexp = "\\.(idat|idat\\.gz)$",
    type = "file"
  )
}

#--- IDENTIFY COMPLETE IDAT PREFIxES -------------------------------------------

find_idat_prefixes <- function(gse_dir) {
  
  idats <- find_idats(gse_dir)
  if (length(idats) == 0) return(character())

  x <- as.character(idats)

  prefixes <- unique(stringr::str_replace(
    x, "(_Red|_Grn)\\.idat(\\.gz)?$", ""
  ))

  has_red <- file.exists(paste0(prefixes, "_Red.idat")) |
    file.exists(paste0(prefixes, "_Red.idat.gz"))

  has_grn <- file.exists(paste0(prefixes, "_Grn.idat")) |
    file.exists(paste0(prefixes, "_Grn.idat.gz"))

  prefixes[has_red & has_grn]
}

#--- DETERMINE WHETHER AN ExPRESSION MATRIx LOOKS LIKE BETA VALUES -------------

looks_like_beta <- function(mat) {
  
  if (is.null(mat) || length(mat) == 0) return(FALSE)

  qs <- suppressWarnings(stats::quantile(
    mat,
    probs = c(0.01, 0.99),
    na.rm = TRUE
  ))

  isTRUE(qs[[1]] >= -0.05 && qs[[2]] <= 1.05)
}

#--- DETERMINE WHETHER AN ExPRESSION MATRIx LOOKS LIKE M-VALUES ----------------

looks_like_m_value <- function(mat) {
  
  if (is.null(mat) || length(mat) == 0) return(FALSE)

  qs <- suppressWarnings(stats::quantile(
    mat,
    probs = c(0.01, 0.99),
    na.rm = TRUE
  ))

  isTRUE(qs[[1]] < -0.5 || qs[[2]] > 1.5)
}

#--- CONVERT M-VALUES TO BETA VALUES -------------------------------------------

m_to_beta <- function(m) {
  
  2^m / (1 + 2^m)
}

#--- READ GEO SERIES MATRIx ExPRESSION DATA ------------------------------------

read_expr_from_series_matrix <- function(gse_id, dest = DIR_GEO_PHENO) {
  
  gse <- fetch_gse_object(gse_id, dest = dest)
  expr <- Biobase::exprs(gse)

  if (is.null(expr) || nrow(expr) == 0 || ncol(expr) == 0) return(NULL)

  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"

  expr
}

#--- NORMALIZE GEO FALLBACK ExPRESSION MATRIx TO BETA SCALE --------------------

normalize_expr_fallback_to_beta <- function(expr) {
  
  if (is.null(expr)) return(NULL)

  if (looks_like_beta(expr)) {
    return(clip_beta(expr))
  }

  if (looks_like_m_value(expr)) {
    return(clip_beta(m_to_beta(expr)))
  }

  rng <- range(expr, na.rm = TRUE)

  if (all(is.finite(rng)) && diff(rng) > 0) {
    expr2 <- (expr - rng[1]) / diff(rng)
    return(clip_beta(expr2))
  }

  NULL
}

#--- PARSE AGE VALUES FROM FREE TExT -------------------------------------------

extract_age_from_string <- function(x) {
  x <- tolower(as.character(x))
  
  out <- rep(NA_real_, length(x))
  
  # years / yo / y / age:
  pat_years <- c(
    "age[^0-9]{0,10}([0-9]+\\.?[0-9]*)",
    "([0-9]+\\.?[0-9]*)\\s*(years|year|yrs|yr|yo|y/o|y\\.o\\.)",
    "([0-9]+\\.?[0-9]*)\\s*y\\b"
  )
  
  for (pat in pat_years) {
    hit <- stringr::str_match(x, pat)[, 2]
    val <- suppressWarnings(as.numeric(hit))
    idx <- is.na(out) & !is.na(val)
    out[idx] <- val[idx]
  }
  
  # months -> years
  pat_months <- c(
    "([0-9]+\\.?[0-9]*)\\s*(months|month|mos|mo)\\b"
  )
  for (pat in pat_months) {
    hit <- stringr::str_match(x, pat)[, 2]
    val <- suppressWarnings(as.numeric(hit)) / 12
    idx <- is.na(out) & !is.na(val)
    out[idx] <- val[idx]
  }
  
  # weeks -> years
  pat_weeks <- c(
    "([0-9]+\\.?[0-9]*)\\s*(weeks|week|wks|wk)\\b"
  )
  for (pat in pat_weeks) {
    hit <- stringr::str_match(x, pat)[, 2]
    val <- suppressWarnings(as.numeric(hit)) / 52.25
    idx <- is.na(out) & !is.na(val)
    out[idx] <- val[idx]
  }
  
  # days -> years
  pat_days <- c(
    "([0-9]+\\.?[0-9]*)\\s*(days|day)\\b"
  )
  for (pat in pat_days) {
    hit <- stringr::str_match(x, pat)[, 2]
    val <- suppressWarnings(as.numeric(hit)) / 365.25
    idx <- is.na(out) & !is.na(val)
    out[idx] <- val[idx]
  }
  
  out
}

extract_age_robust <- function(ph) {

  #--- PURPOSE -----------------------------------------------------------------
  # Extract age from phenotype metadata using heuristic parsing.
  # This is part of Step 1 of the preprocessing pipeline.
  #
  # INPUT:
  #   ph - Data frame with phenotype variables (typically from fetch_gse_pheno)
  #
  # OUTPUT:
  #   Numeric vector of ages in years, with NAs for samples where age cannot
  #   be parsed. Integer ages are shifted by +0.5 to match reference pipeline.
  #
  # NOTES:
  #   - Prioritizes columns with age-like names (e.g., "age", "donor_age")
  #   - Fallback: searches all character/factor columns for age patterns
  #   - Supports multiple formats: "25 years", "25y", "25 yo", "300 months", etc.
  #   - For study-specific age extraction, use extract_age() with age_rules

  char_cols <- names(ph)[vapply(ph, function(z) is.character(z) || is.factor(z), logical(1))]

  # prioritize obvious age-like columns first
  priority_cols <- char_cols[stringr::str_detect(tolower(char_cols), "age|gestational|donor_age|chronological")]
  other_cols <- setdiff(char_cols, priority_cols)
  ordered_cols <- c(priority_cols, other_cols)

  if (length(ordered_cols) == 0) return(rep(NA_real_, nrow(ph)))

  parsed_list <- lapply(ordered_cols, function(col) {
    extract_age_from_string(ph[[col]])
  })

  age <- parsed_list[[1]]
  if (length(parsed_list) > 1) {
    for (j in 2:length(parsed_list)) {
      idx <- is.na(age) & !is.na(parsed_list[[j]])
      age[idx] <- parsed_list[[j]][idx]
    }
  }

  # Also try numeric columns with age-like names
  num_cols <- names(ph)[vapply(ph, is.numeric, logical(1))]
  num_age_cols <- num_cols[stringr::str_detect(tolower(num_cols), "age|gestational|donor_age|chronological")]
  if (length(num_age_cols) > 0) {
    for (col in num_age_cols) {
      val <- suppressWarnings(as.numeric(ph[[col]]))
      idx <- is.na(age) & !is.na(val)
      age[idx] <- val[idx]
    }
  }

  add_half_if_integer(age)
}

#--- LOAD CURATED AGE PARSING RULES --------------------------------------------

load_age_rules <- function(path = AGE_RULES_FILE) {
  
  if (!file.exists(path)) {
    warning("Age rules file not found at: ", path)
    return(tibble::tibble(
      gse = character(),
      field = character(),
      parse_type = character(),
      unit = character(),
      pattern = character(),
      notes = character()
    ))
  }

  readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::mutate(
      gse = as.character(gse),
      field = as.character(field),
      parse_type = as.character(parse_type),
      unit = as.character(unit),
      pattern = as.character(pattern),
      notes = as.character(notes)
    )
}

#--- CONVERT AGE UNITS TO YEARS ------------------------------------------------

convert_age_units_to_years <- function(x, unit) {
  
  unit <- tolower(unit %||% "years")

  dplyr::case_when(
    unit %in% c("year", "years", "yr", "yrs") ~ x,
    unit %in% c("month", "months", "mo", "mos") ~ x / 12,
    unit %in% c("week", "weeks", "wk", "wks") ~ x / 52.25,
    unit %in% c("day", "days") ~ x / 365.25,
    TRUE ~ x
  )
}

#--- PARSE AGE FROM ONE RULE ---------------------------------------------------

parse_age_by_rule <- function(ph, rule_row) {
  
  field <- rule_row$field
  parse_type <- rule_row$parse_type
  unit <- rule_row$unit
  pattern <- rule_row$pattern

  if (!field %in% names(ph)) {
    return(rep(NA_real_, nrow(ph)))
  }

  x <- ph[[field]]

  out <- switch(
    parse_type,

    numeric_column = {
      suppressWarnings(as.numeric(x))
    },

    regex_numeric = {
      hit <- stringr::str_match(as.character(x), pattern)[, 2]
      suppressWarnings(as.numeric(hit))
    },

    title_regex = {
      hit <- stringr::str_match(as.character(x), pattern)[, 2]
      suppressWarnings(as.numeric(hit))
    },

    characteristics_regex = {
      hit <- stringr::str_match(as.character(x), pattern)[, 2]
      suppressWarnings(as.numeric(hit))
    },

    {
      warning("Unknown parse_type: ", parse_type)
      rep(NA_real_, nrow(ph))
    }
  )

  convert_age_units_to_years(out, unit = unit)
}

#--- ExTRACT AGE USING CURATED RULES -------------------------------------------

extract_age_from_rules <- function(ph, gse_id, age_rules) {
  
  rules_gse <- age_rules |>
    dplyr::filter(gse == gse_id)

  if (nrow(rules_gse) == 0) {
    return(rep(NA_real_, nrow(ph)))
  }

  parsed <- rep(NA_real_, nrow(ph))

  for (i in seq_len(nrow(rules_gse))) {
    cur <- parse_age_by_rule(ph, rules_gse[i, , drop = FALSE])
    idx <- is.na(parsed) & !is.na(cur)
    parsed[idx] <- cur[idx]
  }

  parsed
}

#--- FALLBACK AGE ExTRACTION ---------------------------------------------------

extract_age_fallback <- function(ph) {
  
  char_cols <- names(ph)[vapply(
    ph,
    function(z) is.character(z) || is.factor(z),
    logical(1)
  )]

  priority_cols <- char_cols[stringr::str_detect(
    tolower(char_cols),
    "age|gestational|donor_age|chronological|characteristics|title"
  )]

  other_cols <- setdiff(char_cols, priority_cols)
  ordered_cols <- c(priority_cols, other_cols)

  if (length(ordered_cols) == 0) return(rep(NA_real_, nrow(ph)))

  age <- rep(NA_real_, nrow(ph))

  for (col in ordered_cols) {
    cur <- extract_age_from_string(ph[[col]])
    idx <- is.na(age) & !is.na(cur)
    age[idx] <- cur[idx]
  }

  num_cols <- names(ph)[vapply(ph, is.numeric, logical(1))]
  num_age_cols <- num_cols[stringr::str_detect(
    tolower(num_cols),
    "age|gestational|donor_age|chronological"
  )]

  for (col in num_age_cols) {
    cur <- suppressWarnings(as.numeric(ph[[col]]))
    idx <- is.na(age) & !is.na(cur)
    age[idx] <- cur[idx]
  }

  age
}

#--- MAIN AGE ExTRACTION WRAPPER -----------------------------------------------

extract_age <- function(ph, gse_id, age_rules) {
  
  age_rule <- extract_age_from_rules(ph, gse_id, age_rules)
  age_fallback <- extract_age_fallback(ph)

  age <- ifelse(!is.na(age_rule), age_rule, age_fallback)
  add_half_if_integer(age)
}

#--- LOAD ILLUMINA ANNOTATION FOR 27K / 450K ARRAYS ----------------------------

get_27k_annotation <- function() {
  
  ann <- minfi::getAnnotation(IlluminaHumanMethylation27kanno.ilmn12.hg19)

  tibble::tibble(
    probe_id = rownames(ann),
    chr = as.character(ann$chr)
  ) |>
    dplyr::distinct()
}

get_450k_annotation <- function() {
  
  ann <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

  tibble::tibble(
    probe_id = rownames(ann),
    chr = as.character(ann$chr)
  ) |>
    dplyr::distinct()
}

#--- LOAD CROSS-REACTIVE / MULTI-MAPPING PROBE IDS -----------------------------

load_cross_reactive_probes <- function(path = CROSS_REACTIVE_FILE) {
  
  if (!file.exists(path)) {
    message("No cross-reactive probe file found at: ", path)
    return(character())
  }

  x <- readr::read_csv(path, show_col_types = FALSE)

  if (ncol(x) == 0) return(character())

  unique(as.character(x[[1]]))
}

drop_non_autosomal_probes <- function(mat) {

  #--- PURPOSE -----------------------------------------------------------------
  # Remove probes on sex chromosomes (x, Y) if annotation available.
  # This is Step 3b of the preprocessing pipeline.
  #
  # INPUT:
  #   mat - Numeric matrix (probes × samples) with rownames = probe IDs
  #
  # OUTPUT:
  #   Currently returns input unchanged (stub function)
  #
  # NOTES:
  #   - This function is a placeholder for future annotation-based filtering
  #   - When chromosome annotation becomes available, will filter by chr annotation
  #   - For now, the global CpG intersection handles most filtering implicitly

  # Conservative rule:
  # If we later add an annotation table of probe->chr, replace this stub
  # with exact chromosome-based filtering.
  #
  # For now we leave probe set unchanged here and rely on global intersection
  # unless we provide an annotation-based filter externally.
  mat
}

drop_cross_reactive_probes <- function(mat, cross_reactive = character()) {

  #--- PURPOSE -----------------------------------------------------------------
  # Remove probes with cross-reactive or multi-mapping issues.
  # This is Step 3c of the preprocessing pipeline.
  #
  # INPUT:
  #   mat - Numeric matrix (probes × samples) with rownames = probe IDs
  #   cross_reactive - Character vector of probe IDs to remove
  #
  # OUTPUT:
  #   Matrix with problematic probes removed. If cross_reactive is empty,
  #   returns input unchanged.
  #
  # NOTES:
  #   - cross_reactive is typically loaded from data/reference/cross_reactive_probes.csv
  #   - These probes have ambiguous mappings or cross-hybridize in methylation assays
  #   - Removing them improves specificity of downstream epigenetic clock predictions

  if (length(cross_reactive) == 0) return(mat)
  keep <- setdiff(rownames(mat), cross_reactive)
  mat[keep, , drop = FALSE]
}

#--- BUILD ExACT COMMON CPG SET ACROSS 27K AND 450K ----------------------------

common_cpgs_27k_450k <- function(drop_sex = TRUE,
                                 drop_cross_reactive = TRUE,
                                 cross_reactive_file = CROSS_REACTIVE_FILE,
                                 save_manifest = TRUE) {
  
  ann27 <- get_27k_annotation()
  ann450 <- get_450k_annotation()

  common <- dplyr::inner_join(
    ann27 |> dplyr::rename(chr_27k = chr),
    ann450 |> dplyr::rename(chr_450k = chr),
    by = "probe_id"
  ) |>
    dplyr::filter(!is.na(chr_27k), !is.na(chr_450k)) |>
    dplyr::filter(chr_27k == chr_450k) |>
    dplyr::mutate(chr = chr_27k)

  if (drop_sex) {
    common <- common |>
      dplyr::filter(!toupper(chr) %in% c("CHRx", "x", "CHRY", "Y"))
  }

  if (drop_cross_reactive) {
    xr <- load_cross_reactive_probes(cross_reactive_file)
    if (length(xr) > 0) {
      common <- common |>
        dplyr::filter(!(probe_id %in% xr))
    }
  }

  common_ids <- sort(unique(common$probe_id))

  if (save_manifest) {
    manifest <- common |>
      dplyr::select(probe_id, chr) |>
      dplyr::arrange(probe_id)

    readr::write_csv(
      manifest,
      fs::path(DIR_OUT, "common_cpg_manifest.csv")
    )
    qs2::qs_save(
      manifest,
      fs::path(DIR_OUT, "common_cpg_manifest.qs")
    )
  }

  common_ids
}

#--- RESTRICT MATRIx TO COMMON CPGS --------------------------------------------

restrict_to_common_cpgs <- function(mat, common_cpgs) {
  
  if (is.null(mat)) return(NULL)

  rn <- rownames(mat)

  if (is.null(rn)) {
    warning("Matrix has no rownames; cannot restrict to CpGs.")
    return(NULL)
  }

  keep <- common_cpgs[common_cpgs %in% rn]

  if (length(keep) == 0) {
    warning("No common CpGs found in matrix.")
    return(NULL)
  }

  mat[keep, , drop = FALSE]
}

#--- DEDUPLICATE PROBE ROWS ----------------------------------------------------

dedup_probe_rows <- function(mat) {
  
  if (is.null(mat)) return(NULL)

  rn <- rownames(mat)
  if (anyDuplicated(rn) == 0) return(mat)

  dt <- data.table::as.data.table(mat, keep.rownames = "probe_id")
  dt <- dt[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = probe_id]

  out <- as.matrix(dt[, -"probe_id"])
  rownames(out) <- dt$probe_id

  out
}

#--- BUILD SAMPLE MATCHING KEYS FROM PHENOTYPE TABLE ---------------------------

build_pheno_sample_keys <- function(ph) {
  candidate_cols <- c(
    "geo_accession", "title", "supplementary_file", "supplementary_file_1",
    "supplementary_file_2", "source_name_ch1", "source_name_ch1.1"
  )
  candidate_cols <- intersect(candidate_cols, names(ph))
  
  keys <- rep(NA_character_, nrow(ph))
  
  for (col in candidate_cols) {
    cur <- sample_key_normalize(ph[[col]])
    idx <- is.na(keys) & !is.na(cur) & nzchar(cur)
    keys[idx] <- cur[idx]
  }
  
  # always include geo_accession-based fallback
  if ("geo_accession" %in% names(ph)) {
    cur <- sample_key_normalize(ph$geo_accession)
    idx <- is.na(keys) & !is.na(cur) & nzchar(cur)
    keys[idx] <- cur[idx]
  }
  
  keys
}

#--- ALIGN PHENOTYPE ROWS TO MATRIx COLUMNS ------------------------------------
#
# PURPOSE:
#   Match phenotype records to methylation matrix columns by sample keys.
#   This is Step 4 of the preprocessing pipeline.
#
# INPUT:
#   ph  - Phenotype data frame (typically from fetch_gse_pheno + age extraction)
#   mat - Methylation matrix (probes × samples) from IDAT or series matrix
#
# OUTPUT:
#   List with two elements:
#     $ph  - Phenotype table with matched rows only, with geo_accession as stable IDs
#     $mat - Methylation matrix with matched columns only, using geo_accession as colnames
#
# NOTES:
#   - Matches by normalized sample keys (uppercase, punctuation removed, etc.)
#   - Removes samples not present in both phenotype and matrix
#   - Deduplicates by sample ID if needed
#   - Returns (0-row, 0-column) result if no samples match
#   - Validates that ph and mat are perfectly aligned after matching

align_pheno_to_matrix <- function(ph, mat) {
  if (is.null(mat) || is.null(ph) || ncol(mat) == 0 || nrow(ph) == 0) {
    return(list(ph = ph[0, , drop = FALSE], mat = mat[, 0, drop = FALSE]))
  }

  ph <- ph |> dplyr::mutate(sample_key = build_pheno_sample_keys(ph))

  col_keys <- sample_key_normalize(colnames(mat))
  names(col_keys) <- colnames(mat)

  idx <- match(col_keys, ph$sample_key)

  keep_cols <- !is.na(idx)
  mat2 <- mat[, keep_cols, drop = FALSE]
  idx2 <- idx[keep_cols]
  ph2 <- ph[idx2, , drop = FALSE]

  # replace matrix colnames with stable sample IDs
  stable_ids <- if ("geo_accession" %in% names(ph2)) ph2$geo_accession else ph2$sample_key
  colnames(mat2) <- stable_ids

  # de-duplicate sample IDs if needed
  if (anyDuplicated(colnames(mat2)) > 0) {
    dup_id <- duplicated(colnames(mat2))
    mat2 <- mat2[, !dup_id, drop = FALSE]
    ph2 <- ph2[!dup_id, , drop = FALSE]
  }

  stopifnot(identical(colnames(mat2), if ("geo_accession" %in% names(ph2)) ph2$geo_accession else ph2$sample_key))

  list(ph = ph2, mat = mat2)
}

#--- Read and Preprocess IDATs with sesame -----------------------------------
#
# PURPOSE:
#   Read methylation beta values from raw IDAT files using sesame.
#   This is the primary method in Step 2 of the preprocessing pipeline.
#   Automatically decompresses .idat.gz files if needed.
#
# INPUT:
#   gse_dir - Directory containing downloaded IDAT files for a study
#
# OUTPUT:
#   Numeric matrix (probes × samples) with beta values in [0, 1].
#   Returns NULL if no complete IDAT pairs are found.
#
# NOTES:
#   - Searches recursively for paired IDAT files (_Red.idat + _Grn.idat)
#   - Handles compressed .idat.gz files by automatically decompressing
#   - Uses sesame::openSesame() with default settings for preprocessing
#   - Applies beta value clipping to ensure values stay in [0, 1]
#   - If IDAT processing fails, returns NULL (triggering fallback to series matrix)

read_beta_from_idats <- function(gse_dir) {

  # First, decompress any .idat.gz files to .idat
  gz_idats <- fs::dir_ls(
    gse_dir,
    recurse = TRUE,
    regexp = "\\.idat\\.gz$",
    type = "file"
  )

  if (length(gz_idats) > 0) {
    for (gz_file in gz_idats) {
      uncompressed_file <- stringr::str_replace(gz_file, "\\.gz$", "")
      # Only decompress if uncompressed file doesn't already exist
      if (!file.exists(uncompressed_file)) {
        tryCatch({
          R.utils::gunzip(gz_file, destname = uncompressed_file, remove = FALSE)
        }, error = function(e) {
          message("Failed to decompress ", basename(gz_file), ": ", conditionMessage(e))
        })
      }
    }
  }

  prefixes <- find_idat_prefixes(gse_dir)
  if (length(prefixes) == 0) return(NULL)

  # Disable parallelization to ensure sesame can access cached data
  tryCatch({
    BiocParallel::register(BiocParallel::SerialParam())
    beta <- sesame::openSesame(prefixes, func = sesame::getBetas)
  }, error = function(e) {
    message("sesame preprocessing failed in ", gse_dir, ": ", conditionMessage(e))
    return(NULL)
  })

  if (is.null(beta)) return(NULL)

  beta <- as.matrix(beta)
  storage.mode(beta) <- "numeric"
  clip_beta(beta)
}

#-------------------------------------------------------------------------------
# Compute Global Shared CpG Panel
#-------------------------------------------------------------------------------
global_common_cpgs <- function(proc_list) {
  mats <- purrr::keep(proc_list, ~ !is.null(.x$mat) && nrow(.x$mat) > 0 && ncol(.x$mat) > 0)
  if (length(mats) == 0) return(character())
  Reduce(intersect, lapply(mats, function(x) rownames(x$mat))) |> sort()
}

#-------------------------------------------------------------------------------
# Restrict Processed Study Objects to a Shared CpG Set
#-------------------------------------------------------------------------------
subset_proc_list_to_cpgs <- function(proc_list, common_cpgs) {
  lapply(proc_list, function(x) {
    if (is.null(x$mat)) return(x)
    x$mat <- x$mat[common_cpgs, , drop = FALSE]
    x
  })
}

#-------------------------------------------------------------------------------
# Combine Study-Level Matrices
#-------------------------------------------------------------------------------
combine_mats <- function(proc_list) {
  mats <- purrr::keep(proc_list, ~ !is.null(.x$mat) && ncol(.x$mat) > 0)
  if (length(mats) == 0) return(NULL)
  do.call(cbind, lapply(mats, `[[`, "mat"))
}

#-------------------------------------------------------------------------------
# Combine Study-Level Phenotype Tables
#-------------------------------------------------------------------------------
combine_pheno <- function(proc_list) {
  phs <- purrr::keep(proc_list, ~ !is.null(.x$ph) && nrow(.x$ph) > 0)
  if (length(phs) == 0) return(tibble::tibble())
  dplyr::bind_rows(lapply(phs, `[[`, "ph"))
}

#-------------------------------------------------------------------------------
# Impute Missing Values Using Training Probe Medians
#-------------------------------------------------------------------------------
impute_by_train_probe_median <- function(train_cpg_x_sample, test_cpg_x_sample = NULL) {
  train_medians <- apply(train_cpg_x_sample, 1, function(x) stats::median(x, na.rm = TRUE))
  train_imp <- train_cpg_x_sample
  for (i in seq_len(nrow(train_imp))) {
    miss <- is.na(train_imp[i, ])
    if (any(miss)) train_imp[i, miss] <- train_medians[i]
  }
  
  test_imp <- NULL
  if (!is.null(test_cpg_x_sample)) {
    test_imp <- test_cpg_x_sample
    stopifnot(identical(rownames(test_imp), rownames(train_imp)))
    for (i in seq_len(nrow(test_imp))) {
      miss <- is.na(test_imp[i, ])
      if (any(miss)) test_imp[i, miss] <- train_medians[i]
    }
  }
  
  list(
    train = train_imp,
    test = test_imp,
    train_probe_medians = train_medians
  )
}

#-------------------------------------------------------------------------------
# Convert to Scorcher-Ready x/Y Objects
#-------------------------------------------------------------------------------
to_scorcher_xy <- function(beta_cpg_x_sample, ph) {
  stopifnot(identical(colnames(beta_cpg_x_sample), ph$sample_id))
  
  x <- t(beta_cpg_x_sample)
  storage.mode(x) <- "numeric"
  
  y <- ph$age
  names(y) <- ph$sample_id
  
  list(
    x = x,
    y = y,
    sample_ids = ph$sample_id,
    feature_names = colnames(x)
  )
}

#-------------------------------------------------------------------------------
# Save Objects with qs
#-------------------------------------------------------------------------------
save_qs <- function(object, path) {
  qs2::qs_save(object, path)
}

#-------------------------------------------------------------------------------
# Build Age-Parsing QC Summary
#-------------------------------------------------------------------------------
build_age_qc <- function(proc_list) {
  safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  safe_med <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
  safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

  tibble::tibble(
    gse = vapply(proc_list, `[[`, character(1), "gse"),
    method = vapply(proc_list, `[[`, character(1), "method"),
    n_pheno = vapply(proc_list, function(x) nrow(x$ph), integer(1)),
    n_age_nonmissing = vapply(proc_list, function(x) sum(!is.na(x$ph$age)), integer(1)),
    pct_age_nonmissing = 100 * n_age_nonmissing / pmax(n_pheno, 1),
    age_min = vapply(proc_list, function(x) safe_min(x$ph$age), numeric(1)),
    age_median = vapply(proc_list, function(x) safe_med(x$ph$age), numeric(1)),
    age_max = vapply(proc_list, function(x) safe_max(x$ph$age), numeric(1))
  )
}

#-------------------------------------------------------------------------------
# Build Final Summary
#-------------------------------------------------------------------------------
build_final_summary <- function(train_xy, test_xy, common_cpgs) {
  tibble::tibble(
    dataset = c("train", "test"),
    n_samples = c(nrow(train_xy$x), nrow(test_xy$x)),
    n_cpg = c(ncol(train_xy$x), ncol(test_xy$x)),
    age_min = c(min(train_xy$y, na.rm = TRUE), min(test_xy$y, na.rm = TRUE)),
    age_median = c(stats::median(train_xy$y, na.rm = TRUE),
                   stats::median(test_xy$y, na.rm = TRUE)),
    age_max = c(max(train_xy$y, na.rm = TRUE), max(test_xy$y, na.rm = TRUE)),
    target_common_cpgs = length(common_cpgs)
  )
}

#=== CHECKPOINT SAFETY HELPERS ==================================================
#
# These functions wrap qs2 save/load to handle corrupted checkpoint files
# and provide informative error messages
#

safe_checkpoint_read <- function(checkpoint_path) {
  
  if (!file.exists(checkpoint_path)) {
    return(NULL)
  }
  
  tryCatch({
    qs2::qs_read(checkpoint_path)
  }, error = function(e) {
    warning("\n x Corrupted checkpoint at: ", checkpoint_path,
            "\n  Error: ", conditionMessage(e),
            "\n  Attempting recovery by deleting and restarting step...")
    
    tryCatch({
      unlink(checkpoint_path)
      message("   + Corrupted file deleted")
    }, error = function(e2) {
      message("   x Could not delete corrupted file: ", conditionMessage(e2))
    })
    
    return(NULL)
  })
}

safe_checkpoint_save <- function(object, checkpoint_path) {
  
  tryCatch({
    qs2::qs_save(object, checkpoint_path)
    
    # Verify write succeeded
    if (!file.exists(checkpoint_path)) {
      stop("Checkpoint file not created at: ", checkpoint_path,
           " (Check file permissions)")
    }
    
    invisible(NULL)
    
  }, error = function(e) {
    stop(" x Failed to save checkpoint to: ", checkpoint_path, "\n",
         "  Error: ", conditionMessage(e), "\n",
         "  Likely causes:\n",
         "    - Insufficient disk space\n",
         "    - File permissions issue\n",
         "    - Invalid filename\n",
         "  Please check and try again.")
  })
}

