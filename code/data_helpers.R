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
#             - Age parsing from phenotype metadata columns
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
#             - Installed Illumina annotation packages
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
#                 - For publication-quality reproduction, parsed age columns
#                   should be reviewed against study documentation.
#
#             3. Platform harmonization:
#
#                 - The intended preprocessing uses CpGs shared across the
#                   Illumina 27K and 450K arrays. The HM450 sesameData manifest
#                   is used to remove sex chromosome probes and duplicate
#                   genomic coordinates before the final cross-study intersection.
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

progress_message <- function(..., appendLF = TRUE) {
  message("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., appendLF = appendLF)
  flush.console()
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

  key <- x |>
    as.character() |>
    basename() |>
    stringr::str_replace("\\.gz$", "") |>
    stringr::str_replace("\\.idat$", "") |>
    stringr::str_replace("(_Red|_Grn)$", "") |>
    stringr::str_replace("\\s+", "_") |>
    stringr::str_replace_all("[^A-Za-z0-9._-]", "_") |>
    toupper()

  gsm_key <- stringr::str_extract(key, "GSM[0-9]+")
  key[!is.na(gsm_key)] <- gsm_key[!is.na(gsm_key)]
  key
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

extract_geo_raw <- function(gse_id, src = DIR_RAW_DOWNLOADS, dest = DIR_RAW_EXTRACTED) {

  #--- PURPOSE -----------------------------------------------------------------
  # STEP 2: Extract top-level archives from downloaded files.
  #
  # This function extracts tar/tar.gz/tgz files from the download directory.
  # Extracted files go to a SEPARATE directory to keep original archives intact.
  #
  # INPUT:
  #   gse_id - GEO Series accession
  #   src    - Source directory where archives were downloaded (DIR_RAW_DOWNLOADS)
  #   dest   - Destination directory for extracted files (DIR_RAW_EXTRACTED)
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

extract_nested_tar_files <- function(gse_id, src = DIR_RAW_EXTRACTED, dest = DIR_RAW_TAR) {

  #--- PURPOSE -----------------------------------------------------------------
  # STEP 3: Extract secondary/nested TAR archives (e.g., _RAW.tar).
  #
  # Some GEO studies have nested tar files that need secondary extraction.
  # These typically contain raw platform files (_RAW.tar) or additional archives.
  #
  # INPUT:
  #   gse_id - GEO Series accession
  #   src    - Source directory with already-extracted files (DIR_RAW_EXTRACTED)
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
                                  src_dirs = c(DIR_RAW_EXTRACTED, DIR_RAW_TAR),
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

#--- READ PROCESSED SUPPLEMENTARY BETA MATRICES --------------------------------

read_beta_from_supplementary_matrix <- function(gse_id, src = DIR_RAW_DOWNLOADS) {
  gse_dir <- fs::path(src, gse_id)
  if (!dir.exists(gse_dir)) return(NULL)

  candidates <- fs::dir_ls(
    gse_dir,
    regexp = "(MatrixProcessed|Beta).*\\.txt(\\.gz)?$",
    type = "file"
  )

  if (length(candidates) == 0) return(NULL)

  priority <- dplyr::case_when(
    stringr::str_detect(basename(candidates), "MatrixProcessed") ~ 1L,
    stringr::str_detect(basename(candidates), "Normalized") ~ 2L,
    stringr::str_detect(basename(candidates), "Beta") ~ 3L,
    TRUE ~ 4L
  )
  candidates <- candidates[order(priority, basename(candidates))]

  for (path in candidates) {
    beta <- tryCatch({
      progress_message("Reading supplementary matrix: ", path)
      t_start <- Sys.time()
      dt <- data.table::fread(
        path,
        skip = "ID_REF",
        data.table = FALSE,
        check.names = FALSE,
        showProgress = FALSE
      )
      progress_message(
        "Finished supplementary matrix read: ", basename(path),
        " (", nrow(dt), " rows x ", ncol(dt), " cols, ",
        round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1),
        " sec)"
      )

      if (nrow(dt) == 0 || ncol(dt) < 2) {
        NULL
      } else {
        probe_id <- as.character(dt[[1]])
        value_idx <- seq_len(ncol(dt))[-1]
        value_names <- names(dt)[value_idx]
        beta_keep <- !stringr::str_detect(
          tolower(value_names),
          "detect|pval|p\\.val|pvalue"
        )
        beta_idx <- value_idx[beta_keep]

        if (length(beta_idx) == 0) {
          NULL
        } else {
          mat <- as.matrix(dt[, beta_idx, drop = FALSE])
          storage.mode(mat) <- "numeric"
          rownames(mat) <- probe_id
          colnames(mat) <- value_names[beta_keep]

          normalize_expr_fallback_to_beta(mat)
        }
      }
    }, error = function(e) {
      message("supplementary matrix read failed for ", basename(path), ": ",
              conditionMessage(e))
      NULL
    })

    if (!is.null(beta)) {
      message(" (supplementary matrix: ", basename(path), ")", appendLF = FALSE)
      return(beta)
    }
  }

  NULL
}

#--- PARSE AGE VALUES FROM FREE TExT -------------------------------------------

extract_age_from_string <- function(x) {
  x <- tolower(as.character(x))
  
  out <- rep(NA_real_, length(x))

  newborn_idx <- stringr::str_detect(x, "\\b(newborn|neonate|neonatal)\\b")
  out[newborn_idx] <- 0
  
  # years / yo / y / age:
  pat_years <- c(
    "age\\s*(at)?\\s*(recruitment|diagnosis|draw|sample|bl)?\\s*(\\([^)]*\\))?\\s*[:=]?\\s*([0-9]+\\.?[0-9]*)",
    "age[^0-9]{0,10}([0-9]+\\.?[0-9]*)",
    "([0-9]+\\.?[0-9]*)\\s*(years|year|yrs|yr|yo|y/o|y\\.o\\.)",
    "([0-9]+\\.?[0-9]*)\\s*y\\b"
  )
  
  for (pat in pat_years) {
    m <- stringr::str_match(x, pat)
    hit <- m[, ncol(m)]
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

is_age_metadata_column <- function(cols) {
  cols <- tolower(as.character(cols))

  stringr::str_detect(
    cols,
    "(^age)|([^a-z]age([^a-z]|$))|ageat|age_at|age at|agebl|donor_age|chronological|gestational"
  )
}

infer_age_unit_from_column <- function(col) {
  col <- tolower(as.character(col))

  if (stringr::str_detect(col, "month|\\bmo\\b|_mo\\b|months")) {
    return("months")
  }

  if (stringr::str_detect(col, "week|\\bwk\\b|_wk\\b|weeks")) {
    return("weeks")
  }

  if (stringr::str_detect(col, "day|days")) {
    return("days")
  }

  if (stringr::str_detect(col, "gestational")) {
    return("weeks")
  }

  "years"
}

parse_age_column <- function(x, col) {
  unit <- infer_age_unit_from_column(col)

  if (is.numeric(x)) {
    return(convert_age_units_to_years(suppressWarnings(as.numeric(x)), unit))
  }

  parsed <- extract_age_from_string(x)
  plain_numeric <- suppressWarnings(as.numeric(
    stringr::str_replace(as.character(x), "^\\s*[<>]=?\\s*", "")
  ))
  idx <- is.na(parsed) & !is.na(plain_numeric)
  parsed[idx] <- convert_age_units_to_years(plain_numeric[idx], unit)

  parsed
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
  #   - Prioritizes direct age-like columns present in each study
  #     (e.g., "age", "donor_age", "age_months")
  #   - Fallback: searches other character/factor columns for age patterns
  #   - Supports multiple formats: "25 years", "25y", "25 yo", "300 months", etc.
  #   - Plain numeric strings in age-like columns are treated as years unless
  #     the column name indicates months, weeks, or days

  if (nrow(ph) == 0) return(numeric())

  age <- rep(NA_real_, nrow(ph))

  numeric_cols <- names(ph)[vapply(ph, is.numeric, logical(1))]
  numeric_age_cols <- numeric_cols[is_age_metadata_column(numeric_cols)]

  for (col in numeric_age_cols) {
    cur <- parse_age_column(ph[[col]], col)
    idx <- is.na(age) & !is.na(cur)
    age[idx] <- cur[idx]
  }

  char_cols <- names(ph)[vapply(ph, function(z) is.character(z) || is.factor(z), logical(1))]
  char_age_cols <- char_cols[is_age_metadata_column(char_cols)]

  for (col in char_age_cols) {
    cur <- parse_age_column(ph[[col]], col)
    idx <- is.na(age) & !is.na(cur)
    age[idx] <- cur[idx]
  }

  fallback_cols <- setdiff(char_cols, char_age_cols)

  for (col in fallback_cols) {
    cur <- extract_age_from_string(ph[[col]])
    idx <- is.na(age) & !is.na(cur)
    age[idx] <- cur[idx]
  }

  add_half_if_integer(age)
}

#--- INFER CONTROL / CASE STATUS FROM PHENOTYPE METADATA -----------------------

infer_control_status <- function(ph) {

  if (nrow(ph) == 0) return(logical())

  cols <- names(ph)[vapply(ph, function(z) is.character(z) || is.factor(z), logical(1))]
  cols <- cols[stringr::str_detect(
    tolower(cols),
    "case|control|disease|diagnos|group|status|hiv|ards|sample type|used_in_analysis"
  )]

  status <- rep(NA, nrow(ph))

  for (col in cols) {
    x <- tolower(as.character(ph[[col]]))
    x <- stringr::str_squish(x)

    case_hit <- stringr::str_detect(
      x,
      "case|patient|disease|diagnos|cancer|tumou?r|carcinoma|multiple sclerosis|\\bms\\b|\\bibd\\b|crohn|colitis|arthritis|hiv\\s*positive|\\bhiv\\s*1\\b|ards|pneumonia|sepsis|affected"
    )
    control_hit <- stringr::str_detect(
      x,
      "control|healthy|normal|unaffected|non[- ]?case|hiv\\s*negative|\\bhiv\\s*0\\b|no disease|no diagnosis|none"
    )

    idx <- is.na(status) & control_hit
    status[idx] <- TRUE

    idx <- is.na(status) & case_hit
    status[idx] <- FALSE
  }

  status
}

#--- CONVERT AGE UNITS TO YEARS ------------------------------------------------

convert_age_units_to_years <- function(x, unit) {
  
  unit <- tolower(unit %||% "years")

  if (unit %in% c("year", "years", "yr", "yrs")) {
    return(x)
  }

  if (unit %in% c("month", "months", "mo", "mos")) {
    return(x / 12)
  }

  if (unit %in% c("week", "weeks", "wk", "wks")) {
    return(x / 52.25)
  }

  if (unit %in% c("day", "days")) {
    return(x / 365.25)
  }

  x
}

#--- LOAD SESAME MANIFESTS FOR 27K / 450K ARRAYS -------------------------------

get_sesame_manifest <- function(platform) {
  gr <- sesameData::sesameData_getManifestGRanges(platform)

  tibble::tibble(
    probe_id = names(gr),
    chr = as.character(GenomeInfoDb::seqnames(gr)),
    pos = GenomicRanges::start(gr)
  ) |>
    dplyr::filter(!is.na(probe_id), nzchar(probe_id)) |>
    dplyr::distinct()
}

normalize_chr <- function(chr) {
  chr |>
    as.character() |>
    toupper() |>
    stringr::str_replace("^CHR", "")
}

is_autosomal_chr <- function(chr) {
  normalize_chr(chr) %in% as.character(seq_len(22))
}

#--- BUILD AUTOSOMAL CPG PANEL FROM SESAME HM450 MANIFEST ----------------------

autosomal_cpg_panel_hm450 <- function(drop_sex = TRUE,
                                      drop_ambiguous = TRUE,
                                      save_manifest = TRUE) {
  
  panel <- get_sesame_manifest("HM450") |>
    dplyr::mutate(
      chr = normalize_chr(chr)
    ) |>
    dplyr::filter(!is.na(chr), !is.na(pos))

  if (drop_sex) {
    panel <- panel |>
      dplyr::filter(is_autosomal_chr(chr))
  }

  if (drop_ambiguous) {
    panel <- panel |>
      dplyr::add_count(chr, pos, name = "n_probes_at_coordinate") |>
      dplyr::filter(n_probes_at_coordinate == 1)
  } else {
    panel <- panel |>
      dplyr::mutate(n_probes_at_coordinate = NA_integer_)
  }

  panel_ids <- sort(unique(panel$probe_id))

  if (save_manifest) {
    manifest <- panel |>
      dplyr::select(
        probe_id,
        chr,
        pos,
        n_probes_at_coordinate
      ) |>
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

  panel_ids
}

#--- RESTRICT MATRIx TO COMMON CPGS --------------------------------------------

restrict_to_common_cpgs <- function(mat, common_cpgs) {
  
  if (is.null(mat)) return(NULL)

  if (is.null(dim(mat)) || length(dim(mat)) != 2) {
    warning("Methylation object is not a 2D matrix; cannot restrict to CpGs.")
    return(NULL)
  }

  mat <- as.matrix(mat)

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
#   mat - Methylation matrix (probes x samples) from IDAT or series matrix
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
  empty_alignment <- function(ph, mat) {
    ph0 <- if (is.null(ph)) {
      data.frame()
    } else {
      as.data.frame(ph)[0, , drop = FALSE]
    }

    mat0 <- if (!is.null(mat) && !is.null(dim(mat)) && length(dim(mat)) == 2) {
      as.matrix(mat)[, 0, drop = FALSE]
    } else {
      matrix(nrow = 0, ncol = 0)
    }

    list(ph = ph0, mat = mat0)
  }

  if (is.null(mat) || is.null(ph)) {
    return(empty_alignment(ph, mat))
  }

  if (is.null(dim(mat)) || length(dim(mat)) != 2) {
    message("    [DEBUG] Methylation object is not a 2D matrix; skipping alignment.")
    return(empty_alignment(ph, mat))
  }

  mat <- as.matrix(mat)
  ph <- as.data.frame(ph)

  mat_cols <- ncol(mat)
  ph_rows <- nrow(ph)

  if (is.na(mat_cols) || is.na(ph_rows) || mat_cols == 0 || ph_rows == 0) {
    return(empty_alignment(ph, mat))
  }

  ph <- ph |> dplyr::mutate(sample_key = build_pheno_sample_keys(ph))

  col_keys <- sample_key_normalize(colnames(mat))
  names(col_keys) <- colnames(mat)

  idx <- match(col_keys, ph$sample_key)

  keep_cols <- !is.na(idx)
  
  # Debug: Check if any samples matched
  n_matched <- sum(keep_cols)
  if (n_matched == 0) {
    message("    [DEBUG] No samples matched during alignment!")
    message("      Matrix column names sample: ", paste(head(colnames(mat), 3), collapse = ", "))
    message("      Phenotype sample keys sample: ", paste(head(ph$sample_key, 3), collapse = ", "))
    message("      This likely means phenotype and matrix use different sample identifiers.")

    sample_num <- suppressWarnings(as.integer(stringr::str_match(
      colnames(mat),
      stringr::regex("^sample[ ._-]*([0-9]+)$", ignore_case = TRUE)
    )[, 2]))

    if (length(sample_num) == ncol(mat) &&
        all(!is.na(sample_num)) &&
        setequal(sample_num, seq_len(ncol(mat))) &&
        ncol(mat) == nrow(ph)) {
      message("      Falling back to ordinal Sample N matching.")

      ord <- order(sample_num)
      mat2 <- mat[, ord, drop = FALSE]
      ph2 <- ph[seq_len(ncol(mat2)), , drop = FALSE]
      stable_ids <- if ("geo_accession" %in% names(ph2)) ph2$geo_accession else ph2$sample_key
      colnames(mat2) <- as.character(stable_ids)

      return(list(ph = ph2, mat = mat2))
    }

    return(empty_alignment(ph, mat))
  }
  
  mat2 <- mat[, keep_cols, drop = FALSE]
  idx2 <- idx[keep_cols]
  ph2 <- ph[idx2, , drop = FALSE]

  # replace matrix colnames with stable sample IDs
  stable_ids <- if ("geo_accession" %in% names(ph2)) ph2$geo_accession else ph2$sample_key
  colnames(mat2) <- as.character(stable_ids)

  # de-duplicate sample IDs if needed
  # Handle duplicates in both matrix colnames AND phenotype
  if (anyDuplicated(colnames(mat2)) > 0 || anyDuplicated(stable_ids) > 0) {
    # Keep first occurrence of each unique ID
    dup_id <- duplicated(colnames(mat2))
    
    keep_idx <- !dup_id
    mat2 <- mat2[, keep_idx, drop = FALSE]
    ph2 <- ph2[keep_idx, , drop = FALSE]
    
    # Re-compute stable_ids after filtering
    stable_ids <- if ("geo_accession" %in% names(ph2)) ph2$geo_accession else ph2$sample_key
    colnames(mat2) <- as.character(stable_ids)
  }

  # Final check - both lengths and values should match
  expected_ids <- if ("geo_accession" %in% names(ph2)) ph2$geo_accession else ph2$sample_key
  expected_ids <- as.character(expected_ids)
  mat_colnames <- as.character(colnames(mat2))
  
  if (!identical(mat_colnames, expected_ids)) {
    mat_names <- mat_colnames
    ph_names <- expected_ids
    
    message("    [ERROR] Sample ID mismatch during alignment:")
    message("      Matrix columns: ", length(mat_names), " samples")
    message("      Phenotype rows: ", length(ph_names), " samples")
    
    if (length(mat_names) != length(ph_names)) {
      message("      Length mismatch! Matrix and phenotype have different sizes.")
    } else {
      # Same length but different order or values
      mismatch_idx <- which(mat_names != ph_names)
      if (length(mismatch_idx) > 0) {
        message("      Order or value mismatch at positions: ", 
          paste(head(mismatch_idx, 5), collapse = ", "))
        message("      First few matrix names: ", paste(head(mat_names, 3), collapse = ", "))
        message("      First few pheno names: ", paste(head(ph_names, 3), collapse = ", "))
      }
    }
    
    stop("Sample ID mismatch in alignment for this study. Check phenotype duplicates or sample_key generation.")
  }

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
#   Numeric matrix (probes x samples) with beta values in [0, 1].
#   Returns NULL if no complete IDAT pairs are found.
#
# NOTES:
#   - Searches recursively for paired IDAT files (_Red.idat + _Grn.idat)
#   - Handles compressed .idat.gz files by automatically decompressing
#   - Uses sesame::openSesame() with default settings for preprocessing
#   - Applies beta value clipping to ensure values stay in [0, 1]
#   - If IDAT processing fails, returns NULL (triggering fallback to series matrix)

read_beta_from_idats <- function(gse_dir, force_samplewise = FALSE) {
  t0 <- Sys.time()
  progress_message("      IDAT scan started: ", gse_dir)

  # First, decompress any .idat.gz files to .idat
  gz_idats <- fs::dir_ls(
    gse_dir,
    recurse = TRUE,
    regexp = "\\.idat\\.gz$",
    type = "file"
  )

  progress_message("      Found ", length(gz_idats), " compressed IDAT file(s)")

  if (length(gz_idats) > 0) {
    for (gz_file in gz_idats) {
      uncompressed_file <- stringr::str_replace(gz_file, "\\.gz$", "")
      # Only decompress if uncompressed file doesn't already exist
      if (!file.exists(uncompressed_file)) {
        progress_message("      Decompressing ", basename(gz_file))
        tryCatch({
          R.utils::gunzip(gz_file, destname = uncompressed_file, remove = FALSE)
        }, error = function(e) {
          message("Failed to decompress ", basename(gz_file), ": ", conditionMessage(e))
        })
      }
    }
  }

  progress_message("      Finding complete Red/Green IDAT pairs")
  prefixes <- find_idat_prefixes(gse_dir)
  progress_message("      Found ", length(prefixes), " complete IDAT pair(s)")
  if (length(prefixes) == 0) return(NULL)

  # Disable parallelization to ensure sesame can access cached data.
  if (force_samplewise) {
    progress_message("      Skipping batch IDAT preprocessing; using sample-by-sample recovery")
    beta <- NULL
  } else {
    progress_message("      Starting batch sesame::openSesame() for ", length(prefixes), " sample(s)")
    beta <- tryCatch({
      BiocParallel::register(BiocParallel::SerialParam())
      sesame::openSesame(prefixes, func = sesame::getBetas)
    }, error = function(e) {
      message("sesame preprocessing failed in ", gse_dir, ": ", conditionMessage(e))
      NULL
    })
    progress_message(
      "      Batch IDAT preprocessing finished in ",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
      " min"
    )
  }

  if (is.null(beta) && length(prefixes) > 1) {
    message("  Retrying IDAT preprocessing sample-by-sample...")

    beta_list <- list()
    failed_prefixes <- character()

    for (i in seq_along(prefixes)) {
      prefix <- prefixes[[i]]

      if (i == 1 || i %% 25 == 0 || i == length(prefixes)) {
        progress_message(
          "  IDAT recovery progress: ", i, "/", length(prefixes),
          " (", basename(prefix), ")"
        )
      }

      cur <- tryCatch({
        sesame::openSesame(prefix, func = sesame::getBetas)
      }, error = function(e) {
        failed_prefixes <<- c(failed_prefixes, basename(prefix))
        NULL
      })

      if (is.null(cur)) next

      if (is.null(dim(cur))) {
        cur_names <- names(cur)
        cur <- matrix(cur, ncol = 1)
        rownames(cur) <- cur_names
        colnames(cur) <- basename(prefix)
      } else {
        cur <- as.matrix(cur)
      }

      storage.mode(cur) <- "numeric"
      beta_list[[basename(prefix)]] <- cur
    }

    progress_message(
      "      Sample-by-sample IDAT recovery finished in ",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
      " min; recovered ", length(beta_list), "/", length(prefixes), " sample(s)"
    )

    if (length(failed_prefixes) > 0) {
      message("  Skipped ", length(failed_prefixes), " unreadable IDAT pair(s): ",
              paste(head(failed_prefixes, 5), collapse = ", "),
              if (length(failed_prefixes) > 5) ", ..." else "")
    }

    if (length(beta_list) > 0) {
      common_probes <- Reduce(intersect, lapply(beta_list, rownames))
      beta <- do.call(cbind, lapply(beta_list, function(x) {
        x[common_probes, , drop = FALSE]
      }))
    }
  }

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
  
  size_mb <- file.size(checkpoint_path) / 1024^2
  progress_message(
    "Reading checkpoint: ", checkpoint_path,
    " (", round(size_mb, 1), " MB)"
  )
  t_start <- Sys.time()

  tryCatch({
    out <- qs2::qs_read(checkpoint_path)
    progress_message(
      "Finished checkpoint read: ", checkpoint_path,
      " (", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1),
      " sec)"
    )
    out
  }, error = function(e) {
    warning("\n x Corrupted checkpoint at: ", checkpoint_path,
            "\n  Error: ", conditionMessage(e),
            "\n  Leaving checkpoint in place and returning NULL.")
    
    return(NULL)
  })
}

safe_checkpoint_save <- function(object, checkpoint_path) {
  
  progress_message("Saving checkpoint: ", checkpoint_path)
  t_start <- Sys.time()

  tryCatch({
    qs2::qs_save(object, checkpoint_path)
    
    # Verify write succeeded
    if (!file.exists(checkpoint_path)) {
      stop("Checkpoint file not created at: ", checkpoint_path,
           " (Check file permissions)")
    }

    size_mb <- file.size(checkpoint_path) / 1024^2
    progress_message(
      "Finished checkpoint save: ", checkpoint_path,
      " (", round(size_mb, 1), " MB, ",
      round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1),
      " sec)"
    )
    
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
