# =============================================================================
# Extract BVSR Results from Fine-Mapping RDS Files
# =============================================================================
# Step 1 of the AlphaGenome scoring pipeline. Walks the per-gene fine-mapping
# *.rds outputs (raw univariate_bvsr-style nested lists) for each of the
# three SuSiE methods (standard, ash, inf) and flattens them into tables.
#
# For each method, four files are written to output_dir:
#   <method>_top_loci.rds              per-gene-tissue top loci with PIPs
#   <method>_credible_sets.rds         CS-level summary (one row per CS)
#   <method>_gene_tissue_summary.rds   per gene-tissue stats (n_cs, etc.)
#   <method>_variant_cs_membership.rds variant-level CS membership
#
# The middle two (`<method>_credible_sets.rds` and
# `<method>_gene_tissue_summary.rds`) — 6 files across the 3 methods — are the
# Figure 2 panel A inputs and the inputs to alphagenome_cs_scoring.py.
# The other two are byproducts.
#
# Usage:
#   Rscript extract_bvsr_results.R                                # all 3 methods
#   Rscript extract_bvsr_results.R method='susie_standard'        # one method
#   Rscript extract_bvsr_results.R input_dir='/path' output_dir='/path'
# =============================================================================

library(dplyr)
library(tidyr)


# =============================================================================
# Helper functions: Variant ID + tissue parsing
# =============================================================================

#' Parse variant ID into components
#' @param variant_id Character vector of variant IDs in chr:pos:ref:alt format
#' @return Data frame with chr, pos, ref, alt columns
parse_variant_id <- function(variant_id) {
  parts <- strsplit(variant_id, ":")
  data.frame(
    chr = sapply(parts, `[`, 1),
    pos = as.integer(sapply(parts, `[`, 2)),
    ref = sapply(parts, `[`, 3),
    alt = sapply(parts, `[`, 4),
    stringsAsFactors = FALSE
  )
}

#' Extract tissue name from nested list key
#' @param tissue_key Character, e.g., "AC_DeJager_eQTL_ENSG00000047056"
#' @return Character, e.g., "AC"
extract_tissue_name <- function(tissue_key) {
  # Extract first part before underscore (AC, DLPFC, PCC)
  strsplit(tissue_key, "_")[[1]][1]
}

# =============================================================================
# Helper functions: Per-table extractors (one helper per output table)
# =============================================================================

#' Extract top_loci table from a single gene's data
#' @param gene_data List containing tissue data for one gene
#' @param gene_id Character, Ensembl gene ID
#' @return Data frame with top_loci for all tissues
extract_top_loci <- function(gene_data, gene_id) {
  results <- list()

  for (tissue_key in names(gene_data)) {
    tissue_data <- gene_data[[tissue_key]]
    tissue_name <- extract_tissue_name(tissue_key)

    # Check if top_loci exists
    if (is.null(tissue_data$top_loci) || nrow(tissue_data$top_loci) == 0) {
      next
    }

    top_loci <- tissue_data$top_loci

    # Parse variant IDs
    variant_parts <- parse_variant_id(top_loci$variant_id)

    # Build result data frame
    result <- data.frame(
      gene_id = gene_id,
      tissue = tissue_name,
      variant_id = top_loci$variant_id,
      chr = variant_parts$chr,
      pos = variant_parts$pos,
      ref = variant_parts$ref,
      alt = variant_parts$alt,
      betahat = top_loci$betahat,
      sebetahat = top_loci$sebetahat,
      z = top_loci$z,
      maf = top_loci$maf,
      pip = top_loci$pip,
      stringsAsFactors = FALSE
    )

    # Add CS membership columns (handle different naming conventions)
    if ("cs_coverage_0.95" %in% names(top_loci)) {
      result$cs_id_0.95 <- top_loci$`cs_coverage_0.95`
    } else if ("cs_coverage_0_95" %in% names(top_loci)) {
      result$cs_id_0.95 <- top_loci$cs_coverage_0_95
    } else {
      result$cs_id_0.95 <- NA_integer_
    }

    if ("cs_coverage_0.7" %in% names(top_loci)) {
      result$cs_id_0.70 <- top_loci$`cs_coverage_0.7`
    } else if ("cs_coverage_0_7" %in% names(top_loci)) {
      result$cs_id_0.70 <- top_loci$cs_coverage_0_7
    } else {
      result$cs_id_0.70 <- NA_integer_
    }

    if ("cs_coverage_0.5" %in% names(top_loci)) {
      result$cs_id_0.50 <- top_loci$`cs_coverage_0.5`
    } else if ("cs_coverage_0_5" %in% names(top_loci)) {
      result$cs_id_0.50 <- top_loci$cs_coverage_0_5
    } else {
      result$cs_id_0.50 <- NA_integer_
    }

    results[[tissue_key]] <- result
  }

  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

#' Extract credible set summary from a single gene's data
#' @param gene_data List containing tissue data for one gene
#' @param gene_id Character, Ensembl gene ID
#' @return Data frame with credible set summaries
extract_credible_sets <- function(gene_data, gene_id) {
  results <- list()

  for (tissue_key in names(gene_data)) {
    tissue_data <- gene_data[[tissue_key]]
    tissue_name <- extract_tissue_name(tissue_key)

    # Check if susie_result_trimmed exists
    if (is.null(tissue_data$susie_result_trimmed)) {
      next
    }

    susie_result <- tissue_data$susie_result_trimmed
    variant_names <- tissue_data$variant_names

    # Process primary credible sets (0.95 coverage)
    if (!is.null(susie_result$sets) && !is.null(susie_result$sets$cs)) {
      cs_list <- susie_result$sets$cs
      purity <- susie_result$sets$purity
      coverage <- susie_result$sets$coverage

      for (i in seq_along(cs_list)) {
        cs_name <- names(cs_list)[i]
        cs_indices <- cs_list[[i]]
        cs_variants <- variant_names[cs_indices]

        # Get PIPs for variants in this CS
        cs_pips <- susie_result$pip[cs_indices]
        top_idx <- which.max(cs_pips)

        result <- data.frame(
          gene_id = gene_id,
          tissue = tissue_name,
          cs_id = as.integer(gsub("L", "", cs_name)),
          coverage_level = "0.95",
          n_variants = length(cs_indices),
          coverage = if (!is.null(coverage)) coverage[i] else NA_real_,
          min_abs_corr = if (!is.null(purity)) purity$min.abs.corr[i] else NA_real_,
          mean_abs_corr = if (!is.null(purity)) purity$mean.abs.corr[i] else NA_real_,
          median_abs_corr = if (!is.null(purity)) purity$median.abs.corr[i] else NA_real_,
          top_variant = cs_variants[top_idx],
          top_pip = cs_pips[top_idx],
          variants = paste(cs_variants, collapse = ","),
          stringsAsFactors = FALSE
        )

        results[[paste(tissue_key, cs_name, "0.95", sep = "_")]] <- result
      }
    }

    # Process secondary credible sets (0.7 and 0.5 coverage)
    if (!is.null(susie_result$sets_secondary)) {
      for (cov_level in names(susie_result$sets_secondary)) {
        cov_label <- gsub("coverage_", "", cov_level)
        secondary <- susie_result$sets_secondary[[cov_level]]

        if (!is.null(secondary$sets) && !is.null(secondary$sets$cs)) {
          cs_list <- secondary$sets$cs
          purity <- secondary$sets$purity
          coverage <- secondary$sets$coverage

          for (i in seq_along(cs_list)) {
            cs_name <- names(cs_list)[i]
            cs_indices <- cs_list[[i]]
            cs_variants <- variant_names[cs_indices]

            cs_pips <- susie_result$pip[cs_indices]
            top_idx <- which.max(cs_pips)

            result <- data.frame(
              gene_id = gene_id,
              tissue = tissue_name,
              cs_id = as.integer(gsub("L", "", cs_name)),
              coverage_level = cov_label,
              n_variants = length(cs_indices),
              coverage = if (!is.null(coverage)) coverage[i] else NA_real_,
              min_abs_corr = if (!is.null(purity)) purity$min.abs.corr[i] else NA_real_,
              mean_abs_corr = if (!is.null(purity)) purity$mean.abs.corr[i] else NA_real_,
              median_abs_corr = if (!is.null(purity)) purity$median.abs.corr[i] else NA_real_,
              top_variant = cs_variants[top_idx],
              top_pip = cs_pips[top_idx],
              variants = paste(cs_variants, collapse = ","),
              stringsAsFactors = FALSE
            )

            results[[paste(tissue_key, cs_name, cov_label, sep = "_")]] <- result
          }
        }
      }
    }
  }

  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

#' Extract gene-tissue summary from a single gene's data
#' @param gene_data List containing tissue data for one gene
#' @param gene_id Character, Ensembl gene ID
#' @return Data frame with gene-tissue summaries
extract_gene_summary <- function(gene_data, gene_id) {
  results <- list()

  for (tissue_key in names(gene_data)) {
    tissue_data <- gene_data[[tissue_key]]
    tissue_name <- extract_tissue_name(tissue_key)

    # Count credible sets at different coverage levels
    n_cs_95 <- 0
    n_cs_70 <- 0
    n_cs_50 <- 0

    has_susie <- !is.null(tissue_data$susie_result_trimmed)

    if (has_susie) {
      susie_result <- tissue_data$susie_result_trimmed

      if (!is.null(susie_result$sets) && !is.null(susie_result$sets$cs)) {
        n_cs_95 <- length(susie_result$sets$cs)
      }

      if (!is.null(susie_result$sets_secondary)) {
        if (!is.null(susie_result$sets_secondary$coverage_0.7$sets$cs)) {
          n_cs_70 <- length(susie_result$sets_secondary$coverage_0.7$sets$cs)
        }
        if (!is.null(susie_result$sets_secondary$coverage_0.5$sets$cs)) {
          n_cs_50 <- length(susie_result$sets_secondary$coverage_0.5$sets$cs)
        }
      }
    }

    # Extract region info
    region_chrom <- NA_character_
    region_start <- NA_integer_
    region_end <- NA_integer_

    if (!is.null(tissue_data$region_info) && !is.null(tissue_data$region_info$grange)) {
      region_chrom <- tissue_data$region_info$grange$chrom[1]
      region_start <- tissue_data$region_info$grange$start[1]
      region_end <- tissue_data$region_info$grange$end[1]
    }

    # Extract elapsed time
    elapsed_time <- NA_real_
    if (!is.null(tissue_data$total_time_elapsed)) {
      elapsed_time <- tissue_data$total_time_elapsed["elapsed"]
    }

    result <- data.frame(
      gene_id = gene_id,
      tissue = tissue_name,
      n_variants = length(tissue_data$variant_names),
      n_samples = length(tissue_data$sample_names),
      has_susie_result = has_susie,
      n_credible_sets_0.95 = n_cs_95,
      n_credible_sets_0.70 = n_cs_70,
      n_credible_sets_0.50 = n_cs_50,
      n_top_loci = if (!is.null(tissue_data$top_loci)) nrow(tissue_data$top_loci) else 0L,
      region_chrom = region_chrom,
      region_start = region_start,
      region_end = region_end,
      elapsed_time = as.numeric(elapsed_time),
      stringsAsFactors = FALSE
    )

    results[[tissue_key]] <- result
  }

  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

#' Extract variant CS membership details from a single gene's data
#' @param gene_data List containing tissue data for one gene
#' @param gene_id Character, Ensembl gene ID
#' @return Data frame with variant-level CS membership
extract_cs_variants <- function(gene_data, gene_id) {
  results <- list()

  for (tissue_key in names(gene_data)) {
    tissue_data <- gene_data[[tissue_key]]
    tissue_name <- extract_tissue_name(tissue_key)

    if (is.null(tissue_data$susie_result_trimmed)) {
      next
    }

    susie_result <- tissue_data$susie_result_trimmed
    variant_names <- tissue_data$variant_names

    # Helper function to extract CS variants
    extract_cs_level <- function(cs_list, pips, alpha_matrix, coverage_level) {
      level_results <- list()

      for (i in seq_along(cs_list)) {
        cs_name <- names(cs_list)[i]
        cs_id <- as.integer(gsub("L", "", cs_name))
        cs_indices <- cs_list[[i]]

        # Get alpha values for this CS (row corresponds to L number)
        alpha_row <- NULL
        if (!is.null(alpha_matrix) && nrow(alpha_matrix) >= cs_id) {
          alpha_row <- alpha_matrix[cs_id, cs_indices]
        }

        for (j in seq_along(cs_indices)) {
          idx <- cs_indices[j]

          level_results[[paste(cs_name, idx, sep = "_")]] <- data.frame(
            gene_id = gene_id,
            tissue = tissue_name,
            cs_id = cs_id,
            coverage_level = coverage_level,
            variant_id = variant_names[idx],
            variant_index = idx,
            pip = pips[idx],
            alpha = if (!is.null(alpha_row)) alpha_row[j] else NA_real_,
            stringsAsFactors = FALSE
          )
        }
      }

      level_results
    }

    # Primary CS (0.95)
    if (!is.null(susie_result$sets) && !is.null(susie_result$sets$cs)) {
      primary_results <- extract_cs_level(
        susie_result$sets$cs,
        susie_result$pip,
        susie_result$alpha,
        "0.95"
      )
      results <- c(results, primary_results)
    }

    # Secondary CS (0.7, 0.5)
    if (!is.null(susie_result$sets_secondary)) {
      for (cov_level in names(susie_result$sets_secondary)) {
        cov_label <- gsub("coverage_", "", cov_level)
        secondary <- susie_result$sets_secondary[[cov_level]]

        if (!is.null(secondary$sets) && !is.null(secondary$sets$cs)) {
          secondary_results <- extract_cs_level(
            secondary$sets$cs,
            susie_result$pip,
            susie_result$alpha,
            cov_label
          )
          results <- c(results, secondary_results)
        }
      }
    }
  }

  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

# =============================================================================
# Helper functions: Single-file dispatcher
# =============================================================================

#' Process a single RDS file and extract all tables
#' @param file_path Path to RDS file
#' @return List with 4 data frames (top_loci, credible_sets, gene_summary, cs_variants)
process_single_file <- function(file_path) {
  data <- readRDS(file_path)

  # Handle nested structure - data is a list with gene_id as key
  gene_id <- names(data)[1]
  gene_data <- data[[gene_id]]

  list(
    top_loci = extract_top_loci(gene_data, gene_id),
    credible_sets = extract_credible_sets(gene_data, gene_id),
    gene_summary = extract_gene_summary(gene_data, gene_id),
    cs_variants = extract_cs_variants(gene_data, gene_id)
  )
}

# =============================================================================
# Main Processing Function
# =============================================================================

#' Main extraction function
#' @param input_dir Path to mounted S3 directory containing method folders
#' @param output_dir Path to output directory
#' @param method Which method to process ("all", "susie_standard", "susie_ash", "susie_inf")
extract_bvsr_results <- function(input_dir = NULL,
                                  output_dir = NULL,
                                  method = NULL) {

  # Set defaults
  # Mount: statfungen/ftp_fgc_xqtl/interactive_sessions/apm2217 -> /home/apm2217/data
  if (is.null(input_dir)) input_dir <- "/home/apm2217/data"
  if (is.null(output_dir)) output_dir <- "/home/apm2217/output"
  if (is.null(method)) method <- "all"

  # Parse command-line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    for (arg in args) {
      split_arg <- strsplit(arg, "=")[[1]]
      if (length(split_arg) == 2) {
        key <- split_arg[1]
        value <- gsub("^['\"]|['\"]$", "", split_arg[2])

        if (key == "input_dir") {
          input_dir <- value
        } else if (key == "output_dir") {
          output_dir <- value
        } else if (key == "method") {
          method <- value
        }
      }
    }
  }

  # Validate inputs
  if (is.null(input_dir)) {
    stop("input_dir must be provided")
  }

  if (!dir.exists(input_dir)) {
    stop("input_dir does not exist: ", input_dir)
  }

  # Define method paths (relative to input_dir, which should be the mount point)
  # Mount: statfungen/ftp_fgc_xqtl/interactive_sessions/apm2217 -> /home/apm2217/data
  # So input_dir should be /home/apm2217/data
  method_paths <- list(
    susie_standard = "real_data_analysis/susie_standard/fine_mapping",
    susie_ash = "real_data_analysis/updated_ash_inf_02012026/susie_ash/susie_twas/fine_mapping",
    susie_inf = "real_data_analysis/updated_ash_inf_02012026/susie_inf/susie_twas/fine_mapping"
  )

  # Determine which methods to process
  if (method == "all") {
    methods_to_process <- names(method_paths)
  } else if (method %in% names(method_paths)) {
    methods_to_process <- method
  } else {
    stop("Invalid method. Must be 'all', 'susie_standard', 'susie_ash', or 'susie_inf'")
  }

  # Process each method
  for (method_name in methods_to_process) {
    cat("\n", paste(rep("=", 60), collapse = ""), "\n")
    cat("Processing method:", method_name, "\n")
    cat(paste(rep("=", 60), collapse = ""), "\n\n")

    method_dir <- file.path(input_dir, method_paths[[method_name]])

    if (!dir.exists(method_dir)) {
      warning("Method directory does not exist: ", method_dir, ". Skipping.")
      next
    }

    # Find all RDS files
    rds_files <- list.files(method_dir, pattern = "\\.rds$", full.names = TRUE)

    if (length(rds_files) == 0) {
      warning("No RDS files found in: ", method_dir)
      next
    }

    cat("Found", length(rds_files), "RDS files\n\n")

    # Initialize result collectors
    all_top_loci <- list()
    all_credible_sets <- list()
    all_gene_summary <- list()
    all_cs_variants <- list()
    failed_files <- character()

    # Process each file
    start_time <- Sys.time()

    for (i in seq_along(rds_files)) {
      file_path <- rds_files[i]

      # Progress report every 100 files
      if (i %% 100 == 0 || i == 1) {
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        rate <- i / elapsed
        eta <- (length(rds_files) - i) / rate
        cat(sprintf("Processing file %d/%d (%.1f files/sec, ETA: %.1f min)...\n",
                    i, length(rds_files), rate, eta / 60))
      }

      # Process file with error handling
      result <- tryCatch({
        process_single_file(file_path)
      }, error = function(e) {
        warning("Failed to process: ", basename(file_path), " - ", e$message)
        failed_files <<- c(failed_files, file_path)
        NULL
      })

      if (!is.null(result)) {
        if (!is.null(result$top_loci)) {
          all_top_loci[[i]] <- result$top_loci
        }
        if (!is.null(result$credible_sets)) {
          all_credible_sets[[i]] <- result$credible_sets
        }
        if (!is.null(result$gene_summary)) {
          all_gene_summary[[i]] <- result$gene_summary
        }
        if (!is.null(result$cs_variants)) {
          all_cs_variants[[i]] <- result$cs_variants
        }
      }

      # Periodic garbage collection
      if (i %% 500 == 0) {
        gc()
      }
    }

    # Combine results
    cat("\nCombining results...\n")

    top_loci_combined <- bind_rows(all_top_loci)
    credible_sets_combined <- bind_rows(all_credible_sets)
    gene_summary_combined <- bind_rows(all_gene_summary)
    cs_variants_combined <- bind_rows(all_cs_variants)

    # Create output directory (single folder, method name in filenames)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    # Save results with method name prefix to avoid overwriting
    cat("Saving results to:", output_dir, "\n")

    saveRDS(top_loci_combined, file.path(output_dir, paste0(method_name, "_top_loci.rds")))
    cat("  -", paste0(method_name, "_top_loci.rds:"), nrow(top_loci_combined), "rows\n")

    saveRDS(credible_sets_combined, file.path(output_dir, paste0(method_name, "_credible_sets.rds")))
    cat("  -", paste0(method_name, "_credible_sets.rds:"), nrow(credible_sets_combined), "rows\n")

    saveRDS(gene_summary_combined, file.path(output_dir, paste0(method_name, "_gene_tissue_summary.rds")))
    cat("  -", paste0(method_name, "_gene_tissue_summary.rds:"), nrow(gene_summary_combined), "rows\n")

    saveRDS(cs_variants_combined, file.path(output_dir, paste0(method_name, "_variant_cs_membership.rds")))
    cat("  -", paste0(method_name, "_variant_cs_membership.rds:"), nrow(cs_variants_combined), "rows\n")

    # Save failed files log if any
    if (length(failed_files) > 0) {
      failed_log_path <- file.path(output_dir, paste0(method_name, "_failed_files.txt"))
      writeLines(failed_files, failed_log_path)
      cat("  -", paste0(method_name, "_failed_files.txt:"), length(failed_files), "files failed\n")
    }

    elapsed_total <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat("\nMethod", method_name, "completed in", round(elapsed_total, 1), "minutes\n")

    # Clean up
    rm(all_top_loci, all_credible_sets, all_gene_summary, all_cs_variants)
    rm(top_loci_combined, credible_sets_combined, gene_summary_combined, cs_variants_combined)
    gc()
  }

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("All processing complete!\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
}

# Run the extraction if called directly
if (!interactive() && !exists(".testing_mode")) {
  extract_bvsr_results()
}
