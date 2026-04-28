#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S3: Extract TWAS prediction (CV R²) data
# =============================================================================
# Reads prediction RDS files and extracts per-replicate mean R² and winner.
#
# Input:  prediction_cv_results/sparse/pred_cv_nrep300_sparse_*K{1-5}*.rds
#         prediction_cv_results/complex/pred_cv_nrep300_oligo_*.rds  (3 scenarios)
# Output: data/prediction_data.rds
# =============================================================================

library(dplyr)

s3_dir       <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/supplemental/S3"
sparse_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/prediction_cv_results/sparse"
complex_dir  <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/prediction_cv_results/complex"
output_dir   <- file.path(s3_dir, "data")

methods <- c("SuSiE", "SuSiE_inf", "SuSiE_ash")
method_labels <- c("SuSiE", "SuSiE-inf", "SuSiE-ash")

# -----------------------------------------------------------------------------
# Helper: extract per-replicate R² from one RDS file
# -----------------------------------------------------------------------------
extract_r2 <- function(file_path) {
  d <- readRDS(file_path)

  rows <- lapply(seq_along(d$replicates), function(i) {
    rep <- d$replicates[[i]]
    r2_vals <- sapply(methods, function(m) rep$cv_metrics[[m]]$mean_r2)
    data.frame(
      replicate = i,
      Method = method_labels,
      R2 = r2_vals,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(rows)
}

# Helper: compute proportion of replicates won per method
compute_prop <- function(data, group_vars = NULL) {
  if (is.null(group_vars)) {
    winners <- data %>%
      group_by(replicate) %>%
      slice_max(R2, n = 1, with_ties = FALSE) %>%
      ungroup()
  } else {
    winners <- data %>%
      group_by(across(all_of(c(group_vars, "replicate")))) %>%
      slice_max(R2, n = 1, with_ties = FALSE) %>%
      ungroup()
  }
  winners %>%
    count(Method) %>%
    mutate(Proportion = n / sum(n))
}

# -----------------------------------------------------------------------------
# Sparse: K = 1, 2, 3, 4, 5
# -----------------------------------------------------------------------------
cat("Extracting sparse prediction data (K = 1-5)...\n")

sparse_list <- list()
for (K in 1:5) {
  fname <- sprintf("pred_cv_nrep300_sparse_h2persnp0.03_K%d_n1000_nfolds5.rds", K)
  fpath <- file.path(sparse_dir, fname)
  cat(sprintf("  Processing: %s\n", fname))
  df <- extract_r2(fpath)
  df$K <- K
  df$Total_PVE <- 0.03 * K
  sparse_list[[K]] <- df
}
sparse_data <- bind_rows(sparse_list)
sparse_prop <- compute_prop(sparse_data, group_vars = "K")

# -----------------------------------------------------------------------------
# Complex 1: Oligogenic on a Polygenic Background (nInf = 15 SNPs)
# -----------------------------------------------------------------------------
cat("Extracting Complex 1 prediction data (polygenic background)...\n")

complex1_file <- file.path(complex_dir,
  "pred_cv_nrep300_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInf15_n1000_nfolds5.rds")
complex1_data <- extract_r2(complex1_file)
complex1_prop <- compute_prop(complex1_data)

# -----------------------------------------------------------------------------
# Complex 2: Oligogenic on a Moderate Infinitesimal Background
# -----------------------------------------------------------------------------
cat("Extracting Complex 2 prediction data (moderate infinitesimal)...\n")

complex2_file <- file.path(complex_dir,
  "pred_cv_nrep300_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInfall_n1000_nfolds5.rds")
complex2_data <- extract_r2(complex2_file)
complex2_prop <- compute_prop(complex2_data)

# -----------------------------------------------------------------------------
# Complex 3: Oligogenic on an Extensive Infinitesimal Background
# -----------------------------------------------------------------------------
cat("Extracting Complex 3 prediction data (extensive infinitesimal)...\n")

complex3_file <- file.path(complex_dir,
  "pred_cv_nrep300_oligo_h2g0.25_K3_pSparse0.5_pOligo0.15_pInf0.35_nOligo10_nInfall_n1000_nfolds5.rds")
complex3_data <- extract_r2(complex3_file)
complex3_prop <- compute_prop(complex3_data)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
out <- list(
  sparse_data  = sparse_data,
  sparse_prop  = sparse_prop,
  complex1_data = complex1_data,
  complex1_prop = complex1_prop,
  complex2_data = complex2_data,
  complex2_prop = complex2_prop,
  complex3_data = complex3_data,
  complex3_prop = complex3_prop
)

output_file <- file.path(output_dir, "prediction_data.rds")
saveRDS(out, output_file)

cat(sprintf("\nExtraction complete! Saved to: %s\n", output_file))
cat("\n== Sparse proportion of replicates won (pooled K=1-5) ==\n")
print(as.data.frame(sparse_prop))
cat("\n== Complex 1 proportion of replicates won ==\n")
print(as.data.frame(complex1_prop))
cat("\n== Complex 2 proportion of replicates won ==\n")
print(as.data.frame(complex2_prop))
cat("\n== Complex 3 proportion of replicates won ==\n")
print(as.data.frame(complex3_prop))
