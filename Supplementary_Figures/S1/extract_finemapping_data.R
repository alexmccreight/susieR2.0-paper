#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S1: Extract sparse fine-mapping data (WITH LD extension = 0.99)
# =============================================================================
# Reads the per-K *_ldext.rds files (500 replicates each), takes the credible-set
# Power/FDR metrics for the LD-extension construction (version = "xcorr_ext"),
# and saves a summary RDS for fast plotting (panels S1A / S1B).
#
# Only S1A/S1B (credible-set level) use this file. The PIP-level panels
# (S1C-S1J) are produced by extract_pip_roc.R / extract_pip_calibration.R and
# are unaffected by LD extension.
#
# Input:  500_rep_results/sparse/sim_nrep500_sparse_h2persnp0.03_K{1-5}_n1000_ldext.rds
# Output: data/finemapping_data.rds   (in this production S1 folder)
# =============================================================================

library(dplyr)

# Paths
s1_dir     <- "/Users/alexmccreight/StatFunGen/susieR2.0-paper/Supplementary_Figures/S1"
data_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/500_rep_results/sparse"
output_dir <- file.path(s1_dir, "data")

# Parameters
h2_value <- 0.03
K_values <- 1:5
VERSION  <- "xcorr_ext"   # LD-extension (0.99) credible sets

# -----------------------------------------------------------------------------
# Extract metrics from a single *_ldext.rds file
# -----------------------------------------------------------------------------
extract_metrics <- function(file_path, h2_per_snp, K) {
  data <- readRDS(file_path)

  metrics_list <- lapply(data$replicates, function(rep) rep$metrics)
  all_metrics  <- bind_rows(metrics_list)

  # Keep only the LD-extension version and the three main methods
  all_metrics <- all_metrics %>%
    filter(version == VERSION,
           Method %in% c("SuSiE", "SuSiE-inf", "SuSiE.ash"))

  summary_metrics <- all_metrics %>%
    group_by(Method) %>%
    summarise(
      Power_mean = mean(Power, na.rm = TRUE),
      Power_se   = sd(Power, na.rm = TRUE) / sqrt(n()),
      FDR_mean   = mean(FDR, na.rm = TRUE),
      FDR_se     = sd(FDR, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      h2_per_snp = h2_per_snp,
      K = K,
      Total_PVE = h2_per_snp * K
    )

  return(summary_metrics)
}

# -----------------------------------------------------------------------------
# Process all K values
# -----------------------------------------------------------------------------
cat(sprintf("Extracting sparse data WITH LD extension (version = %s, h2_per_snp = %.2f, 500 reps)\n",
            VERSION, h2_value))

all_data <- data.frame()

for (K in K_values) {
  file_name <- sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000_ldext.rds", K)
  file_path <- file.path(data_dir, file_name)

  if (file.exists(file_path)) {
    cat(sprintf("  Processing: %s\n", file_name))
    metrics <- extract_metrics(file_path, h2_value, K)
    all_data <- bind_rows(all_data, metrics)
  } else {
    warning(sprintf("  File not found: %s", file_name))
  }
}

# Standardize method name
all_data$Method <- gsub("SuSiE\\.ash", "SuSiE-ash", all_data$Method)
all_data$Method <- factor(all_data$Method,
                          levels = c("SuSiE", "SuSiE-inf", "SuSiE-ash"))

# Save
output_file <- file.path(output_dir, "finemapping_data.rds")
saveRDS(all_data, output_file)

cat(sprintf("\nExtraction complete!\n"))
cat(sprintf("Saved to: %s\n", output_file))
cat(sprintf("Dimensions: %d rows x %d cols\n", nrow(all_data), ncol(all_data)))
cat(sprintf("Methods: %s\n", paste(levels(all_data$Method), collapse = ", ")))
cat(sprintf("K values: %s\n", paste(sort(unique(all_data$K)), collapse = ", ")))
print(all_data[order(all_data$K, all_data$Method), c("Method","K","Power_mean","FDR_mean")])
