#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S1: Extract sparse fine-mapping data
# =============================================================================
# Reads raw per-K RDS files (500 replicates each), extracts Power/FDR metrics,
# and saves a summary RDS for fast plotting.
#
# Input:  sim_nrep500_sparse_h2persnp0.03_K{1-5}_n1000.rds
# Output: data/finemapping_data.rds
# =============================================================================

library(dplyr)

# Paths
s1_dir     <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/supplemental/S1"
data_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/500_rep_results/sparse"
output_dir <- file.path(s1_dir, "data")

# Parameters
h2_value <- 0.03
K_values <- 1:5

# -----------------------------------------------------------------------------
# Extract metrics from a single RDS file
# -----------------------------------------------------------------------------
extract_metrics <- function(file_path, h2_per_snp, K) {
  data <- readRDS(file_path)

  metrics_list <- lapply(data$replicates, function(rep) rep$metrics)
  all_metrics <- bind_rows(metrics_list)

  # Keep only the three main methods
  all_metrics <- all_metrics %>%
    filter(Method %in% c("SuSiE", "SuSiE-inf", "SuSiE.ash"))

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
cat(sprintf("Extracting sparse data (h2_per_snp = %.2f, 500 reps)\n", h2_value))

all_data <- data.frame()

for (K in K_values) {
  file_name <- sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000.rds", K)
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
