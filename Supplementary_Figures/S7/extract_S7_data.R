#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S7: Extract sparse CS-level metrics for the
#                  no-extension vs LD-extension (0.99) comparison
# =============================================================================
# Same sparse data as S1A/S1B, but keeps BOTH credible-set constructions:
#   xcorr_noext -> "SuSiE" / "SuSiE-inf" / "SuSiE-ash"               (no extension)
#   xcorr_ext   -> "... w/ Extension"                                (LD ext 0.99)
# so S7 can show 2 bars per method per K (6 bars per K).
#
# Both versions use the exact (all-pairs) purity filter; they differ ONLY by the
# LD extension, isolating the extension's effect on Power/FDR.
#
# Input:  500_rep_results/sparse/sim_nrep500_sparse_h2persnp0.03_K{1-5}_n1000_ldext.rds
# Output: data/s7_finemapping_data.rds
# =============================================================================

library(dplyr)

s7_dir     <- "/Users/alexmccreight/StatFunGen/susieR2.0-paper/Supplementary_Figures/S7"
data_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/500_rep_results/sparse"
output_dir <- file.path(s7_dir, "data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

h2_value <- 0.03
K_values <- 1:5
versions <- c("xcorr_noext", "xcorr_ext")

series_levels <- c("SuSiE", "SuSiE w/ Extension",
                   "SuSiE-ash", "SuSiE-ash w/ Extension",
                   "SuSiE-inf", "SuSiE-inf w/ Extension")

# -----------------------------------------------------------------------------
# Aggregate Power/FDR by Method x version for one *_ldext.rds file
# -----------------------------------------------------------------------------
extract_metrics <- function(file_path, K) {
  data <- readRDS(file_path)
  all_metrics <- bind_rows(lapply(data$replicates, function(rep) rep$metrics))

  all_metrics %>%
    filter(version %in% versions,
           Method %in% c("SuSiE", "SuSiE-inf", "SuSiE.ash")) %>%
    group_by(Method, version) %>%
    summarise(
      Power_mean = mean(Power, na.rm = TRUE),
      Power_se   = sd(Power, na.rm = TRUE) / sqrt(n()),
      FDR_mean   = mean(FDR, na.rm = TRUE),
      FDR_se     = sd(FDR, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(K = K, Total_PVE = h2_value * K)
}

cat("Extracting sparse no-ext vs LD-ext metrics...\n")
all_data <- bind_rows(lapply(K_values, function(K) {
  fn <- sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000_ldext.rds", K)
  fp <- file.path(data_dir, fn)
  if (!file.exists(fp)) { warning("File not found: ", fn); return(NULL) }
  cat(sprintf("  Processing: %s\n", fn))
  extract_metrics(fp, K)
}))

# Method display name + 6-level series factor (method-major, version-minor)
all_data$Method <- gsub("SuSiE\\.ash", "SuSiE-ash", all_data$Method)
all_data$series <- ifelse(all_data$version == "xcorr_ext",
                          paste0(all_data$Method, " w/ Extension"),
                          all_data$Method)
all_data$series <- factor(all_data$series, levels = series_levels)
all_data$Method <- factor(all_data$Method, levels = c("SuSiE", "SuSiE-inf", "SuSiE-ash"))
all_data <- all_data[order(all_data$K, all_data$series), ]

output_file <- file.path(output_dir, "s7_finemapping_data.rds")
saveRDS(all_data, output_file)

cat(sprintf("\nSaved: %s\n", output_file))
cat(sprintf("Dimensions: %d rows x %d cols  (expect 30 = 6 series x 5 K)\n",
            nrow(all_data), ncol(all_data)))
print(all_data[all_data$K %in% c(1, 5),
               c("series", "K", "Power_mean", "FDR_mean")])
