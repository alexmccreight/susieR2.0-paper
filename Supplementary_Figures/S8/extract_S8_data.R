#!/usr/bin/env Rscript

# Supplementary Figure S8: 95% CS power and FDR in the sparse simulations with
# and without the correlation-based credible set extension, summarised per
# method and number of causal variants (K).
#
# Both versions use 95% coverage and purity >= 0.5 from the in-sample
# correlations, and differ only in the extension:
#   xcorr_noext   no extension                   series "<method>"
#   xcorr_ext     extension at correlation 0.99  series "<method> w/ Extension"
# The per-replicate metrics come precomputed in the *_ldext.rds files; a
# replicate with no credible set contributes FDR = 0.
#
# Input:  <benchmark>/final_scripts/500_rep_results/sparse/
#           sim_nrep500_sparse_h2persnp0.03_K{1-5}_n1000_ldext.rds
# Output: data/s8_finemapping_data.rds

suppressPackageStartupMessages(library(dplyr))

source("R/paths.R")

input_dir   <- file.path(benchmark_root(), "final_scripts", "500_rep_results", "sparse")
output_file <- file.path(supp_dir(8), "data", "s8_finemapping_data.rds")

h2_per_snp <- 0.03

series_levels <- c("SuSiE", "SuSiE w/ Extension",
                   "SuSiE-ash", "SuSiE-ash w/ Extension",
                   "SuSiE-inf", "SuSiE-inf w/ Extension")

summarise_k <- function(K) {
  sim <- readRDS(file.path(
    input_dir, sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000_ldext.rds", K)
  ))

  bind_rows(lapply(sim$replicates, function(rep) rep$metrics)) %>%
    filter(version %in% c("xcorr_noext", "xcorr_ext"),
           Method %in% c("SuSiE", "SuSiE-inf", "SuSiE.ash")) %>%
    group_by(Method, version) %>%
    summarise(
      Power_mean = mean(Power, na.rm = TRUE),
      Power_se   = sd(Power, na.rm = TRUE) / sqrt(n()),
      FDR_mean   = mean(FDR, na.rm = TRUE),
      FDR_se     = sd(FDR, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(K = K, Total_PVE = h2_per_snp * K)
}

all_data <- bind_rows(lapply(1:5, summarise_k))

all_data$Method <- gsub("SuSiE\\.ash", "SuSiE-ash", all_data$Method)
all_data$series <- ifelse(all_data$version == "xcorr_ext",
                          paste0(all_data$Method, " w/ Extension"),
                          all_data$Method)
all_data$series <- factor(all_data$series, levels = series_levels)
all_data$Method <- factor(all_data$Method, levels = c("SuSiE", "SuSiE-inf", "SuSiE-ash"))
all_data <- all_data[order(all_data$K, all_data$series), ]

saveRDS(all_data, output_file)
cat("Saved:", output_file, "\n")
