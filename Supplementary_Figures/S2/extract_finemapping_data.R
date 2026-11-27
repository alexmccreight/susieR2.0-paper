#!/usr/bin/env Rscript

# Supplementary Figure S2, panels A and B: 95% CS power and FDR in the sparse
# simulations, summarised per method and number of causal variants (K).
#
# Uses the un-extended credible sets (version "xcorr_noext": 95% coverage,
# purity >= 0.5 from the in-sample correlations, no correlation extension).
# The per-replicate metrics come precomputed in the *_ldext.rds files; a
# replicate with no credible set contributes FDR = 0.
#
# Input:  <benchmark>/final_scripts/500_rep_results/sparse/
#           sim_nrep500_sparse_h2persnp0.03_K{1-5}_n1000_ldext.rds
# Output: data/finemapping_data.rds

suppressPackageStartupMessages(library(dplyr))

source("R/paths.R")

input_dir   <- file.path(benchmark_root(), "final_scripts", "500_rep_results", "sparse")
output_file <- file.path(supp_dir(2), "data", "finemapping_data.rds")

h2_per_snp <- 0.03

summarise_k <- function(K) {
  sim <- readRDS(file.path(
    input_dir, sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000_ldext.rds", K)
  ))

  bind_rows(lapply(sim$replicates, function(rep) rep$metrics)) %>%
    filter(version == "xcorr_noext",
           Method %in% c("SuSiE", "SuSiE-inf", "SuSiE.ash")) %>%
    group_by(Method) %>%
    summarise(
      Power_mean = mean(Power, na.rm = TRUE),
      Power_se   = sd(Power, na.rm = TRUE) / sqrt(n()),
      FDR_mean   = mean(FDR, na.rm = TRUE),
      FDR_se     = sd(FDR, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(h2_per_snp = h2_per_snp, K = K, Total_PVE = h2_per_snp * K)
}

all_data <- as.data.frame(bind_rows(lapply(1:5, summarise_k)))
all_data$Method <- factor(sub("SuSiE.ash", "SuSiE-ash", all_data$Method, fixed = TRUE),
                          levels = c("SuSiE", "SuSiE-inf", "SuSiE-ash"))

saveRDS(all_data, output_file)
cat("Saved:", output_file, "\n")
