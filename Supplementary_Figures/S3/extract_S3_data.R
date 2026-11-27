#!/usr/bin/env Rscript

# Supplementary Figure S3: one row per 95% credible set, for every method,
# replicate and simulation scenario.
#
# Uses the un-extended credible sets (version "xcorr_noext": 95% coverage,
# purity >= 0.5 from the in-sample correlations, no correlation extension).
# cs_size is the number of variants in the set and min_abs_corr its purity,
# the smallest absolute correlation between any two of its variants.
# Replicates whose status is not "ok" are skipped.
#
# Scenario labels, as stored in the output and used by plot_S3.R:
#   Sparse      K = 1-5 causal variants, pooled            (panel A)
#   Complex     oligogenic effects, polygenic background   (panel B)
#   Complex S1  moderate infinitesimal background          (panel C)
#   Complex S2  extensive infinitesimal background         (panel D)
#
# Input:  <benchmark>/final_scripts/500_rep_results/
#           sparse/sim_nrep500_sparse_h2persnp0.03_K{1-5}_n1000_ldext.rds
#           complex/sim_nrep500_oligo_*_nOligo5_nInf15_n1000_ldext.rds
#           complex_S1/sim_nrep500_oligo_*_nOligo5_nInfall_n1000_ldext.rds
#           complex_S2/sim_nrep500_oligo_*_nOligo10_nInfall_n1000_ldext.rds
# Output: data/s3_cs_data.rds  (list: $cs_data, $meta)

source("R/paths.R")

rep_results_dir <- file.path(benchmark_root(), "final_scripts", "500_rep_results")
output_file     <- file.path(supp_dir(3), "data", "s3_cs_data.rds")

VERSION <- "xcorr_noext"

scenarios <- list(
  list(
    label = "Sparse",
    files = file.path(rep_results_dir, "sparse",
                      sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000_ldext.rds", 1:5))
  ),
  list(
    label = "Complex",
    files = file.path(rep_results_dir, "complex",
      "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInf15_n1000_ldext.rds")
  ),
  list(
    label = "Complex S1",
    files = file.path(rep_results_dir, "complex_S1",
      "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInfall_n1000_ldext.rds")
  ),
  list(
    label = "Complex S2",
    files = file.path(rep_results_dir, "complex_S2",
      "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.15_pInf0.35_nOligo10_nInfall_n1000_ldext.rds")
  )
)

method_names <- c(susie = "SuSiE", susie_ash = "SuSiE-ash", susie_inf = "SuSiE-inf")

# ---- One row per credible set ----

extract_cs_rows <- function(rds_path, scenario) {
  rows <- list()
  for (rep in readRDS(rds_path)$replicates) {
    if (is.null(rep) || !identical(rep$status, "ok")) next
    for (fit_name in names(method_names)) {
      cs_rec <- rep$credible_sets[[fit_name]][[VERSION]]
      if (length(cs_rec$cs) == 0) next
      rows[[length(rows) + 1]] <- data.frame(
        replicate    = rep$replicate_id,
        method       = method_names[[fit_name]],
        cs_size      = lengths(cs_rec$cs, use.names = FALSE),
        min_abs_corr = unname(cs_rec$purity),
        scenario     = scenario
      )
    }
  }
  do.call(rbind, rows)
}

all_cs <- do.call(rbind, lapply(scenarios, function(sc) {
  do.call(rbind, lapply(sc$files, extract_cs_rows, scenario = sc$label))
}))

# ---- Factor levels and size bins ----

method_levels <- c("SuSiE", "SuSiE-ash", "SuSiE-inf")
all_cs$method <- factor(all_cs$method, levels = method_levels)

scenario_levels <- c("Sparse", "Complex", "Complex S1", "Complex S2")
all_cs$scenario <- factor(all_cs$scenario, levels = scenario_levels)

size_bin_levels <- c("1", "2", "3-5", "6-10", ">10")
all_cs$cs_size_bin <- cut(
  all_cs$cs_size,
  breaks = c(0, 1, 2, 5, 10, Inf),
  labels = size_bin_levels
)

# ---- Save ----

out <- list(
  cs_data = all_cs,
  meta = list(
    method_levels   = method_levels,
    scenario_levels = scenario_levels,
    size_bin_levels = size_bin_levels,
    ld_extend       = NA,      # no correlation extension in this construction
    version         = VERSION
  )
)

saveRDS(out, output_file)
cat("Saved:", output_file, "\n")
