#!/usr/bin/env Rscript

# Figure 1, panels C-F: 95% CS power, FDR and size in the two complex
# simulation scenarios, summarised per method and number N of top variants
# treated as causal.
#   complex S1   oligogenic effects on a moderate infinitesimal background
#   complex S2   oligogenic effects on an extensive infinitesimal background
# SuSiE-ash is summarised too: Supplementary Figure S6 draws the same
# summaries with all three methods.
#
# Uses the un-extended credible sets (version "xcorr_noext": 95% coverage,
# purity >= 0.5 from the in-sample correlations, no correlation extension).
# For each replicate and each N in {3, 5, 10, 15, 20, 25}:
#   causal set   the N variants with the largest |beta|, with N capped at the
#                number of nonzero effects
#   Power        fraction of the causal set captured by at least one CS
#   FDR          fraction of credible sets that contain no variant of the
#                causal set
#   Size         mean number of variants per credible set
# A replicate with no credible set contributes Power = 0, FDR = 0 and
# Size = 0. The summary is the mean and standard error over replicates.
#
# Input:  <benchmark>/final_scripts/500_rep_results/complex_S{1,2}/sim_*_ldext.rds
# Output: data/complex_S{1,2}_finemapping_summary.rds

source("R/paths.R")
source("R/topn_panels.R")

ldext_root <- file.path(benchmark_root(), "final_scripts", "500_rep_results")
fig1_data  <- file.path(fig_dir(1), "data")

VERSION      <- "xcorr_noext"
N_values     <- c(3, 5, 10, 15, 20, 25)
method_names <- c("susie" = "SuSiE", "susie_inf" = "SuSiE-inf", "susie_ash" = "SuSiE-ash")

scenarios <- list(
  list(out  = "complex_S1_finemapping_summary.rds",
       file = file.path(ldext_root, "complex_S1",
         "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInfall_n1000_ldext.rds")),
  list(out  = "complex_S2_finemapping_summary.rds",
       file = file.path(ldext_root, "complex_S2",
         "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.15_pInf0.35_nOligo10_nInfall_n1000_ldext.rds"))
)

# Mean and standard error over replicates, per method and N.
aggregate_summary <- function(rep_data) {
  out <- do.call(rbind, lapply(split(rep_data, list(rep_data$Method, rep_data$N)), function(df) {
    data.frame(
      Method     = df$Method[1],
      K          = df$K[1],
      N          = df$N[1],
      Power_mean = mean(df$Power, na.rm = TRUE),
      Power_se   = sd(df$Power, na.rm = TRUE) / sqrt(sum(!is.na(df$Power))),
      FDR_mean   = mean(df$FDR, na.rm = TRUE),
      FDR_se     = sd(df$FDR, na.rm = TRUE) / sqrt(sum(!is.na(df$FDR))),
      Size_mean  = mean(df$Size, na.rm = TRUE),
      Size_se    = sd(df$Size, na.rm = TRUE) / sqrt(sum(!is.na(df$Size)))
    )
  }))
  rownames(out) <- NULL
  out
}

for (sc in scenarios) {
  if (!file.exists(sc$file)) stop("Missing ldext file: ", sc$file)

  sim <- readRDS(sc$file)
  K   <- sim$parameters$K

  all_rows <- list()
  for (rep in sim$replicates) {
    if (is.null(rep) || !identical(rep$status, "ok")) next

    beta       <- rep$beta
    ranked_idx <- order(abs(beta), decreasing = TRUE)
    n_nonzero  <- sum(beta != 0)

    for (N in N_values) {
      causal_set <- ranked_idx[seq_len(min(N, n_nonzero))]

      for (mkey in names(method_names)) {
        if (is.null(rep$credible_sets[[mkey]])) next   # method had no fit
        cs_list <- rep$credible_sets[[mkey]][[VERSION]]$cs

        all_rows[[length(all_rows) + 1]] <- data.frame(
          Method = method_names[[mkey]], K = K, N = N,
          topn_cs_metrics(cs_list, causal_set)
        )
      }
    }
  }

  out_path <- file.path(fig1_data, sc$out)
  saveRDS(aggregate_summary(do.call(rbind, all_rows)), out_path)
  cat("Saved:", out_path, "\n")
}
