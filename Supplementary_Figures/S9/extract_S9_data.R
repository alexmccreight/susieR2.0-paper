#!/usr/bin/env Rscript

# Supplementary Figure S9: 95% CS power and FDR in the three complex simulation
# scenarios with and without the correlation-based credible set extension,
# summarised per method and number N of top variants treated as causal.
#   Complex      oligogenic effects on a polygenic background
#   Complex S1   oligogenic effects on a moderate infinitesimal background
#   Complex S2   oligogenic effects on an extensive infinitesimal background
#
# Both credible set versions use 95% coverage and purity >= 0.5 from the
# in-sample correlations, and differ only in the extension:
#   xcorr_noext   no extension                   series "<method>"
#   xcorr_ext     extension at correlation 0.99  series "<method> w/ Extension"
# For each replicate and each N in {3, 5, 10, 15, 20, 25}:
#   causal set   the N variants with the largest |beta|, with N capped at the
#                number of nonzero effects (23 in the polygenic background
#                scenario, so its N = 25 point uses all 23 causal variants)
#   Power        fraction of the causal set captured by at least one CS
#   FDR          fraction of credible sets that contain no variant of the
#                causal set
# A replicate with no credible set contributes Power = 0 and FDR = 0. The
# summary is the mean and standard error over replicates.
#
# Input:  <benchmark>/final_scripts/500_rep_results/
#           {complex,complex_S1,complex_S2}/sim_*_ldext.rds
# Output: data/s9_finemapping_data.rds

source("R/paths.R")
source("R/topn_panels.R")

ldext_root  <- file.path(benchmark_root(), "final_scripts", "500_rep_results")
output_file <- file.path(supp_dir(9), "data", "s9_finemapping_data.rds")

N_values     <- c(3, 5, 10, 15, 20, 25)
versions     <- c("xcorr_noext", "xcorr_ext")
method_names <- c("susie" = "SuSiE", "susie_inf" = "SuSiE-inf", "susie_ash" = "SuSiE-ash")

scenario_files <- c(
  "Complex" =
    "complex/sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInf15_n1000_ldext.rds",
  "Complex S1" =
    "complex_S1/sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInfall_n1000_ldext.rds",
  "Complex S2" =
    "complex_S2/sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.15_pInf0.35_nOligo10_nInfall_n1000_ldext.rds"
)

all_rows <- list()

for (label in names(scenario_files)) {
  ldext_file <- file.path(ldext_root, scenario_files[[label]])
  if (!file.exists(ldext_file)) stop("Missing ldext file: ", ldext_file)

  for (rep in readRDS(ldext_file)$replicates) {
    if (is.null(rep) || !identical(rep$status, "ok")) next

    beta       <- rep$beta
    ranked_idx <- order(abs(beta), decreasing = TRUE)
    n_nonzero  <- sum(beta != 0)

    for (N in N_values) {
      causal_set <- ranked_idx[seq_len(min(N, n_nonzero))]
      for (mkey in names(method_names)) {
        if (is.null(rep$credible_sets[[mkey]])) next   # method had no fit
        for (v in versions) {
          m <- topn_cs_metrics(rep$credible_sets[[mkey]][[v]]$cs, causal_set)
          all_rows[[length(all_rows) + 1]] <- data.frame(
            scenario = label,
            Method   = method_names[[mkey]],
            version  = v,
            N        = N,
            Power    = m$Power,
            FDR      = m$FDR
          )
        }
      }
    }
  }
}

rep_data <- do.call(rbind, all_rows)

# Mean and standard error over replicates, per scenario, method, version and N.
agg <- do.call(rbind, lapply(
  split(rep_data, list(rep_data$scenario, rep_data$Method, rep_data$version, rep_data$N), drop = TRUE),
  function(df) {
    data.frame(
      scenario   = df$scenario[1],
      Method     = df$Method[1],
      version    = df$version[1],
      N          = df$N[1],
      Power_mean = mean(df$Power, na.rm = TRUE),
      Power_se   = sd(df$Power, na.rm = TRUE) / sqrt(sum(!is.na(df$Power))),
      FDR_mean   = mean(df$FDR, na.rm = TRUE),
      FDR_se     = sd(df$FDR, na.rm = TRUE) / sqrt(sum(!is.na(df$FDR)))
    )
  }))
rownames(agg) <- NULL

# Series labels and order as in S8.
series_levels <- c("SuSiE", "SuSiE w/ Extension",
                   "SuSiE-ash", "SuSiE-ash w/ Extension",
                   "SuSiE-inf", "SuSiE-inf w/ Extension")
agg$series <- ifelse(agg$version == "xcorr_ext",
                     paste0(agg$Method, " w/ Extension"), agg$Method)
agg$series <- factor(agg$series, levels = series_levels)

# Scenario levels run from the top row of the figure to the bottom row.
agg$scenario <- factor(agg$scenario, levels = c("Complex S2", "Complex S1", "Complex"))

saveRDS(agg, output_file)
cat("Saved:", output_file, "\n")
