#!/usr/bin/env Rscript

# =============================================================================
# Figure 1 (panels C-F): complex S1/S2 fine-mapping summaries WITH LD extension
# =============================================================================
# Recomputes Top-N Power/FDR/Size for the two complex scenarios using the
# LD-extension (0.99) credible sets stored in the *_ldext.rds files
# (version = "xcorr_ext"), replacing the previous un-extended fit$sets$cs
# construction.
#
# This mirrors the original benchmark extract_complex_metrics.R Top-N logic
# EXACTLY -- same N vector, same Power/FDR/Size definitions, same mean/SE
# aggregation -- but reads
#   rep$credible_sets[[method]][["xcorr_ext"]]$cs
# instead of
#   rep$fits[[method]]$fit$sets$cs.
#
# For each replicate and each N in {3, 5, 10, 15, 20, 25}:
#   - "Top N" causal set = the N variants with largest |beta|
#   - Power = fraction of Top-N variants captured by at least one CS
#   - FDR   = fraction of credible sets that contain no Top-N variant
#
# Inputs  (read-only, benchmark):  500_rep_results/complex_S{1,2}/sim_*_ldext.rds
# Outputs (production):            figure_1/data/complex_S{1,2}_finemapping_summary.rds
# =============================================================================

# ---- Paths ----
source("R/paths.R")

# Raw 500-replicate simulation output lives outside the repo (~22 GB).
# Override its location with SUSIER2_BENCHMARK_ROOT if needed.
ldext_root <- file.path(benchmark_root(), "final_scripts", "500_rep_results")
fig1_data  <- file.path(fig_dir(1), "data")

VERSION      <- "xcorr_ext"
N_values     <- c(3, 5, 10, 15, 20, 25)
method_names <- c("susie" = "SuSiE", "susie_inf" = "SuSiE-inf", "susie_ash" = "SuSiE-ash")

scenarios <- list(
  list(label = "Complex S1",
       out   = "complex_S1_finemapping_summary.rds",
       file  = file.path(ldext_root, "complex_S1",
         "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInfall_n1000_ldext.rds")),
  list(label = "Complex S2",
       out   = "complex_S2_finemapping_summary.rds",
       file  = file.path(ldext_root, "complex_S2",
         "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.15_pInf0.35_nOligo10_nInfall_n1000_ldext.rds"))
)

# -----------------------------------------------------------------------------
# Helper: Top-N power / FDR / mean size for a list of credible sets
# (identical formulas to benchmark extract_complex_metrics.R::compute_metrics,
#  refactored to take the credible-set list directly)
# -----------------------------------------------------------------------------
compute_metrics <- function(cs_list, causal_set) {
  if (is.null(cs_list) || length(cs_list) == 0) {
    return(data.frame(Power = 0, FDR = NA, Size = 0))
  }

  n_causal <- length(causal_set)
  n_cs     <- length(cs_list)

  causal_covered <- sapply(causal_set, function(v) {
    any(sapply(cs_list, function(cs) v %in% cs))
  })
  power <- sum(causal_covered) / n_causal

  cs_has_causal <- sapply(cs_list, function(cs) any(cs %in% causal_set))
  fdr <- sum(!cs_has_causal) / n_cs

  avg_size <- mean(sapply(cs_list, length))

  data.frame(Power = power, FDR = fdr, Size = avg_size)
}

# -----------------------------------------------------------------------------
# Helper: aggregate per-rep rows to mean +/- SE by Method x N
# (identical to benchmark aggregation)
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# Process each scenario
# -----------------------------------------------------------------------------
for (sc in scenarios) {
  cat(sprintf("\n=== %s (LD extension = 0.99) ===\n", sc$label))
  if (!file.exists(sc$file)) stop("Missing ldext file: ", sc$file)

  sim  <- readRDS(sc$file)
  reps <- sim$replicates
  K    <- sim$parameters$K
  cat(sprintf("  %d replicates, K=%d\n", length(reps), K))

  all_rows <- list()
  for (r in seq_along(reps)) {
    rep <- reps[[r]]
    if (is.null(rep) || !identical(rep$status, "ok")) next

    beta       <- rep$beta
    ranked_idx <- order(abs(beta), decreasing = TRUE)

    for (N in N_values) {
      causal_set <- ranked_idx[1:N]

      for (mkey in names(method_names)) {
        if (is.null(rep$credible_sets[[mkey]])) next   # method had no fit
        cs_list <- rep$credible_sets[[mkey]][[VERSION]]$cs

        m <- compute_metrics(cs_list, causal_set)
        all_rows[[length(all_rows) + 1]] <- data.frame(
          Method    = method_names[[mkey]],
          K         = K,
          N         = N,
          Power     = m$Power,
          FDR       = m$FDR,
          Size      = m$Size,
          replicate = r
        )
      }
    }

    if (r %% 100 == 0) cat(sprintf("  processed %d / %d replicates\n", r, length(reps)))
  }

  rep_data <- do.call(rbind, all_rows)
  summ     <- aggregate_summary(rep_data)

  out_path <- file.path(fig1_data, sc$out)
  saveRDS(summ, out_path)
  cat(sprintf("  Saved: %s (%d rows)\n", out_path, nrow(summ)))
  print(summ[summ$N %in% c(3, 25),
             c("Method", "N", "Power_mean", "FDR_mean", "Size_mean")])
}

cat("\nDone!\n")
