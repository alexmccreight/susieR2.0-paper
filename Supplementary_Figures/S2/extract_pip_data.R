#!/usr/bin/env Rscript

# Supplementary Figure S2, panels C-J: variant-level ROC curves and PIP
# calibration in the sparse simulations.
#
# PIPs and causal indicators are pooled over the 500 replicates of each K.
#   ROC (C-G):         per K, TPR and FPR as the PIP threshold falls from 1 to
#                      0, thinned to about 500 points per curve; AUROC is taken
#                      from the full curve.
#   Calibration (H-J): pooled over K = 1-5, the observed causal frequency in
#                      ten PIP bins of width 0.1.
#
# Input:  <benchmark>/final_scripts/500_rep_results/sparse/
#           sim_nrep500_sparse_h2persnp0.03_K{1-5}_n1000.rds
# Output: data/pip_roc_data.rds          (list: $curves, $auroc)
#         data/pip_calibration_data.rds

source("R/paths.R")

input_dir  <- file.path(benchmark_root(), "final_scripts", "500_rep_results", "sparse")
output_dir <- file.path(supp_dir(2), "data")

method_names <- c("susie" = "SuSiE", "susie_inf" = "SuSiE-inf", "susie_ash" = "SuSiE-ash")
methods      <- names(method_names)

roc_curves <- list()
roc_auroc  <- list()
pooled     <- lapply(method_names, function(x) list(pip = list(), is_causal = list()))

for (K in 1:5) {
  reps <- readRDS(file.path(
    input_dir, sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000.rds", K)
  ))$replicates

  for (m in methods) {
    fitted <- reps[!vapply(reps, function(rep) is.null(rep$fits[[m]]), logical(1))]
    pip <- unlist(lapply(fitted, function(rep) rep$fits[[m]]$fit$pip), use.names = FALSE)
    is_causal <- unlist(lapply(fitted, function(rep) {
      indicator <- rep(0L, length(rep$beta))
      indicator[rep$causal_indices] <- 1L
      indicator
    }))

    pooled[[m]]$pip[[K]]       <- pip
    pooled[[m]]$is_causal[[K]] <- is_causal

    # ROC curve: walk down the variants in decreasing order of PIP.
    is_causal_sorted <- is_causal[order(pip, decreasing = TRUE)]
    tpr <- cumsum(is_causal_sorted == 1) / sum(is_causal == 1)
    fpr <- cumsum(is_causal_sorted == 0) / sum(is_causal == 0)

    n   <- length(tpr)
    idx <- unique(c(1, seq(1, n, length.out = 500), n))

    # Trapezoidal rule on the full curve.
    tpr_full <- c(0, tpr)
    fpr_full <- c(0, fpr)
    auroc <- sum(diff(fpr_full) * (tpr_full[-1] + tpr_full[-length(tpr_full)]) / 2)

    roc_curves[[length(roc_curves) + 1]] <- data.frame(
      Method = method_names[[m]], K = K, FPR = c(0, fpr[idx]), TPR = c(0, tpr[idx])
    )
    roc_auroc[[length(roc_auroc) + 1]] <- data.frame(
      Method = method_names[[m]], K = K, AUROC = auroc
    )
  }
}

curve_data <- do.call(rbind, roc_curves)
auroc_data <- do.call(rbind, roc_auroc)
rownames(curve_data) <- NULL
rownames(auroc_data) <- NULL

# Calibration bins: [0, 0.1), [0.1, 0.2), ..., [0.9, 1].
bin_breaks <- seq(0, 1, by = 0.1)
bin_mids   <- (bin_breaks[-length(bin_breaks)] + bin_breaks[-1]) / 2

cal_rows <- list()
for (m in methods) {
  pip       <- unlist(pooled[[m]]$pip)
  is_causal <- unlist(pooled[[m]]$is_causal)
  bin_idx   <- findInterval(pip, bin_breaks, rightmost.closed = TRUE)
  bin_idx[bin_idx == 0] <- 1  # a PIP marginally below 0 goes in the first bin

  for (b in 1:10) {
    in_bin   <- bin_idx == b
    n_in_bin <- sum(in_bin)
    obs_freq <- if (n_in_bin == 0) NA else mean(is_causal[in_bin])

    cal_rows[[length(cal_rows) + 1]] <- data.frame(
      Method      = method_names[[m]],
      bin_mid     = bin_mids[b],
      bin_lower   = bin_breaks[b],
      bin_upper   = bin_breaks[b + 1],
      n_variants  = n_in_bin,
      observed    = obs_freq,
      observed_se = sqrt(obs_freq * (1 - obs_freq) / n_in_bin)
    )
  }
}

cal_data <- do.call(rbind, cal_rows)
rownames(cal_data) <- NULL

roc_file <- file.path(output_dir, "pip_roc_data.rds")
cal_file <- file.path(output_dir, "pip_calibration_data.rds")
saveRDS(list(curves = curve_data, auroc = auroc_data), roc_file)
saveRDS(cal_data, cal_file)
cat("Saved:", roc_file, "\n")
cat("Saved:", cal_file, "\n")
