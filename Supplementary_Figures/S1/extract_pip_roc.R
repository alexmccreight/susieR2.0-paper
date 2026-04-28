#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S1: Extract PIP ROC curve data (per K)
# =============================================================================
# For each K in {1, 2, 3, 4, 5}, pools all (PIP, is_causal) pairs across 500 reps,
# sweeps PIP threshold from 1 -> 0, and computes TPR/FPR at each threshold.
# Also computes AUROC per method per K.
#
# Input:  sim_nrep500_sparse_h2persnp0.03_K{1,2,3,4,5}_n1000.rds
# Output: data/pip_roc_data.rds  (list with $curves and $auroc)
# =============================================================================

s1_dir     <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/supplemental/S1"
data_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/500_rep_results/sparse"
output_dir <- file.path(s1_dir, "data")

K_values <- 1:5
methods <- c("susie", "susie_inf", "susie_ash")
method_names <- c("susie" = "SuSiE", "susie_inf" = "SuSiE-inf", "susie_ash" = "SuSiE-ash")

all_curves <- list()
all_auroc  <- list()

for (K in K_values) {
  fname <- sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000.rds", K)
  fpath <- file.path(data_dir, fname)
  if (!file.exists(fpath)) stop("Not found: ", fpath)

  cat(sprintf("Loading K=%d ...\n", K))
  sim <- readRDS(fpath)
  reps <- sim$replicates

  # Pool PIPs per method for this K
  pooled <- list()
  for (m in methods) pooled[[m]] <- list(pip = c(), is_causal = c())

  for (r in seq_along(reps)) {
    rep <- reps[[r]]
    causal_idx <- rep$causal_indices
    p <- length(rep$beta)
    is_causal_vec <- rep(0L, p)
    is_causal_vec[causal_idx] <- 1L

    for (m in methods) {
      fit_obj <- rep$fits[[m]]
      if (is.null(fit_obj)) next
      pooled[[m]]$pip <- c(pooled[[m]]$pip, fit_obj$fit$pip)
      pooled[[m]]$is_causal <- c(pooled[[m]]$is_causal, is_causal_vec)
    }
  }

  cat(sprintf("  %d reps pooled, %d variants per method\n",
              length(reps), length(pooled[[methods[1]]]$pip)))

  # Compute ROC curve and AUROC for each method
  for (m in methods) {
    pip <- pooled[[m]]$pip
    is_causal <- pooled[[m]]$is_causal
    n_pos <- sum(is_causal == 1)
    n_neg <- sum(is_causal == 0)

    # Sort descending by PIP
    ord <- order(pip, decreasing = TRUE)
    is_causal_sorted <- is_causal[ord]

    # Full ROC curve
    tp <- cumsum(is_causal_sorted == 1)
    fp <- cumsum(is_causal_sorted == 0)
    tpr <- tp / n_pos
    fpr <- fp / n_neg

    # Subsample to ~500 evenly spaced points for plotting
    n <- length(tpr)
    idx <- unique(c(1, seq(1, n, length.out = 500), n))
    tpr_sub <- c(0, tpr[idx])
    fpr_sub <- c(0, fpr[idx])

    # AUROC via trapezoidal rule on full curve
    tpr_full <- c(0, tpr)
    fpr_full <- c(0, fpr)
    auroc <- sum(diff(fpr_full) * (tpr_full[-1] + tpr_full[-length(tpr_full)]) / 2)

    all_curves[[length(all_curves) + 1]] <- data.frame(
      Method = method_names[m], K = K, FPR = fpr_sub, TPR = tpr_sub
    )

    all_auroc[[length(all_auroc) + 1]] <- data.frame(
      Method = method_names[m], K = K, AUROC = auroc
    )
  }
}

curve_data <- do.call(rbind, all_curves)
auroc_data <- do.call(rbind, all_auroc)
rownames(curve_data) <- NULL
rownames(auroc_data) <- NULL

cat("\nAUROC values:\n")
print(auroc_data)

out_path <- file.path(output_dir, "pip_roc_data.rds")
saveRDS(list(curves = curve_data, auroc = auroc_data), out_path)
cat(sprintf("\nSaved: %s\n", out_path))
cat("Done!\n")
