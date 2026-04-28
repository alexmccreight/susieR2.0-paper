#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S1: Extract PIP calibration data (pooled across K=1-5)
# =============================================================================
# Pools all (PIP, is_causal) pairs across K=1-5 (500 reps each = 2500 reps),
# bins PIPs into 10 equal-width bins, and computes observed causal frequency
# per bin per method.
#
# Input:  sim_nrep500_sparse_h2persnp0.03_K{1,2,3,4,5}_n1000.rds
# Output: data/pip_calibration_data.rds
# =============================================================================

s1_dir     <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/supplemental/S1"
data_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/500_rep_results/sparse"
output_dir <- file.path(s1_dir, "data")

K_values <- 1:5
methods <- c("susie", "susie_inf", "susie_ash")
method_names <- c("susie" = "SuSiE", "susie_inf" = "SuSiE-inf", "susie_ash" = "SuSiE-ash")

# Pool PIPs across all K values
pooled <- list()
for (m in methods) pooled[[m]] <- list(pip = c(), is_causal = c())

for (K in K_values) {
  fname <- sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000.rds", K)
  fpath <- file.path(data_dir, fname)
  if (!file.exists(fpath)) stop("Not found: ", fpath)

  cat(sprintf("Loading K=%d ...\n", K))
  sim <- readRDS(fpath)
  reps <- sim$replicates

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
}

cat(sprintf("\nPooled %d variants per method\n", length(pooled[[methods[1]]]$pip)))

# Bin into 10 bins: [0, 0.1), [0.1, 0.2), ..., [0.9, 1.0]
bin_breaks <- seq(0, 1, by = 0.1)
bin_mids   <- (bin_breaks[-length(bin_breaks)] + bin_breaks[-1]) / 2

all_rows <- list()

for (m in methods) {
  pip <- pooled[[m]]$pip
  is_causal <- pooled[[m]]$is_causal

  # Assign bins (use right=FALSE so [0,0.1), [0.1,0.2), ...; last bin is [0.9,1.0])
  bin_idx <- findInterval(pip, bin_breaks, rightmost.closed = TRUE)
  bin_idx[bin_idx == 0] <- 1  # edge case

  for (b in 1:10) {
    in_bin <- bin_idx == b
    n_in_bin <- sum(in_bin)
    if (n_in_bin == 0) {
      obs_freq <- NA
      obs_se   <- NA
    } else {
      obs_freq <- mean(is_causal[in_bin])
      obs_se   <- sqrt(obs_freq * (1 - obs_freq) / n_in_bin)
    }

    all_rows[[length(all_rows) + 1]] <- data.frame(
      Method       = method_names[m],
      bin_mid      = bin_mids[b],
      bin_lower    = bin_breaks[b],
      bin_upper    = bin_breaks[b + 1],
      n_variants   = n_in_bin,
      observed     = obs_freq,
      observed_se  = obs_se
    )
  }
}

cal_data <- do.call(rbind, all_rows)
rownames(cal_data) <- NULL

cat("\nCalibration summary:\n")
print(cal_data)

out_path <- file.path(output_dir, "pip_calibration_data.rds")
saveRDS(cal_data, out_path)
cat(sprintf("\nSaved: %s\n", out_path))
cat("Done!\n")
