#!/usr/bin/env Rscript

# Supplementary Figure S7, panels A-D: per-replicate power, FDR, credible set
# size and purity for five ways of building 95% credible sets from the SuSiE
# fits of the sparse simulations (K = 1-5 causal variants, 500 replicates each).
#   SuSiE(0.5)         susie_get_cs with purity >= 0.5 (the susieR default)
#   SuSiE(0.8)         the same sets, keeping those with purity >= 0.8
#   SuSiE+attainable   attainable-coverage sets (ethres = 500), built from
#                      alpha alone with no purity filter
#   SuSiE+BLiP(0.5)    BLiP groups at q = 0.05 with no size cap, candidate
#                      groups filtered at purity 0.5
#   SuSiE+BLiP(0.8)    the same with the filter at purity 0.8
# The BLiP sets come from Simulation_Studies/scripts/fine_mapping/blip_rerun
# (pyblip 0.4.2; arms D_purity_uncapped and D_purity08_uncapped, written in
# this script's input format by combine_blip_to_s7.R).
#
# Purity is the smallest absolute correlation within a set, from R = cor(X).
# X is rebuilt for each replicate from the LD-block genotype file named in
# rep$ld_block_name (first 1000 samples, mean imputation, MAF >= 0.01,
# non-constant columns, first 5000 variants, scaled). A replicate with no
# credible set contributes power 0 and FDR 0.
#
# Replicates run in parallel in forked workers. Each K is saved to
# data/s7_metrics_K<k>.rds when it finishes and is reused on a rerun.
# Environment variables:
#   S7_N_CORES         number of workers (default: all cores but two)
#   S7_FORCE           if set, recompute K values that already have a saved file
#   S7_TEST_K          smoke test: comma-separated K values to run
#   S7_TEST_MAX_REPS   smoke test: number of replicates per K
# A smoke test neither reads nor writes the per-K files and saves
# data/s7_metrics_test.rds instead of the real output.
#
# Input:  <benchmark>/final_scripts/500_rep_results/
#           sparse/sim_nrep500_sparse_h2persnp0.03_K{1-5}_n1000.rds
#           blip_purity{0.5,0.8}/sim_nrep500_sparse_h2persnp0.03_K{1-5}_n1000_blip.rds
#           ld_blocks/<ld_block_name>
# Output: data/s7_metrics.rds  (list: $rep_metrics, $cs_purity, $cs_size, $meta)

# One BLAS thread per worker, so that the workers do not oversubscribe the
# cores. Must be set before any library initialises BLAS.
Sys.setenv(OPENBLAS_NUM_THREADS    = "1",
           VECLIB_MAXIMUM_THREADS  = "1",
           OMP_NUM_THREADS         = "1",
           MKL_NUM_THREADS         = "1",
           BLIS_NUM_THREADS        = "1")

suppressMessages({
  library(susieR)
  library(parallel)
})

source("R/paths.R")

s7_dir          <- supp_dir(7)
rep_results_dir <- file.path(benchmark_root(), "final_scripts", "500_rep_results")

sim_dir        <- file.path(rep_results_dir, "sparse")    # raw SuSiE-fit sim files
blip_dirs      <- c("SuSiE+BLiP(0.5)" = file.path(rep_results_dir, "blip_purity0.5"),
                    "SuSiE+BLiP(0.8)" = file.path(rep_results_dir, "blip_purity0.8"))
LD_BLOCKS_DIR  <- file.path(rep_results_dir, "ld_blocks") # genotype matrices (per LD block)

output_dir     <- file.path(s7_dir, "data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Configuration ----

K_values   <- 1:5
h2_per_snp <- 0.03
n_target   <- 1000

.s7_cores_env <- suppressWarnings(as.integer(Sys.getenv("S7_N_CORES", unset = "")))
N_CORES <- if (!is.na(.s7_cores_env)) .s7_cores_env else max(1L, parallel::detectCores() - 2L)

.s7_test_k <- Sys.getenv("S7_TEST_K", unset = "")
.s7_test_max_reps <- suppressWarnings(as.integer(Sys.getenv("S7_TEST_MAX_REPS", unset = "")))
if (nzchar(.s7_test_k)) {
  K_values <- as.integer(strsplit(.s7_test_k, ",")[[1]])
  cat(sprintf("[TEST MODE] K_values overridden to: %s\n",
              paste(K_values, collapse = ",")))
}
if (!is.na(.s7_test_max_reps)) {
  cat(sprintf("[TEST MODE] max reps per K capped at %d\n", .s7_test_max_reps))
}

cat(sprintf("Using N_CORES = %d (system has %d cores)\n",
            N_CORES, parallel::detectCores()))

.s7_force      <- nzchar(Sys.getenv("S7_FORCE", unset = ""))
.test_mode     <- nzchar(.s7_test_k) || !is.na(.s7_test_max_reps)

method_levels <- c("SuSiE(0.5)", "SuSiE(0.8)", "SuSiE+attainable",
                   "SuSiE+BLiP(0.5)", "SuSiE+BLiP(0.8)")

# ---- Genotype matrix of one LD block ----

reconstruct_genotype_matrix <- function(ld_block_path, n = NULL) {
  ld_block <- readRDS(ld_block_path)
  X <- ld_block$genotypes

  if (!is.null(n) && n < nrow(X)) {
    X <- X[1:n, , drop = FALSE]
  }

  col_means <- colMeans(X, na.rm = TRUE)
  na_idx <- which(is.na(X), arr.ind = TRUE)
  if (nrow(na_idx) > 0) {
    X[na_idx] <- col_means[na_idx[, 2]]
  }

  maf <- colMeans(X) / 2
  maf <- pmin(maf, 1 - maf)
  X <- X[, maf >= 0.01, drop = FALSE]

  col_vars <- apply(X, 2, var)
  X <- X[, col_vars > 0, drop = FALSE]

  if (ncol(X) > 5000) {
    X <- X[, 1:5000, drop = FALSE]
  }

  X <- scale(X)
  return(X)
}

# ---- Attainable-coverage credible sets ----

entropy <- function(y) {
  if (length(which(y > 0)) == 0) return(Inf)
  y <- y / sum(y)
  H <- -sum(y[y > 0] * log(y[y > 0]))
  return(H)
}

# Each variant is assigned to the single effect with the largest alpha for it;
# an effect's attainable coverage is the alpha mass of its assigned variants.
# Credible sets are built for the effects whose attainable coverage exceeds
# `coverage` and whose entropy is below log(ethres).
susie_get_cs_attainable <- function(res, coverage = 0.95, ethres = 500, ...) {
  alpha_attainable <- do.call(cbind, apply(res$alpha, 2,
                                           function(x) ifelse(x == max(x), x, 0),
                                           simplify = FALSE))
  coverage_attained <- rowSums(alpha_attainable)
  entropy_vals <- apply(alpha_attainable, 1, entropy)
  keep_effects <- which(coverage_attained > coverage & entropy_vals < log(ethres))

  if (length(keep_effects) == 0) {
    return(list(cs = NULL, coverage = NULL, requested_coverage = coverage))
  }

  res_filtered <- res
  res_filtered$alpha <- res$alpha[keep_effects, , drop = FALSE]
  res_filtered$V <- res$V[keep_effects]

  dot_args <- list(...)
  if (!is.null(dot_args$X)) dot_args$X <- NULL

  return(do.call(susie_get_cs, c(list(res_filtered, coverage = coverage), dot_args)))
}

# ---- Credible set metrics ----

# Power: fraction of causal variants inside some credible set. FDR: fraction of
# credible sets with no causal variant. Size: mean number of variants per set.
# All three are 0 when no set is reported.
calculate_cs_metrics <- function(test.cs, causal_indices) {
  cs_size <- 0
  cs_fdr <- 0
  cs_power <- 0

  if (length(test.cs) > 0) {
    cs_size <- length(unlist(test.cs)) / length(test.cs)
    TP_fdr <- sum(sapply(test.cs, function(cs) any(cs %in% causal_indices)))
    FP_fdr <- length(test.cs) - TP_fdr
    cs_fdr <- if ((TP_fdr + FP_fdr) > 0) FP_fdr / (TP_fdr + FP_fdr) else 0
    TP_power <- sum(causal_indices %in% unlist(test.cs))
    FN_power <- length(causal_indices) - TP_power
    cs_power <- TP_power / (TP_power + FN_power)
  }

  return(list(power = cs_power, fdr = cs_fdr, size = cs_size))
}

# Purity of each set: the smallest absolute correlation among its variants
# (1 for a single-variant set).
compute_purity_per_cs <- function(cs_list, R) {
  if (length(cs_list) == 0) return(numeric(0))
  vapply(cs_list, function(idx) {
    if (length(idx) <= 1) return(1.0)
    sub <- abs(R[idx, idx, drop = FALSE])
    diag(sub) <- NA_real_
    min(sub, na.rm = TRUE)
  }, numeric(1))
}

# Rows for one variant in one replicate: one summary row, plus one row per
# credible set for size and for purity (NULL when no set is reported).
build_variant_rows <- function(method, cs_list, purity_vec, causal,
                                k, total_pve, rep_id) {
  m <- calculate_cs_metrics(cs_list, causal)

  rep_row <- data.frame(
    Method       = method,
    K            = k,
    Total_PVE    = total_pve,
    replicate_id = rep_id,
    Power        = m$power,
    FDR          = m$fdr,
    Size         = if (length(cs_list) > 0) m$size else NA_real_,
    n_cs         = length(cs_list),
    stringsAsFactors = FALSE
  )

  if (length(cs_list) == 0) {
    return(list(rep_row = rep_row, size_rows = NULL, purity_rows = NULL))
  }

  size_rows <- data.frame(
    Method       = method,
    K            = k,
    Total_PVE    = total_pve,
    replicate_id = rep_id,
    size         = vapply(cs_list, length, integer(1)),
    stringsAsFactors = FALSE
  )

  purity_rows <- data.frame(
    Method       = method,
    K            = k,
    Total_PVE    = total_pve,
    replicate_id = rep_id,
    purity       = purity_vec,
    stringsAsFactors = FALSE
  )

  list(rep_row = rep_row, size_rows = size_rows, purity_rows = purity_rows)
}

# ---- One replicate ----

process_rep <- function(i, rep_sim, rep_blips, k, total_pve) {
  if (is.null(rep_sim) || any(vapply(rep_blips, is.null, logical(1)))) {
    return(list(status = "null_input"))
  }

  rep_id <- if (!is.null(rep_sim$replicate_id)) rep_sim$replicate_id else i

  fit    <- rep_sim$fits$susie$fit
  causal <- rep_sim$causal_indices

  ld_path <- file.path(LD_BLOCKS_DIR, rep_sim$ld_block_name)
  if (!file.exists(ld_path)) {
    return(list(status = "missing_ld_block", rep_id = rep_id,
                ld_block_name = rep_sim$ld_block_name))
  }

  t0 <- Sys.time()

  X <- tryCatch(reconstruct_genotype_matrix(ld_path, n = n_target),
                error = function(e) NULL)
  if (is.null(X)) {
    return(list(status = "X_recon_failed", rep_id = rep_id))
  }
  R <- cor(X)

  s05    <- susie_get_cs(fit, Xcorr = R, coverage = 0.95, min_abs_corr = 0.5)
  cs_05  <- if (is.null(s05$cs)) list() else s05$cs
  pur_05 <- if (length(cs_05) > 0) s05$purity$min.abs.corr else numeric(0)

  keep08 <- which(pur_05 >= 0.8)
  cs_08  <- cs_05[keep08]
  pur_08 <- pur_05[keep08]

  # BLiP at q = 0.05 pairs with the 95% credible sets; one list per purity arm.
  cs_b05  <- rep_blips[["SuSiE+BLiP(0.5)"]]$credible_sets$susie_blip_q0.05$cs
  if (is.null(cs_b05)) cs_b05 <- list()
  pur_b05 <- compute_purity_per_cs(cs_b05, R)

  cs_b08  <- rep_blips[["SuSiE+BLiP(0.8)"]]$credible_sets$susie_blip_q0.05$cs
  if (is.null(cs_b08)) cs_b08 <- list()
  pur_b08 <- compute_purity_per_cs(cs_b08, R)

  # The attainable sets use alpha only. No X or Xcorr is passed, because either
  # would switch on the default purity filter of susie_get_cs
  # (min_abs_corr = 0.5). Their purity is computed afterwards from R, for
  # panel C only.
  sat    <- susie_get_cs_attainable(fit, coverage = 0.95, ethres = 500)
  cs_at  <- if (is.null(sat$cs)) list() else sat$cs
  pur_at <- compute_purity_per_cs(cs_at, R)

  variant_specs <- list(
    list(method = "SuSiE(0.5)",       cs = cs_05,  pur = pur_05),
    list(method = "SuSiE(0.8)",       cs = cs_08,  pur = pur_08),
    list(method = "SuSiE+attainable", cs = cs_at,  pur = pur_at),
    list(method = "SuSiE+BLiP(0.5)",  cs = cs_b05, pur = pur_b05),
    list(method = "SuSiE+BLiP(0.8)",  cs = cs_b08, pur = pur_b08)
  )

  rep_rows    <- vector("list", length(variant_specs))
  size_rows   <- vector("list", length(variant_specs))
  purity_rows <- vector("list", length(variant_specs))
  for (j in seq_along(variant_specs)) {
    v <- variant_specs[[j]]
    bv <- build_variant_rows(v$method, v$cs, v$pur, causal, k, total_pve, rep_id)
    rep_rows[[j]]    <- bv$rep_row
    size_rows[[j]]   <- bv$size_rows
    purity_rows[[j]] <- bv$purity_rows
  }

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("    [pid=%d] K=%d rep=%d done (%.1fs)\n",
              Sys.getpid(), k, rep_id, elapsed))

  list(
    status      = "ok",
    rep_id      = rep_id,
    elapsed     = elapsed,
    rep_rows    = do.call(rbind, rep_rows),
    size_rows   = do.call(rbind, Filter(Negate(is.null), size_rows)),
    purity_rows = do.call(rbind, Filter(Negate(is.null), purity_rows))
  )
}

# ---- All replicates of each K ----

job_start <- Sys.time()
total_reps_done <- 0L
all_rep_metrics <- list()
all_cs_purity   <- list()
all_cs_size     <- list()

for (k in K_values) {
  total_pve <- h2_per_snp * k

  per_k_path <- file.path(output_dir, sprintf("s7_metrics_K%d.rds", k))

  # Reuse the saved file for this K, unless in a smoke test or S7_FORCE is set.
  if (!.test_mode && !.s7_force && file.exists(per_k_path)) {
    cat(sprintf("\n=== K = %d (PVE = %.2f) — RESUMING from %s ===\n",
                k, total_pve, per_k_path))
    cached <- readRDS(per_k_path)
    all_rep_metrics[[length(all_rep_metrics) + 1]] <- cached$rep_metrics
    if (!is.null(cached$cs_purity) && nrow(cached$cs_purity) > 0) {
      all_cs_purity[[length(all_cs_purity) + 1]] <- cached$cs_purity
    }
    if (!is.null(cached$cs_size) && nrow(cached$cs_size) > 0) {
      all_cs_size[[length(all_cs_size) + 1]] <- cached$cs_size
    }
    total_reps_done <- total_reps_done + length(unique(cached$rep_metrics$replicate_id))
    next
  }

  sim_path   <- file.path(sim_dir, sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000.rds", k))
  blip_paths <- file.path(blip_dirs, sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000_blip.rds", k))
  names(blip_paths) <- names(blip_dirs)

  cat(sprintf("\n=== K = %d (PVE = %.2f) ===\n", k, total_pve))
  cat("Loading sim file:   ", sim_path, "\n")
  for (nm in names(blip_paths)) cat(sprintf("Loading BLiP file:   %s (%s)\n", blip_paths[[nm]], nm))

  sim   <- readRDS(sim_path)
  blips <- lapply(blip_paths, readRDS)

  n_reps <- length(sim$replicates)
  if (!is.na(.s7_test_max_reps)) n_reps <- min(n_reps, .s7_test_max_reps)

  cat(sprintf("Dispatching %d reps across %d cores...\n", n_reps, N_CORES))
  k_start <- Sys.time()

  results <- mclapply(seq_len(n_reps), function(i) {
    rep_blips <- lapply(blips, function(b) b$replicates[[i]])
    process_rep(i, sim$replicates[[i]], rep_blips, k, total_pve)
  }, mc.cores = N_CORES, mc.preschedule = FALSE)

  k_elapsed <- as.numeric(difftime(Sys.time(), k_start, units = "secs"))
  cat(sprintf("  K=%d: %d reps in %.1f min (mean %.1fs/rep wall, %.1fs total/N_CORES)\n",
              k, n_reps, k_elapsed / 60, k_elapsed, k_elapsed))

  k_rep_rows <- list()
  k_size_rows <- list()
  k_purity_rows <- list()
  for (r in results) {
    if (inherits(r, "try-error")) {
      cat("  WARN: worker errored:", as.character(r), "\n")
      next
    }
    if (is.null(r) || is.null(r$status) || r$status != "ok") {
      if (!is.null(r$status)) {
        cat(sprintf("  skipped rep id=%s (status=%s)\n",
                    if (!is.null(r$rep_id)) r$rep_id else "?", r$status))
      }
      next
    }
    k_rep_rows[[length(k_rep_rows) + 1]] <- r$rep_rows
    if (!is.null(r$size_rows))   k_size_rows[[length(k_size_rows) + 1]]     <- r$size_rows
    if (!is.null(r$purity_rows)) k_purity_rows[[length(k_purity_rows) + 1]] <- r$purity_rows
    total_reps_done <- total_reps_done + 1L
  }

  k_rep_df    <- do.call(rbind, k_rep_rows)
  k_size_df   <- if (length(k_size_rows) > 0)   do.call(rbind, k_size_rows)   else NULL
  k_purity_df <- if (length(k_purity_rows) > 0) do.call(rbind, k_purity_rows) else NULL

  # A smoke test does not save per-K files, so it cannot replace real results.
  if (!.test_mode) {
    per_k_out <- list(
      rep_metrics = k_rep_df,
      cs_purity   = k_purity_df,
      cs_size     = k_size_df,
      meta = list(
        K              = k,
        total_pve      = total_pve,
        n_cores        = N_CORES,
        elapsed_seconds = k_elapsed,
        timestamp      = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      )
    )
    saveRDS(per_k_out, per_k_path)
    cat(sprintf("  Saved per-K file: %s\n", per_k_path))
  }

  all_rep_metrics[[length(all_rep_metrics) + 1]] <- k_rep_df
  if (!is.null(k_size_df))   all_cs_size[[length(all_cs_size) + 1]]     <- k_size_df
  if (!is.null(k_purity_df)) all_cs_purity[[length(all_cs_purity) + 1]] <- k_purity_df

  rm(sim, blips, results, k_rep_rows, k_size_rows, k_purity_rows,
     k_rep_df, k_size_df, k_purity_df); gc(verbose = FALSE)
}

# ---- Assemble and save ----

cat("\nAssembling output data frames...\n")

rep_metrics <- do.call(rbind, all_rep_metrics)
cs_purity   <- do.call(rbind, all_cs_purity)
cs_size     <- do.call(rbind, all_cs_size)

for (df_name in c("rep_metrics", "cs_purity", "cs_size")) {
  df <- get(df_name)
  df$Method <- factor(df$Method, levels = method_levels)
  rownames(df) <- NULL
  assign(df_name, df)
}

out <- list(
  rep_metrics = rep_metrics,
  cs_purity   = cs_purity,
  cs_size     = cs_size,
  meta = list(
    K_values             = K_values,
    h2_per_snp           = h2_per_snp,
    n_target             = n_target,
    method_levels        = method_levels,
    n_cores              = N_CORES,
    total_reps_processed = total_reps_done,
    elapsed_seconds      = as.numeric(difftime(Sys.time(), job_start, units = "secs"))
  )
)

# A smoke test writes beside the real file, so that it never replaces it.
output_file <- file.path(output_dir, if (.test_mode) "s7_metrics_test.rds" else "s7_metrics.rds")
saveRDS(out, output_file)

cat(sprintf("\nDone. %d reps processed in %.1f minutes (N_CORES = %d).\n",
            total_reps_done, out$meta$elapsed_seconds / 60, N_CORES))
cat(sprintf("Saved: %s\n", output_file))
cat(sprintf("rep_metrics: %d rows  |  cs_purity: %d rows  |  cs_size: %d rows\n",
            nrow(rep_metrics), nrow(cs_purity), nrow(cs_size)))
