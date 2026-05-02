#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S4: Extract per-rep metrics for 5 SuSiE CS variants
# =============================================================================
# Variants compared:
#   1. SuSiE (purity >= 0.5)               - default susieR fit
#   2. SuSiE (purity >= 0.8)               - same fit, post-filtered
#   3. SuSiE+BLiP (q = 0.05)               - BLiP CS list from combined file
#   4. SuSiE+BLiP (q = 0.10)               - BLiP CS list from combined file
#   5. SuSiE+attainable (ethres = 500)     - attainable-coverage CS from fit
#
# For each replicate of each K in 1..5:
#   - Reconstruct genotype matrix X from the local LD-block file referenced
#     by rep$ld_block_name (apply MAF/var filter, cap 5000, scale).
#   - Compute correlation matrix R = cor(X).
#   - Build the 5 CS lists. Compute power, FDR, mean CS size, and per-CS
#     purity (min |corr|) for each.
#
# Inputs (all under final_scripts/500_rep_results/):
#   - sparse/sim_nrep500_sparse_h2persnp0.03_K{1..5}_n1000.rds   (raw SuSiE fits)
#   - blip/sim_nrep500_sparse_h2persnp0.03_K{1..5}_n1000_blip.rds (combined BLiP)
#   - ld_blocks/<ld_block_name>.rds                                (genotype matrices)
#
# Output:
#   - {s4_dir}/data/s4_metrics.rds  (list of 3 long-format data frames; only
#     file user needs to reproduce S4.{pdf,png} via plot_S4.R)
#
# Parallelism:
#   - Replicates within each K are processed in parallel via parallel::mclapply
#     (forked workers).
#   - BLAS threads are pinned to 1 per worker so workers don't oversubscribe
#     cores. Set N_CORES below.
# =============================================================================

# --- Pin BLAS threads BEFORE loading any library that initializes BLAS ---
Sys.setenv(OPENBLAS_NUM_THREADS    = "1",
           VECLIB_MAXIMUM_THREADS  = "1",
           OMP_NUM_THREADS         = "1",
           MKL_NUM_THREADS         = "1",
           BLIS_NUM_THREADS        = "1")

suppressMessages({
  library(susieR)
  library(parallel)
})

# =============================================================================
# Paths
# =============================================================================

s4_dir          <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/supplemental/S4"
rep_results_dir <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/500_rep_results"

# Raw inputs live under final_scripts/500_rep_results/ (not user-facing in S4/).
sim_dir        <- file.path(rep_results_dir, "sparse")    # raw SuSiE-fit sim files
blip_dir       <- file.path(rep_results_dir, "blip")      # combined BLiP outputs
LD_BLOCKS_DIR  <- file.path(rep_results_dir, "ld_blocks") # genotype matrices (per LD block)

output_dir     <- file.path(s4_dir, "data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Configuration
# =============================================================================

K_values   <- 1:5
h2_per_snp <- 0.03
n_target   <- 1000

# Cores for parallel rep processing. Default: leave 2 for system + I/O.
# Override via env var S4_N_CORES.
.s4_cores_env <- suppressWarnings(as.integer(Sys.getenv("S4_N_CORES", unset = "")))
N_CORES <- if (!is.na(.s4_cores_env)) .s4_cores_env else max(1L, parallel::detectCores() - 2L)

# --- Opt-in smoke-test overrides via environment variables ---
# S4_TEST_K        : comma-separated K values to run (e.g. "1" or "1,2")
# S4_TEST_MAX_REPS : integer cap on reps processed per K
.s4_test_k <- Sys.getenv("S4_TEST_K", unset = "")
.s4_test_max_reps <- suppressWarnings(as.integer(Sys.getenv("S4_TEST_MAX_REPS", unset = "")))
if (nzchar(.s4_test_k)) {
  K_values <- as.integer(strsplit(.s4_test_k, ",")[[1]])
  cat(sprintf("[TEST MODE] K_values overridden to: %s\n",
              paste(K_values, collapse = ",")))
}
if (!is.na(.s4_test_max_reps)) {
  cat(sprintf("[TEST MODE] max reps per K capped at %d\n", .s4_test_max_reps))
}

cat(sprintf("Using N_CORES = %d (system has %d cores)\n",
            N_CORES, parallel::detectCores()))

# In normal mode, per-K results are written incrementally to data/s4_metrics_K{k}.rds.
# On a re-run, an existing per-K file is reused (skipping its mclapply) unless
# S4_FORCE=1 is set. Test mode (S4_TEST_K / S4_TEST_MAX_REPS) bypasses both
# saving and reuse so it never corrupts real results.
.s4_force      <- nzchar(Sys.getenv("S4_FORCE", unset = ""))
.test_mode     <- nzchar(.s4_test_k) || !is.na(.s4_test_max_reps)

method_levels <- c("SuSiE(0.5)", "SuSiE(0.8)",
                   "SuSiE+BLiP(q0.05)", "SuSiE+BLiP(q0.10)",
                   "SuSiE+attainable")

# =============================================================================
# Helper: Reconstruct preprocessed genotype matrix from an LD block file
# =============================================================================

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

# =============================================================================
# Helpers (verbatim from user): entropy + susie_get_cs_attainable +
# calculate_cs_metrics
# =============================================================================

entropy <- function(y) {
  if (length(which(y > 0)) == 0) return(Inf)
  y <- y / sum(y)
  H <- -sum(y[y > 0] * log(y[y > 0]))
  return(H)
}

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

# =============================================================================
# Per-CS purity given correlation matrix R
# =============================================================================

compute_purity_per_cs <- function(cs_list, R) {
  if (length(cs_list) == 0) return(numeric(0))
  vapply(cs_list, function(idx) {
    if (length(idx) <= 1) return(1.0)
    sub <- abs(R[idx, idx, drop = FALSE])
    diag(sub) <- NA_real_
    min(sub, na.rm = TRUE)
  }, numeric(1))
}

# =============================================================================
# Pure function: build all rows for one variant of one rep
# Returns a list with three optional row-data-frames.
# =============================================================================

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

# =============================================================================
# Per-rep worker: takes one rep's sim and BLiP slices and returns row tables
# =============================================================================

process_rep <- function(i, rep_sim, rep_blip, k, total_pve) {
  if (is.null(rep_sim) || is.null(rep_blip)) {
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

  # 5 CS lists + purities
  s05    <- susie_get_cs(fit, Xcorr = R, coverage = 0.95, min_abs_corr = 0.5)
  cs_05  <- if (is.null(s05$cs)) list() else s05$cs
  pur_05 <- if (length(cs_05) > 0) s05$purity$min.abs.corr else numeric(0)

  keep08 <- which(pur_05 >= 0.8)
  cs_08  <- cs_05[keep08]
  pur_08 <- pur_05[keep08]

  cs_b05  <- rep_blip$credible_sets$susie_blip_q0.05$cs
  if (is.null(cs_b05)) cs_b05 <- list()
  pur_b05 <- compute_purity_per_cs(cs_b05, R)

  cs_b10  <- rep_blip$credible_sets$susie_blip_q0.1$cs
  if (is.null(cs_b10)) cs_b10 <- list()
  pur_b10 <- compute_purity_per_cs(cs_b10, R)

  # Method: alpha-only, truly post-hoc — no X / Xcorr passed.
  # (Passing Xcorr would silently invoke susie_get_cs's default
  # min_abs_corr = 0.5 filter, which contradicts the post-hoc framing.)
  sat    <- susie_get_cs_attainable(fit, coverage = 0.95, ethres = 500)
  cs_at  <- if (is.null(sat$cs)) list() else sat$cs
  # Purity is a visualization-only post-hoc computation using our R.
  pur_at <- compute_purity_per_cs(cs_at, R)

  variant_specs <- list(
    list(method = "SuSiE(0.5)",        cs = cs_05,  pur = pur_05),
    list(method = "SuSiE(0.8)",        cs = cs_08,  pur = pur_08),
    list(method = "SuSiE+BLiP(q0.05)", cs = cs_b05, pur = pur_b05),
    list(method = "SuSiE+BLiP(q0.10)", cs = cs_b10, pur = pur_b10),
    list(method = "SuSiE+attainable",  cs = cs_at,  pur = pur_at)
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

# =============================================================================
# Main loop (per K, parallel over replicates)
# =============================================================================

job_start <- Sys.time()
total_reps_done <- 0L
all_rep_metrics <- list()
all_cs_purity   <- list()
all_cs_size     <- list()

for (k in K_values) {
  total_pve <- h2_per_snp * k

  per_k_path <- file.path(output_dir, sprintf("s4_metrics_K%d.rds", k))

  # Resume: reuse existing per-K file in normal mode unless S4_FORCE=1
  if (!.test_mode && !.s4_force && file.exists(per_k_path)) {
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

  sim_path  <- file.path(sim_dir, sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000.rds", k))
  blip_path <- file.path(blip_dir, sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000_blip.rds", k))

  cat(sprintf("\n=== K = %d (PVE = %.2f) ===\n", k, total_pve))
  cat("Loading sim file:  ", sim_path,  "\n")
  cat("Loading BLiP file: ", blip_path, "\n")

  sim  <- readRDS(sim_path)
  blip <- readRDS(blip_path)

  n_reps <- length(sim$replicates)
  if (!is.na(.s4_test_max_reps)) n_reps <- min(n_reps, .s4_test_max_reps)

  cat(sprintf("Dispatching %d reps across %d cores...\n", n_reps, N_CORES))
  k_start <- Sys.time()

  results <- mclapply(seq_len(n_reps), function(i) {
    process_rep(i, sim$replicates[[i]], blip$replicates[[i]], k, total_pve)
  }, mc.cores = N_CORES, mc.preschedule = FALSE)

  k_elapsed <- as.numeric(difftime(Sys.time(), k_start, units = "secs"))
  cat(sprintf("  K=%d: %d reps in %.1f min (mean %.1fs/rep wall, %.1fs total/N_CORES)\n",
              k, n_reps, k_elapsed / 60, k_elapsed, k_elapsed))

  # Accumulate this K's rows
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

  # Per-K incremental save (skipped in test mode to avoid corrupting real data)
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

  rm(sim, blip, results, k_rep_rows, k_size_rows, k_purity_rows,
     k_rep_df, k_size_df, k_purity_df); gc(verbose = FALSE)
}

# =============================================================================
# Assemble output data frames
# =============================================================================

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

output_file <- file.path(output_dir, "s4_metrics.rds")
saveRDS(out, output_file)

cat(sprintf("\nDone. %d reps processed in %.1f minutes (N_CORES = %d).\n",
            total_reps_done, out$meta$elapsed_seconds / 60, N_CORES))
cat(sprintf("Saved: %s\n", output_file))
cat(sprintf("rep_metrics: %d rows  |  cs_purity: %d rows  |  cs_size: %d rows\n",
            nrow(rep_metrics), nrow(cs_purity), nrow(cs_size)))
