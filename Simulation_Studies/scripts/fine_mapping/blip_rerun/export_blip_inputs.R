#!/usr/bin/env Rscript

# Supplementary Figure S7, SuSiE+BLiP re-run, step 1: export what BLiP needs
# from one raw sparse simulation file, for one chunk of replicates.
#
# Per replicate:
#   1. Rebuild the genotype matrix X from the replicate's LD block as
#      extract_S7_data.R does (first n rows, mean-impute, MAF >= 0.01, drop
#      zero variance, first 5000 columns, scale) and compute R = cor(X).
#   2. Compute the SuSiE credible sets at 95% and 90% coverage with
#      susie_get_cs(fit, Xcorr = R, coverage, min_abs_corr = 0), so that the
#      purity-free comparators in blip_metrics.R see every set, and record
#      each set's purity. The sets with purity >= 0.5 are checked against the
#      Figure S7 call (min_abs_corr = 0.5).
#   3. Write alpha.npy (L x p SuSiE alpha, float64), corr_sub.npy (float32)
#      and meta.json (replicate id, K, n, p, L, causal indices, PIPs, prior
#      variances, LD block name, the sets with their purity, export settings).
#      Indices are 1-based; meta.json records this as index_base.
#
# corr_sub.npy is the block of R over meta$corr_index: the SNPs with
# PIP >= pip_floor (0.01 is the prefilter run_blip_arms.py applies), the
# members of every exported set with purity >= corr_set_purity (0.5 is the
# lowest purity filter used) and the causal SNPs. Members of impure sets are
# left out on purpose: near-null components yield sets spanning most of the
# region, which would make the block the full p x p matrix (~100 MB per
# replicate), and no purity-matched arm can report such a set.
# run_blip_arms.py stops when its PIP threshold retains a SNP outside the
# block.
#
# A replicate whose meta.json exists is skipped unless overwrite=TRUE, so a
# chunk can be re-run after an interruption.
#
# Input:  sim_file       sim_nrep500_sparse_h2persnp0.03_K<k>_n1000.rds
#         ld_blocks_dir  the LD block rds files the replicates name
#         A relative path that does not exist from the working directory is
#         taken relative to <benchmark>/final_scripts/500_rep_results.
# Output: <output_dir>/K<k>/sim_parameters.json
#         <output_dir>/K<k>/rep_<id>/{alpha.npy, corr_sub.npy, meta.json}
#
# Usage: Rscript export_blip_inputs.R sim_file=<rds> ld_blocks_dir=<dir> \
#          output_dir=<root> [start_replicate=1] [end_replicate=<last>] \
#          [n_cores=1] [n=1000] [pip_floor=0.01] [corr_set_purity=0.5] \
#          [overwrite=FALSE]

# One BLAS thread per forked worker; set before any package initialises BLAS.
Sys.setenv(OPENBLAS_NUM_THREADS   = "1",
           VECLIB_MAXIMUM_THREADS = "1",
           OMP_NUM_THREADS        = "1",
           MKL_NUM_THREADS        = "1",
           BLIS_NUM_THREADS       = "1")

suppressMessages({
  library(susieR)
  library(parallel)
  library(jsonlite)
})

.script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else getwd()
})
source(file.path(.script_dir, "blip_common.R"))
locate_paper_paths(.script_dir)

COVERAGE_LEVELS <- c(0.95, 0.90)

# ---- Per-replicate worker ----

export_replicate <- function(rep_sim, k, ld_blocks_dir, k_dir, n_target,
                             pip_floor, corr_set_purity, overwrite = FALSE,
                             position = NA) {
  if (is.null(rep_sim)) return(list(status = "null_input", rep_id = position))

  rep_id <- if (!is.null(rep_sim$replicate_id)) rep_sim$replicate_id else position
  rep_dir <- file.path(k_dir, sprintf("rep_%04d", rep_id))
  meta_path <- file.path(rep_dir, "meta.json")
  if (file.exists(meta_path) && !overwrite) {
    return(list(status = "exists", rep_id = rep_id))
  }

  fit <- rep_sim$fits$susie$fit
  if (is.null(fit) || is.null(fit$alpha)) return(list(status = "no_fit", rep_id = rep_id))

  ld_path <- file.path(ld_blocks_dir, rep_sim$ld_block_name)
  if (!file.exists(ld_path)) {
    return(list(status = "missing_ld_block", rep_id = rep_id,
                ld_block_name = rep_sim$ld_block_name))
  }

  t0 <- Sys.time()
  X <- reconstruct_genotype_matrix(ld_path, n = n_target)
  R <- cor(X)
  p <- ncol(X)
  if (ncol(fit$alpha) != p) {
    stop(sprintf("alpha has %d columns but the reconstructed X has %d",
                 ncol(fit$alpha), p))
  }
  causal <- as.integer(rep_sim$causal_indices)
  pip    <- unname(susie_get_pip(fit))

  cs_by_cov <- lapply(COVERAGE_LEVELS, function(cov) {
    # With Xcorr the development susieR computes purity over the whole CS;
    # CRAN susieR subsamples CSs above n_purity SNPs, so seed both calls alike.
    set.seed(rep_id)
    s_all <- susie_get_cs(fit, Xcorr = R, coverage = cov, min_abs_corr = 0)
    set.seed(rep_id)
    s_05  <- susie_get_cs(fit, Xcorr = R, coverage = cov, min_abs_corr = 0.5)
    if (is.null(s_all$cs)) return(list())

    entries <- lapply(seq_along(s_all$cs), function(i) list(
      component       = as.integer(s_all$cs_index[i]),
      cs              = I(as.integer(s_all$cs[[i]])),
      min_abs_corr    = unname(s_all$purity$min.abs.corr[i]),
      mean_abs_corr   = unname(s_all$purity$mean.abs.corr[i]),
      median_abs_corr = unname(s_all$purity$median.abs.corr[i]),
      coverage        = unname(s_all$coverage[i])
    ))

    # The min_abs_corr >= 0.5 subset must reproduce the Figure S7 call.
    pure     <- vapply(entries, function(e) e$min_abs_corr >= 0.5, logical(1))
    from_all <- vapply(entries[pure], function(e) paste(sort(as.integer(e$cs)), collapse = ","), "")
    direct   <- if (is.null(s_05$cs)) character(0) else
      vapply(s_05$cs, function(cs) paste(sort(cs), collapse = ","), "")
    if (!setequal(from_all, direct)) {
      warning(sprintf("rep %d coverage %.2f: purity-0.5 subset differs from susie_get_cs(min_abs_corr = 0.5)",
                      rep_id, cov))
    }
    entries
  })
  names(cs_by_cov) <- as.character(COVERAGE_LEVELS)

  cs_members <- unlist(lapply(cs_by_cov, function(l) {
    lapply(l, function(e) if (e$min_abs_corr >= corr_set_purity) as.integer(e$cs) else integer(0))
  }))
  corr_index <- sort(unique(c(which(pip >= pip_floor), cs_members, causal)))

  dir.create(rep_dir, showWarnings = FALSE, recursive = TRUE)
  write_npy(unname(fit$alpha), file.path(rep_dir, "alpha.npy"), "float64")
  write_npy(R[corr_index, corr_index, drop = FALSE],
            file.path(rep_dir, "corr_sub.npy"), "float32")

  meta <- list(
    replicate_id  = rep_id,
    K             = k,
    n             = nrow(X),
    p             = p,
    L             = nrow(fit$alpha),
    index_base    = 1L,
    ld_block_name = rep_sim$ld_block_name,
    causal        = I(causal),
    pip           = I(pip),
    V             = I(if (is.null(fit$V)) numeric(0) else as.numeric(fit$V)),
    corr_index    = I(corr_index),
    susie_cs      = cs_by_cov,
    export = list(
      susieR_version  = as.character(packageVersion("susieR")),
      susieR_path     = find.package("susieR"),
      coverage_levels = I(COVERAGE_LEVELS),
      min_abs_corr    = 0.5,
      n_purity        = as.character(formals(susie_get_cs)$n_purity),
      purity_seed     = rep_id,
      pip_floor       = pip_floor,
      corr_set_purity = corr_set_purity,
      n_target        = n_target,
      timestamp       = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
  )
  # meta.json is written last: its presence marks the replicate as complete.
  write_json(meta, meta_path, auto_unbox = TRUE, digits = NA)

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("    [pid=%d] K=%d rep=%d exported (%d SNPs in corr_sub, %.1fs)\n",
              Sys.getpid(), k, rep_id, length(corr_index), elapsed))
  list(status = "ok", rep_id = rep_id, elapsed = elapsed, n_corr = length(corr_index))
}

# ---- Main ----

main <- function() {
  a <- parse_key_value_args(list(
    sim_file        = NULL,
    ld_blocks_dir   = NULL,
    output_dir      = NULL,
    start_replicate = 1L,
    end_replicate   = NA_integer_,
    n_cores         = 1L,
    n               = 1000L,
    pip_floor       = 0.01,
    corr_set_purity = 0.5,
    overwrite       = FALSE
  ))
  for (key in c("sim_file", "ld_blocks_dir", "output_dir")) {
    if (is.null(a[[key]])) stop(key, " must be specified", call. = FALSE)
  }

  sim_file      <- resolve_benchmark_path(a$sim_file)
  ld_blocks_dir <- resolve_benchmark_path(a$ld_blocks_dir)
  if (!file.exists(sim_file))     stop("sim_file does not exist: ", sim_file, call. = FALSE)
  if (!dir.exists(ld_blocks_dir)) stop("ld_blocks_dir does not exist: ", ld_blocks_dir, call. = FALSE)

  cat("\n========================================\n")
  cat("BLiP re-run: export SuSiE inputs\n")
  cat("========================================\n")
  cat("Sim file:      ", sim_file, "\n")
  cat("LD blocks:     ", ld_blocks_dir, "\n")
  cat("Output dir:    ", a$output_dir, "\n")

  sim <- readRDS(sim_file)
  k <- sim$parameters$K
  if (is.null(k)) {
    m <- regmatches(basename(sim_file), regexpr("_K[0-9]+", basename(sim_file)))
    if (length(m) == 0) stop("Cannot determine K from sim$parameters or the file name", call. = FALSE)
    k <- as.integer(sub("_K", "", m))
  }
  k <- as.integer(k)

  k_dir <- file.path(a$output_dir, sprintf("K%d", k))
  dir.create(k_dir, recursive = TRUE, showWarnings = FALSE)

  params_path <- file.path(k_dir, "sim_parameters.json")
  if (!file.exists(params_path)) {
    tryCatch(
      write_json(c(sim$parameters, list(sim_file = basename(sim_file))),
                 params_path, auto_unbox = TRUE, digits = NA, pretty = TRUE),
      error = function(e) warning("Could not write sim_parameters.json: ", conditionMessage(e)))
  }

  n_total <- length(sim$replicates)
  start   <- a$start_replicate
  end     <- if (is.na(a$end_replicate)) n_total else min(a$end_replicate, n_total)
  if (start < 1 || end < start) stop("Invalid replicate range ", start, "-", end, call. = FALSE)
  idx <- start:end

  cat(sprintf("K = %d, replicates %d-%d of %d, n = %d, n_cores = %d, pip_floor = %g, corr_set_purity = %g\n",
              k, start, end, n_total, a$n, a$n_cores, a$pip_floor, a$corr_set_purity))
  job_start <- Sys.time()

  results <- mclapply(idx, function(i) {
    tryCatch(
      export_replicate(sim$replicates[[i]], k, ld_blocks_dir, k_dir,
                       n_target = a$n, pip_floor = a$pip_floor,
                       corr_set_purity = a$corr_set_purity,
                       overwrite = a$overwrite, position = i),
      error = function(e) list(status = "error", rep_id = i, message = conditionMessage(e)))
  }, mc.cores = a$n_cores, mc.preschedule = FALSE)

  statuses <- vapply(results, function(r) {
    if (inherits(r, "try-error") || is.null(r$status)) "error" else r$status
  }, character(1))

  cat("\nStatus counts:\n")
  print(table(statuses))
  for (r in results[statuses == "error"]) {
    msg <- if (inherits(r, "try-error")) as.character(r) else r$message
    cat(sprintf("  ERROR rep %s: %s\n", r$rep_id, msg))
  }
  for (r in results[statuses == "missing_ld_block"]) {
    cat(sprintf("  MISSING LD BLOCK rep %s: %s\n", r$rep_id, r$ld_block_name))
  }

  cat(sprintf("\nDone in %.1f min.\n",
              as.numeric(difftime(Sys.time(), job_start, units = "mins"))))
  if (any(statuses %in% c("error", "missing_ld_block"))) quit(status = 1)
}

if (sys.nframe() == 0L) main()
