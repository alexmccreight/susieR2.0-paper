#!/usr/bin/env Rscript
# =============================================================================
# 3-Way Benchmark: susie_rss (lambda > 0)
# Comparison: susieR 1.0 (R) vs susieR 2.0 (R) vs susieR 2.0 (X)
# =============================================================================
#
# This script compares the runtime of three methods:
#   - susieR 1.0 (R): susie_rss_lambda(z, R, lambda) from stephenslab/susieR
#   - susieR 2.0 (R): susie_rss(z, R, lambda) from StatFunGen/susieR
#   - susieR 2.0 (X): susie_rss(z, X, lambda) from StatFunGen/susieR
#
# Uses sparse phenotype simulation with LD blocks.
# Supports parallel execution via start_replicate/end_replicate args.

# Load required libraries
suppressPackageStartupMessages({
  library(pkgload)
  library(MASS)
})

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Simulate multiple traits from genotype and effect sizes
sim_multi_traits <- function(G, B, h2g, is_h2g_total = TRUE, max_h2g = 1,
                             residual_corr = NULL, null_sigma = sqrt(0.1)) {
  if (!is_h2g_total) {
    max_causal <- max(apply(B, 2, function(x) sum(x != 0)))
    h2g <- min(h2g, max_h2g / max_causal)
  }

  P <- matrix(0, nrow = ncol(B), ncol = nrow(G))
  mu <- G %*% B
  sigma <- numeric(length = ncol(B))

  for (i in 1:ncol(mu)) {
    if (is_h2g_total) {
      if (sum(abs(mu[, i])) == 0) {
        sigma[i] <- null_sigma
      } else {
        sigma[i] <- sqrt(var(mu[, i]) * (1 - h2g) / h2g)
      }
    } else {
      if (sum(abs(mu[, i])) == 0) {
        sigma[i] <- null_sigma
      } else {
        first_index <- min(which(B[, i] == 1))
        if (var(G[, first_index]) / h2g - var(mu[, i]) >= 0) {
          sigma[i] <- sqrt(var(G[, first_index]) / h2g - var(mu[, i]))
        } else {
          stop("Per SNP heritability too large, residual variance will be less than 0.")
        }
      }
    }
  }

  if (is.null(residual_corr)) {
    residual_corr <- diag(length(sigma))
  }

  residual_var <- sweep(sweep(residual_corr, 2, sigma, "*"), 1, sigma, "*")
  P <- mu + mvrnorm(n = nrow(G), mu = rep(0, ncol(residual_var)), Sigma = residual_var)
  colnames(P) <- paste0("Trait_", 1:ncol(P))

  return(list(P = P, residual_var = residual_var))
}

#' Simulate sparse phenotype
simulate_phenotype_sparse <- function(X, n_causal, h2_per_snp, seed) {
  set.seed(seed)

  causal_idx <- select_causal_variants_simple(X, n_causal)

  h2_total <- h2_per_snp * n_causal

  beta <- rep(0, ncol(X))
  beta[causal_idx] <- 1

  pheno <- sim_multi_traits(
    G = X,
    B = as.matrix(beta),
    h2g = h2_total,
    is_h2g_total = TRUE
  )

  y <- drop(pheno$P)
  var_y <- var(y)
  var_residual <- var_y * (1 - h2_total)
  y_scaled <- as.vector(scale(y))

  list(
    y = y_scaled,
    beta = beta,
    causal = causal_idx,
    h2_total = h2_total,
    h2_per_snp = h2_per_snp,
    residual_variance = var_residual
  )
}

#' Select approximately independent causal variants
select_causal_variants_simple <- function(X, n_causal, ld_threshold = 0.15) {
  if (n_causal > ncol(X)) stop("n_causal > ncol(X)")
  if (n_causal == 1) return(sample(1:ncol(X), 1))

  LD_vars <- 1
  max_attempts <- 1000
  attempt <- 0

  while (length(LD_vars) != 0 && attempt < max_attempts) {
    vars <- sample(1:ncol(X), size = n_causal)
    cor_mat <- cor(X[, vars])
    LD_vars <- which(colSums(abs(cor_mat) > ld_threshold) > 1)
    attempt <- attempt + 1
  }

  if (attempt >= max_attempts) {
    warning("Could not find fully independent variants. Using best attempt.")
  }

  return(vars)
}

# =============================================================================
# MAIN SIMULATION FUNCTION
# =============================================================================

simulation <- function(total_replicates = NULL,
                      start_replicate = NULL,
                      end_replicate = NULL,
                      h2_per_snp = NULL,
                      K = NULL,
                      n = NULL,
                      p = NULL,
                      LD_blocks_dir = NULL,
                      lambda = 1e-5,
                      L = 10,
                      checkpoint_every = 5,
                      resume = TRUE,
                      tolerance = 1e-5) {

  # ---- Set default values ----
  if (is.null(total_replicates)) total_replicates <- 100
  if (is.null(h2_per_snp))      h2_per_snp <- 0.03
  if (is.null(K))               K <- 3
  if (is.null(n))               n <- 500
  if (is.null(p))               p <- 5000
  if (is.null(LD_blocks_dir))   LD_blocks_dir <- "oligo_LD_blocks"

  # ---- Parse command-line arguments ----
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    for (arg in args) {
      split_arg <- strsplit(arg, "=")[[1]]
      key <- split_arg[1]
      value <- split_arg[2]

      # Remove quotes if present
      value <- gsub("^'|'$", "", value)
      value <- gsub('^"|"$', "", value)

      if (key == "total_replicates") {
        total_replicates <- as.integer(value)
      } else if (key == "start_replicate") {
        start_replicate <- as.integer(value)
      } else if (key == "end_replicate") {
        end_replicate <- as.integer(value)
      } else if (key == "h2_per_snp") {
        h2_per_snp <- as.numeric(value)
      } else if (key == "K") {
        K <- as.integer(value)
      } else if (key == "n") {
        n <- as.integer(value)
      } else if (key == "p") {
        p <- as.integer(value)
      } else if (key == "LD_blocks_dir") {
        LD_blocks_dir <- value
      } else if (key == "lambda") {
        lambda <- as.numeric(value)
      } else if (key == "L") {
        L <- as.integer(value)
      } else if (key == "checkpoint_every") {
        checkpoint_every <- as.integer(value)
      } else if (key == "resume") {
        resume <- as.logical(value)
      } else if (key == "tolerance") {
        tolerance <- as.numeric(value)
      }
    }
  }

  # ---- Set replicate range ----
  if (is.null(start_replicate)) start_replicate <- 1
  if (is.null(end_replicate))   end_replicate <- total_replicates
  num_replicates <- end_replicate - start_replicate + 1

  # ---- Validate inputs ----
  if (!dir.exists(LD_blocks_dir)) {
    stop("LD_blocks_dir does not exist: ", LD_blocks_dir)
  }
  if (start_replicate > end_replicate) {
    stop("start_replicate (", start_replicate, ") > end_replicate (", end_replicate, ")")
  }

  # ---- Load susieR packages ----
  cat("\n========================================\n")
  cat("LOADING SUSIER PACKAGES\n")
  cat("========================================\n")

  # Use PID-specific paths so parallel jobs on the same node don't collide
  pid <- Sys.getpid()

  # Helper to clone and checkout a susieR package
  setup_package <- function(path, repo_url, commit, label) {
    if (dir.exists(path)) {
      # Check if it's a valid git repo with the correct commit
      current <- tryCatch(
        system(sprintf("cd %s && git rev-parse HEAD 2>/dev/null", path), intern = TRUE),
        warning = function(w) character(0)
      )

      if (length(current) == 1 && grepl(paste0("^", commit), current)) {
        cat(label, "already at correct commit\n")
        return(invisible(NULL))
      }

      # Directory exists but is corrupted or wrong commit — wipe it
      cat(label, "directory exists but needs fresh clone. Removing...\n")
      system(sprintf("rm -rf %s", path))
    }

    cat("Cloning", label, "...\n")
    result <- system(sprintf("git clone -q %s %s", repo_url, path))
    if (result != 0) stop("Failed to clone ", label)
    result <- system(sprintf("cd %s && git checkout -q %s", path, commit))
    if (result != 0) stop("Failed to checkout ", label, " commit ", commit)
  }

  # Reference package (susieR 1.0)
  REF_PATH <- sprintf("/tmp/susieR_reference_%d", pid)
  REF_COMMIT <- "1f9166c"
  setup_package(REF_PATH, "https://github.com/stephenslab/susieR.git", REF_COMMIT, "Reference (susieR 1.0)")

  ref_env <- pkgload::load_all(REF_PATH, export_all = FALSE, quiet = TRUE)
  cat("Reference package loaded (susieR 1.0 @", REF_COMMIT, ")\n")

  # Development package (susieR 2.0)
  DEV_PATH <- sprintf("/tmp/susieR_development_%d", pid)
  DEV_COMMIT <- "84616edd3469f741d75be7a64f10530ce6c87976"
  setup_package(DEV_PATH, "https://github.com/stephenslab/susieR.git", DEV_COMMIT, "Development (susieR 2.0)")

  dev_env <- pkgload::load_all(DEV_PATH, export_all = FALSE, quiet = TRUE)
  cat("Development package loaded (susieR 2.0 @", DEV_COMMIT, ")\n\n")

  # ---- Setup output paths ----
  output_dir <- "/home/apm2217/output"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Create base filename
  base_name <- sprintf("benchmark_n%d_p%d_lambda%s_K%d_h2persnp%s_nrep%d",
                       n, p, lambda, K, h2_per_snp, total_replicates)

  checkpoint_dir <- file.path(output_dir, "checkpoints")
  if (!dir.exists(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE)
  }

  # ---- Find LD block files ----
  ld_block_files <- sort(list.files(path = LD_blocks_dir, pattern = "\\.rds$", full.names = TRUE))

  if (length(ld_block_files) == 0) {
    stop("No .rds files found in LD_blocks_dir: ", LD_blocks_dir)
  }

  cat("Found", length(ld_block_files), "LD block files\n")

  # ---- Check for existing checkpoints and resume ----
  all_results <- vector("list", num_replicates)
  actual_start <- start_replicate

  if (resume) {
    # Build checkpoint pattern for this specific job range
    job_base <- sprintf("%s_rep%d-%d", base_name, start_replicate, end_replicate)
    checkpoint_pattern <- paste0(job_base, "_checkpoint_rep")
    checkpoint_files <- list.files(checkpoint_dir,
                                  pattern = paste0(gsub("\\.", "\\\\.", checkpoint_pattern), ".*\\.rds$"),
                                  full.names = TRUE)

    if (length(checkpoint_files) > 0) {
      cat("\n========================================\n")
      cat("CHECKPOINT FOUND - Resuming simulation\n")
      cat("========================================\n")

      for (ckpt_file in checkpoint_files) {
        ckpt_data <- readRDS(ckpt_file)
        for (rep_result in ckpt_data$replicates) {
          if (!is.null(rep_result)) {
            rep_id <- rep_result$replicate_id
            result_index <- rep_id - start_replicate + 1
            all_results[[result_index]] <- rep_result
            cat("Loaded checkpoint for replicate", rep_id, "\n")
          }
        }
      }

      completed <- sapply(all_results, function(x) !is.null(x))
      if (all(completed)) {
        cat("All replicates already completed!\n")
        actual_start <- end_replicate + 1
      } else {
        # Find first incomplete in global replicate numbering
        first_incomplete_idx <- which(!completed)[1]
        actual_start <- start_replicate + first_incomplete_idx - 1
        cat("Resuming from replicate", actual_start, "\n")
      }
      cat("========================================\n\n")
    }
  }

  # ---- Print simulation configuration ----
  if (actual_start <= end_replicate) {
    cat("\n========================================\n")
    cat("3-WAY BENCHMARK CONFIGURATION\n")
    cat("========================================\n")
    cat("Replicate range:", start_replicate, "to", end_replicate,
        "(", num_replicates, "replicates)\n")
    cat("Per-SNP heritability:", h2_per_snp, "\n")
    cat("Number of causal variants (K):", K, "\n")
    cat("Sample size (n):", n, "\n")
    cat("Number of variants (p):", p, "\n")
    cat("Lambda:", lambda, "\n")
    cat("Number of single effects (L):", L, "\n")
    cat("LD blocks directory:", LD_blocks_dir, "\n")
    cat("Methods:\n")
    cat("  1. susieR 1.0 (R): susie_rss_lambda(z, R, lambda)\n")
    cat("  2. susieR 2.0 (R): susie_rss(z, R, lambda)\n")
    cat("  3. susieR 2.0 (X): susie_rss(z, X, lambda)\n")
    cat("Starting from replicate:", actual_start, "\n")
    cat("========================================\n\n")
  }

  # ---- Main simulation loop ----
  job_start_time <- Sys.time()

  for (i in actual_start:end_replicate) {

    replicate_start <- Sys.time()
    result_index <- i - start_replicate + 1

    cat("========================================\n")
    cat("Replicate", i, "(", result_index, "of", num_replicates, ")\n")

    # Cycle through LD blocks
    ld_block_index <- ((i - 1) %% length(ld_block_files)) + 1
    cat("Using LD block:", basename(ld_block_files[ld_block_index]), "\n")
    cat("========================================\n")

    # Load LD block
    ld_block <- readRDS(ld_block_files[ld_block_index])
    X_working <- ld_block$genotypes

    # Subset samples
    if (n < nrow(X_working)) {
      X_working <- X_working[1:n, , drop = FALSE]
    }

    # Mean imputation for NAs
    col_means <- colMeans(X_working, na.rm = TRUE)
    na_idx <- which(is.na(X_working), arr.ind = TRUE)
    if (nrow(na_idx) > 0) {
      X_working[na_idx] <- col_means[na_idx[, 2]]
    }

    # MAF filter (>= 1%)
    maf <- colMeans(X_working) / 2
    maf <- pmin(maf, 1 - maf)
    X_working <- X_working[, maf >= 0.01, drop = FALSE]

    # Remove zero-variance columns
    col_vars <- apply(X_working, 2, var)
    X_working <- X_working[, col_vars > 0, drop = FALSE]

    # Cap at p variants
    if (ncol(X_working) > p) {
      X_working <- X_working[, 1:p, drop = FALSE]
    }

    n_use <- nrow(X_working)
    cat("Using genotype matrix: n=", n_use, ", p=", ncol(X_working), "\n")

    # ---- Generate phenotype data ----
    seed <- i + 1000
    cat("Generating phenotype (seed=", seed, ")...\n")

    data <- simulate_phenotype_sparse(
      X = X_working,
      n_causal = K,
      h2_per_snp = h2_per_snp,
      seed = seed
    )

    cat("Phenotype generated. Causal variants:", paste(data$causal, collapse = ", "), "\n")

    # ---- Compute summary statistics (outside timing) ----
    cat("Computing summary statistics...\n")
    input_ss <- dev_env$env$compute_suff_stat(X_working, data$y, standardize = TRUE)
    ss <- dev_env$env$univariate_regression(X_working, data$y)

    # Validate summary statistics
    if (any(is.na(input_ss$XtX)) || any(is.infinite(input_ss$XtX))) {
      stop("XtX contains NA or Inf values")
    }
    if (any(is.na(ss$betahat)) || any(is.infinite(ss$betahat))) {
      stop("betahat contains NA or Inf values")
    }
    if (any(ss$sebetahat == 0)) {
      stop("sebetahat contains zero values")
    }

    R <- with(input_ss, cov2cor(XtX))
    R <- (R + t(R)) / 2  # Ensure symmetry

    if (any(is.na(R)) || any(is.infinite(R))) {
      stop("R matrix contains NA or Inf values")
    }

    z <- with(ss, betahat / sebetahat)

    # Scale X for the X-path
    # scale() gives colSums(X^2) = n-1; multiply by sqrt(n/(n-1)) to get colSums = n
    # susie_rss internally rescales by sqrt((n-1)/n), giving csd=1 matching the R-path
    X_scaled <- scale(X_working) * sqrt(n_use / (n_use - 1))

    # ---- Run the three methods ----

    # Method 1: susieR 1.0 (R)
    cat("Running susieR 1.0 (R)...\n")
    time_start <- Sys.time()
    fit_1.0R <- tryCatch(
      suppressWarnings(ref_env$env$susie_rss_lambda(
        z = z, R = R, L = L, lambda = lambda,
        estimate_residual_variance = TRUE,
        estimate_prior_method = "optim",
        verbose = FALSE
      )),
      error = function(e) {
        cat("  ERROR:", e$message, "\n")
        NULL
      }
    )
    time_1.0R <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
    cat("  susieR 1.0 (R) completed in", round(time_1.0R, 3), "seconds\n")

    # Method 2: susieR 2.0 (R)
    cat("Running susieR 2.0 (R)...\n")
    time_start <- Sys.time()
    fit_2.0R <- tryCatch(
      suppressWarnings(dev_env$env$susie_rss(
        z = z, R = R, L = L, lambda = lambda,
        estimate_residual_variance = TRUE,
        estimate_prior_method = "optim",
        verbose = FALSE
      )),
      error = function(e) {
        cat("  ERROR:", e$message, "\n")
        NULL
      }
    )
    time_2.0R <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
    cat("  susieR 2.0 (R) completed in", round(time_2.0R, 3), "seconds\n")

    # Method 3: susieR 2.0 (X)
    cat("Running susieR 2.0 (X)...\n")
    time_start <- Sys.time()
    fit_2.0X <- tryCatch(
      suppressWarnings(dev_env$env$susie_rss(
        z = z, X = X_scaled, n = n_use, L = L, lambda = lambda,
        estimate_residual_variance = TRUE,
        estimate_prior_method = "optim",
        verbose = FALSE
      )),
      error = function(e) {
        cat("  ERROR:", e$message, "\n")
        NULL
      }
    )
    time_2.0X <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
    cat("  susieR 2.0 (X) completed in", round(time_2.0X, 3), "seconds\n")

    # ---- Sanity checks ----
    pip_diff_1.0R_2.0R <- NA
    pip_diff_2.0R_2.0X <- NA

    if (!is.null(fit_1.0R) && !is.null(fit_2.0R)) {
      pip_diff_1.0R_2.0R <- max(abs(fit_1.0R$pip - fit_2.0R$pip))
    }
    if (!is.null(fit_2.0R) && !is.null(fit_2.0X)) {
      pip_diff_2.0R_2.0X <- max(abs(fit_2.0R$pip - fit_2.0X$pip))
    }

    # ---- Store results ----
    replicate_elapsed <- as.numeric(difftime(Sys.time(), replicate_start, units = "secs"))

    all_results[[result_index]] <- list(
      replicate_id = i,
      n = nrow(X_working),
      p = ncol(X_working),
      ld_block_name = basename(ld_block_files[ld_block_index]),
      seed = seed,
      time_1.0_R = time_1.0R,
      time_2.0_R = time_2.0R,
      time_2.0_X = time_2.0X,
      speedup_2.0R_vs_1.0R = time_1.0R / time_2.0R,
      speedup_2.0X_vs_1.0R = time_1.0R / time_2.0X,
      speedup_2.0X_vs_2.0R = time_2.0R / time_2.0X,
      pip_diff_1.0R_2.0R = pip_diff_1.0R_2.0R,
      pip_diff_2.0R_2.0X = pip_diff_2.0R_2.0X,
      niter_1.0R = if (!is.null(fit_1.0R)) fit_1.0R$niter else NA,
      niter_2.0R = if (!is.null(fit_2.0R)) fit_2.0R$niter else NA,
      niter_2.0X = if (!is.null(fit_2.0X)) fit_2.0X$niter else NA,
      converged_1.0R = if (!is.null(fit_1.0R)) fit_1.0R$converged else NA,
      converged_2.0R = if (!is.null(fit_2.0R)) fit_2.0R$converged else NA,
      converged_2.0X = if (!is.null(fit_2.0X)) fit_2.0X$converged else NA,
      elapsed_time = replicate_elapsed
    )

    # Print replicate summary
    cat(sprintf("\nReplicate %d summary:\n", i))
    cat(sprintf("  1.0(R): %.2fs  |  2.0(R): %.2fs  |  2.0(X): %.2fs\n",
                time_1.0R, time_2.0R, time_2.0X))
    cat(sprintf("  Speedup 2.0(R) vs 1.0(R): %.1fx\n", time_1.0R / time_2.0R))
    cat(sprintf("  Speedup 2.0(X) vs 1.0(R): %.1fx\n", time_1.0R / time_2.0X))
    if (!is.na(pip_diff_2.0R_2.0X)) {
      cat(sprintf("  PIP diff 2.0(R) vs 2.0(X): %.6f\n", pip_diff_2.0R_2.0X))
    }
    cat("\n")

    # Clean up memory
    rm(ld_block, X_working, X_scaled, data,
       fit_1.0R, fit_2.0R, fit_2.0X, R, z, ss, input_ss)
    gc(verbose = FALSE)

    # ---- Save checkpoint if needed ----
    if (result_index %% checkpoint_every == 0 || i == end_replicate) {
      job_base <- sprintf("%s_rep%d-%d", base_name, start_replicate, end_replicate)
      checkpoint_file <- file.path(checkpoint_dir,
                                  paste0(job_base, "_checkpoint_rep", i, ".rds"))

      checkpoint_data <- list(
        replicates = all_results[1:result_index],
        parameters = list(
          total_replicates = total_replicates,
          start_replicate = start_replicate,
          end_replicate = end_replicate,
          h2_per_snp = h2_per_snp,
          K = K,
          n = n,
          p = p,
          lambda = lambda,
          L = L,
          tolerance = tolerance,
          LD_blocks_dir = LD_blocks_dir,
          checkpoint_replicate = i
        )
      )

      cat(">>> Saving checkpoint at replicate", i, "\n")
      saveRDS(checkpoint_data, checkpoint_file)

      # Delete old checkpoints for this job
      if (file.exists(checkpoint_file)) {
        old_checkpoints <- list.files(checkpoint_dir,
                                      pattern = paste0(gsub("\\.", "\\\\.", job_base),
                                                       "_checkpoint_rep.*\\.rds$"),
                                      full.names = TRUE)
        old_checkpoints <- setdiff(old_checkpoints, checkpoint_file)
        if (length(old_checkpoints) > 0) {
          file.remove(old_checkpoints)
        }
      }
    }
  }

  # ---- Compile final results ----
  job_elapsed <- as.numeric(difftime(Sys.time(), job_start_time, units = "secs"))

  results <- list(
    replicates = all_results,
    parameters = list(
      total_replicates = total_replicates,
      start_replicate = start_replicate,
      end_replicate = end_replicate,
      h2_per_snp = h2_per_snp,
      K = K,
      n = n,
      p = p,
      lambda = lambda,
      L = L,
      tolerance = tolerance,
      LD_blocks_dir = LD_blocks_dir,
      ref_commit = REF_COMMIT,
      dev_commit = DEV_COMMIT,
      job_elapsed_time = job_elapsed
    )
  )

  # ---- Save final results ----
  final_file <- sprintf("%s_rep%d-%d.rds", base_name, start_replicate, end_replicate)
  output_path <- file.path(output_dir, final_file)

  cat("========================================\n")
  cat("Saving final results to:", output_path, "\n")
  cat("========================================\n")

  saveRDS(results, output_path)

  # ---- Clean up checkpoint files for this job ----
  job_base <- sprintf("%s_rep%d-%d", base_name, start_replicate, end_replicate)
  checkpoint_files <- list.files(checkpoint_dir,
                                pattern = paste0(gsub("\\.", "\\\\.", job_base),
                                                 "_checkpoint_rep.*\\.rds$"),
                                full.names = TRUE)
  if (length(checkpoint_files) > 0) {
    file.remove(checkpoint_files)
  }

  # ---- Print summary ----
  cat("\n========================================\n")
  cat("3-WAY BENCHMARK SUMMARY\n")
  cat("========================================\n")

  times_1.0R <- sapply(all_results, function(x) if (!is.null(x)) x$time_1.0_R else NA)
  times_2.0R <- sapply(all_results, function(x) if (!is.null(x)) x$time_2.0_R else NA)
  times_2.0X <- sapply(all_results, function(x) if (!is.null(x)) x$time_2.0_X else NA)

  cat("Replicates:", start_replicate, "to", end_replicate, "\n")
  cat("n =", n, ", p =", p, "\n\n")

  cat("Runtime susieR 1.0 (R):\n")
  cat("  Mean:", round(mean(times_1.0R, na.rm = TRUE), 3), "s\n")
  cat("  Median:", round(median(times_1.0R, na.rm = TRUE), 3), "s\n\n")

  cat("Runtime susieR 2.0 (R):\n")
  cat("  Mean:", round(mean(times_2.0R, na.rm = TRUE), 3), "s\n")
  cat("  Median:", round(median(times_2.0R, na.rm = TRUE), 3), "s\n\n")

  cat("Runtime susieR 2.0 (X):\n")
  cat("  Mean:", round(mean(times_2.0X, na.rm = TRUE), 3), "s\n")
  cat("  Median:", round(median(times_2.0X, na.rm = TRUE), 3), "s\n\n")

  cat("Speedup 2.0(R) vs 1.0(R): ", round(mean(times_1.0R / times_2.0R, na.rm = TRUE), 2), "x\n")
  cat("Speedup 2.0(X) vs 1.0(R): ", round(mean(times_1.0R / times_2.0X, na.rm = TRUE), 2), "x\n")
  cat("Speedup 2.0(X) vs 2.0(R): ", round(mean(times_2.0R / times_2.0X, na.rm = TRUE), 2), "x\n")

  cat("\nJob elapsed time:", round(job_elapsed / 60, 2), "minutes\n")
  cat("========================================\n")

  cat("\nBenchmark completed successfully!\n")

  return(results)
}

# =============================================================================
# RUN SIMULATION
# =============================================================================

results <- simulation()
