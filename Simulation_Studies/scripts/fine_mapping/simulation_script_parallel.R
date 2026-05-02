# =============================================================================
# Simulation Script — Fine-Mapping Method Evaluation (Parallel Version)
# =============================================================================
# Per-rep worker script invoked by command_generator_parallel.R via Rscript
# with command-line arguments. For each replicate index in
# [start_replicate, end_replicate], this script:
#
#   1. Loads an LD block (recycles through LD_blocks_dir),
#   2. Preprocesses the genotype matrix (mean imputation, MAF, variance
#      filtering, scaling, cap at 5000 variants),
#   3. Simulates the phenotype under one of two regimes:
#        - "sparse"  K independent causal variants with per-SNP h2
#        - "oligo"   sparse + oligogenic + infinitesimal h2 partition
#   4. Fits SuSiE, SuSiE-inf, and SuSiE-ash on the resulting (X, y),
#   5. Records per-method credible-set metrics (Power / FDR / Size).
#
# Output: one rds per invocation, written to /home/apm2217/output/, named:
#   sim_nrep<N>_<data_type>_<params>_rep<start>-<end>.rds
#
# Example HPC invocation:
#   Rscript simulation_script_parallel.R \
#     total_replicates=500 start_replicate=1 end_replicate=10 \
#     data_type='sparse' h2_per_snp=0.03 K=3 n=1000 \
#     LD_blocks_dir='oligo_LD_blocks'
# =============================================================================

library(susieR)
library(dplyr)
library(magrittr)
library(MASS)


# =============================================================================
# Helper functions: Phenotype simulation (multi-trait + sparse-K)
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

#' Select causal variants from genotype matrix
select_causal_variants <- function(X, n_causal, independent = FALSE, ld_threshold = 0.10) {
  if (!is.matrix(X)) stop("Input X must be a numeric matrix.")
  if (n_causal > ncol(X)) stop("n_causal cannot exceed number of SNPs in X.")
  if (n_causal == 1) return(sample(1:ncol(X), size = 1))

  LD_vars <- 1

  if (independent) {
    while (length(LD_vars) != 0) {
      vars <- sample(1:ncol(X), size = n_causal)
      cor_mat <- cor(X[, vars])
      LD_vars <- which(colSums(abs(cor_mat) > ld_threshold) > 1)
    }
  } else {
    while (length(LD_vars) != 0) {
      vars <- sample(1:ncol(X), size = n_causal)
      cor_mat <- cor(X[, vars])
      LD_vars <- which(colSums(abs(cor_mat) == 1) > 1)
    }
  }

  return(vars)
}

#' Simulate phenotype with controllable per-SNP heritability
simulate_phenotype <- function(X, n_causal = 5, h2_per_snp = 0.01, independent = TRUE) {
  causal_idx <- select_causal_variants(X, n_causal = n_causal, independent = independent)
  h2_total <- h2_per_snp * n_causal

  beta <- rep(0, ncol(X))
  beta[causal_idx] <- 1

  pheno <- sim_multi_traits(G = X, B = as.matrix(beta), h2g = h2_total, is_h2g_total = TRUE)
  y <- drop(pheno$P)

  var_y <- var(y)
  var_genetic <- var_y * h2_total
  var_residual <- var_y * (1 - h2_total)

  list(
    G = X,
    y = y,
    beta = beta,
    causal = causal_idx,
    h2_total = h2_total,
    h2_per_snp = h2_per_snp,
    residual_variance = var_residual
  )
}

# =============================================================================
# Helper functions: Effect-component simulators (sparse / oligo / infinitesimal)
# =============================================================================

#' Simulate sparse effects for eQTL data
simulate_sparse_effects <- function(G, h2_sparse, n_sparse, effect_sd = 0.5, valid_cols = NULL) {
  n_features <- ncol(G)
  beta <- rep(0, n_features)

  if (is.null(valid_cols)) {
    valid_cols <- 1:n_features
  }

  if (n_sparse > 0 && length(valid_cols) >= n_sparse) {
    sparse_indices <- sample(valid_cols, n_sparse)
    beta[sparse_indices] <- rnorm(n_sparse, mean = 0, sd = effect_sd)
  } else {
    sparse_indices <- integer(0)
  }

  if (length(sparse_indices) > 0) {
    sparse_effects <- as.vector(G[, sparse_indices, drop = FALSE] %*% beta[sparse_indices])
    current_var <- var(sparse_effects)
    if (current_var > 0) {
      scaling_factor <- sqrt(h2_sparse / current_var)
      beta[sparse_indices] <- beta[sparse_indices] * scaling_factor
    }
  }

  return(list(beta = beta, sparse_indices = sparse_indices))
}

#' Simulate oligogenic effects for eQTL data
simulate_oligogenic_effects <- function(G, h2_oligogenic, n_oligogenic, mixture_props,
                                        non_sparse_indices, effect_sds = c(0.05, 0.15)) {
  if (abs(sum(mixture_props) - 1) > 1e-6) {
    stop("mixture_props must sum to 1.")
  }

  n_features <- ncol(G)
  beta <- rep(0, n_features)

  n_available <- length(non_sparse_indices)
  n_oligogenic <- min(n_oligogenic, n_available)
  oligogenic_indices <- sample(non_sparse_indices, n_oligogenic, replace = FALSE)

  mixture_assignments <- sample(1:length(mixture_props), n_oligogenic, replace = TRUE, prob = mixture_props)
  beta[oligogenic_indices] <- rnorm(n_oligogenic, mean = 0, sd = effect_sds[mixture_assignments])

  oligogenic_effects <- as.vector(G[, oligogenic_indices] %*% beta[oligogenic_indices])
  current_var <- var(oligogenic_effects)

  if (current_var > 0) {
    scaling_factor <- sqrt(h2_oligogenic / current_var)
    beta[oligogenic_indices] <- beta[oligogenic_indices] * scaling_factor
  }

  mixture_assignments_full <- rep(NA, n_features)
  mixture_assignments_full[oligogenic_indices] <- mixture_assignments

  return(list(beta = beta, oligogenic_indices = oligogenic_indices,
              mixture_assignments = mixture_assignments_full))
}

#' Simulate infinitesimal effects for eQTL data
simulate_infinitesimal_effects <- function(G, h2_infinitesimal, infinitesimal_indices,
                                           effect_sd = 0.01) {
  n_features <- ncol(G)
  beta <- rep(0, n_features)
  n_inf <- length(infinitesimal_indices)

  if (n_inf > 0) {
    beta[infinitesimal_indices] <- rnorm(n_inf, mean = 0, sd = effect_sd)

    inf_effects <- as.vector(G[, infinitesimal_indices, drop = FALSE] %*% beta[infinitesimal_indices])
    current_var <- var(inf_effects)

    if (current_var > 0 && !is.na(current_var)) {
      scaling_factor <- sqrt(h2_infinitesimal / current_var)
      beta[infinitesimal_indices] <- beta[infinitesimal_indices] * scaling_factor
    }
  }
  return(beta)
}

# =============================================================================
# Helper functions: Causal-variant definition + oligogenic eQTL pipeline
# =============================================================================

#' Identify causal SNPs based on statistical power
is_causal_power <- function(G, beta, residual_variance, power = 0.80) {
  n <- nrow(G)
  p <- ncol(G)
  alpha <- 0.05 / p

  var_snp <- apply(G, 2, var)
  ncp <- n * (beta^2) * var_snp / residual_variance

  crit_val <- qchisq(alpha, df = 1, lower.tail = FALSE)
  power_per_snp <- pchisq(crit_val, df = 1, ncp = ncp, lower.tail = FALSE)

  causal_idx <- which(power_per_snp >= power)
  return(causal_idx)
}

#' Generate eQTL data with multiple genetic architecture components
generate_cis_qtl_data <- function(G, h2g = 0.25, prop_h2_sparse = 0.50,
                                  prop_h2_oligogenic = 0.35, prop_h2_infinitesimal = 0.15,
                                  n_sparse = 2, n_oligogenic = 5, n_inf = 15,
                                  mixture_props = c(0.75, 0.25), sparse_sd = 0.5,
                                  oligo_sds = c(0.05, 0.15), inf_sd = 0.01,
                                  standardize = TRUE, independent = TRUE,
                                  ld_threshold = 0.10, ld_step = 0.10, ld_max = 0.50,
                                  max_attempts = 200, seed = NULL) {
  if (abs(prop_h2_sparse + prop_h2_oligogenic + prop_h2_infinitesimal - 1) > 1e-6) {
    stop("The sum of prop_h2_sparse, prop_h2_oligogenic, and prop_h2_infinitesimal must equal 1.")
  }
  if (abs(sum(mixture_props) - 1) > 1e-6) {
    stop("mixture_props must sum to 1.")
  }

  n_samples <- nrow(G)
  n_features <- ncol(G)

  h2_sparse <- h2g * prop_h2_sparse
  h2_oligogenic <- h2g * prop_h2_oligogenic
  h2_infinitesimal <- h2g * prop_h2_infinitesimal

  valid_cols <- 1:n_features
  ld_satisfied <- FALSE
  current_ld_threshold <- ld_threshold
  total_attempts <- 0

  if (!is.null(seed)) set.seed(seed)

  while (!ld_satisfied && current_ld_threshold <= ld_max) {
    attempt <- 0

    while (!ld_satisfied && attempt < max_attempts) {
      attempt <- attempt + 1
      total_attempts <- total_attempts + 1

      if (!is.null(seed) && total_attempts > 1) {
        set.seed(seed * 1000 + total_attempts)
      }

      sparse_res <- simulate_sparse_effects(G, h2_sparse, n_sparse, sparse_sd, valid_cols)
      beta_sparse <- sparse_res$beta
      sparse_indices <- sparse_res$sparse_indices

      non_sparse_indices <- setdiff(valid_cols, sparse_indices)
      oligo_res <- simulate_oligogenic_effects(G, h2_oligogenic, n_oligogenic,
                                               mixture_props, non_sparse_indices, oligo_sds)
      beta_oligo <- oligo_res$beta
      oligogenic_indices <- oligo_res$oligogenic_indices

      infinitesimal_pool <- setdiff(non_sparse_indices, oligogenic_indices)

      if (!is.null(n_inf)) {
        n_inf_actual <- min(n_inf, length(infinitesimal_pool))
        infinitesimal_indices <- sample(infinitesimal_pool, n_inf_actual)
      } else {
        infinitesimal_indices <- infinitesimal_pool
      }

      beta_inf <- simulate_infinitesimal_effects(G, h2_infinitesimal, infinitesimal_indices, inf_sd)

      beta <- beta_sparse + beta_oligo + beta_inf

      g <- as.vector(G %*% beta)
      var_g <- var(g)
      var_epsilon <- var_g * (1 - h2g) / h2g
      epsilon <- rnorm(n_samples, 0, sqrt(var_epsilon))
      y <- g + epsilon

      if (independent) {
        ld_satisfied <- TRUE

        if (length(sparse_indices) > 1) {
          sparse_cor <- cor(G[, sparse_indices, drop = FALSE])
          high_ld_sparse <- which(abs(sparse_cor) >= current_ld_threshold &
                                    upper.tri(sparse_cor, diag = FALSE), arr.ind = TRUE)
          if (nrow(high_ld_sparse) > 0) ld_satisfied <- FALSE
        }

        if (ld_satisfied && length(oligogenic_indices) > 1) {
          oligo_cor <- cor(G[, oligogenic_indices, drop = FALSE])
          high_ld_oligo <- which(abs(oligo_cor) >= current_ld_threshold &
                                   upper.tri(oligo_cor, diag = FALSE), arr.ind = TRUE)
          if (nrow(high_ld_oligo) > 0) ld_satisfied <- FALSE
        }

        if (ld_satisfied && length(sparse_indices) > 0 && length(oligogenic_indices) > 0) {
          cross_cor <- cor(G[, sparse_indices, drop = FALSE],
                           G[, oligogenic_indices, drop = FALSE])
          high_ld_cross <- which(abs(cross_cor) >= current_ld_threshold, arr.ind = TRUE)
          if (nrow(high_ld_cross) > 0) ld_satisfied <- FALSE
        }
      } else {
        ld_satisfied <- TRUE
      }
    }

    if (!ld_satisfied && independent) {
      current_ld_threshold <- current_ld_threshold + ld_step
    }
  }

  if (independent && !ld_satisfied) {
    warning(paste0("Failed to satisfy LD constraints after ", total_attempts,
                   " total attempts (up to ld_threshold = ", ld_max,
                   "). Returning last generated data with LD violations."))
  }

  causal_indices <- is_causal_power(G, beta, var_epsilon, power = 0.80)

  var_y <- var(y)
  h2_sparse_actual <- var(as.vector(G[, sparse_indices] %*% beta[sparse_indices])) / var_y
  h2_oligogenic_actual <- var(as.vector(G[, oligogenic_indices] %*% beta[oligogenic_indices])) / var_y
  h2_infinitesimal_actual <- var(as.vector(G[, infinitesimal_indices] %*% beta_inf[infinitesimal_indices])) / var_y
  h2g_actual <- var(as.vector(G %*% beta)) / var_y

  return(list(
    G = G, y = y, beta = beta,
    h2g = h2g_actual, h2_sparse = h2_sparse_actual,
    h2_oligogenic = h2_oligogenic_actual, h2_infinitesimal = h2_infinitesimal_actual,
    sparse_indices = sparse_indices, oligogenic_indices = oligogenic_indices,
    infinitesimal_indices = infinitesimal_indices, residual_variance = var_epsilon,
    causal = causal_indices, ld_threshold_used = if (independent) current_ld_threshold else NA
  ))
}

# =============================================================================
# Helper functions: Data dispatcher (sparse vs oligo)
# =============================================================================

#' Generate simulation data based on data type
generate_data <- function(X, data_type, h2_per_snp = NULL, h2g = NULL, seed, K = 10, n = NULL,
                          prop_h2_sparse = 0.5, prop_h2_oligogenic = 0.35,
                          prop_h2_infinitesimal = 0.15, n_oligogenic = 5, n_inf = 15) {
  set.seed(seed)

  if (!is.null(n)) {
    n_available <- nrow(X)
    if (n > n_available) {
      warning("Requested n (", n, ") exceeds available samples (", n_available, "). Using all available samples.")
      n <- n_available
    }
    X <- X[1:n, , drop = FALSE]
  }

  if (data_type == "sparse") {
    if (is.null(h2_per_snp)) stop("h2_per_snp must be provided for sparse data type")
    data <- simulate_phenotype(X = X, n_causal = K, h2_per_snp = h2_per_snp, independent = TRUE)
  } else if (data_type == "oligo") {
    if (is.null(h2g)) stop("h2g must be provided for oligo data type")
    n_inf_use <- if (identical(n_inf, "all")) NULL else n_inf
    data <- generate_cis_qtl_data(G = X, h2g = h2g, n_sparse = K, seed = seed,
                                  prop_h2_sparse = prop_h2_sparse,
                                  prop_h2_oligogenic = prop_h2_oligogenic,
                                  prop_h2_infinitesimal = prop_h2_infinitesimal,
                                  n_oligogenic = n_oligogenic, n_inf = n_inf_use,
                                  standardize = FALSE, max_attempts = 500,
                                  ld_threshold = 0.1)
  } else {
    stop("Invalid data_type. Must be 'sparse' or 'oligo'.")
  }

  return(data)
}

# =============================================================================
# Helper functions: Method execution & evaluation
# =============================================================================

#' Run all fine-mapping methods on the data
run_all_methods <- function(data) {
  results <- list()

  susie_start <- Sys.time()
  susie_fit <- tryCatch({
    susie(X = data$G, y = data$y, unmappable_effects = 'none')
  }, error = function(e) { warning("SuSiE failed: ", e$message); NULL })
  results$susie <- list(fit = susie_fit,
                        elapsed_time = as.numeric(difftime(Sys.time(), susie_start, units = "secs")))

  susie_inf_start <- Sys.time()
  susie_inf_fit <- tryCatch({
    susie(X = data$G, y = data$y, unmappable_effects = 'inf')
  }, error = function(e) { warning("SuSiE-inf failed: ", e$message); NULL })
  results$susie_inf <- list(fit = susie_inf_fit,
                            elapsed_time = as.numeric(difftime(Sys.time(), susie_inf_start, units = "secs")))

  susie_ash_start <- Sys.time()
  susie_ash_fit <- tryCatch({
    susie(X = data$G, y = data$y, unmappable_effects = 'ash')
  }, error = function(e) { warning("SuSiE.ash failed: ", e$message); NULL })
  results$susie_ash <- list(fit = susie_ash_fit,
                            elapsed_time = as.numeric(difftime(Sys.time(), susie_ash_start, units = "secs")))

  return(results)
}

#' Calculate credible set metrics
calculate_cs_metrics <- function(fit, causal_indices, method) {
  if (method %in% c("susie", "susie_inf", "susie_ash")) {
    test.cs <- fit$sets$cs
  } else {
    stop("Unknown method specified for metric calculation: ", method)
  }

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

#' Evaluate performance metrics for all methods
evaluate_all_methods <- function(method_fits, causal_indices) {
  methods <- c("SuSiE", "SuSiE-inf", "SuSiE.ash")
  method_keys <- c("susie", "susie_inf", "susie_ash")

  metrics_list <- list()

  for (i in seq_along(methods)) {
    method_name <- methods[i]
    method_key <- method_keys[i]
    fit_obj <- method_fits[[method_key]]$fit

    if (is.null(fit_obj)) {
      metrics_list[[i]] <- data.frame(Method = method_name, Power = NA_real_,
                                      FDR = NA_real_, Size = NA_real_, stringsAsFactors = FALSE)
      next
    }

    tryCatch({
      cs_results <- calculate_cs_metrics(fit_obj, causal_indices, method = method_key)
      metrics_list[[i]] <- data.frame(Method = method_name, Power = cs_results$power,
                                      FDR = cs_results$fdr, Size = cs_results$size,
                                      stringsAsFactors = FALSE)
    }, error = function(e) {
      warning("Metric calculation failed for ", method_name, ": ", e$message)
      metrics_list[[i]] <- data.frame(Method = method_name, Power = NA_real_,
                                      FDR = NA_real_, Size = NA_real_, stringsAsFactors = FALSE)
    })
  }

  metrics_df <- do.call(rbind, metrics_list)
  return(metrics_df)
}

# =============================================================================
# Main Simulation Function
# =============================================================================

#' Main simulation function for parallel execution
#'
#' @param total_replicates Integer, total number of replicates in the full study (for filename)
#' @param start_replicate Integer, first replicate index to run in this job
#' @param end_replicate Integer, last replicate index to run in this job (if NULL, calculated from num_replicates)
#' @param num_replicates Integer, number of replicates to run in this job (alternative to end_replicate)
#' @param h2_per_snp Numeric, per-SNP heritability (only used for sparse data)
#' @param h2g Numeric, total SNP heritability (only used for oligo data)
#' @param data_type Character, either "sparse" or "oligo"
#' @param LD_blocks_dir Character, path to directory containing LD block RDS files
#' @param K Integer, number of sparse causal variants:
#'   - For sparse data: total number of independent causal variants
#'   - For oligo data: number of sparse effects (default = 10)
#' @param n Integer, number of samples to use from genotype matrix (default = NULL uses all samples)
#' @param prop_h2_sparse Numeric, proportion of h2g from sparse effects (only used for oligo data, default = 0.5)
#' @param prop_h2_oligogenic Numeric, proportion of h2g from oligogenic effects (only used for oligo data, default = 0.35)
#' @param prop_h2_infinitesimal Numeric, proportion of h2g from infinitesimal effects (only used for oligo data, default = 0.15)
#' @param n_oligogenic Integer, number of oligogenic effects (only used for oligo data, default = 5)
#' @param n_inf Integer or "all", number of infinitesimal effects (only used for oligo data, default = 15). Set to "all" for full infinitesimal background.
#'
#' @return List containing replicate results and simulation parameters
simulation_parallel <- function(total_replicates = NULL,
                               start_replicate = NULL,
                               end_replicate = NULL,
                               num_replicates = NULL,
                               h2_per_snp = NULL,
                               h2g = NULL,
                               data_type = NULL,
                               LD_blocks_dir = NULL,
                               K = NULL,
                               n = NULL,
                               prop_h2_sparse = NULL,
                               prop_h2_oligogenic = NULL,
                               prop_h2_infinitesimal = NULL,
                               n_oligogenic = NULL,
                               n_inf = NULL) {

  # Set default values
  if (is.null(start_replicate))   start_replicate <- 1
  if (is.null(h2_per_snp))        h2_per_snp <- 0.01
  if (is.null(h2g))               h2g <- 0.3
  if (is.null(data_type))         data_type <- "sparse"
  if (is.null(LD_blocks_dir))     LD_blocks_dir <- "LD_blocks"
  if (is.null(K))                 K <- 10
  if (is.null(prop_h2_sparse))         prop_h2_sparse <- 0.5
  if (is.null(prop_h2_oligogenic))     prop_h2_oligogenic <- 0.35
  if (is.null(prop_h2_infinitesimal))  prop_h2_infinitesimal <- 0.15
  if (is.null(n_oligogenic))           n_oligogenic <- 5
  if (is.null(n_inf))                  n_inf <- 15

  # Parse command-line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    for (arg in args) {
      split_arg <- strsplit(arg, "=")[[1]]
      key <- split_arg[1]
      value <- split_arg[2]
      value <- gsub("^'|'$", "", value)
      value <- gsub('^"|"$', "", value)

      if (key == "total_replicates") {
        total_replicates <- as.integer(value)
      } else if (key == "start_replicate") {
        start_replicate <- as.integer(value)
      } else if (key == "end_replicate") {
        end_replicate <- as.integer(value)
      } else if (key == "num_replicates") {
        num_replicates <- as.integer(value)
      } else if (key == "h2_per_snp") {
        h2_per_snp <- as.numeric(value)
      } else if (key == "h2g") {
        h2g <- as.numeric(value)
      } else if (key == "data_type") {
        data_type <- value
      } else if (key == "LD_blocks_dir") {
        LD_blocks_dir <- value
      } else if (key == "K") {
        K <- as.integer(value)
      } else if (key == "n") {
        n <- as.integer(value)
      } else if (key == "prop_h2_sparse") {
        prop_h2_sparse <- as.numeric(value)
      } else if (key == "prop_h2_oligogenic") {
        prop_h2_oligogenic <- as.numeric(value)
      } else if (key == "prop_h2_infinitesimal") {
        prop_h2_infinitesimal <- as.numeric(value)
      } else if (key == "n_oligogenic") {
        n_oligogenic <- as.integer(value)
      } else if (key == "n_inf") {
        if (tolower(value) == "all") {
          n_inf <- "all"
        } else {
          n_inf <- as.integer(value)
        }
      }
    }
  }

  # Calculate end_replicate if not provided
  if (is.null(end_replicate)) {
    if (is.null(num_replicates)) {
      stop("Either end_replicate or num_replicates must be provided")
    }
    end_replicate <- start_replicate + num_replicates - 1
  }

  if (is.null(num_replicates)) {
    num_replicates <- end_replicate - start_replicate + 1
  }

  if (is.null(total_replicates)) {
    total_replicates <- end_replicate
  }

  # Validate inputs
  if (!data_type %in% c("sparse", "oligo")) {
    stop("data_type must be either 'sparse' or 'oligo', got: ", data_type)
  }

  if (start_replicate < 1) {
    stop("start_replicate must be >= 1, got: ", start_replicate)
  }

  if (end_replicate < start_replicate) {
    stop("end_replicate must be >= start_replicate")
  }

  if (!dir.exists(LD_blocks_dir)) {
    stop("LD_blocks_dir does not exist: ", LD_blocks_dir)
  }

  # Setup output paths
  output_dir <- "/home/apm2217/output"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Create output filename
  base_name <- paste0("sim_nrep", total_replicates, "_", data_type)
  if (data_type == "sparse") {
    base_name <- paste0(base_name, "_h2persnp", h2_per_snp, "_K", K)
  } else if (data_type == "oligo") {
    base_name <- paste0(base_name, "_h2g", h2g, "_K", K,
                       "_pSparse", prop_h2_sparse,
                       "_pOligo", prop_h2_oligogenic,
                       "_pInf", prop_h2_infinitesimal,
                       "_nOligo", n_oligogenic,
                       "_nInf", n_inf)
  }

  if (!is.null(n)) {
    base_name <- paste0(base_name, "_n", n)
  }
  replicate_suffix <- paste0("_rep", start_replicate, "-", end_replicate)
  final_name <- paste0(base_name, replicate_suffix)

  # Find LD block files
  ld_block_files <- list.files(path = LD_blocks_dir, pattern = "\\.rds$", full.names = TRUE)

  if (length(ld_block_files) == 0) {
    stop("No .rds files found in LD_blocks_dir: ", LD_blocks_dir)
  }

  cat("Found", length(ld_block_files), "LD block files (will recycle if replicates exceed this).\n")

  cat("Starting simulation: replicates", start_replicate, "-", end_replicate,
      "(", data_type, ", K=", K, ")\n", sep = "")

  # Main simulation loop
  all_results <- vector("list", num_replicates)
  job_start_time <- Sys.time()

  for (i in start_replicate:end_replicate) {

    replicate_start <- Sys.time()
    result_index <- i - start_replicate + 1

    cat("Replicate", i, "/", total_replicates, "... ")

    ld_block_index <- ((i - 1) %% length(ld_block_files)) + 1
    ld_block <- readRDS(ld_block_files[ld_block_index])
    X_working <- ld_block$genotypes

    # Subset samples
    if (!is.null(n) && n < nrow(X_working)) {
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

    # Cap at 5000 variants
    if (ncol(X_working) > 5000) {
      X_working <- X_working[, 1:5000, drop = FALSE]
    }

    # Scale
    X_scaled <- scale(X_working)

    # Generate Data
    seed <- i + 1000
    data <- generate_data(X = X_scaled,
                          data_type = data_type,
                          h2_per_snp = h2_per_snp,
                          h2g = h2g,
                          seed = seed,
                          K = K,
                          n = n,
                          prop_h2_sparse = prop_h2_sparse,
                          prop_h2_oligogenic = prop_h2_oligogenic,
                          prop_h2_infinitesimal = prop_h2_infinitesimal,
                          n_oligogenic = n_oligogenic,
                          n_inf = n_inf)

    # Run methods and calculate metrics
    method_fits <- run_all_methods(data = data)
    metrics <- evaluate_all_methods(method_fits = method_fits,
                                    causal_indices = data$causal)

    replicate_elapsed <- as.numeric(difftime(Sys.time(), replicate_start, units = "secs"))
    all_results[[result_index]] <- list(
      replicate_id = i,
      ld_block_name = basename(ld_block_files[ld_block_index]),
      seed = seed,
      causal_indices = data$causal,
      beta = data$beta,
      residual_variance = data$residual_variance,
      metrics = metrics,
      elapsed_time = replicate_elapsed,
      fits = method_fits
    )

    cat("done (", round(replicate_elapsed, 1), "s)\n", sep = "")
    rm(ld_block, data, method_fits)
    gc()
  }

  job_elapsed <- as.numeric(difftime(Sys.time(), job_start_time, units = "secs"))
  results <- list(
    replicates = all_results,
    parameters = list(
      total_replicates = total_replicates,
      start_replicate = start_replicate,
      end_replicate = end_replicate,
      num_replicates = num_replicates,
      h2_per_snp = h2_per_snp,
      h2g = h2g,
      data_type = data_type,
      K = K,
      n = n,
      prop_h2_sparse = prop_h2_sparse,
      prop_h2_oligogenic = prop_h2_oligogenic,
      prop_h2_infinitesimal = prop_h2_infinitesimal,
      n_oligogenic = n_oligogenic,
      n_inf = n_inf,
      LD_blocks_dir = LD_blocks_dir,
      job_elapsed_time = job_elapsed
    )
  )

  final_file <- paste0(final_name, ".rds")
  output_path <- file.path(output_dir, final_file)

  saveRDS(results, output_path)
  cat("Completed in", round(job_elapsed / 60, 1), "min. Saved to:", output_path, "\n")

  return(results)
}

# Run the simulation (parameters will be overridden by command-line args if provided)
results <- simulation_parallel(total_replicates = NULL,
                              start_replicate = NULL,
                              end_replicate = NULL,
                              num_replicates = NULL,
                              h2_per_snp = NULL,
                              h2g = NULL,
                              data_type = NULL,
                              LD_blocks_dir = NULL,
                              K = NULL,
                              n = NULL,
                              n_inf = NULL)
