#!/usr/bin/env Rscript

# Supplementary Figure S4: worker for the cross-validated prediction benchmark.
# For each replicate i in [start_replicate, end_replicate]:
#   1. take LD block i (the sorted list of blocks is recycled), keep the first
#      n samples, mean-impute missing genotypes, drop variants with MAF < 0.01
#      or zero variance, and keep the first 5000 of the rest;
#   2. simulate the phenotype with seed i + 1000 under data_type "sparse"
#      (K causal variants, each with the same h2) or "oligo" (sparse,
#      oligogenic and infinitesimal components sharing h2g);
#   3. assign the samples to n_folds folds and fit SuSiE, SuSiE-inf and
#      SuSiE-ash (L = 10, max_iter = 100) on each training fold;
#   4. predict the held-out fold from coef(fit) and record RMSE and
#      R^2 = cor(y, y_hat)^2; the replicate's value is the mean over folds.
#
# The phenotype is simulated on the standardized genotype matrix, which defines
# the truth. The models are fit on the unscaled filtered genotypes of the
# training fold only: susie() standardizes within the fold and coef() returns
# effects on the input scale, so no test-fold statistic enters a fit. The MAF
# filter, the imputation and the 5000-variant cap use all samples; they do not
# involve the phenotype. Block order, preprocessing, seed and simulation
# functions are those of fine_mapping/simulation_script_parallel.R, so
# replicate i of a scenario has the same genotypes and phenotype in both
# benchmarks.
#
# Usage: Rscript prediction_cv_script_parallel.R total_replicates=300 \
#          start_replicate=1 end_replicate=2 data_type='sparse' h2_per_snp=0.03 \
#          K=1 n=1000 n_folds=5 convergence_method='pip' \
#          LD_blocks_dir='oligo_LD_blocks' output_dir='/home/apm2217/output'
#        Every argument is key=value; the full list is above
#        prediction_cv_parallel().
#
# Input:  <LD_blocks_dir>/*.rds, each a list with a $genotypes matrix
# Output: <output_dir>/pred_cv_nrep<N>_<data_type>_<parameters>[_n<n>]_nfolds<F>
#           _rep<start>-<end>.rds  (list: $replicates, $parameters)

library(susieR)
library(dplyr)
library(magrittr)
library(MASS)

# ---- Phenotype simulation ----

# The functions of this section are identical to those of the same names in
# fine_mapping/simulation_script_parallel.R.

# Phenotypes G %*% B + Gaussian noise, one column of B per trait. With
# is_h2g_total = TRUE (the only mode used here) the residual variance is
# var(G b) * (1 - h2g) / h2g, so that the genetic share of the variance is h2g.
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

# Draws n_causal variants. With independent = TRUE, redraws until no pair has
# |r| > ld_threshold; otherwise only until no pair is perfectly correlated.
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

# Sparse architecture: n_causal variants with effect 1 on the scale of X, and
# total h2 = h2_per_snp * n_causal.
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

# The three component simulators of the oligo architecture draw normal effects
# and rescale them so that the variance of the component's genetic value
# equals its h2 argument.

# n_sparse effects with standard deviation effect_sd, on variants drawn from
# valid_cols.
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

# n_oligogenic effects on variants drawn from non_sparse_indices, from a
# normal mixture with standard deviations effect_sds and weights mixture_props.
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

# One effect with standard deviation effect_sd on every variant in
# infinitesimal_indices.
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

# Causal variants of the oligo architecture: those a marginal test at the
# Bonferroni level 0.05 / p detects with probability >= `power`, given the
# effect, the variance of the variant and the residual variance.
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

# Oligo architecture: beta is the sum of sparse, oligogenic and infinitesimal
# effects, with h2g split by the prop_h2_* arguments, and y = G beta + noise
# with residual variance var(G beta) * (1 - h2g) / h2g. n_inf = NULL gives an
# infinitesimal effect to every variant that is neither sparse nor oligogenic.
# With independent = TRUE the draw is repeated until every pair among the
# sparse and oligogenic variants has |r| < ld_threshold; after max_attempts
# failed draws the threshold rises by ld_step, up to ld_max. Each draw after
# the first is seeded with seed * 1000 + its number. `standardize` is not used.
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

# One simulated data set: set.seed(seed), then the sparse or the oligo
# architecture. n_inf = "all" is passed on as NULL.
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

# ---- Prediction ----

# Test-set predictions come from coef.susie(), which handles
# X_column_scale_factors and theta (the unmappable effects of SuSiE-inf and
# SuSiE-ash); predict.susie() returns fitted values on the training set and
# cannot be used here.
# This function is not called: the cross-validation loop uses
# get_predictions_with_padding(), defined there.
get_predictions <- function(fit, X_test) {
  beta_with_intercept <- coef(fit)
  intercept <- beta_with_intercept[1]
  beta <- beta_with_intercept[-1]

  y_pred <- as.vector(intercept + X_test %*% beta)

  return(y_pred)
}

calculate_rmse <- function(y_true, y_pred) {
  sqrt(mean((y_true - y_pred)^2))
}

# ---- Worker ----

# Runs the replicates of one job and saves them. The arguments are read from
# the command line as key=value (defaults in brackets):
#   total_replicates    replicates in the whole study; goes into the output
#                       file name [end_replicate]
#   start_replicate     first replicate of this job [1]
#   end_replicate       last replicate of this job (or give num_replicates)
#   data_type           "sparse" or "oligo" ["sparse"]
#   h2_per_snp          sparse: h2 of each causal variant [0.03]
#   h2g                 oligo: total h2 [0.3]
#   K                   sparse: number of causal variants; oligo: number of
#                       sparse effects [10]
#   n                   number of samples, taken from the top of the genotype
#                       matrix [all]
#   prop_h2_sparse, prop_h2_oligogenic, prop_h2_infinitesimal
#                       oligo: shares of h2g [0.5, 0.35, 0.15]
#   n_oligogenic        oligo: number of oligogenic effects [5]
#   n_inf               oligo: number of infinitesimal effects, or "all" [15]
#   n_folds             cross-validation folds [5]
#   convergence_method  for SuSiE-inf and SuSiE-ash ["pip"]
#   LD_blocks_dir       ["LD_blocks"]
#   output_dir          ["/home/apm2217/output"]
#   skip_existing       TRUE: exit at once if the output file exists [FALSE]
prediction_cv_parallel <- function(total_replicates = NULL,
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
                                   n_inf = NULL,
                                   n_folds = 5,
                                   convergence_method = "pip",
                                   output_dir = NULL,
                                   skip_existing = FALSE) {

  # ---- Defaults ----
  if (is.null(start_replicate))   start_replicate <- 1
  if (is.null(output_dir))        output_dir <- "/home/apm2217/output"
  if (is.null(h2_per_snp))        h2_per_snp <- 0.03
  if (is.null(h2g))               h2g <- 0.3
  if (is.null(data_type))         data_type <- "sparse"
  if (is.null(LD_blocks_dir))     LD_blocks_dir <- "LD_blocks"
  if (is.null(K))                 K <- 10
  if (is.null(prop_h2_sparse))         prop_h2_sparse <- 0.5
  if (is.null(prop_h2_oligogenic))     prop_h2_oligogenic <- 0.35
  if (is.null(prop_h2_infinitesimal))  prop_h2_infinitesimal <- 0.15
  if (is.null(n_oligogenic))           n_oligogenic <- 5
  if (is.null(n_inf))                  n_inf <- 15
  # n stays NULL unless given: all samples are used.

  # ---- Command-line arguments ----
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
      } else if (key == "n_folds") {
        n_folds <- as.integer(value)
      } else if (key == "convergence_method") {
        convergence_method <- value
      } else if (key == "output_dir") {
        output_dir <- value
      } else if (key == "skip_existing") {
        skip_existing <- as.logical(value)
      }
    }
  }

  # ---- Replicate range ----
  if (is.null(end_replicate)) {
    if (is.null(num_replicates)) {
      stop("Either end_replicate or num_replicates must be provided")
    }
    end_replicate <- start_replicate + num_replicates - 1
  }

  if (is.null(num_replicates)) {
    num_replicates <- end_replicate - start_replicate + 1
  }

  # Without total_replicates, this job is taken to be the whole study.
  if (is.null(total_replicates)) {
    total_replicates <- end_replicate
  }

  # ---- Checks ----
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

  # ---- Output ----
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # base_name identifies the parameter set; command_generator_parallel_cv.R
  # builds the same string for the combine commands.
  base_name <- paste0("pred_cv_nrep", total_replicates, "_", data_type)
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

  base_name <- paste0(base_name, "_nfolds", n_folds)

  replicate_suffix <- paste0("_rep", start_replicate, "-", end_replicate)
  final_name <- paste0(base_name, replicate_suffix)

  # Resume support: a finished job is not repeated.
  if (skip_existing && file.exists(file.path(output_dir, paste0(final_name, ".rds")))) {
    cat("Output already exists, skipping:", file.path(output_dir, paste0(final_name, ".rds")), "\n")
    return(invisible(NULL))
  }

  # ---- LD blocks ----
  # Byte-order sort (the cluster's C-locale order), so replicate i gets the same
  # block on every machine regardless of the session locale.
  ld_block_files <- sort(list.files(path = LD_blocks_dir, pattern = "\\.rds$", full.names = TRUE),
                         method = "radix")

  if (length(ld_block_files) == 0) {
    stop("No .rds files found in LD_blocks_dir: ", LD_blocks_dir)
  }

  cat("Found", length(ld_block_files), "LD block files (will recycle if replicates exceed this).\n")

  cat("\n========================================\n")
  cat("Cross-Validation Prediction (Parallel)\n")
  cat("========================================\n")
  cat("Total replicates in study:", total_replicates, "\n")
  cat("Replicates for this job:", start_replicate, "to", end_replicate, "(", num_replicates, "replicates)\n")
  cat("Data type:", data_type, "\n")
  if (data_type == "sparse") {
    cat("Per-SNP heritability (h2_per_snp):", h2_per_snp, "\n")
    cat("Number of causal variants (K):", K, "\n")
    cat("Total heritability:", h2_per_snp * K, "\n")
  } else if (data_type == "oligo") {
    cat("Total SNP heritability (h2g):", h2g, "\n")
    cat("Number of sparse effects (K):", K, "\n")
    cat("Proportion h2 from sparse:", prop_h2_sparse, "\n")
    cat("Proportion h2 from oligogenic:", prop_h2_oligogenic, "\n")
    cat("Proportion h2 from infinitesimal:", prop_h2_infinitesimal, "\n")
    cat("Number of oligogenic effects:", n_oligogenic, "\n")
    cat("Number of infinitesimal effects (n_inf):", n_inf, "\n")
  }
  if (!is.null(n)) {
    cat("Number of samples (n):", n, "\n")
  }
  cat("Number of CV folds:", n_folds, "\n")
  cat("LD blocks directory:", LD_blocks_dir, "\n")
  cat("Output file:", final_name, ".rds\n")
  cat("========================================\n\n")

  # ---- Replicates ----
  all_results <- vector("list", num_replicates)

  job_start_time <- Sys.time()

  for (i in start_replicate:end_replicate) {

    replicate_start <- Sys.time()
    result_index <- i - start_replicate + 1  # Index in all_results vector

    cat("========================================\n")
    cat("Running replicate", i, "of", total_replicates,
        "(job replicate", result_index, "of", num_replicates, ")\n")
    ld_block_index <- ((i - 1) %% length(ld_block_files)) + 1
    cat("Using LD block file:", basename(ld_block_files[ld_block_index]), "\n")
    cat("========================================\n")

    tryCatch({
      ld_block <- readRDS(ld_block_files[ld_block_index])
      cat("Successfully loaded LD block:", basename(ld_block_files[ld_block_index]), "\n")
    }, error = function(e) {
      cat("\n!!! ERROR loading LD block file:", ld_block_files[ld_block_index], "\n")
      cat("Error message:", e$message, "\n")
      stop("Failed to load LD block file")
    })

    # Subset, impute, filter and cap as the fine-mapping benchmark does.
    tryCatch({
      X_working <- ld_block$genotypes
      original_n <- nrow(X_working)
      original_p <- ncol(X_working)
      cat("Original genotype matrix: n=", original_n, ", p=", original_p, "\n")

      if (!is.null(n)) {
        if (n > original_n) {
          cat("WARNING: Requested n=", n, " but only", original_n, "samples available. Using all.\n")
          n_use <- original_n
        } else {
          cat("Subsetting to first", n, "samples\n")
          X_working <- X_working[1:n, , drop = FALSE]
          n_use <- n
        }
      } else {
        n_use <- original_n
      }

      col_means <- colMeans(X_working, na.rm = TRUE)
      na_idx <- which(is.na(X_working), arr.ind = TRUE)
      if (nrow(na_idx) > 0) {
        cat("Mean-imputing", nrow(na_idx), "missing value(s)\n")
        X_working[na_idx] <- col_means[na_idx[, 2]]
      }

      maf <- colMeans(X_working) / 2
      maf <- pmin(maf, 1 - maf)
      maf_keep <- maf >= 0.01
      n_maf_removed <- sum(!maf_keep)
      if (n_maf_removed > 0) {
        cat("MAF filter: removing", n_maf_removed, "variant(s) with MAF < 1%\n")
      }
      X_working <- X_working[, maf_keep, drop = FALSE]

      col_vars <- apply(X_working, 2, var)
      non_zero_var_cols <- which(col_vars > 0)
      n_zero_var <- ncol(X_working) - length(non_zero_var_cols)
      if (n_zero_var > 0) {
        cat("Filtering out", n_zero_var, "zero-variance variant(s)\n")
        X_working <- X_working[, non_zero_var_cols, drop = FALSE]
      }

      # The cap of 5000 variants applies after the filters.
      if (ncol(X_working) > 5000) {
        cat("Subsetting to first 5000 variants (had", ncol(X_working), ")\n")
        X_working <- X_working[, 1:5000, drop = FALSE]
      } else {
        cat("Using all", ncol(X_working), "variants (<=5000)\n")
      }

      # Scale X for the phenotype simulation; the fits below see X_working.
      X_scaled <- scale(X_working)

      cat("Final genotype matrix dimensions: n=", nrow(X_scaled), ", p=", ncol(X_scaled), "\n")
    }, error = function(e) {
      cat("\n!!! ERROR subsetting/filtering genotype matrix\n")
      cat("Error message:", e$message, "\n")
      stop("Failed to preprocess genotype matrix")
    })

    # The phenotype is simulated on the standardized matrix, so that the
    # heritability of a causal variant does not depend on its MAF.
    cat("Generating data using scaled X...\n")
    seed <- i + 1000
    data <- tryCatch({
      generate_data(X = X_scaled,
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
    }, error = function(e) {
      cat("\n!!! ERROR generating data\n")
      cat("Data type:", data_type, "\n")
      if (data_type == "sparse") {
        cat("h2_per_snp:", h2_per_snp, ", K:", K, "\n")
      } else {
        cat("h2g:", h2g, ", K:", K, "\n")
      }
      if (!is.null(n)) cat("n (samples):", n, "\n")
      cat("Seed:", seed, "\n")
      cat("Error message:", e$message, "\n")
      cat("Traceback:\n")
      print(sys.calls())
      stop("Failed to generate data")
    })
    cat("Data generated successfully. Causal variants:", length(data$causal), "\n")

    # Cross-validation works on the unscaled filtered genotypes: susie()
    # standardizes each training fold internally and coef() maps back to this
    # scale, so the test fold never contributes to centering or scaling.
    data$G <- X_working
    cat("Data dimensions: n=", nrow(data$G), ", p=", ncol(data$G), "\n\n")

    # Balanced random folds, reseeded so that they depend only on the replicate.
    set.seed(seed)
    n_samples <- nrow(data$G)
    fold_ids <- sample(rep(1:n_folds, length.out = n_samples))

    # Drop the columns that are constant over all samples and remap beta and
    # the causal indices to the columns kept. X_working has no such column
    # (they were filtered out above), so this removes nothing.
    cat("  Pre-filtering globally constant columns...\n")
    original_p <- ncol(data$G)
    col_vars <- apply(data$G, 2, var)
    cols_to_keep <- which(col_vars > 0)
    n_removed <- original_p - length(cols_to_keep)

    original_causal <- data$causal
    original_beta <- data$beta

    if (n_removed > 0) {
      cat("  Removing", n_removed, "globally constant column(s)\n")

      causal_removed <- original_causal[!original_causal %in% cols_to_keep]
      if (length(causal_removed) > 0) {
        cat("  WARNING:", length(causal_removed), "causal variant(s) were constant and removed\n")
      }

      data$beta <- original_beta[cols_to_keep]

      causal_kept <- original_causal[original_causal %in% cols_to_keep]
      data$causal <- match(causal_kept, cols_to_keep)

      data$G <- data$G[, cols_to_keep, drop = FALSE]
    }

    data$kept_variant_indices <- cols_to_keep
    data$original_causal <- original_causal
    data$original_beta <- original_beta

    cat("  Final dimensions for CV: n=", nrow(data$G), ", p=", ncol(data$G), "\n")
    cat("  Causal variants in filtered matrix:", length(data$causal), "of", length(original_causal), "\n\n")

    # ---- Cross-validation ----
    fold_results <- list()

    for (fold_k in 1:n_folds) {

      cat("  Fold", fold_k, "of", n_folds, "\n")

      train_idx <- fold_ids != fold_k
      test_idx <- fold_ids == fold_k

      X_train_full <- data$G[train_idx, , drop = FALSE]
      y_train <- data$y[train_idx]
      X_test_full <- data$G[test_idx, , drop = FALSE]
      y_test <- data$y[test_idx]

      # Columns that are constant within the training fold are dropped for
      # fitting and get a zero coefficient at prediction.
      train_col_vars <- apply(X_train_full, 2, var)
      nonzero_cols <- which(train_col_vars > 1e-10)
      zero_cols <- which(train_col_vars <= 1e-10)
      p_full <- ncol(X_train_full)

      if (length(zero_cols) > 0) {
        cat("    Found", length(zero_cols), "constant column(s) - will filter for fitting\n")
      }

      X_train <- X_train_full[, nonzero_cols, drop = FALSE]
      X_test <- X_test_full[, nonzero_cols, drop = FALSE]

      cat("    Train samples:", sum(train_idx), ", Test samples:", sum(test_idx),
          ", p_full:", p_full, ", p_fit:", ncol(X_train), "\n")

      # get_predictions() with zero coefficients for the dropped columns.
      get_predictions_with_padding <- function(fit, X_test_full, nonzero_cols, p_full) {
        beta_with_intercept <- coef(fit)
        intercept <- beta_with_intercept[1]
        beta_filtered <- beta_with_intercept[-1]

        beta_full <- rep(0, p_full)
        beta_full[nonzero_cols] <- beta_filtered

        as.vector(intercept + X_test_full %*% beta_full)
      }

      fold_result <- list()

      # Fit one method and score it on the test fold. An error in the fit or
      # in the prediction gives NA metrics and a recorded message, so a single
      # bad fold cannot abort the job.
      fit_and_score <- function(label, fit_call) {
        tryCatch({
          fit <- fit_call()
          y_pred <- get_predictions_with_padding(fit, X_test_full, nonzero_cols, p_full)
          list(rmse = calculate_rmse(y_test, y_pred), r2 = cor(y_test, y_pred)^2, error = NULL)
        }, error = function(e) {
          cat("    WARNING:", label, "failed:", conditionMessage(e), "\n")
          list(rmse = NA_real_, r2 = NA_real_, error = conditionMessage(e))
        })
      }

      cat("    Running SuSiE...\n")
      fold_result$SuSiE <- fit_and_score("SuSiE", function()
        susie(X = X_train, y = y_train, L = 10, max_iter = 100, verbose = FALSE))

      cat("    Running SuSiE-inf...\n")
      fold_result$SuSiE_inf <- fit_and_score("SuSiE-inf", function()
        susie(X = X_train, y = y_train, L = 10, unmappable_effects = "inf",
              convergence_method = convergence_method, max_iter = 100, verbose = FALSE))

      cat("    Running SuSiE-ash...\n")
      fold_result$SuSiE_ash <- fit_and_score("SuSiE-ash", function()
        susie(X = X_train, y = y_train, L = 10, unmappable_effects = "ash",
              convergence_method = convergence_method, max_iter = 100, verbose = FALSE))

      gc(verbose = FALSE)

      fold_results[[paste0("fold_", fold_k)]] <- fold_result

      cat("    Fold", fold_k, "complete\n\n")
    }

    # Replicate metrics: mean and SD over folds, leaving out failed folds (NA).
    cat("Calculating averaged CV metrics...\n")

    methods <- c("SuSiE", "SuSiE_inf", "SuSiE_ash")
    cv_metrics <- list()

    for (method in methods) {
      rmse_values <- sapply(fold_results, function(f) f[[method]]$rmse)
      r2_values <- sapply(fold_results, function(f) f[[method]]$r2)

      cv_metrics[[method]] <- list(
        mean_rmse = mean(rmse_values, na.rm = TRUE),
        sd_rmse = sd(rmse_values, na.rm = TRUE),
        mean_r2 = mean(r2_values, na.rm = TRUE),
        sd_r2 = sd(r2_values, na.rm = TRUE)
      )

      cat("  ", method, ": RMSE =", round(cv_metrics[[method]]$mean_rmse, 4),
          "(SD:", round(cv_metrics[[method]]$sd_rmse, 4), "), R² =",
          round(cv_metrics[[method]]$mean_r2, 4),
          "(SD:", round(cv_metrics[[method]]$sd_r2, 4), ")\n")
    }

    n_fit_errors <- sum(vapply(fold_results, function(f)
      sum(vapply(f, function(m) !is.null(m$error), logical(1))), integer(1)))
    if (n_fit_errors > 0) {
      cat("  WARNING:", n_fit_errors, "fit(s) failed in this replicate; those folds are NA\n")
    }

    replicate_elapsed <- as.numeric(difftime(Sys.time(), replicate_start, units = "secs"))

    all_results[[result_index]] <- list(
      replicate_id = i,
      ld_block_name = basename(ld_block_files[ld_block_index]),
      seed = seed,
      # Truth in the columns of the matrix used for cross-validation.
      causal_indices = data$causal,
      beta = data$beta,
      # Truth in the columns of the simulated matrix, and the positions of the
      # cross-validation columns in it. The two sets of columns are the same
      # unless the constant-column step above dropped some.
      original_causal = data$original_causal,
      original_beta = data$original_beta,
      kept_variant_indices = data$kept_variant_indices,
      fold_results = fold_results,
      cv_metrics = cv_metrics,
      n_fit_errors = n_fit_errors,
      elapsed_time = replicate_elapsed
    )

    cat("Replicate", i, "completed in", round(replicate_elapsed, 2), "seconds\n\n")

    rm(ld_block, data, fold_results, cv_metrics)
    gc()
  }

  # ---- Save ----
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
      n_folds = n_folds,
      convergence_method = convergence_method,
      susieR_version = as.character(packageVersion("susieR")),
      job_elapsed_time = job_elapsed
    )
  )

  final_file <- paste0(final_name, ".rds")
  output_path <- file.path(output_dir, final_file)

  cat("========================================\n")
  cat("Saving results to:", output_path, "\n")
  cat("Job completed in", round(job_elapsed / 60, 2), "minutes\n")
  cat("========================================\n")

  # Write to a temporary name and rename, so an interrupted job never leaves
  # a truncated file that a resumed run would take for a finished one.
  tmp_path <- paste0(output_path, ".tmp")
  saveRDS(results, tmp_path)
  file.rename(tmp_path, output_path)

  cat("\nCross-validation job completed successfully!\n")

  return(results)
}

# Every setting comes from the command line, parsed inside the function.
results <- prediction_cv_parallel(
  total_replicates = NULL,
  start_replicate = NULL,
  end_replicate = NULL,
  num_replicates = NULL,
  h2_per_snp = NULL,
  h2g = NULL,
  data_type = NULL,
  LD_blocks_dir = NULL,
  K = NULL,
  n = NULL,
  n_inf = NULL
)
