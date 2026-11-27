#!/usr/bin/env Rscript

# Supplementary Figure S7, SuSiE+BLiP re-run: synthetic data for the smoke
# test, a small stand-in for the raw sparse benchmark.
#
# LD blocks are Gaussian AR(1) matrices with rho = 0.99, so that 95% credible
# sets span 20-40 variants and several exceed 25. Phenotypes follow the sparse
# regime of simulation_script_parallel.R: K causal variants with pairwise
# |r| <= 0.1, unit effects, total h2 = K x h2_per_snp, replicate seed i + 1000,
# LD blocks recycled in order. Each replicate stores replicate_id,
# causal_indices, beta, ld_block_name and fits$susie$fit, as the
# sim_nrep500_sparse_h2persnp0.03_K<k>_n1000.rds files do.
#
# Output: <output_dir>/ld_blocks/smoke_block_<j>.rds
#         <output_dir>/sparse/sim_nrep<n_reps>_sparse_h2persnp<h2>_K<K>_n<n>.rds
#
# Usage: Rscript simulate_smoke_data.R output_dir=<dir> [n_reps=20] [n=800] \
#          [p=1500] [rho=0.99] [K=2] [h2_per_snp=0.03] [n_blocks=4] [seed=1]

suppressMessages(library(susieR))

.script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else getwd()
})
source(file.path(.script_dir, "..", "blip_common.R"))

# Column j has correlation rho^|j - k| with column k. Dichotomised haplotypes
# with random allele frequencies would attenuate that to ~0.7 between
# neighbours; the Gaussian block keeps the intended LD. Values are centred at
# 1 so the MAF filter in reconstruct_genotype_matrix() keeps every column.
simulate_block <- function(n, p, rho, seed) {
  set.seed(seed)
  z <- matrix(rnorm(n * p), n, p)
  for (j in 2:p) z[, j] <- rho * z[, j - 1] + sqrt(1 - rho^2) * z[, j]
  G <- 1 + 0.5 * z
  G[sample(length(G), round(0.001 * length(G)))] <- NA   # exercises the mean imputation
  G
}

# The rule of select_causal_variants(independent = TRUE) in
# simulation_script_parallel.R.
select_causal <- function(X, K, ld_threshold = 0.1) {
  if (K == 1) return(sample(ncol(X), 1))
  repeat {
    vars <- sample(ncol(X), K)
    r <- cor(X[, vars])
    if (all(abs(r[upper.tri(r)]) <= ld_threshold)) return(vars)
  }
}

main <- function() {
  a <- parse_key_value_args(list(
    output_dir = NULL, n_reps = 20L, n = 800L, p = 1500L, rho = 0.99,
    K = 2L, h2_per_snp = 0.03, n_blocks = 4L, seed = 1L
  ))
  if (is.null(a$output_dir)) stop("output_dir must be specified", call. = FALSE)
  ld_dir  <- file.path(a$output_dir, "ld_blocks")
  sim_dir <- file.path(a$output_dir, "sparse")
  dir.create(ld_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)

  block_files <- character(a$n_blocks)
  for (j in seq_len(a$n_blocks)) {
    path <- file.path(ld_dir, sprintf("smoke_block_%d.rds", j))
    if (!file.exists(path)) {
      saveRDS(list(genotypes = simulate_block(a$n, a$p, a$rho, a$seed * 100 + j),
                   name = basename(path)), path)
    }
    block_files[j] <- path
  }

  reps  <- vector("list", a$n_reps)
  n_big <- 0L
  for (i in seq_len(a$n_reps)) {
    ld_path <- block_files[((i - 1) %% a$n_blocks) + 1]
    X <- reconstruct_genotype_matrix(ld_path, n = a$n)

    set.seed(i + 1000)
    causal <- select_causal(X, a$K)
    beta <- numeric(ncol(X))
    beta[causal] <- 1
    h2_total <- a$h2_per_snp * a$K
    mu <- drop(X %*% beta)
    y  <- mu + rnorm(nrow(X), 0, sqrt(var(mu) * (1 - h2_total) / h2_total))

    t0  <- Sys.time()
    fit <- susie(X, y)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    cs <- susie_get_cs(fit, Xcorr = cor(X), coverage = 0.95, min_abs_corr = 0.5)
    max_cs <- if (is.null(cs$cs)) 0L else max(lengths(cs$cs))
    if (max_cs > 25) n_big <- n_big + 1L
    cat(sprintf("rep %2d: %s, %d CS(s), largest %d, susie %.1fs\n",
                i, basename(ld_path), length(cs$cs), max_cs, elapsed))

    reps[[i]] <- list(
      replicate_id   = i,
      causal_indices = causal,
      beta           = beta,
      ld_block_name  = basename(ld_path),
      fits           = list(susie = list(fit = fit, elapsed_time = elapsed)),
      metrics        = NULL
    )
  }

  sim <- list(
    replicates = reps,
    parameters = list(num_replicates = a$n_reps, h2_per_snp = a$h2_per_snp,
                      data_type = "sparse", K = a$K, n = a$n,
                      LD_blocks_dir = "ld_blocks", rho = a$rho, smoke_test = TRUE)
  )
  out <- file.path(sim_dir, sprintf("sim_nrep%d_sparse_h2persnp%s_K%d_n%d.rds",
                                    a$n_reps, format(a$h2_per_snp), a$K, a$n))
  saveRDS(sim, out)
  cat(sprintf("\nWrote %s\n%d of %d replicates have a SuSiE(0.5) 95%% CS above 25 variants.\n",
              out, n_big, a$n_reps))
}

if (sys.nframe() == 0L) main()
