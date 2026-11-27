#!/usr/bin/env Rscript

# Supplementary Figure S4: cross-validated prediction R2 of SuSiE, SuSiE-inf
# and SuSiE-ash in the sparse and complex simulations.
#
# A replicate's R2 for a method is its mean over the five cross-validation
# folds (cv_metrics$<method>$mean_r2), and the best method of a replicate is
# the one with the highest R2. The *_prop tables give the proportion of
# replicates in which each method was best; the sparse table pools K = 1-5.
#
# Scenario names in the output:
#   sparse    K = 1-5 causal variants
#   complex1  oligogenic effects on a polygenic background
#   complex2  oligogenic effects on a moderate infinitesimal background
#   complex3  oligogenic effects on an extensive infinitesimal background
#
# Input:  <benchmark>/final_scripts/prediction_cv_results/
#           sparse/pred_cv_nrep300_sparse_h2persnp0.03_K{1-5}_n1000_nfolds5.rds
#           complex/pred_cv_nrep300_oligo_*_n1000_nfolds5.rds  (three scenarios)
# Output: data/prediction_data.rds  (list: <scenario>_data, <scenario>_prop)

suppressPackageStartupMessages(library(dplyr))

source("R/paths.R")

input_dir   <- file.path(benchmark_root(), "final_scripts", "prediction_cv_results")
output_file <- file.path(supp_dir(4), "data", "prediction_data.rds")

methods <- c("SuSiE", "SuSiE_inf", "SuSiE_ash")
method_labels <- c("SuSiE", "SuSiE-inf", "SuSiE-ash")

complex_files <- c(
  complex1 = "pred_cv_nrep300_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInf15_n1000_nfolds5.rds",
  complex2 = "pred_cv_nrep300_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInfall_n1000_nfolds5.rds",
  complex3 = "pred_cv_nrep300_oligo_h2g0.25_K3_pSparse0.5_pOligo0.15_pInf0.35_nOligo10_nInfall_n1000_nfolds5.rds"
)

# Per-replicate R2 of each method in one results file.
extract_r2 <- function(file_path) {
  d <- readRDS(file_path)

  rows <- lapply(seq_along(d$replicates), function(i) {
    rep <- d$replicates[[i]]
    # sapply() names the values by method; the saved tables keep these as row names.
    r2_vals <- sapply(methods, function(m) rep$cv_metrics[[m]]$mean_r2)
    data.frame(
      replicate = i,
      Method = method_labels,
      R2 = r2_vals
    )
  })

  bind_rows(rows)
}

# Proportion of replicates (within each level of group_vars, if given) in which
# each method has the highest R2. A tie would go to the method listed first.
compute_prop <- function(data, group_vars = NULL) {
  data %>%
    group_by(across(all_of(c(group_vars, "replicate")))) %>%
    slice_max(R2, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    count(Method) %>%
    mutate(Proportion = n / sum(n))
}

# ---- Sparse: K = 1-5 ----

sparse_data <- bind_rows(lapply(1:5, function(K) {
  df <- extract_r2(file.path(
    input_dir, "sparse",
    sprintf("pred_cv_nrep300_sparse_h2persnp0.03_K%d_n1000_nfolds5.rds", K)
  ))
  df$K <- K
  df$Total_PVE <- 0.03 * K
  df
}))

out <- list(
  sparse_data = sparse_data,
  sparse_prop = compute_prop(sparse_data, group_vars = "K")
)

# ---- Complex ----

for (scenario in names(complex_files)) {
  r2 <- extract_r2(file.path(input_dir, "complex", complex_files[[scenario]]))
  out[[paste0(scenario, "_data")]] <- r2
  out[[paste0(scenario, "_prop")]] <- compute_prop(r2)
}

saveRDS(out, output_file)
cat("Saved:", output_file, "\n")
