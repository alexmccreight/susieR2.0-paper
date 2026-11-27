#!/usr/bin/env Rscript

# Supplementary Figure S11: runtime per method. Collects the wall-clock time of
# every SuSiE, SuSiE-ash and SuSiE-inf fit (rep$fits$<method>$elapsed_time, in
# seconds) over the 500 replicates of the runtime simulation, which uses the
# settings of the complex S2 scenario (oligogenic effects on an extensive
# infinitesimal background), and summarises it per method as mean, standard
# deviation, standard error and median.
#
# Input:  <benchmark>/final_scripts/500_rep_results/runtime/
#           sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.15_pInf0.35_nOligo10_nInfall_n1000.rds
# Output: data/s11_runtime_data.rds           (per-method summary, for plot_S11.R)
#         data/s11_runtime_per_replicate.rds  (one row per replicate and method)

source("R/paths.R")

s11_dir     <- supp_dir(11)
input_file <- file.path(
  file.path(benchmark_root(), "final_scripts", "500_rep_results", "runtime"),
  "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.15_pInf0.35_nOligo10_nInfall_n1000.rds")
output_dir <- file.path(s11_dir, "data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Key in rep$fits -> display name, in plotting order.
method_names <- c("susie" = "SuSiE", "susie_ash" = "SuSiE-ash", "susie_inf" = "SuSiE-inf")
method_levels <- c("SuSiE", "SuSiE-ash", "SuSiE-inf")

if (!file.exists(input_file)) stop("Missing runtime file: ", input_file)

cat("Loading runtime simulation (this is a large file)...\n")
sim  <- readRDS(input_file)
reps <- sim$replicates
cat(sprintf("  %d replicates found\n", length(reps)))

# ---- Per-replicate times ----
# Replicates and methods with no fit or no finite time are skipped.

rows <- list()
for (r in seq_along(reps)) {
  rep <- reps[[r]]
  if (is.null(rep) || is.null(rep$fits)) next

  for (mkey in names(method_names)) {
    fit_m <- rep$fits[[mkey]]
    if (is.null(fit_m)) next
    et <- fit_m$elapsed_time
    if (is.null(et) || !is.finite(et)) next
    rows[[length(rows) + 1]] <- data.frame(
      replicate    = r,
      Method       = method_names[[mkey]],
      elapsed_time = as.numeric(et),
      stringsAsFactors = FALSE
    )
  }
}
per_rep <- do.call(rbind, rows)
per_rep$Method <- factor(per_rep$Method, levels = method_levels)

rm(sim, reps); gc(verbose = FALSE)

# ---- Summary per method ----

agg <- do.call(rbind, lapply(split(per_rep, per_rep$Method, drop = TRUE), function(df) {
  n <- sum(is.finite(df$elapsed_time))
  data.frame(
    Method      = df$Method[1],
    n           = n,
    mean_time   = mean(df$elapsed_time, na.rm = TRUE),
    sd_time     = sd(df$elapsed_time, na.rm = TRUE),
    se_time     = sd(df$elapsed_time, na.rm = TRUE) / sqrt(n),
    median_time = median(df$elapsed_time, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
agg$Method <- factor(agg$Method, levels = method_levels)
agg <- agg[order(agg$Method), ]
rownames(agg) <- NULL

saveRDS(agg,     file.path(output_dir, "s11_runtime_data.rds"))
saveRDS(per_rep, file.path(output_dir, "s11_runtime_per_replicate.rds"))

cat(sprintf("\nSaved: %s\n", file.path(output_dir, "s11_runtime_data.rds")))
cat(sprintf("Saved: %s\n", file.path(output_dir, "s11_runtime_per_replicate.rds")))
cat(sprintf("Per-replicate rows: %d  (expect ~%d = 3 methods x valid reps)\n",
            nrow(per_rep), 3 * length(unique(per_rep$replicate))))
cat("\n=== Runtime summary (seconds) ===\n")
print(agg, row.names = FALSE)
