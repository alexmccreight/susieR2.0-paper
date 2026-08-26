#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S9: Extract per-method runtime for SuSiE, SuSiE-ash, SuSiE-inf
# =============================================================================
# Each replicate stores per-method wall-clock fit times under
#   rep$fits$<method>$elapsed_time   (seconds; from Sys.time() differences)
# for the three methods susie / susie_ash / susie_inf. This script pulls those
# scalars across all 500 replicates and aggregates to mean +/- SE per method.
#
# Input (read-only):
#   500_rep_results/runtime/sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_
#     pOligo0.15_pInf0.35_nOligo10_nInfall_n1000.rds
# Outputs:
#   data/s9_runtime_data.rds            (aggregated summary, used by plot_S9.R)
#   data/s9_runtime_per_replicate.rds   (raw long frame, for transparency)
# =============================================================================

source("R/paths.R")

s9_dir     <- supp_dir(9)
input_file <- file.path(
  file.path(benchmark_root(), "final_scripts", "500_rep_results", "runtime"),
  "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.15_pInf0.35_nOligo10_nInfall_n1000.rds")
output_dir <- file.path(s9_dir, "data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# fits key -> display name (x-axis / legend order)
method_names <- c("susie" = "SuSiE", "susie_ash" = "SuSiE-ash", "susie_inf" = "SuSiE-inf")
method_levels <- c("SuSiE", "SuSiE-ash", "SuSiE-inf")

if (!file.exists(input_file)) stop("Missing runtime file: ", input_file)

cat("Loading runtime simulation (this is a large file)...\n")
sim  <- readRDS(input_file)
reps <- sim$replicates
cat(sprintf("  %d replicates found\n", length(reps)))

# -----------------------------------------------------------------------------
# Pull per-method elapsed_time (seconds) from every valid replicate
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# Aggregate: mean +/- SE (and sd / median) per method across replicates
# -----------------------------------------------------------------------------
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

saveRDS(agg,     file.path(output_dir, "s9_runtime_data.rds"))
saveRDS(per_rep, file.path(output_dir, "s9_runtime_per_replicate.rds"))

cat(sprintf("\nSaved: %s\n", file.path(output_dir, "s9_runtime_data.rds")))
cat(sprintf("Saved: %s\n", file.path(output_dir, "s9_runtime_per_replicate.rds")))
cat(sprintf("Per-replicate rows: %d  (expect ~%d = 3 methods x valid reps)\n",
            nrow(per_rep), 3 * length(unique(per_rep$replicate))))
cat("\n=== Runtime summary (seconds) ===\n")
print(agg, row.names = FALSE)
