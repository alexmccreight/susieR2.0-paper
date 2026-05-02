#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S6: Extract per-replicate fit runtime by method × scenario
# =============================================================================
# For each scenario (Sparse, Complex, Complex S1, Complex S2), pulls the wall
# time recorded for each of the three fine-mapping fits per replicate. Wall
# times come from rep$fits$<method>$elapsed_time (set in
# simulation_script_parallel.R's run_all_methods()).
#
# Inputs (under final_scripts/500_rep_results/):
#   - sparse/sim_nrep500_sparse_h2persnp0.03_K{1..5}_n1000.rds
#   - complex/sim_nrep500_oligo_..._nOligo5_nInf15_n1000.rds        -> "Complex"
#   - complex_S1/sim_nrep500_oligo_..._nOligo5_nInfall_n1000.rds    -> "Complex S1"
#   - complex_S2/sim_nrep500_oligo_..._nOligo10_nInfall_n1000.rds   -> "Complex S2"
#
# Output:
#   - {s6_dir}/data/s6_runtime_data.rds  (long-format data frame, one row
#                                          per replicate x method)
# =============================================================================

# =============================================================================
# Paths
# =============================================================================

s6_dir          <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/supplemental/S6"
rep_results_dir <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/500_rep_results"

output_dir <- file.path(s6_dir, "data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Each scenario points at one or more raw sim files
scenarios <- list(
  list(
    label = "Sparse",
    files = file.path(rep_results_dir, "sparse",
                      sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000.rds", 1:5))
  ),
  list(
    label = "Complex",
    files = file.path(rep_results_dir, "complex",
      "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInf15_n1000.rds")
  ),
  list(
    label = "Complex S1",
    files = file.path(rep_results_dir, "complex_S1",
      "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInfall_n1000.rds")
  ),
  list(
    label = "Complex S2",
    files = file.path(rep_results_dir, "complex_S2",
      "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.15_pInf0.35_nOligo10_nInfall_n1000.rds")
  )
)

# =============================================================================
# Helper: extract per-method elapsed_time for every replicate in one sim file
# =============================================================================

method_map <- c(susie = "SuSiE",
                susie_ash = "SuSiE-ash",
                susie_inf = "SuSiE-inf")

extract_runtime_rows <- function(rds_path) {
  d <- readRDS(rds_path)
  rows <- list()

  for (rep in d$replicates) {
    if (is.null(rep)) next
    rep_id <- rep$replicate_id
    for (fit_name in names(method_map)) {
      et <- rep$fits[[fit_name]]$elapsed_time
      if (is.null(et)) next
      rows[[length(rows) + 1]] <- data.frame(
        replicate    = rep_id,
        method       = method_map[[fit_name]],
        elapsed_time = as.numeric(et),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

# =============================================================================
# Main loop
# =============================================================================

job_start <- Sys.time()
all_rows <- list()

for (sc in scenarios) {
  cat(sprintf("\n=== Scenario: %s (%d file%s) ===\n",
              sc$label, length(sc$files),
              if (length(sc$files) == 1) "" else "s"))

  for (f in sc$files) {
    if (!file.exists(f)) {
      warning("Missing input: ", f)
      next
    }
    cat(sprintf("  Reading: %s\n", basename(f)))
    rows <- extract_runtime_rows(f)
    if (is.null(rows)) next
    rows$scenario <- sc$label
    all_rows[[length(all_rows) + 1]] <- rows
  }
}

runtime_data <- do.call(rbind, all_rows)
rownames(runtime_data) <- NULL

# =============================================================================
# Set factor levels for downstream plotting
# =============================================================================

method_levels   <- c("SuSiE", "SuSiE-ash", "SuSiE-inf")
scenario_levels <- c("Sparse", "Complex", "Complex S1", "Complex S2")

runtime_data$method   <- factor(runtime_data$method,   levels = method_levels)
runtime_data$scenario <- factor(runtime_data$scenario, levels = scenario_levels)

# =============================================================================
# Save
# =============================================================================

out <- list(
  runtime_data = runtime_data,
  meta = list(
    method_levels   = method_levels,
    scenario_levels = scenario_levels,
    elapsed_seconds = as.numeric(difftime(Sys.time(), job_start, units = "secs"))
  )
)

output_file <- file.path(output_dir, "s6_runtime_data.rds")
saveRDS(out, output_file)

cat(sprintf("\nExtraction complete in %.1f s\n", out$meta$elapsed_seconds))
cat(sprintf("Saved: %s\n", output_file))
cat(sprintf("Total rows: %d (= reps x 3 methods)\n", nrow(runtime_data)))

cat("\n=== Mean elapsed_time (s) by scenario x method ===\n")
agg <- aggregate(elapsed_time ~ scenario + method, data = runtime_data,
                 FUN = function(v) c(mean = round(mean(v, na.rm = TRUE), 1),
                                     sd   = round(sd(v,   na.rm = TRUE), 1)))
print(agg)
