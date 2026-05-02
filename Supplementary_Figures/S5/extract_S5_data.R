#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S5: Extract CS-level data for size & purity figure
# =============================================================================
# Reads raw SuSiE-fit simulation files for 4 scenarios and extracts one row
# per credible set per method per replicate, with cs_size + min_abs_corr.
#
# Inputs (under final_scripts/500_rep_results/):
#   - sparse/sim_nrep500_sparse_h2persnp0.03_K{1..5}_n1000.rds
#   - complex/sim_nrep500_oligo_..._nOligo5_nInf15_n1000.rds        -> "Complex"
#   - complex_S1/sim_nrep500_oligo_..._nOligo5_nInfall_n1000.rds    -> "Complex S1"
#   - complex_S2/sim_nrep500_oligo_..._nOligo10_nInfall_n1000.rds   -> "Complex S2"
#
# Output:
#   - {s5_dir}/data/s5_cs_data.rds  (long-format data frame, one row per CS)
# =============================================================================

# =============================================================================
# Paths
# =============================================================================

s5_dir          <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/supplemental/S5"
rep_results_dir <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/500_rep_results"

output_dir <- file.path(s5_dir, "data")
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
# Helper: extract one row per CS from a single sim RDS file
# =============================================================================

method_map <- c(susie = "SuSiE",
                susie_ash = "SuSiE-ash",
                susie_inf = "SuSiE-inf")

extract_cs_rows <- function(rds_path) {
  d <- readRDS(rds_path)
  rows <- list()

  for (rep in d$replicates) {
    if (is.null(rep)) next
    rep_id <- rep$replicate_id
    for (fit_name in names(method_map)) {
      sets <- rep$fits[[fit_name]]$fit$sets
      if (is.null(sets) || is.null(sets$cs) || length(sets$cs) == 0) next
      for (j in seq_along(sets$cs)) {
        rows[[length(rows) + 1]] <- data.frame(
          replicate    = rep_id,
          method       = method_map[[fit_name]],
          cs_size      = length(sets$cs[[j]]),
          min_abs_corr = sets$purity$min.abs.corr[j],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

# =============================================================================
# Main loop
# =============================================================================

job_start <- Sys.time()
all_cs_list <- list()

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
    rows <- extract_cs_rows(f)
    if (is.null(rows)) next
    rows$scenario <- sc$label
    all_cs_list[[length(all_cs_list) + 1]] <- rows
  }
}

all_cs <- do.call(rbind, all_cs_list)
rownames(all_cs) <- NULL

# =============================================================================
# Set factor levels for downstream plotting
# =============================================================================

method_levels <- c("SuSiE", "SuSiE-ash", "SuSiE-inf")
all_cs$method <- factor(all_cs$method, levels = method_levels)

scenario_levels <- c("Sparse", "Complex", "Complex S1", "Complex S2")
all_cs$scenario <- factor(all_cs$scenario, levels = scenario_levels)

size_bin_levels <- c("1", "2", "3-5", "6-10", ">10")
all_cs$cs_size_bin <- cut(
  all_cs$cs_size,
  breaks = c(0, 1, 2, 5, 10, Inf),
  labels = size_bin_levels,
  right  = TRUE
)

# =============================================================================
# Save
# =============================================================================

out <- list(
  cs_data = all_cs,
  meta = list(
    method_levels   = method_levels,
    scenario_levels = scenario_levels,
    size_bin_levels = size_bin_levels,
    elapsed_seconds = as.numeric(difftime(Sys.time(), job_start, units = "secs"))
  )
)

output_file <- file.path(output_dir, "s5_cs_data.rds")
saveRDS(out, output_file)

cat(sprintf("\nExtraction complete in %.1f s\n", out$meta$elapsed_seconds))
cat(sprintf("Saved: %s\n", output_file))
cat(sprintf("Total CSs: %d  (across %d scenarios x 3 methods)\n",
            nrow(all_cs), length(scenario_levels)))

cat("\n=== CS counts by scenario x method ===\n")
print(table(all_cs$scenario, all_cs$method))
