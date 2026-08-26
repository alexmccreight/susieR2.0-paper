#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S5: Extract CS-level data for size & purity figure
#                  (WITH LD extension = 0.99)
# =============================================================================
# Reads the *_ldext.rds files for 4 scenarios and extracts one row per credible
# set per method per replicate, using the LD-extension construction
# (version = "xcorr_ext"): cs_size = length of the extended CS, and
# min_abs_corr = the exact all-pairs minimum |r| within that CS (the purity the
# ldext pipeline stored).
#
# NOTE: vs the previous (un-extended) S5, extension grows credible sets, so the
# CS-size bins shift toward larger sizes and median purity may drop; the total
# number of CS is essentially unchanged (extension does not add/remove CS).
#
# Inputs (read-only, under final_scripts/500_rep_results/):
#   - sparse/sim_nrep500_sparse_h2persnp0.03_K{1..5}_n1000_ldext.rds
#   - complex/sim_nrep500_oligo_..._nOligo5_nInf15_n1000_ldext.rds      -> "Complex"
#   - complex_S1/sim_nrep500_oligo_..._nOligo5_nInfall_n1000_ldext.rds  -> "Complex S1"
#   - complex_S2/sim_nrep500_oligo_..._nOligo10_nInfall_n1000_ldext.rds -> "Complex S2"
#
# Output:
#   - {s5_dir}/data/s5_cs_data.rds  (long-format data frame, one row per CS)
# =============================================================================

# =============================================================================
# Paths
# =============================================================================

source("R/paths.R")

s5_dir          <- supp_dir(5)
rep_results_dir <- file.path(benchmark_root(), "final_scripts", "500_rep_results")

output_dir <- file.path(s5_dir, "data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

VERSION <- "xcorr_ext"   # LD-extension (0.99) credible sets

# Each scenario points at one or more *_ldext.rds files
scenarios <- list(
  list(
    label = "Sparse",
    files = file.path(rep_results_dir, "sparse",
                      sprintf("sim_nrep500_sparse_h2persnp0.03_K%d_n1000_ldext.rds", 1:5))
  ),
  list(
    label = "Complex",
    files = file.path(rep_results_dir, "complex",
      "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInf15_n1000_ldext.rds")
  ),
  list(
    label = "Complex S1",
    files = file.path(rep_results_dir, "complex_S1",
      "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.35_pInf0.15_nOligo5_nInfall_n1000_ldext.rds")
  ),
  list(
    label = "Complex S2",
    files = file.path(rep_results_dir, "complex_S2",
      "sim_nrep500_oligo_h2g0.25_K3_pSparse0.5_pOligo0.15_pInf0.35_nOligo10_nInfall_n1000_ldext.rds")
  )
)

# =============================================================================
# Helper: extract one row per CS from a single *_ldext.rds file
# =============================================================================

method_map <- c(susie = "SuSiE",
                susie_ash = "SuSiE-ash",
                susie_inf = "SuSiE-inf")

extract_cs_rows <- function(rds_path) {
  d <- readRDS(rds_path)
  rows <- list()

  for (rep in d$replicates) {
    if (is.null(rep) || !identical(rep$status, "ok")) next
    rep_id <- rep$replicate_id
    for (fit_name in names(method_map)) {
      cs_rec <- rep$credible_sets[[fit_name]][[VERSION]]
      if (is.null(cs_rec) || is.null(cs_rec$cs) || length(cs_rec$cs) == 0) next
      cs_list <- cs_rec$cs
      purity  <- cs_rec$purity
      for (j in seq_along(cs_list)) {
        rows[[length(rows) + 1]] <- data.frame(
          replicate    = rep_id,
          method       = method_map[[fit_name]],
          cs_size      = length(cs_list[[j]]),
          min_abs_corr = purity[j],
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
    ld_extend       = 0.99,
    version         = VERSION,
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

cat("\n=== Median purity by scenario x method ===\n")
print(round(tapply(all_cs$min_abs_corr,
                   list(all_cs$scenario, all_cs$method),
                   median, na.rm = TRUE), 3))
