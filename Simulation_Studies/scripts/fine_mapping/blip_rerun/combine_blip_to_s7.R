#!/usr/bin/env Rscript

# Supplementary Figure S7, SuSiE+BLiP re-run, step 4: write one BLiP arm in the
# format Supplementary_Figures/S7/extract_S7_data.R reads.
#
# Gathers the per-replicate blip_arms.json files of each K into one rds:
#   replicates[[i]]$replicate_id
#   replicates[[i]]$causal_indices
#   replicates[[i]]$credible_sets$susie_blip_q0.05$cs   list of 1-based index vectors
#   replicates[[i]]$credible_sets$susie_blip_q0.1$cs
#   parameters                                          sim parameters, BLiP settings and versions
# Slot i of `replicates` is replicate i (NULL when it was not processed).
#
# extract_S7_data.R reads arm D_purity_uncapped from blip_purity0.5/ and arm
# D_purity08_uncapped from blip_purity0.8/, both under
# <benchmark>/final_scripts/500_rep_results/.
#
# Input:  <input_dir>/K<k>/sim_parameters.json
#         <input_dir>/K<k>/rep_<id>/{meta.json, blip_arms.json}
# Output: <output_dir>/<sim file base name>_blip.rds, one per K
#
# Usage: Rscript combine_blip_to_s7.R input_dir=<export root> [arm=A_capped] \
#          [output_dir=<input_dir>/s7_<arm>] [k_values=<all found>]

suppressMessages(library(jsonlite))

.script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else getwd()
})
source(file.path(.script_dir, "blip_common.R"))
locate_paper_paths(.script_dir)

read_json_lists <- function(path) {
  fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
}

combine_k <- function(k_dir, arm, output_dir) {
  params_path <- file.path(k_dir, "sim_parameters.json")
  params <- if (file.exists(params_path)) read_json_lists(params_path) else list()

  rep_dirs <- list.files(k_dir, pattern = "^rep_[0-9]+$", full.names = TRUE)
  reps <- list()
  settings <- NULL
  versions <- NULL
  n_missing <- 0L
  for (d in rep_dirs) {
    arms_path <- file.path(d, "blip_arms.json")
    if (!file.exists(arms_path)) { n_missing <- n_missing + 1L; next }
    arms <- read_json_lists(arms_path)
    meta <- read_json_lists(file.path(d, "meta.json"))
    res <- arms$arms[[arm]]
    if (is.null(res)) stop("Arm '", arm, "' not found in ", arms_path, call. = FALSE)

    credible_sets <- list()
    for (qk in names(res)) {
      credible_sets[[paste0("susie_blip_", qk)]] <- list(cs = lapply(res[[qk]]$groups, as.integer))
    }
    rid <- as.integer(meta$replicate_id)
    reps[[rid]] <- list(
      replicate_id   = rid,
      causal_indices = as.integer(meta$causal),
      credible_sets  = credible_sets
    )
    if (is.null(settings)) { settings <- arms$settings; versions <- arms$versions }
  }
  n_present <- length(Filter(Negate(is.null), reps))
  if (n_present == 0) {
    cat(sprintf("  %s: no replicate has blip_arms.json, skipped\n", basename(k_dir)))
    return(invisible(NULL))
  }

  n_total <- if (!is.null(params$num_replicates)) as.integer(params$num_replicates) else length(reps)
  if (length(reps) < n_total) reps[n_total] <- list(NULL)

  base <- if (!is.null(params$sim_file)) sub("\\.rds$", "", params$sim_file) else basename(k_dir)
  out <- list(
    replicates = reps,
    parameters = c(
      params[setdiff(names(params), "sim_file")],
      list(blip_arm      = arm,
           blip_q_values = settings$q_values,
           blip_settings = settings,
           blip_versions = versions,
           blip_source   = "run_blip_arms.py",
           blip_n_replicates_present = n_present)
    )
  )
  out_path <- file.path(output_dir, paste0(base, "_blip.rds"))
  saveRDS(out, out_path)
  cat(sprintf("  %s: %d of %d replicates (%d folders without arms output) -> %s\n",
              basename(k_dir), n_present, n_total, n_missing, out_path))
  invisible(out_path)
}

main <- function() {
  a <- parse_key_value_args(list(input_dir = NULL, arm = "A_capped", output_dir = "", k_values = ""))
  if (is.null(a$input_dir)) stop("input_dir must be specified", call. = FALSE)
  input_dir  <- path.expand(a$input_dir)
  output_dir <- if (nzchar(a$output_dir)) path.expand(a$output_dir) else file.path(input_dir, paste0("s7_", a$arm))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  k_dirs <- list.files(input_dir, pattern = "^K[0-9]+$", full.names = TRUE)
  if (nzchar(a$k_values)) {
    keep <- sprintf("K%s", trimws(strsplit(a$k_values, ",")[[1]]))
    k_dirs <- k_dirs[basename(k_dirs) %in% keep]
  }
  if (length(k_dirs) == 0) stop("No K folders under ", input_dir, call. = FALSE)

  cat("\n========================================\n")
  cat("BLiP re-run: combine arm '", a$arm, "' into the S7 schema\n", sep = "")
  cat("========================================\n")
  for (k_dir in k_dirs) combine_k(k_dir, a$arm, output_dir)
  cat("Done.\n")
}

if (sys.nframe() == 0L) main()
