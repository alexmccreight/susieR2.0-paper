#!/usr/bin/env Rscript

# Supplementary Figure S4: command generator for the cross-validated prediction
# benchmark. Walks the parameter grids and writes two command files:
#   commands_cv_to_submit.txt  one worker command (prediction_cv_script_parallel.R)
#                              per parameter set and replicate chunk
#   commands_cv_combine.txt    one combine command (combine_parallel_cv_results.R)
#                              per parameter set; run after all worker jobs finish
#
# Scenarios (all emitted by default) and their panels in S4:
#   sparse      K = 1-5 causal variants, h2 = 0.03 per variant          A, E
#   complex     oligogenic, on a polygenic background                   B, F
#   complex_s1  oligogenic, on a moderate infinitesimal background      C, G
#   complex_s2  oligogenic, on an extensive infinitesimal background    D, H
#
# Usage: Rscript command_generator_parallel_cv.R [key=value ...]
#   scripts_dir         directory holding the worker and the combiner
#   ld_blocks_dir       directory of LD block rds files
#   output_dir          where the worker jobs write their rds files
#   scenarios           comma-separated subset of the four scenarios
#   total_replicates    replicates per parameter set (300; the fine-mapping
#                       benchmark uses 500)
#   replicates_per_job  replicates per worker job (2); every replicate is seeded
#                       by its own index, so the chunking does not change results
#   convergence_method  for SuSiE-inf and SuSiE-ash ("pip"; the package default
#                       is "elbo")
#   skip_existing       TRUE makes a worker job exit at once if its output exists
#   commands_dir        where the two command files are written (".")
# The path defaults are the cluster layout (see `opts`).
#
# Input:  none
# Output: <commands_dir>/commands_cv_to_submit.txt
#         <commands_dir>/commands_cv_combine.txt

args <- commandArgs(trailingOnly = TRUE)
opts <- list(
  scripts_dir        = "/home/apm2217/data/benchmark",
  ld_blocks_dir      = "oligo_LD_blocks",
  output_dir         = "/home/apm2217/output",
  scenarios          = "sparse,complex,complex_s1,complex_s2",
  total_replicates   = 300L,
  replicates_per_job = 2L,
  convergence_method = "pip",
  skip_existing      = FALSE,
  commands_dir       = "."
)
for (arg in args) {
  kv <- regmatches(arg, regexpr("=", arg), invert = TRUE)[[1]]
  if (length(kv) != 2 || !kv[1] %in% names(opts)) stop("Unknown argument: ", arg)
  value <- gsub("^'|'$", "", gsub('^"|"$', "", kv[2]))
  opts[[kv[1]]] <- switch(class(opts[[kv[1]]]),
                          integer = as.integer(value),
                          logical = as.logical(value),
                          value)
}

# ---- Parameter grids ----

# n_inf is the number of variants with an infinitesimal effect; "all" means
# every variant that is neither sparse nor oligogenic. It is a string so that
# both forms fit in one column.

grids <- list(
  sparse = expand.grid(
    total_replicates = opts$total_replicates,
    h2_per_snp       = c(0.03),
    data_type        = c("sparse"),
    K                = c(1, 2, 3, 4, 5),
    n                = c(1000),
    n_folds          = c(5),
    stringsAsFactors = FALSE
  ),
  complex = expand.grid(
    total_replicates      = opts$total_replicates,
    h2g                   = c(0.25),
    data_type             = c("oligo"),
    K                     = c(3),
    n                     = c(1000),
    prop_h2_sparse        = c(0.5),
    prop_h2_oligogenic    = c(0.35),
    prop_h2_infinitesimal = c(0.15),
    n_oligogenic          = c(5),
    n_inf                 = c("15"),
    n_folds               = c(5),
    stringsAsFactors      = FALSE
  ),
  complex_s1 = expand.grid(
    total_replicates      = opts$total_replicates,
    h2g                   = c(0.25),
    data_type             = c("oligo"),
    K                     = c(3),
    n                     = c(1000),
    prop_h2_sparse        = c(0.5),
    prop_h2_oligogenic    = c(0.35),
    prop_h2_infinitesimal = c(0.15),
    n_oligogenic          = c(5),
    n_inf                 = c("all"),
    n_folds               = c(5),
    stringsAsFactors      = FALSE
  ),
  complex_s2 = expand.grid(
    total_replicates      = opts$total_replicates,
    h2g                   = c(0.25),
    data_type             = c("oligo"),
    K                     = c(3),
    n                     = c(1000),
    prop_h2_sparse        = c(0.5),
    prop_h2_oligogenic    = c(0.15),
    prop_h2_infinitesimal = c(0.35),
    n_oligogenic          = c(10),
    n_inf                 = c("all"),
    n_folds               = c(5),
    stringsAsFactors      = FALSE
  )
)

scenarios <- trimws(strsplit(opts$scenarios, ",")[[1]])
unknown <- setdiff(scenarios, names(grids))
if (length(unknown) > 0) stop("Unknown scenario(s): ", paste(unknown, collapse = ", "))
# Stack the selected grids; a column that a scenario does not use is NA.
parameter_grid <- do.call(rbind, lapply(scenarios, function(s) {
  g <- grids[[s]]
  for (col in setdiff(unique(unlist(lapply(grids, names))), names(g))) g[[col]] <- NA
  g$scenario <- s
  g
}))

worker   <- file.path(opts$scripts_dir, "prediction_cv_script_parallel.R")
combiner <- file.path(opts$scripts_dir, "combine_parallel_cv_results.R")
dir.create(opts$commands_dir, recursive = TRUE, showWarnings = FALSE)

# Scenario-specific worker arguments and output base name for one grid row.
# The base name must match the one the worker builds: the combine command uses
# it to find the per-job files.
scenario_args <- function(params) {
  if (params[["data_type"]] == "sparse") {
    list(args = paste0(" h2_per_snp=", params[["h2_per_snp"]], " K=", params[["K"]]),
         pattern = paste0("pred_cv_nrep", params[["total_replicates"]], "_sparse",
                          "_h2persnp", params[["h2_per_snp"]], "_K", params[["K"]]))
  } else {
    list(args = paste0(" h2g=", params[["h2g"]], " K=", params[["K"]],
                       " prop_h2_sparse=", params[["prop_h2_sparse"]],
                       " prop_h2_oligogenic=", params[["prop_h2_oligogenic"]],
                       " prop_h2_infinitesimal=", params[["prop_h2_infinitesimal"]],
                       " n_oligogenic=", params[["n_oligogenic"]],
                       " n_inf=", params[["n_inf"]]),
         pattern = paste0("pred_cv_nrep", params[["total_replicates"]], "_oligo",
                          "_h2g", params[["h2g"]], "_K", params[["K"]],
                          "_pSparse", params[["prop_h2_sparse"]],
                          "_pOligo", params[["prop_h2_oligogenic"]],
                          "_pInf", params[["prop_h2_infinitesimal"]],
                          "_nOligo", params[["n_oligogenic"]],
                          "_nInf", params[["n_inf"]]))
  }
}

# ---- Worker commands ----

commands_file <- file.path(opts$commands_dir, "commands_cv_to_submit.txt")
file_conn <- file(commands_file, open = "w")
total_commands <- 0

for (i in seq_len(nrow(parameter_grid))) {
  params <- parameter_grid[i, ]
  total_replicates <- params[["total_replicates"]]
  num_jobs <- ceiling(total_replicates / opts$replicates_per_job)
  sa <- scenario_args(params)

  cat("\n========================================\n")
  cat("Parameter set", i, "of", nrow(parameter_grid), "(", params[["scenario"]], ")\n")
  cat("Data type:", params[["data_type"]], " K:", params[["K"]], " CV folds:", params[["n_folds"]], "\n")
  cat("Total replicates:", total_replicates, " Replicates per job:", opts$replicates_per_job,
      " Jobs:", num_jobs, "\n")
  cat("========================================\n")

  for (job in seq_len(num_jobs)) {
    start_rep <- (job - 1) * opts$replicates_per_job + 1
    end_rep   <- min(job * opts$replicates_per_job, total_replicates)
    command <- paste0("Rscript ", worker,
                      " total_replicates=", total_replicates,
                      " start_replicate=", start_rep,
                      " end_replicate=", end_rep,
                      " data_type='", params[["data_type"]], "'",
                      sa$args,
                      " n=", params[["n"]],
                      " n_folds=", params[["n_folds"]],
                      " convergence_method='", opts$convergence_method, "'",
                      " LD_blocks_dir='", opts$ld_blocks_dir, "'",
                      " output_dir='", opts$output_dir, "'",
                      if (opts$skip_existing) " skip_existing=TRUE" else "")
    writeLines(command, file_conn)
    total_commands <- total_commands + 1
  }
}
close(file_conn)

cat("\n========================================\n")
cat("Commands file", commands_file, "created successfully.\n")
cat("Total commands generated:", total_commands, "\n")
cat("Parameter sets:", nrow(parameter_grid), "\n")
cat("========================================\n\n")

# ---- Combine commands ----

combine_commands_file <- file.path(opts$commands_dir, "commands_cv_combine.txt")
combine_conn <- file(combine_commands_file, open = "w")
for (i in seq_len(nrow(parameter_grid))) {
  params <- parameter_grid[i, ]
  pattern <- paste0(scenario_args(params)$pattern, "_n", params[["n"]], "_nfolds", params[["n_folds"]])
  writeLines(paste0("Rscript ", combiner,
                    " pattern='", pattern, "'",
                    " input_dir='", opts$output_dir, "'"), combine_conn)
}
close(combine_conn)

cat("Combine commands file", combine_commands_file, "created successfully.\n")
cat("Total combine commands:", nrow(parameter_grid), "\n")
cat("\n========================================\n")
cat("USAGE:\n")
cat("1. Submit all CV prediction jobs: submit commands from 'commands_cv_to_submit.txt'\n")
cat("2. After all jobs complete, combine results: run commands from 'commands_cv_combine.txt'\n")
cat("========================================\n\n")
