# =============================================================================
# Command Generator — Prediction CV Pipeline
# =============================================================================
# Walks a parameter grid for one scenario and writes two HPC command files:
#
#   commands_cv_to_submit.txt   one Rscript per (parameter set x rep batch);
#                               submit these jobs in parallel.
#   commands_cv_combine.txt     one Rscript per parameter set; run AFTER all
#                               prediction jobs from the submit file finish.
#
# Each per-rep job runs 5-fold cross-validated prediction R^2 for SuSiE,
# SuSiE-ash, and SuSiE-inf on a single simulated (X, y).
#
# Four scenarios are pre-defined as expand.grid blocks below. To switch
# scenarios, comment out the currently active block and uncomment another:
#
#   - Sparse                 (h2_per_snp x K, K = 1..5)
#   - Complex                (oligogenic + polygenic background)
#   - Complex S1             (oligogenic + moderate infinitesimal background)
#   - Complex S2             (oligogenic + extensive infinitesimal background)
#
# Outputs (`pred_cv_*.rds` files) land in /home/apm2217/output/ on the HPC.
# Note: prediction CV uses 300 replicates (fine-mapping uses 500).
# =============================================================================


# =============================================================================
# Parameter grids (one active at a time; comment / uncomment to switch)
# =============================================================================

# --- Complex: Polygenic Background      (Figure S3, panels B & F) ---------- #
# parameter_grid <- expand.grid(
#   total_replicates      = c(300),
#   h2g                   = c(0.25),
#   data_type             = c('oligo'),
#   K                     = c(3),
#   n                     = c(1000),
#   prop_h2_sparse        = c(0.5),
#   prop_h2_oligogenic    = c(0.35),
#   prop_h2_infinitesimal = c(0.15),
#   n_oligogenic          = c(5),
#   n_inf                 = c(15),     # use 'all' for full infinitesimal background
#   n_folds               = c(5),
#   LD_blocks_dir         = c('oligo_LD_blocks'),
#   stringsAsFactors      = FALSE
# )
# replicates_per_job <- 2

# --- Complex S1: Moderate Infinitesimal (Figure S3, panels C & G) ---------- #
# parameter_grid <- expand.grid(
#   total_replicates      = c(300),
#   h2g                   = c(0.25),
#   data_type             = c('oligo'),
#   K                     = c(3),
#   n                     = c(1000),
#   prop_h2_sparse        = c(0.5),
#   prop_h2_oligogenic    = c(0.35),
#   prop_h2_infinitesimal = c(0.15),
#   n_oligogenic          = c(5),
#   n_inf                 = c('all'),
#   n_folds               = c(5),
#   LD_blocks_dir         = c('oligo_LD_blocks'),
#   stringsAsFactors      = FALSE
# )
# replicates_per_job <- 2

# --- Complex S2: Extensive Infinitesimal (Figure S3, panels D & H) --------- #
# parameter_grid <- expand.grid(
#   total_replicates      = c(300),
#   h2g                   = c(0.25),
#   data_type             = c('oligo'),
#   K                     = c(3),
#   n                     = c(1000),
#   prop_h2_sparse        = c(0.5),
#   prop_h2_oligogenic    = c(0.15),
#   prop_h2_infinitesimal = c(0.35),
#   n_oligogenic          = c(10),
#   n_inf                 = c('all'),
#   n_folds               = c(5),
#   LD_blocks_dir         = c('oligo_LD_blocks'),
#   stringsAsFactors      = FALSE
# )
# replicates_per_job <- 2

# --- Sparse  [ACTIVE]                   (Figure S3, panels A & E) ---------- #
parameter_grid <- expand.grid(
  total_replicates = c(300),
  h2_per_snp       = c(0.03),
  data_type        = c('sparse'),
  K                = c(1, 2, 3, 4, 5),
  n                = c(1000),
  n_folds          = c(5),
  LD_blocks_dir    = c('oligo_LD_blocks'),
  stringsAsFactors = FALSE
)
replicates_per_job <- 10


# =============================================================================
# Generate per-job CV commands -> commands_cv_to_submit.txt
# =============================================================================

commands_file <- "commands_cv_to_submit.txt"
file_conn <- file(commands_file, open = "w")

total_commands <- 0

for (i in 1:nrow(parameter_grid)) {
  params <- parameter_grid[i, ]

  # Extract parameter values
  total_replicates <- params[["total_replicates"]]
  data_type        <- params[["data_type"]]
  K                <- params[["K"]]
  n                <- params[["n"]]
  n_folds          <- params[["n_folds"]]
  LD_blocks_dir    <- params[["LD_blocks_dir"]]

  # Calculate how many jobs needed for this parameter set
  num_jobs <- ceiling(total_replicates / replicates_per_job)

  cat("\n========================================\n")
  cat("Parameter set", i, "of", nrow(parameter_grid), "\n")
  cat("Data type:", data_type, "\n")
  if (data_type == 'sparse') {
    cat("h2_per_snp:", params[["h2_per_snp"]], "\n")
  } else {
    cat("h2g:", params[["h2g"]], "\n")
  }
  cat("K:", K, "\n")
  cat("CV folds:", n_folds, "\n")
  cat("Total replicates:", total_replicates, "\n")
  cat("Replicates per job:", replicates_per_job, "\n")
  cat("Number of parallel jobs:", num_jobs, "\n")
  cat("========================================\n")

  # Generate commands for parallel jobs
  for (job in 1:num_jobs) {
    # Calculate replicate range for this job
    start_rep <- (job - 1) * replicates_per_job + 1
    end_rep <- min(job * replicates_per_job, total_replicates)
    num_reps_this_job <- end_rep - start_rep + 1

    # Create the base command string
    command <- paste0("Rscript /home/apm2217/data/benchmark/prediction_cv_script_parallel.R",
                      " total_replicates=", total_replicates,
                      " start_replicate=", start_rep,
                      " end_replicate=", end_rep,
                      " data_type='", data_type, "'")

    # Add heritability parameter based on data_type
    if (data_type == 'sparse') {
      h2_per_snp <- params[["h2_per_snp"]]
      command <- paste0(command, " h2_per_snp=", h2_per_snp)
      command <- paste0(command, " K=", K)
    } else if (data_type == 'oligo') {
      h2g <- params[["h2g"]]
      command <- paste0(command, " h2g=", h2g)
      command <- paste0(command, " K=", K)

      # Add oligogenic-specific parameters
      prop_h2_sparse <- params[["prop_h2_sparse"]]
      prop_h2_oligogenic <- params[["prop_h2_oligogenic"]]
      prop_h2_infinitesimal <- params[["prop_h2_infinitesimal"]]
      n_oligogenic <- params[["n_oligogenic"]]
      n_inf <- params[["n_inf"]]

      command <- paste0(command, " prop_h2_sparse=", prop_h2_sparse)
      command <- paste0(command, " prop_h2_oligogenic=", prop_h2_oligogenic)
      command <- paste0(command, " prop_h2_infinitesimal=", prop_h2_infinitesimal)
      command <- paste0(command, " n_oligogenic=", n_oligogenic)
      command <- paste0(command, " n_inf=", n_inf)
    }

    # Add n parameter if specified (not NA)
    if (!is.na(n)) {
      command <- paste0(command, " n=", n)
    }

    # Add CV-specific parameters
    command <- paste0(command, " n_folds=", n_folds)

    # Add LD_blocks_dir
    command <- paste0(command, " LD_blocks_dir='", LD_blocks_dir, "'")

    writeLines(command, file_conn)
    total_commands <- total_commands + 1

    cat("  Job", job, "of", num_jobs, ": replicates", start_rep, "to", end_rep,
        "(", num_reps_this_job, "replicates)\n")
  }
}

close(file_conn)

cat("\n========================================\n")
cat("Commands file 'commands_cv_to_submit.txt' created successfully.\n")
cat("Total commands generated:", total_commands, "\n")
cat("Parameter sets:", nrow(parameter_grid), "\n")
cat("Average jobs per parameter set:", round(total_commands / nrow(parameter_grid), 1), "\n")
cat("========================================\n\n")


# =============================================================================
# Generate per-scenario combine commands -> commands_cv_combine.txt
# =============================================================================

cat("Generating combine commands...\n")

combine_commands_file <- "commands_cv_combine.txt"
combine_conn <- file(combine_commands_file, open = "w")

for (i in 1:nrow(parameter_grid)) {
  params <- parameter_grid[i, ]

  # Extract parameter values
  total_replicates <- params[["total_replicates"]]
  data_type        <- params[["data_type"]]
  K                <- params[["K"]]
  n                <- params[["n"]]
  n_folds          <- params[["n_folds"]]

  # Build the pattern to match result files
  pattern <- paste0("pred_cv_nrep", total_replicates, "_", data_type)

  if (data_type == 'sparse') {
    h2_per_snp <- params[["h2_per_snp"]]
    pattern <- paste0(pattern, "_h2persnp", h2_per_snp, "_K", K)
  } else if (data_type == 'oligo') {
    h2g <- params[["h2g"]]
    prop_h2_sparse <- params[["prop_h2_sparse"]]
    prop_h2_oligogenic <- params[["prop_h2_oligogenic"]]
    prop_h2_infinitesimal <- params[["prop_h2_infinitesimal"]]
    n_oligogenic <- params[["n_oligogenic"]]
    n_inf <- params[["n_inf"]]

    pattern <- paste0(pattern, "_h2g", h2g, "_K", K,
                     "_pSparse", prop_h2_sparse,
                     "_pOligo", prop_h2_oligogenic,
                     "_pInf", prop_h2_infinitesimal,
                     "_nOligo", n_oligogenic,
                     "_nInf", n_inf)
  }

  # Add sample size if specified
  if (!is.na(n)) {
    pattern <- paste0(pattern, "_n", n)
  }

  # Add number of folds
  pattern <- paste0(pattern, "_nfolds", n_folds)

  # Create combine command
  combine_cmd <- paste0("Rscript /home/apm2217/data/benchmark/combine_parallel_cv_results.R",
                       " pattern='", pattern, "'",
                       " input_dir='/home/apm2217/output'")

  writeLines(combine_cmd, combine_conn)
}

close(combine_conn)

cat("Combine commands file 'commands_cv_combine.txt' created successfully.\n")
cat("Total combine commands:", nrow(parameter_grid), "\n")
cat("\n========================================\n")
cat("USAGE:\n")
cat("1. Submit all CV prediction jobs: submit commands from 'commands_cv_to_submit.txt'\n")
cat("2. After all jobs complete, combine results: run commands from 'commands_cv_combine.txt'\n")
cat("========================================\n\n")
