# =============================================================================
# Combine Parallel Simulation Results — Fine-Mapping Pipeline
# =============================================================================
# Final stage of the 3-step fine-mapping pipeline (generator -> worker ->
# combiner). Run AFTER all per-rep batches written by simulation_script_parallel.R
# have completed.
#
# Reads every per-batch rds in `input_dir` whose filename matches
# `<pattern>_rep<X>-<Y>.rds`, validates structural and parameter consistency
# across jobs, merges the per-rep entries by `replicate_id` (warning on any
# gaps), and writes a single combined rds named `<pattern>.rds` in the same
# directory. Optionally deletes the per-batch source files afterwards.
#
# Example HPC invocation (one per scenario):
#   Rscript combine_parallel_results.R \
#     pattern='sim_nrep500_sparse_h2persnp0.03_K1_n1000' \
#     input_dir='/home/apm2217/output'
# =============================================================================

#' Combine parallel simulation results
#'
#' @param input_dir Character, directory containing per-batch result files
#'   (default = "/home/apm2217/output").
#' @param pattern Character, basename pattern that matches the scenario, e.g.
#'   "sim_nrep500_sparse_h2persnp0.03_K1_n1000". The script appends
#'   `_rep[0-9]+-[0-9]+\.rds$` to find the per-batch files.
#' @param output_file Character, output filename (if NULL, auto-generated as
#'   `<pattern>.rds`).
#' @param remove_individual_files Logical, whether to delete per-batch source
#'   files after combining (default = FALSE).
#'
#' @return List with two elements: `replicates` (a length-N list ordered by
#'   replicate_id, NULL slots for any missing reps) and `parameters` (the
#'   shared simulation parameters plus combiner provenance).
combine_parallel_results <- function(input_dir = NULL,
                                    pattern = NULL,
                                    output_file = NULL,
                                    remove_individual_files = FALSE) {

  # ---- Parse command-line arguments ----
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    for (arg in args) {
      split_arg <- strsplit(arg, "=")[[1]]
      key <- split_arg[1]
      value <- split_arg[2]

      # Remove quotes if present
      value <- gsub("^'|'$", "", value)
      value <- gsub('^"|"$', "", value)

      if (key == "input_dir") {
        input_dir <- value
      } else if (key == "pattern") {
        pattern <- value
      } else if (key == "output_file") {
        output_file <- value
      } else if (key == "remove_individual_files") {
        remove_individual_files <- as.logical(value)
      }
    }
  }

  # ---- Set defaults ----
  if (is.null(input_dir)) {
    input_dir <- "/home/apm2217/output"
  }

  if (is.null(pattern)) {
    stop("Pattern must be specified to identify result files to combine")
  }

  # ---- Find matching result files ----
  cat("\n========================================\n")
  cat("Combining Parallel Simulation Results\n")
  cat("========================================\n")
  cat("Input directory:", input_dir, "\n")
  cat("Pattern:", pattern, "\n\n")

  # Find files matching pattern with _rep{X}-{Y}.rds suffix
  all_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

  # Filter for files matching the pattern and having replicate range suffix
  pattern_regex <- paste0(gsub("\\.", "\\\\.", pattern), "_rep[0-9]+-[0-9]+\\.rds$")
  result_files <- grep(pattern_regex, all_files, value = TRUE)

  if (length(result_files) == 0) {
    stop("No result files found matching pattern: ", pattern, "_rep*-*.rds")
  }

  cat("Found", length(result_files), "result files to combine:\n")
  for (f in basename(result_files)) {
    cat("  -", f, "\n")
  }
  cat("\n")

  # ---- Load and validate all results ----
  all_job_results <- list()
  all_replicates <- list()
  replicate_ids <- integer(0)

  for (i in seq_along(result_files)) {
    cat("Loading:", basename(result_files[i]), "...\n")

    job_results <- tryCatch({
      readRDS(result_files[i])
    }, error = function(e) {
      cat("ERROR loading file:", result_files[i], "\n")
      cat("Error message:", e$message, "\n")
      stop("Failed to load result file")
    })

    # Validate structure
    if (!all(c("replicates", "parameters") %in% names(job_results))) {
      stop("Invalid result file structure in: ", basename(result_files[i]))
    }

    # Store job results
    all_job_results[[i]] <- job_results

    # Extract replicates and their IDs
    for (rep in job_results$replicates) {
      if (!is.null(rep)) {
        rep_id <- rep$replicate_id

        # Check for duplicates
        if (rep_id %in% replicate_ids) {
          stop("Duplicate replicate ID found: ", rep_id, " in file: ", basename(result_files[i]))
        }

        replicate_ids <- c(replicate_ids, rep_id)
        all_replicates[[rep_id]] <- rep
      }
    }
  }

  cat("\nLoaded", length(result_files), "job files\n")
  cat("Total replicates collected:", length(replicate_ids), "\n")

  # ---- Validate replicate completeness ----
  replicate_ids_sorted <- sort(replicate_ids)
  expected_total <- max(replicate_ids_sorted)

  cat("\nReplicate range:", min(replicate_ids_sorted), "to", max(replicate_ids_sorted), "\n")

  # Check for missing replicates
  missing_replicates <- setdiff(1:expected_total, replicate_ids_sorted)

  if (length(missing_replicates) > 0) {
    warning("Missing replicates detected: ", paste(missing_replicates, collapse = ", "))
    cat("WARNING: Not all replicates are present. Expected 1 to", expected_total, "\n")
    cat("Missing replicates:", paste(missing_replicates, collapse = ", "), "\n")
  } else {
    cat("All replicates present (1 to", expected_total, ")\n")
  }

  # ---- Validate parameter consistency ----
  cat("\nValidating parameter consistency across jobs...\n")

  # Get parameters from first job
  ref_params <- all_job_results[[1]]$parameters
  param_names <- setdiff(names(ref_params), c("start_replicate", "end_replicate",
                                               "num_replicates", "job_elapsed_time"))

  # Check all other jobs have same parameters
  for (i in 2:length(all_job_results)) {
    job_params <- all_job_results[[i]]$parameters

    for (param in param_names) {
      if (!identical(ref_params[[param]], job_params[[param]])) {
        warning("Parameter mismatch in job ", i, ": ", param)
        cat("  Job 1:", ref_params[[param]], "\n")
        cat("  Job", i, ":", job_params[[param]], "\n")
      }
    }
  }

  cat("Parameter validation complete\n")

  # ---- Combine results ----
  cat("\nCombining results...\n")

  # Create combined results list ordered by replicate ID
  combined_replicates <- vector("list", expected_total)
  for (i in 1:expected_total) {
    if (i %in% replicate_ids) {
      combined_replicates[[i]] <- all_replicates[[i]]
    } else {
      combined_replicates[[i]] <- NULL  # Missing replicate
    }
  }

  # Calculate total elapsed time across all jobs
  total_elapsed_time <- sum(sapply(all_job_results, function(x) {
    if ("job_elapsed_time" %in% names(x$parameters)) {
      return(x$parameters$job_elapsed_time)
    } else {
      return(0)
    }
  }))

  # Create combined results object
  combined_results <- list(
    replicates = combined_replicates,
    parameters = list(
      num_replicates = expected_total,
      h2_per_snp = ref_params$h2_per_snp,
      h2g = ref_params$h2g,
      data_type = ref_params$data_type,
      K = ref_params$K,
      n = ref_params$n,
      prop_h2_sparse = ref_params$prop_h2_sparse,
      prop_h2_oligogenic = ref_params$prop_h2_oligogenic,
      prop_h2_infinitesimal = ref_params$prop_h2_infinitesimal,
      n_oligogenic = ref_params$n_oligogenic,
      n_inf = ref_params$n_inf,
      LD_blocks_dir = ref_params$LD_blocks_dir,
      combined_from_parallel = TRUE,
      num_parallel_jobs = length(result_files),
      total_elapsed_time = total_elapsed_time
    )
  )

  # ---- Save combined results ----
  if (is.null(output_file)) {
    # Auto-generate output filename from pattern
    output_file <- paste0(pattern, ".rds")
  }

  output_path <- file.path(input_dir, output_file)

  cat("\n========================================\n")
  cat("Saving combined results to:", output_path, "\n")
  cat("Total replicates:", expected_total, "\n")
  cat("Parallel jobs combined:", length(result_files), "\n")
  cat("Total time across all jobs:", round(total_elapsed_time / 60, 2), "minutes\n")
  cat("========================================\n")

  saveRDS(combined_results, output_path)

  # ---- Remove individual files if requested ----
  if (remove_individual_files) {
    cat("\nRemoving individual job files...\n")
    removed_count <- 0
    for (f in result_files) {
      if (file.remove(f)) {
        removed_count <- removed_count + 1
        cat("  Removed:", basename(f), "\n")
      } else {
        warning("Failed to remove:", basename(f))
      }
    }
    cat("Removed", removed_count, "of", length(result_files), "files\n")
  }

  cat("\nCombining completed successfully!\n\n")

  return(combined_results)
}

# Run the combine function
results <- combine_parallel_results()
