#!/usr/bin/env Rscript

# Fine-mapping simulations (Figure 1 C-F; Supplementary Figures S2, S3, S6-S9
# and S11), step 3 of 3: combiner. Run once every worker job has finished.
#
# Merges the per-job files of one parameter set by replicate_id. Stops on a
# duplicate replicate; warns when a replicate between 1 and the largest id is
# missing (its slot stays NULL) or when a parameter differs between jobs.
#
# Input:  <input_dir>/<pattern>_rep<start>-<end>.rds
#         (simulation_script_parallel.R)
# Output: <input_dir>/<pattern>.rds, a list with
#           replicates  one slot per replicate, in replicate_id order
#           parameters  the shared simulation parameters, the number of jobs
#                       combined and their total run time
#
# Usage: Rscript combine_parallel_results.R \
#          pattern='sim_nrep500_sparse_h2persnp0.03_K1_n1000' \
#          [input_dir='/home/apm2217/output'] [output_file=<pattern>.rds] \
#          [remove_individual_files=FALSE]

combine_parallel_results <- function(input_dir = NULL,
                                    pattern = NULL,
                                    output_file = NULL,
                                    remove_individual_files = FALSE) {

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    for (arg in args) {
      split_arg <- strsplit(arg, "=")[[1]]
      key <- split_arg[1]
      value <- split_arg[2]

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

  if (is.null(input_dir)) {
    input_dir <- "/home/apm2217/output"
  }

  if (is.null(pattern)) {
    stop("Pattern must be specified to identify result files to combine")
  }

  # ---- Find the per-job files ----
  cat("\n========================================\n")
  cat("Combining Parallel Simulation Results\n")
  cat("========================================\n")
  cat("Input directory:", input_dir, "\n")
  cat("Pattern:", pattern, "\n\n")

  all_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

  # Dots in the pattern (as in h2persnp0.03) are matched literally.
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

  # ---- Load and check ----
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

    if (!all(c("replicates", "parameters") %in% names(job_results))) {
      stop("Invalid result file structure in: ", basename(result_files[i]))
    }

    all_job_results[[i]] <- job_results

    for (rep in job_results$replicates) {
      if (!is.null(rep)) {
        rep_id <- rep$replicate_id

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

  replicate_ids_sorted <- sort(replicate_ids)
  expected_total <- max(replicate_ids_sorted)

  cat("\nReplicate range:", min(replicate_ids_sorted), "to", max(replicate_ids_sorted), "\n")

  missing_replicates <- setdiff(1:expected_total, replicate_ids_sorted)

  if (length(missing_replicates) > 0) {
    warning("Missing replicates detected: ", paste(missing_replicates, collapse = ", "))
    cat("WARNING: Not all replicates are present. Expected 1 to", expected_total, "\n")
    cat("Missing replicates:", paste(missing_replicates, collapse = ", "), "\n")
  } else {
    cat("All replicates present (1 to", expected_total, ")\n")
  }

  cat("\nValidating parameter consistency across jobs...\n")

  # Fields that legitimately differ between jobs are left out of the comparison.
  ref_params <- all_job_results[[1]]$parameters
  param_names <- setdiff(names(ref_params), c("start_replicate", "end_replicate",
                                               "num_replicates", "job_elapsed_time"))

  for (i in seq_len(length(all_job_results))[-1]) {
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

  # ---- Combine ----
  cat("\nCombining results...\n")

  # Assigning NULL removes a slot rather than filling it; the indexed
  # assignments that follow restore the length, so a missing replicate still
  # ends up as a NULL slot.
  combined_replicates <- vector("list", expected_total)
  for (i in 1:expected_total) {
    if (i %in% replicate_ids) {
      combined_replicates[[i]] <- all_replicates[[i]]
    } else {
      combined_replicates[[i]] <- NULL
    }
  }

  total_elapsed_time <- sum(sapply(all_job_results, function(x) {
    if ("job_elapsed_time" %in% names(x$parameters)) {
      return(x$parameters$job_elapsed_time)
    } else {
      return(0)
    }
  }))

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

  # ---- Save ----
  if (is.null(output_file)) {
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

results <- combine_parallel_results()
