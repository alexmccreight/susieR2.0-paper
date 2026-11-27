#!/usr/bin/env Rscript

# Supplementary Figure S7, SuSiE+BLiP re-run: command generator for the HPC
# route, for K = 1-5 of the sparse scenario.
#
# Output: <commands_dir>/
#   commands_blip_export.txt   export_blip_inputs.R, one per K x replicate
#                              chunk; submit in parallel
#   commands_blip_arms.txt     run_blip_arms.py, one per K x replicate chunk;
#                              submit once every export job has finished
#   commands_blip_combine.txt  blip_metrics.R, then combine_blip_to_s7.R for
#                              each arm in convert_arms; run in order once
#                              every arms job has finished
#
# Usage: Rscript command_generator_blip.R [key=value ...]
#   scripts_dir         this folder on the cluster
#   sim_dir             raw sparse simulation files
#   ld_blocks_dir       LD block files
#   output_dir          export root (replicate folders, metrics/, s7_<arm>/)
#   python              interpreter with pyblip 0.4.2
#   total_replicates    500
#   replicates_per_job  25
#   k_values            "1,2,3,4,5"
#   n                   1000, LD block rows to keep
#   n_cores             1, forked workers per export job
#   sim_pattern         simulation file name with %d for K
#   arms_args           "", extra run_blip_arms.py options
#   convert_arms        "A_capped,B_uncapped", arms written in the S7 format
#   commands_dir        ".", where the command files are written
# The path defaults below are the cluster layout.

.script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else getwd()
})
source(file.path(.script_dir, "blip_common.R"))

a <- parse_key_value_args(list(
  scripts_dir        = "/home/apm2217/data/benchmark/blip_rerun",
  sim_dir            = "/home/apm2217/data/benchmark/final_scripts/500_rep_results/sparse",
  ld_blocks_dir      = "/home/apm2217/data/benchmark/final_scripts/500_rep_results/ld_blocks",
  output_dir         = "/home/apm2217/output/blip_rerun",
  python             = "/home/apm2217/blip_env/venv/bin/python",
  total_replicates   = 500L,
  replicates_per_job = 25L,
  k_values           = "1,2,3,4,5",
  n                  = 1000L,
  n_cores            = 1L,
  sim_pattern        = "sim_nrep500_sparse_h2persnp0.03_K%d_n1000.rds",
  arms_args          = "",
  convert_arms       = "A_capped,B_uncapped",
  commands_dir       = "."
))

k_values <- as.integer(strsplit(a$k_values, ",")[[1]])
num_jobs <- ceiling(a$total_replicates / a$replicates_per_job)
dir.create(a$commands_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Export and arms commands ----

export_file <- file.path(a$commands_dir, "commands_blip_export.txt")
arms_file   <- file.path(a$commands_dir, "commands_blip_arms.txt")
export_conn <- file(export_file, open = "w")
arms_conn   <- file(arms_file,   open = "w")

total_commands <- 0
for (k in k_values) {
  sim_file <- file.path(a$sim_dir, sprintf(a$sim_pattern, k))
  cat("\n========================================\n")
  cat("K:", k, "\n")
  cat("Sim file:", sim_file, "\n")
  cat("Total replicates:", a$total_replicates, "\n")
  cat("Replicates per job:", a$replicates_per_job, "\n")
  cat("Number of parallel jobs:", num_jobs, "\n")
  cat("========================================\n")

  for (job in seq_len(num_jobs)) {
    start_rep <- (job - 1) * a$replicates_per_job + 1
    end_rep   <- min(job * a$replicates_per_job, a$total_replicates)

    export_cmd <- paste0("Rscript ", file.path(a$scripts_dir, "export_blip_inputs.R"),
                         " sim_file='", sim_file, "'",
                         " ld_blocks_dir='", a$ld_blocks_dir, "'",
                         " output_dir='", a$output_dir, "'",
                         " start_replicate=", start_rep,
                         " end_replicate=", end_rep,
                         " n=", a$n,
                         " n_cores=", a$n_cores)
    arms_cmd <- paste0(a$python, " ", file.path(a$scripts_dir, "run_blip_arms.py"),
                       " --input-dir '", a$output_dir, "'",
                       " --k ", k,
                       " --start-replicate ", start_rep,
                       " --end-replicate ", end_rep,
                       if (nzchar(a$arms_args)) paste0(" ", a$arms_args) else "")

    writeLines(export_cmd, export_conn)
    writeLines(arms_cmd,   arms_conn)
    total_commands <- total_commands + 1
    cat("  Job", job, "of", num_jobs, ": replicates", start_rep, "to", end_rep, "\n")
  }
}
close(export_conn)
close(arms_conn)

# ---- Combine commands ----

convert_arms <- trimws(strsplit(a$convert_arms, ",")[[1]])
combine_file <- file.path(a$commands_dir, "commands_blip_combine.txt")
combine_conn <- file(combine_file, open = "w")
writeLines(paste0("Rscript ", file.path(a$scripts_dir, "blip_metrics.R"),
                  " input_dir='", a$output_dir, "'"), combine_conn)
for (arm in convert_arms) {
  writeLines(paste0("Rscript ", file.path(a$scripts_dir, "combine_blip_to_s7.R"),
                    " input_dir='", a$output_dir, "'",
                    " arm='", arm, "'"), combine_conn)
}
close(combine_conn)

cat("\n========================================\n")
cat("Command files written to", normalizePath(a$commands_dir), "\n")
cat("  commands_blip_export.txt  :", total_commands, "jobs\n")
cat("  commands_blip_arms.txt    :", total_commands, "jobs\n")
cat("  commands_blip_combine.txt :", 1 + length(convert_arms), "commands\n")
cat("USAGE:\n")
cat("1. Submit every line of commands_blip_export.txt\n")
cat("2. When all have finished, submit every line of commands_blip_arms.txt\n")
cat("3. When all have finished, run commands_blip_combine.txt in order\n")
cat("========================================\n\n")
