#!/usr/bin/env Rscript

# Supplementary Figure S1: command generator for the runtime benchmark of the
# lambda-regularized RSS model. Writes one worker command
# (run_rss_lambda_benchmark.R) per (n, p) setting and replicate chunk, and one
# combine command (combine_rss_lambda_results.R) per setting, with the cluster
# paths.
#
# Settings: p = 5000 with n = 100, 250, 500, 1000, the four feature-sample
# ratios drawn in S1; and n = 500 with p = 500, 1000, 2500, which are run but
# not drawn. 100 replicates per setting.
#
# Input:  none
# Output: commands_to_submit_3way.txt, commands_combine_3way.txt
#         (in the working directory)

# ---- Settings ----

parameter_grid <- data.frame(
  n = c(500, 500, 500, 500, 100, 250, 1000),
  p = c(500, 1000, 2500, 5000, 5000, 5000, 5000),
  stringsAsFactors = FALSE
)

total_replicates <- 100
replicates_per_job <- 2    # 50 jobs per setting, 350 in all
lambda <- 1e-5
L <- 5
K <- 3
h2_per_snp <- 0.03
LD_blocks_dir <- "oligo_LD_blocks"
script_path <- "/home/apm2217/data/benchmark/run_rss_lambda_benchmark.R"

# ---- Worker commands ----

commands_file <- "commands_to_submit_3way.txt"
file_conn <- file(commands_file, open = "w")

total_commands <- 0

for (i in 1:nrow(parameter_grid)) {
  n_val <- parameter_grid$n[i]
  p_val <- parameter_grid$p[i]

  num_jobs <- ceiling(total_replicates / replicates_per_job)

  cat("\n========================================\n")
  cat("Setting", i, "of", nrow(parameter_grid), "\n")
  cat("n =", n_val, ", p =", p_val, " (p/n =", p_val / n_val, ")\n")
  cat("Number of parallel jobs:", num_jobs, "\n")
  cat("========================================\n")

  for (job in 1:num_jobs) {
    start_rep <- (job - 1) * replicates_per_job + 1
    end_rep <- min(job * replicates_per_job, total_replicates)

    command <- paste0(
      "Rscript ", script_path,
      " n=", n_val,
      " p=", p_val,
      " lambda=", lambda,
      " L=", L,
      " K=", K,
      " h2_per_snp=", h2_per_snp,
      " total_replicates=", total_replicates,
      " start_replicate=", start_rep,
      " end_replicate=", end_rep,
      " LD_blocks_dir='", LD_blocks_dir, "'"
    )

    writeLines(command, file_conn)
    total_commands <- total_commands + 1

    cat("  Job", job, "of", num_jobs, ": replicates", start_rep, "to", end_rep, "\n")
  }
}

close(file_conn)

cat("\n========================================\n")
cat("Simulation commands file created:", commands_file, "\n")
cat("Total commands:", total_commands, "\n")
cat("Settings:", nrow(parameter_grid), "\n")
cat("========================================\n\n")

# ---- Combine commands ----

cat("Generating combine commands...\n")

combine_commands_file <- "commands_combine_3way.txt"
combine_conn <- file(combine_commands_file, open = "w")

for (i in 1:nrow(parameter_grid)) {
  n_val <- parameter_grid$n[i]
  p_val <- parameter_grid$p[i]

  # Must match base_name in run_rss_lambda_benchmark.R.
  pattern <- sprintf("benchmark_n%d_p%d_lambda%s_K%d_h2persnp%s_nrep%d",
                     n_val, p_val, lambda, K, h2_per_snp, total_replicates)

  combine_cmd <- paste0(
    "Rscript /home/apm2217/data/benchmark/combine_rss_lambda_results.R",
    " pattern='", pattern, "'",
    " input_dir='/home/apm2217/output'"
  )

  writeLines(combine_cmd, combine_conn)
}

close(combine_conn)

cat("Combine commands file created:", combine_commands_file, "\n")
cat("Total combine commands:", nrow(parameter_grid), "\n")
cat("\n========================================\n")
cat("USAGE:\n")
cat("1. Submit all simulation jobs from:", commands_file, "\n")
cat("2. After all jobs complete, run combine commands from:", combine_commands_file, "\n")
cat("========================================\n\n")
