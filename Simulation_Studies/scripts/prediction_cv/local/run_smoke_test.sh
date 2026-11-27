#!/usr/bin/env bash

# Supplementary Figure S4: smoke test of the prediction CV pipeline. Runs one
# replicate of the sparse K = 1 scenario and one of the extensive infinitesimal
# scenario (two worker jobs in parallel, 5 folds x 3 fits each), combines them,
# and prints the per-method cross-validated R^2 and the parameters recorded in
# the output. Fails if a mean R^2 is not finite. Takes about 5 minutes.
#
# Usage: bash run_smoke_test.sh [work dir]   (default: a fresh temp directory)
# Needs the environment from fine_mapping/blip_rerun/local/setup_local_env.sh
# and the LD blocks (LD_DIR, default as in run_prediction_cv_local.sh).
#
# Input:  $LD_DIR/*.rds
# Output: <work dir>/out/*.rds, <work dir>/commands/, <work dir>/generate.log

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$(dirname "$HERE")"
REPO="$(cd "$SCRIPTS/../../.." && pwd)"
WORK="${1:-$(mktemp -d -t predcv_smoke)}"
SUSIER2_ENV_DIR="${SUSIER2_ENV_DIR:-$HOME/.susieR2_blip_env}"
export R_LIBS_USER="${R_LIBS_USER:-$SUSIER2_ENV_DIR/Rlib}"
BENCH="${SUSIER2_BENCHMARK_ROOT:-$(dirname "$REPO")/susieR2.0-benchmark}"
LD_DIR="${LD_DIR:-$BENCH/final_scripts/500_rep_results/ld_blocks}"
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$WORK"
echo "Smoke test work dir: $WORK"

Rscript "$SCRIPTS/command_generator_parallel_cv.R" scripts_dir="$SCRIPTS" ld_blocks_dir="$LD_DIR" \
  output_dir="$WORK/out" scenarios=sparse,complex_s2 total_replicates=1 replicates_per_job=1 skip_existing=TRUE \
  commands_dir="$WORK/commands" > "$WORK/generate.log"
grep -E "data_type='sparse' h2_per_snp=0.03 K=1 |n_oligogenic=10 " "$WORK/commands/commands_cv_to_submit.txt" > "$WORK/commands/smoke_jobs.txt"
echo "Running $(wc -l < "$WORK/commands/smoke_jobs.txt") worker jobs in parallel..."
export LOG_DIR="$WORK"
tr '\n' '\0' < "$WORK/commands/smoke_jobs.txt" | xargs -0 -P 2 -n 1 bash "$HERE/run_one_job.sh"

while IFS= read -r cmd; do eval "$cmd" > /dev/null; done < <(grep -E "K1_|nOligo10" "$WORK/commands/commands_cv_combine.txt")

Rscript -e '
files <- list.files("'"$WORK"'/out", pattern = "nfolds5[.]rds$", full.names = TRUE)
for (f in files) {
  d <- readRDS(f); r <- d$replicates[[1]]
  cat("\n", basename(f), "\n  susieR", d$parameters$susieR_version, "| convergence", d$parameters$convergence_method,
      "| LD block", r$ld_block_name, "| causal", length(r$causal_indices), "| elapsed", round(r$elapsed_time), "s\n")
  for (m in names(r$cv_metrics)) cat(sprintf("  %-10s mean R2 %.4f (sd %.4f)  mean RMSE %.4f\n", m, r$cv_metrics[[m]]$mean_r2, r$cv_metrics[[m]]$sd_r2, r$cv_metrics[[m]]$mean_rmse))
  stopifnot(all(vapply(r$cv_metrics, function(m) is.finite(m$mean_r2), logical(1))))
}
cat("\nSMOKE TEST PASSED\n")'
