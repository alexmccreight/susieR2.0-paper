#!/usr/bin/env bash

# Supplementary Figure S4: runs the cross-validated prediction benchmark on one
# machine, with the scripts the cluster command files run. Steps (STEPS):
#   generate  command_generator_parallel_cv.R writes the worker and combine
#             commands with local paths: the sparse scenarios go to
#             $OUT_DIR/sparse and the three complex scenarios to
#             $OUT_DIR/complex, the layout
#             Supplementary_Figures/S4/extract_prediction_data.R reads;
#   run       the worker commands run N_JOBS at a time through xargs, each a
#             separate R process on REPS_PER_JOB replicates with one BLAS
#             thread;
#   combine   combine_parallel_cv_results.R merges each parameter set; then one
#             summary line per combined file (replicates present, failed fits,
#             non-finite mean R2);
#   figure    extract_prediction_data.R and plot_S4.R rebuild the figure.
#
# A worker job writes its result file through a rename when its replicates are
# done, and exits at once on a re-run if that file exists. An interrupted run
# therefore resumes with the same command and loses, per worker, only the job
# in progress (one replicate with the default REPS_PER_JOB).
# logs/progress.log gets one line per finished job; logs/failed_commands.txt
# lists the failed jobs, whose logs are kept as logs/job.*; the driver stops
# before combining if any job failed.
#
# Expected cost with the current susieR: about 3.5 to 6 minutes per replicate
# (5 folds x 3 fits), 2400 replicates, so roughly 150 CPU-hours. A worker peaks
# at about 4.3 GB, so a 48 GB machine runs 8 jobs at once and the whole study
# takes about 19 hours.
#
# Usage: caffeinate -ims bash run_prediction_cv_local.sh
#        (on a Mac, plugged in and with the lid open: caffeinate does not
#        prevent sleep on lid close)
# Environment variables (defaults in brackets):
#   SUSIER2_ENV_DIR   environment from fine_mapping/blip_rerun/local/setup_local_env.sh
#                     [$HOME/.susieR2_blip_env]
#   R_LIBS_USER       R library with the development susieR [$SUSIER2_ENV_DIR/Rlib]
#   SUSIER2_BENCHMARK_ROOT  benchmark tree [../susieR2.0-benchmark]
#   LD_DIR            LD block files
#                     [$SUSIER2_BENCHMARK_ROOT/final_scripts/500_rep_results/ld_blocks]
#   OUT_DIR           [$SUSIER2_BENCHMARK_ROOT/final_scripts/prediction_cv_results]
#   N_JOBS            parallel worker processes
#                     [min(cores - 2, (RAM in GB - 8) / MEM_PER_JOB_GB), at least 1]
#   MEM_PER_JOB_GB    memory budget per worker used for that bound [5]
#   REPS_PER_JOB      replicates per worker job [1]
#   TOTAL_REPLICATES  replicates per parameter set [300]
#   CONVERGENCE       convergence_method for SuSiE-inf and SuSiE-ash [pip]
#   STEPS             [generate run combine figure]
#
# Input:  $LD_DIR/*.rds
# Output: $OUT_DIR/{sparse,complex}/pred_cv_*.rds, $OUT_DIR/commands/,
#         $OUT_DIR/logs/, and the files the two S4 scripts write

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$(dirname "$HERE")"
REPO="$(cd "$SCRIPTS/../../.." && pwd)"

SUSIER2_ENV_DIR="${SUSIER2_ENV_DIR:-$HOME/.susieR2_blip_env}"
export R_LIBS_USER="${R_LIBS_USER:-$SUSIER2_ENV_DIR/Rlib}"
BENCH="${SUSIER2_BENCHMARK_ROOT:-$(dirname "$REPO")/susieR2.0-benchmark}"
export SUSIER2_BENCHMARK_ROOT="$BENCH"
LD_DIR="${LD_DIR:-$BENCH/final_scripts/500_rep_results/ld_blocks}"
OUT_DIR="${OUT_DIR:-$BENCH/final_scripts/prediction_cv_results}"
MEM_PER_JOB_GB="${MEM_PER_JOB_GB:-5}"
if [ -z "${N_JOBS:-}" ]; then
  # A worker peaks at ~4.3 GB resident (SuSiE-ash on p = 5000), so the job
  # count is bounded by memory as well as by cores, leaving 8 GB for the rest.
  cores="$(sysctl -n hw.ncpu 2>/dev/null || nproc)"
  mem_gb="$(( $(sysctl -n hw.memsize 2>/dev/null || grep MemTotal /proc/meminfo | awk '{print $2 * 1024}') / 1073741824 ))"
  by_cores=$(( cores > 3 ? cores - 2 : 1 ))
  by_mem=$(( (mem_gb - 8) / MEM_PER_JOB_GB ))
  N_JOBS=$(( by_mem < by_cores ? by_mem : by_cores ))
  [ "$N_JOBS" -ge 1 ] || N_JOBS=1
fi
REPS_PER_JOB="${REPS_PER_JOB:-1}"
TOTAL_REPLICATES="${TOTAL_REPLICATES:-300}"
CONVERGENCE="${CONVERGENCE:-pip}"
STEPS="${STEPS:-generate run combine figure}"

# One BLAS thread per worker process; the parallelism is across processes.
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

LOG_DIR="$OUT_DIR/logs"
CMD_DIR="$OUT_DIR/commands"
mkdir -p "$OUT_DIR/sparse" "$OUT_DIR/complex" "$LOG_DIR" "$CMD_DIR"
stamp() { date "+%Y-%m-%d %H:%M:%S"; }
has_step() { case " $STEPS " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

[ -d "$LD_DIR" ] || { echo "LD_DIR does not exist: $LD_DIR" >&2; exit 1; }
Rscript -e 'invisible(library(susieR)); cat("susieR", as.character(packageVersion("susieR")), "from", find.package("susieR"), "\n")'

echo "[$(stamp)] prediction CV (local): N_JOBS = $N_JOBS, $REPS_PER_JOB replicates per job, $TOTAL_REPLICATES replicates per scenario"
echo "  LD dir:  $LD_DIR"
echo "  out dir: $OUT_DIR"
echo "  R libs:  $R_LIBS_USER"

if has_step generate; then
  echo "[$(stamp)] generating commands"
  Rscript "$SCRIPTS/command_generator_parallel_cv.R" scripts_dir="$SCRIPTS" ld_blocks_dir="$LD_DIR" \
    output_dir="$OUT_DIR/sparse" scenarios=sparse total_replicates="$TOTAL_REPLICATES" \
    replicates_per_job="$REPS_PER_JOB" convergence_method="$CONVERGENCE" skip_existing=TRUE \
    commands_dir="$CMD_DIR/sparse" > "$LOG_DIR/generate_sparse.log"
  Rscript "$SCRIPTS/command_generator_parallel_cv.R" scripts_dir="$SCRIPTS" ld_blocks_dir="$LD_DIR" \
    output_dir="$OUT_DIR/complex" scenarios=complex,complex_s1,complex_s2 total_replicates="$TOTAL_REPLICATES" \
    replicates_per_job="$REPS_PER_JOB" convergence_method="$CONVERGENCE" skip_existing=TRUE \
    commands_dir="$CMD_DIR/complex" > "$LOG_DIR/generate_complex.log"
  cat "$CMD_DIR/sparse/commands_cv_to_submit.txt" "$CMD_DIR/complex/commands_cv_to_submit.txt" > "$CMD_DIR/all_jobs.txt"
  cat "$CMD_DIR/sparse/commands_cv_combine.txt" "$CMD_DIR/complex/commands_cv_combine.txt" > "$CMD_DIR/all_combine.txt"
  echo "  $(wc -l < "$CMD_DIR/all_jobs.txt") worker jobs, $(wc -l < "$CMD_DIR/all_combine.txt") combine commands"
fi

if has_step run; then
  echo "[$(stamp)] running worker jobs ($N_JOBS in parallel; finished jobs exit at once)"
  rm -f "$LOG_DIR/failed_commands.txt"
  export LOG_DIR
  # Each line is one full command, handed to run_one_job.sh unparsed.
  tr '\n' '\0' < "$CMD_DIR/all_jobs.txt" | xargs -0 -P "$N_JOBS" -n 1 bash "$HERE/run_one_job.sh" || true
  n_done=$(find "$OUT_DIR/sparse" "$OUT_DIR/complex" -name 'pred_cv_*_rep*-*.rds' | wc -l | tr -d ' ')
  echo "[$(stamp)] worker jobs finished: $n_done output files of $(wc -l < "$CMD_DIR/all_jobs.txt") expected"
  if [ -s "$LOG_DIR/failed_commands.txt" ]; then
    echo "  $(wc -l < "$LOG_DIR/failed_commands.txt") jobs failed; see $LOG_DIR/failed_commands.txt and the job.*.log files" >&2
    exit 1
  fi
fi

if has_step combine; then
  echo "[$(stamp)] combining"
  while IFS= read -r cmd; do eval "$cmd" > /dev/null; done < "$CMD_DIR/all_combine.txt"
  # Completeness and fit-failure summary of every combined scenario file.
  Rscript -e '
    for (f in Sys.glob(c("'"$OUT_DIR"'/sparse/pred_cv_*nfolds*.rds", "'"$OUT_DIR"'/complex/pred_cv_*nfolds*.rds"))) {
      if (grepl("_rep[0-9]+-[0-9]+[.]rds$", f)) next
      d <- readRDS(f); reps <- Filter(Negate(is.null), d$replicates)
      n_err <- sum(vapply(reps, function(r) if (is.null(r$n_fit_errors)) 0L else r$n_fit_errors, integer(1)))
      n_na  <- sum(vapply(reps, function(r) sum(!is.finite(vapply(r$cv_metrics, function(m) m$mean_r2, 1))), integer(1)))
      cat(sprintf("  %-95s %3d/%d replicates, %d failed fits, %d non-finite mean R2\n", basename(f), length(reps), length(d$replicates), n_err, n_na))
    }'
fi

if has_step figure; then
  echo "[$(stamp)] extracting S4 data and plotting"
  cd "$REPO"
  Rscript Supplementary_Figures/S4/extract_prediction_data.R 2>&1 | tee "$LOG_DIR/extract_S4.log"
  Rscript Supplementary_Figures/S4/plot_S4.R 2>&1 | tee "$LOG_DIR/plot_S4.log"
fi
echo "[$(stamp)] Done."
