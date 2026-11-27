#!/usr/bin/env bash
# Supplementary Figure S7, SuSiE+BLiP re-run: local route. Runs the pipeline
# on one machine: for each K, export_blip_inputs.R (simulation file loaded
# once, replicates forked over N_CORES) and run_blip_arms.py (process pool of
# N_CORES); then blip_metrics.R, and combine_blip_to_s7.R for each arm in
# CONVERT_ARMS. The per-replicate files are identical to those of the HPC
# route. Every step skips replicates that are already done, so the script can
# be re-run after an interruption.
#
# Usage: bash run_blip_local.sh
#        (caffeinate -i bash run_blip_local.sh keeps a Mac awake)
# Environment variables [default]:
#   BLIP_ENV_DIR   environment from setup_local_env.sh     [$HOME/.susieR2_blip_env]
#   PYTHON         python with pyblip                      [$BLIP_ENV_DIR/venv/bin/python]
#   R_LIBS_USER    R library with susieR, jsonlite, dplyr  [$BLIP_ENV_DIR/Rlib]
#   SUSIER2_BENCHMARK_ROOT  raw simulation output ($BENCH) [../susieR2.0-benchmark]
#   SIM_DIR        sparse sim files    [$BENCH/final_scripts/500_rep_results/sparse]
#   LD_DIR         LD block files      [$BENCH/final_scripts/500_rep_results/ld_blocks]
#   OUT_DIR        export root         [$BENCH/final_scripts/500_rep_results/blip_rerun]
#   SIM_PATTERN    file name with %d for K [sim_nrep500_sparse_h2persnp0.03_K%d_n1000.rds]
#   K_VALUES       [1 2 3 4 5]
#   N_TARGET       LD block rows to keep [1000]
#   N_CORES        [cores - 2]
#   ARMS_ARGS      extra run_blip_arms.py options []
#   CONVERT_ARMS   arms written in the S7 format [A_capped B_uncapped]
#   STEPS          [export arms metrics combine]
#
# Output: $OUT_DIR/K<k>/, metrics/, s7_<arm>/ and logs/
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$(dirname "$HERE")"
REPO="$(cd "$SCRIPTS/../../../.." && pwd)"

BLIP_ENV_DIR="${BLIP_ENV_DIR:-$HOME/.susieR2_blip_env}"
PYTHON="${PYTHON:-$BLIP_ENV_DIR/venv/bin/python}"
export R_LIBS_USER="${R_LIBS_USER:-$BLIP_ENV_DIR/Rlib}"
BENCH="${SUSIER2_BENCHMARK_ROOT:-$(dirname "$REPO")/susieR2.0-benchmark}"
SIM_DIR="${SIM_DIR:-$BENCH/final_scripts/500_rep_results/sparse}"
LD_DIR="${LD_DIR:-$BENCH/final_scripts/500_rep_results/ld_blocks}"
OUT_DIR="${OUT_DIR:-$BENCH/final_scripts/500_rep_results/blip_rerun}"
SIM_PATTERN="${SIM_PATTERN:-sim_nrep500_sparse_h2persnp0.03_K%d_n1000.rds}"
K_VALUES="${K_VALUES:-1 2 3 4 5}"
N_TARGET="${N_TARGET:-1000}"
if [ -z "${N_CORES:-}" ]; then
  total="$(sysctl -n hw.ncpu 2>/dev/null || nproc)"
  N_CORES=$(( total > 3 ? total - 2 : 1 ))
fi
ARMS_ARGS="${ARMS_ARGS:-}"
CONVERT_ARMS="${CONVERT_ARMS:-A_capped B_uncapped}"
STEPS="${STEPS:-export arms metrics combine}"

LOG_DIR="$OUT_DIR/logs"
mkdir -p "$LOG_DIR"
stamp() { date "+%Y-%m-%d %H:%M:%S"; }
has_step() { case " $STEPS " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

[ -x "$PYTHON" ] || { echo "Python not found at $PYTHON; run local/setup_local_env.sh first" >&2; exit 1; }
[ -d "$SIM_DIR" ] || { echo "SIM_DIR does not exist: $SIM_DIR" >&2; exit 1; }
[ -d "$LD_DIR" ]  || { echo "LD_DIR does not exist: $LD_DIR" >&2; exit 1; }

echo "[$(stamp)] BLiP re-run (local): K = $K_VALUES, N_CORES = $N_CORES, steps = $STEPS"
echo "  sim dir: $SIM_DIR"
echo "  LD dir:  $LD_DIR"
echo "  out dir: $OUT_DIR"
echo "  python:  $PYTHON"
echo "  R libs:  $R_LIBS_USER"

for k in $K_VALUES; do
  # shellcheck disable=SC2059
  sim_file="$SIM_DIR/$(printf "$SIM_PATTERN" "$k")"
  if has_step export; then
    echo "[$(stamp)] K=$k export"
    Rscript "$SCRIPTS/export_blip_inputs.R" sim_file="$sim_file" ld_blocks_dir="$LD_DIR" \
      output_dir="$OUT_DIR" n="$N_TARGET" n_cores="$N_CORES" 2>&1 | tee "$LOG_DIR/export_K$k.log"
  fi
  if has_step arms; then
    echo "[$(stamp)] K=$k arms"
    # shellcheck disable=SC2086
    "$PYTHON" "$SCRIPTS/run_blip_arms.py" --input-dir "$OUT_DIR" --k "$k" --n-jobs "$N_CORES" \
      --quiet $ARMS_ARGS 2>&1 | tee "$LOG_DIR/arms_K$k.log"
  fi
done

if has_step metrics; then
  echo "[$(stamp)] metrics"
  Rscript "$SCRIPTS/blip_metrics.R" input_dir="$OUT_DIR" n_cores="$N_CORES" 2>&1 | tee "$LOG_DIR/metrics.log"
fi
if has_step combine; then
  for arm in $CONVERT_ARMS; do
    echo "[$(stamp)] combine $arm"
    Rscript "$SCRIPTS/combine_blip_to_s7.R" input_dir="$OUT_DIR" arm="$arm" 2>&1 | tee "$LOG_DIR/combine_$arm.log"
  done
fi
echo "[$(stamp)] Done. Outputs under $OUT_DIR"
