#!/usr/bin/env bash
# Supplementary Figure S7, SuSiE+BLiP re-run: end-to-end smoke test on
# synthetic data.
#   1. simulate_smoke_data.R builds a small sparse benchmark: 20 replicates,
#      n = 800, p = 1500, AR(1) rho = 0.99, K = 2, h2 = 0.03 per SNP.
#   2. local/run_blip_local.sh runs the local route on it.
#   3. command_generator_blip.R and its three command files, run in order, are
#      the HPC route on the same data (7 replicates per job).
#   4. check_smoke_result.R verifies that the two routes wrote the same
#      per-replicate results, prints the stratified table and checks that the
#      capped arm loses power on replicates whose SuSiE CS exceeds 25 variants
#      while the uncapped arm does not.
#
# Usage: bash run_smoke_test.sh [work dir]   (default: a fresh temp directory)
# Needs the environment from local/setup_local_env.sh, or PYTHON and
# R_LIBS_USER pointing at equivalents. N_REPS [20], RHO [0.99] and N_CORES [4]
# can be set in the environment.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$(dirname "$HERE")"
WORK="${1:-$(mktemp -d -t blip_smoke)}"
N_REPS="${N_REPS:-20}"
RHO="${RHO:-0.99}"
N_CORES="${N_CORES:-4}"
BLIP_ENV_DIR="${BLIP_ENV_DIR:-$HOME/.susieR2_blip_env}"
export PYTHON="${PYTHON:-$BLIP_ENV_DIR/venv/bin/python}"
export R_LIBS_USER="${R_LIBS_USER:-$BLIP_ENV_DIR/Rlib}"
SIM_PATTERN="sim_nrep${N_REPS}_sparse_h2persnp0.03_K%d_n800.rds"
mkdir -p "$WORK"
echo "Smoke test work dir: $WORK"

echo
echo "== 1. simulate =="
Rscript "$HERE/simulate_smoke_data.R" output_dir="$WORK/data" n_reps="$N_REPS" rho="$RHO"

echo
echo "== 2. local route =="
SIM_DIR="$WORK/data/sparse" LD_DIR="$WORK/data/ld_blocks" OUT_DIR="$WORK/out_local" \
  SIM_PATTERN="$SIM_PATTERN" K_VALUES=2 N_TARGET=800 N_CORES="$N_CORES" \
  bash "$SCRIPTS/local/run_blip_local.sh"

echo
echo "== 3. HPC route =="
mkdir -p "$WORK/hpc_commands"
Rscript "$SCRIPTS/command_generator_blip.R" scripts_dir="$SCRIPTS" \
  sim_dir="$WORK/data/sparse" ld_blocks_dir="$WORK/data/ld_blocks" \
  output_dir="$WORK/out_hpc" python="$PYTHON" total_replicates="$N_REPS" \
  replicates_per_job=7 k_values=2 n=800 sim_pattern="$SIM_PATTERN" \
  commands_dir="$WORK/hpc_commands"
for f in commands_blip_export.txt commands_blip_arms.txt commands_blip_combine.txt; do
  echo "-- $f"
  bash -e "$WORK/hpc_commands/$f"
done

echo
echo "== 4. checks =="
Rscript "$HERE/check_smoke_result.R" local_dir="$WORK/out_local" hpc_dir="$WORK/out_hpc"
