#!/usr/bin/env bash

# Supplementary Figure S4, prediction CV on one machine: runs one worker
# command and keeps its log only if the command fails. Called by
# run_prediction_cv_local.sh and run_smoke_test.sh through `xargs -0 -n 1`,
# which hands over each command line unparsed.
#
# Usage:  bash run_one_job.sh "<command>"   (LOG_DIR: log directory, default .)
# Output: $LOG_DIR/progress.log, one line per successful command (a job that
#         exits at once because its output exists counts); for a failed
#         command, a line in $LOG_DIR/failed_commands.txt and its log in
#         $LOG_DIR/job.XXXXXX

cmd="$1"
log_dir="${LOG_DIR:-.}"
log="$(mktemp "$log_dir/job.XXXXXX")"
if eval "$cmd" > "$log" 2>&1; then
  rm -f "$log"
  # The line holds the scenario arguments and the replicate range.
  printf '%s ok %s\n' "$(date '+%F %T')" \
    "$(echo "$cmd" | sed -E "s/.* (start_replicate=[0-9]+) (end_replicate=[0-9]+) (data_type=.*) n=[0-9]+ n_folds.*/\3 \1 \2/")" \
    >> "$log_dir/progress.log"
else
  echo "$cmd" >> "$log_dir/failed_commands.txt"
  echo "FAILED (log: $log): $cmd" >&2
  exit 1
fi
