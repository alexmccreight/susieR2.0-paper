#!/usr/bin/env Rscript

# =============================================================================
# Figure 2, Panel B — Data Preparation
# =============================================================================
# Loads the CS group assignments (from the Jaccard >= 0.75 concordance
# analysis), counts signals per cross-method intersection group, and outputs a
# summary for the UpSet-style plot.
#
# Method-agnostic: the set of methods is read from the assignments file, so the
# same script serves main Figure 2 (SuSiE + SuSiE-inf -> 3 groups) and
# Supplementary S11 (all three methods -> 7 groups).
#
# Input:  data/alphagenome_cs_group_assignments__standard_inf.rds
#           (from alphagenome_cs_group_comparison.R; pass a different file as
#            the first command-line argument to rebuild the S11 version)
# Output: panel_B_data.rds
# =============================================================================

# =============================================================================
# Paths
# =============================================================================

source("R/paths.R")
source("R/aesthetics.R")

fig2_dir   <- fig_dir(2)
output_dir <- file.path(fig2_dir, "figure_2_panel_B")

.cli <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(.cli) > 0 && nzchar(.cli[1])) {
  .cli[1]
} else {
  file.path(fig2_dir, "data", "alphagenome_cs_group_assignments__standard_inf.rds")
}
output_path <- if (length(.cli) > 1 && nzchar(.cli[2])) {
  .cli[2]
} else {
  file.path(output_dir, "panel_B_data.rds")
}

# =============================================================================
# Load data
# =============================================================================

cat("Loading group assignments...\n")
cat(sprintf("  %s\n", basename(input_file)))
if (!file.exists(input_file)) stop("Missing input: ", input_file)
signals <- readRDS(input_file)
cat(sprintf("  %d CS entries across %d signals\n",
            nrow(signals), length(unique(signals$signal_id))))

METHOD_ORDER <- c("standard", "ash", "inf")
DISPLAY      <- c(standard = "SuSiE", ash = "SuSiE-ash", inf = "SuSiE-inf")
SHORT        <- c(standard = "SuSiE", ash = "ash",       inf = "inf")

methods <- METHOD_ORDER[METHOD_ORDER %in% unique(signals$method)]
cat(sprintf("  methods present: %s\n", paste(methods, collapse = ", ")))

# =============================================================================
# Signal-level method membership
# =============================================================================

cat("Classifying signals into concordance groups...\n")

sig_ids <- sort(unique(signals$signal_id))
present <- vapply(methods, function(m) {
  sig_ids %in% unique(signals$signal_id[signals$method == m])
}, logical(length(sig_ids)))
present <- matrix(present, nrow = length(sig_ids),
                  dimnames = list(NULL, methods))

# Display label for a method set: all -> "Consensus"; one -> its display name;
# otherwise the short names joined ("SuSiE&inf"), matching the palette keys in
# R/aesthetics.R.
label_for <- function(sel) {
  ms <- methods[sel]
  if (length(ms) == length(methods)) return("Consensus")
  if (length(ms) == 1) return(unname(DISPLAY[ms]))
  paste(unname(SHORT[ms]), collapse = "&")
}

signal_group <- apply(present, 1, label_for)

# =============================================================================
# Group ordering: Consensus, then pairwise, then single-method
# =============================================================================

combos <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(methods))))
combos <- combos[rowSums(combos) > 0, , drop = FALSE]
combos <- combos[order(-rowSums(combos),
                       apply(combos, 1, function(r) which(r)[1])), , drop = FALSE]
colnames(combos) <- methods

group_labels     <- apply(combos, 1, label_for)
group_membership <- combos
rownames(group_membership) <- group_labels

signal_counts <- table(factor(signal_group, levels = group_labels))

concordance_df <- data.frame(
  group     = factor(group_labels, levels = group_labels),
  n_signals = as.integer(signal_counts),
  stringsAsFactors = FALSE
)

method_totals <- data.frame(
  method = unname(DISPLAY[methods]),
  total  = vapply(methods, function(m) sum(signals$method == m), integer(1)),
  stringsAsFactors = FALSE
)

# =============================================================================
# Verification
# =============================================================================

cat("\n--- Signal counts per concordance group ---\n")
for (i in seq_len(nrow(concordance_df))) {
  cat(sprintf("  %-15s %5d signals\n",
              as.character(concordance_df$group[i]), concordance_df$n_signals[i]))
}
cat(sprintf("  %-15s %5d\n", "TOTAL", sum(concordance_df$n_signals)))
stopifnot(sum(concordance_df$n_signals) == length(sig_ids))

cat("\n--- Total CSs per method ---\n")
for (i in seq_len(nrow(method_totals))) {
  cat(sprintf("  %-10s %5d CSs\n", method_totals$method[i], method_totals$total[i]))
}
stopifnot(sum(method_totals$total) == nrow(signals))

# =============================================================================
# Save
# =============================================================================

output <- list(
  concordance      = concordance_df,
  group_membership = group_membership,          # groups x methods, logical
  methods          = unname(DISPLAY[methods]),  # display names, plot order
  method_totals    = method_totals
)

saveRDS(output, output_path)
cat(sprintf("\nSaved: %s\n", output_path))
