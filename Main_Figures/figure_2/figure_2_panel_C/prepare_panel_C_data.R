#!/usr/bin/env Rscript

# Figure 2, panel C (data): number of signals recovered by both methods or by
# one method only.
#
# A signal is a group of credible sets matched across methods (Jaccard >= 0.75)
# by alphagenome_cs_group_comparison.R. Each signal is assigned to the set of
# methods that report it, and the signals in each set are counted. The methods
# are read from the input, so three-method assignments give the seven groups
# of Supplementary Figure S10.
#
# Usage:  prepare_panel_C_data.R [assignments.rds] [output.rds]
# Input:  data/alphagenome_cs_group_assignments__standard_inf.rds
# Output: figure_2_panel_C/panel_C_data.rds

source("R/paths.R")
source("R/aesthetics.R")

fig2_dir    <- fig_dir(2)
input_file  <- cli_arg(1, file.path(fig2_dir, "data",
                                    "alphagenome_cs_group_assignments__standard_inf.rds"))
output_path <- cli_arg(2, file.path(fig2_dir, "figure_2_panel_C", "panel_C_data.rds"))

signals <- readRDS(input_file)

# Method keys used in the input -> display names, and the short names used in
# the labels of pairwise groups.
display <- setNames(names(METHOD_KEYS), METHOD_KEYS)
short   <- c(standard = "SuSiE", ash = "ash", inf = "inf")
methods <- names(display)[names(display) %in% signals$method]

# ---- Methods reporting each signal ----

sig_ids <- sort(unique(signals$signal_id))
present <- vapply(methods, function(m) {
  sig_ids %in% signals$signal_id[signals$method == m]
}, logical(length(sig_ids)))
present <- matrix(present, nrow = length(sig_ids),
                  dimnames = list(NULL, methods))

# Label of a set of methods: "Consensus" for all of them, the display name for
# one, otherwise the short names joined ("SuSiE&inf"). These are the keys of
# concordance_colors in R/aesthetics.R.
label_for <- function(sel) {
  ms <- methods[sel]
  if (length(ms) == length(methods)) return("Consensus")
  if (length(ms) == 1) return(unname(display[ms]))
  paste(short[ms], collapse = "&")
}

signal_group <- apply(present, 1, label_for)

# ---- Groups, ordered Consensus, pairwise, single-method ----

# Every non-empty set of methods, one per row.
membership <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(methods))))
membership <- membership[rowSums(membership) > 0, , drop = FALSE]
membership <- membership[order(-rowSums(membership),
                               apply(membership, 1, function(r) which(r)[1])), , drop = FALSE]
colnames(membership) <- methods

group_labels <- apply(membership, 1, label_for)
rownames(membership) <- group_labels

concordance_df <- data.frame(
  group     = factor(group_labels, levels = group_labels),
  n_signals = as.integer(table(factor(signal_group, levels = group_labels)))
)

method_totals <- data.frame(
  method = unname(display[methods]),
  total  = vapply(methods, function(m) sum(signals$method == m), integer(1))
)

# Every signal and every credible set belongs to a known method.
stopifnot(sum(concordance_df$n_signals) == length(sig_ids),
          sum(method_totals$total) == nrow(signals))

saveRDS(list(
  concordance      = concordance_df,
  group_membership = membership,                # groups x methods, logical
  methods          = unname(display[methods]),  # display names, plot order
  method_totals    = method_totals
), output_path)
cat("Saved:", output_path, "\n")
