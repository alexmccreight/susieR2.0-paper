#!/usr/bin/env Rscript

# Figure 2, panel A (data): the 95% credible sets of SuSiE and SuSiE-inf in
# the ROSMAP eQTL analysis, each with its size bin.
#
# Only gene-tissue pairs present for all three methods (SuSiE, SuSiE-ash,
# SuSiE-inf) are kept, so that Figure 2 and the three-method Supplementary
# Figure S10 count credible sets over the same pairs.
#
# Usage:  prepare_panel_A_data.R [--include-ash]
#         --include-ash adds SuSiE-ash and writes to
#         Supplementary_Figures/S10/data instead (the three-method panel of
#         Supplementary Figure S10).
# Input:  data/panel_A/susie_{standard,ash,inf}_gene_tissue_summary.rds
#         data/panel_A/susie_{standard,ash,inf}_credible_sets.rds
# Output: figure_2_panel_A/panel_A_data.rds

source("R/paths.R")
source("R/aesthetics.R")

args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% "--include-ash")) {
  stop("Unknown argument(s): ", paste(args[!args %in% "--include-ash"], collapse = ", "))
}
include_ash <- "--include-ash" %in% args

fig2_dir   <- fig_dir(2)
data_dir   <- file.path(fig2_dir, "data", "panel_A")
output_dir <- if (include_ash) {
  file.path(supp_dir(10), "data")
} else {
  file.path(fig2_dir, "figure_2_panel_A")
}
output_path <- file.path(output_dir, "panel_A_data.rds")
methods     <- if (include_ash) ALL_METHODS else MAIN_METHODS

read_method <- function(method, what) {
  readRDS(file.path(data_dir, sprintf("susie_%s_%s.rds", METHOD_KEYS[[method]], what)))
}
make_pair <- function(df) paste(df$gene_id, df$tissue, sep = ":")

common_pairs <- Reduce(intersect, lapply(ALL_METHODS, function(method) {
  make_pair(read_method(method, "gene_tissue_summary"))
}))

cs_combined <- do.call(rbind, lapply(methods, function(method) {
  cs <- read_method(method, "credible_sets")
  cs <- cs[make_pair(cs) %in% common_pairs & cs$coverage_level == "0.95", ]
  cs$method <- method
  cs
}))
rownames(cs_combined) <- NULL

cs_combined$cs_size_bin <- cut(
  cs_combined$n_variants,
  breaks = c(0, 1, 2, 5, 10, Inf),
  labels = c("1", "2", "3-5", "6-10", ">10")
)
cs_combined$method <- factor(cs_combined$method, levels = methods)

saveRDS(cs_combined, output_path)
cat("Saved:", output_path, "\n")
