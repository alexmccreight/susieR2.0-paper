#!/usr/bin/env Rscript

# Figure 2, panel E (data): PIPs and credible sets of SuSiE and SuSiE-inf at
# the LMOD1 locus (ENSG00000163431, chr1) in DLPFC.
#
# For each method, takes the PIP, position and credible-set membership of every
# variant from the fitted model, and one summary row per credible set. Credible
# sets are numbered in the order in which the fit lists them.
#
# Usage:  prepare_panel_E_data.R [--include-ash]
#         --include-ash adds SuSiE-ash and writes to
#         Supplementary_Figures/S10/data instead (the three-method panel of
#         Supplementary Figure S10).
# Input:  data/panel_E/{STANDARD,INF,ASH}-ROSMAP_DeJager_eQTL_twas.chr1_
#           ENSG00000163431.univariate_bvsr.rds
#           (the pipeline's output for this gene, without its per-sample fields)
# Output: figure_2_panel_E/panel_E_data.rds
#         (list: $pip_data, $cs_summary, $region_info, $gene_id, $tissue)

source("R/paths.R")
source("R/aesthetics.R")

args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% "--include-ash")) {
  stop("Unknown argument(s): ", paste(args[!args %in% "--include-ash"], collapse = ", "))
}
include_ash <- "--include-ash" %in% args

fig2_dir   <- fig_dir(2)
data_dir   <- file.path(fig2_dir, "data", "panel_E")
output_dir <- if (include_ash) {
  file.path(supp_dir(10), "data")
} else {
  file.path(fig2_dir, "figure_2_panel_E")
}

gene_id    <- "ENSG00000163431"
tissue_key <- "DLPFC_DeJager_eQTL_ENSG00000163431"
methods    <- if (include_ash) ALL_METHODS else MAIN_METHODS

# The entry of one method's fit file for this gene and tissue.
read_locus <- function(method) {
  fit_file <- sprintf("%s-ROSMAP_DeJager_eQTL_twas.chr1_%s.univariate_bvsr.rds",
                      toupper(METHOD_KEYS[[method]]), gene_id)
  readRDS(file.path(data_dir, fit_file))[[gene_id]][[tissue_key]]
}
loci <- lapply(setNames(methods, methods), read_locus)

all_variants <- list()
cs_summaries <- list()

for (method in methods) {
  fit       <- loci[[method]]$susie_result_trimmed
  var_names <- loci[[method]]$variant_names

  df <- data.frame(
    variant_name = var_names,
    # Variant names are chr:pos:ref:alt.
    position     = as.numeric(sapply(strsplit(var_names, ":"), `[`, 2)),
    pip          = fit$pip,
    in_cs        = FALSE,
    cs_id        = NA_integer_,
    method       = method
  )

  for (cs_idx in seq_along(fit$sets$cs)) {
    cs_vars <- fit$sets$cs[[cs_idx]]
    cs_pips <- fit$pip[cs_vars]
    df$in_cs[cs_vars] <- TRUE
    df$cs_id[cs_vars] <- cs_idx

    cs_summaries[[length(cs_summaries) + 1]] <- data.frame(
      method       = method,
      cs_id        = cs_idx,
      n_variants   = length(cs_vars),
      lead_variant = var_names[cs_vars[which.max(cs_pips)]],
      lead_pip     = max(cs_pips),
      coverage     = fit$sets$coverage[cs_idx],
      min_abs_corr = fit$sets$purity[cs_idx, 1]
    )
  }

  all_variants[[method]] <- df
}

pip_data <- do.call(rbind, all_variants)
pip_data$method <- factor(pip_data$method, levels = methods)
rownames(pip_data) <- NULL

cs_summary <- do.call(rbind, cs_summaries)
cs_summary$method <- factor(cs_summary$method, levels = methods)
rownames(cs_summary) <- NULL

output_path <- file.path(output_dir, "panel_E_data.rds")
saveRDS(list(
  pip_data    = pip_data,
  cs_summary  = cs_summary,
  region_info = loci[["SuSiE"]]$region_info,
  gene_id     = gene_id,
  tissue      = "DLPFC"
), output_path)
cat("Saved:", output_path, "\n")
