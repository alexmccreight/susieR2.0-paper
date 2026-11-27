#!/usr/bin/env Rscript

# Figure 2, panel E (gene models): the protein-coding genes overlapping
# chr1:201,750,000-201,970,000 (GRCh38) and the exons of their canonical
# transcripts, from the Ensembl REST API (https://rest.ensembl.org). Genes are
# put on two alternating rows, in order of start position, so that their
# labels do not overlap in the gene track.
#
# Needs network access, so run_all.R does not run it; its output ships with
# the repository.
#
# Output: gene_models.rds  (list: $genes, $exons; read by plot_panel_E.R)

library(httr)
library(jsonlite)

source("R/paths.R")

script_dir <- file.path(fig_dir(2), "figure_2_panel_E")

chrom      <- "1"
region_start <- 201750000
region_end   <- 201970000

ensembl_base <- "https://rest.ensembl.org"

# GET an endpoint as parsed JSON. A rate-limited request (HTTP 429) is retried
# once, after the delay the server asks for.
ensembl_get <- function(endpoint) {
  url <- paste0(ensembl_base, endpoint)
  cat(sprintf("  GET %s\n", endpoint))

  resp <- GET(url, content_type("application/json"))

  if (status_code(resp) == 429) {
    retry_after <- as.numeric(headers(resp)[["retry-after"]])
    if (is.na(retry_after)) retry_after <- 1
    cat(sprintf("  Rate limited — waiting %g s\n", retry_after))
    Sys.sleep(retry_after)
    resp <- GET(url, content_type("application/json"))
  }

  if (status_code(resp) != 200) {
    stop(sprintf("Ensembl API error %d for %s", status_code(resp), endpoint))
  }

  content(resp, as = "parsed", simplifyVector = TRUE)
}

# ---- Protein-coding genes in the region ----

cat("Fetching genes in chr1:", region_start, "-", region_end, "...\n")

overlap <- ensembl_get(sprintf(
  "/overlap/region/homo_sapiens/%s:%d-%d?feature=gene",
  chrom, region_start, region_end
))

pc_genes <- overlap[overlap$biotype == "protein_coding", ]
cat(sprintf("Found %d protein-coding genes.\n\n", nrow(pc_genes)))

# ---- Exons of each gene's canonical transcript ----

cat("Fetching exon structures for canonical transcripts...\n")

genes_list <- list()
exons_list <- list()

for (i in seq_len(nrow(pc_genes))) {
  gene_id   <- pc_genes$gene_id[i]
  gene_name <- pc_genes$external_name[i]

  cat(sprintf("  [%d/%d] %s (%s)\n", i, nrow(pc_genes), gene_name, gene_id))

  gene_info <- ensembl_get(sprintf("/lookup/id/%s?expand=1", gene_id))

  # The canonical transcript, or the longest one if none is flagged.
  transcripts <- gene_info$Transcript
  canonical_idx <- which(transcripts$is_canonical == 1)
  if (length(canonical_idx) == 0) {
    canonical_idx <- which.max(transcripts$end - transcripts$start)
    cat(sprintf("    No canonical flag — using longest transcript\n"))
  }
  tx <- transcripts[canonical_idx[1], ]

  cat(sprintf("    Transcript: %s (%d exons)\n", tx$id, nrow(tx$Exon[[1]])))

  genes_list[[i]] <- data.frame(
    name        = gene_name,
    start       = gene_info$start,
    end         = gene_info$end,
    strand      = gene_info$strand,
    ensembl_id  = gene_id,
    transcript_id = tx$id,
    stringsAsFactors = FALSE
  )

  exon_data <- tx$Exon[[1]]
  exons_list[[i]] <- data.frame(
    gene  = gene_name,
    start = exon_data$start,
    end   = exon_data$end,
    stringsAsFactors = FALSE
  )

  # Be polite to the API
  Sys.sleep(0.2)
}

genes_df <- do.call(rbind, genes_list)
exons_df <- do.call(rbind, exons_list)

genes_df <- genes_df[order(genes_df$start), ]

# ---- Rows of the gene track ----

cat("\nAssigning row positions for gene stacking...\n")

# Genes alternate between row 1 (1st, 3rd, 5th) and row 2 (2nd, 4th).
genes_df$row <- ifelse(seq_len(nrow(genes_df)) %% 2 == 1, 1, 2)

cat(sprintf("  Assigned %d genes across %d rows.\n",
            nrow(genes_df), max(genes_df$row)))

# ---- Summary and save ----

cat("\n=== Gene Model Summary ===\n")
for (i in seq_len(nrow(genes_df))) {
  g <- genes_df[i, ]
  n_ex <- sum(exons_df$gene == g$name)
  strand_char <- ifelse(g$strand == 1, "+", "-")
  cat(sprintf("  %s (%s) %s:%d-%d [%s] row=%d, %d exons, tx=%s\n",
              g$name, g$ensembl_id, chrom, g$start, g$end,
              strand_char, g$row, n_ex, g$transcript_id))
}

output_path <- file.path(script_dir, "gene_models.rds")
saveRDS(list(genes = genes_df, exons = exons_df), output_path)
cat(sprintf("\nSaved: %s\n", output_path))

cat("Done!\n")
