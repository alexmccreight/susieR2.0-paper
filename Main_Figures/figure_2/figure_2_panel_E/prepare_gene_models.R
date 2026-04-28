#!/usr/bin/env Rscript

# =============================================================================
# Figure 2, Panel E — Gene Model Preparation
# =============================================================================
# Fetches protein-coding gene coordinates and exon structures from the Ensembl
# REST API (GRCh38) for the Panel E zoom region, then saves them as an RDS file
# for use by the plotting script.
#
# Region:  chr1:201750000-201970000
# Source:  Ensembl REST API (https://rest.ensembl.org)
# Output:  gene_models.rds
# =============================================================================

library(httr)
library(jsonlite)

# =============================================================================
# Configuration
# =============================================================================

script_dir <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/figure_2/figure_2_panel_E"

# Zoom region (bp)
chrom      <- "1"
region_start <- 201750000
region_end   <- 201970000

ensembl_base <- "https://rest.ensembl.org"

# =============================================================================
# Helper: GET from Ensembl REST API with rate-limit handling
# =============================================================================

ensembl_get <- function(endpoint) {
  url <- paste0(ensembl_base, endpoint)
  cat(sprintf("  GET %s\n", endpoint))

  resp <- GET(url, content_type("application/json"))

  # Respect rate limits
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

# =============================================================================
# Step 1: Get all genes overlapping the region
# =============================================================================

cat("Fetching genes in chr1:", region_start, "-", region_end, "...\n")

overlap <- ensembl_get(sprintf(
  "/overlap/region/homo_sapiens/%s:%d-%d?feature=gene",
  chrom, region_start, region_end
))

# Filter to protein-coding genes only
pc_genes <- overlap[overlap$biotype == "protein_coding", ]
cat(sprintf("Found %d protein-coding genes.\n\n", nrow(pc_genes)))

# =============================================================================
# Step 2: For each gene, fetch canonical transcript exons
# =============================================================================

cat("Fetching exon structures for canonical transcripts...\n")

genes_list <- list()
exons_list <- list()

for (i in seq_len(nrow(pc_genes))) {
  gene_id   <- pc_genes$gene_id[i]
  gene_name <- pc_genes$external_name[i]

  cat(sprintf("  [%d/%d] %s (%s)\n", i, nrow(pc_genes), gene_name, gene_id))

  # Fetch gene with expanded transcripts
  gene_info <- ensembl_get(sprintf("/lookup/id/%s?expand=1", gene_id))

  # Find canonical transcript
  transcripts <- gene_info$Transcript
  canonical_idx <- which(transcripts$is_canonical == 1)
  if (length(canonical_idx) == 0) {
    # Fallback: pick the longest transcript
    canonical_idx <- which.max(transcripts$end - transcripts$start)
    cat(sprintf("    No canonical flag — using longest transcript\n"))
  }
  tx <- transcripts[canonical_idx[1], ]

  cat(sprintf("    Transcript: %s (%d exons)\n", tx$id, nrow(tx$Exon[[1]])))

  # Store gene-level info
  genes_list[[i]] <- data.frame(
    name        = gene_name,
    start       = gene_info$start,
    end         = gene_info$end,
    strand      = gene_info$strand,
    ensembl_id  = gene_id,
    transcript_id = tx$id,
    stringsAsFactors = FALSE
  )

  # Store exon coordinates
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

# Sort genes by start position
genes_df <- genes_df[order(genes_df$start), ]

# =============================================================================
# Step 3: Assign row positions to avoid label overlaps
# =============================================================================

cat("\nAssigning row positions for gene stacking...\n")

# Alternating rows: odd genes (1st, 3rd, 5th) on row 1, even genes on row 2
genes_df$row <- ifelse(seq_len(nrow(genes_df)) %% 2 == 1, 1, 2)

cat(sprintf("  Assigned %d genes across %d rows.\n",
            nrow(genes_df), max(genes_df$row)))

# =============================================================================
# Step 4: Print summary and save
# =============================================================================

cat("\n=== Gene Model Summary ===\n")
for (i in seq_len(nrow(genes_df))) {
  g <- genes_df[i, ]
  n_ex <- sum(exons_df$gene == g$name)
  strand_char <- ifelse(g$strand == 1, "+", "-")
  cat(sprintf("  %s (%s) %s:%d-%d [%s] row=%d, %d exons, tx=%s\n",
              g$name, g$ensembl_id, chrom, g$start, g$end,
              strand_char, g$row, n_ex, g$transcript_id))
}

# Save
output_path <- file.path(script_dir, "gene_models.rds")
saveRDS(list(genes = genes_df, exons = exons_df), output_path)
cat(sprintf("\nSaved: %s\n", output_path))

cat("Done!\n")
