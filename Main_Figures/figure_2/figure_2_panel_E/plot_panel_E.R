#!/usr/bin/env Rscript

# =============================================================================
# Figure 2, Panel E — Publication-ready slimmed version
# =============================================================================
# Tracks selected to tell the key story:
#   Red CS sentinel (201806618) overlaps Exc. Neuron H3K27ac (invisible in bulk)
#   Blue CS ASH/INF (201917641) overlaps bulk H3K27ac,
#     astrocyte DNase, astrocyte H3K4me1, astrocyte H3K4me3
#   Blue CS SuSiE (201915540) overlaps Exc. Neuron H3K27me3 (REPRESSIVE)
#
# Input:  panel_E_data.rds, celltype_tracks/*.bed, bulk BED files
# Output: panel_E.pdf, panel_E.png, panel_E_plot.rds
# =============================================================================

library(ggplot2)
library(cowplot)

# =============================================================================
# Paths
# =============================================================================

fig2_dir   <- "/Users/alexmccreight/StatFunGen/susieR2.0-benchmark/final_scripts/figure_2"
script_dir <- file.path(fig2_dir, "figure_2_panel_E")
data_dir   <- file.path(fig2_dir, "data", "panel_E")
data_path  <- file.path(script_dir, "panel_E_data.rds")
ct_dir     <- file.path(data_dir, "celltype_tracks")

# =============================================================================
# Signal-based CS color mapping
# =============================================================================

signal_map <- data.frame(
  method = c("SuSiE",  "SuSiE",  "SuSiE-ash", "SuSiE-ash", "SuSiE-inf"),
  cs_id  = c(1,        2,        1,            2,            1),
  signal = c("A",      "C",      "B",          "A",          "A"),
  stringsAsFactors = FALSE
)

signal_colors <- c(
  "A" = "#1E88E5",   # blue   — shared Standard CS1/INF CS1/ASH CS2
  "B" = "#E53935",   # red    — ASH-unique novel CS1
  "C" = "#43A047"    # green  — Standard diffuse CS2
)

# =============================================================================
# Load data
# =============================================================================

cat("Loading panel G data...\n")
panel_data <- readRDS(data_path)
pip_data   <- panel_data$pip_data
cs_summary <- panel_data$cs_summary

pip_data$pos_mb <- pip_data$position / 1e6

pip_data$signal <- NA_character_
for (i in seq_len(nrow(signal_map))) {
  mask <- pip_data$method == signal_map$method[i] &
          !is.na(pip_data$cs_id) &
          pip_data$cs_id == signal_map$cs_id[i]
  pip_data$signal[mask] <- signal_map$signal[i]
}

# =============================================================================
# Zoom window & key positions
# =============================================================================

zoom_start <- 201.78   # Mb
zoom_end   <- 201.97   # Mb

ash_lead_pos_mb <- 201806618 / 1e6
ash_cs1_positions_mb <- c(201799661, 201799712, 201801928,
                          201805081, 201805083, 201806618) / 1e6
shared_sentinel_pos_mb <- 201917641 / 1e6
susie_cs1_pos_mb <- 201915540 / 1e6

# =============================================================================
# Helper: build one PIP track
# =============================================================================

make_pip_track <- function(df, method_name, show_x_axis = FALSE) {

  d <- df[df$method == method_name, ]
  d_bg <- d[!d$in_cs, ]
  d_cs <- d[d$in_cs, ]

  p <- ggplot() +
    geom_vline(xintercept = ash_lead_pos_mb, linetype = "dashed",
               color = "#E53935", linewidth = 0.8, alpha = 0.7) +
    geom_vline(xintercept = shared_sentinel_pos_mb, linetype = "dashed",
               color = "#1E88E5", linewidth = 0.8, alpha = 0.7) +
    geom_vline(xintercept = susie_cs1_pos_mb, linetype = "dotdash",
               color = "#90CAF9", linewidth = 0.8, alpha = 0.8) +
    geom_point(data = d_bg, aes(x = pos_mb, y = pip),
               size = 1.2, color = "gray70", alpha = 0.5) +
    geom_point(data = d_cs, aes(x = pos_mb, y = pip, color = signal),
               shape = 21, fill = NA, size = 8, stroke = 2.0,
               alpha = 0.9) +
    geom_point(data = d_cs, aes(x = pos_mb, y = pip),
               size = 4, color = "black", alpha = 0.9) +
    scale_color_manual(values = signal_colors, guide = "none") +
    scale_y_continuous(limits = c(0, 1.08), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_continuous(limits = c(zoom_start, zoom_end)) +
    labs(y = "PIP") +
    theme_classic(base_size = 20) +
    theme(
      axis.title.y = element_text(face = "bold", size = 19),
      axis.text.y  = element_text(color = "black", size = 17),
      plot.title   = element_text(face = "bold", size = 20, hjust = 0),
      plot.margin  = margin(2, 10, 2, 5)
    )

  p <- p + ggtitle(method_name)

  if (show_x_axis) {
    p <- p + labs(x = "Position on chr1 (Mb)") +
      theme(
        axis.title.x = element_text(face = "bold", size = 19),
        axis.text.x  = element_text(color = "black", size = 17)
      )
  } else {
    p <- p + labs(x = NULL) +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank()
      )
  }

  return(p)
}

# =============================================================================
# Build 3 PIP tracks
# =============================================================================

cat("Creating PIP tracks...\n")

track_susie <- make_pip_track(pip_data, "SuSiE", show_x_axis = FALSE)

track_susie <- track_susie +
  annotate("segment",
           x = 201.88, xend = susie_cs1_pos_mb + 0.001,
           y = 0.80, yend = 0.996,
           linewidth = 0.4, color = "gray30") +
  annotate("label",
           x = 201.88, y = 0.80,
           label = 'atop("chr1:201915540:A:G", max~bgroup("|", Q[RNA], "|") == 0.125)',
           parse = TRUE,
           hjust = 0.5, vjust = 1,
           size = 4.2, fontface = "bold",
           fill = "white", label.r = unit(0.1, "lines"),
           label.padding = unit(0.15, "lines"),
           color = "gray20")

track_ash <- make_pip_track(pip_data, "SuSiE-ash", show_x_axis = FALSE)

track_ash <- track_ash +
  annotate("segment",
           x = 201.85, xend = ash_lead_pos_mb + 0.003,
           y = 0.40, yend = 0.37,
           linewidth = 0.4, color = "gray30") +
  annotate("label",
           x = 201.85, y = 0.40,
           label = 'atop("chr1:201806618:G:A", max~bgroup("|", Q[RNA], "|") == 0.910)',
           parse = TRUE,
           hjust = 0.5, vjust = 0,
           size = 5, fontface = "bold",
           fill = "white", label.r = unit(0.1, "lines"),
           label.padding = unit(0.15, "lines"),
           color = "gray20")

track_inf <- make_pip_track(pip_data, "SuSiE-inf", show_x_axis = FALSE)

track_inf <- track_inf +
  annotate("segment",
           x = 201.88, xend = shared_sentinel_pos_mb - 0.003,
           y = 0.80, yend = 1.0,
           linewidth = 0.4, color = "gray30") +
  annotate("label",
           x = 201.88, y = 0.80,
           label = 'atop("chr1:201917641:TA:T", max~bgroup("|", Q[RNA], "|") == 0.999)',
           parse = TRUE,
           hjust = 0.5, vjust = 1,
           size = 4.2, fontface = "bold",
           fill = "white", label.r = unit(0.1, "lines"),
           label.padding = unit(0.15, "lines"),
           color = "gray20")

# Shared vertical reference lines for annotation tracks
vline_layer <- list(
  geom_vline(xintercept = ash_lead_pos_mb, linetype = "dashed",
             color = "#E53935", linewidth = 0.7, alpha = 0.7),
  geom_vline(xintercept = shared_sentinel_pos_mb, linetype = "dashed",
             color = "#1E88E5", linewidth = 0.7, alpha = 0.7),
  geom_vline(xintercept = susie_cs1_pos_mb, linetype = "dotted",
             color = "#90CAF9", linewidth = 0.7, alpha = 0.8)
)

# =============================================================================
# Gene track — fishbone style (exon blocks + intron lines + strand chevrons)
# =============================================================================

cat("Building gene track...\n")

# Load gene models (generated by prepare_gene_models.R from Ensembl REST API)
gene_models <- readRDS(file.path(script_dir, "gene_models.rds"))
genes_df    <- gene_models$genes
exons_df    <- gene_models$exons

# Convert bp to Mb for exons and add row assignment
exons_df$start <- exons_df$start / 1e6
exons_df$end   <- exons_df$end / 1e6
exons_df$row   <- genes_df$row[match(exons_df$gene, genes_df$name)]

# Assign y positions: row 1 = upper, row 2 = lower
exon_h <- 0.14  # half-height of exon blocks
genes_df$y <- ifelse(genes_df$row == 1, 0.62, 0.22)
exons_df$y <- ifelse(exons_df$row == 1, 0.62, 0.22)

# Clip gene spans to zoom window (in Mb)
genes_df$start_mb <- pmax(genes_df$start / 1e6, zoom_start)
genes_df$end_mb   <- pmin(genes_df$end / 1e6, zoom_end)
genes_df$label_x  <- pmin((genes_df$start_mb + genes_df$end_mb) / 2,
                          zoom_end - 0.015)

# Clip exons to zoom window
exons_df$start <- pmax(exons_df$start, zoom_start)
exons_df$end   <- pmin(exons_df$end, zoom_end)
exons_df <- exons_df[exons_df$start < zoom_end & exons_df$end > zoom_start, ]

# Gene labels: italic name + strand indicator
genes_df$label <- paste0(genes_df$name, " (", ifelse(genes_df$strand == 1, "+", "\u2212"), ")")

track_genes <- ggplot() +
  annotate("rect", xmin = zoom_start, xmax = zoom_end,
           ymin = 0, ymax = 1, fill = "gray95") +
  # Intron lines (thin)
  geom_segment(data = genes_df,
               aes(x = start_mb, xend = end_mb, y = y, yend = y),
               linewidth = 0.5, color = "gray30") +
  # Exon blocks (thick rectangles)
  geom_rect(data = exons_df,
            aes(xmin = start, xmax = end,
                ymin = y - exon_h, ymax = y + exon_h),
            fill = "gray30", color = NA) +
  # Gene labels with strand
  geom_text(data = genes_df,
            aes(x = label_x, y = y + 0.22, label = label),
            size = 5, fontface = "bold.italic", color = "gray20") +
  vline_layer +
  scale_x_continuous(limits = c(zoom_start, zoom_end)) +
  scale_y_continuous(breaks = NULL, limits = c(0, 1)) +
  labs(y = "Genes", x = NULL) +
  theme_classic(base_size = 20) +
  theme(
    axis.title.y = element_text(face = "bold", size = 18, color = "gray30",
                                angle = 0, hjust = 1, vjust = 0.5),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    plot.margin  = margin(1, 10, 1, 5)
  )

# =============================================================================
# Helper: load BED + build annotation track
# =============================================================================

load_peak_bed <- function(filepath, min_width_bp = 500) {
  if (!file.exists(filepath) || file.info(filepath)$size == 0) {
    return(data.frame(chr = character(), start_mb = numeric(), end_mb = numeric()))
  }
  d <- read.delim(filepath, header = FALSE, stringsAsFactors = FALSE)
  df <- data.frame(
    chr      = d[, 1],
    start_mb = d[, 2] / 1e6,
    end_mb   = d[, 3] / 1e6
  )
  width_mb <- df$end_mb - df$start_mb
  min_width_mb <- min_width_bp / 1e6
  too_narrow <- width_mb < min_width_mb
  if (any(too_narrow)) {
    midpoints <- (df$start_mb[too_narrow] + df$end_mb[too_narrow]) / 2
    df$start_mb[too_narrow] <- midpoints - min_width_mb / 2
    df$end_mb[too_narrow]   <- midpoints + min_width_mb / 2
  }
  return(df)
}

make_annot_track <- function(peak_df, fill_color, track_label,
                             show_x_axis = FALSE) {
  p <- ggplot() +
    annotate("rect", xmin = zoom_start, xmax = zoom_end,
             ymin = 0, ymax = 1, fill = "gray95") +
    scale_x_continuous(limits = c(zoom_start, zoom_end)) +
    scale_y_continuous(breaks = NULL, limits = c(0, 1)) +
    labs(y = track_label) +
    theme_classic(base_size = 20) +
    theme(
      axis.title.y = element_text(face = "bold", size = 18, color = fill_color,
                                  angle = 0, hjust = 1, vjust = 0.5),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      plot.margin  = margin(1, 10, 1, 5)
    )

  if (nrow(peak_df) > 0) {
    p <- p + geom_rect(data = peak_df,
                       aes(xmin = start_mb, xmax = end_mb, ymin = 0, ymax = 1),
                       fill = fill_color, color = NA, alpha = 0.85)
  }

  # Add vertical reference lines on top of annotations
  p <- p + vline_layer

  if (show_x_axis) {
    p <- p + labs(x = "Position on chr1 (Mb)") +
      theme(
        axis.title.x = element_text(face = "bold", size = 19),
        axis.text.x  = element_text(color = "black", size = 17)
      )
  } else {
    p <- p + labs(x = NULL) +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank()
      )
  }

  return(p)
}

# =============================================================================
# Section header helper
# =============================================================================

make_section_header <- function(label_text) {
  ggplot() +
    annotate("rect", xmin = zoom_start, xmax = zoom_end,
             ymin = 0, ymax = 1, fill = "white") +
    annotate("text", x = (zoom_start + zoom_end) / 2, y = 0.5,
             label = label_text,
             size = 6, fontface = "bold.italic", color = "gray40") +
    scale_x_continuous(limits = c(zoom_start, zoom_end)) +
    scale_y_continuous(breaks = NULL, limits = c(0, 1)) +
    labs(y = NULL, x = NULL) +
    theme_void() +
    theme(plot.margin = margin(1, 10, 0, 5))
}

# =============================================================================
# BULK DLPFC — 1 key track
# =============================================================================

cat("Building bulk DLPFC tracks...\n")

peaks_h3k27ac    <- load_peak_bed(file.path(data_dir, "H3K27ac_DLPFC.bed"))

cat(sprintf("  Bulk H3K27ac:     %d peaks\n", nrow(peaks_h3k27ac)))

track_h3k27ac    <- make_annot_track(peaks_h3k27ac,    "#FF9800", "H3K27ac")

# =============================================================================
# EXCITATORY NEURONS — 2 key tracks
# =============================================================================

cat("Building excitatory neuron tracks...\n")

ct_exc_h3k27ac  <- load_peak_bed(file.path(ct_dir, "ExcNeuron_H3K27ac.bed"))
ct_exc_h3k27me3 <- load_peak_bed(file.path(ct_dir, "ExcNeuron_H3K27me3.bed"))

cat(sprintf("  Exc. neuron H3K27ac:  %d peaks\n", nrow(ct_exc_h3k27ac)))
cat(sprintf("  Exc. neuron H3K27me3: %d peaks\n", nrow(ct_exc_h3k27me3)))

track_ct_exc_h3k27ac  <- make_annot_track(ct_exc_h3k27ac,  "#FF9800", "H3K27ac")
track_ct_exc_h3k27me3 <- make_annot_track(ct_exc_h3k27me3, "#7B1FA2", "H3K27me3")

# =============================================================================
# ASTROCYTES — 3 key tracks
# =============================================================================

cat("Building astrocyte tracks...\n")

# Merge DNase narrowpeak replicates
ct_astro_dnase_files <- list.files(ct_dir, pattern = "Astrocyte_DNase_narrowpeak_rep",
                                    full.names = TRUE)
ct_astro_dnase_list <- lapply(ct_astro_dnase_files, load_peak_bed)
ct_astro_dnase <- do.call(rbind, ct_astro_dnase_list)
if (nrow(ct_astro_dnase) > 0) {
  ct_astro_dnase <- ct_astro_dnase[!duplicated(ct_astro_dnase), ]
}

ct_astro_h3k4me1 <- load_peak_bed(file.path(ct_dir, "Astrocyte_H3K4me1.bed"))
ct_astro_h3k4me3 <- load_peak_bed(file.path(ct_dir, "Astrocyte_H3K4me3.bed"))

cat(sprintf("  Astrocyte DNase:   %d peaks (merged)\n", nrow(ct_astro_dnase)))
cat(sprintf("  Astrocyte H3K4me1: %d peaks\n", nrow(ct_astro_h3k4me1)))
cat(sprintf("  Astrocyte H3K4me3: %d peaks\n", nrow(ct_astro_h3k4me3)))

track_ct_astro_dnase   <- make_annot_track(ct_astro_dnase,   "#4CAF50", "DNase")
track_ct_astro_h3k4me1 <- make_annot_track(ct_astro_h3k4me1, "#FFC107", "H3K4me1")
track_ct_astro_h3k4me3 <- make_annot_track(ct_astro_h3k4me3, "#F44336", "H3K4me3",
                                            show_x_axis = TRUE)

# =============================================================================
# Section headers
# =============================================================================

header_bulk  <- make_section_header("Bulk DLPFC (E073)")
header_exc   <- make_section_header("Excitatory Neurons (ENCODE Mint-ChIP)")
header_astro <- make_section_header("Astrocytes (ENCODE Mint-ChIP)")

# =============================================================================
# Combine all tracks
# =============================================================================

cat("Combining tracks...\n")

panel_E <- plot_grid(
  # --- PIP tracks (3) ---
  track_susie, track_ash, track_inf,
  # --- Gene track (1) ---
  track_genes,
  # --- Bulk DLPFC (1) ---
  header_bulk,
  track_h3k27ac,
  # --- Excitatory Neurons (2) ---
  header_exc,
  track_ct_exc_h3k27ac, track_ct_exc_h3k27me3,
  # --- Astrocytes (3) ---
  header_astro,
  track_ct_astro_dnase, track_ct_astro_h3k4me1, track_ct_astro_h3k4me3,
  ncol = 1, align = "v", axis = "lr",
  rel_heights = c(
    # PIP tracks (3)
    0.7, 0.7, 0.7,
    # Gene track
    0.4,
    # Bulk header + 1 track
    0.12,
    0.18,
    # Excitatory header + 2 tracks
    0.12,
    0.18, 0.18,
    # Astrocyte header + 3 tracks (last has x-axis)
    0.12,
    0.18, 0.18, 0.45
  )
)

# =============================================================================
# Save
# =============================================================================

output_pdf <- file.path(script_dir, "panel_E.pdf")
ggsave(output_pdf, panel_E, width = 6, height = 13, units = "in", bg = "white")
cat(sprintf("Saved: %s\n", output_pdf))

output_png <- file.path(script_dir, "panel_E.png")
ggsave(output_png, panel_E, width = 6, height = 13, units = "in",
       dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", output_png))

plot_path <- file.path(script_dir, "panel_E_plot.rds")
saveRDS(panel_E, plot_path)
cat(sprintf("Saved: %s\n", plot_path))

cat("\nDone!\n")
