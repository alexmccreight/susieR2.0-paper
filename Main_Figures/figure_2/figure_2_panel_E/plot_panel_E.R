#!/usr/bin/env Rscript

# Figure 2, panel E: the LMOD1 locus (ENSG00000163431) in DLPFC, zoomed to
# chr1:201.78-201.97 Mb. One PIP track per method, then a gene track and six
# epigenomic peak tracks (bulk DLPFC, excitatory neurons, astrocytes).
#
# Credible-set variants are circled in the colour of their signal. The main
# figure has three signals, none shared by the two methods:
#   SuSiE CS1      chr1:201915540:A:G           blue (the SuSiE colour)
#   SuSiE CS2      six variants across ~150 kb  magenta
#   SuSiE-inf CS1  chr1:201917641:TA:T          green (the SuSiE-inf colour)
#
# Usage:  plot_panel_E.R [--include-ash]
#         --include-ash draws the three-method panel of Supplementary Figure
#         S10 from Supplementary_Figures/S10/data/panel_E_data.rds and writes
#         only panel_E_plot.rds, in that folder.
# Input:  panel_E_data.rds  (prepare_panel_E_data.R)
#         gene_models.rds   (prepare_gene_models.R)
#         data/panel_E/H3K27ac_DLPFC.bed, data/panel_E/celltype_tracks/*.bed
# Output: panel_E.pdf, panel_E.png, panel_E_plot.rds

library(ggplot2)
library(cowplot)

source("R/paths.R")
source("R/aesthetics.R")

args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% "--include-ash")) {
  stop("Unknown argument(s): ", paste(args[!args %in% "--include-ash"], collapse = ", "))
}
include_ash <- "--include-ash" %in% args

fig2_dir   <- fig_dir(2)
script_dir <- file.path(fig2_dir, "figure_2_panel_E")
data_dir   <- file.path(fig2_dir, "data", "panel_E")
ct_dir     <- file.path(data_dir, "celltype_tracks")
output_dir <- if (include_ash) file.path(supp_dir(10), "data") else script_dir

zoom_start <- 201.78   # Mb
zoom_end   <- 201.97   # Mb

# Lead variants of the credible sets.
susie_cs1_pos_mb <- 201915540 / 1e6
susie_cs2_pos_mb <- 201836798 / 1e6
inf_cs1_pos_mb   <- 201917641 / 1e6   # also the lead of SuSiE-ash CS2
ash_lead_pos_mb  <- 201806618 / 1e6   # SuSiE-ash CS1

# ---- Signals and their colours ----

# SuSiE CS1 and SuSiE-inf CS1 lie 2.1 kb apart but are different signals: the
# two variants correlate at r = -0.60 in the DLPFC genotypes, and the
# concordance analysis behind panels C and D classifies one as SuSiE-specific
# and the other as SuSiE-inf-specific. Each therefore takes its method colour.
# SuSiE's diffuse CS2 takes magenta, a hue that none of the peak tracks uses.
signal_map <- data.frame(
  method = c("SuSiE",  "SuSiE",  "SuSiE-inf"),
  cs_id  = c(1,        2,        1),
  signal = c("A",      "B",      "C")
)

signal_colors <- c(
  "A" = method_colors[["SuSiE"]],
  "B" = "#C2185B",
  "C" = method_colors[["SuSiE-inf"]]
)

if (include_ash) {
  # SuSiE-ash CS2 has the same lead variant as the SuSiE-inf singleton and
  # shares its green, although it holds three variants. SuSiE-ash CS1 is orange.
  signal_map <- rbind(signal_map, data.frame(
    method = c("SuSiE-ash", "SuSiE-ash"),
    cs_id = c(1, 2), signal = c("D", "C")
  ))
  signal_colors <- c(signal_colors, "D" = "#FF9800")
}

# ---- PIP data ----

pip_data <- readRDS(file.path(output_dir, "panel_E_data.rds"))$pip_data
pip_data$pos_mb <- pip_data$position / 1e6

# The fitted region spans about 5.5 Mb; only the zoom window is drawn.
pip_data <- pip_data[pip_data$pos_mb >= zoom_start & pip_data$pos_mb <= zoom_end, ]

pip_data$signal <- signal_map$signal[match(
  paste(pip_data$method, pip_data$cs_id),
  paste(signal_map$method, signal_map$cs_id)
)]

# ---- Elements shared by the tracks ----

# Vertical lines at the lead variants: SuSiE-inf CS1 (green), SuSiE CS1 (blue)
# and, with --include-ash, SuSiE-ash CS1 (orange).
lead_lines <- function(linewidth, susie_linetype) {
  c(
    if (include_ash) {
      list(geom_vline(xintercept = ash_lead_pos_mb, linetype = "dashed",
                      color = signal_colors[["D"]], linewidth = linewidth, alpha = 0.7))
    },
    list(
      geom_vline(xintercept = inf_cs1_pos_mb, linetype = "dashed",
                 color = signal_colors[["C"]], linewidth = linewidth, alpha = 0.7),
      geom_vline(xintercept = susie_cs1_pos_mb, linetype = susie_linetype,
                 color = signal_colors[["A"]], linewidth = linewidth, alpha = 0.8)
    )
  )
}

# Only the bottom track shows the x-axis.
no_x_axis <- theme(
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank(),
  axis.line.x  = element_blank()
)

# ---- PIP tracks ----

make_pip_track <- function(method_name) {
  d <- pip_data[pip_data$method == method_name, ]
  d_bg <- d[!d$in_cs, ]
  d_cs <- d[d$in_cs, ]

  ggplot() +
    lead_lines(0.8, "dotdash") +
    geom_point(data = d_bg, aes(x = pos_mb, y = pip),
               size = 1.2, color = "gray70", alpha = 0.5) +
    geom_point(data = d_cs, aes(x = pos_mb, y = pip, color = signal),
               shape = 21, fill = NA, size = 8, stroke = 2.0,
               alpha = 0.9) +
    geom_point(data = d_cs, aes(x = pos_mb, y = pip),
               size = 4, color = "black", alpha = 0.9) +
    scale_color_manual(values = signal_colors, guide = "none") +
    # The upper limit of 1.08 keeps a circle at PIP = 1 inside the panel.
    scale_y_continuous(limits = c(0, 1.08), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_continuous(limits = c(zoom_start, zoom_end)) +
    labs(x = NULL, y = "PIP", title = method_name) +
    theme_classic(base_size = 20) +
    theme(
      axis.title.y = element_text(face = "bold", size = 19),
      axis.text.y  = element_text(color = "black", size = 17),
      plot.title   = element_text(face = "bold", size = 20, hjust = 0),
      plot.margin  = margin(2, 10, 2, 5)
    ) +
    no_x_axis
}

# A variant callout: a box with the variant and its max |Q_RNA| (AlphaGenome
# RNA-seq quantile score, entered by hand), and a leader line from the box
# anchor (box_x, box_y) to the variant. The anchor is the top centre of the box
# (the bottom centre with vjust = 0).
callout <- function(variant, score, box_x, box_y, to_x, to_y,
                    size = 8.5, padding = 0.22, vjust = 1) {
  list(
    annotate("segment", x = box_x, xend = to_x, y = box_y, yend = to_y,
             linewidth = 0.4, color = "gray30"),
    annotate("label", x = box_x, y = box_y,
             label = sprintf('atop("%s", max~bgroup("|", Q[RNA], "|") == %s)',
                             variant, score),
             parse = TRUE, hjust = 0.5, vjust = vjust,
             size = size, fontface = "bold", fill = "white",
             label.r = unit(0.1, "lines"), label.padding = unit(padding, "lines"),
             color = "gray20")
  )
}

if (include_ash) {
  # Three shorter tracks with smaller boxes.
  pip_tracks <- list(
    make_pip_track("SuSiE") +
      callout("chr1:201915540:A:G", "0.125", 201.88, 0.80,
              susie_cs1_pos_mb + 0.001, 0.996, size = 4.2, padding = 0.15),
    make_pip_track("SuSiE-ash") +
      callout("chr1:201806618:G:A", "0.910", 201.85, 0.40,
              ash_lead_pos_mb + 0.003, 0.37, size = 5, padding = 0.15, vjust = 0),
    make_pip_track("SuSiE-inf") +
      callout("chr1:201917641:TA:T", "0.999", 201.88, 0.80,
              inf_cs1_pos_mb - 0.003, 1.0, size = 4.2, padding = 0.15)
  )
  pip_heights <- c(0.7, 0.7, 0.7)
} else {
  # Two taller tracks. The PIP axis ends just above 1, so the boxes sit inside
  # the data area, where they cover no credible-set variant.
  pip_tracks <- list(
    make_pip_track("SuSiE") +
      callout("chr1:201915540:A:G", "0.125", 201.9445, 0.82,
              susie_cs1_pos_mb + 0.002, 0.990) +
      callout("chr1:201836798:T:C", "0.913", 201.8112, 1.00,
              susie_cs2_pos_mb, 0.585),
    make_pip_track("SuSiE-inf") +
      callout("chr1:201917641:TA:T", "0.999", 201.8869, 0.82,
              inf_cs1_pos_mb - 0.002, 0.992)
  )
  pip_heights <- c(1.1, 1.1)
}

# ---- Gene and peak tracks ----

# A grey band across the window with the track name to its left.
track_base <- function(track_label, label_color) {
  ggplot() +
    annotate("rect", xmin = zoom_start, xmax = zoom_end,
             ymin = 0, ymax = 1, fill = "gray95") +
    scale_x_continuous(limits = c(zoom_start, zoom_end)) +
    scale_y_continuous(breaks = NULL, limits = c(0, 1)) +
    labs(x = NULL, y = track_label) +
    theme_classic(base_size = 20) +
    theme(
      axis.title.y = element_text(face = "bold", size = 18, color = label_color,
                                  angle = 0, hjust = 1, vjust = 0.5),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      plot.margin  = margin(1, 10, 1, 5)
    )
}

# Gene track: exon blocks on a line spanning each gene, on the two rows
# assigned by prepare_gene_models.R.
gene_models <- readRDS(file.path(script_dir, "gene_models.rds"))
genes_df    <- gene_models$genes
exons_df    <- gene_models$exons

exon_h <- 0.14  # half-height of exon blocks
genes_df$y <- ifelse(genes_df$row == 1, 0.62, 0.22)
exons_df$y <- genes_df$y[match(exons_df$gene, genes_df$name)]

# Gene spans and exons are clipped to the zoom window.
genes_df$start_mb <- pmax(genes_df$start / 1e6, zoom_start)
genes_df$end_mb   <- pmin(genes_df$end / 1e6, zoom_end)
genes_df$label_x  <- pmin((genes_df$start_mb + genes_df$end_mb) / 2,
                          zoom_end - 0.015)

exons_df$start <- pmax(exons_df$start / 1e6, zoom_start)
exons_df$end   <- pmin(exons_df$end / 1e6, zoom_end)
exons_df <- exons_df[exons_df$start < zoom_end & exons_df$end > zoom_start, ]

# ASCII "-" for the minus strand: the pdf device cannot encode the Unicode
# minus sign (U+2212) and prints a dot in its place.
genes_df$label <- paste0(genes_df$name, " (", ifelse(genes_df$strand == 1, "+", "-"), ")")

track_genes <- track_base("Genes", "gray30") +
  geom_segment(data = genes_df,
               aes(x = start_mb, xend = end_mb, y = y, yend = y),
               linewidth = 0.5, color = "gray30") +
  geom_rect(data = exons_df,
            aes(xmin = start, xmax = end,
                ymin = y - exon_h, ymax = y + exon_h),
            fill = "gray30", color = NA) +
  geom_text(data = genes_df,
            aes(x = label_x, y = y + 0.22, label = label),
            size = 5, fontface = "bold.italic", color = "gray20") +
  lead_lines(0.7, "dotted") +
  no_x_axis

# Peaks of a BED file, in Mb. Peaks narrower than 500 bp are widened to 500 bp
# around their midpoint so that they stay visible.
load_peak_bed <- function(filepath) {
  bed   <- read.delim(filepath, header = FALSE)
  peaks <- data.frame(start_mb = bed[, 2] / 1e6, end_mb = bed[, 3] / 1e6)

  min_width_mb <- 500 / 1e6
  too_narrow   <- peaks$end_mb - peaks$start_mb < min_width_mb
  midpoints    <- (peaks$start_mb[too_narrow] + peaks$end_mb[too_narrow]) / 2
  peaks$start_mb[too_narrow] <- midpoints - min_width_mb / 2
  peaks$end_mb[too_narrow]   <- midpoints + min_width_mb / 2
  peaks
}

# A peak that extends beyond the zoom window is not drawn (the x scale turns
# its out-of-range edge into NA); it is not clipped to the window.
make_annot_track <- function(peak_df, fill_color, track_label,
                             show_x_axis = FALSE) {
  p <- track_base(track_label, fill_color) +
    geom_rect(data = peak_df,
              aes(xmin = start_mb, xmax = end_mb, ymin = 0, ymax = 1),
              fill = fill_color, color = NA, alpha = 0.85) +
    lead_lines(0.7, "dotted")

  if (show_x_axis) {
    p + labs(x = "Position on chr1 (Mb)") +
      theme(
        axis.title.x = element_text(face = "bold", size = 19),
        axis.text.x  = element_text(color = "black", size = 17)
      )
  } else {
    p + no_x_axis
  }
}

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

# The astrocyte DNase track is the union of the three replicates' narrow peaks.
astro_dnase <- unique(do.call(rbind, lapply(
  list.files(ct_dir, pattern = "Astrocyte_DNase_narrowpeak_rep", full.names = TRUE),
  load_peak_bed
)))

lower_tracks <- list(
  track_genes,
  make_section_header("Bulk DLPFC (E073)"),
  make_annot_track(load_peak_bed(file.path(data_dir, "H3K27ac_DLPFC.bed")),
                   "#FF9800", "H3K27ac"),
  make_section_header("Excitatory Neurons (ENCODE Mint-ChIP)"),
  make_annot_track(load_peak_bed(file.path(ct_dir, "ExcNeuron_H3K27ac.bed")),
                   "#FF9800", "H3K27ac"),
  make_annot_track(load_peak_bed(file.path(ct_dir, "ExcNeuron_H3K27me3.bed")),
                   "#7B1FA2", "H3K27me3"),
  make_section_header("Astrocytes (ENCODE Mint-ChIP)"),
  make_annot_track(astro_dnase, "#4CAF50", "DNase"),
  make_annot_track(load_peak_bed(file.path(ct_dir, "Astrocyte_H3K4me1.bed")),
                   "#FFC107", "H3K4me1"),
  make_annot_track(load_peak_bed(file.path(ct_dir, "Astrocyte_H3K4me3.bed")),
                   "#F44336", "H3K4me3", show_x_axis = TRUE)
)
# Genes 0.4, section headers 0.12, peak tracks 0.18; the last track is taller
# because it carries the x-axis.
lower_heights <- c(0.4, 0.12, 0.18, 0.12, 0.18, 0.18, 0.12, 0.18, 0.18, 0.45)

# ---- Assemble ----

panel_E <- plot_grid(
  plotlist = c(pip_tracks, lower_tracks),
  ncol = 1, align = "v", axis = "lr",
  rel_heights = c(pip_heights, lower_heights)
)

# The standalone panel is saved at the size it has in Figure 2 (the full-width
# bottom row, about 18 x 10 in), so that its text has the same relative size.
if (!include_ash) {
  save_figure(panel_E, output_dir, "panel_E", width = 18, height = 10.1)
}

plot_path <- file.path(output_dir, "panel_E_plot.rds")
saveRDS(panel_E, plot_path)
cat("Saved:", plot_path, "\n")
