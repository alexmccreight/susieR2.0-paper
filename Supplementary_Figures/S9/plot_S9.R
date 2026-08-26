#!/usr/bin/env Rscript

# =============================================================================
# Supplemental S9: Runtime comparison -- SuSiE vs SuSiE-ash vs SuSiE-inf
# =============================================================================
# One bar per method (mean per-method wall-clock fit time over 500 replicates)
# with mean +/- SE error bars. Same method colors as S1/S5/S7/S8:
#   SuSiE = blue (#4A90E2), SuSiE-ash = red (#E53935), SuSiE-inf = green (#7CB342)
#
# Input:  data/s9_runtime_data.rds  (built by extract_S9_data.R)
# Output: S9.pdf, S9.png
# =============================================================================

suppressMessages(library(ggplot2))

source("R/paths.R")

s9_dir    <- supp_dir(9)
data_path <- file.path(s9_dir, "data", "s9_runtime_data.rds")
if (!file.exists(data_path)) stop("Run extract_S9_data.R first.")

dat <- readRDS(data_path)

# =============================================================================
# Aesthetics (shared method palette)
# =============================================================================

method_levels <- c("SuSiE", "SuSiE-ash", "SuSiE-inf")
dat$Method <- factor(dat$Method, levels = method_levels)

method_colors <- c(
  "SuSiE"     = "#4A90E2",
  "SuSiE-ash" = "#E53935",
  "SuSiE-inf" = "#7CB342"
)

theme_bar <- theme_classic(base_size = 14) +
  theme(
    axis.title      = element_text(face = "bold", size = 18),
    axis.text       = element_text(color = "black", size = 15),
    axis.text.x     = element_text(face = "bold", size = 16),
    legend.position = "none",
    plot.margin     = margin(12, 14, 12, 12)
  )

# =============================================================================
# Bar panel: mean runtime +/- SE
# =============================================================================

y_top <- max(dat$mean_time + dat$se_time) * 1.12

p <- ggplot(dat, aes(x = Method, y = mean_time, fill = Method)) +
  geom_bar(stat = "identity", width = 0.65) +
  geom_errorbar(
    aes(ymin = mean_time - se_time, ymax = mean_time + se_time),
    width = 0.18, linewidth = 0.5
  ) +
  geom_text(
    aes(label = sprintf("%.1f s", mean_time)),
    vjust = -0.6, nudge_y = max(dat$se_time) * 1.1,
    fontface = "bold", size = 5
  ) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, y_top)) +
  labs(x = NULL, y = "Runtime (seconds)") +
  theme_bar

out_pdf <- file.path(s9_dir, "S9.pdf")
ggsave(out_pdf, p, width = 7, height = 6, units = "in", bg = "white")
cat(sprintf("Saved: %s\n", out_pdf))

out_png <- file.path(s9_dir, "S9.png")
ggsave(out_png, p, width = 7, height = 6, units = "in", dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", out_png))

cat("\nDone!\n")
