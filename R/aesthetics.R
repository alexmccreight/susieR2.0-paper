# Method names, palettes and the panel theme shared by the figure scripts.
# Usage: source("R/aesthetics.R"), then e.g.
#     scale_color_manual(values = method_colors, breaks = MAIN_METHODS)

# ---- Methods ----
# The main figures show SuSiE and SuSiE-inf; the supplement keeps all three.
MAIN_METHODS <- c("SuSiE", "SuSiE-inf")
ALL_METHODS  <- c("SuSiE", "SuSiE-ash", "SuSiE-inf")

# Internal method keys as they appear in the real-data pipeline outputs.
METHOD_KEYS <- c(SuSiE = "standard", `SuSiE-ash` = "ash", `SuSiE-inf` = "inf")

method_colors <- c(
  "SuSiE"     = "#4A90E2",  # blue
  "SuSiE-inf" = "#7CB342",  # green
  "SuSiE-ash" = "#E53935"   # red (supplement only)
)

# ---- Other prediction methods (Figure 2B) ----
# SuSiE and SuSiE-inf take their colours from method_colors.
prediction_colors <- c(
  "Elastic Net" = "#FF9800",
  "LASSO"       = "#D81B60",
  "BayesR"      = "#6A1B9A",
  "BayesL"      = "#78909C"
)

# ---- Concordance groups (Figure 2C, 2D and S10) ----
concordance_colors <- c(
  "Consensus"   = "#888888",
  "SuSiE&ash"   = "#9C27B0",  # blue + red   -> purple
  "SuSiE&inf"   = "#2E9D8F",  # blue + green -> teal
  "ash&inf"     = "#FF9800",  # red + green  -> orange
  "SuSiE"       = "#4A90E2",
  "SuSiE-ash"   = "#E53935",
  "SuSiE-inf"   = "#7CB342"
)

# Group labels written by alphagenome_cs_group_comparison.R -> display labels.
# The order of this vector is the plotting order.
GROUP_DISPLAY <- c(
  "Consensus"         = "Consensus",
  "Standard+ASH"      = "SuSiE&ash",
  "Standard+Inf"      = "SuSiE&inf",
  "ASH+Inf"           = "ash&inf",
  "Standard-specific" = "SuSiE",
  "ASH-specific"      = "SuSiE-ash",
  "Inf-specific"      = "SuSiE-inf"
)

# ---- Panel theme ----
theme_panel <- ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    axis.title  = ggplot2::element_text(face = "bold", size = 10),
    axis.text   = ggplot2::element_text(color = "black", size = 9),
    legend.position = "none",
    plot.margin = ggplot2::margin(5, 5, 5, 5)
  )
