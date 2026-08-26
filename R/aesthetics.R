# =============================================================================
# Shared aesthetics for susieR2.0-paper figures
# =============================================================================
# Single source of truth for the method palette and panel theme. These used to
# be copy-pasted into every panel script, which meant a palette change required
# editing five or more files.
#
# Usage:
#     source("R/aesthetics.R")
#     scale_color_manual(values = method_colors, breaks = MAIN_METHODS)
# =============================================================================

# --- Method identity -----------------------------------------------------
# Methods shown in the MAIN figures. SuSiE-ash was moved to the supplement,
# so it is deliberately absent here.
MAIN_METHODS <- c("SuSiE", "SuSiE-inf")

# Methods shown in the SUPPLEMENTARY figures (S10, S11), which retain the
# full three-method comparison.
ALL_METHODS <- c("SuSiE", "SuSiE-ash", "SuSiE-inf")

# Internal method keys as they appear in the real-data pipeline outputs.
METHOD_KEYS <- c(SuSiE = "standard", `SuSiE-ash` = "ash", `SuSiE-inf` = "inf")

method_colors <- c(
  "SuSiE"     = "#4A90E2",  # blue
  "SuSiE-inf" = "#7CB342",  # green
  "SuSiE-ash" = "#E53935"   # red (supplement only)
)

# --- Prediction methods (Figure 2 panel C) --------------------------------
# SuSiE and SuSiE-inf are taken from method_colors above so this panel can never
# drift from the rest of the figure. BayesR is a deep purple rather than the
# teal it used to be: teal (#2E9D8F) is the "SuSiE&inf" concordance colour and
# read as another green next to SuSiE-inf. #6A1B9A is absent from every other
# palette here, so it collides with nothing.
prediction_colors <- c(
  "Elastic Net" = "#FF9800",
  "LASSO"       = "#D81B60",
  "BayesR"      = "#6A1B9A",
  "BayesL"      = "#78909C"
)

# --- Concordance groups (Figure 2 panels B and D) ------------------------
# Three-method vocabulary, retained for the supplement.
concordance_colors <- c(
  "Consensus"   = "#888888",
  "SuSiE&ash"   = "#9C27B0",  # blue + red   -> purple
  "SuSiE&inf"   = "#2E9D8F",  # blue + green -> teal
  "ash&inf"     = "#FF9800",  # red + green  -> orange
  "SuSiE"       = "#4A90E2",
  "SuSiE-ash"   = "#E53935",
  "SuSiE-inf"   = "#7CB342"
)

# Canonical mapping from the internal group labels produced by
# alphagenome_cs_group_comparison.R to the display labels used in panels B and
# D. The ORDER of this vector is the canonical plotting order:
# Consensus -> pairwise -> single-method.
GROUP_DISPLAY <- c(
  "Consensus"         = "Consensus",
  "Standard+ASH"      = "SuSiE&ash",
  "Standard+Inf"      = "SuSiE&inf",
  "ASH+Inf"           = "ash&inf",
  "Standard-specific" = "SuSiE",
  "ASH-specific"      = "SuSiE-ash",
  "Inf-specific"      = "SuSiE-inf"
)

# --- Shared panel theme --------------------------------------------------
theme_panel <- ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    axis.title  = ggplot2::element_text(face = "bold", size = 10),
    axis.text   = ggplot2::element_text(color = "black", size = 9),
    legend.position = "none",
    plot.margin = ggplot2::margin(5, 5, 5, 5)
  )
