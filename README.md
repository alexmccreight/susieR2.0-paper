# susieR 2.0 Manuscript Resources

Code and data to reproduce figures in susieR 2.0 manuscript.


## About susieR 2.0

susieR 2.0 is a complete architectural redesign that addresses code duplication and fragmented architecture of the original SuSiE implementation while adding new features including unmappable effects modeling and variuos performance optimizations, all while maintaing full backward compatibility.

susieR 2.0 eliminates this duplicative architecture through a unified framework built on modular design principles. The user-facing interface remains unchanged, but the implementation now uses generic functions with data-type specific backends through R’s S3 dispatch system. This architecture enables the integration of new SuSiE extensions while maintaining identical results to their original versions. Beyond architectural improvements, this release introduces substantial algorithmic advances including support for unmappable effects modeling, enhanced computational speed for regularized LD matrices, new convergence criteria, and improved refinement procedures.

## Repository Structure

This repository contains the code and data used to generate all figures in the susieR 2.0 manuscript: <https://github.com/alexmccreight/susieR2.0-paper>

A rendered version of this resource is published with GitHub Pages at <https://alexmccreight.github.io/susieR2.0-paper/>.

- **Simulation Studies**: This directory contains simulation designs and implementation codes used in the paper. 

- **Data Applications**: The ROSMAP bulk-brain eQTL fine-mapping run behind Figure 2, and the AlphaGenome scoring applied to its credible sets.

- **Main Figures**: Assembly code for Figures 1 and 2, plus the panel-level scripts and vendored inputs each one needs.

- **Supplementary Figures**: Assembly code for Supplementary Figures S1-S11, including the three-method (SuSiE / SuSiE-ash / SuSiE-inf) versions of the main figures.

## Getting Started

To navigate this resource, use the table of contents in the left sidebar. Each
section documents one part of the manuscript:

1. The code used to generate the analyses
2. The data associated with each figure
3. How to reproduce every figure (see **Reproducing the figures** below)

## Reproducing the figures

**All scripts are run from the repository root**, so that `source("R/paths.R")` resolves:

```
Rscript Main_Figures/figure_1/figure_1_panel_C/plot_panel_C.R
Rscript Main_Figures/figure_1/create_figure1.R
```

Or rebuild everything that does not require external data:

```
Rscript run_all.R
```

### Paths

Scripts contain no absolute paths. Two shared helpers replace them:

- `R/paths.R` provides `paper_root()`, `fig_dir(n)`, `supp_dir(n)`, `benchmark_root()`
  and `susier_root()`. The repository root is located by walking up for the
  `.susieR2-paper-root` marker file, so scripts work under `Rscript`, `source()`
  and RStudio.
- `R/aesthetics.R` provides the shared method palette (`method_colors`,
  `MAIN_METHODS`, `ALL_METHODS`, `concordance_colors`, `GROUP_DISPLAY`) and
  `theme_panel`, so a palette change is a one-line edit rather than the same
  block repeated across panel scripts.

Every figure input needed to rebuild the main and supplementary figures is
vendored in this repository. Two categories of script reach outside it:

| Environment variable | Default | Needed by |
|---|---|---|
| `SUSIER2_BENCHMARK_ROOT` | `../susieR2.0-benchmark` | the `extract_*` scripts, which read raw 500-replicate simulation output (~22 GB, far too large to vendor) |
| `SUSIER_ROOT` | `../susieR` | `Supplementary_Figures/S2/fit_models.R`, which fits models with the in-development susieR package |

Set them if your checkouts live elsewhere:

```
export SUSIER2_BENCHMARK_ROOT=/path/to/susieR2.0-benchmark
```

The HPC scripts under `Simulation_Studies/` still contain cluster paths
(`/home/apm2217/...`); they are run on the cluster and are not needed to rebuild
any figure.

### Methods shown

Main Figures 1 and 2 show **SuSiE** and **SuSiE-inf**. The three-method
comparisons that additionally include **SuSiE-ash** are preserved as
Supplementary Figures **S10** (was Figure 1C–1F) and **S11** (was Figure 2).

## Computational Requirements

The analyses in this book were performed using:
- R version 4.1 or higher
- Key R packages: data.table, ggplot2, dplyr

## susieR 2.0 Tutorial Website

Learn how to perform colocalization analysis with step-by-step examples. For detailed tutorials and use cases in [Tutorials](https://statfungen.github.io/susieR/articles/index.html).

## Citation

If you use susieR 2.0 in your research, please cite:
