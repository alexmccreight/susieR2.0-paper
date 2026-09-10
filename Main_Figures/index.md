# Main Figures

This section contains the scripts used to assemble each main figure from
pre-computed results. The upstream simulation and analysis code lives under
**Simulation Studies** and **Data Applications**.

Both figures show **SuSiE** and **SuSiE-inf**. The three-method comparisons that
additionally include **SuSiE-ash** are preserved as Supplementary Figures
**S10** (Figure 1C–1F) and **S11** (Figure 2).

## Figure 1

SuSiE 2.0 software architecture overview, the runtime speedup of its new
SuSiE-inf implementation over the original, and simulation benchmarks comparing
fine-mapping power and FDR across two complex (oligogenic-plus-infinitesimal)
simulation scenarios.

![Figure 1](figure_1/susieR2_figure1.png)

Built by `figure_1/create_figure1.R` from the per-panel objects in
`figure_1/figure_1_panel_{A..F}/`. Panels C–F are drawn by the shared function
in `R/topn_panels.R` from `figure_1/data/complex_S{1,2}_finemapping_summary.rds`
(95% credible sets, exact purity ≥ 0.5, no correlation extension; a replicate
with no credible set contributes FDR = 0). Supplementary Figure S10 draws the
same tables with SuSiE-ash included.

## Figure 2

ROSMAP eQTL fine-mapping application of SuSiE 2.0 across three brain tissues,
covering credible set discovery, TWAS prediction, cross-method signal
concordance, AlphaGenome functional enrichment, and a locus example.

![Figure 2](figure_2/susieR2_figure2.png)

Built by `figure_2/create_figure2.R`. Panel D tests one observation per matched
signal and modality (member credible sets of a consensus signal are averaged
before the t-test). Panels are laid out and labelled in reading
order **A, B, C, D** across the top row with **E** spanning the second row; note
that the source directories keep their original names, so `figure_2_panel_C`
supplies the panel labelled **B** and `figure_2_panel_B` supplies the panel
labelled **C**.

## Reproducing

From the repository root:

```bash
Rscript run_all.R
```

This rebuilds every main and supplementary figure from data vendored in this
repository. See the root `README.md` for the path conventions and the two
environment variables that only the `extract_*` scripts need.
