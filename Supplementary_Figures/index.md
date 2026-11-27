# Supplementary Figures

Scripts that build Supplementary Figures S1-S11. Folder numbers are the manuscript's
figure numbers, and each `S<n>/S<n>.pdf` is the file used in the manuscript.

| Figure | Shows | Drawn by | Saved data written by |
|---|---|---|---|
| S1 | Fold speedup of the SuSiE 2.0 summary-statistics backend for the regularized RSS model, with an LD matrix or a genotype matrix, at p = 5,000 | `plot_S1.R` | `data_generation/` (cluster benchmark) |
| S2 | Sparse simulations, K = 1–5 causal variants: 95% CS power and FDR, variant-level ROC curves, PIP calibration | `plot_S2.R` | `extract_finemapping_data.R`, `extract_pip_data.R` |
| S3 | Credible set size distributions in four simulation scenarios, and median purity | `plot_S3.R` | `extract_S3_data.R` |
| S4 | Gene expression prediction in simulations: cross-validated R² and how often each method is best | `plot_S4.R` | `extract_prediction_data.R` |
| S5 | Illustrative locus (susieR's `unmappable_data` example): PIPs and credible sets of the three methods | `plot_S5.R` | `fit_models.R` |
| S6 | Complex simulations with all three methods: the version of Figure 1C–F that includes SuSiE-ash | `build_S6_panels.R`, `plot_S6.R` | `Main_Figures/figure_1/extract_complex_ldext.R` |
| S7 | Credible set construction: SuSiE and SuSiE+BLiP at purity 0.5 and 0.8, and attainable coverage (A–D); posterior mass of the sets reported by SparsePro and susieR (E) | `plot_S7.R` | `extract_S7_data.R`, `simulate_S7E.R` |
| S8 | Correlation-based credible set extension in the sparse simulations | `plot_S8.R` | `extract_S8_data.R` |
| S9 | Correlation-based credible set extension in the complex simulations | `plot_S9.R` | `extract_S9_data.R` |
| S10 | ROSMAP analysis with all three methods: the version of Figure 2 that includes SuSiE-ash | `plot_S10.R` | the Figure 2 panel scripts (see below) |
| S11 | Mean runtime per method | `plot_S11.R` | `extract_S11_data.R` |

`Rscript run_all.R` runs every script in the "Drawn by" column. The scripts in the
last column need inputs from outside this repository; the root `README.md` lists
what each needs.

## Credible set conventions

Every simulation figure reports 95% credible sets with purity of at least 0.5,
computed from the in-sample correlation matrix, and no correlation extension
(`xcorr_noext` in the `extract_*` scripts). This matches the susieR default and the
real-data analysis. Figures S8 and S9 add the correlation extension at 0.99 for
comparison. A replicate that reports no credible set contributes power 0 and FDR 0.
In the complex scenarios the N variants with the largest true effects are treated
as causal, with N capped at the number of nonzero effects.

## Panel E of Figure S7

Panel E is drawn by `S7/plot_S7.R` from the saved table
`S7/data/paired-posterior-mass.csv`. To regenerate that table:

```bash
Rscript Supplementary_Figures/S7/simulate_S7E.R
```

The simulation uses the public 121-variant ADSP sketch example in `susieR`, with
seed 83011: 150 replicates for each combination of one or two causal variants
and total signal strength `nh²` of 24 or 80. Both reporting rules receive the
same fitted SuSiE posteriors. The script downloads the released 2023 SparsePro
source (commit `e13ac87d6e650d7682ff1f33f6f9898d66932d64`) and calls its original
`get_effect()` method; the ADSP example and the `susieR` reporting code come from
commit `ed6e491e3d0235ea5fa2e8545ea206b5d759eb6e`. It needs Python 3 with numpy,
pandas and scipy, and the pinned `susieR` 0.16.6 fitting build:

```r
remotes::install_github(
  "stephenslab/susieR@ef213feed2cb82419677661a8c986e1504df2c73",
  upgrade = "never"
)
```

The expected result is 599 converged fits and 567 paired sets, with 152
SparsePro sets and 0 `susieR` sets below 95% of the full single-effect posterior
mass.

## The three-method panels of Figure S10

The panels of Figure S10 are drawn by the Figure 2 panel scripts run in
three-method mode, which write into `S10/data/`. Panel B is the panel of
Figure 2 itself, because SuSiE-ash is not among the prediction methods it
compares.

```bash
# cross-method matching of the credible sets of all three methods
Rscript Main_Figures/figure_2/data/data_generation/alphagenome_scoring/alphagenome_cs_group_comparison.R \
  standard,ash,inf
# panel A
Rscript Main_Figures/figure_2/figure_2_panel_A/prepare_panel_A_data.R --include-ash
Rscript Main_Figures/figure_2/figure_2_panel_A/plot_panel_A.R --include-ash
# panel C
Rscript Main_Figures/figure_2/figure_2_panel_C/prepare_panel_C_data.R \
  Main_Figures/figure_2/data/alphagenome_cs_group_assignments.rds \
  Supplementary_Figures/S10/data/panel_C_data.rds
Rscript Main_Figures/figure_2/figure_2_panel_C/plot_panel_C.R \
  Supplementary_Figures/S10/data/panel_C_data.rds Supplementary_Figures/S10/data/panel_C
# panel D
Rscript Main_Figures/figure_2/figure_2_panel_D/prepare_panel_D_data.R \
  Main_Figures/figure_2/data/alphagenome_cs_group_scores.csv \
  Main_Figures/figure_2/data/alphagenome_cs_group_assignments.rds \
  Supplementary_Figures/S10/data/panel_D_data.rds
Rscript Main_Figures/figure_2/figure_2_panel_D/plot_panel_D.R \
  Supplementary_Figures/S10/data/panel_D_data.rds Supplementary_Figures/S10/data
# panel E
Rscript Main_Figures/figure_2/figure_2_panel_E/prepare_panel_E_data.R --include-ash
Rscript Main_Figures/figure_2/figure_2_panel_E/plot_panel_E.R --include-ash
# the figure
Rscript Supplementary_Figures/S10/plot_S10.R
```

`run_all.R` runs these steps, after the Figure 2 steps that draw panel B.
