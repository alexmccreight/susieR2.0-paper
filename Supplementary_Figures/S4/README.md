# Supplementary Figure S7 in the manuscript

The repository calls this figure S4. Panels A–D use the saved results in
`data/s4_metrics.rds`; panel E uses `data/paired-posterior-mass.csv`.

From the repository root, redraw the complete figure with

```sh
Rscript --vanilla Supplementary_Figures/S4/plot_S4.R
```

This writes `S4.pdf` and `S4.png` here. R packages `ggplot2`, `cowplot`,
and `dplyr` are required. Drawing does not refit models. To rebuild the
A–D input from the original benchmark files, see `extract_S4_data.R`.

## Reproduce the panel E comparison

The comparison uses the public 121-variant ADSP sketch example in `susieR`,
with seed 83011. It runs 150 replicates for each combination of one or two
causal variants and total signal strength `nh²` of 24 or 80. Both reporting
rules receive the same fitted SuSiE posteriors.

Use R with the pinned fitting version of `susieR`. Install it in a separate
R library if you already use another version.

```r
remotes::install_github(
  "stephenslab/susieR@ef213feed2cb82419677661a8c986e1504df2c73",
  upgrade = "never"
)
```

Python 3 with `numpy`, `pandas`, and `scipy`, plus `curl` and `shasum`, are
also required. Then run from the repository root

```sh
sh Supplementary_Figures/S4/reproducibility/reproduce.sh /tmp/susie2-S7
```

The script downloads the original 2023 SparsePro `get_effect()` implementation
at commit `e13ac87d6e650d7682ff1f33f6f9898d66932d64`. The ADSP example and
`susieR` reporting source come from
`ed6e491e3d0235ea5fa2e8545ea206b5d759eb6e`. Downloads are checked by SHA-256.

The regenerated table is checked against the saved table before drawing
all five panels together. The output directory receives `S4.pdf`, `S4.png`,
and `paired-posterior-mass.csv`. The expected result is 599 converged fits
and 567 paired sets. SparsePro reports 152 sets below 95% of the full
single effect posterior mass; all corresponding `susieR` sets reach 95%.
