# SuSiE 2.0 manuscript: figure code and data

Code and data behind the figures of

> McCreight, A., Cho, Y., Li, R., Carbonetto, P., Nachun, D., Gan, H-Y., Kunkel, D.,
> Stephens, M., Denault, W.R.P. & Wang, G. (2025). SuSiE 2.0:
> improved methods and implementations for genetic fine-mapping and
> phenotype prediction. *bioRxiv*.
> https://doi.org/10.1101/2025.11.25.690514

The software itself is the [`susieR`](https://github.com/stephenslab/susieR) package
([documentation](https://stephenslab.github.io/susieR)). A rendered version of this repository is at
<https://alexmccreight.github.io/susieR2.0-paper/>.

## What is here

| Folder | Contents |
|---|---|
| `Main_Figures/` | Figures 1 and 2: one folder per panel, the data each panel reads, and the scripts that assemble the figures |
| `Supplementary_Figures/` | Supplementary Figures S1-S11, numbered as in the manuscript |
| `Simulation_Studies/` | The simulation pipelines behind Figure 1 and the simulation supplements (run on a cluster) |
| `Data_Applications/` | The ROSMAP eQTL fine-mapping runs behind Figure 2 (run on a cluster) |
| `R/` | Helpers shared by the figure scripts: paths, palettes and the top-N power/FDR panel |
| `run_all.R` | Rebuilds every figure from the data in this repository |

## Rebuilding the figures

You need R with the packages `ggplot2`, `cowplot`, `dplyr`, `tidyr`, `png`, `gridExtra` and `gtable`.
From the repository root:

```bash
Rscript run_all.R
```

This takes about a minute and rewrites every panel and figure under `Main_Figures/` and
`Supplementary_Figures/`. A single script is run the same way, always from the repository root:

```bash
Rscript Supplementary_Figures/S2/plot_S2.R
```

Last checked with R 4.6.1, ggplot2 4.0.3, cowplot 1.2.0, dplyr 1.2.1, tidyr 1.3.2, png 0.1.9,
gridExtra 2.3.1 and gtable 0.3.6.

## What each step needs

Everything `run_all.R` does uses only the data in this repository. The steps upstream of it need more:

| Step | Scripts | Needs |
|---|---|---|
| Figures from saved data | every step of `run_all.R` | the R packages above |
| Figure 1A diagram | `Main_Figures/figure_1/figure_1_panel_A/panel_A.tex` | XeLaTeX with TikZ, the Arial font, and a PDF rasterizer |
| Saved data from raw simulation output | the `extract_*` scripts | the raw simulation output, which is too large to include here |
| Model fits for Figure S5 | `Supplementary_Figures/S5/fit_models.R` | the `susieR` source tree and the `pkgload` package |
| Table behind Figure S7E | `Supplementary_Figures/S7/simulate_S7E.R` | a pinned `susieR` build and Python 3 (see the Supplementary Figures page) |
| Simulations | `Simulation_Studies/` | a cluster (see that page) |
| Real-data fine-mapping and AlphaGenome scoring | `Data_Applications/`, `Main_Figures/figure_2/data/data_generation/` | ROSMAP data access, a cluster, and an AlphaGenome API key |

Two environment variables say where the outside inputs are:

| Variable | Default | Used by |
|---|---|---|
| `SUSIER2_BENCHMARK_ROOT` | `../susieR2.0-benchmark` | the `extract_*` scripts |
| `SUSIER_ROOT` | `../susieR` | `Supplementary_Figures/S5/fit_models.R` |

To rebuild the Figure 1A diagram and then Figure 1:

```bash
cd Main_Figures/figure_1/figure_1_panel_A
xelatex panel_A.tex
pdftoppm -png -r 700 -singlefile panel_A.pdf panel_A
cd ../../..
Rscript Main_Figures/figure_1/create_figure1.R
```

## Conventions

- Main Figures 1 and 2 show SuSiE and SuSiE-inf. Supplementary Figures S6 and S10 are the versions
  with SuSiE-ash included.
- Simulation figures report 95% credible sets with purity of at least 0.5 and no correlation extension;
  Supplementary Figures S8 and S9 add the extension for comparison. The Supplementary Figures page has
  the details.
- The figure scripts contain no absolute paths: `R/paths.R` finds the repository root from the
  `.susieR2-paper-root` marker file. The cluster and pipeline scripts (`Simulation_Studies/`,
  `Data_Applications/`, steps 1 and 2 of the AlphaGenome scoring, and the track download scripts)
  are kept as they were run and still carry their original paths.

## Citing this work

If you find the `susieR` package or any of the source code in this
repository useful for your work, please cite both:

> Wang, G., Sarkar, A., Carbonetto, P. & Stephens, M. (2020). A
> simple new approach to variable selection in regression, with
> application to genetic fine mapping. *Journal of the Royal
> Statistical Society, Series B* **82**, 1273–1300.
> https://doi.org/10.1111/rssb.12388

> McCreight, A., Cho, Y., Li, R., Carbonetto, P., Nachun, D., Gan, H-Y., Kunkel, D.,
> Stephens, M., Denault, W.R.P. & Wang, G. (2025). SuSiE 2.0:
> improved methods and implementations for genetic fine-mapping and
> phenotype prediction. *bioRxiv*.
> https://doi.org/10.1101/2025.11.25.690514

If you use any of the summary data methods such as `susie_ss` or 
`susie_rss`, please also cite:

> Zou, Y., Carbonetto, P., Wang, G. & Stephens, M. (2022). Fine-mapping
> from summary data with the "Sum of Single Effects" model. *PLoS
> Genetics* **18**, e1010299. https://doi.org/10.1371/journal.pgen.1010299

If you use the Normal-Inverse-Gamma prior on residual variance estimates
(`estimate_residual_method = "NIG"`), please also cite:

> Denault, W.R.P., Carbonetto, P., Li, R., Alzheimer's Disease Functional
> Genomics Consortium, Wang, G. & Stephens, M. (2025). Accounting for
> uncertainty in residual variances improves calibration of the "Sum of
> Single Effects" model for small sample sizes. *bioRxiv*, 2025-05.
> Under review for *Nature Methods*.

If you use infinitesimal effects modeling (`unmappable_effects = "inf"`), 
please also cite:

> Cui, R., Elzur, R.A., Kanai, M. et al. (2024). Improving fine-mapping
> by modeling infinitesimal effects. *Nature Genetics* **56**, 162–169.
> https://doi.org/10.1038/s41588-023-01597-3

## License

MIT; see `LICENSE`.
