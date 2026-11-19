# Supplementary Figures

This section contains the scripts used to assemble each supplementary figure from pre-computed results. The upstream simulation and analysis code lives under **Simulation Studies** and **Data Applications**.

- **Figure S1**: Sparse fine-mapping benchmarks comparing SuSiE, SuSiE-ash, and SuSiE-inf on power, FDR, PIP ROC curves, and PIP calibration across K = 1–5 causal variants (credible sets: 95% coverage, exact purity ≥ 0.5, no correlation extension)
- **Figure S2**: Unmappable-effects worked example showing PIPs and credible sets from SuSiE (L=10), SuSiE-inf, and SuSiE-ash on a single simulated dataset
- **Figure S3**: TWAS cross-validation prediction R² and per-replicate winner proportions for SuSiE, SuSiE-ash, and SuSiE-inf across sparse and complex (oligogenic-plus-infinitesimal) scenarios
- **Figure S4**: Sparse fine-mapping comparison of SuSiE credible-set variants: default (purity ≥ 0.5), high-purity filter (≥ 0.8), SuSiE+BLiP, and SuSiE+attainable coverage on power, FDR, purity, and CS size across K = 1–5 causal variants, with an additional panel showing the full posterior mass of CSs reported by SparsePro (2023) and susieR
- **Figure S5**: Credible-set count by size bin and median purity for SuSiE, SuSiE-ash, and SuSiE-inf across sparse and three complex (oligogenic-plus-infinitesimal) scenarios (un-extended credible sets, as in Figure 1)
- **Figure S6**: Runtime speedup of susieR 2.0 over susieR 1.0 for RSS fine-mapping with regularization (`susie_rss_lambda`), shown as fold reduction in median runtime across feature/sample (p/n) ratios at p = 5000, for both LD-matrix and genotype-matrix inputs
- **Figure S7**: Sparse fine-mapping power and FDR for SuSiE, SuSiE-ash, and SuSiE-inf comparing credible sets without LD extension against credible sets with LD extension (0.99), shown as two bars per method across K = 1–5 causal variants
- **Figure S8**: Complex fine-mapping power and FDR versus the top-N causal-variant definition, comparing credible sets without LD extension against credible sets with LD extension (0.99) for each method (six lines per panel), arranged by simulation scenario (rows: Complex S2, Complex S1, Complex) and metric (columns: power, FDR)
- **Figure S9**: Per-method runtime comparison of SuSiE, SuSiE-ash, and SuSiE-inf, shown as mean wall-clock fit time (± SE) per replicate across 500 replicates of the complex (oligogenic-plus-infinitesimal) simulation scenario
- **Figure S10**: Complex fine-mapping power and FDR versus the top-N causal-variant definition for SuSiE, SuSiE-ash, and SuSiE-inf, across the Complex S1 and Complex S2 scenarios. Built by `S10/build_S10_panels.R` from the same summary tables as Figure 1 panels C–F (which show SuSiE and SuSiE-inf only), using the shared panel function in `R/topn_panels.R`
- **Figure S11**: Real-data comparison of SuSiE, SuSiE-ash, and SuSiE-inf on ROSMAP eQTL data — credible-set counts by size bin, cross-method signal concordance (Jaccard ≥ 0.75), TWAS cross-validation R², AlphaGenome CS scores by concordance group, and a locus example. This was Figure 2; the main figure now shows SuSiE and SuSiE-inf only. Panel D is rebuilt from the three-method scores by `run_all.R`

Credible-set conventions: every simulation figure reports 95% credible sets with exact purity ≥ 0.5 computed from the in-sample correlation matrix and **no** correlation extension (`version = "xcorr_noext"` in the `extract_*` scripts), matching the susieR default and the real-data analysis. Figures S7 and S8 are the only ones that add the optional 0.99 correlation extension, for comparison. A replicate that reports no credible set contributes power 0 and FDR 0. In the complex scenarios the top-N causal set is capped at the number of nonzero effects.

## Manuscript numbering

The manuscript numbers supplementary figures by order of first mention, which differs from the directory numbers used here. The `figures/supplementary/figS*.pdf` files of the manuscript repository are copies of the `S*/S*.pdf` outputs below.

| Manuscript figure | Repository directory | Content |
|---|---|---|
| S1 | S6 | regularized RSS backend speedup |
| S2 | S1 | fine-mapping under sparse settings |
| S3 | S5 | credible set size and purity |
| S4 | S3 | prediction R² in simulation |
| S5 | S2 | illustrative locus |
| S6 | S10 | complex settings, three methods |
| S7 | S4 | alternative credible set coverage strategies |
| S8 | S7 | correlation-based extension, sparse settings |
| S9 | S8 | correlation-based extension, complex settings |
| S10 | S11 | real data application, three methods |
| S11 | S9 | runtime comparison |

Main figures: `Main_Figures/figure_1/susieR2_figure1.pdf` and `Main_Figures/figure_2/susieR2_figure2.pdf`.

## Rebuilding the three-method locus panel

The repository's S11 is **Figure S10 in the manuscript**. Its locus panel is
constructed by the same production scripts as Figure 2E, with `--include-ash`
to include all three methods:

```bash
Rscript Main_Figures/figure_2/figure_2_panel_E/prepare_panel_E_data.R --include-ash
Rscript Main_Figures/figure_2/figure_2_panel_E/plot_panel_E.R --include-ash
Rscript Supplementary_Figures/S11/plot_S11.R
```

These steps use the vendored fit objects, gene models, and annotation BED files.
They write the three-method plotting data and plot object to `S11/data/`, then
assemble `S11/S11.pdf` and `S11/S11.png`. Colors are specified in
`Main_Figures/figure_2/figure_2_panel_E/plot_panel_E.R` when the plot is constructed.
Running the panel scripts without `--include-ash` builds the two-method main figure.
