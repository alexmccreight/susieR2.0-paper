# Supplementary Figures

This section contains the scripts used to assemble each supplementary figure from pre-computed results. The upstream simulation and analysis code lives under **Simulation Studies** and **Data Applications**.

- **Figure S1**: Sparse fine-mapping benchmarks comparing SuSiE, SuSiE-ash, and SuSiE-inf on power, FDR, PIP ROC curves, and PIP calibration across K = 1–5 causal variants
- **Figure S2**: Unmappable-effects worked example showing PIPs and credible sets from SuSiE (L=10), SuSiE+BB, SuSiE-inf, and SuSiE-ash on a single simulated dataset
- **Figure S3**: TWAS cross-validation prediction R² and per-replicate winner proportions for SuSiE, SuSiE-ash, and SuSiE-inf across sparse and complex (oligogenic-plus-infinitesimal) scenarios
- **Figure S4**: Sparse fine-mapping comparison of SuSiE credible-set variants: default (purity ≥ 0.5), high-purity filter (≥ 0.8), SuSiE+BLiP, and SuSiE+attainable coverage on power, FDR, purity, and CS size across K = 1–5 causal variants
- **Figure S5**: Credible-set count by size bin and median purity for SuSiE, SuSiE-ash, and SuSiE-inf across sparse and three complex (oligogenic-plus-infinitesimal) scenarios
- **Figure S6**: Runtime speedup of susieR 2.0 over susieR 1.0 for RSS fine-mapping with regularization (`susie_rss_lambda`), shown as fold reduction in median runtime across feature/sample (p/n) ratios at p = 5000, for both LD-matrix and genotype-matrix inputs
- **Figure S7**: Sparse fine-mapping power and FDR for SuSiE, SuSiE-ash, and SuSiE-inf comparing credible sets without LD extension against credible sets with LD extension (0.99), shown as two bars per method across K = 1–5 causal variants
- **Figure S8**: Complex fine-mapping power and FDR versus the top-N causal-variant definition, comparing credible sets without LD extension against credible sets with LD extension (0.99) for each method (six lines per panel), arranged by simulation scenario (rows: Complex S2, Complex S1, Complex) and metric (columns: power, FDR)
- **Figure S9**: Per-method runtime comparison of SuSiE, SuSiE-ash, and SuSiE-inf, shown as mean wall-clock fit time (± SE) per replicate across 500 replicates of the complex (oligogenic-plus-infinitesimal) simulation scenario
- **Figure S10**: Complex fine-mapping power and FDR versus the top-N causal-variant definition for SuSiE, SuSiE-ash, and SuSiE-inf, across the Complex S1 and Complex S2 scenarios. These were panels C–F of Figure 1; the main figure now shows SuSiE and SuSiE-inf only
- **Figure S11**: Real-data comparison of SuSiE, SuSiE-ash, and SuSiE-inf on ROSMAP eQTL data — credible-set counts by size bin, cross-method signal concordance (Jaccard ≥ 0.75), TWAS cross-validation R², AlphaGenome CS scores by concordance group, and a locus example. This was Figure 2; the main figure now shows SuSiE and SuSiE-inf only
