# Data Application

This section contains the data applications discussed in the ColocBoost paper. Below, we outline key components and protocols for analyzing data using the `colocboost` R package and related methodologies.


## 1. Vignettes for `colocboost`

This subsection provides a vignette demonstrating how to use the `colocboost::colocboost()` function and other related functions to analyze data. The [tutorials](https://statfungen.github.io/colocboost/) contains:

- **Input data format**: Standard input data formats to perform `colocboost`.
- **Running `colocboost`**: Example workflows for multi-trait colocalization analysis using individual level data and summary statistics.
- **Interpreting results**: Guidance on interpreting and visualization the output of `colocboost`.


---

## 2. xQTL Protocols

This subsection outlines protocols for xQTL analysis, including:

### a. Bioinformatics pipeline for ColocBoost
- Tutorial of bioinformatics pipeline for ColocBoost [tutorial link](https://statfungen.github.io/colocboost/articles/ColocBoost_Wrapper_Pipeline.html).
- Steps to perform multi-trait colocalization using ColocBoost and data preparation [protocol link](https://statfungen.github.io/xqtl-protocol/code/mnm_analysis/mnm_methods/colocboost.html).

### b. Comparison with fine-mapping using SuSiE
- Protocol of fine-mapping analysis using SuSiE [protocol link](https://statfungen.github.io/xqtl-protocol/code/mnm_analysis/univariate_fine_mapping_twas_vignette.html).
- Example workflows for identifying causal variants.

### c. Enrichment Analysis
- Steps to perform enrichment analysis using xQTL data [protocol link](https://statfungen.github.io/xqtl-protocol/code/enrichment/eoo_enrichment.html).
- Integration with external datasets for functional annotation.

### d. SLDSC Analysis
- Protocol for stratified linkage disequilibrium score regression (SLDSC) [protocol link](https://statfungen.github.io/xqtl-protocol/code/enrichment/sldsc_enrichment.html#workflow-steps).
- Example use cases for partitioning heritability.

---
