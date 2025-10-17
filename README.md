# susieR 2.0 Manuscript Resources

Code and data to reproduce figures in susieR 2.0 manuscript.


## About susieR 2.0

susieR 2.0 is a complete architectural redesign that addresses code duplication and fragmented architecture of the original SuSiE implementation while adding new features including unmappable effects modeling and variuos performance optimizations, all while maintaing full backward compatibility.

susieR 2.0 eliminates this duplicative architecture through a unified framework built on modular design principles. The user-facing interface remains unchanged, but the implementation now uses generic functions with data-type specific backends through R’s S3 dispatch system. This architecture enables the integration of new SuSiE extensions while maintaining identical results to their original versions. Beyond architectural improvements, this release introduces substantial algorithmic advances including support for unmappable effects modeling, enhanced computational speed for regularized LD matrices, new convergence criteria, and improved refinement procedures.

## Repository Structure

This repository contains the codes and data used to generate all figures from our manuscript, available at: https://github.com/StatFunGen/susieR2.0-paper

- **Simulation Studies**: This directory contains simulation designs and implementation codes used in the paper. 

- **Data Applications**: This section contains the data applications discussed in the susieR 2.0 paper.

- **Main Figures**: Each notebook is fully executable and documented to ensure reproducibility of our results.

- **Supplementary Figures**: Each notebook is fully executable and documented to ensure reproducibility of our results shown in Supplementary Figures. 

## Getting Started

To navigate this resource, use the table of contents in the left sidebar. Each figure section contains interactive notebooks that allow you to:

1. View the code used to generate analyses
2. Examine data associated with the figures
3. Reproduce visualizations

## Computational Requirements

The analyses in this book were performed using:
- R version 4.1 or higher
- Key R packages: data.table, ggplot2, dplyr

## susieR 2.0 Tutorial Website

Learn how to perform colocalization analysis with step-by-step examples. For detailed tutorials and use cases in [Tutorials](https://statfungen.github.io/susieR/articles/index.html).

## Citation

If you use susieR 2.0 in your research, please cite:
