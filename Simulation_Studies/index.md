# Simulation Studies

This directory contains simulation designs and implementation codes used in the paper.


## 1. Simulation Design Overview

This section provides a summary of the simulation notebooks that implement a comprehensive benchmark of susieR 2.0 in both a sparse and oligogenic setting.

- **Simulation Schematics**: Comprehensive overview of the simulation framework used to evaluate and benchmark SuSiE with competing methods.

- **Phenotype Data Simulation**: Establishes the fundamental simulation framework for generating synthetic phenotype data. Simulates phenotype data (y vector) for based on real genotype data (X matrix) using varying levels of total heritability approaches. Configurable for different numbers of sparse effects with controllable narrow-sense heritability.

- **Run susieR2.0**: Executes the SuSiE algorithm on simulated datasets to identify causal variants. Processes and standardizes results for performance evaluation with key output metrics.

- **Result Summary**: Calculates performance metrics including power and false discovery rates, AUROC/AUPRC curves, and prediction metrics from method results. Generates comparison tables summarizing method effectiveness across simulation scenarios.

## 2. References

[1] [simxQTL](https://github.com/StatFunGen/simxQTL): In house simulation R package to support investigations of various QTL association methods.

[2] [susieR]()

[3] [susieR for SS/RSS data]()

[4] [ash methodology]()

[5] [mr.ash]()

[6] [susie-inf]()
