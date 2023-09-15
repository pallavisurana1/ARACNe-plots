# ARACNe Network Plots

This repository contains an R script for constructing network plots based on ARACNe output files. ARACNe (Algorithm for the Reconstruction of Accurate Cellular Networks) is an information-theoretic algorithm designed to infer gene regulatory networks.

## Reference

Margolin, A.A., Nemenman, I., Basso, K. et al. ARACNe: An Algorithm for the Reconstruction of Gene Regulatory Networks in a Mammalian Cellular Context. BMC Bioinformatics 7, S7 (2006). [DOI: 10.1186/1471-2105-7-S1-S7](https://doi.org/10.1186/1471-2105-7-S1-S7)

## Dependencies

The following R packages are required for running the script:

- `BiocManager`
- `reshape2`
- `EnsDb.Hsapiens.v86`
- `dplyr`
- `tidyverse`
- `igraph`
- `networkD3`
- `magrittr`
- `linkcomm`

You can install these packages using the `pacman` package manager for R.

## Instructions

### 1. Data Import

The script reads ARACNe output and a list of genes of interest from specified file paths.

### 2. Annotation

It maps Ensembl IDs for regulators and targets to their corresponding gene symbols using the `EnsDb.Hsapiens.v86` package.

### 3. Data Processing

The script subsets the data to focus only on genes of interest and then processes it to prepare for different levels of the network.

### 4. Network Construction

It uses the `linkcomm` and `networkD3` libraries to generate network plots for different hierarchical levels of the gene network.

### Helper Functions

The `process_data_level()` function is introduced to make the code reusable and more readable by handling data manipulation at various levels.

## Running the Code

To execute the script, simply source it in your R environment after installing the required dependencies and setting the file paths appropriately for the ARACNe output and genes of interest.

## Plotting

The script includes functionality to plot the networks both as static PDFs and as interactive HTML files.

