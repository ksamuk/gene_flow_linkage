
[![License](http://img.shields.io/:license-mit-blue.svg)](http://doge.mit-license.org)

# recombination_bias

This repository contains all the scripts and processed data files used to perform the data processing and analyses described in:

**Gene flow favors linkage between adaptive alleles in a globally distributed species**

*Kieran Samuk, Greg Owens, Diana Rennison, Kira Delmore, Sara Miller and Dolph Schluter.*

### Requirements

- **R version 3.2.2**. Certainly dependencies will only work with this version of R.
- **[Rstudio](https://www.rstudio.com/) (any version that supports packrat)**. This is streamlines the use of [packrat](https://rstudio.github.io/packrat/).

### How to use this repo

Reccomended method:

1. Clone this repository using Rstudio
2. Open the provide .Rproj file
3. packrat will initialize the provided dependency library.
4. Run/browse code etc.

The following top-level scripts will run 'as-is':

* 4_analyze_75k_permute_S2.R
* 5_analyze_clustering_output_S3.R
* 6_plot_figure_1.R
* 7_plot_figures_2.R
* 8_plot_figures_S4_S5.R
* 9_plot_figure_S1.R
* 10_prepare_supplemental_data_files.R

The remaining below require precursors to function (these will be made available on Dryad, accession number to follow):

* 1_process_75k_fit_models.R
* 2_process_snp_clustering_fst.R
* 3_process_clustering_output.R

