
[![License](http://img.shields.io/:license-mit-blue.svg)](http://doge.mit-license.org)

# gene_flow_linkage

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

### Repository contents, by folder

#### analysis_ready
Processed divergence data for populations used in the project. This is the 'raw'(ish) data used for the analysis. These files are generated by scripts 1-3 (the 'process' scripts). 

#### ev_prep_scripts
The scripts used to generate the explanatory variables ('evs;) used in our regressions. Recombination rate, mutation rate, and genet density calculations are found here.

#### evs
The explanatory variables used in the regressions (recombination, mutation, gene density). These are generated by the ev_prep_scripts

#### figures
These are the figures from the manuscript. The 'plot' scripts output figures to the `figures/raw` folder. The figures in the root `figures` folder are the good copies of the figures, post formating in Inkscape.

#### meta_data
Miscellaneous metadata, including population locations, overall FST estimates, p-values, etc.

#### process_raw_data
The scripts used to do all analyses upstream of the ones outlined here. These scripts move from raw .fastq data (e.g. from a database), to population genetic summary stastistics. A separate readme is provided for this section.

#### shared_functions
The library of functions used by the top-level scripts.
