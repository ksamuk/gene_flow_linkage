
[![License](http://img.shields.io/:license-mit-blue.svg)](http://doge.mit-license.org)

# recombination_bias

This repository contains all the scripts and processed data files used to perform the data processing and analyses described in:

**Gene flow favors linkage between adaptive alleles in a globally distributed species**

*Kieran Samuk, Greg Owens, Diana Rennison, Kira Delmore, Sara Miller and Dolph Schluter.*

### How to use this repo

The scripts hosted here fall into two broad categories: 

1. Prepartory scripts that rely on raw DNA sequence data - authored largely by Greg Owens 
2. Analytical scripts that rely on the output of the prepartory scripts - authored largely by Kieran Samuk

The scripts in category (1) above require one to aquire the raw DNA sequence data described in Table S1 of our paper, install the dependancies listed in the supplemental material, and have a general proficiency with the UNIX shell. In other words, they require some outside preparation and will not work 'out-of-the-box'.

The scripts in category (2) include all necessary data can be run 'as-is' by cloning this repository and running the individual scripts in R. The easiest way to do this would be to open the included .Rproj file using Rstudio. Note you will need to install several necessary R packages in order for these scripts to function.
