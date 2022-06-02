# Traits4


## Overview

All code is in `scripts`. R code assumes the working directory is the project root. 

## Directory structures ##

```
Main directory (git root)
  README.md
  data
    Vascular_Plants_rooted.tre
    RawData.csv
    Germination(sd_se).xlsx
  results (directory contents not in version control)
    .gitignore (with two lines, see below)
    Figure1_Phylogenetic_compara.pdf
    Figure2_Phylogenetic_Position.pdf
    table2_phylogenetic_signal.csv
    appendix_model_selection.csv
    resample_k_distr_Height.pdf
    resample_k_distr_Mass.pdf
    resample_k_distr_Area.pdf
    resample_k_distr_Germination.pdf
    resample_Height.pdf
    resample_Mass.pdf
    resample_Area.pdf
    resample_Germination.pdf
  scripts
  	 lib_dat_prep.R
  	 Part1_PhyloCompara.R
  	 Part2_Sim_Para.R
    func_prun_replac.R
    func_bind_sister_tip.R
    func_phylosig.R
    func_model_aicc.R
    func_phylosig_parallel.R
```
