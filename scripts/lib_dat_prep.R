## lib_dat_prep.R

## Reads in data and does simple summarizing. Code assumes working directory is
## root of the git repository.

## load required libraries
library(ape)
library(nlme)
#library(geiger)
#library(corpcor)
#library(nloptr)
library(RColorBrewer)
#library(OUwie)
library(readxl)
#library(plyr) ## Dangerous to mix these!
#library(tidyverse)  ## Best not to dump this in.
library(dplyr)
library(phytools)
library(MASS)
library(ggplot2)
library(ggrepel)
library(AICcmodavg)
library(broom)
library(patchwork)
#library(ggtree)  ## Where is this from?

## DWS: ah, bioconductor! This is non standard so make sure to document:

## if (!require("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")

## BiocManager::install("ggtree")


## DWS: This is a large number of dependencies! Some these are also redundant
## with one another. In production code, I'd prefer to see explicit library
## calls rather than having so much dumped into the global environment. Some of
## these are never even used.


###############################################################################
####################### 0. Import Data: Data Prep #############################
###############################################################################

######################## 0.1 Master Tree Import ###############################
treeVascularPlants <- read.tree("./data/Vascular_Plants_rooted.tre")

######################## 0.2 Trait Value Data #################################
InputMatrix <- read.csv("./data/RawData.csv") # Import Raw Measurement Data

Germination = read_xlsx("./data/Germination(sd_se).xlsx") %>%
  mutate(mean_Germination = mean) %>%
  mutate(se_Germination = se) %>%
  dplyr::select(Species, mean_Germination, se_Germination)  # Import
                                                            # Germination Data

Rep = 10 # Number of replicates in measurements
## DWS: why is this hard coded?

dat = InputMatrix %>% group_by(Family, Species, Classification) %>% 
  dplyr::summarise (mean_Area = mean(Area), se_Area = sd(Area)/sqrt(Rep), 
                    mean_Height = mean(Height), se_Height = sd(Height)/sqrt(Rep),
                    mean_Mass = mean(Mass), 
                    se_Mass = sd(Mass/sqrt(Rep))) %>% 
  as.data.frame()
