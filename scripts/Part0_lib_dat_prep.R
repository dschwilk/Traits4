## Part0 lib_dat_prep.R
## load require library and data transformation.
## need to run at the beginning

# library for part1
library(ape)
library(tidyverse)
library(phytools)
library(ggrepel)
library(ggtree)

source("scripts/func_prun_replac.R")
source("scripts/func_bind_sister_tip.R")
source("scripts/func_phylosig.R")
source("scripts/func_model_aicc.R")

#source("scripts/func_phylosig_parallel.R")

###############################################################################
####################### 0. Import Data: Data Prep #############################
###############################################################################

######################## 0.1 Master Tree Import ###############################
treeVascularPlants <- read.tree("./data/Vascular_Plants_rooted.tre")

######################## 0.2 Trait Value Data #################################
# The number of replicate in measurements is ten
# Each morphological traits were measured from ten difference seeds
REP = 10 


InputMatrix <- read.csv("./data/RawData.csv") # Import Raw Measurement Data

# Germination data
Germination = readxl::read_xlsx("./data/Germination(sd_se).xlsx") %>%
  mutate(mean_Germination = mean) %>%
  mutate(se_Germination = se) %>%
  dplyr::select(Species, mean_Germination, se_Germination)  # Import

# Trait measurements 
dat = InputMatrix %>% group_by(Family, Species, Classification) %>% 
  dplyr::summarise (mean_Area = mean(Area), se_Area = sd(Area)/sqrt(REP), 
                    mean_Height = mean(Height), se_Height = sd(Height)/sqrt(REP),
                    mean_Mass = mean(Mass), 
                    se_Mass = sd(Mass/sqrt(REP))) %>% as.data.frame()


dat_G = merge(dat, Germination)

dat_G_scaled = dat_G %>%
  mutate(scaled_Area = scale(mean_Area)) %>%
  mutate(scaled_Height = scale(mean_Height)) %>%
  mutate(scaled_Mass = scale(mean_Mass)) %>%
  mutate(scaled_Germination = scale(mean_Germination)) %>%
  dplyr::select(scaled_Area, scaled_Height, scaled_Mass, scaled_Germination)

