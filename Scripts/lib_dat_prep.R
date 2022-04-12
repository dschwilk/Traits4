## set working directory
setwd("~/Desktop/Traits_Validated") 

## load required libraries
library(ape)
library(nlme)
library(geiger)
library(corpcor)
library(nloptr)
library(RColorBrewer)
library(OUwie)
library(readxl)
library(plyr)
library(tidyverse)
library(phytools)
library(MASS)
library(ggplot2)
library(ggrepel)
library(AICcmodavg)
library(MASS)
library(broom)
library(patchwork)
library(ggtree)

##########################################################################################
####################### 0. Import Data: Data Prep ########################################
##########################################################################################

######################## 0.1 Master Tree Import ##########################################
treeVascularPlants <- read.tree("Input/Vascular_Plants_rooted.tre")

######################## 0.2 Trait Value Data #############################################
InputMatrix <- read.csv("Input/RawData.csv") # Import Raw Measurement Data

Germination = read_xlsx("Input/Germination(sd_se).xlsx") %>%
  mutate(mean_Germination = mean) %>%
  mutate(se_Germination = se) %>%
  dplyr::select(Species, mean_Germination, se_Germination)  # Import Germination Data

Rep = 10 # Number of replicates in measurements

dat = InputMatrix %>% group_by(Family, Species, Classification) %>% 
  dplyr::summarise (mean_Area = mean(Area), se_Area = sd(Area)/sqrt(Rep), 
                    mean_Height = mean(Height), se_Height = sd(Height)/sqrt(Rep),
                    mean_Mass = mean(Mass), 
                    se_Mass = sd(Mass/sqrt(Rep))) %>% 
  as.data.frame() 