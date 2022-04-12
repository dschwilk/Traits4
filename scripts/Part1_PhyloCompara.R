#

#setwd("~/Desktop/Research Projects/Traits/GitRepo/Traits4") 

## Never use absolute paths in code and never set the working directory in
## code. Also, this is a different working directory than was set in
## `lib_dat_prep.R`!


##########################################################################################
############################## 0. Import Raw Measurement Data ############################
##########################################################################################
##library(ape)
## library(nlme)
## library(geiger)
## library(corpcor)
## library(nloptr)
## library(RColorBrewer)
## library(OUwie)
## library(readxl)
## library(plyr)
## library(tidyverse)
## library(phytools)
## library(MASS)
## library(ggplot2)
## library(ggrepel)
## library(AICcmodavg)
## library(MASS)
## library(broom)
## library(patchwork)
## library(ggtree)

## Trait Value Data
InputMatrix <- read.csv("./data/RawData.csv")  # Import Raw Measurement Data

Rep = 10 # Number of replicates in measurements

## DWS: hard coded data?

appendix = InputMatrix %>% group_by(Family, Species, Classification) %>% 
  dplyr::summarise (mean_Area = mean(Area), se_Area = sd(Area)/sqrt(Rep), 
             mean_Height = mean(Height), se_Height = sd(Height)/sqrt(Rep),
             mean_Mass = mean(Mass), 
             se_Mass = sd(Mass/sqrt(Rep))) 

Germination = read_xlsx("data/Germination(sd_se).xlsx") %>%
  mutate(mean_Germination = mean) %>%
  mutate(se_Germination = se) %>%
  dplyr::select(Species, mean_Germination, se_Germination) 

##########################################################################################
############################## 1. Building Phylogeny #####################################
##########################################################################################
## Prepare Phylogenetic Tree
treeVascularPlants <- read.tree("data/Vascular_Plants_rooted.tre")
tips <-treeVascularPlants$tip.label
Genera<-unique(sapply(strsplit(tips,"_"),function(x) x[1]))

#PruneTree Function
PruneTree <- function(x){
  ## x, a list of desired species in "xx_xx" format
  DesiredSpecies <- unique(x)
  ## Desired Genera
  DesiredGenera <- sapply(strsplit(DesiredSpecies, "_"),function(x) x[1])
  
  ## TreeSpecies Data Frame
  SpeciesGeneraSpecies <- tibble(TreeSpecies = treeVascularPlants$tip.label,
                                     TreeGenera = sapply(strsplit(tips,"_"),function(x) x[1]),
                                     DesiringGenera = TreeGenera %in% intersect(Genera,DesiredGenera) ,
                                     DesiringSpecies = TreeSpecies %in% DesiredSpecies) 
  
  SpeciesListSpecies <- filter(SpeciesGeneraSpecies, DesiringSpecies == "TRUE")
  
  GeneraSpecies <- filter(SpeciesGeneraSpecies, SpeciesGeneraSpecies$TreeGenera %in% 
                            setdiff(DesiredGenera, SpeciesListSpecies$TreeGenera))
  SpeciesListGenera <- group_by(GeneraSpecies, TreeGenera) %>% group_modify(~ head(.x, 1L))
  
  LISTALLSPECIES <- rbind.data.frame(SpeciesListSpecies, SpeciesListGenera)
  
  treeTestedSpecies <- keep.tip(treeVascularPlants, LISTALLSPECIES$TreeSpecies)
  
  ## Replacing tip label
  
  aaa <- data_frame(DesiredSpecies = DesiredSpecies, 
                    TreeGenera = sapply(strsplit(DesiredSpecies,"_"),function(x) x[1]))
  bbb <- merge(aaa, SpeciesListGenera)

   ## DWS: where is mapvalues function from? Ah, old plyr! Rewrite!
  treeTestedSpecies$tip.label <- plyr::mapvalues(treeTestedSpecies$tip.label,
                                           c(bbb$TreeSpecies), c(bbb$DesiredSpecies))
  tree_x <- treeTestedSpecies
  plotTree(tree_x, ftype="i")
  return(tree_x)
}

TreeInputMatrix <- PruneTree(unique(InputMatrix$Species))

b = setdiff(unique(InputMatrix$Species), TreeInputMatrix$tip.label)


## DWS: why is there hard coded data here? Why is this not from a data file?

## Adding "Callirhoe_leiocarpa"
tip1 <- "Callirhoe_leiocarpa"
sister1 <- "Callirhoe_involucrata"
tree1 <- bind.tip(TreeInputMatrix,tip1,
                  edge.length = 0.5*TreeInputMatrix$edge.length[which(TreeInputMatrix$edge[,2]==
                                                                        which(TreeInputMatrix$tip.label==sister1))],
                  where=which(TreeInputMatrix$tip.label==sister1),
                  position=0.5*TreeInputMatrix$edge.length[which(TreeInputMatrix$edge[,2]==
                                                                   which(TreeInputMatrix$tip.label==sister1))])
## Adding "Digitaria_californic"
tip2 <- "Digitaria_californica"
sister2 <- "Digitaria_ciliaris"
TreeAllMatrix <- bind.tip(tree1,tip2,
                          edge.length = 0.5*tree1$edge.length[which(tree1$edge[,2]==
                                                                      which(tree1$tip.label==sister2))],
                          where = which(tree1$tip.label==sister2),
                          position=0.5*tree1$edge.length[which(tree1$edge[,2]==
                                                                 which(tree1$tip.label==sister2))])

plotTree(TreeAllMatrix, ftype="i")
##########################################################################################
############################## 2. Calculate Phylogenetic Positiion #######################
##########################################################################################
## Phylogenetic Position
par(mar=c(5,5,5,1))
dist <- cophenetic.phylo(TreeAllMatrix)
phyloposi <-isoMDS(dist) %>% as.data.frame()

phyloposi_species = phyloposi %>% 
  mutate(Species = row.names(phyloposi)) %>% 
  mutate(phy1 = round(scale(points.1),digits = 1)) %>% 
  mutate(phy2 = round(scale(points.2), digits = 1))

AllMatrix = merge (appendix, phyloposi_species)

phyloposi_family = AllMatrix %>% group_by(Family) %>% 
  summarise(Family, Classification, phy1, phy2) %>%
  distinct() 

p_family = ggplot(phyloposi_family, aes(x=phy1, y=phy2, color=Classification)) +
  geom_point() + labs(x="phy1", y="phy2") + 
  geom_text_repel(aes(label = Family), size =3.5) + 
  theme_bw()

print(p_family) 

##########################################################################################
################### 3. Array Trait Value with Phylogenetic Tree ##########################
##########################################################################################
# Phylogenetic Tree: TreeAllMatrix


# Trait Data: appendix
A = TreeAllMatrix$tip.label
#names(A) = c(1: length(A))

appendix_G = merge(appendix, Germination)

B1 = appendix_G %>% mutate (Species1 = factor(Species, level = rev(A)),
                  Classification1 = factor(Classification, levels = c("Asteracea", "Other_Dicot", "Monoct"), 
                                                   ordered = TRUE)) 

G1 = B1[with(B1, order(Classification1, Species1)),] 
#write_csv(G1, "Input/Reordered_AllMatrix.csv")
G2 = G1 %>% mutate(Species1 = NULL, Classification1 = NULL, SpeciesNum = 1:length(TreeAllMatrix$tip.label))


#G4 = G1 %>% mutate(Species1 = NULL, Classification1 = NULL, SpeciesNum = 1:length(TreeAllMatrix$tip.label))

#3.1 Trait Value Figures (Seed Mass and Seed Height, Colored)
Mass_Figure <- ggplot(G2, aes(y = mean_Mass, x = reorder(Species, -SpeciesNum))) + 
  geom_point(aes(shape = Classification, color = Classification)) + 
  geom_errorbar(aes(ymin=mean_Mass - se_Mass, ymax=mean_Mass + se_Mass, color = Classification))+
  labs(y = "Mass(log)") +
  theme_bw() + 
  theme (axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         ) +
  coord_flip()
#Mass_Figure

Height_Figure <- ggplot(G2, aes(y = mean_Height, x = reorder(Species, -SpeciesNum))) + 
  geom_point(aes(shape = Classification, color = Classification)) + 
  geom_errorbar(aes(ymin=mean_Height - se_Height, ymax=mean_Height + se_Height, color = Classification))+
  labs(y = "Height (mm)") +
  theme_bw() +
  theme (axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
        ) +
  coord_flip()
#Height_Figure

Area_Figure <- ggplot(G2, aes(y = mean_Area, x = reorder(Species, -SpeciesNum))) + 
  geom_point(aes(shape = Classification, color = Classification)) + 
  geom_errorbar(aes(ymin=mean_Area - se_Area, ymax=mean_Area + se_Area, color = Classification))+
  labs(y = "Area (sq mm)") +
  theme_bw() +
  theme (axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
  ) +
  coord_flip()
#Area_Figure

Germination_Figure <- ggplot(G2, aes(y = mean_Germination, x = reorder(Species, -SpeciesNum))) + 
  geom_point(aes(shape = Classification, color = Classification)) + 
  geom_errorbar(aes(ymin=mean_Germination - se_Germination, ymax=mean_Germination + se_Germination, color = Classification))+
  labs(y = "Germination (%)") +
  theme_bw() +
  theme (axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
  ) +
  coord_flip()
#Germination_Figure


#3.2 Plotting/Coloring Phylogenetic Tree, Colored
# http://bioconductor.riken.jp/packages/3.7/bioc/vignettes/ggtree/inst/doc/treeAnnotation.html#labelling-associated-taxa-monophyletic-polyphyletic-or-paraphyletic
# https://yulab-smu.top/treedata-book/chapter12.html 
# https://groups.google.com/g/bioc-ggtree/c/BA7g-npY1BM

classification = list(asteracea = G2$Species[1:8],
                      other_dicoct = G2$Species[9:30],
                      monocot = G2$Species[31:45])

t = ggtree(TreeAllMatrix, branch.length ='none')

t2 = groupOTU(t, classification, 'Classification') + aes(color = Classification) + geom_treescale(x=20, y =1, offset = 2, color = "white") + geom_tiplab(size = 3.8, fontface = 'italic') + theme(legend.position = "none")

t2

#3.3 Patched Figures

(t2 + guides(colour = "none")) + Mass_Figure + Height_Figure + Area_Figure + Germination_Figure + plot_layout(widths = c(2,1,1,1,1), guides = "collect") & theme(legend.position = 'bottom')





##########################################################################################
############################## 4. Calculate Phylogenetic Signal ##########################
##########################################################################################

AllMatrix_G = merge(AllMatrix, Germination)

### Format data before test phylogenetic signal
Mass.All <- AllMatrix_G$mean_Mass
names(Mass.All) <- AllMatrix$Species

Height.All <- AllMatrix_G$mean_Height
names(Height.All) <- AllMatrix$Species

Area.All <- AllMatrix_G$mean_Area
names(Area.All) <- AllMatrix$Species

Germination.All <- AllMatrix_G$mean_Germination
names(Germination.All) <- AllMatrix_G$Species

### Test phylogenetic signal use Blomberg's K

phyloSigMass.k <- phylosig(TreeAllMatrix, Mass.All, method="K", test=TRUE, nsim=999)

phyloSigHeight.k <- phylosig(TreeAllMatrix, Height.All, method="K", test=TRUE, nsim=999)

phyloSigArea.k <- phylosig(TreeAllMatrix, Area.All, method="K", test=TRUE, nsim=999)

phyloSigGermination.k <- phylosig(TreeAllMatrix, Germination.All, method="K", test=TRUE, nsim=999)

### Phylogenetic signal display
# table display results
tab_phyloSig = matrix(c(phyloSigMass.k$K, phyloSigHeight.k$K, phyloSigArea.k$K, 
                        phyloSigGermination.k$K, phyloSigMass.k$P, phyloSigHeight.k$P, phyloSigArea.k$P, 
                        phyloSigGermination.k$P), ncol = 2)

colnames(tab_phyloSig) = c("Blomberg's K", "P-value")

Variables = c("Mass", "Height", "Area", "Germination")
tab_phyloSig1 = cbind(Variables,tab_phyloSig)
tab_phyloSig1



##########################################################################################
############################## 4.Model Selection  ########################################
##########################################################################################


AllMatrix_G_scaled = AllMatrix_G %>% 
  mutate(scaled_Area = scale(mean_Area)) %>%
  mutate(scaled_Height = scale(mean_Height)) %>%
  mutate(scaled_Mass = scale(mean_Mass)) %>%
  mutate(scaled_Germination = scale(mean_Germination)) %>%
  dplyr::select(scaled_Area, scaled_Height, scaled_Mass, scaled_Germination, phy1, phy2)

vars = names(AllMatrix_G_scaled[-4])
models = list()
for (i in seq_along(vars)){
  vc = combn(vars, i)
  for (j in 1:ncol(vc)){
    model = as.formula(paste0('scaled_Germination ~', paste0(vc[, j], collapse = '+')))
    models = c(models, model)
  }
}

glmmodels = lapply(models, function(x) glm(x, data = AllMatrix_G_scaled))

# AICc calculation

AICc = sapply(models, function(x) AICc(glm(x, data = AllMatrix_G_scaled), 
                                       return.K = FALSE))
names(AICc) = sapply(glmmodels, function(x) x$formula)
AICc
min(AICc)
max(AICc)

## tablePlot
tab_models = data_frame(names(AICc),AICc)
colnames(tab_models) = c("Model", "AICc")

tab_models

# Correlation between mass and height 
cor.test(AllMatrix_G_scaled$scaled_Height, AllMatrix_G_scaled$scaled_Mass, 
         method = "pearson")
cor.test(AllMatrix_G_scaled$scaled_Mass, AllMatrix_G_scaled$scaled_Area, 
         method = "pearson")
cor.test(AllMatrix_G_scaled$scaled_Height, AllMatrix_G_scaled$scaled_Area, 
         method = "pearson")



