##############################################################################################################
############################# Function of Calculation Phylogenetic Signal ####################################
##############################################################################################################

#### Requirement inputs include phylogenetic tree and its comparable trait variable

PhyloSig_para = function(TreeAllMatrix, dat, resample_spp){
  ## TreeAllMatrix are the species tree after PruneTree function
  ## dat is a list of trait value of different species
  
  
  dist = cophenetic.phylo(TreeAllMatrix)
  phyloposi = isoMDS(dist, trace = F) %>% as.data.frame()
  
  phyloposi_species = phyloposi %>% mutate(Species = row.names(phyloposi)) %>% 
    mutate(phy1 = round(scale(points.1),digits = 1)) %>% 
    mutate(phy2 = round(scale(points.2), digits = 1))
  
  AllMatrix = merge(dat, phyloposi_species)
  
  AllMatrix_G = merge(AllMatrix, Germination)
  
  Mass.All <- AllMatrix_G$mean_Mass
  names(Mass.All) <- AllMatrix_G$Species
  
  Height.All <- AllMatrix_G$mean_Height
  names(Height.All) <- AllMatrix$Species
  
  Area.All <- AllMatrix_G$mean_Area
  names(Area.All) <- AllMatrix$Species
  
  Germination.All <- AllMatrix_G$mean_Germination
  names(Germination.All) <- AllMatrix_G$Species
  
  # calculate the phylogenetic signals
  
  phyloSigMass.k <- phylosig(TreeAllMatrix, Mass.All, method="K", test=TRUE, nsim=999)
  
  phyloSigHeight.k <- phylosig(TreeAllMatrix, Height.All, method="K", test=TRUE, nsim=999)
  
  phyloSigArea.k <- phylosig(TreeAllMatrix, Area.All, method="K", test=TRUE, nsim=999)
  
  phyloSigGermination.k <- phylosig(TreeAllMatrix, Germination.All, method="K", test=TRUE, nsim=999)
  
  
  # table display results
  tab_phyloSig = matrix(c(phyloSigMass.k$K, phyloSigHeight.k$K, phyloSigArea.k$K, 
                          phyloSigGermination.k$K, phyloSigMass.k$P, phyloSigHeight.k$P, phyloSigArea.k$P, 
                          phyloSigGermination.k$P), ncol = 2)
  
  colnames(tab_phyloSig) = c("Blomberg's K", "P-value")
  
  Variables = c("Mass", "Height", "Area", "Germination")
  
  tab_phyloSig1 = data.frame(cbind(Variables,tab_phyloSig)) 
  
  #%>% 
  #  mutate(resample_spp = resample_spp) %>%
  #  mutate(resample_rpt = resample_rpt)
  
  tab_phyloSig1
}